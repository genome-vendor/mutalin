#include "commande.h"
#include "basproc1.h"
#include "ma.h"

/* Sequence description methods */
void InitSequenceDescription (t_pDescripteurSequence DS)
{
	DS->Sequence=NULL;
	DS->NoSeq=NULL;
	DS->NbSeq=0;
}

void LibereDescripteurSequence( t_pDescripteurSequence DS )
/* libere les tableaux de sequences */
{
	t_indN i;
	t_pLigneSeq tmp;

	/* libere les sequences */
	for( i=0,tmp=DS->Sequence; i<=DS->NbSeq; i++,tmp++ )
		if( tmp->Seq )
			free( tmp->Seq );

	/* libere le tableau de sequences */
	free( DS->Sequence );

	/* libere le tableau d'ordre des sequences */
	if (DS->NoSeq) free (DS->NoSeq);
}

/*-------------------------------------------------------------------------*/
/* local functions for MultipleAlignment */

static int AlloueTabArbre(t_pTabClassHier TabArbre, t_indN NbSeq )
/* allocation de NbSeq+1 branches aux arbres 0,1 et 2 */
{
	int i;

	/* allocation des 3 elements Arbre */
	memset(TabArbre,0,3*sizeof(t_ClassHier));
	for( i=0; i<3; i++ )
		if((TabArbre[i]=(t_ClassHier)calloc(NbSeq+1, sizeof(t_Branche))) ==NULL)
			return faux;
	return vrai;
}

static int AlloueTabDessin(t_DessinArb *TabDessin, t_indN NbSeq )
/* allocation de NbSeq+1 branches aux dessins 0,1 et 2 */
{
	int i;

	/* allocation des 3 elements Dessin */
	memset(TabDessin,0,3*sizeof(t_DessinArb));
	for( i=0; i<3; i++ )
		if((TabDessin[i]=(t_DessinArb)calloc(NbSeq+1,sizeof(t_Dessin)))==NULL )
			return faux;
	return vrai;
}

static int InitDescripteurGroupe(t_pDescripteurGroupe DG,t_indN NbSeq)
/* allocation et initialisation des groupes de sequences */
{
	t_indN i;
	/* allocation et initialisation des tableaux Groupe et No Groupe */
	DG->NoGroupe=NULL;
	if( (DG->Groupe=(t_indN*)calloc(NbSeq+1,sizeof(t_indN)))==NULL )
		return faux;
	if( (DG->NoGroupe=(t_indN*)calloc(NbSeq+1,sizeof(t_indN))) ==NULL)
		return faux;

	for(i=0; i<=NbSeq; i++)
	{
		DG->Groupe[i]= 1;
		DG->NoGroupe[i]= i;
	}

	/* initialisations */
	DG->G1=0;
	DG->G2=0;
	DG->M=0;
	DG->M1=0;
	DG->N=0;
	DG->N1=0;
	DG->Nb1=0;
	DG->Nb2=0;

	return vrai;
}

/*-------------------------------------------------------------------------*/
static void LibereTabArbre( t_pTabClassHier TabArbre )
/* liberation du tableau de 3 arbres  */
{
	int i;

	for( i=0; i<3; i++ )
		if (TabArbre[i]) free( TabArbre[i] );
}

static void LibereTabDessin( t_DessinArb *TabDessin )
/* liberation du tableau de 3 dessins */
{
	int i;

	for( i=0; i<3; i++ )
		if ( TabDessin[i] ) free( TabDessin[i] );
}
static void LibereDescripteurGroupe( t_pDescripteurGroupe DG )
/* libere les tableaux Groupe et NoGroupe */
{
	if (DG->Groupe) free( DG->Groupe );
	if (DG->NoGroupe) free( DG->NoGroupe );
}

/*-------------------------------------------------------------------------*/
static void MetOrdre(t_ClassHier Arbre, t_pDescripteurSequence DS )
/* met a jour NoSeq avec le nouvel ordre des sequences donne par Arbre */
{
	t_indN i, Seq, *n;
	for( i=0,Seq=0, n=DS->NoSeq; i<=DS->NbSeq; i++,n++,Seq= Arbre[Seq].SeqSuiv)
		*n= Seq;
}

/*-------------------------------------------------------------------------*/
#ifndef clus2dom
/* Si on a peur que Div deborde, utilser un float pour A */
#define Div(A,B,C) (((A)*(B) + ((C)>>1)) /(C))
static int SaveWeights (t_Weight **Save, t_pTabSeq Sequence, t_indN NbSeq, t_indN SeqDejaAligne)
{
	int i;
	t_Weight *w, w1, w2;
	t_pTabSeq Seq;
	if ((*Save = (t_Weight*)calloc(sizeof(t_Weight), NbSeq))==NULL) return 0;
	for (i=1,Seq = Sequence+1,w1=0;i<SeqDejaAligne;i++,Seq++) w1+=Seq->Weight;
	for (w2=0;i<=NbSeq;i++,Seq++) w2+=Seq->Weight;
	if (w1 > 0) for (i=1, Seq=Sequence+1, w=*Save;i<SeqDejaAligne;i++,w++,Seq++)
		*w= Div(Seq->Weight,100*(SeqDejaAligne-1),w1);
	if (w2 > 0) for (i=SeqDejaAligne, Seq=Sequence+i,w=*Save+i-1;i<=NbSeq;i++,w++,Seq++)
		*w= Div(Seq->Weight,100*(NbSeq-SeqDejaAligne+1),w2);
	return 1;
}

static void RestoreWeights (t_Weight *Save, t_pTabSeq Sequence, t_indN NbSeq)
{
	int i;
	t_Weight *w;
	t_pTabSeq S = Sequence+1;
	for (i=0, w=Save;i<NbSeq;i++,w++,S++) S->Weight=*w;
	MajWeightConsens(Sequence,NbSeq);
}

static void GapInsert (t_pLigneSeq NewSeq, t_pLigneSeq OldSeq)
{
	NewSeq->PremCol=OldSeq->PremCol;
	NewSeq->DerCol=OldSeq->DerCol;
	EditSizeLine(NewSeq, OldSeq->Taille);
	NewSeq->Taille=OldSeq->Taille;
	memcpy(NewSeq->Seq, OldSeq->Seq, OldSeq->BuffLen*sizeof(t_Seq));
}

static int SpecialWeights2Groups( t_indN N2, t_pDescripteurSequence DS,
			t_pParamAlign Param, t_indN *Neighbour)
{
	t_indN Seq1, Seq2, *Seq, Nb1 = N2-1, Nb2=DS->NbSeq-Nb1, N1maxi, N2maxi;
	t_Weight Fact;    /* scoring factors */
	t_score sc, maxi;
	int special = 0;
	t_pLigneSeq s1;
	if ((Nb1 > 1) && (Nb2 > 1)) return 1;
	if (Nb1==1) {*Neighbour = 1;Seq=&Seq2;}
	else if (Nb2==1) {*Neighbour = DS->NbSeq;Seq=&Seq1;}
	Fact=0;
	for( Seq1=1,s1=DS->Sequence+1, maxi=0; (Seq1 < N2); Seq1++,s1++)
	{
		for( Seq2=N2; (Seq2 <= DS->NbSeq); Seq2++)
		{
			Id2a2(Param,DS->Sequence,&sc,Seq1,Seq2);
			if (sc > maxi)
			{
				maxi = sc;
				N1maxi = Seq1;
				N2maxi = Seq2;
				if ((maxi >= 1000.0)&&(DS->Sequence[Seq1].NbResid==DS->Sequence[Seq2].NbResid))
				{
					special = 1;
					break;
				}
			}
			Fact += DS->Sequence[*Seq].Weight = Div(1.0*DS->Sequence[*Seq].Weight,1000,1000-sc);
		}
	}
	if (Nb1==1) *Neighbour = N2maxi;
	else if (Nb2==1) *Neighbour = N1maxi;
	if (special) return 2;

	if (Nb1 > 1)
	{
		if (Fact==0) for (Seq1=1;Seq1< N2;Seq1++)
			DS->Sequence[Seq1].Weight =1;
		else for (Seq1=1;Seq1< N2;Seq1++)
			DS->Sequence[Seq1].Weight = Div(1.0*DS->Sequence[Seq1].Weight,100L*Nb1,Fact);
	}
	if (Nb2 > 1)
	{
		if (Fact==0) for (Seq2=N2;Seq2<=DS->NbSeq;Seq2++)
			DS->Sequence[Seq2].Weight = 1;
		else for (Seq2=N2;Seq2<=DS->NbSeq;Seq2++)
			DS->Sequence[Seq2].Weight = Div(1.0*DS->Sequence[Seq2].Weight,100L*Nb2,Fact);
	}
	MajWeightConsens(DS->Sequence,DS->NbSeq);

	return 1;
}

static int ProfileAlignment (t_pDescripteurSequence Sequences, t_pParamAlign Param,
	int SeqDejaAligne, 	t_ClassHier TabArbre, 	t_pDescripteurGroupe DG)
{
	t_Weight *Save=NULL;
	int OK = 0;
	int special;
	t_indN Neighbour;

	if (!Align2Groups(SeqDejaAligne,TabArbre,Sequences,DG,Param))
		return 0;
	if ((SeqDejaAligne!=2)&&(SeqDejaAligne!=Sequences->NbSeq)) return 1;

	if (!SaveWeights(&Save,Sequences->Sequence,Sequences->NbSeq,SeqDejaAligne)||
	!(special=SpecialWeights2Groups (SeqDejaAligne,Sequences,Param,&Neighbour)) ||
	((special<2)&&!Align2Groups(SeqDejaAligne,TabArbre,Sequences,DG,Param)))
		goto fin;
	else if ((SeqDejaAligne==2)&&(Neighbour > 1))
	{
		if (special == 2) GapInsert(&Sequences->Sequence[1],&Sequences->Sequence[Neighbour]);
		PESetOutputOrder();
		TabArbre[0].SeqSuiv = TabArbre[2].Pere = 2;
		TabArbre[1].Pere = Neighbour;
		TabArbre[1].SeqSuiv = TabArbre[Neighbour].SeqSuiv;
		TabArbre[Neighbour].SeqSuiv = 1;
	}
	else if ((SeqDejaAligne==Sequences->NbSeq)&&(Neighbour < Sequences->NbSeq))
	{
		if (special == 2) GapInsert(&Sequences->Sequence[Sequences->NbSeq],&Sequences->Sequence[Neighbour]);
		PESetOutputOrder();
		TabArbre[Sequences->NbSeq-1].SeqSuiv=0;
		TabArbre[Sequences->NbSeq].Pere = Neighbour;
		TabArbre[Sequences->NbSeq].SeqSuiv = TabArbre[Neighbour].SeqSuiv;
		TabArbre[Neighbour].SeqSuiv = Sequences->NbSeq;
	}
	RestoreWeights(Save,Sequences->Sequence,Sequences->NbSeq);
	OK = 1;
fin:
	if (Save) free(Save);
	return OK;
}
#endif
/*----------------------------------------------------------------------------*/

static void InitArbre1( t_ClassHier Arbre, t_indN NbSeq )
{
	register t_indN i;
	t_ClassHier A;

	for( i=1, A= Arbre; i<=NbSeq; i++,A++ )
	{
		A->Pere=0;
		A->SeqSuiv=i;
	}
	A->Pere=0;
	A->SeqSuiv=0;
}

int MultipleAlignment (t_pDescripteurSequence Sequences, t_pParamAlign Param,
	int SeqDejaAligne)
{
	t_TabClassHier TabArbre;
	t_DessinArb TabDessin[3];
	t_DescripteurGroupe DG;
	int OK=0;
	int special;
#ifndef clus2dom
	int drawtree = PEGetFormatArbre();
#endif

/*-- allocation et initialisation d'un tableau de 3 arbres ---------*/
	if( !AlloueTabArbre(TabArbre, Sequences->NbSeq) ) goto fin;
	if( !AlloueTabDessin(TabDessin, Sequences->NbSeq )) goto fin;


/*-- allocation et initialisation des tableaux Groupe et NoGroupe ----*/
	if( !InitDescripteurGroupe(&DG,Sequences->NbSeq) ) goto fin;

	InitArbre1( TabArbre[1], Sequences->NbSeq );

	if (SeqDejaAligne == 1)
	{
		TraducLetToByte (Sequences, Param);
		OK = ArbreAlignement(Param,Sequences,TabArbre[1],TabDessin+1,
			Sequences->NbSeq);
		TraducByteToLet (Sequences,Param);
		if (!OK) goto fin;
/* unused
		if (Phylip > 1)
		{
			updateweight(Sequences->Sequence,TabArbre[1],TabDessin[1],
				Sequences->NbSeq);
			MajWeightConsens (Sequences->Sequence,Sequences->NbSeq);
		}
*/
		goto affichseul;
	}
#ifndef clus2dom
	else if (SeqDejaAligne)
	{
		TraducLetToByte (Sequences, Param);
		OK = ProfileAlignment(Sequences,Param,SeqDejaAligne,TabArbre[1],&DG);
		TraducByteToLet(Sequences,Param);
		if (!OK) goto fin;
		MajWeightConsens(Sequences->Sequence,Sequences->NbSeq);
		OK = 1;
		goto affichseul;
	}
#endif
/*-- Clustering */
	if (! (special = GetArbre (Sequences,TabArbre[2], TabDessin[2]))) goto fin;
	else if (special < 0)/* on doit recalculer l'arbre */
	{
		if  (Sequences->NbSeq==2)
			OK = CalculArbre( NULL, TabArbre[2], TabDessin[2],Sequences->NbSeq, 0);
		else
		{
			OK = FastAction2( Sequences, Param, TabArbre[2], TabDessin[2] );
#ifndef clus2dom
			if (OK) SauvArbre( TabArbre[2], TabDessin[2],Sequences,".clu");
#endif
		}
		if (!OK) goto fin;
		OK = 0;
	}
/*-- alignement */
	OK = AligneAction(TabArbre,TabDessin,Sequences,Param,&DG);
	MajWeightConsens(Sequences->Sequence, Sequences->NbSeq);

affichseul:
	if(PEIsOutputOrder() &&
		((Sequences->NoSeq=(t_indN*)calloc(Sequences->NbSeq+1,sizeof(t_indN)))
			!=NULL )) MetOrdre( TabArbre[1], Sequences );

#ifndef clus2dom
	if (!SeqDejaAligne || drawtree)
		SauvArbre(TabArbre[1], TabDessin[1], Sequences, ".cl2");
#endif

	OK = 1;

/* liberer tous les tableaux */
fin:
	LibereDescripteurGroupe( &DG );
	LibereTabArbre( TabArbre );
	LibereTabDessin (TabDessin);
	return OK;
}

