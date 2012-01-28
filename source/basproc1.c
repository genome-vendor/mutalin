/*------------------------------------------------------------------------
	fichier BASPROC1.C

				fonctions utilisees pour l'alignement
2 juin 1994 : dans EditSizeLine il faut Buffer > NlleTaille a cause de
	l'indice 0
-------------------------------------------------------------------------*/
/*

7 juillet 1994 : utilisation de TabArbre[0] pour AncArbre

9 fevrier 1995 : AgregeNoeuds doit tester Score!=NULL a cause du cas NbSeq=2
*/
#include "basproc1.h"

#include "msgerr.h"

int TailleGraphe= (LargeurGraphe - DebutGraphe);
/*int Phylip=0; unused */

typedef struct
{
	t_ClassHier AncArbre;
	int* Bon;
	t_indN* NoGroupe;
}* t_pData;

/*++++++++++++++++++++++++++++ PROTOTYPES ++++++++++++++++++++++++++++++++*/
static void AgregeNoeuds( 	t_indN* Groupe,t_indN NbSeq,t_DessinArb Dessin,t_ClassHier Arbre,
			t_indN Noeud1,t_indN Noeud2,t_score** Score,t_indN *Som,t_score Sup );

static void CalculScore (  t_pParamAlign Param,t_pDescripteurSequence DS,t_score **Score,
				float *ScoreTotal);

static t_score ScoreIdentity (  t_pParamAlign Param, t_pDescripteurSequence DS)
;

static void ChercheScoreMaxi(t_score** Score,t_indN Som,t_score *Sup,t_indN *Noeud1,t_indN *Noeud2 );

static void Depart( t_indN NbSeq,	t_indN* Groupe, t_ClassHier Arbre, t_DessinArb Dessin,
		t_score *Sup);

void Sc2a2( t_pParamAlign Param,t_pLigneSeq Sequence,
	t_score *x ,t_indN Prot,t_indN Prot2 );

void SupprimeGroupe( t_ClassHier Arbre,t_indN* Groupe,t_indN Groupe1,t_indN Groupe2 );

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Si on a peur que Div deborde, utilser un float pour A */
#define Div(A,B,C) (((A)*(B) + ((C)>>1)) /(C))



/*=======================================================================*/
/*unused
int maketree(t_score** Score,	t_ClassHier Arbre, t_DessinArb *Dessin,
	t_indN NbSeq, int NJ);
*/
int ArbreAlignement(t_pParamAlign Param,
				t_pDescripteurSequence DS,
				t_ClassHier Arbre,
				t_DessinArb *Dessin,
				t_indN NbSeq)
{
	int OK;
	t_indN i;
	t_score** Score, **sz;
	float ScoreTotal;
	if (Score=(t_score**)calloc(NbSeq-1,sizeof(t_score*)))
	for (i = 1, sz = Score; i < NbSeq; i++, sz++)
			if (!(*sz=(t_score*)calloc(i,sizeof(t_score))))
				break;

	memset( Arbre, 0, (NbSeq+1)*sizeof(t_Branche));
	memset (*Dessin, 0,(NbSeq+1)*sizeof(t_Dessin));

	if( OK= (Score && (i==NbSeq)))
	{

		CalculScore( Param, DS, Score, &ScoreTotal);
/* unused
#ifdef _domainer
		if (Phylip) OK = maketree(Score,Arbre, Dessin, NbSeq,Phylip-1);
		else
#endif
*/
			OK = CalculArbre(Score, Arbre, *Dessin, NbSeq,ScoreIdentity(Param,DS));
	}

	if (Score)
	{
		for (i=1, sz = Score; (i<NbSeq)&&(sz); i++, sz++)
		 free (*sz);
		free(Score);
	}
	return OK;
}

/*----------------------------------------------------------------------*/

static void Depart( t_indN NbSeq,  			/* in */
		t_indN* Groupe,     	        /* out */
		t_ClassHier Arbre,  		/* out */
		t_DessinArb Dessin,		     	/* out */
		t_score *Sup)               	/* out */
/* initialisation des groupes, de Bon, de Sup */
{
	t_indN i, *g;
	t_Dessin *d;

	*Sup= 0;

	for( i=0, g=Groupe; i<=NbSeq; i++,g++ ) *g=1;

	memset(Arbre, 0, (NbSeq+1)*sizeof(t_Branche));
	memset (Dessin,0,(NbSeq+1)*sizeof(t_Dessin));
	/* initialiser Noeud pour y ranger les scores */ /* modif */
		for (i=1,d=Dessin+1; i<=NbSeq; i++,d++) d->Noeud = i;

}

/*----------------------------------------------------------------------*/
static void ChercheScoreMaxi( 	t_score** Score,		/* in */
			t_indN Som,		/* in */
			t_score *Sup,		/* out */
			t_indN *Noeud1,		/* out */
			t_indN *Noeud2 )		/* out */
{
	t_indN i, j;
	t_score *k,**sz;

	if( Score )
		*Sup= **Score;
	else
		*Sup= SCO_MIN;
	*Noeud1= 1;
	*Noeud2= 2;

	for( i= 2, sz= Score+1; i< Som; i++,sz++ )
		for( j= 0, k = *sz; j<i; j++, k++ )
			if( *k > *Sup)
			{
				*Sup= *k;
				*Noeud1= j+1;
				*Noeud2= i+1;
			}
}		/* de ChercheScoreMaxi */

/*----------------------------------------------------------------------*/

static void Distance( t_score* Score1, t_score Score2,
			t_indN Poids1, t_indN Poids2, t_indN PoidsT )
{
	*Score1 = (Poids1 * (*Score1) + Poids2 * Score2 + (PoidsT>>1))/ PoidsT;
}


/*----------------------------------------------------------------------*/

void SupprimeGroupe( t_ClassHier Arbre,
			t_indN* Groupe,
			t_indN Groupe1,
			t_indN Groupe2 )
/* Accroche Groupe2 a Groupe1 et le supprime
  Entree : Arbre et Groupe
		 contiennent les infos sur les groupes deja formes
		 NbSeq est le nombre de sequences

  Sortie : Arbre et Groupe ont avances d'un niveau
		 dans la hierarchie

  Arbre.SeqSuiv peut etre a complete (SeqSuiv <> 0 en fin de Groupe). */
{
	t_indN Seq1,i;

	Seq1= Groupe1;
	i= 1;
	while(i < Groupe[Groupe1])
	{
		Seq1= Arbre[Seq1].SeqSuiv;
		i++;
	}
	Arbre[Seq1].SeqSuiv= Groupe2;
	Arbre[Groupe2].Pere= Groupe1;
	Groupe[Groupe1]+= Groupe[Groupe2];
	Groupe[Groupe2]= 0;
} 		/* de SupprimeGroupe */

/*----------------------------------------------------------------------*/
static void AgregeNoeuds( 	t_indN* Groupe, t_indN NbSeq,
			t_DessinArb Dessin,
			t_ClassHier Arbre,
			t_indN Noeud1,
			t_indN Noeud2,
			t_score** Score,
			t_indN *Som,
			t_score Sup )


{
	t_indN Poids1, Poids2, PoidsT;
	t_score *Som1, *Som2;
	t_indN Groupe1, Groupe2,i;

	/* le Noeud 2 doit etre le plus grand */
	if( Dessin[Noeud1].Noeud > Dessin[Noeud2].Noeud )
	{
		i= Noeud1;
		Noeud1= Noeud2;
		Noeud2= i;
	}
	Groupe1= Dessin[Noeud1].Noeud;
	Groupe2= Dessin[Noeud2].Noeud;
	Poids1= Groupe[Groupe1];
	Poids2= Groupe[Groupe2];
	PoidsT= Poids1 + Poids2;

	if (Score)
	{
		for( i=1; i<=*Som; i++ )
			if( (i!=Noeud1) && (i!=Noeud2) )
				Distance( Place (i,Noeud1,Score) , *Place(i,Noeud2,Score),
						Poids1, Poids2,PoidsT);

		for( i=1, Som1=Score[*Som-2],Som2=Score[Noeud2-2]; i< Noeud2; i++,Som1++,Som2++ )
			*Som2 = *Som1;

		for( i=Noeud2-1, Som1=Score[*Som-2]+i+1; i<*Som-2; i++, Som1++ )
			Score[i][Noeud2-1]= *Som1;

		Score[NbSeq-2][*Som-2]= Sup;  	/* on conserve des renseignements sur cette */
		i= Dessin[*Som].Noeud;   		/* agregation pour la sortie */
		Dessin[*Som].Noeud= Dessin[Noeud2].Noeud;
		Dessin[Noeud2].Noeud= i;
	}

	/* suppression du groupe2 */
	SupprimeGroupe( Arbre, Groupe, Groupe1, Groupe2);
	(*Som)--;
} 		/* de AgregeNoeuds */



int CalculArbre( t_score** Score, 		/* scores puis resultats */
		t_ClassHier Arbre,	/* nouvel arbre */
		t_DessinArb Dessin,
		t_indN NbSeq,
		t_score Maxi )
/* Entree : Score pointe vers le tableau des scores 2 a 2
		 NbSeq est le nombre de sequences
		 Maxi est une valeur maximale du score (optionel)
	Sortie : Arbre  et Dessin decrivent le nouvel arbre	*/
{
	t_score Sup, SupSup, InfSup, *s, **ss, Correction;
	t_indN* Groupe;       		/* Nb de sequences ds le groupe */
	t_indN Noeud1, Noeud2 ;	  	/* numero des noeuds regroupes */
	t_indN Som;              	/* nombre de groupes restant a regrouper */
	t_DessinArb d;

	/* allocation du tableau Groupe */
	if( !(Groupe= (t_indN*)malloc((NbSeq+1)*sizeof(t_indN))) )
	{
		TraiteErr(8);
		return faux;
	}
	Correction = 0;
	if (Score)
	{
		for (Noeud1=1,ss=Score;Noeud1<NbSeq;Noeud1++,ss++)
			for (Noeud2=0,s=*ss;Noeud2<Noeud1;Noeud2++,s++)
				if (*s<Correction) Correction=*s;
		if (Correction < 0)	for (Noeud1=1,ss=Score;Noeud1<NbSeq;Noeud1++,ss++)
			for (Noeud2=0,s=*ss;Noeud2<Noeud1;Noeud2++,s++) *s-=Correction;
	}

	Depart(NbSeq, Groupe, Arbre, Dessin, &Sup);

	Som= NbSeq;
	while( Som > 1)
	{
		ChercheScoreMaxi( Score, Som, &Sup, &Noeud1, &Noeud2);
		if (Som==NbSeq) SupSup = Sup;
		InfSup = Sup;
		AgregeNoeuds( Groupe,NbSeq,Dessin, Arbre, Noeud1, Noeud2,
					Score, &Som, Sup);
	}
	Arbre[0].Pere= 1;
	Arbre[0].SeqSuiv= 1;
	InfSup = SupSup - InfSup;
	if (!InfSup) InfSup = 1;
	if (TailleGraphe==0) TailleGraphe=InfSup=1;
	Dessin[1].Col=DebutGraphe+SupSup+Correction;
	if (NbSeq==2) Dessin[2].Col = DebutGraphe;
	else if (Score) for (Som=1, s=Score[NbSeq-2],d=Dessin+2; Som <NbSeq; Som++,s++,d++)
		d->Col = DebutGraphe +
			TailleGraphe * (SupSup - *s)/ InfSup;
	if (Maxi >= SupSup+Correction)
	Dessin->Col = DebutGraphe + TailleGraphe * (SupSup - Maxi+Correction) / InfSup; /* si les
	scores sont bornes, cette valeur permet de placer le Maxi  01/07/95*/
	else Dessin->Col = DebutGraphe - 1;
	free( Groupe );

	return vrai;
} 		/* de CalculArbre */



/*========================================================================*/
void Sc2a2( t_pParamAlign Param,
			t_pLigneSeq Sequence,
			t_score *x ,
			t_indN Prot,
			t_indN Prot2 )
{
	int Gap1, Gap2;
	int Res1, Res2;
	t_indL Debut, Fin, i;
	t_Seq *s1, *s2;

/* 	attention EnsLet - EnsSymb contient des symboles non comparables
	il faut les sauter */
	Debut = max(Sequence[Prot].PremCol, Sequence[Prot2].PremCol );

	Fin = min (Sequence[Prot].DerCol, Sequence[Prot2].DerCol );

	*x= 0;
	Gap1= faux;
	Gap2= faux;
	for( i=Debut,s1=Sequence[Prot].Seq+i,s2=Sequence[Prot2].Seq+i; i<=Fin;
		i++,s1++,s2++ )
	{
		Res1= s1->Symb;
		Res2= s2->Symb;
		if(Res1 >= Param->PremLettre)
		{
			Gap1= faux;
			if( Res2 >= Param->PremLettre )
			{
				Gap2= faux;
				*x+= Param->Coeff[Res1][Res2];
			}
			else
			{
				*x -= Param->Gap2;
				if( !Gap2 ) *x-= Param->Gap;
				Gap2= vrai;
			}
		}
		else
			if(Res2 >= Param->PremLettre)
			{
				Gap2= faux;
				if( !Gap1 )	*x-= Param->Gap;
				*x -= Param->Gap2;
				Gap1= vrai;
			}
	}
	if ((Param->GapExtT >>1) && (Debut > 1) && (Debut<=Fin))
		*x-= Param->Gap + (Debut-1)*Param->Gap2;

	if ((Param->GapExtT % 2)
		&& ((Fin=max(Sequence[Prot].DerCol, Sequence[Prot2].DerCol )- Fin) > 0))
		*x-= Param->Gap + Fin * Param->Gap2;

/* Normalisation */
	if ((Param->ScMethod > ABSOLU) && (Param->ScMethod < NORMALISED))
		*x = float_to_score(1000.0 * *x/
			min(Sequence[Prot].NbResid,Sequence[Prot2].NbResid));

}		/* de Sc2a2 */

void Id2a2  (t_pParamAlign Param,
			t_pLigneSeq Sequence,
			t_score *x ,
			t_indN Prot,
			t_indN Prot2 )
{
	t_indL Debut, Fin, i,Count = 0;
	t_Seq *s1, *s2;
	byte PremLettre = Param->PremLettre;

	Debut = max(Sequence[Prot].PremCol, Sequence[Prot2].PremCol );

	Fin = min (Sequence[Prot].DerCol, Sequence[Prot2].DerCol );

	for( i=Debut,s1=Sequence[Prot].Seq+i,s2=Sequence[Prot2].Seq+i; i<=Fin;
		i++,s1++,s2++ )
		if ((s1->Symb == s2->Symb) && (s1->Symb>=PremLettre)) Count ++;
	*x = float_to_score (1000.0 * Count /
		min (Sequence[Prot].NbResid,Sequence[Prot2].NbResid));
}
/* unused
void pamdist(t_pParamAlign Param,t_pLigneSeq Sequence,t_score *x,
	t_indN Prot,t_indN Prot2);
*/
static void CalculScore (  t_pParamAlign Param,
				t_pDescripteurSequence DS,
				t_score** Score,
				float *ScoreTotal)
{

/* Calcule les scores des alignements 2 a 2 contenus dans un align. multiple
  Entree : Sequence contient les sequences alignees sous forme Byte
		 Score pointe sur un tableau de 0 ou on pourra ranger les scores
		 NbSeq est le nombre de sequences
		 [PremLettre..DerLettre] symboles possible
		 Coeff est le tableau des comparaisons de symboles
		 Gap est la penalite pour un trou

  Sortie : Score^[(i-1)*(i-2) div 2 + j] = Score des sequences i et j (j<i)
		 ScoreTotal = somme de tous les scores/NbSeq/pred(NbSeq)
*/
/*	FILE *f = fopen("scores","wt");*/
/*	FILE *f = fopen("scores","rt");*/
	t_score ScTot, sci, scj;
	t_indN Prot, Prot2;
	t_score *sc, **sz;
	float Norm;
	void (*Score2a2)(t_pParamAlign Param,
			t_pLigneSeq Sequence,
			t_score *x ,
			t_indN Prot,
			t_indN Prot2);

	if( DS->NbSeq < 2 )
	{
		*ScoreTotal= 0;
		exit(1);
	}

	*ScoreTotal= 0;
	ScTot= 0;
	Norm= 2.0/DS->NbSeq/(DS->NbSeq-1);
	switch (Param->ScMethod)
	{
		case IDENTITY: Score2a2 = Id2a2;break;
/* unused
#ifdef _domainer
		case 4: Score2a2 = pamdist; break;
#endif
*/
		default: Score2a2 = Sc2a2;
	}

	for( Prot=2,sz=Score; Prot<=DS->NbSeq; Prot++,sz++ )
	{
		if (Param->ScMethod==3) Sc2a2(Param,DS->Sequence,&sci,Prot,Prot);
		for( Prot2=1,sc=*sz; Prot2<Prot; Prot2++,sc++)
		{
			Score2a2(Param, DS->Sequence, sc, Prot, Prot2 );
			if (Param->ScMethod==3)
			{
				Sc2a2(Param,DS->Sequence,&scj,Prot2,Prot2);
				*sc = float_to_score (2000.0* *sc /(sci+scj));
			}
			ScTot= ScTot +  *sc;
/*			if (f) {fprintf(f,"%6d",*sc);fflush(f);}*/
/*			if (f) {fscanf(f,"%6d",sc);}*/
		}
/*		if (f) fputc('\n',f);*/
/*		if (f) fscanf(f,"%*[^\n]");*/
	}
	*ScoreTotal= Norm*ScTot;
/*	if (f) fclose(f);*/
} 		/* de calculScore */

static t_score ScoreIdentity (  t_pParamAlign Param, t_pDescripteurSequence DS)
{
	t_indN Prot;
	t_score sc, scmax=0;

	for (Prot = 1; Prot<=DS->NbSeq;Prot++)
	{
		Sc2a2(Param, DS->Sequence, &sc, Prot, Prot);/*/
		sc = Id2a2 (DS->Sequence,Prot,Prot,Param->PremLettre); */
		if (sc > scmax) scmax = sc;
	}
	return scmax;
}

int EditSizeLine(t_pLigneSeq p, t_indL NlleTaille )
/* realloue la chaine de caractere p->Seq a la taille NlleTaille +1.

	Retourne faux:	-si NlleTaille n'est pas valide
			-ou si la reallocation est impossible
		Dans ces 2 cas, p->Seq est inchange.

	Retourne vrai:	-si la reallocation n'est pas necessaire
			-ou si la reallocation a ete realisee( dans ce cas,
				p->Seq a pu etre modifie */
{
	t_indL Buffer;
	t_Seq* Ancien;

	if( p->BuffLen > NlleTaille ) /*02/06/94 > et non >= */
		return vrai;	/* la zone memoire reservee est deja assez longue */

	if( (NlleTaille<0) || (NlleTaille>LongMax ))
		return faux;	/* erreur */

	Buffer= (NlleTaille|63)+1;/* allocation par blocs de 32 */
	/*02/06/94 NlleTaille et non NlleTaille -1*/
	/* on doit allonger la sequence */
	if (Ancien= p->Seq) p->Seq=(t_pSeq)realloc(p->Seq,Buffer*sizeof(t_Seq));
	else p->Seq=(t_pSeq)malloc(Buffer*sizeof(t_Seq));
	if (!(p->Seq))
	{
		p->Seq= Ancien;	/* on restaure l'ancien pointeur */
		return faux;
	}

	/* initialiser a 0 la zone ajoutee */
	memset( p->Seq+p->BuffLen, 0,
		(Buffer - p->BuffLen)*sizeof(t_Seq));

	/* maj de BufffLen et Taille de la sequence */
	p->BuffLen= Buffer;

	return vrai;
} 		/*	de EditSizeline 	*/

/*=========================================================================*/

void MajTailleConsens(t_pLigneSeq Sequence, t_indN NbSeq, char complet )
/* Alloue ou realloue la sequence 0 avec la taille de
	la plus longue sequence */
/* 1/7/94 if (complet) all sequences are edited to the same length */
{
	t_indL Long;
	t_indN Seq,i;
	int libere;
	t_pLigneSeq tmp;
	t_Seq *Ancien;

	Long= 0;
	/* on recherche la taille de la plus longue sequence */
	for( Seq=0, tmp=Sequence+1; Seq<NbSeq; Seq++,tmp++ )
		if( tmp->DerCol > Long ) Long= tmp->DerCol;
	libere = (Sequence[0].Taille > Long);

/*	if( !Long )
	{
		TraiteErr( 8, "All sequences have null length" );
		return faux;
	}
*/
	Sequence[0].DerCol= Long;
	Sequence[0].Taille= Long;

	if (!complet) return ;/*1/7/94*/

	if (libere)
	{
		Long = (Long | 63) +1;
		for( i=0,tmp=Sequence; i<= NbSeq; i++,tmp++ )
			if( tmp->BuffLen > Long)
			{
				Ancien = tmp->Seq;
				if ((tmp->Seq=realloc (tmp->Seq, Long*sizeof(t_Seq)))==NULL)
					tmp->Seq=Ancien;
				else tmp->BuffLen = Long;
			}
	}
	if (!EditSizeLine( Sequence, Long ) )
  {
		free (Sequence[0].Seq);
		Sequence[0].Seq=NULL;
		Sequence[0].BuffLen=0;
  }
}


t_score *Place(t_indN i,t_indN j,t_score**Score)
{
	if( i < j )
		return (Score[j-2] + i -1);
	else
		return (Score[i-2] + j - 1);
}

void MajWeightConsens (t_pLigneSeq Sequence, t_indN NbSeq)
{
	t_Weight w;
	t_indN i;
	t_pLigneSeq tmp;

	for (i=0, w=0, tmp = Sequence+1; i<NbSeq; i++,tmp++) w+= tmp->Weight;
	if (w==0) {
		for (i=0, tmp=Sequence+1;i<NbSeq;i++,tmp++) tmp->Weight=1;
		Sequence->Weight = NbSeq;
	}
	else Sequence->Weight = w;
}

/*=======================================================================*/
void TraducByteToLet( t_pDescripteurSequence DS, t_pParamAlign Param )
/* mettre les sequences sous forme lettre */
{
	t_pLigneSeq tmp;
	t_indN n;
	t_indL i;
	t_Seq *pos;

	for( n= 0, tmp=DS->Sequence+1; n< DS->NbSeq; n++,tmp++ )
		for( i= 0,pos=tmp->Seq+1; i<tmp->Taille; i++,pos++ )
				pos->Car= Param->Symb[pos->Symb];
}


void TraducLetToByte( t_pDescripteurSequence DS, t_pParamAlign Param )
/* mettre les sequences sous forme symbole */
{
	t_pLigneSeq tmp;
	t_indN n;
	int j;
	t_indL i;
	t_Seq *pos;

	for( n=0,tmp=DS->Sequence+1 ; n<DS->NbSeq; n++,tmp++ )
	{
		for( i=0,pos=tmp->Seq+1; i<tmp->Taille; i++,pos++ )
		{
			if( (j= Param->NumSymb[pos->Car]) > Param->DerLettre )
				pos->Symb= Param->PremLettre - 1;
			else
				pos->Symb= j;
		}
	}
}


int TraiteBranche(  t_ClassHier Arbre, t_indN NbSeq,
				t_indN *Groupe ,
				void* V,    /* pointeur sur n'importe quoi */
				ProcAction Action)
{


	t_indN PremPile, Courant, Pere, PPere;

	PremPile= Arbre[0].SeqSuiv;
	Courant= Arbre[PremPile].SeqSuiv;

	while( PremPile != 0 )
	{
		Pere= Arbre[Courant].Pere;
		while(Pere != PremPile)
		{
			PPere= Arbre[PremPile].Pere;
			if (!Action (PPere,PremPile,NbSeq,Groupe,Arbre,V)) return faux;
			PremPile= PPere;
		}
		PremPile= Courant;
		Courant= Arbre[Courant].SeqSuiv;
	}
	return vrai;
}
/*=========================================================================*/

/*=======================================================================*/






/*=======================================================================*/

