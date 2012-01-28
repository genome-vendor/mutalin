/*-------------------------------------------------------------------------
	ALIGNE.C
					Alignement multiple des sequences

mai-juillet 1993 :
	traduction du programme pascal de Florence Corpet par Pascal Delobbe
mai-juillet 1994 :
	corrections et ameliorations de Florence Corpet
derniere MAJ avril 1995
-------------------------------------------------------------------------------
	CORRECTIONS

2 juin 1994 : as scores can be negative, values must be initialized to - infini

2 juin 1994 : size of Adresse is Nb2+1 and not MaxPointMu

2 juin 1994 : TabTotCoeff is a static variable no free!

1 juillet 1994 : in NettoyerSequences :
	Prot and Nb were inverted
	loops were badly nested

1 juillet 1994 : MajTailleConsens has a third parameter "complet" that must be
	false when the sequences are not yet aligned

1 juillet 1994 : in CalculeCoeff, initial. of TabTotCoeff from 0 instead of 1

1 juillet 1994 : NoGroupe must not be initial. after NettoyerSequences

13 fevrier 1995 : add IterMax to stop iteration when NbIter > IterMax

4 octobre 1995 : reservation des lignes de Direction au fur et a mesure,
	d'ou grands changements dans swap
-------------------------------------------------------------------------------
	OPTIMISATIONS

2 juin 1994 : line counter (ii) in Direction array is better initialized in the
	for declaration than computed inside the loop

16 juin 1994 : variables intermediaires pour les calculs de Direction[ii]+j

1 juillet 1994 : variable intermediaire tmp dans CalculeCoeff

1 juillet 1994 : les indices de Direction partent de 0 meme qd on ne swappe pas
	et, dans la descente, on range les Saut dans PMu[s-1] au lieu de PMu[s],
	cela utilise PMu[0] et evite les debordements

1 juillet 1994 : variable intermediaire DerLigne pour savoir quand on swappe a
	la descente

7 juillet 1994 : utilisation de TabArbre[0] pour AncArbre

7 juillet 1994 : variable intermediaire PremLigne a la place de xx (Descendre1)

-------------------------------------------------------------------------------
	AMELIORATIONS

9 mai 1994 : possibilite d'aligner des sequences partielles :
	si DG->G1 ==-1, on supprime l'initialisation de M, M1, N, N1
	dans ActAligner, et l'alignement se fait de M-M1 a M et de N-N1 a N.
	En sortie, DG->M et DG->N contiennent les valeurs dans la nouvelle
	numerotation des premiers residus non traites (nouvelle taille +1 pour un
	alignement total).

10 mai 1994 : passage de int Coeff a t_score Coeff

2 juin 1994 : Gap can be fonction of gap length
	GapLocal = Nb1*Nb2*(Gap1 + length * gap2)
10 avril 1995 : the similarity coefficient at a position is still the mean of
	all pairwise coefficients at this position, BUT only the sequences for which
	the position is internal are counted.
8 decembre 1998 : Scores are rescaled so that SCO_MAX is not reached
-------------------------------------------------------------------------*/


#include <string.h>
#include <math.h>

#include "ma.h"

#include "basproc1.h"
#include "swapd.h"
#include "msgerr.h"
#include "afichseq.h"

#define IterMax 10

static const char savename[] = "sauve.tmp";

static int ToutBlanc(	t_indL Col,t_indN  Nb,t_indN Prot,t_pTabSeq TabSequence,
					t_ClassHier Arbre,t_pParamAlign ParamAlign );
					/*1/7/94 Nb before Prot */

static int NettoyerSequences( t_pDescripteurSequence DS,t_indN* Groupe,
				t_ClassHier Arbre,t_pParamAlign ParamAlign );

static t_score Sup( t_score X1,t_score X2, t_score X3, t_indL P1, t_indL P2, t_indL *P);

static void CalculeCoeff(t_pLigneSeq Sequence,t_pParamAlign Param,t_ClassHier Arbre,
				t_score* TabTotCoeff,t_indL Rang,t_indN Nb1,t_indN Groupe1,t_Weight WeightM);
/*10/05/94 t_score*/

typedef void Ajoutefunc
(t_score Score,t_indN Nb2,t_score* TabTotCoeff,t_Seq** Adresse,
	t_Weight** Weight,t_Weight TWeight, t_score *Res);

static Ajoutefunc Ajoute, Ajoute2;

static void Ajoute (t_score Score,t_indN Nb2,t_score* TabTotCoeff,t_Seq** Adresse,
	t_Weight** Weight,t_Weight TWeight, t_score *Res);
							/*10/05/94 t_score*/
static void CalculMu (t_score DGap, t_score X1, t_score *X2, t_indL *P );
/*02/06/94 Gap fonc(gap length)*/


static int Monter( t_indN Seq21,t_pLigneSeq Sequence,t_ClassHier Arbre,
			t_pParamAlign Param,t_pGrdTab GrdTab,t_pDescripteurGroupe DG,
			t_Seq** Adresse, t_Weight** Profile, t_Weight* Weight,/*t_Weight* RDGap,*/
			t_score*TabTotCoeff,t_indL** Direction,
			t_indL* Saut,
			t_indL* NbEnrg);

static int Monter2( t_indN Seq21,t_pLigneSeq Sequence,t_ClassHier Arbre,
			t_pParamAlign Param,t_pGrdTab GrdTab,t_pDescripteurGroupe DG,
			t_Seq** Adresse, t_Weight* Weight,
			t_score*TabTotCoeff,t_indL** Direction,
			t_indL* Saut,
			t_indL* NbEnrg);

static int Descendre1( t_indL** Direction,t_indL Saut,
			t_pDescripteurGroupe DG,t_pGrdTab GrdTab, t_indL NbEnrg);

static int MettreGap (t_indL x,t_indL y,t_indL *DerPos,unsigned Ext,
			t_pDescripteurGroupe DG,t_pLigneSeq Sequence,t_ClassHier Arbre );

static int Descendre2( t_pDescripteurGroupe DG,t_pLigneSeq Sequence,
		t_ClassHier Arbre,t_pGrdTab GrdTab,unsigned Ext,t_indL Saut );

static t_indL AlloueTabDirection( t_indL*** Direction, t_indL NbL, t_indL NbC );

static void LibererDirection( t_indL** Direction, t_indL NbL );

static int AlloueTabLong( t_score* *pTab, t_indL taille );

static int AlloueTabInt( t_indL* *pTab, t_indL taille );

static int AlloueTabPointSeq( t_Seq** *pTab, t_indL taille);

static int Initialiser(int Gros, t_pDescripteurGroupe DG,t_pGrdTab GrdTab,
			t_Seq** *pAdresse, t_Weight** *pProfile, t_Weight* *pWeight,
			/*t_Weight* *pRealWeight, */	t_indL** *pDirection);

int Aligner (t_pDescripteurGroupe DG,t_pDescripteurSequence DS,
			t_ClassHier Arbre,t_pGrdTab GrdTab,t_pParamAlign ParamAlign);


static int ActAligner(t_indN Groupe1,t_indN Groupe2,t_indN NbSeq, t_indN *Groupe,
	t_ClassHier Arbre, void *V);

static int Alignement ( t_ClassHier Arbre, t_DessinArb Dessin,t_pDescripteurGroupe DG,
			t_pDescripteurSequence DS, t_pParamAlign Param );

int AligneAction( t_TabClassHier TabArbre,t_DessinArb *Dessin,
		t_pDescripteurSequence DS,t_pParamAlign Param,t_pDescripteurGroupe DG );


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

static int ToutBlanc(	t_indL Col,		/* numero de colonne */
		t_indN  Nb,            	/* nombre de sequences */
		t_indN Prot,           	/* premiere sequence *//*1/7/94 Nb before Prot */
		t_pTabSeq TabSequence,
		t_ClassHier Arbre,
		t_pParamAlign ParamAlign )/* parametres d'alignement */

/*	 Regarde s'il y a 1 symbole en colonne L dans les Nb sequences commencant
  par Prot et se suivant selon Arbre
  En Entree : L no de la colonne a regarder
		 Lmax longueur des sequences
		 Nb nombre de sequences a regarder
		 Prot premiere sequence a regarder
		 Arbre doit decrire l'enchainement des sequences
  En Sortie : ToutBlanc s'il n'y a pas de symb en colonne L     */

{
	int PremLettre= ParamAlign->PremLettre;

	while( (Nb > 0) && (TabSequence[Prot].Seq[Col].Symb < PremLettre) )
	{
		Nb--;
		Prot= Arbre[Prot].SeqSuiv;
	}
	return( Nb== 0 );
}


static int NettoyerSequences( t_pDescripteurSequence DS,
				t_indN* Groupe,
				t_ClassHier Arbre,
				t_pParamAlign ParamAlign )
/* remet les sequences dans les groupes definis par Groupe
  Entree:
		TabSequence donne les sequences alignees
		Groupe donne les groupes initiaux
  Sortie:
	   TabSequence donne les sequences alignees par groupes initiaux   */
{
	register t_indN  i,j,Prot;

	t_indN Nb, Trou;
	t_indL Longueur, DerPos, Lon, Mini, Maxi, Reste;

	t_LigneSeq* tmp;

	for(i=1; i; i=Prot )
	{
		Nb= Groupe[i];	/* nbre de sequences dans le groupe i */

		tmp= DS->Sequence+i;

		Longueur = tmp->Taille;	/* Longueur de la sequence i */

		/* suppression des blancs de debut : Mini = inf PremCol
					Maxi = sup DerCol; Lon = sup NbResid */
		Mini= tmp->PremCol;
		Maxi= tmp->DerCol;
		Lon= tmp->NbResid;
		Prot= i;
		for( j=2; j<=Nb; j++)
		{
			Prot = Arbre[Prot].SeqSuiv;
			tmp= DS->Sequence+Prot;
			if( tmp->PremCol < Mini )
				Mini= tmp->PremCol;
			if( Maxi < tmp->DerCol )
				Maxi= tmp->DerCol;
			if( Lon < tmp->NbResid )
				Lon= tmp->NbResid;
		}

		/*1/7/94 this is still in the i loop */
		Prot= i;
		Mini--;
		Longueur= Maxi - Mini;
		if( Mini>0 )
			for (j=1; j<=Nb; j++ )
			{

				tmp= DS->Sequence+Prot;
				tmp->DerCol-= Mini;
				memmove( tmp->Seq+1, tmp->Seq+Mini+1,
						Longueur*sizeof(t_Seq));
				tmp->PremCol-= Mini;
				Prot= Arbre[Prot].SeqSuiv;
			}

	/* suppression des trous communs a toutes les sequences du groupe
			Taille mini = sup NbResid = Lon */
		DerPos= 1;
		while( (Longueur > Lon) && (DerPos < Longueur) )
		{
		/* DerPos passe du premier non blanc au premier blanc */
			DerPos++;

			while( (DerPos < Longueur)
			  && !ToutBlanc (DerPos,Nb,i,DS->Sequence,Arbre,ParamAlign) )
				DerPos++;

			if( DerPos < Longueur )
			{
		/* Mini passe du premier blanc au premier non blanc */
				Mini = DerPos + 1;
				while( ToutBlanc (Mini,Nb,i,DS->Sequence,Arbre,ParamAlign) )
					Mini++;
		/* on supprime le trou */
				Trou= Mini - DerPos;
				Longueur-= Trou;
				Reste= Longueur - DerPos + 1;
				Prot= i;
				for( j=1; j<=Nb; j++ )
				{
					tmp= DS->Sequence+Prot;
					memmove( tmp->Seq+DerPos, tmp->Seq+Mini, Reste*sizeof(t_Seq));
					if( tmp->PremCol >= Mini )
						tmp->PremCol-= Trou;
					if( tmp->DerCol >= Mini )
						tmp->DerCol-= Trou;
					Prot= Arbre[Prot].SeqSuiv;
				}
			}
		}

	/* suppression des blancs de fin */
		Prot= i;
		for( j=1; j<=Nb; j++)
		{
			tmp= DS->Sequence+Prot;
			tmp->Taille= Longueur;
			Prot= Arbre[Prot].SeqSuiv;
		}

/* passage au groupe suivant */
	}

	MajTailleConsens(DS->Sequence, DS->NbSeq, faux );

	return vrai;
}


/*------------------------------------------------------------------------*/
static t_score Sup( t_score X1,t_score X2, t_score X3, t_indL P1, t_indL P2, t_indL *P)
/* calcul du sup de X1, X2, X3 */
{
	t_score X;

	if(X1 >= X2)
	{
		if(X1 >= X3)
		{
			X= X1;
			*P= 0;
		}
		else
		{
			X= X3;
			*P= P2;
		}
	}
	else
		if( (X2>X3) || ((X2==X3) && (P1<P2)) )
		{
			X= X2;
			*P= -P1;
		}
		else
		{
			X= X3;
			*P= P2;
		}
	return X;
}

/*------------------------------------------------------------------------*/
/*t_Weight Div (t_Weight A, t_Weight B, t_Weight C)
{
	if (1.0*A*B + (C >>1) > SCO_MAX) printf ("trop gros %ld %ld %ld\n",A,B,C);
	return (1.0*A*B + (C >>1) ) / C;
}*/
/* Si on a peur que Div deborde, utiliser un float pour A ou bien float_to_score*/
#define Div(A,B,C) (((A)*(B) + ((C)>>1)) /(C))

static void CalculeCoeff(t_pLigneSeq Sequence,t_pParamAlign Param,t_ClassHier Arbre,
				t_score* TabTotCoeff,t_indL Rang,t_indN Nb1,t_indN Groupe1,t_Weight WeightM)
/*10/05/94 t_score*/
/* remplit TabTotCoeff avec les coefficients correspondant a toutes
	les associations possibles */
/*20/09/94 Param[0,x]=Param[x,O]=-Gap2 */
/*10/04/95 le coeff est moyenne seulement sur les sequences non terminees (Res1!=1)*/
/* 22/10/96 un gap externe est code par 0 et non 1 */
/* 01/12/97 TabTotCoeff[0] est le poids des sequences non gap (Res1 >= PremLettre */
/* 16/07/98 TabTotCoeff[0] est pondéré par WeightM/InSeq */
/* 26/11/98 suppression des 2 dernieres modifs */
/* 08/12/98 ponderation par WeightM */
{
	t_indN Seq1,j;
	int Res1, Res2, r;
	t_score *t, *p;
	t_pLigneSeq tmp;
	t_Weight w, InSeq=0;

	for(r=0, t=TabTotCoeff; r<=Param->DerLettre; r++,t++)	*t= 0;

	for( j=1,Seq1= Groupe1,tmp=Sequence+Seq1;j<=Nb1;
				 j++,Seq1= Arbre[Seq1].SeqSuiv,tmp=Sequence+Seq1)
	{
				if ((w = tmp->Weight)<=0) w = 1;
				if ((Res1= tmp->Seq[Rang].Symb)!=0) InSeq+=w;
/*				if (Res1>=Param->PremLettre) TabTotCoeff[0]+= w;*/
				for( Res2= 1,t=TabTotCoeff+Res2,
					p=Param->Coeff[Res1]+Res2;
					Res2<= Param->DerLettre; Res2++, t++, p++ )
					*t+= *p * w;
	}
	if (WeightM!=InSeq) for(r=1, t=TabTotCoeff+1; r<=Param->DerLettre; r++,t++)
		*t = Div(*t, WeightM ,InSeq);
/*	TabTotCoeff[0]=Div(TabTotCoeff[0],WeightM,InSeq);*/
} /* de CalculeCoeff */

/*-----------------------------------------------------------------------*/
static void Ajoute (t_score Score,t_indN Nb2,t_score* TabTotCoeff,t_Seq** Adresse ,
	t_Weight** Weight,t_Weight TWeight, t_score *Res )
/* En sortie, Res recoit le nouveau score */
/*10/04/95 le coeff est moyenne seulement sur les sequences non terminees (Res!=1)*/
/* 22/10/96 un gap externe est code par 0 et non 1 */
{
	register t_indN Ns;
	t_Seq** A;
	t_score S;
	t_Weight InSeq=0;
	t_Weight *w;
	for( Ns=1,A=Adresse+1,w=*Weight,S=0; Ns<=Nb2; (*A)--,A++,Ns++,w++ )
	 if ((*A)->Symb!=0)
	{
		InSeq+= *w;
		S+= TabTotCoeff[(*A)->Symb]* *w;
	}
	S = Div (S , TWeight ,InSeq);
	*Res= Score+S;
}
/*-----------------------------------------------------------------------*/

static void Ajoute2 (t_score Score,t_indN N,t_score* TabTotCoeff,t_Seq** dummy,
	t_Weight** Profile,t_Weight TWeight, t_score *Res)
{
	register int Ns;
	t_Weight** A;
	t_score S, *C;
	t_Weight InSeq=0;
	InSeq = *(*Profile)--;
	for( Ns=1,A=Profile+1,C=TabTotCoeff+1,S=0; Ns<N; (*A)--,A++,Ns++,C++ )
		S+= **A * *C;
	S = Div(S , TWeight, InSeq);
	*Res= Score+S;
}

static t_Weight RealWeight (t_indL j, t_Weight **Profile,t_indN Seq2,t_indN Nb,
	t_Weight *w,t_pLigneSeq Sequence,t_ClassHier Arbre, int First)
{
	int Ns,c;
	t_Weight ww;
	for (Ns=0,ww=0;Ns<Nb;Ns++,w++)
	{
		if ((c=Sequence[Seq2].Seq[j].Symb)>0)
			{
/*				if (c>First) ww+=*w;*/
				*(Profile[0])+=*w;
				*(Profile[c])+=*w;
			}
		Seq2= Arbre[Seq2].SeqSuiv;
	}
	return ww;
}

static void CalculMu (t_score DGap, t_score X1, t_score *X2, t_indL *P )
/*02/06/94 Gap fonc(gap length)*/
{
	if( (*X2-=DGap) <= X1 )/*02/06/94 Gap fonc(gap length)*/

	{
		*X2= X1;
		*P= 0;
	}
	else
		(*P)++;
}

/*-----------------------------------------------------------------------*/
static int Monter( t_indN Seq21,t_pLigneSeq Sequence,t_ClassHier Arbre,
			t_pParamAlign Param,t_pGrdTab GrdTab,t_pDescripteurGroupe DG,
			t_Seq** Adresse, t_Weight** Profile, t_Weight* Weight,/*t_Weight* RDGap,*/
			t_score*TabTotCoeff,t_indL** Direction,
			t_indL* Saut,
			t_indL* NbEnrg)
{
	register t_indL j;
	register t_indL i;
	size_t LineSize;
	t_score Nu, *mu;
	t_indL PNu, *pmu;
	t_score GapLocal, DGap, DGap1, DGap2, /*DGapi,*/Temp, T1, GapExt, *T;
	t_indN Ns;
	t_indL ii, Pos;
	t_indN Seq2;
	t_indL ** Direci, *Direcij;
	t_Weight** P;
	t_Seq** A;
	char ** D, **Dc/*, c*/;
	int Dinc;
	t_Weight *w/*, *rw*/;
	int Last, First=Param->PremLettre-1;
	t_Weight WeightM, WeightN, ww;
	Ajoutefunc *Add;
	if (CtrlBreakHit) return faux;

	/* Show Message */
	AffichShow ((long)(DG->M1+1));
	Temp = min (DG->M1+1,DG->N1+1);
	T1 = SCO_MAX /Param->CoeffMax -1;
	Temp= T1/Temp -1;
/*	fprintf(stderr,"Temp=%ld T1=%ld N1=%d N2=%d ",Temp, T1,DG->Nb1,DG->Nb2);*/
	if ((Temp < 1)||(T1/DG->Nb1/DG->Nb2 <2)) {fprintf(stderr,"Temp\n");TraiteErr(9); return faux;}
/* calcul du poids total du groupe 1*/
	for ( Seq2 = DG->G1, Ns = 0, WeightM=0; Ns < DG->Nb1;
		Seq2=Arbre[Seq2].SeqSuiv,Ns++)
	{
		if ((ww =Sequence[Seq2].Weight) <=0) ww=1;
		WeightM += ww;
	}
/* calcul du poids total du groupe 2 et tableau des poids Weight*/
	for ( Seq2 = Seq21, Ns = 0, w= Weight, WeightN=0; Ns < DG->Nb2;
		Seq2 = Arbre[Seq2].SeqSuiv, Ns++, w++)
	{
		if ((*w = Sequence[Seq2].Weight)<=0) *w = 1;
		WeightN+= *w;
	}
/*	fprintf(stderr,"WM=%ld WN=%ld\n",WeightM,WeightN);*/
	WeightM = min (WeightM, T1/WeightM);
	WeightM = min (WeightM, Temp);
	WeightM = min (WeightM, T1/WeightN);
	if (WeightM < 1){fprintf(stderr,"wM\n");TraiteErr(9); return faux;}
	WeightN = min (WeightN, T1/WeightN/WeightM);
	WeightN = min (WeightN, Temp/WeightM);
	if (WeightN < 1){fprintf(stderr,"wN\n");TraiteErr(9); return faux;}
/*	WeightM=DG->Nb1;WeightN=DG->Nb2;*/
	DGap = WeightM * WeightN;
	GapLocal= DGap * Param->Gap;
	DGap *= Param->Gap2;
	DGap1 = (Param->GapExtT >> 1) ? DGap : 0;
	DGap2 = (Param->GapExtT & 1) ? DGap : 0;

/* initialisation de Mu a - l'infini */
	Temp= SCO_MIN+DGap;
	for( i=0, T=GrdTab->Mu; i<= DG->N1+1; i++,T++ ) *T= Temp;

/* ligne M : calcul du profile et du poids des sequences non vides*/
	CalculeCoeff(Sequence,Param,Arbre,TabTotCoeff,DG->M,DG->Nb1,DG->G1,WeightM);

/* colonne N: memoriser les adresses des sequences
	poids ses sequences non vides
	cellule (M,N)
*/
	j = DG->N1+1;
	T=GrdTab->X+j;
/*	rw = RDGap+j;*/
	if (Adresse)
	{
		Last = DG->Nb2;
		Add = Ajoute;
		D = (char**)(Adresse+1);
		Dinc = j*sizeof(t_Seq);
		Profile = &Weight;
		for(Seq2=Seq21, Ns= 0,A=Adresse+1/*,w=Weight,ww=0,InSeq=0*/; Ns<DG->Nb2; Ns++/*,w++*/,A++)
		{
			*A= Sequence[Seq2].Seq+DG->N;
/*			c = (*A)->Symb;
			if (c>0) {
				InSeq+=*w;
				if (c>First) ww += *w;
			}*/
			Seq2= Arbre[Seq2].SeqSuiv;
		}
/*		*rw--=Div(DGap,ww,InSeq);*/
		Ajoute(0,Last,TabTotCoeff,Adresse,&Weight,WeightN,T);
		for( j--,T--,
			GapExt = (Param->GapExtT & 1) ? GapLocal+DGap : 0;
/*			GapExt = (Param->GapExtT & 1) ? GapLocal+RDGap[j+1] : 0;*/
			j>0; j--,T--,GapExt += DGap2 /*? RDGap[j+1]:0*/)
		{
/*			for (Ns=0,A=Adresse+1,w=Weight,ww=0,InSeq=0;Ns<DG->Nb2;Ns++,w++)
			{
				c= (*A++)->Symb;
				if (c>0) {
					InSeq+=*w;
					if (c>First) ww+=*w;
				}
			}
			*rw--=Div(DGap,ww,InSeq);*/
			Ajoute(-GapExt,Last,TabTotCoeff,Adresse,&Weight,WeightN,T);
		}
	}
	else
	{
		Last=Param->DerLettre+1;
		Add=Ajoute2;
		D = (char**)Profile;
		Dinc=j*sizeof(t_Weight);
		for (Ns=0,P=Profile;Ns<Last;Ns++,P++) *P+=DG->N1+1;
		Pos=DG->N;
		T=GrdTab->X+j;
		/*ww= */RealWeight(Pos,Profile,Seq21,DG->Nb2,Weight,Sequence,
			Arbre,First);
/*		*rw=Div(DGap,ww,**Profile);*/
		Ajoute2(0,Last,TabTotCoeff,NULL,Profile,WeightN,T);
		for( j--,T--,Pos--,
				GapExt = (Param->GapExtT & 1)? GapLocal+ DGap2 : 0;
/*			 GapExt =(Param->GapExtT & 1)? GapLocal+*rw : 0,rw--;*/
			 j>0; j--,T--,Pos--,GapExt += DGap2/* ? *rw:0, rw--*/ )
		{
			/*ww= */RealWeight(Pos,Profile,Seq21,DG->Nb2,Weight,
				Sequence,Arbre,First);
/*			*rw=Div(DGap,ww,**Profile);*/
			Ajoute2(-GapExt,Last,TabTotCoeff,NULL,Profile,WeightN,T);
		}
	}

/*	GapExt = (Param->GapExtT & 1) ? GapLocal+Div(DGap,TabTotCoeff[0],WeightM) : 0;*/
	GapExt = (Param->GapExtT & 1) ? GapLocal+DGap : 0;
	LineSize = (DG->N1+1)*sizeof(t_indL);
	for( i=DG->M1,Pos=DG->M-1,ii= i-1,
	Direci= Direction+ii;
	i>0; i--,ii--,Direci--,GapExt += DGap2/*Div(DGap2,TabTotCoeff[0],WeightM)*/)
	{
		Temp= GrdTab->X[DG->N1+1];
		Nu=SCO_MIN+DGap;/* 02/06/94 scores can be < 0 */
		PNu=0;
/*		DGapi=Div(DGap,TabTotCoeff[0],WeightM);
		DGap1 = (Param->GapExtT >> 1) ? DGapi : 0;*/
		CalculeCoeff(Sequence,Param,Arbre,TabTotCoeff,Pos--,DG->Nb1,DG->G1,WeightM);
/* 4/10/95 on alloue ici les nouveaux tableaux */
		if(!*Direci && !(*Direci=(t_indL*)malloc(LineSize)))
		{
			*NbEnrg = DG->M1-1-ii;
			if (InitSwapD (*NbEnrg,DG->M1,LineSize)
			 ==SWAPIMPOSSIBLE) { TraiteErr(8); return faux; }
			AffichUpDate ((long)i+1);
			if (CtrlBreakHit) return faux;
			if (SwapD((void**)Direci+1,i)==SWAPIMPOSSIBLE)
				 { TraiteErr(8); return faux; }
			ii = *NbEnrg-1;
			Direci = Direction  + DG->M1-1;
		}

		for( Ns= 0,Dc=D; Ns<Last; Ns++,Dc++ )
			*Dc+=Dinc;/* tableau de pointeurs  */
		Add( -GapExt,Last,TabTotCoeff,Adresse,Profile,WeightN,GrdTab->X+DG->N1+1);
		/* pour chaque residu */

		for(j=DG->N1,T=GrdTab->X+j,mu=GrdTab->Mu+j,pmu=GrdTab->PMu+j,/*rw=RDGap+j+1,*/
			Direcij=*Direci+j ;
			j>0; j--,T--,mu--,pmu--,Direcij-- )
		{
			T1= *T;
/*			CalculMu(DGapi,Temp,mu,pmu);
			CalculMu(*rw--,Temp,&Nu,&PNu);*/
			CalculMu(DGap,Temp,mu,pmu);
			CalculMu(DGap,Temp,&Nu,&PNu);


			Add (Sup (Temp, *mu - GapLocal, Nu - GapLocal, *pmu, PNu, Direcij),
						Last,TabTotCoeff,Adresse,Profile,WeightN, T);
			Temp= T1;
		}

		CalculMu (DGap1, Temp, GrdTab->Mu+0, GrdTab->PMu+0);

		if ( !(ii%100))
		{
			AffichUpDate ((long)i);
			if (CtrlBreakHit) return faux;
		}

		if( ( !ii ) && (i>1))	/* on swappe */
		{
/* 4/10/95 c'est la fin de Direction qui est remplie
 Direction devient Direci (c'etait deja vrai)
 et Direci devient Direction + DG->M1*/
			if (SwapD((void**) Direci, i-1 )==SWAPIMPOSSIBLE) {
				TraiteErr (8);
				return faux;
			}
			ii =  *NbEnrg;/* 02/06/94 updating of ii */
			Direci = Direction + DG->M1;
		}
	}

	for( j=DG->N1,T=GrdTab->X+j,Nu=*(T+1),PNu=0; j>0; j--,T--  )
		CalculMu(DGap1 /*? RDGap[j]:0*/,*T,&Nu,&PNu);

/*	CalculMu(Div(DGap1,TabTotCoeff[0],WeightM),GrdTab->X[1],GrdTab->Mu+0,GrdTab->PMu+0);*/
	CalculMu(DGap1,GrdTab->X[1],GrdTab->Mu+0,GrdTab->PMu+0);

	GapExt = (Param->GapExtT >> 1) ? GapLocal : 0;
	Temp= Sup(GrdTab->X[1],GrdTab->Mu[0]-GapExt,Nu-GapExt,GrdTab->PMu[0],PNu,Saut);
	return vrai;
} 	/*de Monter*/


/*------------------------------------------------------------------------*/

static int Descendre1( t_indL** Direction,
		t_indL Saut,
		t_pDescripteurGroupe DG,
		t_pGrdTab GrdTab,
		t_indL NbEnrg )
/* si on swappe, les indices de Direction vont de 0 a NbEnrg,
	sinon de 0 a DG->M1 - 1*/
/* 4/10/95 les blocs de Swap partent de la ligne M1 au lieu de 1 */
{
	t_indL r,s,x,y;
	t_indL DerLigne = /*min (NbEnrg, DG->M1);*/
		/* 4/10/95*/   DG->M1 - ((DG->M1-1)/NbEnrg)*NbEnrg;
	t_indL PremLigne = DerLigne - NbEnrg +1;/* 4/10/95 */

	r=1;
	s=1;
	while( (r<=DG->M1) && (s<=DG->N1) )/*09/05/94 pour seq.partielles N!=N1+1*/
	{
		GrdTab->PMu[s-1]= Saut;
		if( Saut >= 0 )
		{
			x= r;
			y= s + Saut;
		}
		else
		{
			x= r - Saut;
			y= s;
		}

		/* a priori la paire suivante est sans gap */
		r= x + 1;
		s= y + 1;


		if (x > DerLigne )
		{
			/* on recalcule chaque fois le numero car on peut sauter
				des fichiers de swap */
			DerLigne = /*((x-1)/NbEnrg + 1) * NbEnrg;*/
			/* 4/10/95 */ DG->M1 - ((DG->M1-x)/NbEnrg)*NbEnrg;
			if /* 4/10/95 */(x > DG->M1)
				PremLigne = x;
			else
			{
			 PremLigne = DerLigne - NbEnrg + 1;
			AffichUpDate ((long) PremLigne);
			if (CtrlBreakHit) return faux;
			 if (/*(PremLigne <= DG->M1)
				&& */(RecupSwapD((void**)Direction, x-1 )==SWAPIMPOSSIBLE))
			{
				TraiteErr (16);
				return faux;
			};
			}

		}

		Saut= Direction[x - PremLigne][y];

	}
	GrdTab->PMu[s-1]= Saut;/*1/7/94 s-1 au lieu de s -> de 0 a N1+1 */

	return vrai;
} 		/* de Descendre1 */

/*------------------------------------------------------------------------*/

static int Insere(    t_pLigneSeq Sequence,
			t_ClassHier Arbre,
			unsigned Ext,
			t_indN Seq1,
			t_indN Nb,
			t_indL DerPos,
			t_indL Long )
{
	t_pLigneSeq tmp;
	t_indN j;
	long code= Sequence[Seq1].Taille + Long;

	if(code > LongMax)
	{
		TraiteErr( 11, Sequence[Seq1].NomSeq );
		return faux;
	}
	for( j=0; j<Nb; j++ )
	{
		tmp= Sequence+Seq1;

		if( !EditSizeLine( tmp, code ) )
		{
			TraiteErr (18);
			return faux;
		}

		memmove(tmp->Seq+DerPos+Long+1,tmp->Seq+DerPos+1,code-Long-DerPos);
		memset(tmp->Seq+DerPos+1, 1, Long*sizeof(t_Seq));
		if(DerPos < tmp->PremCol)
		{
			tmp->PremCol+= Long;  	/* blanc en tete */
			if( Ext < 2 )
				memset(tmp->Seq+DerPos+1,0,Long*sizeof(t_Seq));
		}
		if( DerPos < tmp->DerCol )
			tmp->DerCol+= Long;
		else
			if( Ext%2 == 0 )
				memset(tmp->Seq+DerPos+1,0,Long*sizeof(t_Seq));
		tmp->Taille= code;
		Seq1= Arbre[Seq1].SeqSuiv;
	}
	return vrai;
}		/* de insere */

static int MettreGap( t_indL x, t_indL y, 	/* x ou y est nul */
		t_indL *DerPos,
		unsigned Ext,
		t_pDescripteurGroupe DG,
		t_pLigneSeq Sequence,
		t_ClassHier Arbre )
{
	if( (y > 0) && !Insere (Sequence,Arbre,Ext,DG->G1, DG->Nb1, DG->M, y) )
		return faux;									/*09/05/94 M pour seq.partielles*/

	if( (x > 0) && !Insere (Sequence,Arbre,Ext,DG->G2, DG->Nb2, DG->N, x) )
		return faux;                           /*09/05/94 N pour seq.partielles*/

/*09/05/94 pour seq.partielles : updating M and N */
	x += y+1;
	*DerPos+= x;
	DG->M += x;
	DG->N += x;


	return vrai;
} 		/*de MettreGap */

/*------------------------------------------------------------------------*/

static int Descendre2( t_pDescripteurGroupe DG,
		t_pLigneSeq Sequence,
		t_ClassHier Arbre,
		t_pGrdTab GrdTab,
		unsigned Ext,
		t_indL Saut )

{
	t_indL r,s,DerPos,x,y;
/*09/05/94 pour seq.partielles : init. M and N to the first aligned position */
	DG->M -= DG->M1 + 1;
	DG->N -= DG->N1 + 1;
	DerPos = 0;
	r= 1;
	s= 1;
	while( (r<=DG->M1+1) && (s<=DG->N1+1) )
	/*09/05/94 pour seq.partielles N!=N1+1*/
	{
		Saut= GrdTab->PMu[s-1];
		if( Saut >= 0 )
		{
			x= r;
			y= s + Saut;
			if( !MettreGap(0, Saut,&DerPos,Ext,DG,Sequence,Arbre) )
				return faux;
		}
		else
		{
			x= r - Saut;
			y= s;
			if( !MettreGap(- Saut, 0,&DerPos,Ext,DG,Sequence,Arbre) )
				return faux;
		}
		/* a priori la paire suivante est sans gap */
		r= x + 1;
		s= y + 1;
	}
	return MettreGap(DG->M1-x+1,DG->N1-y+1,&DerPos,Ext,DG,Sequence,Arbre);
	/*09/05/94 pour seq.partielles N!=N1+1 */
} 		/* de Descendre2 */


/************************************************************************/

static t_indL AlloueTabDirection( t_indL*** Direction, t_indL NbL, t_indL NbC )
/* allocation d'un tableau de NbL lignes sur NbC colonnes.
Retourne SWAPIMPOSSIBLE si l'allocation minimale necessaire est impossible,
 ou l'indice de la derniere ligne allouee si l'allocation totale est impossible. */
{
/*	register t_indL i;
	t_indL** D;*/

	/* allocation d'un tableau de NbL pointeurs de ligne */
	if( ! (*Direction= (t_indL**)calloc(NbL,sizeof(t_indL *))) )
		return SWAPIMPOSSIBLE;
	/* sur chaque pointeur de ligne, allocation d'un tableau de NbC colonnes */
	/* 4/10/95 : allocation d'un seul tableau */
/*	for( i=0,D=(*Direction); i<NbL; i++,D++ )*/
		if( ! ((*Direction)[NbL-1]= (t_indL*)calloc(NbC,sizeof(t_indL))) )
			return 0;
	return NbL;		/* toute la matrice a pu etre allouee */
}


static void LibererDirection( t_indL** Direction, t_indL NbL )
/* liberation d'un tableau a 2 dimensions de NbL lignes */
{
	register t_indL i;
	t_indL** D;
	/* liberation des lignes */
	for( i=0, D=Direction; i<NbL; i++,D++ )
		if (*D) free(*D);

	/* liberation d'un tableau de pointeurs */
	free(Direction);
}
/*------------------------------------------------------------------------*/

/* allocation des grands tableaux */
static int AlloueTabLong( t_score* *pTab, t_indL taille )
{
	if(! (*pTab= (t_score *)calloc(taille,sizeof(t_score))) )
		return faux;
	return vrai;
}

static int AlloueTabInt( t_indL* *pTab, t_indL taille )
{
	if(! (*pTab= (t_indL *)calloc( taille, sizeof(t_indL))) )
		return faux;
	return vrai;
}

static int AlloueTabPointSeq( t_Seq** *pTab, t_indL taille)
{
	if(! (*pTab= (t_Seq**)calloc( taille, sizeof(t_Seq*))) )
		return faux;
	return vrai;
}

/*------------------------------------------------------------------------*/

static int AlloueProfile (t_Weight** *pProfile, int NbL, int NbC)
{
	int i;
	t_Weight **p;
	if ((*pProfile= (t_Weight**)calloc(NbL,sizeof(t_Weight*)))==NULL)
		return 0;
	for (i=0,p=*pProfile;i<NbL;i++,p++)
		if ((*p=(t_Weight*)calloc(NbC,sizeof(t_Weight)))==NULL) return 0;
	return 1;
}

static int Initialiser(int Gros, t_pDescripteurGroupe DG,t_pGrdTab GrdTab,
			t_Seq** *pAdresse, t_Weight** *pProfile, t_Weight* *pWeight,
			/*t_Weight* *pRealWeight,*/	t_indL** *pDirection)
/*allocation des tableaux de GrdTab( Mu, PMu, X), Adresse/Profile et Direction*/
{
	t_indL MaxTab;

	/* allocation des grands tableaux */
	MaxTab= (DG->N1+2);

	if( !AlloueTabLong( &GrdTab->Mu, MaxTab ) )
		return SWAPIMPOSSIBLE;

	if( !AlloueTabLong( &GrdTab->X, MaxTab ) )
		return SWAPIMPOSSIBLE;

	if( !AlloueTabInt( &GrdTab->PMu, MaxTab) )
		return SWAPIMPOSSIBLE;

	if (Gros)
	{
		if ( !AlloueProfile( pProfile, Gros+1, MaxTab) )
			return SWAPIMPOSSIBLE;
	}
	else
	{
		if( !AlloueTabPointSeq( pAdresse, DG->Nb2+1) )
			return SWAPIMPOSSIBLE;
	}

	if ( (*pWeight = calloc (DG->Nb2, sizeof(t_Weight)))==NULL)
		return SWAPIMPOSSIBLE;

/*	if ( (*pRealWeight = calloc (MaxTab, sizeof(t_Weight)))==NULL)
		return SWAPIMPOSSIBLE;*/
/* 4/10/95 Direction n'a que DG->M1 lignes */

	 return AlloueTabDirection( pDirection, DG->M1, DG->N1+1 );

}  /* de initialiser */


/*----------------------------------------------------------------------*/
static void Eclairage (int i,int N,int Col,t_ClassHier Arbre)
{
	int j;
	for (j=0;j<N;j++,i=Arbre[i].SeqSuiv) AffichAllume(i,Col);
}

int Aligner( t_pDescripteurGroupe DG,
		t_pDescripteurSequence DS,
		t_ClassHier Arbre,
		t_pGrdTab GrdTab,
		t_pParamAlign Param)
/* Aligne les Groupe[Groupe1] sequences commencant par Groupe1 en suivant Arbre
avec les Groupe[Groupe2] sequences commencant par Groupe2 en suivant Arbre
*/

{
	t_indL** Direction= NULL;
	t_Seq** Adresse= NULL;
	t_Weight** Profile= NULL;
	t_score TabTotCoeff[DERLETTREMAX+1];
	int OK = vrai;
	t_indL Saut;
	t_indL NbEnrg;
	t_indN i;
	t_Weight *WeightN= NULL/*, *RDGap=NULL*/;
	int Gros =(DG->Nb2 > Param->DerLettre)? Param->DerLettre:0;
	int Mult = (DG->Nb1>1)||(DG->Nb2>1);
/* Eclairage des sequences concernees */
	Eclairage(DG->G1,DG->Nb1,2,Arbre);
	Eclairage(DG->G2,DG->Nb2,3,Arbre);

	/* allocation de tous les tableaux et selection du swapping */
		NbEnrg = Initialiser( Gros,DG,GrdTab,&Adresse,&Profile,&WeightN,/*&RDGap,*/
			&Direction );
		if (NbEnrg <=0)
		{
				TraiteErr (8);
				OK= faux;
		}
		if( OK )
		{
			if (Mult) OK= Monter(DG->G2,DS->Sequence,Arbre,Param,GrdTab,DG,
				Adresse,Profile,WeightN,/*RDGap,*/TabTotCoeff,Direction,&Saut,&NbEnrg);
			else OK= Monter2(DG->G2,DS->Sequence,Arbre,Param,GrdTab,DG,
				Adresse,WeightN,TabTotCoeff,Direction,&Saut,&NbEnrg);
		}
	/* on libere les grands tableaux */
	if( GrdTab->Mu != NULL )
		free(GrdTab->Mu);
	if( GrdTab->X != NULL )
		free(GrdTab->X);
	if( Adresse != NULL )
		free(Adresse);
	if ( Profile != NULL )
	{
		for (i=0;i<=Param->DerLettre;i++) if (Profile[i]) free(Profile[i]);
		free(Profile);
	}
	if ( WeightN != NULL)
		free(WeightN);
/*	if (RDGap != NULL)
		free(RDGap);*/
/* 4/10/95 on utilise la fin de Direction */
	OK= OK && Descendre1( Direction + DG->M1 - NbEnrg,Saut, DG, GrdTab, NbEnrg );
	LibererDirection( Direction,DG->M1 );

	DetruitSwapD();

	OK= OK && Descendre2( DG, DS->Sequence, Arbre, GrdTab, Param->GapExtT, Saut );

	if( GrdTab->PMu != NULL )
		free(GrdTab->PMu);

	AffichHide ();
	Eclairage(DG->G1,DG->Nb1,1,Arbre);
	Eclairage(DG->G2,DG->Nb2,1,Arbre);

	MajTailleConsens(DS->Sequence,DS->NbSeq,faux);
	return( OK );
}

/*------------------------------------------------------------------------*/

static int ActAligner(t_indN Groupe1,t_indN Groupe2,t_indN NbSeq, t_indN *Groupe,
	t_ClassHier Arbre, void *V)

  /*	- Aligne Groupe1 et Groupe2 (df par NoGroupe et Groupe) si necessaire
	- Supprime Groupe2 dans Groupe  */
{
	struct _Aligner {
		t_pDescripteurGroupe DG;
		t_pDescripteurSequence DS;
		t_pParamAlign ParamAlign;} *Data = (struct _Aligner *)V;

		t_pDescripteurGroupe DG = Data->DG;
		t_pDescripteurSequence DS = Data->DS;
		t_pParamAlign ParamAlign = Data->ParamAlign;

	t_GrdTab GrdTab;
	int init;
	t_indL tmp;
	t_indN t;

	/*----------- initialisation des tableaux de GrdTab ------------*/
	GrdTab.Mu= NULL;
	GrdTab.PMu= NULL;
	GrdTab.X= NULL;

	/*----------- initialisation de DG -------------*/
	init = (DG->G1==-1);/*09/05/94 pour seq.partielles*/
	DG->G1= DG->NoGroupe[Groupe1];
	DG->G2= DG->NoGroupe[Groupe2];

	/* G1 doit etre le plus gros */
	if( DG->Groupe[DG->G1] < DG->Groupe[DG->G2] )
	{/* on permute */
		t = DG->G1; DG->G1= DG->G2; DG->G2= t;
		tmp = DG->M; DG->M= DG->N; DG->N= tmp;
		tmp = DG->M1; DG->M1= DG->N1; DG->N1= tmp;
	}
	/* Localisation des descriptions de groupes si pas encore fait */
	if (init)/*09/05/94 pour seq.partielles*/
	{
		ParamAlign->GapExtT = 0;
		if ((DG->M-DG->M1!=1) && (DG->N-DG->N1!=1)) ParamAlign->GapExtT++;
		if ((DG->M <DS->Sequence[DG->G1].Taille)
			&& (DG->N <DS->Sequence[DG->G2].Taille)) ParamAlign->GapExtT += 2;
	}
	else
	{
		DG->M= DS->Sequence[DG->G1].Taille;
		DG->M1= DG->M-1;

		DG->N= DS->Sequence[DG->G2].Taille;
		DG->N1= DG->N-1;
	}
	DG->Nb1= DG->Groupe[DG->G1];
	DG->Nb2= DG->Groupe[DG->G2];

	if( (DG->G1!=DG->G2) && !Aligner(DG, DS, Arbre,&GrdTab,ParamAlign))
		return faux;
	/* on regroupe */
	DG->Groupe[Groupe1]+= DG->Groupe[Groupe2];
	DG->Groupe[Groupe2]= 0;

	return vrai;
}

/*----------------------------------------------------------------------*/

static void CalculWeighti (t_indN N,t_indN i,t_indN D,t_ClassHier Arbre,
	t_pTabSeq Sequence)
{
		long tot;
		t_indN k, n;
		float s;
		for (n=0,k=i,tot=0;n<N;k=Arbre[k].SeqSuiv,n++) tot += Sequence[k].Weight;
		if (!tot) for (n=0,k=i;n<N;k=Arbre[k].SeqSuiv,n++) Sequence[k].Weight+= D;
		else
		{
			s = 1.0 + (1.0*D)/tot;
			for (n=0,k=i;n<N;k=Arbre[k].SeqSuiv,n++)
			Sequence[k].Weight = float_to_score(s*Sequence[k].Weight);
		}
}

static int Alignement ( t_ClassHier Arbre,
		t_DessinArb Dessin,
		t_pDescripteurGroupe DG,
		t_pDescripteurSequence DS,
		t_pParamAlign Param )

/* cree un alignement
Entree : Sequence donne les sequences alignees dans les groupes initiaux
		 Arbre est l'arbre a suivre pour aligner
		 NoGroupe donne les groupes initiaux
		 NbSeq est le nombre de sequences
		 Coeff donne le tableau des comparaisons de symbole
		 DerLettre donne le nombre de symboles
		 Gap + length * Gap2 donne la penalite pour un trou

Sortie : Sequence donne l'alignement selon Arbre
		 Alignement si on est alle jusqu'au bout
*/
{
	struct _Aligner {
		t_pDescripteurGroupe DG;
		t_pDescripteurSequence DS;
		t_pParamAlign ParamAlign;} Data;

	t_indN i,j,n, N=DS->NbSeq;
	t_indN *g, *noeud, Ni, Nj;
	t_pLigneSeq tmp;
	t_Dessin *des;
	t_score T2,T1 = (SCO_MAX / Param->CoeffMax) * 4 -1;
	Data.DG = DG;
	Data.DS = DS;
	Data.ParamAlign = Param;
	DG->G1 = 0; /* global alignment */
	for(i=0,g=DG->Groupe; i<=DS->NbSeq; i++,g++) *g=1;

/* sans les Weight on utilise TraiteBranche */
if (!Param->Weighted)
{
	if (N * N > T1){fprintf(stderr,"ss poids\n");TraiteErr(9); return faux;}
	for (n=0,tmp = DS->Sequence;n<= N;n++,tmp++) tmp->Weight = 1;
	return TraiteBranche(Arbre, DS->NbSeq,DG->Groupe,&Data,ActAligner);
}

/* avec les poids : boucle sur les noeuds*/
	g= DG->Groupe;
	if ((N*N*(Dessin[2].Col-Dessin[0].Col) > T1)&&(Dessin[2].Col>Dessin[0].Col))
	{
		T1 -= T2 = 1L * N * N;
		T2 *= (Dessin[2].Col-Dessin[0].Col);
		if (T1 > 0) for (n=0;n<=N;n++) Dessin[n].Col = Div(1.0*Dessin[n].Col,T1,T2);
		else for (n=0;n<=N;n++) Dessin[n].Col = 1;
/*      fprintf(stderr,"Dess/2\n");*/
/*	if (Dessin[2].Col<=Dessin[0].Col){fprintf(stderr,"Dess\n");TraiteErr(9); return faux;}*/
	}
	if ((noeud=(t_indN*)(calloc(N+1,sizeof(t_indN))))==NULL) return faux;
	for (n=0,i=Dessin[0].Col;n<=N;n++) noeud[n] = i;
	for (n=0,tmp = DS->Sequence;n<= N;n++,tmp++) tmp->Weight = 0;

	for (n=N,des=Dessin+n;n>1;n--,des--)
	{
		j = des->Noeud;
		i = Arbre[j].Pere;
		Ni = g[i];
		Nj = g[j];
		if (!ActAligner(i,j,N,g,Arbre,&Data)) { free (noeud); return faux;}

		CalculWeighti (Ni,i,des->Col-noeud[i],Arbre,DS->Sequence);
		CalculWeighti (Nj,j,des->Col-noeud[j],Arbre,DS->Sequence);
		noeud[i]=Dessin[n].Col;;
	}
	free(noeud);
	return vrai;
}	/* de Alignement */

/*-----------------------------------------------------------------------*/
static int ActRacine(	t_indN Groupe1, t_indN Groupe2, t_indN NbSeq, t_indN* Groupe,
			t_ClassHier Arbre, void *V)

/* On sait que sur Arbre, Groupe1, Groupe2 et Groupe3 se suivent avec
  Groupe1 pere de Groupe2 et Groupe2 non pere de Groupe3.
  Groupe2 est inclus dans Groupe1 si  c'est la meme chose sur AncArbre:
  Arbre       AncArbre
  G1 +-         G1 +-
  G2 + |        G2 + |
  G3 --+       AG3 --+	*/
{
	t_indN AncGroupe3;
	t_indN i;
	struct _RacineData {
	 t_ClassHier AncArbre;
	 int *Bon;
	 t_indN *NoGroupe;
	 } *Data = (struct _RacineData *)V;

	AncGroupe3= Data->AncArbre[Groupe2].SeqSuiv;
	Data->Bon[Groupe1]= Data->Bon[Groupe1] && Data->Bon[Groupe2]
			&& (Data->AncArbre[Groupe1].SeqSuiv == Groupe2)
			&& (Data->AncArbre[Groupe2].Pere == Groupe1)
			&& (Data->AncArbre[AncGroupe3].Pere != Groupe2);
	if( Groupe1 != 0 )
	{
		if( Data->Bon[Groupe1] )
		{
			for( i= 1; i<=NbSeq; i++ )
				if( Data->NoGroupe[i] == Groupe2 )
					Data->NoGroupe[i]= Groupe1;
			Groupe[Groupe1]+= Groupe[Groupe2];
			Groupe[Groupe2]= 0;
		}
		Arbre[Groupe1].SeqSuiv= Arbre[Groupe2].SeqSuiv;
		Data->AncArbre[Groupe1].SeqSuiv= Data->AncArbre[Groupe2].SeqSuiv;
	}
	return vrai;
} 	/* de ActRacine */

static int ArbEgal(t_ClassHier A1, t_ClassHier A2, t_indN NbSeq )
{
	t_indN j;

	for(j= 1; (j <= NbSeq)&&(A1[j].Pere==A2[j].Pere)
				&&(A1[j].SeqSuiv==A2[j].SeqSuiv);j++);
	return (j > NbSeq);
}

static int ArbCommun( t_ClassHier Arbre,
			t_ClassHier AncArbre,
			t_pDescripteurGroupe DG,
			t_indN NbSeq)
{
	struct
	{
		t_ClassHier AncArbre;
		int* Bon;
		t_indN* NoGroupe;
	} Data;
	t_ClassHier NvArbre; /* 7/7/94 local copy */
	t_indN i;

	/* initialisation et allocation des champs de Data */
	if(! (Data.Bon=(int*)calloc(NbSeq+1,sizeof(int))) )
	{
		TraiteErr(8);
		return faux;
	}
	for( i=0; i<=NbSeq; i++ )
		Data.Bon[i]= vrai;

	/*7/7/94 local copy */
	if(! (Data.AncArbre=(t_Branche*)calloc(NbSeq+1,sizeof(t_Branche))) )
	{
		TraiteErr(8);
		free (Data.Bon);
		return faux;
	}
	memmove (Data.AncArbre, AncArbre, (NbSeq+1)*sizeof(t_Branche));

	Data.NoGroupe= DG->NoGroupe;

	/* copie locale de Arbre */
	if(! (NvArbre=(t_Branche*)calloc(NbSeq+1,sizeof(t_Branche))) )
	{
		free (Data.AncArbre);
		free (Data.Bon);
		TraiteErr(8);
		return faux;
	}
	memmove (NvArbre, Arbre, (NbSeq+1)*sizeof(t_Branche));


	/* initialisation de NoGroupe et Groupe */
	for( i= 1; i<= NbSeq; i++ )
	{
		DG->NoGroupe[i]= i;
		DG->Groupe[i]= 1;
	}
	TraiteBranche( NvArbre, NbSeq, DG->Groupe,(void *)&Data, ActRacine);

	/* 7/7/94 free local copies */
	free(NvArbre);
	free(Data.AncArbre);
	free( Data.Bon );

	return vrai;
} 		/* de ArbCommun */


static int ActClusMark(t_indN Groupe1, t_indN Groupe2, t_indN NbSeq, t_indN* Groupe,
				t_ClassHier Arbre, void *V)
{
	struct _ClusMark {
		t_pLigneSeq Sequence;
		t_pParamAlign Param;
		float score;
		float norm;
	} *Data = (struct _ClusMark *)V;
	t_score sc, s;
	t_indN Seq1, Seq2, i, j;
	t_indN Nb1, Nb2;
	t_indL Pos;
	t_score GapLoc, g1, g2;
	int Ch, Ch2;
	int Ferme, Vide, Vide2;
	int Lettre;
	t_score TabCoeff[DERLETTREMAX+1];
	t_pLigneSeq tmp;


	Nb1= Groupe[Groupe1];
	Nb2= Groupe[Groupe2];

	if( Nb1 && Nb2 )
	{

	GapLoc= Nb1*Nb2;
	sc= 0;
	g1= 0;
	g2 =0;
	Ferme= (Data->Param->GapExtT >> 1) ? vrai : faux;

	for( Pos=1; Pos<=Data->Sequence[0].DerCol; Pos++ )
	{
		s= 0;
		Vide= vrai;
		Seq1= Groupe1;
		memset(TabCoeff, 0, (DERLETTREMAX+1)*sizeof(t_score));
		for( i=0; i< Nb1; i++ )
		{
			tmp= Data->Sequence+Seq1;
			Ch= tmp->Seq[Pos].Car;
			if( Ch && Pos<=tmp->DerCol )
			{
				if (Ch >= Data->Param->PremLettre) Vide= faux;
				for( Lettre=1; Lettre<=Data->Param->DerLettre;
							Lettre++ )
					TabCoeff[Lettre]+= Data->Param->Coeff[Ch][Lettre];
			}
			Seq1= Arbre[Seq1].SeqSuiv;
		}
		Vide2= vrai;
		Seq2= Groupe2;
		for( j=0; j<Nb2; j++ )
		{
			tmp= Data->Sequence+Seq2;
			Ch2= tmp->Seq[Pos].Car;
			if( Ch2 && Pos<=tmp->DerCol )
			{
				if (Ch2 >= Data->Param->PremLettre) Vide2= faux;
				s+= TabCoeff[Ch2];
			}
			Seq2= Arbre[Seq2].SeqSuiv;
		}
		if( !(Vide && Vide2) )
		{
			if (!Vide && g1 )
			{
				if (Ferme) sc -= (Data->Param->Gap + g1* Data->Param->Gap2) * GapLoc;
				Ferme= vrai;
				g1= 0;
			}
			else if (!Vide2 && g2)
			{
				if (Ferme) sc -= (Data->Param->Gap + g2* Data->Param->Gap2) * GapLoc;
				Ferme= vrai;
				g2= 0;
			}
			Ferme = Ferme || !(Vide || Vide2);
			if (Vide) g1++;
			else if (Vide2) g2++;
		}
	}
	if( Ferme && (Data->Param->GapExtT%1) && (g1+=g2))
		sc-= (Data->Param->Gap + g1*Data->Param->Gap2) * GapLoc;
	Data->score += (float)(sc) / GapLoc;


	Groupe[Groupe1]+= Nb2;
	Groupe[Groupe2]= 0;
	}
	return vrai;
}
static int CalculClusMark( t_pDescripteurSequence DS, t_pParamAlign Param,
				t_ClassHier Arbre, float *Score )
{
	t_indN i, *Groupe;
	struct _ClusMark {
		t_pLigneSeq Sequence;
		t_pParamAlign Param;
		float score;
		float norm;
	} Data;

	if( DS->NbSeq<2 )
	{
		*Score= 0;
		return vrai;
	}

	if( !(Groupe=(t_indN*)calloc(DS->NbSeq+1,sizeof(t_indN))) )
	{
		TraiteErr( 8);
		return faux;
	}

	for( i=0; i<=DS->NbSeq; i++ )
		Groupe[i]= 1;

	Data.Sequence = DS->Sequence;
	Data.Param = Param;
	Data.score= 0;
	Data.norm= 1.0/DS->NbSeq-1;

	TraiteBranche(Arbre, DS->NbSeq, Groupe, &Data, ActClusMark);

	*Score= Data.score * Data.norm;
	free (Groupe);
	return vrai;
}

static int CompareArbre (int Iter,
			int *OK,
			t_TabClassHier TabArbre,
			t_DessinArb *TabDessin,
			t_pDescripteurGroupe DG,
			t_pDescripteurSequence DS,
			t_pParamAlign Param )
/*	- regarde la nouveaute de Arbre[2]
	- s'il est nouveau, garde Arbre[1] dans AncArb et Arbre[2] dans Arbre[1]
	et decrit les groupes communs aux 2 arbres dans NoGroupe et Groupe*/
/*7/7/94 AncArb = Arbre[0]
 en Entree : OK si tout va bien
		Arbre[0] = Arbre[1] de l'etape precedente
		Arbre[1] arbre qui a servi a faire l'alignement
		Arbre[2] arbre deduit de l'alignement

 en Sortie :
		si Arbre[1] = Arbre[2], fin
		si Arbre[2] = AncArb, il y a un cycle.
		si pas de cycle
		 AncArb = Arbre[1]; Arbre[1] = Arbre[2];
		 NoGroupe et Groupe decrivent les groupes communs aux 2 arbres

 Retourne faux si les arbres sont identiques, vrai sinon.

 */
{
	register int i;
	int Egaux;
	float ClusMark[3];
	void *Ancien;

	if( ArbEgal(TabArbre[1], TabArbre[2], DS->NbSeq ))
	{
		if( Iter )
			CalculClusMark( DS, Param, TabArbre[1], ClusMark+1 );
		ClusMark[2]= ClusMark[1];
		return faux;
	}

	if( Iter ) /* au premier passage, AncArb==NULL */
	{
		Egaux=  (Iter>IterMax) || ArbEgal(TabArbre[0], TabArbre[2], DS->NbSeq);
		/*7/7/94 TabArbre[0]=AncArbre */
		if( !Egaux )
			*OK= !(CtrlBreakHit || Param->UNEITER);
		if( Egaux || ! *OK )
		{
			CalculClusMark( DS, Param, TabArbre[1], ClusMark+1 );
			CalculClusMark( DS, Param, TabArbre[2], ClusMark+2 );
			return faux;
		}
	}

	/* permutation des arbres  et recherche de la partie commune */
	if (!ArbCommun( TabArbre[1], TabArbre[2], DG, DS->NbSeq ))
	{
		OK = faux;
		return faux;
	}


	Ancien = TabArbre[0];
	for (i=0;i<2;i++) TabArbre[i]= TabArbre[i+1];
	TabArbre[2] = Ancien;
	memset (TabArbre[2],0,(DS->NbSeq+1)*sizeof(t_Branche));
	Ancien = TabDessin[0];
	for (i=0;i<2;i++) TabDessin[i]= TabDessin[i+1];
	TabDessin[2] = Ancien;
	memset (TabDessin[2],0,(DS->NbSeq+1)*sizeof(t_Dessin));

	return vrai;
}/*  de CompareArbre   */

/*-----------------------------------------------------------------------*/
static void CompleteArbre( t_ClassHier Arbre, t_indN NbSeq )
{
	t_indN PremPere, i;
	t_ClassHier A;

	PremPere= Arbre[0].SeqSuiv;
	if( PremPere == 0 )
		PremPere= 1;
	for(i=0, A=Arbre; i<PremPere; i++,A++)
		if( A->Pere == 0 )
			A->Pere= PremPere;
	for(i=PremPere+1,A=Arbre+i; i<=NbSeq; i++,A++ )
		if( A->Pere == 0)
			A->Pere= PremPere;
	Arbre[0].SeqSuiv= PremPere;
}



static void Sauvegarde(t_pDescripteurSequence DS)
{
	FILE *f;
	t_indN i;
	t_pLigneSeq Seq;
	int OK;

		f = fopen(savename,"wb");
		if (f == NULL) {unlink (savename); return;}
		OK = vrai;
		for (i = 0, Seq = DS->Sequence ; (i <= DS->NbSeq) && OK ; i++, Seq++) {
			if ((OK=fwrite(Seq,sizeof(t_LigneSeq),1,f)) &&
			 (Seq->BuffLen > 0))
				OK=fwrite(Seq->Seq+1,1,Seq->Taille,f);
		}
		if ((fclose(f) != 0) || !OK) unlink (savename);
}

static int Recupere(t_pDescripteurSequence DS,
			t_TabClassHier TabArbre,
			t_DessinArb *TabDessin)
{
	FILE *f;
	t_indN i;
	t_indL NewBufflen,OldBufflen;
	t_pLigneSeq Seq;
	t_pSeq NewSeq;
	int OK;
	void *tmp;

		MessageAction (96, NULL);
		OK = vrai;
		f = fopen(savename,"rb");
		if (f == NULL) return faux;
		for (i = 0, Seq= DS->Sequence ; (i <= DS->NbSeq) && OK ; i++, Seq++) {
			NewBufflen = Seq->BuffLen;
			NewSeq  = Seq->Seq;
			fread(Seq,sizeof(t_LigneSeq),1,f);

			if (Seq->BuffLen == 0) OldBufflen = 0;
			else OldBufflen = (Seq->Taille | 63)+1;

			Seq->BuffLen = NewBufflen;
			Seq->Seq = NewSeq;
			if (NewSeq && OldBufflen) {
				if ((OldBufflen==NewBufflen)
					||(Seq->Seq = realloc (NewSeq, OldBufflen)))
						Seq->BuffLen= OldBufflen;
				else Seq->Seq= NewSeq;
			}
			else if (OldBufflen) {
				if (Seq->Seq = malloc(OldBufflen)) Seq->BuffLen= OldBufflen;
				else Seq->BuffLen = OldBufflen = 0;
			}
			else {
				if (NewSeq) free (NewSeq);
				Seq->Seq = NULL;
				Seq->BuffLen = 0;
			}
			if (OldBufflen && ((Seq->Taille > Seq->BuffLen) ||
			 (!fread(Seq->Seq+1,1,Seq->Taille,f)))) OK = faux;
		}
		fclose(f);
		if (!OK) TraiteErr (16);
		tmp= TabArbre[2];
		for (i=1;i>=0;i--) TabArbre[i+1] = TabArbre[i];
		TabArbre[0] = tmp;
		tmp= TabDessin[2];
		for (i=1;i>=0;i--) TabDessin[i+1] = TabDessin[i];
		TabDessin[0] = tmp;

		return OK;
}



int AligneAction( t_TabClassHier TabArbre,
			t_DessinArb *Dessin,
			t_pDescripteurSequence DS,
			t_pParamAlign Param,
			t_pDescripteurGroupe DG )

/* calcule un alignement multiple avec iteration jusqu'a convergence ou cycle
Entree : Sequence contient les sequences alignees selon arbre[1]
	 Arbre[2] est l'arbre propose pour l'alignement
	 Consensus pointe vers le consensus
	 NbSeq est le nombre de sequences
	 Coeff donne le tableau des comparaisons de symbole
	 DerLettre donne le nombre de symboles
	 Gap donne la penalite pour un trou

Sortie :
	 Sequence donne l'alignement
	 Arbre[1] est l'arbre donnant naissance a l'alignement
	 Arbre[2] est l'arbre deduit de l'alignement		*/
{
	int OK, recup;
	t_indN i;
	int NbIter;
	t_pLigneSeq tmp;
	char M[20];
/*	FILE *debug=fopen("48.deb","wt");*/

	if( Param->DerLettre == 0 )
	{
		TraiteErr( 17);
		return faux;
	}
	TraducLetToByte(DS, Param);

/*
				for (i=1;i<=DS->NbSeq;i++)
				{
					fscanf(debug,"%4d %4d %4d %4d %4ld %*d %*d\n",
					&TabArbre[2][i].Pere,&TabArbre[2][i].SeqSuiv,&Dessin[2][i].Noeud,&Dessin[2][i].Col,
					&DS->Sequence[i].Weight);
				}  */
	CompleteArbre (TabArbre[2],DS->NbSeq);

	NbIter= 0;
	OK= vrai;
	recup= faux;

	AffichStart ("Aligning", DS->Sequence);
	AffichInit ("Position #", 5);
	tmp=DS->Sequence;
	if (tmp->Seq != NULL)
	{
		free (tmp->Seq);
		tmp->Seq = NULL;
		tmp->BuffLen = 0;
	};
	for (i=0,tmp=DS->Sequence+1;
		(i<DS->NbSeq)&&EditSizeLine(tmp,DS->Sequence->Taille); i++,tmp++);
	if (i<DS->NbSeq)
	{
		TraiteErr (8);
		OK = faux;
	}

	/* boucle des iterations */
	while( OK && CompareArbre(NbIter, &OK, TabArbre, Dessin, DG, DS, Param) )
	/*7/7/94 TabArbre[2]->TabArbre[1]->TabArbre[0] */
	{
		Sauvegarde (DS);
		NettoyerSequences( DS, DG->Groupe,TabArbre[1], Param);
/*				for (i=1;i<=DS->NbSeq;i++)
				{
					fprintf(debug,"%4d %4d %4d %4d %4ld %4d %hd\n",
					TabArbre[1][i].Pere,TabArbre[1][i].SeqSuiv,Dessin[1][i].Noeud,Dessin[1][i].Col,
					DS->Sequence[i].Weight,
					DS->Sequence[i].name,
					DS->Sequence[i].Seq[1].Symb);
				}
				fprintf(debug,"-----------------------\n");
*/
		if(OK= Alignement(TabArbre[1], Dessin[1], DG, DS, Param))
		{
			ArbreAlignement(Param,DS,TabArbre[2],Dessin+2,DS->NbSeq);
			NbIter++;
		}
		else if (NbIter>0) recup=Recupere (DS, TabArbre, Dessin);
	}

	TraducByteToLet(DS, Param);
	unlink (savename);
	CtrlBreakHit = faux;
	OK = OK || (Param->UNEITER && (NbIter>0)) || recup;
	sprintf (M," %d iteration(s)",NbIter);
	AffichEnd (OK,vrai,M);
/*	fclose(debug);*/
	return OK;
}		/* de AligneAction */

#ifndef clus2dom
int Align2Groups (int N2,t_ClassHier Arbre, t_pDescripteurSequence DS,
	t_pDescripteurGroupe DG, t_pParamAlign Param)
{
	struct _Aligner {
		t_pDescripteurGroupe DG;
		t_pDescripteurSequence DS;
		t_pParamAlign ParamAlign;} Data;
	int i;
	int N = DS->NbSeq;
	t_score wm,wn,T1 = SCO_MAX / Param->CoeffMax -1;
	float s;
	if ((N2-1) * (N-N2+1) > T1) {fprintf(stderr,"2gr debut\n");TraiteErr(9); return faux;}
	for (i=1;i<N2;i++) DG->NoGroupe[i]=1;
	for (;i<=N;i++) DG->NoGroupe[i]=N2;
	for (i=0;i<=N;i++) DG->Groupe[i]=0;
	DG->Groupe[1]=N2-1;
	DG->Groupe[N2]=N-N2+1;
	for (i=0;i<=N2;i++)
	{
		Arbre[i].Pere=1;
		Arbre[i].SeqSuiv = i+1;
	}
	for (;i<=N;i++)
	{
		Arbre[i].Pere=N2;
		Arbre[i].SeqSuiv=i+1;
	}
	Arbre[N].SeqSuiv=0;
	NettoyerSequences (DS,DG->Groupe,Arbre,Param);
	Data.DG = DG;
	Data.DS = DS;
	Data.ParamAlign = Param;
	DG->G1=1;

	for (i=1,wm=0;i<=N2;i++) wm += DS->Sequence[i].Weight;
	wm = max (N2-1, wm);
	for (wn=0;i<=N;i++) wn += DS->Sequence[i].Weight;
	wn = max (N-N2+1, wn);
	if (wm * wn > T1)
	{
		s = sqrt (1.0*T1/wm/wn);
		for (i=1;i<=N;i++) DS->Sequence[i].Weight =
			float_to_score(s*DS->Sequence[i].Weight);
/*			fprintf(stderr,"W/2\n");*/
	}
	return ActAligner(1,N2,N,DG->Groupe,Arbre,&Data);
}
#endif

static int Monter2( t_indN Seq21,t_pLigneSeq Sequence,t_ClassHier Arbre,
			t_pParamAlign Param,t_pGrdTab GrdTab,t_pDescripteurGroupe DG,
			t_Seq** Adresse, t_Weight* Weight,
			t_score*TabTotCoeff,t_indL** Direction,
			t_indL* Saut,
			t_indL* NbEnrg)
{
	register t_indL j;
	register t_indL i;
	size_t LineSize;
	t_score Nu, *mu;
	t_indL PNu, *pmu;
	t_score GapLocal, DGap, DGap1, DGap2, Temp, T1, GapExt, *T;
	t_indL ii, Pos;
	t_indL ** Direci, *Direcij;
	int Dinc;
	if (CtrlBreakHit) return faux;

	/* Show Message */
	AffichShow ((long)(DG->M1+1));

	Weight[0]=1;
	DGap = (t_score)1;
	GapLocal= DGap * Param->Gap;
	DGap *= Param->Gap2;
	DGap1 = (Param->GapExtT >> 1) ? DGap : 0;
	DGap2 = (Param->GapExtT & 1) ? DGap : 0;

/* initialisation de Mu a - l'infini */
	Temp= SCO_MIN+DGap;
	for( i=0, T=GrdTab->Mu; i<= DG->N1+1; i++,T++ ) *T= Temp;

/* ligne M : calcul du profile et du poids des sequences non vides*/
	CalculeCoeff(Sequence,Param,Arbre,TabTotCoeff,DG->M,1,DG->G1,1);

/* colonne N: memoriser les adresses des sequences
	poids ses sequences non vides
	cellule (M,N)
*/
	j = DG->N1+1;
	T=GrdTab->X+j;
	Dinc = j;
	Adresse[1]= Sequence[Seq21].Seq+DG->N;
	Ajoute(0,1,TabTotCoeff,Adresse,&Weight,1,T);
	for( j--,T--,
		GapExt = (Param->GapExtT & 1) ? GapLocal+DGap: 0;
		j>0; j--,T--,GapExt += DGap2)
		Ajoute(-GapExt,1,TabTotCoeff,Adresse,&Weight,1,T);

	GapExt = (Param->GapExtT & 1) ? GapLocal+DGap: 0;
	LineSize = (DG->N1+1)*sizeof(t_indL);
	for( i=DG->M1,Pos=DG->M-1,ii= i-1,
	Direci= Direction+ii;
	i>0; i--,ii--,Direci--,GapExt += DGap2)
	{
		Temp= GrdTab->X[DG->N1+1];
		Nu=SCO_MIN+DGap;/* 02/06/94 scores can be < 0 */
		PNu=0;
		CalculeCoeff(Sequence,Param,Arbre,TabTotCoeff,Pos--,1,DG->G1,1);
/* 4/10/95 on alloue ici les nouveaux tableaux */
		if(!*Direci && !(*Direci=(t_indL*)malloc(LineSize)))
		{
			*NbEnrg = DG->M1-1-ii;
			if (InitSwapD (*NbEnrg,DG->M1,LineSize)
			 ==SWAPIMPOSSIBLE) { TraiteErr(8); return faux; }
			AffichUpDate ((long)i+1);
			if (CtrlBreakHit) return faux;
			if (SwapD((void**)Direci+1,i)==SWAPIMPOSSIBLE)
				 { TraiteErr(8); return faux; }
			ii = *NbEnrg-1;
			Direci = Direction  + DG->M1-1;
		}

		Adresse[1]+=Dinc;/* tableau de pointeurs  */
		Ajoute( -GapExt,1,TabTotCoeff,Adresse,&Weight,1,GrdTab->X+DG->N1+1);
		/* pour chaque residu */

		for(j=DG->N1,T=GrdTab->X+j,mu=GrdTab->Mu+j,pmu=GrdTab->PMu+j,
			Direcij=*Direci+j ;
			j>0; j--,T--,Direcij-- )
		{
			T1= *T;
			CalculMu(DGap,Temp,mu,pmu);
			CalculMu(DGap,Temp,&Nu,&PNu);


			Ajoute (Sup (Temp, *mu-- - GapLocal, Nu - GapLocal, *pmu--, PNu, Direcij),
						1,TabTotCoeff,Adresse,&Weight,1, T);
			Temp= T1;
		}

		CalculMu (DGap1, Temp, GrdTab->Mu+0, GrdTab->PMu+0);

		if ( !(ii%100))
		{
			AffichUpDate ((long)i);
			if (CtrlBreakHit) return faux;
		}

		if( ( !ii ) && (i>1))	/* on swappe */
		{
/* 4/10/95 c'est la fin de Direction qui est remplie
 Direction devient Direci (c'etait deja vrai)
 et Direci devient Direction + DG->M1*/
			if (SwapD((void**) Direci, i-1 )==SWAPIMPOSSIBLE) {
				TraiteErr (8);
				return faux;
			}
			ii =  *NbEnrg;/* 02/06/94 updating of ii */
			Direci = Direction + DG->M1;
		}
	}

	for( j=DG->N1,T=GrdTab->X+j,Nu=*(T+1),PNu=0; j>0; j--,T--  )
		CalculMu(DGap1,*T,&Nu,&PNu);

	CalculMu(DGap1,GrdTab->X[1],GrdTab->Mu+0,GrdTab->PMu+0);

	GapExt = (Param->GapExtT >> 1) ? GapLocal : 0;
	Temp= Sup(GrdTab->X[1],GrdTab->Mu[0]-GapExt,Nu-GapExt,GrdTab->PMu[0],PNu,Saut);
	return vrai;
} 	/*de Monter2*/



