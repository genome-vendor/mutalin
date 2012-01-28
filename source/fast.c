/*------------------------------------------------------------------------
	fichier FAST.C

				Alignement rapide de toutes les sequences 2 a 2
					N*(N-1) comparaisons.

------------------------------------------------------------------------*/
/* Aligne les Sequences 2 a 2 par Fast
 Entree : Sequence contient les sequences (Taille, SeqS)
		NbSeq nombre de Sequences
		LongueurMax taille de la plus grande sequence
		[PremLettre..DerLettre] no des symboles reconnus
		Coeff table de comparaison de symboles
		Gap penalite pour un trou
		Score pointe vers un tableau de 0

 Sortie : Score pointe vers le tableau des scores

*/
#include "ma.h"
#include "msgerr.h"
#include "basproc1.h"
#include "afichseq.h"

#define MAXWINDOW 16
#define NB_ELEM_LIGNESCORE 2*MAXWINDOW
#define NB_DIAG 5
#define SCOMAX 800
#define ScoreMax(x) ((x)*SCOMAX + 500)/1000
/* types specifiques a FAST */
/*typedef t_indL* TabHash;
typedef t_indL* TabSuite;*/

typedef struct
{
	int Hmax;           	/* =2**(hshift*ktup)-1 <= 2**12 */
	int Hshift;             	/* 2**(hshift-1)<= DerLettre<2**hshift */
	int Ktup;		    	  	/* codage par ktuple */
	int Kt1;	   		  	/* ktup - 1 */
	t_indL* HashTab;
	t_indL* Suite;
	t_indL* DebutDiag;	    	/* debut de diagonale */
	t_indL* FinDiag;     	    	/* fin de diagonale */
	t_score* ScoreDiag;		/* score de diagonale */
	t_score** Score;
	t_score* ScoreId;
}t_DescripteurFast;
typedef t_DescripteurFast* t_pDescripteurFast;

typedef struct
{
	t_indL NoDiag;
	t_indL Debut;
	t_indL Fin;
	t_score Score;
} t_MaxElm; 	/* <-> TabMaxElm */
typedef t_MaxElm* t_pMaxElm;
typedef t_MaxElm t_TabMax[NB_DIAG+1];
typedef t_MaxElm* t_pTabMax;

typedef struct
{
	t_pSeq SeqT;   /* sequences */
	t_pSeq SeqC;
	t_indL LSeqT; 	/* longueur des sequences */
	t_indL LSeqC;
} t_FastSeq;
typedef t_FastSeq* t_pFastSeq;


typedef t_score t_LigneScore[NB_ELEM_LIGNESCORE];
typedef t_indL t_LigneIdent[NB_ELEM_LIGNESCORE];

/*************************   Prototypes  *********************************/
/*extern t_Weight Div (t_Weight A, t_Weight B, t_Weight C);*/
/* Si on a peur que Div deborde, utilser un float pour A */
#define Div(A,B,C) (((A)*(B) + ((C)>>1)) /(C))

static int FastAllocMat(t_pDescripteurFast DF,	t_pDescripteurSequence DS,
			t_pParamAlign Param );

static void FastFreeMat( t_pDescripteurFast DF, t_indN NbSeq );

static void HashTest( t_pDescripteurFast DF, t_pFastSeq FS, t_indL Hmask );

static void SauveDiag( t_pDescripteurFast DF, t_indL Diag, t_pTabMax DiagMax,
			t_score* ScoreMini, int* iScoreMini );

static t_score Recalcule( t_pMaxElm Dmax, t_pParamAlign Param , t_pFastSeq FS );

static t_indL MeilleureDiag( t_pDescripteurFast DF, t_pFastSeq FS, t_pParamAlign Param,
			/*t_indN nSeqC, t_pLigneSeq Sequence,*/ t_indL Hmask, t_indL nd0,
			t_indL Diag0, t_score Kfact, t_score Fact );

static t_score Cellule( t_indL i, t_indL j, t_pFastSeq FS, t_pParamAlign Param ,
	t_indL *Ident );

static t_score ScoreAlign( t_indL Decal, t_pFastSeq FS, t_pParamAlign Param ,
	t_indL *Ident );

/*int FastAction( t_pDescripteurSequence DS, t_pParamAlign Param,
			t_ClassHier TabArbre, t_DessinArb Dessin  );
*/
/*************************************************************************/


static int FastAllocMat(t_pDescripteurFast DF,
			t_pDescripteurSequence DS,
			t_pParamAlign Param )
/* allocation des tableaux et matrices pour le calcul des scores */
{
	t_indL NDiag, Nd2,j;
	int i0,c;


	/* calcul de hshift pour que  2**(hshift-1)<= DerLettre<2**hshift */

	for( c= 1,DF->Hshift= 0; c<=Param->DerLettre; c<<=1,DF->Hshift++);


	/* calcul de ktup pour que hmax = 2**(ktup*hshift)-1 < 2**12 */
	DF->Ktup= 12 / DF->Hshift;
	DF->Kt1= DF->Ktup-1;


	for( i0=0,j= 1; i0<=DF->Kt1; i0++,j<<=DF->Hshift);

	DF->Hmax= j-1;

	DF->HashTab= NULL;
	DF->Suite= NULL;
	DF->DebutDiag= NULL;
	DF->FinDiag= NULL;
	DF->ScoreDiag= NULL;
	DF->Score= NULL;

	/* determination des tailles de tableaux */
	NDiag= DS->Sequence[0].Taille;
	Nd2= 2*NDiag-1;

	/* allocation des tableaux */
	if ((DF->ScoreId=(t_score*)calloc(DS->NbSeq,sizeof(t_score)))==NULL)
		return faux;
	if( !(DF->Score=(t_score**)calloc(DS->NbSeq-1,sizeof(t_score*))) )
		return faux;
	if( !(DF->HashTab=(t_indL*)calloc(j,sizeof(t_indL))) )
		return faux;
	if( !(DF->Suite=(t_indL*)calloc(NDiag+1,sizeof(t_indL))) )
		return faux;
	if( !(DF->DebutDiag=(t_indL*)calloc(Nd2+1,sizeof(t_indL))) )
		return faux;
	if( !(DF->FinDiag=(t_indL*)calloc(Nd2+1,sizeof(t_indL))) )
		return faux;
	if( !(DF->ScoreDiag=(t_score*)calloc(Nd2+1,sizeof(t_score))) )
		return faux;

	return vrai;
} /* de AllocMat */


static void FastFreeMat( t_pDescripteurFast DF, t_indN NbSeq )
/* liberation de tous les tableaux et matrices de DF */
{
	t_indN i;
	t_score ** sz;

	if( DF->ScoreDiag != NULL )
		free(DF->ScoreDiag);
	if( DF->FinDiag != NULL )
		free(DF->FinDiag);
	if( DF->DebutDiag != NULL )
		free(DF->DebutDiag);
	if( DF->Suite != NULL )
		free(DF->Suite);
	if( DF->HashTab != NULL )
		free(DF->HashTab);
	if( DF->Score != NULL )
	{
		for (i=1,sz=DF->Score;(i<NbSeq)&&(*sz!=NULL);i++,sz++) free (*sz);
		free(DF->Score);
	}
	if (DF->ScoreId != NULL) free(DF->ScoreId);
}



static void HashTest( t_pDescripteurFast DF, t_pFastSeq FS, t_indL Hmask )
{
	t_indL i;
	t_indL *h,j;
	t_Seq *s;

	for( i=0,h=DF->HashTab; i<=DF->Hmax; i++,h++ )
		*h= -1;

	memset( DF->Suite, 0,(FS->LSeqT+1)*sizeof(t_indL) );

	for( i=1, j=0, s=FS->SeqT+1; i<=DF->Kt1; i++,s++ )
		j= (j << DF->Hshift) + s->Symb;

	for( i=DF->Ktup,s=FS->SeqT+i,h=DF->Suite+i; i<=FS->LSeqT; i++,s++,h++ )
	{
		j=((j & Hmask) << DF->Hshift) + s->Symb;
		*h= DF->HashTab[j];
		DF->HashTab[j]= i;
	}
} /* de HashTest */



static void SauveDiag( t_pDescripteurFast DF, t_indL Diag,	t_pTabMax DiagMax,
			t_score* ScoreMini, int* iScoreMini )
{
	t_MaxElm* tmp;
	int i,ii;

	for (i= 1,tmp=DiagMax+1;(i<=NB_DIAG) && ((tmp->NoDiag != Diag)
				|| (tmp->Debut != DF->DebutDiag[Diag]));
			i++,tmp++);

	if (i <= NB_DIAG)
	{
		tmp->Fin= DF->FinDiag[Diag];
		if( DF->ScoreDiag[Diag] > tmp->Score )
			tmp->Score= DF->ScoreDiag[Diag];
		if( i == *iScoreMini )
		{
			*ScoreMini= LONG_MAX;
			for( ii=1,tmp=DiagMax+1; ii<=NB_DIAG; ii++, tmp++ )
				if( tmp->Score < *ScoreMini )
				{
					*iScoreMini= ii;
					*ScoreMini= tmp->Score;
				}
		}
	}
	else
	{
		tmp= DiagMax+(*iScoreMini);
		tmp->Score= DF->ScoreDiag[Diag];
		tmp->NoDiag= Diag;
		tmp->Debut= DF->DebutDiag[Diag];
		tmp->Fin= DF->FinDiag[Diag];
		*ScoreMini= tmp->Score;
		for( i=1,tmp=DiagMax + 1; i<=NB_DIAG; i++, tmp++ )
			if( tmp->Score < *ScoreMini )
			{
				*iScoreMini= i;
				*ScoreMini= tmp->Score;
			}
	}
} 	/* de SauveDiag */





static t_score Recalcule( t_pMaxElm Dmax, t_pParamAlign Param , t_pFastSeq FS )
{
	t_indL Pos, PosT;
	t_score Tot;
	t_Seq *sC, *sT;

	for (Tot= 0,Pos=Dmax->Debut, sC=FS->SeqC+Pos,
		PosT= Dmax->Debut + FS->LSeqT - Dmax->NoDiag,sT=FS->SeqT+PosT;
		Pos<=Dmax->Fin; Pos++,sC++,PosT++,sT++ )

		Tot+= Param->Coeff[sT->Symb][sC->Symb];

	return Tot;
} /* de Recalcule */


static t_indL MeilleureDiag( t_pDescripteurFast DF,
			t_pFastSeq FS,
			t_pParamAlign Param,
/*			t_indN nSeqC ,
			t_pLigneSeq Sequence,*/
			t_indL Hmask,
			t_indL nd0,
			t_indL Diag0,
			t_score Kfact,
			t_score Fact	)
{
	t_indL Pos, PosT, Diag, Diag1;
	t_indL Val,Fin, Dist, dd;
	int iScoreMini,   im, ib;
	t_indL nd;
	t_TabMax DiagMax;
	t_MaxElm* tmp;
	t_Seq *sC;
	t_score* sco, ScoreMini,Score, NvScore,ScoreMax,Ajoute;

/*	FS->LSeqC= Sequence[nSeqC].DerCol;
	FS->SeqC= Sequence[nSeqC].Seq;*/

	nd= nd0 + FS->LSeqC;


	memset(DF->DebutDiag, 0, (nd+1)*sizeof(t_indL));
	memset(DF->FinDiag, 0, (nd+1)*sizeof(t_indL));
	memset(DF->ScoreDiag, 0, (nd+1)*sizeof(t_score));
	memset(DiagMax, 0, sizeof(t_TabMax));
	iScoreMini= 1;
	ScoreMini= 0;

	Val= 0;
	for( Pos=1,sC=FS->SeqC+1; Pos<=DF->Kt1; Pos++,sC++ )
		Val= ((Val & Hmask)<<DF->Hshift) + sC->Symb;

	Diag1= Diag0;
	for( Pos=DF->Ktup,sC=FS->SeqC+Pos; Pos<=FS->LSeqC; Pos++,sC++ )
	{
		Diag1++;
		Val= ((Val & Hmask)<<DF->Hshift) + sC->Symb;
		PosT= DF->HashTab[Val];
		while( PosT > 0)
		{
			Diag= Diag1 - PosT;
			Fin= DF->FinDiag[Diag];
			if( Fin == 0 )		/* nouvelle diagonale */
			{
				DF->DebutDiag[Diag]= Pos - DF->Kt1;
				DF->FinDiag[Diag]= Pos;
				DF->ScoreDiag[Diag]= Kfact;
			}
			else 	/* poursuite d'une diagonale */
			{
				Dist= Pos - DF->Ktup - Fin;
				if( Dist < 0 )      /* ktuples chevauchants */
				{
					Ajoute= Fact*(Pos-Fin);
					Dist= 0;
				}
				else
					Ajoute= Kfact;
				Score= DF->ScoreDiag[Diag];
				NvScore= Score + Ajoute;
				if( NvScore > Dist )
					NvScore-= Dist;
				else
					NvScore= 0;
				/* penalite proportionnelle a la longueur du dismatch */
				if( (NvScore<Score) && (Score>ScoreMini) )
					SauveDiag(DF,Diag,DiagMax,&ScoreMini,&iScoreMini);

				if( NvScore < Kfact )	/* nouveau debut de diagonale */
				{
					DF->DebutDiag[Diag]= Pos - DF->Kt1;
					DF->FinDiag[Diag]= Pos;
					DF->ScoreDiag[Diag]= Kfact;
				}
				else /* vraie poursuite de diagonale */
				{
					DF->FinDiag[Diag]= Pos;
					DF->ScoreDiag[Diag]= NvScore;
				}
			} /* end else */
			PosT= DF->Suite[PosT];
		} /* end PosT */
	} /* end Pos */

	for( Diag=1,sco=DF->ScoreDiag+1; Diag<=nd; Diag++,sco++ )
		if( *sco > ScoreMini )
			SauveDiag(DF,Diag,DiagMax,&ScoreMini,&iScoreMini);

	ib= 0;
	ScoreMax= SCO_MIN;
	dd= 0;
	for( im=1,tmp=DiagMax+im; im<=NB_DIAG; im++,tmp++ )
	{
		Score= Recalcule(tmp, Param, FS );
		if( (Score > ScoreMax)
		|| ((Score==ScoreMax) && (abs(tmp->NoDiag - FS->LSeqT) < dd)) )
		{
			ib= im;
			ScoreMax= Score;
			dd= abs(tmp->NoDiag - FS->LSeqT);
		}
	}
	if( ib > 0 )
		return DiagMax[ib].NoDiag;
	else
		return 0;
}	 /* de MeilleureDiag	*/





static t_score Cellule( t_indL i, t_indL j, t_pFastSeq FS, t_pParamAlign Param,
	t_indL *Ident )
{
	t_score Res= 0;
   *Ident=0;

	if( (i >= 1) && (i <= FS->LSeqT) && (j >= 1) && (j <= FS->LSeqC) )
	{
		Res= Param->Coeff[FS->SeqT[i].Symb][FS->SeqC[j].Symb];
		*Ident = (FS->SeqT[i].Symb == FS->SeqC[j].Symb)? 1:0;
	}
	return Res;
}

static t_score MaxWay2 (t_score x1,t_score x2,t_indL p1,t_indL p2,t_indL *p)
{
	if (x1 >= x2)
	{
		if (x1 == x2) *p = max(p1,p2);
		else *p = p1;
		return x1;
	}
	*p = p2;
	return x2;
}

static t_score MaxWay3(t_score x, t_score m, t_score n, t_indL im, t_indL in,
	t_indL *id)
{
	if (x >= m)
	{
		if ((x==n)&&(*id<in)) *id = in;
		if (x >= n) return x;
	}
	else if (m>=n)
	{
		if (m == n)	*id = max (im,in);
		else *id = im;
		return m;
	}
	*id = in;
	return n;
}

static t_score ScoreAlign( t_indL Decal, t_pFastSeq FS, t_pParamAlign Param,
	t_indL *Ident )
{
	t_LigneScore Mu, X;
	t_LigneIdent Id,Idm;
	t_score *x, *m, Nu, Res;
	t_indL i,j,jj,*id,*idm,idn,it;
	int k, ik, Window;
	t_indL DebutT, DebutC, Fin;
	t_score Gap = Param->Gap;
	t_score Gap2 = Param->Gap2;
	t_score mini = SCO_MIN+Gap+Gap2;
	t_indL LC = FS->LSeqC;

	if( Decal >= 0 )
	{
		DebutT= 1;
		DebutC= Decal + 1;
	}
	else
	{
		DebutT= 1 - Decal;
		DebutC= 1;
	}

	Window= min(MAXWINDOW-1,(LC-DebutC)/2);

	for (k=-Window,m=Mu;k<=Window+1;k++,m++) *m = mini;

	Fin= min (LC + Window - Decal,FS->LSeqT);

	if ((i=DebutT-Window-1) <= 0)/* first row */
		for (k=-Window,i=1,j=1+Decal,ik=j+k,x=X,m=Mu,id=Id,idm=Idm;k<=Window;
			k++,ik++,x++,m++,id++,idm++)
		{
			*x= Cellule (i,ik,FS,Param,id);
			*m = *x;
			*idm = *id;
		}

	for (i++,j=i+Decal,jj=1-j+Window;(i<=Fin)&&(j-Window<2);i++,j++,jj--)
	{
		x=X+jj;/* first column */
		m=Mu+jj;
		id=Id+jj;
		idm=Idm+jj;
		*m = *x = Cellule (i,1,FS,Param,id);
		Nu = mini;
		idn = 0;
		for (k=2-j,ik=2,x++,m++,id++,idm++;(k<=Window);k++,ik++,x++,m++,id++,idm++)
		{
			Nu = MaxWay2(Nu-Gap2, *x,idn,*id,&idn);
			*x= Cellule(i,ik,FS,Param,&it) + MaxWay3(*x,*m-Gap,Nu-Gap,*idm,idn,id);
			*id += it;
			*m= MaxWay2(*(m+1)-Gap2,*x,*(idm+1),*id,idm);
		}
	}


	for(Res=X[2*Window],*Ident=Id[2*Window]; i<=Fin; i++,j++ )/* Window left side complete */
	{
		for( k=-Window,ik=j+k,x=X,m=Mu,Nu = mini,id=Id,idm=Idm,idn=0;
		 (k<=Window)&&(ik<=LC); k++,ik++,x++,m++,id++,idm++ )
		{
			Nu = MaxWay2 (Nu - Gap2, *x,idn,*id,&idn);
			*x= Cellule(i,ik,FS,Param,&it) + MaxWay3(*x,*m-Gap,Nu-Gap,*idm,idn,id);
			*id +=it;
			*m= MaxWay2( *x, *(m+1) - Gap2,*id,*(idm+1),idm);
		}
		if(k<=Window) Res = MaxWay2 (Res,*(x-1),*Ident,*(id-1),Ident);
			/* Window right side incomplete */
		else { Res = *(x-1); *Ident = *(id-1);}
	}

	for (k=max(-Window,2-j),ik=j-1+k,x=X+Window+k,id =Id+Window+k;
		(k<Window)&&(ik<LC);k++,ik++,x++,id++)
			Res = MaxWay2 (Res,*x,*Ident,*id,Ident); /* last row */

	return Res;
} /* de ScoreAlign */


static t_score ScoreIdent(t_pParamAlign Param,t_pLigneSeq S)
{
	t_indL i;
	t_pSeq Seq;
	t_score sco;
	t_FormatScores test=Param->ScMethod;
	if(test>PERCENT) return 1000;
	for (i=0,Seq = S->Seq+1, sco=0;i<S->DerCol;i++,Seq++)
		sco+= Param->Coeff[Seq->Symb][Seq->Symb];
	if (test==1) sco = sco*1000/S->DerCol;
	return sco;
}

static void MeanScore( t_score* Score1, t_score Score2,
			t_indN Poids1, t_indN Poids2, t_indN PoidsT )
{
	*Score1 = (Poids1 * (*Score1) + Poids2 * Score2 + (PoidsT>>1))/ PoidsT;
}

int FastAction2( t_pDescripteurSequence DS,
			t_pParamAlign Param,
			t_ClassHier Arbre,
			t_DessinArb Dessin )
{
	t_DescripteurFast DF;
	t_FastSeq FS;
	t_indL Hmask;
	t_indN Seq1, Seq2, i,j,k,l,m, Noeud1, Noeud2,G1,G2,P1,P2,PT;
	t_score Kfact, Fact;    /* scoring factors */
	t_indL nd0, Diag0, Ident;
	t_score *sc, **sz, sco, scomax, Sup, InfSup, SupSup;
	int Res= faux;
	t_pLigneSeq s1;
	t_indN *Groupe, *NoSeq;
/*	FILE *f; * pour sauver la matrice de score*/

	if( !FastAllocMat( &DF, DS, Param ) )
	{
		TraiteErr( 8);
		FastFreeMat (&DF, DS->NbSeq);
		return faux;
	}


	TraducLetToByte( DS, Param );

	Hmask=  ((DF.Hmax+1) >> DF.Hshift)-1;

	Fact= 4*DF.Ktup;
	Kfact= DF.Ktup*Fact;

	AffichStart ("Clustering with fast",DS->Sequence);
	for (Seq1=1,s1=DS->Sequence+1,sc=DF.ScoreId,scomax=0;Seq1 <= DS->NbSeq;
	Seq1++,s1++,sc++)
		if ((*sc = ScoreIdent(Param,s1))>scomax) scomax=*sc;

	/*f=fopen("cytc.sco","wt");*/
	if (((Groupe=calloc(DS->NbSeq+1,sizeof(t_indN)))==NULL)||
		((NoSeq=calloc(DS->NbSeq,sizeof(t_indN)))==NULL))
	{
		TraiteErr( 8);
		if (Groupe) free(Groupe);
		FastFreeMat (&DF, DS->NbSeq);
		return faux;
	}

	for (i=1;i<=DS->NbSeq;i++) Groupe[i]=1;
	memset(Arbre, 0, (DS->NbSeq+1)*sizeof(t_Branche));
	NoSeq[0]=1;
	m=k=DS->NbSeq;
	SupSup = ScoreMax(scomax);
	for( i=0, sz=DF.Score; (k>1) && !CtrlBreakHit; i++, sz++ )
	{
		Seq1=NoSeq[i];
		AffichAllume (Seq1,2);
		s1 = DS->Sequence+Seq1;
/*	fprintf (f,"%8s",s1->NomSeq);*/
		FS.LSeqT= s1->DerCol;
		FS.SeqT= s1->Seq;
		nd0= FS.LSeqT-1;
		Diag0= FS.LSeqT + DF.Kt1;
		HashTest( &DF, &FS, Hmask );

		if (!(*sz=(t_score*)calloc(k-1,sizeof(t_score)))) break;
		k=0;
		l=Seq1;
		for( Seq2=Seq1+1; (Seq2 <= DS->NbSeq) && !CtrlBreakHit; Seq2++)
		if (Groupe[Seq2]>0)
		{
			AffichAllume (Seq2,3);
			FS.LSeqC= DS->Sequence[Seq2].DerCol;
			FS.SeqC= DS->Sequence[Seq2].Seq;
			sco= ScoreAlign( MeilleureDiag(&DF,&FS,Param,/*Seq2,
			DS->Sequence,*/Hmask,nd0,Diag0,Kfact,Fact)- FS.LSeqT,
			&FS, Param, &Ident);
			switch (Param->ScMethod)
			{
				case ABSOLU: break;
				case PERCENT: sco = 1000.0 * sco / min(FS.LSeqC, FS.LSeqT); break;
				case IDENTITY: sco = 1000.0 * Ident / min(FS.LSeqC, FS.LSeqT); break;
				case NORMALISED: sco = 2000.0 * sco /(DF.ScoreId[Seq1-1]+DF.ScoreId[Seq2-1]);
				 break;
			}
/*		fprintf (f,"%5ld",*sc);*/
			if (sco>=SupSup)
			{
				Groupe[Seq1]++;
				Groupe[Seq2]=0;
				Dessin[m--].Noeud=l=Arbre[l].SeqSuiv=Seq2;
				Arbre[Seq2].Pere=Seq1;
			}
			else
			{
				if (!k) NoSeq[i+1]=Seq2;
				(*sz)[k++]=sco;
			}
			AffichAllume (Seq2,1);
		}
		AffichAllume (Seq1,1);
/*	fprintf (f,"\n");*/
	}

/*fclose(f);*/
	Res = !CtrlBreakHit && (k<2);
	AffichEnd (Res,vrai,"");
	CtrlBreakHit = faux;

	if (Res)
	{
		l=m-1;
		for (i=0;i<l;i++) for (j=0,k=0,Seq2=NoSeq[i+1];j<l-i;Seq2++)
		{
			if (Groupe[Seq2]>0) DF.Score[i][j++]=DF.Score[i][k++];
			else if (Arbre[Seq2].Pere>NoSeq[i]) k++;
		}

		while (l>0)
		{
			Sup=DF.Score[l-1][0];
			Noeud1=l-1;
			Noeud2=l;
			for (i=0;i<l-1;i++) for (j=0;j<l-i;j++) if (DF.Score[i][j]>Sup)
			{
				Noeud1=i;
				Noeud2=i+j+1;
				Sup=DF.Score[i][j];
			}
			if (l==DS->NbSeq-1) SupSup=Sup;
			InfSup=Sup;
			G1 = NoSeq[Noeud1];
			G2 = NoSeq[Noeud2];
			P1 = Groupe[G1];
			P2 = Groupe[G2];
			PT = P1+P2;
			/* (i,Noeud1)+(i,Noeud2)->(i,Noeud1) & (i,l)->(i,Noeud2) */
			for (i=0;i<Noeud1;i++)
			{
				MeanScore(DF.Score[i]+Noeud1-i-1,DF.Score[i][Noeud2-1-i],P1,P2,PT);
				DF.Score[i][Noeud2-1-i]=DF.Score[i][l-i-1];
			}
			for (i++;i<Noeud2;i++)
			{
				MeanScore(DF.Score[Noeud1]+abs(i-Noeud1)-1,
					DF.Score[i][Noeud2-i-1],P1,P2,PT);
				DF.Score[i][Noeud2-i-1]=DF.Score[i][l-i-1];
			}
			for (i++;i<l;i++)
			{
				MeanScore(DF.Score[Noeud1]+i-Noeud1-1,DF.Score[Noeud2][i-Noeud2-1],
					P1,P2,PT);
				DF.Score[Noeud2][i-1-Noeud2]=DF.Score[i][l-i-1];
			}
			if ((l!=Noeud1)&&(l!=Noeud2))
			{
				/* (Noeud1,l)+(Noeud2,l)->(Noeud1,l)->(Noeud1,Noeud2) */
				MeanScore(DF.Score[Noeud1]+l-Noeud1-1,DF.Score[Noeud2][l-Noeud2-1],
					P1,P2,PT);
				DF.Score[Noeud1][Noeud2-Noeud1-1]=DF.Score[Noeud1][l-Noeud1-1];
			}

			DF.Score[l-1][0]=Sup;
			NoSeq[Noeud2]=NoSeq[l];
			if (G1>G2) {i=G1;G1=G2;G2=i;}
			NoSeq[Noeud1]=G1;
			Dessin[l+1].Noeud=G2;
			SupprimeGroupe(Arbre,Groupe,G1,G2);
			l--;
		}
	Arbre[0].Pere= 1;
	Arbre[0].SeqSuiv= 1;
	InfSup = SupSup - InfSup;
	if (!InfSup) InfSup = 1;
	if (TailleGraphe==0) TailleGraphe=InfSup=1;
	Dessin[1].Col=DebutGraphe+scomax;
	for (l=DS->NbSeq;l>m;l--) Dessin[l].Col = DebutGraphe;
	for (sz=DF.Score+l-2; l>1; l--,sz--)
		Dessin[l].Col = DebutGraphe + TailleGraphe * (SupSup - **sz)/ InfSup;
	Dessin->Col=DebutGraphe-1;
	}

	TraducByteToLet( DS, Param );

	FastFreeMat( &DF,DS->NbSeq );
	free(Groupe);
	free(NoSeq);
	return Res;
}

