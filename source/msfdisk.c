/*------------------------------------------------------------------------
	Fichier MSFDISK.C
		Sauvegarde des fichiers au format MSF (GCG)
------------------------------------------------------------------------*/

#include <string.h>
#include <ctype.h>
#include "commande.h"
#include "disk.h"
#include "ma.h"

static void Ecrit1LigneSeqDoc(FILE* outf,t_pSeq Seq,t_pSeq Cons,t_pSeq First_Seq,int IsCons,int IsFirst_Seq ,
	t_indL Start,t_indL Finish,t_pParamAlign Param,int BlockSize,t_Style OutputStyle);

static const char MSFformat[56] =
	" Name: %-14s Len: %5ld  Check: %4d  Weight:  %3.2f\n";
#define MSFnameLength 14

char ConsSymb[DERLETTREMAX+1][DERLETTREMAX+1];

int CalculeSymboles (t_pParamAlign Param)
{
	char *s;
	int i,h,m;
	memset (ConsSymb,0,sizeof(ConsSymb));
	m=0;
	for (i=1;i<= DERLETTREMAX;i++)if ((h=Param->Hom[i])!=i)
	{
		s = ConsSymb[h];
		m++;
		if (!s[0]) s[0]=Param->Symb[h];
		s[strlen(s)]=Param->Symb[i];
	}
	return m;
}

static void EcritSymboles (FILE* OUTFILE,t_pParamAlign Param)
{
	int i;
	char *s;
	if (!CalculeSymboles(Param)) return;
	fprintf(OUTFILE,"Consensus symbols:\n");
	for (i=1;i<= DERLETTREMAX;i++)if ((s=ConsSymb[i])[0])
		fprintf(OUTFILE," %c is anyone of %s\n",s[0],s+1);
}

static void Ecrit1LigneSeq( FILE* OUTFILE,t_pSeq Seq,int Start,int Finish,int BlockSize)
/* ecrit BlockSize residus de la sequence p_Seq */
{
	int c;
	t_indL Pos;
	int  Count=0;
	t_Seq *S;

	for( Pos=Start, S=Seq+Pos; Pos<=Finish; Pos++,S++ )
	{
		c= S->Car;
		if( c==' ' || c==Insert )
			c= '.';

		fputc( c, OUTFILE );
		Count++;
		if( Count==BlockSize && Pos!=Finish)
		{
			fputc( ' ', OUTFILE );
			Count= 0;
		}
	}
	fputc( '\n', OUTFILE );
}

static int GCGchar (int c)
{
	if ((c==Insert)||(c==' ')) return '.';
	else return toupper(c);
}

static int GCGCheckSum (t_pLigneSeq Seq)
{
	t_indL i;
	int Count;
	long Check;
	t_pSeq s;

	for (i=1, s= Seq->Seq+1, Count=1, Check=0; i<=Seq->Taille; i++,s++,Count++)
	{
		Check += Count * GCGchar(s->Car);
		if (Count == 57) Count=0;
	}

	return (int)(Check % 10000);
}

static void MsfEcritEnTete( FILE* OUTFILE,t_pDescripteurSequence DS,t_pParamAlign Param,
	long MaxLen)
{
	t_indN i, *n;
	t_pLigneSeq tmp;
	t_StringNom name;
	float Weight0 = ((float)DS->NbSeq) / DS->Sequence[0].Weight;
	PrintCopyRight (OUTFILE);
	fprintf(OUTFILE,"Symbol comparison table: %s\n", Param->NomCoeff);
	fprintf(OUTFILE,"Gap weight: %g\n", (float)Param->Gap);
	fprintf(OUTFILE,"Gap length weight: %g\n", (float)Param->Gap2);
	fprintf(OUTFILE,"Consensus levels: high=%d%% low=%d%%\n",(int)Param->ConsHlevel,
		(int)Param->ConsLlevel);
	EcritSymboles(OUTFILE,Param);
	fprintf(OUTFILE,"\n MSF: %6ld    Check:   0         ..\n",MaxLen);
	for( i=1,n= DS->NoSeq; i<= DS->NbSeq; i++)
	{
		tmp=DS->Sequence+ (n ? *(++n) : i);
		strcpy(name,tmp->NomSeq);
		name[MSFnameLength]=0;
		fprintf( OUTFILE, MSFformat,
			name, MaxLen, GCGCheckSum(tmp), tmp->Weight*Weight0 );

	}
	fprintf (OUTFILE, MSFformat,
		DS->Sequence[0].NomSeq, MaxLen, GCGCheckSum(DS->Sequence),0.0);
	fprintf (OUTFILE, "\n//\n\n");
}

void MsfEcritSeq( FILE* OUTFILE,t_pDescripteurSequence DS,t_pParamAlign Param,
	int HasStyle, int Group2 )
/* sauvegarde au format MSF les sequences de DS dans OUTFILE */
{
	int LineSize= 50;		/* nombre de residus par ligne */
	int BlockSize= 10;     /* decoupage des residus par bloc de 10 */
	t_Style style= Normal;
	t_indL Start, Finish, Limit, Spaces,j;
	long MaxLen;
	t_indN i, *n, ind;
	t_pLigneSeq tmp,FirstSeq;
	int nomf=20;
	t_StringNom name;
	if (HasStyle)
	{
		nomf=10;
		PEGetOutputStyleParam(&style,&LineSize,&BlockSize);
	}
	MaxLen= DS->Sequence[0].DerCol;

	/* ecrire l'entete de fichier */
	MsfEcritEnTete (OUTFILE,DS,Param,MaxLen);

	/* ecrire les sequences */
	Start= 1;
	Limit= 10;
	do
	{
		while( Start>=Limit)
			Limit+= 10;

		Finish= (Start+LineSize)-1;
		if( Finish > MaxLen )
			Finish= MaxLen;
		Spaces= (Finish-Start) / BlockSize;

		/* ecrire les numeros de debut et de fin de ligne */
		fprintf(OUTFILE,"%*c%-5ld",nomf+2,' ',(long)Start);
		if( (Finish-Start+1) >= 10 )
		{
			for( j=1; j<=(Finish-Start+1+Spaces-10); j++ )
				fputc(' ',OUTFILE);
			fprintf(OUTFILE,"%5ld\n", (long)Finish);
		}
		else
			fputc('\n', OUTFILE);

		/* ecrire les NbSeq lignes de sequences */
		for( i=1, n= DS->NoSeq; i<=DS->NbSeq; i++ )
		{
			ind =  (n ? *(++n) : i);
			tmp= DS->Sequence+ ind;
			if (i==1) FirstSeq=tmp;
			if( Start <= tmp->Taille )
			{
				strcpy(name,tmp->NomSeq);
            name[nomf]=0;
				if (!HasStyle)
				{
					fprintf( OUTFILE,"%*s  ",nomf,name);
					Ecrit1LigneSeq( OUTFILE, tmp->Seq, Start, Finish, BlockSize);
				}
				else
				{
					if (ind < Group2)
						fprintf( OUTFILE,"[%*s]  ",nomf,name);
					else
						fprintf( OUTFILE,"%*s  ",nomf,name);
					Ecrit1LigneSeqDoc(OUTFILE,tmp->Seq,DS->Sequence[0].Seq,
						FirstSeq->Seq,0,i==1,Start,Finish,Param,BlockSize,style);
				}
			}
		}
		/* ecrire le consensus */
		fprintf( OUTFILE,"%*s  ",nomf,DS->Sequence[0].NomSeq);
		if (!HasStyle)
		Ecrit1LigneSeq( OUTFILE, DS->Sequence[0].Seq, Start, Finish, BlockSize);
		else
		Ecrit1LigneSeqDoc(OUTFILE,DS->Sequence[0].Seq,DS->Sequence[0].Seq,
					FirstSeq->Seq,1,i==1,Start,Finish,Param,BlockSize,style);

		fputc( '\n', OUTFILE );
		Start+= LineSize;
	} while( Finish < MaxLen);

}

/*----------------------------------------------------------
	Sauvegarde d'une sequence alignee avec un style
------------------------------------------------------------*/

static const char Ouvre[]=" ([", Ferme[]=" )]";
static void Ecrit1LigneSeqDoc(FILE* outf,t_pSeq Seq,t_pSeq Cons,t_pSeq First_Seq,int IsCons,int IsFirst_Seq ,
	t_indL Start,t_indL Finish,t_pParamAlign Param,int BlockSize,t_Style OutputStyle)
/* ecrit BlockSize residus de la sequence p_Seq */
{


	int c,c1,c2,h,h1,Colour,OldColour;
	t_indL Pos;
	int  Count=0;
	t_Seq *S,*C,*F;
	char Gap[4]="  .";

	OldColour = 0;
	Gap[1]=Insert;
	for( Pos=Start, S=Seq+Pos,C=Cons+Pos,F=First_Seq+Pos; Pos<=Finish; Pos++,S++,C++,F++ )
	{
		c= S->Car;
		h= Param->Hom[Param->NumSymb[toupper(c)]];
		c1= C->Car;
		if (strchr(Gap,c1)) h1=c1;
		else h1= Param->Hom[Param->NumSymb[toupper(c1)]];
		c2= F->Car;
		if (!strchr(Gap,c)) /* le caractere a afficher n'est pas le '-' */
		{
			if (h == h1) /* residu = consensus */
				if (!islower(c1)) /* consensus fort */
					Colour=2;
				else /* consensus faible */
					Colour=1;
			else if (IsCons) Colour = 1;
			else/* residu != consensus */
			{
				Colour = 0;
				if (OutputStyle==Case) c = tolower(c);
			}
		}
		else
		{
			Colour=0;
		}
		if ((OutputStyle==Difference) && (!IsFirst_Seq) && (!strchr(Gap,c)))
		{
			if (c == c2)
				c = '.';

		}
		if (Colour!=OldColour)
		{
			if (OldColour) fputc(Ferme[OldColour],outf);
			if (Colour) fputc(Ouvre[Colour],outf);
			OldColour = Colour;
		}
		fputc(c,outf);
		Count++;
		if ( Count==BlockSize && Pos!=Finish)
		{
			if (OldColour) fputc(Ferme[OldColour],outf);
         OldColour = 0;
			fputc(' ',outf);
			Count= 0;
		}
	}
	if (OldColour) fputc(Ferme[OldColour],outf);
	fputc('\n',outf);
}


