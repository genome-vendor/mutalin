/*-------------------------------------------------------------------------
	fichier MULDISK.C

		Lecture et ecriture au format MultAlin

		( 1 seul fichier contenant plusieurs sequences )

Commentaires:  '>' en debut de ligne et en dehors des sequences



-------------------------------------------------------------------------*/
/* CORRECTIONS

1 juillet 1994 : MajTailleConsens has a third parameter "complet" that must be
	false when the sequences are not yet aligned

*/

#include <string.h>
#include <ctype.h>
#include "disk.h"
#include "afichseq.h"
#include "commande.h"
#include "msgerr.h"
#include "parametr.h"


/*=============    Lecture de fichier de sequence   ================*/
static void LireNomSequence( FILE* FICH, t_pDescripteurSequence DS, t_indN NoSeq )
{
	char buf[256], *s;
	char format[6];
	float w;
	sprintf(format,"%%%ds",LONGNOM-1);
	fgets (buf,255,FICH);
	/* lecture du nom */
	sscanf(buf,format,DS->Sequence[NoSeq].NomSeq);
	/* recherche du poids s'il existe */
	if ((s=strstr(buf,"Weight:"))==NULL) DS->Sequence[NoSeq].Weight = 100;
	else
	{
		sscanf(s+7,"%f",&w);
		DS->Sequence[NoSeq].Weight = (int)(w*100.0 + 0.5);
	}
	/* passer sur la fin de la ligne */
}


/*================  Lecture au format MultAlin  ==================*/

int MulLitSeq( FILE* FICH, t_pDescripteurSequence DS, int Partiel,char* EnsLet)
{
	int SeqDejaAligne = PEGetSeqDejaAligne();
	t_pLigneSeq tmp;
	t_Seq *S;
	long DebLect= 0;
	long FinLect= LONG_MAX;
	long Cpt;
	int c;
	t_indN i;
   t_indL l;
	int doublet;

	for(i=DS->NbSeq+1, c=0; !feof(FICH) && (DS->NbSeq <NbSeqMax);
		i++ ,AffichUpDate(++DS->NbSeq))
	{
		/* passer les commentaires */
		if( c != '>' )
			while( ((c=fgetc(FICH))!=EOF) && (c!='>') );

		/*-------------- lire le nom de la sequence ----------------*/
		if( !feof(FICH) )
			/* allouer 4 sequences de plus si necesaire */
			if( !((i-1)%4) && !(DS->Sequence=(t_pLigneSeq)realloc(DS->Sequence,
					(i+4)*sizeof(t_LigneSeq))) )
			{
				TraiteErr(8);
				return faux;
			}
		LireNomSequence( FICH, DS, i );
		/*----- verifier que ce nom de ce sequence n'a pas encore
			ete lu --------------------------------------------*/
		if( doublet = ExisteSequence( DS, i ) )
			TraiteErr( 4, DS->Sequence[i].NomSeq );

		/*----- selectionner la portion de sequence a utiliser ----*/
		else if( Partiel )
			PESelectPartSeq(DS->Sequence[i].NomSeq,&DebLect,&FinLect);

		/*-------------------- lire la sequence -------------------*/
		if( !feof(FICH) )
		{
			tmp= DS->Sequence+i;

			/* allouer 64 positions de la chaine Seq */
			if(! (tmp->Seq=(t_pSeq)calloc(64,sizeof(t_Seq))) )
			{
				TraiteErr(8);
				return faux;
			}
			memset( tmp->Seq, 0, 64*sizeof(t_Seq) );
			S = tmp->Seq+1;
			tmp->BuffLen=64;
			tmp->Taille= 0;

			/* lecture */
			for( Cpt=1; (c=fgetc(FICH))!=EOF && c!='>' && (tmp->Taille < LongMax);
				 Cpt++ )
				if((c>' ') &&  !isspace(c) && !isdigit(c) && (SeqDejaAligne || (c!=Insert)))
				{
					/* reallouer 64 positions de plus si necessaire */
					if(tmp->Taille+1 == tmp->BuffLen)
					{
						if(!(tmp->Seq=(t_pSeq)realloc(tmp->Seq,
								(tmp->BuffLen+64)*sizeof(t_Seq))) )
						{
							TraiteErr(8);
							return faux;
						}
						S= tmp->Seq + tmp->BuffLen;
						memset( S, 0, 64*sizeof(t_Seq) );
						tmp->BuffLen+=64;
					}

					if( Cpt>=DebLect && Cpt<=FinLect )
					{
						c= toupper(c);
						AjouteSymb (EnsLet, c);
						S->Car= c;
						S++;
						tmp->Taille++;
					}
				}

			if (!doublet)
			{
				if( tmp->Taille < 3 )
				{
					TraiteErr( 6, tmp->NomSeq );
					doublet = vrai;
				}
				else if (tmp->Taille >= LongMax)
					TraiteErr (5, (long)LongMax);
			}

			if (doublet )
			{
				free (tmp->Seq);
				DS->NbSeq--;
				i--;
			}
			else
			{
			/* maj des caracteristiques de la sequence */
			if (SeqDejaAligne)
			{
				for (l=1;tmp->Seq[l].Car==Insert;l++) tmp->Seq[l].Car = ' ';
				tmp->PremCol= l;
				for (Cpt=0;l<=tmp->Taille;l++) if (tmp->Seq[l].Car!=Insert)
				{
					Cpt++;
					tmp->DerCol = l;
				}
				for (l= tmp->DerCol+1;l<=tmp->Taille;l++) tmp->Seq[l].Car=' ';
				tmp->NbResid= Cpt;
			}
			else
			{
				tmp->PremCol = 1;
				tmp->DerCol = tmp->NbResid = tmp->Taille;
			}
			}
		}
	}
	if (!feof(FICH)) TraiteErr (7,(long)NbSeqMax);
	return vrai;
}






/*================  Sauvegarde au format MultAlin  ==================*/
void Ecrit1Seq(FILE* f, t_pLigneSeq Sequence, float Weight0)
{
	int Count;
	int Pos;
	int c;

	/* ecrire le nom et le nombres de residus de la sequence */
	fprintf(f,"%c%s%10ld Weight: %3.2f\n",DebutLigne, Sequence->NomSeq,
		(long)Sequence->NbResid, Weight0*Sequence->Weight);

	Count= 0;
	Pos= 1;
#ifndef NOCOMMENT
	fprintf(f,"%6ld ",(long)Pos);
#endif
	/* ecrire eventuellement les caracteres d'insertion en debut de sequence */
	while( Pos < Sequence->PremCol )
	{
		/* 50 par ligne */
		if(Count == 50)
		{
			fputc('\n',f);
#ifndef NOCOMMENT
			fprintf(f,"%6ld ",(long)Pos);
#endif
			Count= 0;
		}
		/* separer les groupes de 10 */
#ifndef NOCOMMENT
		if( !(Count % 10) )
			fputc(' ',f);
#endif
		fputc(Insert,f);
		Count++;
		Pos++;
	}
	/* ecrire la sequence */
	while( Pos <= Sequence->DerCol)
	{
		if(Count == 50)
		{
			fputc('\n',f);
#ifndef NOCOMMENT
			fprintf(f,"%6ld ",(long)Pos);
#endif
			Count= 0;
		}
#ifndef NOCOMMENT
		if( !(Count % 10))
			fputc(' ',f);
#endif
		fputc((c=Sequence->Seq[Pos].Car)==' ' ? Insert : c,f);
		Count++;
		Pos++;
	}
	/* ecrire eventuellement les caracteres d'insertion en fin de sequence */
	while( Pos <= Sequence->Taille )
	{
		if(Count == 50)
		{
			fputc('\n',f);
#ifndef NOCOMMENT
			fprintf(f,"%6ld ",(long)Pos);
#endif
			Count= 0;
		}
#ifndef NOCOMMENT
		if( !(Count % 10))
			fputc(' ',f);
#endif
		fputc(Insert,f);
		Count++;
		Pos++;
	}
	fputs("\n\n",f);

}

void MulEcritSeq( FILE* F, t_pDescripteurSequence DS )
{
	t_indN i,*n;
	float Weight0 = ((float)DS->NbSeq) / DS->Sequence[0].Weight;

	for( i=1, n=DS->NoSeq; i<=DS->NbSeq; i++ )
		Ecrit1Seq(F,DS->Sequence+ (n ? *(++n) : i),Weight0);
}
