/*-----------------------------------------------------------------------
	GCGDISK.H

			Lecture des fichiers de sequences au format GCG

			( 1 sequence par fichier, plusieurs fichiers )

Commentaires d'entete:	terminaison par ".."
Delimiteurs de commentaires dans les sequences:	'<' ...'>' et '$'...'$'

-----------------------------------------------------------------------*/
/* CORRECTIONS

1 juillet 1994 : MajTailleConsens has a third parameter "complet" that must be
	false when the sequences are not yet aligned

*/


#include <ctype.h>

#include "disk.h"
#include "afichseq.h"
#include "commande.h"
#include "msgerr.h"
#include "parametr.h"


/*---------------------------------------------------------------------*/
int GcgLitSeq( FILE* FICH, t_pDescripteurSequence DS,
			char* NomFich , int Partiel, char* EnsLet)
{
	int c;
	t_Seq* S;
	t_pLigneSeq tmp;
	long Cpt;
	long DebLect= 0;
	long FinLect= LONG_MAX;
	int doublet = faux;

	/* position de debut de la ligne courante dans le fichier */
	long PosDebLigne= 0;
	/* position de fin des commentaires */
	long PosDebSeq= 0;

	int FinCommentaires= faux;

	if (DS->NbSeq == NbSeqMax)
	{
		TraiteErr (7,(long)NbSeqMax);
		return vrai;
	}

	/*----- passer les commentaires ------*/
	while( !FinCommentaires )
	{
		/* lire jusqu'a un '.' */
		while( (c= fgetc(FICH))!=EOF && c!='.' )
			if( c=='\n' )
				PosDebLigne= ftell( FICH );

		if( c!=EOF )
			if( (c= fgetc(FICH))==EOF || c=='.' )
				FinCommentaires= vrai;
	}

	if( feof(FICH) )
	{
		TraiteErr( 3, NomFich );
		return faux;
	}

	PosDebSeq= ftell( FICH );

	/*------- allouer 4 sequences de plus si necessaire -------*/
	if( !(DS->NbSeq %4) && !(DS->Sequence=(t_pLigneSeq)realloc(DS->Sequence,
			(DS->NbSeq+5)*sizeof(t_LigneSeq))) )
	{
		TraiteErr(8);
		return faux;
	}
	AffichUpDate (	++DS->NbSeq);


	/*------- lire le nom de la sequence ----------*/
	fseek( FICH, PosDebLigne, 0 );
	fscanf(FICH, "%s", DS->Sequence[DS->NbSeq].NomSeq );

	/*----- verifier que ce nom de ce sequence n'a pas encore ete lu ---*/
	if( doublet=ExisteSequence( DS, DS->NbSeq ) )
		TraiteErr( 4, DS->Sequence[DS->NbSeq].NomSeq );

	/*----- selectionner la portion de sequence a utiliser ----*/
	if( Partiel )
		PESelectPartSeq(DS->Sequence[DS->NbSeq].NomSeq,&DebLect,&FinLect);

	/*------ lire la sequence -------------*/
	fseek( FICH, PosDebSeq, 0 );
	tmp= DS->Sequence+DS->NbSeq;

	/* allouer 64 positions de la chaine Seq */
	if(! (tmp->Seq=(t_pSeq)calloc(64,sizeof(t_Seq))) )
	{
		TraiteErr(8);
		return faux;
	}
	memset( tmp->Seq, 0, 64*sizeof(t_Seq) );
	S=tmp->Seq+1;
	tmp->BuffLen=64;
	tmp->Taille= 0;


	for( Cpt=1; !feof(FICH) && (tmp->Taille < LongMax); Cpt++ )
	{
		/* lecture */
		while( (c=fgetc(FICH))!=EOF && c!='<' && c!='$')
			if( (c>' ') && !isspace(c) && !isdigit(c) && c!='.' )
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
		/* passer les commentaires */
		while( (c=fgetc(FICH))!=EOF && c!='>' && c!='$' );
	}

	if (doublet)
	{
		free (tmp->Seq);
		AffichUpDate (--DS->NbSeq);
		return vrai;
	}

	if( tmp->Taille < 3 )
	{
		TraiteErr( 6, tmp->NomSeq);
		free (tmp->Seq);
		AffichUpDate (--DS->NbSeq);
		return vrai;
	}
	else if ((tmp->Taille == LongMax) && !feof(FICH))
		TraiteErr (5, Cpt);

	/* maj des caracteristiques de la sequence */
	tmp->PremCol= 1;
	tmp->DerCol= tmp->Taille;
	tmp->NbResid= tmp->Taille;
   tmp->Weight = 1;

	return vrai;
}

