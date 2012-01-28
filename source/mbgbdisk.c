/*-----------------------------------------------------------------------
	MBGBDISK.C

		Lecture des fichiers de sequences au format EMBL et GenBank
		( 1 ou plsieurs sequences par fichier, 1 ou plusieurs fichiers )

                                  EMBL                  GenBank
                                 -------------------------------
Nom de sequence precede de:       ID                    LOCUS
Debut de sequence:                SQ                    ORIGIN
Fin de sequence:                  //                    //
-----------------------------------------------------------------------*/

/* CORRECTIONS

1 juillet 1994 : MajTailleConsens has a third parameter "complet" that must be
	false when the sequences are not yet aligned

*/

#include <string.h>
#include <ctype.h>
#include "disk.h"
#include "afichseq.h"
#include "msgerr.h"
#include "parametr.h"
#include "commande.h"


static int ChercheRubrique( FILE* FICH, char* Rubrique );

int EmblGenbankLitSeq( FILE* FICH, t_pDescripteurSequence DS,
			t_FormatFichSeq Format, char* NomFich, int Partiel,char* EnsLet );

/*---------------------------------------------------------------------*/
static int ChercheRubrique( FILE* FICH, char* Rubrique )
/* lit la ligne courante du fichier jusqu'a '\n' inclus,
	Rubrique recoit le premier mot de la ligne */
{
	register int c;

	if( fscanf( FICH, "%s", Rubrique ) == EOF )
		return EOF;
	/* passer la fin de ligne */
	while( (c=fgetc(FICH)) != '\n' )
		if( c == EOF )
			return EOF;
	return vrai;
}


/*---------------------------------------------------------------------*/
int EmblGenbankLitSeq( FILE* FICH, t_pDescripteurSequence DS,
			t_FormatFichSeq Format, char* NomFich, int Partiel, char* EnsLet )
{
	int c;
	t_Seq *S;
	t_pLigneSeq tmp;
	long Cpt;
	long DebLect= 0;
	long FinLect=  LONG_MAX;

	char Rubrique[256];
	char RubIdent[10];
	char RubDebut[10];
	int PremiereSequence, FinLecture, FinSequence, doublet;
	int LgRubIdent, LgRubDebut;

	/* selection des mots cles (rubriques) du fichier selon le format */
	switch( Format )
	{
		case EMBL:
			strcpy( RubIdent, "ID" );
			strcpy( RubDebut, "SQ" );
			LgRubIdent= 2;
			LgRubDebut= 2;
			break;
		case GENBANK:
			strcpy( RubIdent, "LOCUS" );
			strcpy( RubDebut, "ORIGIN" );
			LgRubIdent= 5;
			LgRubDebut= 6;
			break;
	}


	/* lecture des sequences */
	PremiereSequence= vrai;
	FinLecture= faux;
	while( !FinLecture && (DS->NbSeq<NbSeqMax))
	{

		/*----- passer les espaces de debut ----*/
		while( (c= fgetc(FICH))!=EOF && isspace(c) );
		if( feof(FICH) )
			if( PremiereSequence )
			{
				TraiteErr( 3, NomFich );
				return faux;
			}
			else
			{
				FinLecture= vrai;
				break;
			}

		PremiereSequence= faux;
		ungetc( c, FICH );

		/*------- allouer 4 sequences de plus si necessaire -------*/
		if(!(DS->NbSeq %4) && !(DS->Sequence=(t_pLigneSeq)realloc(DS->Sequence,
				(DS->NbSeq+5)*sizeof(t_LigneSeq))) )
		{
			TraiteErr(8);
			return faux;
		}
		AffichUpDate (++DS->NbSeq);

		/*----- lire le nom de la sequence -----*/
		fscanf( FICH, "%s", Rubrique );
		if( strncmp(Rubrique,RubIdent,LgRubIdent) )
		{
			TraiteErr( 28, NomFich );
			return faux;
		}
		else
			fscanf( FICH, "%s", DS->Sequence[DS->NbSeq].NomSeq );

		/*----- verifier que ce nom de ce sequence n'a pas encore
			ete lu --------------------------------------------*/
		if( doublet = ExisteSequence( DS, DS->NbSeq ) )
			TraiteErr( 4, DS->Sequence[DS->NbSeq].NomSeq );

		/*----- selectionner la portion de sequence a utiliser ----*/
		if( Partiel && !doublet)
			PESelectPartSeq(DS->Sequence[DS->NbSeq].NomSeq,
						&DebLect,&FinLect);

		/*----- passer les commentaires ------*/
		while( ChercheRubrique( FICH, Rubrique ) != EOF
			&& strncmp( Rubrique, RubDebut, LgRubDebut ) );
/*		if( feof(FICH) )
		{
			TraiteErr( 28, NomFich );
			return faux;
		} */

		/*------ lire la sequence -------------*/
		tmp= DS->Sequence+DS->NbSeq;

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


		/* lecture de la sequence */
		FinSequence= faux;
		for( Cpt=1; (c=fgetc(FICH))!=EOF && !FinSequence && (tmp->Taille<LongMax);
			 Cpt++ )
			if((c>' ') &&  !isspace(c) && !isdigit(c) && c!='/' )
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
					S=tmp->Seq+tmp->BuffLen;
					memset(S , 0, 64*sizeof(t_Seq) );
					tmp->BuffLen+=64;
				}
				if( Cpt>=DebLect && Cpt<=FinLect )
				{
					c=toupper(c);
					AjouteSymb (EnsLet,c);
					S->Car= c;
					S++;
					tmp->Taille++;
				}
			}
			else
			{
				if( c=='/' ) 
					if( (c=fgetc(FICH))=='/' )
						FinSequence= vrai;
					else
					{
						TraiteErr( 28, NomFich );
						return faux;
					}
			}

		/* maj des caracteristiques de la sequence */
		if( !doublet )
		{
			if (tmp->Taille < 3)
			{
				TraiteErr(6,tmp->NomSeq);
				doublet = vrai;
			}
			else if (tmp->Taille == LongMax)
				TraiteErr (5,Cpt);
		}

		if (doublet)
		{
			free (tmp->Seq);
			AffichUpDate (--DS->NbSeq);
		}
		else
		{
			tmp->PremCol= 1;
			tmp->DerCol= tmp->Taille;
			tmp->NbResid= tmp->Taille;
         tmp->Weight = 1;
		}
	}
	if (!FinLecture) TraiteErr (7,(long)NbSeqMax);
	return vrai;
}

