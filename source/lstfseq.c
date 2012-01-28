/*------------------------------------------------------------------------
	LSTFSEQ.C
			Gestion de la liste chainee des noms
				de fichiers de sequences
-------------------------------------------------------------------------*/

#include <string.h>
#include <ctype.h>

#include "disk.h"
#include "msgerr.h"

/*=============  LISTES CHAINEES  =======================================*/
/* liste des noms de fichiers de sequences */
typedef struct s_ElemListeNomFich
{
	t_CheminFichier Nom;
	struct s_ElemListeNomFich* Suiv;
} t_ElemListeNomFich;
typedef t_ElemListeNomFich* t_pElemListeNomFich;

typedef struct
{
	t_pElemListeNomFich Tete;
	t_pElemListeNomFich Courant;
} t_ListeNomFich;
typedef t_ListeNomFich* t_pListeNomFich;

static t_ListeNomFich ListeFich;






/*-----------------------------------------------------------------------*/
static int MetDansListeNomFich( char* Nom );
int InitListeNomFich( char * NomFichFich )
/* initialise la liste avec les noms de fichiers de sequences
	contenus dans le fichier NomFichFich */
{
	FILE* FICH;
	t_CheminFichier Nom;
	int c;

	ListeFich.Tete= NULL;
	ListeFich.Courant=NULL;
	if((FICH= fopen(NomFichFich,"rt")) ==NULL)
	{
		TraiteErr( 28, NomFichFich );
		fclose( FICH );
		return faux;
	}
	while( !feof(FICH) )
	{
		while( isspace(c=fgetc(FICH)) && c!=EOF );
		if( c==EOF )
			break;
		ungetc( c, FICH );
		fscanf( FICH, "%s", Nom );
		if( !MetDansListeNomFich(Nom ) )
		{
			fclose( FICH );
			return faux;
		}
	}
	fclose( FICH );
	return vrai;
}

/*-----------------------------------------------------------------------*/

static int MetDansListeNomFich( char* Nom )
/* ajoute en tete dans la liste chainee des noms de fichiers de sequences */
{
	t_pElemListeNomFich Nouveau;

	if( (Nouveau= (t_pElemListeNomFich)malloc(sizeof(t_ElemListeNomFich)))==NULL )
		return faux;

	/* initialiser le nouvel element */
	strcpy( Nouveau->Nom, Nom );
	Nouveau->Suiv= ListeFich.Tete;

	/* l'inserer */
	ListeFich.Tete= Nouveau;
	return vrai;
}

/*-----------------------------------------------------------------------*/

void DetruitListeNomFich( void )
/* detruit les elements chaines depuis Tete,
		jusqu'a la fin de la liste*/

{
	t_pElemListeNomFich Suivant;
	t_pElemListeNomFich Courant= ListeFich.Tete;

	while( Courant )
	{
		Suivant= Courant->Suiv;
		free( Courant );
		Courant= Suivant;
	}
	ListeFich.Tete= NULL;
	ListeFich.Courant= NULL;
}

char * GetNextFileinList (void)
{
	t_pElemListeNomFich Fich;
	if (!ListeFich.Courant) ListeFich.Courant = ListeFich.Tete;
	Fich = ListeFich.Tete;
	if (! Fich)
	{
		ListeFich.Tete = ListeFich.Courant;
		return NULL;
	}
	ListeFich.Tete=Fich->Suiv;
	return Fich->Nom;
}
