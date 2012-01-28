/*------------------------------------------------------------------------
	Fichier SWAP.C

			Gestion de swap disque pour une matrice d'entiers
					( gestion independante de MA )
------------------------------------------------------------------------*/
/* 4/10/95 les blocs de Swap se fon a partir de la fin de la matrice et font
	toujours exactemnet NbEnrg lignes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "swapd.h"

#ifndef vrai
#	define vrai 1
#	define faux 0
#endif

typedef struct
{
	int Swap;
	int DerSwap;			/* no du dernier fichier de swap */
	size_t NbEnrg;	 		/* nbre de lignes de Direction par fichier swap */
	size_t NbEnrgTotal;		/* nbre total d'enregistrements de la matrice */
	size_t TailleEnrg;     	/* taille en octets d'un enregistrement */
	char PathTmp[256];  	/* path du repertoire tmp s'il existe */
}t_DescripteurSwapDisk;

static t_DescripteurSwapDisk DSD;


typedef char t_NomFich[15];
typedef char t_Path[256];

/*=======================================================================*/
static void NomFichSwapD( char* PathTmp, int NoSwap, char* Path )
/* determine le nom du NoSwap ieme fichier du repertoire PathTmp.
	En sortie, Path contient le chemin complet d'acces a ce fichier */
{
	t_NomFich NomFich;

	sprintf (NomFich,"\\MA#%05d.!@~",NoSwap);
	strcpy( Path, PathTmp );
	strcat( Path, NomFich );
}


int InitSwapD( size_t NbEnrg, size_t NbEnrgTotal, size_t LgEnrg )
/* initialise les parametres de swapping */
{
	DSD.Swap = SWAPIMPOSSIBLE;
	if( !NbEnrg || !LgEnrg ) return SWAPIMPOSSIBLE;

	DSD.DerSwap= 0;
	DSD.NbEnrgTotal= NbEnrgTotal;
	DSD.NbEnrg= NbEnrg;
	DSD.TailleEnrg= LgEnrg;

	if (NbEnrg >= NbEnrgTotal)
	{
		DSD.Swap = NOSWAP;
		return NOSWAP;
	}

	/* recherche si un repertoire tmp existe,
		sinon PathTmp est vide (= repertoire courant) */
	/* DOS # UNIX */
	if( getenv("TEMP") )
		strcpy( DSD.PathTmp, getenv("TEMP") );
	else
		DSD.PathTmp[0]= '\0';
	DSD.Swap = SWAP;
	return SWAP;
}




int SwapD( void** Matrice , size_t RangPremEnrg)
/* realise un swap disque de DSD.NbEnrg enregistrements de Matrice,
	a partir de l'element d'indice fictif DebSwap */
{
	t_Path Path;

	FILE* F;
	size_t i, imax;
	int NoSwap; 	/* no du fichier de swap a creer */
	void **Ligne;

	if(DSD.Swap != SWAP)	return DSD.Swap;

	/* determination du numero du fichier de swap a creer */
/*	NoSwap= (RangPremEnrg+DSD.NbEnrg-1)/DSD.NbEnrg;*/
/* 4/10/95 on compte dans l'autre sens */
	NoSwap = (DSD.NbEnrgTotal -1 - RangPremEnrg)/DSD.NbEnrg;

	/* creer le nouveau fichier de swap */
	NomFichSwapD( DSD.PathTmp, NoSwap, Path );
	if( !(F=fopen( Path, "wb" )) )
		return SWAPIMPOSSIBLE;

	/* mettre a jour DerSwap */
	if( NoSwap > DSD.DerSwap )
		DSD.DerSwap= NoSwap;

	/* !! ne pas utiliser fwrite en une seule ecriture :
	les lignes de matrice peuvent ne pas etre consecutives en memoire */
/* 4/10/95 les paquets font toujours NbEnrg
	imax = DSD.NbEnrgTotal - NoSwap*DSD.NbEnrg;
	if (imax > DSD.NbEnrg)*/ imax = DSD.NbEnrg;
	for( i=0, Ligne=Matrice; (i<imax); i++,Ligne++ )
		if( fwrite( *Ligne, DSD.TailleEnrg, 1, F ) != 1 )
		{
			fclose(F);
			return SWAPIMPOSSIBLE;
		}
	fclose( F );
	return SWAP;
}


int RecupSwapD( void** Matrice, size_t RangEnrg )
{
	t_Path Path;
	FILE* F;
	size_t i, imax;
	int NoSwap;	/* no du fichier de swap a lire */
	void **Ligne;


	if( DSD.Swap != SWAP) return DSD.Swap;

	/* recherche du fichier de swap correspondant a RangPremEnrg */
/*	NoSwap= RangEnrg / DSD.NbEnrg;*/
/* 4/10/95 on compte dans l'autre sens */
	NoSwap = (DSD.NbEnrgTotal -1 - RangEnrg) / DSD.NbEnrg;

	/* lire le fichier de swap */
	NomFichSwapD( DSD.PathTmp, NoSwap, Path );

	if( !(F=fopen( Path, "rb" )) )
		return SWAPIMPOSSIBLE;

	/* lecture */
	/* !! ne pas utiliser fread en une seule ecriture :
		les lignes de matrice peuvent ne pas etre consecutives en memoire */
/* 4/10/95 les paquets font toujours NbEnrg
	imax = DSD.NbEnrgTotal - NoSwap*DSD.NbEnrg;
	if (imax > DSD.NbEnrg)*/ imax = DSD.NbEnrg;
	for( i= 0, Ligne=Matrice; (i<imax); i++,Ligne++ )
		if( fread( *Ligne, DSD.TailleEnrg, 1, F ) != 1 )
		{
			fclose(F);
			return SWAPIMPOSSIBLE;
		}
	fclose( F );
	return SWAP;
}


int DetruitSwapD( void )
/* supprime les fichiers de swap de numeros 0 a DerSwap,
	du repertoire PathTmp */
{
	int i;
	t_Path Path;
	int Res= SWAP;

	if (DSD.Swap != SWAP) return DSD.Swap;
	for( i=0; i<=DSD.DerSwap; i++ )
	{
		NomFichSwapD( DSD.PathTmp, i,  Path );
		if( remove(Path) )
		{
			Res= SWAPIMPOSSIBLE;
		}
	}
	DSD.NbEnrg = 0;
	return Res;
}





