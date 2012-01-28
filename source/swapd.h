/*------------------------------------------------------------------------
	Fichier SWAP.H

			Gestion de swap disque pour une matrice d'entiers
------------------------------------------------------------------------*/



#ifndef __SWAPD_H

#define __SWAPD_H

#define NOSWAP -1
#define SWAPIMPOSSIBLE 0
#define SWAP 1

/* a virtual array of NbEnrgTotal records of TailleEnrg ints is maintained
	one block of NbEnrg only can be in memory at one time
	each block begins at an index that is an integer multiple of NbEnrg
	the memory block is of type int** */

int InitSwapD( size_t NbEnrg, size_t NbEnrgTotal, size_t LgEnrg );
/* initialize the swapping system for an array of NbEnrgTotal records, each of
	LgEnrg bytes. NbEnrg is the number of records that can be in memory
	return SWAPIMPOSSIBLE, NOSWAP (Swap unnecessary) or SWAP */

int SwapD( void** Matrice , size_t RangPremEnrg);
/*	save NbEnrg records (or less if RangPremEnrg + NbEnrg > NbEnrgTotal) into the
	virtual array beginning with index RangPremEnrg, from Matrice */

int RecupSwapD( void** Matrice, size_t RangEnrg );
/* put NbEnrg records (or less) from the virtual array beginning with index
	RangPremEnrg, into Matrice */

int DetruitSwapD( void );
/* close the swapping system */

#endif

