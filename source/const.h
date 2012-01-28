/*-------------------------------------------------------------------------
	Fichier CONST.H

			Definitions des constantes utilisees
					dans MA

-------------------------------------------------------------------------*/

#ifndef __CONST_H

#define __CONST_H

#include "portab.h"

#ifndef vrai
#	define vrai 1
#	define faux 0
#endif


#define LONGNOM 21       /* nom d'une sequence */

#define DERLETTREMAX 31	/* taille maxi de l'alphabet de codage des sequences */


#define FICH_SELECT "select.ma"	/* fichier contenant la liste des fichiers
				correspondants au masque de selection */

#define DEF_TAB "blosum62.tab"

/*-------- caracteres reserves dans les fichiers --------*/
#define Insert  '-'     /* caractere d'insertion */
#define ResEg '.'       /* caractere mis a la place du residu s'il est egal a la 1e seq*/
#define DebutLigne '>'  /* debut de ligne de commantaire */
#define NoSymb -1
#define FreeSymb -127



/*------- macros min() et max() ---------------*/
#ifndef max
#	undef max
#	undef min
#	define max(a,b)    (((a) > (b)) ? (a) : (b))
#	define min(a,b)    (((a) <= (b)) ? (a) : (b))
#endif


#endif 	/* __CONST_H */
