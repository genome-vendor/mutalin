/*-------------------------------------------------------------------------
fichier UTIL.H :
			Procedures de manipulation d'ensembles
				 de chaines de caracteres

--------------------------------------------------------------------------*/

#ifndef __UTIL_H

#define __UTIL_H

/* Traitements sur les chaines de caractere */
  /*-------- General purpose string manipulation --------*/
char *strupr (char *S);

char *strlwr (char *S);

	/* File name string manipulation */
char *NomSeul (char* NomFich);
char* ChangeExt( char* Fich, char* NlleExtension );

	/* basic file manipulation */
void Backup( char* NomFich );
/* create a backup file for NomFich */


#endif	/* __UTIL_H */
