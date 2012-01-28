/*-------------------------------------------------------------------------
	fichier DISK.H
			Procedures de traitement sur les fichiers

-------------------------------------------------------------------------*/

#ifndef __DISK_H
#	define __DISK_H
#include "deftypes.h"
/* dans disk.c */
int EcritSeq( t_pDescripteurSequence DS, t_pParamAlign Param );
int LitSeqCoeff( t_pDescripteurSequence DS, t_pParamAlign Param);
/* dans msfdisk.c */
extern char ConsSymb[DERLETTREMAX+1][DERLETTREMAX+1];
int CalculeSymboles (t_pParamAlign Param);
void MsfEcritSeq( FILE* OUTFILE ,t_pDescripteurSequence DS,
			t_pParamAlign Param, int style, int Group2 );
/* dans muldisk.c */
int MulLitSeq( FILE* FICH, t_pDescripteurSequence DS, int Partiel,
	char* EnsLet );
void MulEcritSeq( FILE* F, t_pDescripteurSequence DS );
void Ecrit1Seq(FILE* f, t_pLigneSeq Sequence, float Weight0);
/* dans mbgbdisk.c */
int EmblGenbankLitSeq( FILE* FICH, t_pDescripteurSequence DS,
			t_FormatFichSeq Format, char* NomFich, int Partiel, char* EnsLet );
/* dans gcgdisk.c */
int GcgLitSeq( FILE* FICH, t_pDescripteurSequence DS,
					char* NomFich, int Partiel, char* EnsLet );
/* dans lstfseq */
int InitListeNomFich( char * NomFichFich );
char * GetNextFileinList (void);
void DetruitListeNomFich( void );
/* dans coefdisk.c */
int LitFichCoeff (FILE *FICH, char* FileName, t_pParamAlign Param);
/* dans image.c */
void CreerImage(char* FileName ,t_pDescripteurSequence DS,t_pParamAlign Param,
	int complete,int Group2 );

#endif
