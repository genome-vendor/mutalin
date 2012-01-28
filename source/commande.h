#ifndef COMMANDE_H
#define COMMANDE_H
/*=============  PARAMETRES D'EXECUTION  ================================*/
#include "deftypes.h"
typedef enum Liste_Styles {Normal, Case , Difference } t_Style;

/* informations de la ligne de commande */
int PECanSelect (void);
/* return 1 if the user can be asked for partial sequences */

char * PEGetCoeffFileName (void);
/* return the coeffivient table file name */

int PEGetFormatArbre (void);
/* return 1 if the tree must be drawn, 0 for a parantheses format */

int PEGetFormatEntree( FILE* FICH);

int PEGetKindOfTreeInput (void);
/* return 1 if there is a tree file, 2 if there is a score file and -1 if no file */

char * PEGetNomFichSeq (void);
/* return the sequence file name */

int PEGetOutputFormat (void);
/* return the output format */

void PEGetOutputStyleParam(t_Style *Astyle, int *Alinesize, int *Ablocksize);
/* return style info for doc and gif output */

int PEGetSeqDejaAligne (void);
/* return as PEUpDateMul */

char * PEGetTreeFile (void);
/* return the tree or score file name */

int PEGifOutput (void);
/* return 1 if a gif image is wanted */

void PEInit (void);
/* give default values to the different options */

int PEIsList (void);
/* return 1 if the input file is a list of sequence file names */

int PEIsOutputOrder (void);
/* return 1 if output order must follow the tree order */

int PEIsQuiet (void);
/* return 1 if option -q (quiet) is set */

int PEReadSeqCoeff (t_pDescripteurSequence DS, t_pParamAlign Param);
/* Read Sequence file(s) and Coefficient file */

void PEMajParam (t_pParamAlign Param);
/* Set Param with the user defined values */

void PEModeInteract (void);
/* set  the options with an interactive dialog */

int PEModeLigneCmde (int Argc, char * Argv[]);
/* set the options by reading the command line */

void PESauveConfig (void);
/* Write the configuration to disk */

void 	PESaveInputIfMul(void);
/* rename the input file if its extension is mul */

void PESelectPartSeq( char* NomSeq, long *Debut, long *Fin );
/* ask the user for start and end position in a sequence */

void PESetOutputOrder (void);
/* Force the output order to follow the tree */

int PEUpDateMul (int NbSeq);
/* return 0 in the normal case; return 1 if sequences are already aligned,
return the index of the the second group first sequence for the alignment of
two profiles; if this index is invalid, return 1 */

#ifdef _madomain
char * PEGetDomainOutputFile(void);
int PENoOutput (void);
void 	PESaveInputIfMul(void);
#endif
#endif
