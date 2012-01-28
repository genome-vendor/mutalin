#include <stdlib.h>
#include <stdio.h>
#ifdef __BORLANDC__
 #include <io.h> /* non ANSI : for use of isatty */
#endif
#include <stdarg.h>
#include <string.h>
#include "const.h"
#include "msgerr.h"

typedef struct _StringListElem {
	int Key;
	char *Text;
} StringListElem;

typedef struct _TStringList {
	int size;
	StringListElem *Index;
} TStringList;

static StringListElem Errors[16]= {
{0,"Error #%ld."},
{8,"Not enough memory to perform this action. "},
{18,"Not enough memory to make an insertion. "},
{23,"Not enough memory to calculate consensus. "},
{28,"Error reading %s. "},
{29,"Error writing %s. "},
{3,"No sequence loaded from %s . "},
{4,"Two sequences have the same name : %s.\nThe second one is discarded. "},
{5,"This sequence is too long.\nIt is cut at %ld. "},
{6,"Sequence %s is too short.\nIt is discarded.\n "},
{7,"Too many sequences.\nOnly %ld are kept. "},
{17,"No symbol comparison table selected. "},
{12,"%s contains no symbol comparison table. "},
{11,"Cannot make an insertion :\n%s becomes too long. "},
{16,"Cannot retrieve temporary file.\nCancel action. "},
{34,"%s is invalid in the configuration file. "}
};

static StringListElem Warnings = {50,"%c exists in a sequence and\nnot in the symbol comparison table.\nDelete from sequence ? "};

static StringListElem Actions[9] ={
{90,"Reading "},
{91,"Saving "},
{92,"Printing "},
{93,"Printing to "},
{94,"Saving Configuration."},
{95,"Reading Configuration."},
{96,"Restoring previous iteration."}
};

static TStringList Error[4] = {{16,Errors},{1,&Warnings},{0,NULL},{3,Actions}};
static TStringList *Erreur = Error;

int Muet;
int SortieConsole = !0;
/* Declaration of the output file for HTML errors */
#ifdef HTML
	static FILE *ErrorFile = NULL;
#else
#	define ErrorFile stdout
#endif


void Usage()
/* affiche la syntaxe d'utilisation de MA */
{
#ifdef HTML
		  /*open output file */
		  if ((ErrorFile = fopen("error_file.txt", "w"))==NULL) exit(1);
		  fprintf(ErrorFile,"Parameters on command line are incorrect !\n");
		  fclose (ErrorFile);
#else
	printf("\n\n");
	printf("Usage:  (1) ma or ma [<Option> ...] <SequencesFile> [<ClusterFile.xxx>]\n");
	printf("        (2) ma [<Option>...] <AlignmentFile.mul>\n");
	printf("Options are:\n");
	printf("   -r                     Recover last configuration.\n");

	printf("   -i:<Input format>       with Format:  gcg, embl, genbank,\n");
	printf("                           mul (MultAlin=Pearson),auto(matically determined).\n");
	printf("   -p                    select Parts of sequences.\n");

	printf("   -c:<symbol Comparison table>\n");
	printf("   -g:<Gap value> with gap_penalty = gap_value + gap_length_value x length(gap)\n");
	printf("   -l:<Gap length value>\n");
	printf("   -x:<Gap at ext>       0 no penalties, 1 penalty at the end\n");
	printf("                         2 penalty at the beginning, 3 at both extremities\n");

	printf("   -1                    One iteration only.\n");
	printf("   -u                    Unweighted sequences.\n");
	printf("   -s:<Scoring method>   with Scor meth: abs(olute), per(centage), ide(ndity)\n");

	printf("   -o:<Output format>    with Format:  msf (GCG), mul (MultAlin) or doc (Word).\n");

	printf("   -a                    output sequences ordered as input.\n");
	printf("   -A                                                 aligned.\n");
#ifdef _gif
	printf("   -f                    create a GIF image.\n");
#endif
	printf("   -d                    Draw the clustering in the output clustering file.\n");
	printf("   -k:<U>.<L>            U and L are minimal %% for Upcase and Lowcase consensus\n");
	printf("   -q                    quiet : no message during alignment\n");
	printf("Default options: -i:auto -c:%s -g:12 -l:2 -x:0 -s:abs -o:msf -A -k:90.50\n\n",DEF_TAB);
	printf(" if <ClusterFile.xxx> is not present, first clustering with FAST\n");
	printf(" if xxx = clu, first clustering read from <ClusterFile.clu>\n");
	printf(" if xxx = sco, scores read from <ClusterFile.sco> give the first clustering\n");
	printf(" if <SequencesFile> has .mul extension syntax (2) is used\n");
	printf("Syntax (2): ma [<Option> ...] <AlignmentFile.mul>\n");
	printf(" no alignment, new outputs with options -o, -f, -k, -q\n");
	printf("Special option for syntax (2):\n");
	printf("    -2:<n>	            alignment of profile 1 (aligned sequences 1 to n-1) with\n");
	printf("                         profile 2 (aligned sequences n to end of input file).\n");

#endif
	exit(1);
}


void InitMessages(char * ProgName)
{
#ifdef HTML
 /*open output file */
 if ((ErrorFile = fopen("error_file.txt", "w"))==NULL) ErrorFile = stdout;
/* close error file */
 if (ErrorFile!=stdout) fclose (ErrorFile);
#endif
#ifdef __BORLANDC__
 SortieConsole = isatty (fileno(stdout));
#endif
}

void DoneMessages(void)
{
}

void SetMessageMode (int isQuiet)
{
	Muet = isQuiet;
}

static void GetMessage (char *Mes, int NoErr)
{
	int i;
	StringListElem *x;
	Mes[0]=0;;
	for (i=0, x=Erreur->Index; (i<Erreur->size)&&
		(NoErr != x->Key ); i++, x++);
	if (i < Erreur->size)
		strcpy (Mes, x->Text);
}

unsigned int TraiteErr(int NoErr, ...)
{
	char Mes[256],S[256];
	va_list ap;
	TStringList *VraiErreur = Erreur;

#ifdef HTML
		  /* open error file */
if ((ErrorFile = fopen("error_file.txt", "a"))==NULL) ErrorFile = stdout;
#endif
	va_start (ap, NoErr);
	switch((int)(NoErr / 10)) {
		case 0:
		case 1:
		case 2:
		case 3:
		case 4:
			fprintf(ErrorFile,"\nError: ");break;
		case 5 :
			fprintf(ErrorFile,"\nWarning: ");Erreur++;break;
		default:
			fprintf(ErrorFile,"\nInformation: ");Erreur+=2;
	}

	GetMessage (Mes,NoErr);
	if (Mes == "")
	{
		Erreur = VraiErreur;
		GetMessage (Mes,0);
		sprintf (S, Mes, NoErr);
	}
	else vsprintf (S, Mes, ap);
	fprintf(ErrorFile,"%s\n",S);

	va_end (ap);
		  /* close error file */
	if (ErrorFile!=stdout)	  fclose (ErrorFile);
	Erreur= VraiErreur;
	return cmOK;
}

void MessageAction(int NoMess,char *P)
{
	char S[256];

	if (Muet) return;
#ifdef HTML
		  /* open error file */
if ((ErrorFile = fopen("error_file.txt", "a"))==NULL) return;
#endif
	strcpy (S, Erreur[3].Index[NoMess-90].Text);
	if (P != NULL) strcat (S,P);
	fprintf(ErrorFile,"%s\n",S);
#ifdef HTML
		  /* close error file */
		  fclose (ErrorFile);
#endif
}

