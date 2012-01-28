#include <ctype.h>
#include <string.h>

#include "commande.h"
#include "msgerr.h"
#include "util.h"

#define GAPDFLT 5

/*=============  PARAMETRES D'EXECUTION  ================================*/

/* informations de la ligne de commande */

typedef struct
{
	t_CheminFichier NomFichSeq;
	t_CheminFichier NomFichArbre;

	t_FormatFichSeq FormatEntree;
	t_FormatFichSeq FormatSortie;

	t_CheminFichier NomFichCoeff;


	t_score Gap, Gap2; /*14/06/94 Gap=func(length)*/
	unsigned GapExtT:2;

	short FICHIERS;
	short SELECTPART;
	short LIREARBRE;
	t_FormatScores SCMETHOD;
	short UNEITER;
	short WEIGHTED;
	short FormatArbre;
	short OutputGif;
	short OutputOrder;
	short ConsH;
	short ConsL;
	short RECUP;
	short MUET;
	t_indN SeqDejaAligne;
	t_Style Style;
	int BlockSize;
	int LineSize;
} t_ParamExec;

static t_ParamExec ParamExec;

void PEInit(void)
{
	/* parametres par defaut */
	strcpy( ParamExec.NomFichSeq, "");
	strcpy( ParamExec.NomFichArbre, "");

	ParamExec.FormatEntree= NON_DEFINI;/* -i:auto */
#ifdef _madomain
	ParamExec.FormatSortie= MUL;/* -o:mul */
#else
	ParamExec.FormatSortie= MSF;/* -o:msf */
#endif
	strcpy( ParamExec.NomFichCoeff, DEF_TAB);/* -c:DEF_TAB */
	ParamExec.Gap= -1;/* the actual value will be read from DEF_TAB */
	ParamExec.Gap2=-1;/* the actual value will be read from DEF_TAB */
	ParamExec.GapExtT=0;/* -x:0 */

	ParamExec.FICHIERS= faux;
	ParamExec.SELECTPART= faux;
	ParamExec.LIREARBRE= faux;
	ParamExec.SCMETHOD= ABSOLU;/* -s:abs */
	ParamExec.UNEITER= faux;
	ParamExec.WEIGHTED= vrai;
#ifdef HTML
	ParamExec.FormatArbre = vrai;/* -d */
#else
	ParamExec.FormatArbre = faux;
#endif
	ParamExec.OutputGif= faux;
	ParamExec.OutputOrder= vrai;/* -A */
	ParamExec.ConsH = 90;
	ParamExec.ConsL = 50;/* -k:90.50 */
	ParamExec.RECUP= faux;
	ParamExec.MUET= faux;
	ParamExec.SeqDejaAligne = 0;
	ParamExec.Style= Normal;
	ParamExec.BlockSize= 10;
	ParamExec.LineSize= 50;
}

void PESetOutputOrder (void)
{
	ParamExec.OutputOrder = 1;
}

char * PEGetNomFichSeq (void)
{
	return ParamExec.NomFichSeq;
}

char * PEGetTreeFile (void)
{
	return ParamExec.NomFichArbre;
}

int PEGetOutputFormat (void)
{
	return ParamExec.FormatSortie;
}

char* PEGetCoeffFileName (void)
{
	return ParamExec.NomFichCoeff;
}
#ifndef _madomain
int PEIsList (void)
{
	return ParamExec.FICHIERS;
}
#endif

int PECanSelect (void)
{
	return ParamExec.SELECTPART;
}

int PEGetKindOfTreeInput (void)
{
	return ParamExec.LIREARBRE;
}

int PEGetFormatArbre (void)
{
	return ParamExec.FormatArbre;
}
#ifdef _gif
int PEGifOutput (void)
{
	return ParamExec.OutputGif;
}
#endif

int PEIsOutputOrder (void)
{
	return ParamExec.OutputOrder;
}

int PEIsQuiet (void)
{
	return ParamExec.MUET;
}

static char* ExtensionSeule( char* Fich )
{
	char* s = strrchr (Fich,'.');
	if (s) return s;
	else return "";
}

int PEUpDateMul (int NbSeq) /* set and return SeqDejaAligne */
{
	if (!strncmp (ExtensionSeule(ParamExec.NomFichSeq),".mul",4)&& !ParamExec.SeqDejaAligne)
		ParamExec.SeqDejaAligne=1;
	if (ParamExec.SeqDejaAligne > NbSeq) ParamExec.SeqDejaAligne = 1;
	return ParamExec.SeqDejaAligne;
}

int PEGetSeqDejaAligne (void)
{
	return ParamExec.SeqDejaAligne;
}

void PEGetOutputStyleParam(t_Style *Astyle, int *Alinesize, int *Ablocksize)
{
	*Astyle = ParamExec.Style;
	*Alinesize = ParamExec.LineSize;
	*Ablocksize = ParamExec.BlockSize;
}

void PESauveConfig (void)
{
	char * FichCfg;
	FILE* FICH;
	char Format[11];
	float lu;

	if( ParamExec.RECUP )   return;
#ifdef HTML
		FichCfg = ChangeExt (ParamExec.NomFichSeq,".Cfg");
#else
		FichCfg = "ma.cfg";
#endif
	Backup( FichCfg );

	if( (FICH= fopen(FichCfg,"wt"))==NULL )
		TraiteErr( 29, FichCfg );
	else
	{
		MessageAction(hcActWconfig, NULL);

		/* sauvegarde du format d'entree */
		switch( ParamExec.FormatEntree )
		{
			case GCG:
				strcpy( Format, "gcg" );
				break;
			case MUL:
				strcpy( Format, "mul" );
				break;
			case EMBL:
				strcpy( Format, "embl" );
				break;
			case GENBANK:
				strcpy( Format, "genbank" );
				break;
			case NON_DEFINI:
				strcpy( Format, "undefined" );
		}
		fprintf( FICH, "[InputFormat]\n%s\n\n", Format);

		/* sauvegarde du nom du fichier de coeff */
		fprintf( FICH, "[SymbolCompTableFile]\n%s\n\n", ParamExec.NomFichCoeff );

		/* sauvegarde de la valeur de gap */
		lu =ParamExec.Gap;
		fprintf( FICH, "[GapValue]\n%g\n\n", lu);
		lu= ParamExec.Gap2;
		fprintf( FICH, "[Gap2Value]\n%g\n\n", lu);
		/* 14/06/94 Gap2 is user defined*/
		lu = ParamExec.GapExtT;
		fprintf( FICH, "[GapExtValue]\n%g\n\n", lu);

		/* sauvegarde du nombre d'iterations */
		if( ParamExec.UNEITER )
			fprintf( FICH, "[OneIter]\ntrue\n\n" );
		else
			fprintf( FICH, "[OneIter]\nfalse\n\n" );

		/* sauvegarde du poids */
		if( ParamExec.WEIGHTED)
			fprintf( FICH, "[Weighted]\ntrue\n\n" );
		else
			fprintf( FICH, "[Weighted]\nfalse\n\n" );

		/* sauvegarde de la methode de score */
		switch (ParamExec.SCMETHOD)
		{
			case 0: strcpy(Format,"absolute"); break;
			case 1: strcpy(Format,"percentage"); break;
			case 2: strcpy(Format,"identity");break;
			case 3: strcpy(Format,"normalised");break;
		}
		fprintf( FICH, "[ScoringMethod]\n%s\n\n",Format);

		/* sauvegarde du format de sortie */
		switch( ParamExec.FormatSortie )
		{
			case MSF:
				strcpy( Format, "msf" );
				break;
			case MUL:
				strcpy( Format, "mul" );
				break;
			case DOC:
				strcpy( Format, "doc" );
				break;
		}
		fprintf( FICH, "[OutputFormat]\n%s\n\n", Format);

		/* sauvegarde de l'ordre de sortie */
		if (ParamExec.OutputOrder)
			fprintf(FICH,"[OutputOrder]\naligned\n\n");
		else fprintf(FICH,"[OutputOrder]\ninput\n\n");
#ifdef _gif
		/* sauvegarde du choix gif */
		if (ParamExec.OutputGif)
			fprintf(FICH,"[GifImage]\nyes\n\n");
		else fprintf(FICH,"[GifImage]\nno\n\n");
#endif

		/* sauvegarde du format de sortie de clustering */
		if (ParamExec.FormatArbre) strcpy (Format, "drawing");
		else strcpy (Format, "list");
		fprintf( FICH, "[ClusteringOutputFormat]\n%s\n\n", Format);

		/* sauvegarde des niveaux de consensus */
		fprintf(FICH, "[ConsensusLevel]\n%d %d\n\n",ParamExec.ConsH,ParamExec.ConsL);

		/* sauvegarde des parametres de style */
		switch (ParamExec.Style)
		{
			case Normal: strcpy(Format,"Normal"); break;
			case Case: strcpy(Format,"Case"); break;
			case Difference: strcpy(Format,"Difference"); break;
		}
		fprintf(FICH, "[OutputStyle]\n%s\n\n", Format);
		fprintf(FICH, "[LineSize]\n%d\n\n",ParamExec.LineSize);
		fprintf(FICH, "[GraduationStep]\n%d\n\n",ParamExec.BlockSize);

		fclose( FICH );
	}
}

static int RecupParametre(FILE* FICH )
{
	char Ligne[256];
	float GapFloat;
	int c;

	/* passer les espaces de debut */
	while( isspace(c= fgetc(FICH)) );
	if( c==EOF )
		return EOF;
	else
		ungetc( c, FICH );

	/* lire la ligne */
	if( !fgets( Ligne, 255, FICH ) )
		return faux;

	/*--- format d'entree ---*/
	if( !strncmp( Ligne, "[InputFormat]", 13 ) )
	{
		fscanf(FICH,"%s",Ligne);
		if( !strncmp( Ligne, "gcg", 3 ) )
			ParamExec.FormatEntree= GCG;
		else
		if( !strncmp( Ligne, "mul", 3 ) )
			ParamExec.FormatEntree= MUL;
		else
		if( !strncmp( Ligne, "embl", 4 ) )
			ParamExec.FormatEntree= EMBL;
		else
		if( !strncmp( Ligne, "genbank", 7 ) )
			ParamExec.FormatEntree= GENBANK;
		else
		if( !strncmp( Ligne, "undefi", 6 ) || !strncmp( Ligne, "auto", 4))
			ParamExec.FormatEntree= NON_DEFINI;
		else
		{
			TraiteErr( 34, Ligne );
			return faux;
		}
	}
	else

	/*--- fichier de coefficients ---*/
	if( !strncmp( Ligne, "[SymbolCompTableFile]", 21 ) )
	{
		fscanf(FICH,"%s",Ligne);
		strcpy( ParamExec.NomFichCoeff, Ligne );
	}
	else
	/*--- valeur de gap ---*/
	if( !strncmp( Ligne, "[GapValue]", 10 ) )
	{
		fscanf( FICH, "%f", &GapFloat );
		ParamExec.Gap= (t_score)GapFloat;
	}
	else
	/*14/06/94 Gap2 is user-defined */
	if( !strncmp( Ligne, "[Gap2Value]",11) )
	{
		fscanf( FICH, "%f", &GapFloat );
		ParamExec.Gap2= (t_score)GapFloat;
	}
	else
	if ( !strncmp( Ligne,"[GapExtValue]",12) )
	{
		fscanf ( FICH, "%f", &GapFloat );
		ParamExec.GapExtT= (short)GapFloat;
	}
	else

	/*--- nombre d'iterations ---*/
	if( !strncmp( Ligne, "[OneIter]", 9 ) )
	{
		fscanf(FICH,"%s",Ligne);
		if( !strncmp( Ligne, "true", 4 ) )
			ParamExec.UNEITER= vrai;
		else
		if( !strncmp( Ligne, "false", 5 ) )
			ParamExec.UNEITER= faux;
		else
		{
			TraiteErr( 34, Ligne );
			return faux;
		}
	}
	else

	/*--- poids ---*/
	if( !strncmp( Ligne, "[Weighted]", 10 ) )
	{
		fscanf(FICH,"%s",Ligne);
		if( !strncmp( Ligne, "true", 4 ) )
			ParamExec.WEIGHTED= vrai;
		else
		if( !strncmp( Ligne, "false", 5 ) )
			ParamExec.WEIGHTED= faux;
		else
		{
			TraiteErr( 34, Ligne );
			return faux;
		}
	}
	else

	/*------- methode de calcul du score ------------*/
	if ( !strncmp( Ligne, "[ScoringMethod]",15))
	{
		fscanf(FICH,"%s",Ligne);
		if ( !strncmp( Ligne, "abs", 3))
			ParamExec.SCMETHOD = ABSOLU;
		else
		if ( !strncmp( Ligne, "per", 3))
			ParamExec.SCMETHOD = PERCENT;
		else
		if ( !strncmp( Ligne, "ide", 3))
			ParamExec.SCMETHOD = IDENTITY;
		else
		if ( !strncmp( Ligne, "nor", 3))
			ParamExec.SCMETHOD = NORMALISED;
		else
		{
			TraiteErr( 34, Ligne );
			return faux;
		}
	}
	else
	/*--- format de sortie ---*/
	if( !strncmp( Ligne, "[OutputFormat]", 14 ) )
	{
		fscanf(FICH,"%s",Ligne);
		if( !strncmp( Ligne, "msf", 3 ) )
			ParamExec.FormatSortie= MSF;
		else
		if( !strncmp( Ligne, "mul", 3 ) )
			ParamExec.FormatSortie= MUL;
		else
		if( !strncmp( Ligne, "doc", 3 ) )
			ParamExec.FormatSortie= DOC;
		else
		{
			TraiteErr( 34, Ligne );
			return faux;
		}
	}
	else

	/* --- ordre des sequences en sortie ---*/
	if( !strncmp( Ligne, "[OutputOrder]",13))
	{
		fscanf( FICH, "%s", Ligne);
		if( !strncmp( Ligne, "align", 5))
			ParamExec.OutputOrder = vrai;
		else
		if( !strncmp( Ligne, "input", 5))
			ParamExec.OutputOrder = faux;
		else
		{
			TraiteErr( 34, Ligne );
			return faux;
		}
	}
	else

	/* --- choix de gif ---*/
	if( !strncmp( Ligne, "[GifImage]",10))
	{
		fscanf(FICH,"%s",Ligne);
		if ( !strncmp( Ligne, "yes", 3))
			ParamExec.OutputGif = vrai;
		else ParamExec.OutputGif = faux;
	}
	else

	/* --- format de sortie de clustering ---*/
	if ( !strncmp( Ligne, "[ClusteringOutputFormat]",24))
	{
		fscanf(FICH,"%s",Ligne);
		if (!strncmp( Ligne, "list",4)) ParamExec.FormatArbre= faux;
		else
		if (!strncmp( Ligne, "draw",4)) ParamExec.FormatArbre= vrai;
		else
		{
			TraiteErr (34, Ligne);
			return faux;
		}
	}
	else

	/*--- niveaux de consensus ---*/
	if( !strncmp(Ligne, "[ConsensusLevel]",16))
	{
		fscanf (FICH, "%hd %hd", &ParamExec.ConsH, &ParamExec.ConsL);
	}
	else

	/* style des sorties gif et doc */
	if( !strncmp( Ligne, "[OutputStyle]", 13 ) )
	{
		fscanf(FICH,"%s",Ligne);
		if (!strncmp(Ligne,"Normal",6)) ParamExec.Style= Normal;
		else
		if (!strncmp(Ligne,"Case",4))  ParamExec.Style= Case;
		else
		if (!strncmp(Ligne,"Difference",10)) ParamExec.Style= Difference;
		else
		{
			TraiteErr (34, Ligne);
			return faux;
		}
	}
	else
	if( !strncmp( Ligne, "[LineSize]", 10 ) )
	{
		fscanf(FICH,"%s",Ligne);
		ParamExec.LineSize=atoi(Ligne);
		if (ParamExec.LineSize <=0) ParamExec.LineSize = 50;
		if (ParamExec.LineSize < ParamExec.BlockSize)
			ParamExec.BlockSize=ParamExec.LineSize;
	}
	else
	if( !strncmp( Ligne, "[GraduationStep]", 16 ) )
	{
		fscanf(FICH,"%s",Ligne);
		ParamExec.BlockSize = atoi(Ligne);
		if (ParamExec.LineSize <=0) ParamExec.LineSize = 10;
		if (ParamExec.LineSize < ParamExec.BlockSize)
			ParamExec.BlockSize=ParamExec.LineSize;
	}
	else

		/*--- rubrique inconnue ---*/
	{
		TraiteErr( 34, Ligne );
		return faux;
	}

	return vrai;
}


/*------------------------------------------------------------------------*/

static int RecupConfig(void)
{
	FILE* FICH;
	int i;

#ifdef HTML
	if( (FICH= fopen(ChangeExt(ParamExec.NomFichSeq,".Cfg"),"rt")) ==NULL)
#else
	if( (FICH= fopen("ma.cfg","rt")) ==NULL)
#endif
	{
		TraiteErr( 28, "ma.cfg" );
		return faux;
	}

	MessageAction(hcActRconfig, NULL);

	/* recuperation des parametres du fichier de configuration */
	while (((i=RecupParametre(FICH ))!=0) && (i!=EOF));
	fclose( FICH );

	return (i==EOF);
}
/*---------------------------------------------------------------------------*/
#ifndef _madomain
static int TestLectFich( char* NomFich )
/* retourne vrai si NomFich est accessible en lecture, faux sinon */
{
	FILE* FICH;

	if( (FICH= fopen(NomFich,"r") )!=NULL)
	{
		fclose(FICH);
		return vrai;
	}
	return faux;
}

#	ifndef _Windows
/*MasqRech avec system ne fonctionne pas sous Windows*/
static int MasqRech(char* Masque )
{
	char* Redirection= " > ";
	t_CheminFichier Cmde;

	Cmde[0]= '\0';

	strcat( Cmde, CmdeListeFich );
	strcat( Cmde, Masque );
	strcat( Cmde, Redirection );
	strcat( Cmde, FICH_SELECT );

	if( system(Cmde) )
	{
		TraiteErr( 28, Masque );
		return faux;
	}
	return vrai;
}

/*-----------------------------------------------------------------------*/
#	else
#include <dir.h>
static int MasqRech(char* Masque )
{
	struct ffblk ffblk;
	int done;
	FILE* f;
	if ((f=fopen(FICH_SELECT,"wt"))==NULL)
	{
		TraiteErr( 29, FICH_SELECT );
		return faux;
	}
	done = findfirst(Masque,&ffblk,0);
	while (!done)
	{
		fprintf(f,"%s\n", ffblk.ff_name);
		done = findnext(&ffblk);
	}
	fclose(f);

	return vrai;
}
#	endif
#endif
/* analyse du nom de fichier de sequences */

static int AnalSeqFileName (char* Name)
{
	t_FormatFichSeq def;
	ParamExec.FICHIERS = faux;
#ifndef _madomain
	/* si fichier de noms de fichiers de sequences */
	if( Name[0]=='@' )
	{
		strcpy( ParamExec.NomFichSeq, Name+1 );
		ParamExec.FICHIERS= vrai;
		return vrai;
	}

		/* si c'est un masque pour noms de fichiers */
	if( strchr(Name,'*') || strchr(Name,'?') )
	{
		strcpy( ParamExec.NomFichSeq, FICH_SELECT );
		ParamExec.FICHIERS = vrai;
		return MasqRech(Name );
	}
#endif
		/* sinon c'est un nom de fichier */
	strcpy( ParamExec.NomFichSeq, Name );
	/* Si on part du fichier .mul ou les sequences sont deja alignees */
	if (!strncmp (ExtensionSeule(Name),".mul",4))
	{
		def = ParamExec.FormatSortie;
		ParamExec.RECUP = RecupConfig();
		ParamExec.FormatEntree = MUL;
		ParamExec.FormatSortie = def;
		ParamExec.SeqDejaAligne=1;
	}

	return vrai;
}


/*-------------------------------------------------------------------------*/
#ifndef _madomain
static void SaisieNomFich (char *question, char *answer)
{
	int Correct= faux;

	do
	{
		printf("\n%s: ",question);
		gets( answer );
		if ( !TestLectFich(answer) )
			printf("File not found, or read-protected file !");
		else
			Correct= vrai;
	} while( !Correct );
}

void PEModeInteract(void)
{
	char Format[256];
	char OuiNon[256];
	t_CheminFichier NomFich;
	float lu;

	int Correct;

	/* saisie du nom du fichier de sequences */
	Correct= faux;
	do
	{
		printf("\nSequence file:  ");
		gets( NomFich );
		if (!AnalSeqFileName (NomFich) || !TestLectFich(ParamExec.NomFichSeq) )
			printf("File not found, or read-protected file !");
		else
			Correct= vrai;
	} while( !Correct );

	/* saisie du format d'entree des sequences */
	do
	{
		Correct= vrai;
		printf("\nInput format (gcg, mul, embl, genbank, auto): (def = auto)  ");
		gets(Format );
		if (!*Format) strcpy (Format,"auto");
		if( !strcmp( Format, "gcg" ) )
			ParamExec.FormatEntree= GCG;
		else
		if( !strcmp( Format, "mul" ) )
			ParamExec.FormatEntree= MUL;
		else
		if( !strcmp( Format, "embl" ) )
			ParamExec.FormatEntree= EMBL;
		else
		if( !strcmp( Format, "genbank" ) )
			ParamExec.FormatEntree= GENBANK;
		else
		if( !strcmp( Format, "auto" ) )
			ParamExec.FormatEntree= NON_DEFINI;
		else
		{
			printf("Unknown file format!    ");
			Correct= faux;
		}
	}while( !Correct );

	/* autres parametres d'entree */
	printf("\nOther Input parameters (y/n) ?: (def=n) ");
	gets( OuiNon );
	if( tolower(OuiNon[0]) == 'y' )
	{
	/* choix des portions de sequences */
	printf("\nSelect parts of sequences (y/n) ?: (def=n)  ");
	gets( OuiNon );
	if( tolower(OuiNon[0]) == 'y' )
		ParamExec.SELECTPART= vrai;

	/* saisie du nom du fichier contenant l'arbre */
	printf("\nUse a cluster or a score file (y/n) ?: (def = n) ");
	gets( OuiNon );
	if( tolower(OuiNon[0]) == 'y' ) {
		SaisieNomFich ("\tCluster file (.sco for a score file)",ParamExec.NomFichArbre);
		if (strstr(ParamExec.NomFichArbre,".sco")) ParamExec.LIREARBRE= 2;
		}
	}

	/* saisie du nom de la table de comparaison des symboles */
	printf("\nSymbol comparison table: (def = %s) ",DEF_TAB);
	gets (ParamExec.NomFichCoeff);
	if (!*ParamExec.NomFichCoeff) strcpy(ParamExec.NomFichCoeff,DEF_TAB);

	/* saisie de la valeur de gap */
	printf("\nGap value ? (def = value in Symbol comparison file):   ");
	lu = -1;
	gets (Format);
	sscanf(Format,"%f", &lu);
	ParamExec.Gap = lu;
	if (ParamExec.Gap >=0)
	{
	/* 14/06/94 Gap2 is user-defined */
		printf("\nGap value is %g + length x Gap2; Gap2 value ?: (def = 0) ",lu);
		lu =0;
		gets (Format);
		sscanf(Format,"%f", &lu);
		ParamExec.Gap2 = lu;
	}

	/* autres parametres pour l'alignement */
	printf("\nOther Alignment parameters (y/n) ?: (def=n) ");
	gets( OuiNon );
	if( tolower(OuiNon[0]) == 'y' )
	{
	printf("\n Gap penalty at extremities ? \n");
	printf ("0: none, 1: end, 2: beginning, 3: both : (def = 0) ");
	lu = 0;
	gets (Format);
	sscanf(Format,"%f", &lu);
	if (lu < 4) ParamExec.GapExtT = lu;

	/* saisie du nombre d'iterations */
	printf("\nOne iteration only (y/n) ?: (def = n) ");
	gets( OuiNon );
	if( tolower(OuiNon[0]) == 'y' )
		ParamExec.UNEITER= vrai;

	/* saisie de poids */
	printf("\n Weighted sequences ?  (def=y)\n");
	gets(OuiNon);
	if( tolower(OuiNon[0]) == 'n' )
		ParamExec.WEIGHTED= faux;

	/* saisie de la methode de score */
	printf("\nScoring method ? \n");
	printf("0: absolute, 1:percentage, 2:identity, 3:normalised : (def=0) ");
	lu = 0;
	gets(Format);
	sscanf(Format,"%f",&lu);
	if (lu < 4) ParamExec.SCMETHOD= (t_FormatScores)lu;
	}

	/* saisie du format de sortie des sequences */
	do
	{
		Correct= vrai;
		printf("\nOutput format (msf, mul, doc): (def = msf) ");
		gets(Format);
		if( !strcmp( Format, "msf" )
			|| !strcmp( Format, "" ) )
			ParamExec.FormatSortie= MSF;
		else
		if( !strcmp( Format, "mul" ) )
			ParamExec.FormatSortie= MUL;
		else
		if( !strcmp( Format, "doc" ) )
			ParamExec.FormatSortie= DOC;
		else
		{
			printf("Unknown file format!    ");
			Correct= faux;
		}
	}while( !Correct );

	/* autres parametres pour l'alignement */
	printf("\nOther Output parameters (y/n) ?: (def=n) ");
	gets( OuiNon );
	if( tolower(OuiNon[0]) == 'y' )
	{
	/* saisie de l'ordre */
	printf("\nOuput order as (aligned, input)) ?: (def=aligned) ");
	gets(Format);
	ParamExec.OutputOrder=strncmp( Format, "in",2 );

#ifdef _gif
	printf("\nCreation of a GIF image (y/n) ?: (def=%c) ",
		(ParamExec.OutputGif) ? 'y':'n');
	gets( OuiNon );
	switch (tolower(OuiNon[0]))
	{
		case 0: break;
		case 'y': ParamExec.OutputGif = vrai; break;
		case 'n': ParamExec.OutputGif = faux; break;
	}
#endif
	/* saisie du format de sortie pour clustering */
	printf("\nSave cluster as a drawing or as a list (d/l) ?: (def = l) ");
	gets (OuiNon);
	ParamExec.FormatArbre = (tolower(OuiNon[0])=='d');

	/* saisie des niveaux de consensus */
	printf(
	"\n Consensus levels (default = 90 and 50) ? :  ");
	lu = 90;
	gets (Format);
	sscanf(Format,"%f", &lu);
	if ((lu > 0) && (lu <= 100))
	{
		ParamExec.ConsH = lu;
		lu = 50;
		if (*Format)
		{
			gets (Format);
			sscanf(Format,"%f",&lu);
		}
		if ((lu>0)&& (lu<=100)) ParamExec.ConsL = lu;
	}
	}
}

/*-------------------------------------------------------------------------*/
void PESelectPartSeq( char* NomSeq, long *Debut, long *Fin )
/* demande a l'utilisateur les numeros de residus de debut et de fin
	a prendre en compte dans la sequence */
{
	printf("\nCurrent sequence is:  <%s>\n", NomSeq );
	printf("Enter range:    first (1):  ");
	scanf("%d", Debut);
	printf("                last (0=end):  ");
	scanf("%d", Fin);
	if( *Fin==0 )
		*Fin= UINT_MAX;
}



/*-----------------------------------------------------------------------*/

#endif

static int ExtraitOption( char* Arg, char* Option, char* Param )
/* extrait de Arg le caractere d'Option et la chaine de caracteres
	correspondant au parametre.
	Le parametre optionnel doit etre de la forme:
	-<Caractere_d_option>:<Chaine_parametre> */
{
	int i;

	if( Arg[0] != '-' )
		return 0;

	/* extraire l'option */
	*Option= Arg[1];

	/* extraire le parametre */
	for( i=0; (Param[i]= Arg[i+3])!=NULL; i++ );

	return !0;
}

int PEModeLigneCmde(int Argc, char* Argv[] )
{
	t_CheminFichier Parametre;
	int Arg;
	char Option;
	int Fin;
	char Ext[5];

	/*---  Recuperation des noms de fichiers de sequences et d'arbre ---*/
	/* on doit lire un fichier contenant l'arbre */
	strncpy (Ext, strlwr(ExtensionSeule(Argv[Argc-1])),5);
	Fin=faux;
	if(((Fin=!strncmp( Ext,".clu",5))!=0)	|| !strncmp (Ext,".sco",5))
	{
		if( Argc < 3 )
			Usage();
		ParamExec.LIREARBRE= Fin ? 1 : 2;
		strcpy( ParamExec.NomFichArbre, Argv[Argc-1]);

		Argc-= 2;
	}
	else	/* on doit calculer l'arbre */
	{
		ParamExec.LIREARBRE= faux;
		Argc--;
	}

	/*--- analyse du nom de fichier de sequences ---*/
	if (!AnalSeqFileName (Argv[Argc])) Usage();

	/*--- analyse des options de la ligne de commande ---*/
	for( Arg=1; Arg<Argc; Arg++ )
	{
		if( Argv[Arg][0] != '-' )
			Usage();
		else			/* c'est une option */
		{
			ExtraitOption( Argv[Arg], &Option, Parametre );
			ParamExec.RECUP= faux;
			switch( Option)
			{
			case 'r':  /* utiliser les parametres d ma.cfg */
				if (ParamExec.SeqDejaAligne) break;
				if( !RecupConfig() )
					return faux;
				ParamExec.RECUP= vrai;
				break;

			case 'i':  /* format de fichier d'entree */
				if (ParamExec.SeqDejaAligne) break;
				if( !strncmp(Parametre,"gcg",3) )
					ParamExec.FormatEntree= GCG;
				else
				if( !strncmp(Parametre,"mul", 3) )
					ParamExec.FormatEntree= MUL;
				else
				if( !strncmp(Parametre,"embl", 4) )
					ParamExec.FormatEntree= EMBL;
				else
				if( !strncmp(Parametre,"genbank", 7) )
					ParamExec.FormatEntree= GENBANK;
				else
					Usage();
				break;

			case 'c':  /* fichier de coefficients */
				strcpy( ParamExec.NomFichCoeff, Parametre );
				break;

			case 'g':  /* valeur de gap */
				ParamExec.Gap= atof(Parametre);
				break;

			case 'l': /* 14/06/94 Gap2 is user-defined */
				ParamExec.Gap2=atof(Parametre);
				break;

			case '2':/* with a .mul file, 2 blocks to aligned, 2nd one begins with
				the Parametre_th sequence */
				if (!ParamExec.SeqDejaAligne || ((ParamExec.SeqDejaAligne=atoi(Parametre))<2))Usage();
				break;
			case 'o':  /* format de fichier de sauvegarde */
				if (ParamExec.SeqDejaAligne) break;
				if( !strncmp(Parametre,"msf",3) )
					ParamExec.FormatSortie= MSF;
				else
				if( !strncmp(Parametre,"mul", 2) )
					ParamExec.FormatSortie= MUL;
				else
				if( !strncmp(Parametre,"doc",3) )
					ParamExec.FormatSortie= DOC;
				else
					Usage();
				break;

			case 'p':  /* selection des parties de sequences */
				if (ParamExec.SeqDejaAligne) break;
				ParamExec.SELECTPART= vrai;
				break;

			case 'x': /* extremities gap */
				ParamExec.GapExtT= atoi(Parametre);
				break;

			case '1':  /* une seule iteration */
				ParamExec.UNEITER= vrai;
				break;

			case 'u': /* unweighted sequences */
				ParamExec.WEIGHTED = faux;
				break;

			case 's': /* scoring method */
				if( !strncmp(Parametre,"abs",3) )
					ParamExec.SCMETHOD= ABSOLU;
				else
				if( !strncmp(Parametre,"per",3) )
					ParamExec.SCMETHOD= PERCENT;
				else
				if( !strncmp(Parametre,"ide",3) )
					ParamExec.SCMETHOD= IDENTITY;
				else
				if( !strncmp(Parametre,"nor",3) )
					ParamExec.SCMETHOD= NORMALISED;
				else
					Usage();
				break;

			case 'a': /* order as input */
				ParamExec.OutputOrder=faux;
				break;
			case 'A': /* order as aligned */
				ParamExec.OutputOrder=vrai;
				break;

#ifdef _gif
			case 'f': /* image gif */
				ParamExec.OutputGif= vrai;
				break;
#endif
			case 'd': /* drawing */
				ParamExec.FormatArbre = vrai;
				break;

			case 'k': /* consensus */
				sscanf(Parametre,"%hd.%hd",&ParamExec.ConsH,&ParamExec.ConsL);
				break;

			case 'y': /* undocumented option for doc style */
				sscanf(Parametre,"%c_%d_%d",&Option,&ParamExec.LineSize,&ParamExec.BlockSize);
				switch (Option)
				{
					case 'n':
					case 'N': ParamExec.Style=Normal;break;
					case 'c':
					case 'C': ParamExec.Style=Case;break;
					case 'd':
					case 'D': ParamExec.Style=Difference;break;
					default : Usage();
				}
				if (ParamExec.LineSize <=0) ParamExec.LineSize = 50;
				if (ParamExec.BlockSize <=0) ParamExec.BlockSize = 10;
				if (ParamExec.LineSize < ParamExec.BlockSize)
					ParamExec.BlockSize = ParamExec.LineSize;
				break;
			case 'q': /* quiet */
				ParamExec.MUET=vrai;
				break;

			default:
				Usage();
			}
		}
	}

	return vrai;
}

void PEMajParam(t_pParamAlign Param )
/* mise a jour de Param apres lecture de fichier et de la ligne de commande */
{
	short tmp;
	t_score g, *l;
	t_LignCoeff *c;
	int i;

	/* valeur de gap */
	if( ParamExec.Gap >= 0) 	/* si Gap est defini en ligne de commande */
		Param->Gap= ParamExec.Gap;
	else
	{
		if( !Param->Gap )	/* non defini -> valeur par defaut */
			Param->Gap= GAPDFLT;
		ParamExec.Gap= Param->Gap;
	}
	/* 14/06/94 Gap2 is user-defined */
	if (ParamExec.Gap2>=0) Param->Gap2= ParamExec.Gap2;

	if ((g=Param->Gap2)>0)
	{
		ParamExec.Gap2 = Param->Gap2;
		/* 20/09/94  -g is the score of any residue aligned with a gap */
		for (i=Param->PremLettre,c=Param->Coeff+i,l=Param->Coeff[1]+i;
			i<=Param->DerLettre;i++,c++,l++)
		{
			(*c)[1] = -g;
			 *l = -g;
		}
	}
	else Param->Gap2 = 0;

	Param->GapExtT=ParamExec.GapExtT;


	/* niveaux de consensus */
	if (ParamExec.ConsH) Param->ConsHlevel = ParamExec.ConsH;
	if (ParamExec.ConsL) Param->ConsLlevel = ParamExec.ConsL;
	if (Param->ConsHlevel < Param->ConsLlevel)
	{
		tmp = Param->ConsHlevel;
		Param->ConsHlevel = Param->ConsLlevel;
		Param->ConsLlevel = tmp;
	}

	/* nombre d'iterations */
	Param->UNEITER= ParamExec.UNEITER;
	Param->Weighted = ParamExec.WEIGHTED;
	Param->ScMethod= ParamExec.SCMETHOD;
}

#ifndef _madomain
int PEGetFormatEntree( FILE* FICH)
/* determine le format du fichier de sequences FICH */
{
	char Mot[256],Mot2[256],s[256];
	int c;
	int Trouve;
	char First;

	if (ParamExec.FormatEntree != NON_DEFINI) return ParamExec.FormatEntree;
	rewind( FICH );

	/* passer les espaces de debut */
	while( (c=fgetc(FICH))!=EOF && isspace(c) );
	ungetc( c, FICH );

	/* analyser le premier mot du fichier */
	First=' ';
	fscanf( FICH, "%s", Mot );

	if (Mot[0]=='>')  ParamExec.FormatEntree=MUL;
	else
	{
		if( !strncmp(Mot,"LOCUS",5) )	First='g';
		else if ( !strncmp( Mot, "ID", 2 ) ) First='e';

		/* rechercher la sequence ".." de GCG ou le 2e mot de Genbank et Embl*/
		Trouve=faux;
		while (((fgets(s,256,FICH))!=NULL)&&!Trouve)
		{
			if (!strchr(s,'\n')) continue;
			c = strlen(s);
			sscanf(s,"%s %s",Mot,Mot2);
			if ((s[c-2]=='.')&&(s[c-3]=='.')&&(!strncmp(Mot2,"Length:",7)))
			{
				ParamExec.FormatEntree= GCG;
				Trouve= vrai;
			}
			else if ((First=='g')&&(!strncmp(Mot,"ORIGIN",6)))
					ParamExec.FormatEntree= GENBANK;
			else if ((First=='e')&&(!strncmp(Mot,"SQ",2)))
					ParamExec.FormatEntree= EMBL;
		}

		if (ParamExec.FormatEntree == NON_DEFINI)
		{
			/* format non reconnu */
			TraiteErr( 28, ParamExec.NomFichSeq );
			return NON_DEFINI;
		}

	}
	rewind( FICH );

	return ParamExec.FormatEntree;
}

#else
void 	PESaveInputIfMul(void)
{
	char *s;
	if ((ParamExec.SeqDejaAligne>1)&&(ParamExec.FormatSortie==MUL))
	{
		s = ChangeExt(ParamExec.NomFichSeq,".in");
		remove (s);
		rename( ParamExec.NomFichSeq,s );
		strcpy(ParamExec.NomFichSeq,s);
	}
}

char * PEGetDomainOutputFile(void)
{
	if ((ParamExec.SeqDejaAligne!=1)&&(ParamExec.FormatSortie==MUL)) return ParamExec.NomFichSeq;
	else return NULL;/* no domain output file */
}

int PENoOutput (void)
{
	return (ParamExec.SeqDejaAligne==1) && (ParamExec.FormatSortie==MUL);
}

#endif

