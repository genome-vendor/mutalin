/*-------------------------------------------------------------------------
	fichier DISK.C
			Procedures generales de traitement sur les fichiers

-------------------------------------------------------------------------*/
#include <string.h>

#include "disk.h"

#include "afichseq.h"
#include "basproc1.h"
#include "commande.h"
#include "msgerr.h"
#include "parametr.h"
#include "util.h"

#ifdef _madomain
#include "madomain.h"
#endif

/*=====    Lecture de fichier de sequences et de coefficients   ==========*/

static int InitConsensus( t_pDescripteurSequence DS )
{
	t_pLigneSeq tmp;

	/* initialisations des valeurs de consensus */
	tmp= DS->Sequence+0;

	strcpy( tmp->NomSeq, "Consensus" );
	tmp->NbResid= 0;
	tmp->PremCol= 1;
	tmp->DerCol= 0;
	tmp->Taille= 0;
	tmp->BuffLen= 0;
	if( (tmp->Seq= (t_pSeq)calloc(64, sizeof(t_Seq)))==NULL )
		return faux;
	tmp->BuffLen= 64;

	return vrai;
}

static int InitLectureFichSeq( t_pDescripteurSequence DS)
{
	/* alloue et initialise le consensus */
	if( (DS->Sequence= (t_pLigneSeq)malloc(sizeof(t_LigneSeq)))==NULL)
	{
		TraiteErr(8);
		return faux;
	}
	InitConsensus( DS );
	return vrai;
}


/*------------------------------------------------------------------------*/
static int LitUnSeulFichier( t_pDescripteurSequence DS, char *Name, int select,
	char* EnsLet)
{
	FILE* FICH;

#ifdef _madomain
	PESaveInputIfMul();
#endif
	if((FICH= fopen(Name,"rt")) ==NULL)
	{
		TraiteErr(28, Name);
		return faux;
	}

#ifdef _madomain
	if (!MulProdomLitSeq (FICH,DS,PEGetDomainOutputFile(),EnsLet))
	{
		fclose(FICH);
		return faux;
	}
#else

	/* selection selon le format */
	switch( PEGetFormatEntree(FICH) )
	{
		case MUL:
			if( ! MulLitSeq( FICH, DS, select , EnsLet ) )
			{
				fclose(FICH);
				return faux;
			}
			break;
		case EMBL:
			if( ! EmblGenbankLitSeq( FICH, DS, EMBL,
						Name, select, EnsLet ) )
			{
				fclose(FICH);
				return faux;
			}
			break;
		case GENBANK:
			if( ! EmblGenbankLitSeq( FICH, DS, GENBANK,
						Name, select, EnsLet ) )
			{
				fclose(FICH);
				return faux;
			}
			break;
		default :
			TraiteErr( 28, Name );
			fclose(FICH);
			return faux;
	}
#endif
	fclose( FICH );
	MajTailleConsens( DS->Sequence, DS->NbSeq, faux );
	return vrai;
}
#ifndef _madomain
/*------------------------------------------------------------------------*/
static int LitPlsFichiers( t_pDescripteurSequence DS, char *Name, int select,
	char* EnsLet )
{
	FILE* FICH;
	char *Nom;
	if( ! InitListeNomFich(Name ) )
		return faux;

	while ((Nom = GetNextFileinList())!=NULL)
	{
		if((FICH= fopen(Nom,"rt"))==NULL )
		{
			TraiteErr(28, Nom);
			return faux;
		}

		/* selection selon le format */
		switch( PEGetFormatEntree(FICH) )
		{
			case MUL:
				if( ! MulLitSeq( FICH, DS, select , EnsLet ) )
				{
					fclose(FICH);
					return faux;
				}
				break;
			case GCG:
				if( !GcgLitSeq( FICH, DS, Nom, select, EnsLet ) )
				{
					fclose(FICH);
					return faux;
				}
				break;

			case EMBL:
				if( ! EmblGenbankLitSeq( FICH, DS, EMBL,
					Name, select, EnsLet ) )
				{
					fclose(FICH);
					return faux;
				}
				break;

			case GENBANK:
				if( ! EmblGenbankLitSeq( FICH, DS, GENBANK,
					Name, select, EnsLet ) )
				{
					fclose(FICH);
					return faux;
				}
				break;
		}
		fclose(FICH);
	}
	DetruitListeNomFich();
	MajTailleConsens( DS->Sequence, DS->NbSeq, faux );

	return vrai;
}
#endif


/*------------------------------------------------------------------------*/

static int LitSeq(t_pDescripteurSequence DS, char* EnsLet)
/* selection du mode de lecture selon le format */
{
	int select = PECanSelect();
	char* FileName = PEGetNomFichSeq();
	/*---- initialiser les allocations de sequences et le consensus ----*/
	if( ! InitLectureFichSeq( DS) )
		return faux;

	MessageAction (hcActRead, FileName);
	AffichInit ("Sequence #",4);
	AffichShow (0);

	/*----- plusieurs fichiers contenant une seule sequence -----
			(Gcg, Embl, GenBank) */
#ifndef _madomain
	if( PEIsList() )
	{
		if( !LitPlsFichiers( DS, FileName, select, EnsLet ) )
			return faux;
	}

	/*--- un seul fichier contenant plusieurs sequences ---
			(MultAlin et eventuellement Embl et GenBank) */
	else
#endif
	{
		if( !LitUnSeulFichier( DS, FileName, select, EnsLet ) )
			return faux;
	}
	AffichDone ();
	MajTailleConsens (DS->Sequence,DS->NbSeq,faux);
	MajWeightConsens (DS->Sequence,DS->NbSeq);
	return vrai;
}

/*------------------------------------------------------------------------*/
int LitSeqCoeff( t_pDescripteurSequence DS, t_pParamAlign Param)
{
	FILE* FICH;
	char EnsLet[128], Nom[80], *l;
	int OK;
	char * FileName = PEGetCoeffFileName();
	memset (EnsLet, 0,128);
	/*----------	Lecture des sequences  ------------------*/
	if( !LitSeq(DS, EnsLet) )
		return faux;

	/*----------  Lecture de la table de coefficients  ------*/
	if((FICH= fopen(FileName,"rt"))==NULL )
	{
		if ((l=getenv("MULTALIN"))!=0)
		{
			strcpy (Nom,l);
			strcat (Nom,FileName);
			FICH= fopen(Nom,"rt");
		}
		if (!FICH)
		{
			TraiteErr(28, FileName);
			return faux;
		}
	}

	strcpy(Param->NomCoeff,NomSeul(FileName));
	MessageAction (hcActRead, FileName);
	OK = LitFichCoeff (FICH,  FileName,Param);
	fclose(FICH);
	if (!OK) return faux;

	if( Param->DerLettre > 0 )
		VerifiSymb(DS,Param, EnsLet);
	return vrai;
}

/*==================  Sauvegarde des sequences  ==========================*/
int EcritSeq( t_pDescripteurSequence DS, t_pParamAlign Param )
{
	FILE* F;
	char NomFich[256];
	t_ExtensionFichier Extension;
#if defined HTML && !defined _madomain
	int pass;
#endif
	t_FormatFichSeq Format = (t_FormatFichSeq)PEGetOutputFormat();
	char* NomFichOri = PEGetNomFichSeq();
#ifdef _gif
	/* creer l'image gif correspondante */
	if (PEGifOutput()) CreerImage(NomFichOri ,DS, Param,0,PEGetSeqDejaAligne());
#endif
#ifdef _madomain
	if (PENoOutput()) return vrai; /* pas d'ecriture si pas de nouvel alignement */
#elif defined HTML
pass = 3;
Format = MUL;
while (pass)
{
#endif
	switch( Format )
	{
		case MSF:
			strcpy( Extension, ".msf" );
			break;
		case MUL:
			strcpy( Extension, ".mul" );
			break;
		case DOC:
			strcpy( Extension, ".doc" );
			break;
	}
	/* meme nom avec extension selon format */
	strcpy( NomFich, ChangeExt( NomFichOri,Extension) );
#ifdef _madomain
 if (Format==MUL) F= fopen(NomFich,"at");
 else
#endif
 {
	/* si le fichier existe deja, le renommer en .bak */
	Backup( NomFich );
	F= fopen(NomFich,"wt");
 }
	if( !F )
	{
		TraiteErr(29,NomFich);
		return faux;
	}
	MessageAction(hcActWrite, NomFich );
	switch( Format )
	{
		case MSF:
			MsfEcritSeq( F, DS, Param,0,0 );
			break;
		case DOC:
			MsfEcritSeq(F,DS,Param,1,PEGetSeqDejaAligne());
			break;
		case MUL:
#ifdef _madomain
			MulProdomEcritSeq (F,DS);
#else
			MulEcritSeq( F,DS );
			fclose(F);
			strcpy( NomFich, ChangeExt( NomFichOri,".con") );
			if( (F= fopen(NomFich,"wt")) !=NULL)
			{
				MessageAction(hcActWrite, NomFich );
				Ecrit1Seq (F,DS->Sequence,DS->NbSeq);
			}
#endif
			break;
	}
	fclose(F);
#if (defined(HTML)&&!defined(_madomain))
	pass--;
	if (pass>1)Format = MSF;
	else Format = DOC;
}
#endif


	return vrai;
}



