/*-------------------------------------------------------------------------
	DEFTYPES.H
				Definition des types utilises pour MA
-------------------------------------------------------------------------*/
/*
2 juin 1994 : Gap can be fonction of gap length
	Gap penalty = Gap + Gap2 * length

14 juin 1994 : Gap2 is a user defined parameter

8 july 1994 : modification for domainer program (conditional compilation)

18 july 1994 : GapExtT defines the penalty for the extremities : two bits, the
	higher for the beginning, the lower for the end. Bits are set when a penalty
	must be counted for an extremity gap.
*/

#ifndef __DEFTYPES_H

#define __DEFTYPES_H
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "const.h"
/* type de fichier de sequences */
typedef enum { NON_DEFINI=-1,GCG, MUL, EMBL, GENBANK, MSF, DOC } t_FormatFichSeq;
/* type de methodes de score pour l'arbre guide */
typedef enum { ABSOLU, PERCENT, IDENTITY, NORMALISED } t_FormatScores;

/*============== TYPES USUELS ======================================*/
typedef long t_score;/*/
typedef float t_score;*/
#define SCO_MIN LONG_MIN
#define SCO_MAX LONG_MAX
#define float_to_score(F) (t_score)((F) + (((F)>=0) ? 0.5 : -0.5))

typedef char t_StringNom[LONGNOM];
typedef char t_ExtensionFichier[5];
typedef char t_CheminFichier[256]; 	/* nom complet d'un fichier */





/*=============  PARAMETRES D'ALIGNEMENT  ===============================*/
/* matrice des coefficients */
typedef t_score t_LignCoeff[DERLETTREMAX+1];
typedef t_LignCoeff t_TabCoeff[DERLETTREMAX+1];

/* tableau des symboles de substitution */
typedef signed char byte;
typedef byte t_TabHom[DERLETTREMAX+1];
typedef byte* t_pTabHom;

/* tableaux pour la conversion
	NumeroSymbole (0 a 30) <-> CaractereSymbole (0 a 127: code ASCII)*/
typedef byte t_TabNumSymb[128];		/* correspondance Numero <-> Symbole */
typedef byte t_TabSymb[DERLETTREMAX+1]; /* correspondance Symbole <-> Numero */

/* parametres d'alignement */
typedef struct
{
	t_TabCoeff Coeff;
	t_score Gap;
	t_score Gap2; /*2/6/94 gap=func(length)*/
	char NomCoeff[40];
	t_score CoeffMax;
	t_TabNumSymb NumSymb;
	t_TabSymb Symb;
	t_TabHom Hom;			/* homologies */
	byte PremLettre;
	byte DerLettre;
	int UNEITER;
	int Weighted;
	unsigned GapExtT:2;
	t_FormatScores ScMethod;
	short ConsHlevel;
	short ConsLlevel;
} t_ParamAlign;
typedef t_ParamAlign* t_pParamAlign;


/*======================  SEQUENCES  ====================================*/
/* residus */
typedef  union
{
	byte Car;     	/* 0 a 127 */
	byte Symb;     /* 0 a DERLETTREMAX */
} t_Seq;

typedef t_Seq* t_pSeq;

typedef long int t_Weight;/*/
typedef float t_Weight;*/
/* description d'une sequence */
typedef struct
{
	t_StringNom NomSeq;
#ifdef _domainer
	int name;
	int shift;
#endif
	t_indL NbResid;
	t_indL PremCol;
	t_indL DerCol;
	t_indL Taille;
	t_indL BuffLen;
	t_pSeq Seq;
	t_Weight Weight;
} t_LigneSeq;
typedef t_LigneSeq* t_pLigneSeq;
typedef t_pLigneSeq t_pTabSeq;

/* ensemble des sequences */
typedef struct
{
	t_pTabSeq Sequence;
	t_indN* NoSeq;		/* ordre des sequences selon l'arbre */
	t_indN NbSeq;
} t_DescripteurSequence;
typedef t_DescripteurSequence* t_pDescripteurSequence;


/*===============  CLASSIFICATION HIERARCHIQUE  ==========================*/
/* branche */
typedef struct
{
	t_indN SeqSuiv;
	t_indN Pere;
} t_Branche;

/* arbre */
typedef t_Branche* t_ClassHier;

typedef struct
{
	t_indN Noeud;
	t_indN Col;
} t_Dessin;

typedef t_Dessin* t_DessinArb;

/* tableau de 3 arbres */
typedef t_ClassHier t_TabClassHier[3];
typedef t_ClassHier* t_pTabClassHier;

typedef struct
{
	t_indN G1;
	t_indN G2;
	t_indN *NoGroupe;				/* tableau de numero de la premiere sequence de chaque groupe */
	t_indN *Groupe;   			/* tableau de nombres de sequences par groupe */
	t_indL M;                        /* taille de Groupe1 */
	t_indL M1;                       /* M-1 */
	t_indL N;                        /* taille de Groupe2 */
	t_indL N1;                       /* N-1 */
	t_indN Nb1;
	t_indN Nb2;
} t_DescripteurGroupe;

typedef t_DescripteurGroupe* t_pDescripteurGroupe;


/*==========  MATRICES ET TABLEAUX POUR L'ALIGNEMENT ===============*/

typedef struct
{
	t_score* Mu;
	t_score* X;
	t_indL* PMu;
} t_GrdTab;
typedef t_GrdTab* t_pGrdTab;

/* standard prototypes that are not always in the same standard header files */
void *memset(void *s, int c, size_t n);
void *memcpy(void *dest, const void *src, size_t n);

#endif	/* __DEFTYPES_H */




