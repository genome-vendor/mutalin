/*-------------------------------------------------------------------------
	Fichier MA.C
				Programme principal
--------------------------------------------------------------------------*/
/*
14 juin 1994 : Gap penalty is a function of gap length : Gap + length * Gap2
	Gap2 is user-defined (0 by default)
*/
#include <ctype.h>
#include "ma.h"

#include "basproc1.h"
#include "commande.h"
#include "disk.h"
#include "msgerr.h"
#include "parametr.h"
#ifdef _madomain
#	include "madomain.h"
#endif

#ifdef __BORLANDC__
#	include <dos.h>
	/* augmentation de la taille de pile a 32 ko ( 4 ko par defaut ) */
	extern unsigned _stklen= 0x8000;
#endif

static int ChercheConsensus(t_pDescripteurSequence DS, t_pParamAlign Param)
/* recherche la sequence consensus en fonction des lettres et homologues
	les plus frequents dans les sequences */
{
	int C;
	t_Weight Tab[DERLETTREMAX+1]; /*frequence des lettres homologues*/
	t_Weight Tab2[DERLETTREMAX+1];		/* frequence des lettres */

	int LettreMax, LettreMax2, Lettre;
	t_Weight Max, Max2, *t, NivH, NivL;
	t_indL ColCour, Long;
	t_indN Seq1;
	t_pSeq Consens;
	t_pLigneSeq tmp;
	long WeightT=DS->Sequence->Weight;

	MajTailleConsens( DS->Sequence, DS->NbSeq, vrai ) ;/*1/7/94 add vrai */
	if( DS->Sequence[0].BuffLen== 0 )
		return faux;
  /*	MessageAction( ACT_CHERCHE_CONSENS );*/
	tmp = DS->Sequence;
	Long = tmp->DerCol;
	NivH = (Param->ConsHlevel * WeightT +99)/100;
	NivL = (Param->ConsLlevel * WeightT +99)/100;

	/* recherche par colonne */
	for( ColCour=1, Consens= DS->Sequence[0].Seq+1; ColCour<=Long; ColCour++, Consens++)
	{
		memset(Tab,0,(DERLETTREMAX+1)*sizeof(t_Weight));
		memset(Tab2,0,(DERLETTREMAX+1)*sizeof(t_Weight));
		for( Seq1=0,tmp=DS->Sequence+1; Seq1<DS->NbSeq; Seq1++,tmp++)
		{
			C= tmp->Seq[ColCour].Car;
			if( (ColCour>tmp->DerCol) || (C==Insert) || (C==' '))
				Lettre= 0;
			else
				Lettre= Param->NumSymb[C];
			Tab2[Lettre]+=tmp->Weight;
			Tab[Param->Hom[Lettre]]+=tmp->Weight;
		}

		/* recherche de la classe d'homologie la + frequente */
		Max= 0;
		LettreMax= 0;
		for( Lettre=0, t= Tab; Lettre<=DERLETTREMAX; Lettre++,t++ )
			if( *t > Max )
			{
				LettreMax= Lettre;
				Max= *t;
			}
		LettreMax = Param->Symb[LettreMax];

		/* recherche de la lettre la + frequente */
		Max2= 0;
		LettreMax2= 0;
		for( Lettre=0, t= Tab2; Lettre<=Param->DerLettre; Lettre++, t++ )
			if( *t > Max2 )
			{
				LettreMax2= Lettre;
				Max2= *t;
			}
		LettreMax2 = Param->Symb[LettreMax2];

		if( Max2 >= NivH)
		/* lettre presente  dans + de 90% des seq */
			Consens->Symb= LettreMax2;
		else
			if( Max >= NivH)
			/* classe d'homol presente dans + de 90% des seq */
				Consens->Symb= LettreMax;
			else
				if(( Max2 >= NivL ) && isupper(LettreMax2))
				/* lettre presente dans + de 50% des seq */
					Consens->Symb= tolower(LettreMax2);
				else
					/* pas de consensus */
						Consens->Symb= '.';
	}
	DS->Sequence->Seq[DS->Sequence->Taille+1].Car=0;
	return vrai;
} 	/* de ChercheConsensus */

static int Run(void )
/* true main function: return 1 if success and 0 if error */
{
	t_DescripteurSequence DS;
	t_ParamAlign Param;
	int SeqDejaAligne;

	/* init the sequence structure */
	InitSequenceDescription(&DS);
	InitParam(&Param);   	/* in PARAMETR.C */

	/*-- lecture des sequences ---------------------------------------*/
	if( !LitSeqCoeff(&DS, &Param) )
		return 0;

	if(Param.DerLettre < 1)/* empty alphabet */
		return 0;

	if(DS.NbSeq < 2)/* not enough sequences */
		return 0;

	PEMajParam(&Param );

	/*-- on sauvegarde les parametres d'execution dans ma.cfg ---------*/
	PESauveConfig();

	SeqDejaAligne = PEUpDateMul(DS.NbSeq);
	if (SeqDejaAligne==1)
		printf("Sequences are already aligned\n");

	if (!MultipleAlignment(&DS,&Param,SeqDejaAligne))
	{
		TraiteErr(8);
		return 0;
	}
#	ifdef _madomain
	if (SeqDejaAligne != 1)
		EcritProdomConsensus (PEGetNomFichSeq(),&DS,&Param);
#	endif

	if( !ChercheConsensus( &DS, &Param ) )
		TraiteErr( 23);


	/* sauvegarde des sequences */
	EcritSeq(&DS, &Param);

	LibereDescripteurSequence( &DS );
	return vrai;
}

/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
int main( int argc, char* argv[] )
{
	int OK;

	ChargeCopyRight ();
	PrintCopyRight(stdout);
	PEInit();
	InitMessages (argv[0]);

#ifdef _madomain
	if ((argc < 3)&&(argv[argc-1][0]=='-'))
	{
		printf(" madomain needs more argument");
		exit(1);
	}
#else
	if( argc==2 && argv[1][0]=='-' )
		Usage();

	/* on utilise le mode interactif */
	if( argc==1 )
		PEModeInteract();
#endif

	/* on utilise le mode ligne de commande */
	else
		if( !PEModeLigneCmde(argc, argv ) )
			exit(1);

	SetMessageMode (PEIsQuiet());
	InitBreakControl ();
	OK = Run();
	DoneBreakControl ();
	DoneMessages ();

	return OK ? 0 : -1;
}
