/*------------------------------- AFICHSEQ -----------------------------------*/
/*                                                                            */
/* to print the sequence list with different colors / selection               */
/*                                                                            */
/* an updatable message can be added to wait patiently                        */
/*                                                                            */
/* First call : AffichStart (Message_Title, Sequence_List_Pointer)            */
/*   this initialize the display : all sequences in color 1                   */
/*                                                                            */
/* then : AffichAllume (No, Col) display Sequence[No] with color #Col         */
/*   color 0 to erase the name, color 2,3,... for more colors                 */
/*                                                                            */
/* then : AffichInit (Constant_Text, Var_Text_Length) when you want a message */
/*    to be added.                                                            */
/*        AffichShow (long_var) when you want the No message to be present    */
/*        	with the variable text set to long_var                            */
/*			 AffichHide (No) when you want it to be hidden                       */
/*			 AffichUpDate (Long_Integer) when you want the variable text set to  */
/*				Long_Integer in the No message                                    */
/*                                                                            */
/* Last call : AffichEnd (OK,MessFin,Mess)                                    */
/*    if (MessFin) a final message will be displayed : action successful or   */
/*       action cancelled, following OK value.                                */
/*			If OK, Mess is displayed after "action successful"                   */
/*    the display system is deleted                                           */
/*                                                                            */
/*                                                                            */
/*  Affich.c gives default values for these functions :                       */
/*     display messages on the standard output                                */
/*                                                                            */
/*                                                                            */
/* Example :                                                                  */
/*  AffichStart ("Aligning", Sequence);                                       */
/*	 AffichInit ("Position #",5);                                              */
/*  for (iter=0;iter<10;iter++) {                                             */
/*  for (n=1;n<NbSeq;n++)  {                                                  */
/*    AffichAllume (n,2);                                                     */
/*    for (m=0;m<n;m++) {                                                     */
/*       AffichAllume (m,3);                                                  */
/*       AffichShow ((long)PosMax);                                           */
/*       for (i=0; i<PosMax; i++)                                             */
/*          if (i%10) AffichUpDate (i);                                       */
/*       AffichHide ();                                                       */
/*       AffichAllume (m,1);                                                  */
/*   }                                                                        */
/*   AffichAllume (n,1);                                                      */
/*  }                                                                         */
/*  }                                                                         */
/*	 sprintf (Mess," %d iterations",iter);                                     */
/*  AffichEnd (1,1,Mess);                                                     */
/*                                                                            */
/*----------------------------------------------------------------------------*/
#include <string.h>
#include "msgerr.h"
#include "afichseq.h"

typedef char TTitleStr[40];

typedef struct {
	int Len;
	char* Text;
	char* Erase;
	int Visible;
	} TPatience;

typedef TPatience *PPatience;

typedef struct {
	int Color;
	t_pLigneSeq Sequence;
	} TListeSeq;

typedef TListeSeq *PListeSeq;

static PListeSeq AffichGroup = NULL;
static PPatience Pat = NULL;

void AffichInit(char *AText, int MaxLen)
{
	if (Muet || (Pat) ||
		((Pat = (PPatience)malloc(sizeof(TPatience)))==NULL)) return;
	Pat->Len = MaxLen;
	Pat->Text = (char *) malloc(strlen(AText)+1);
	strcpy(Pat->Text,AText);
	Pat->Erase = (char *) malloc(MaxLen+1);
	memset (Pat->Erase,'\b',MaxLen);
	Pat->Erase[MaxLen] = 0;
	Pat->Visible = 0;
}

void AffichDone(void)
{
	if (!Pat) return;
	if (Pat->Visible) printf("\n");
	free(Pat->Text);
	free(Pat->Erase);
	free(Pat);
	Pat = NULL;
}

void AffichShow(long l)
{
	if (!Pat) return;
	if (!(Pat->Visible = SortieConsole)) return;
	printf ("%s%*s",Pat->Text,Pat->Len," ");
	AffichUpDate (l);
}

void AffichHide(void)
{
	if (Pat) Pat->Visible = 0;
}

void AffichUpDate( long l)
{
	char F[16];
	int len;

	if ((!Pat)||(!Pat->Visible)) return;
	len = Pat->Len;
	sprintf(F,"%*ld",len,l);
	if ((int)strlen(F) > len) memmove (F, F+strlen(F) - len,len+1);
	printf("%s%s",Pat->Erase,F);
}

void AffichStart(char *ATitle, t_pLigneSeq ASequence)
{
	if ((AffichGroup = malloc (sizeof(TListeSeq)))==NULL) return;
	AffichGroup->Color = 1;
	AffichGroup->Sequence = ASequence;
	if (!Muet) printf("%s\n",ATitle);
}


void AffichEnd (int OK,int MessFin, char *Mess)
{
	AffichDone ();
	if (!AffichGroup) return;
	if (!Muet && MessFin)
	{
		if (OK) printf("Action successful %s\n",Mess);
		else printf("Action cancelled\n");
	}
	free (AffichGroup);
	AffichGroup = NULL;
}

void AffichAllume(t_indN No, unsigned int Col)
{
	if (!AffichGroup || !AffichGroup->Sequence) return;
	if (Muet | !SortieConsole) return;
	if (Col == AffichGroup->Color) {
		if (Col != 1) printf(".");
		return;
	}
	switch(Col) {
		case 1 :
			printf("%c%79s%c",13," ",13);break;
		case 3 :
			printf("with %s ",AffichGroup->Sequence[No].NomSeq);break;
		case 2 :
			printf("%s ",AffichGroup->Sequence[No].NomSeq);break;
	}
	AffichGroup->Color = Col;
}

