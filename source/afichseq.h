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

#ifndef __AFICHSEQ_H
#define __AFICHSEQ_H
#include "deftypes.h"

void AffichStart(char *ATitle, t_pLigneSeq ASequence);
void AffichEnd (int OK,int MessFin,char *Mess);
void AffichAllume(t_indN No,unsigned int Col);

void AffichInit(char *AText, int MaxLen);
void AffichShow(long l);
void AffichHide(void);
void AffichUpDate(long l);
void AffichDone(void);

#endif
