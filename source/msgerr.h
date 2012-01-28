#ifndef __MSGERR_H
#define __MSGERR_H

#define cmOK 10
#define hcActPrint 92
#define hcActPrintDisk 93
#define hcActRead 90
#define hcActWrite 91
#define hcActWconfig 94
#define hcActRconfig 95

extern int Muet;
extern int SortieConsole;

/* print usage of MultAlin and exit with value 1 */
void Usage(void);

/* init the use of error and progression messages */
void InitMessages(char *ProgName);

/* choose between a quiet or verbose mode for messages */
void SetMessageMode (int isQuiet );

/* print an error message */
unsigned int TraiteErr(int NoErr, ...);

/* print a progression message */
void MessageAction(int NoMess,char *P);

/* end the use of error and progression messages */
void DoneMessages(void);


#endif
