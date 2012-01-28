#include <string.h>
#include "ma.h"
char  CopyRight[5][57];

void ChargeCopyRight (void)
{
const char  TextAbout[57][6]
			  = {"MCPMF", "uouu.", "lpbl ", "tyltC", "ariiO",
				 "lispR", "ighlP", "nheeE", " td T", "v  s,", "eIre ", "r.eq1",
				 "sNsu9", "i.ee8", "oRan8", "n.rc,", " Ace ", "5.h N", ".  au",
				 "0Fulc", ".rsil", "0aig.", "\0nnn ", " cgmA", " e ec", "  tni",
				 " 1htd", " 9i s", " 8sw ", " 9 iR", " ,ste", "  ohs", " 1f .",
				 " 9th,", " 9wi ", " 1ae1", " ,rr6", "  ea ", " 1 r(", " 9sc2",
				 " 9hh2", " 4oi)", " ,uc,", "  la ", " 1dl1", " 9  0", " 9cc8",
				 " 6il8", " \0tu1", "  es-", "  \0t1", "   e0", "   r8", "   i9",
				 "   n0", "   g\0", "\0\0\0\0\0"};

const char version[6] = "5.4.1";
int i,j;

 for (i= 0;i<5;i++)
	for (j=0;j<57;j++)
		CopyRight[i][j] = TextAbout[j][i];
 strncpy (CopyRight[0]+17, version, 5);
}

int PrintCopyRight (FILE *f)
{
	int i;
	if (f!=NULL) for (i=0;i<5;i++) fprintf (f,"%s\n",CopyRight[i]);
   return 5;
}

