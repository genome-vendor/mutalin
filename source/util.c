/*-------------------------------------------------------------------------
fichier UTIL.C :

		Procedures de manipulation d'ensembles
				et de chaines de caracteres
--------------------------------------------------------------------------*/


#include <stdio.h>
#include <string.h>
#include "util.h"
#include "const.h"

  /*-------- General purpose string manipulation --------*/
char *strupr (char *S);

char *strlwr (char *S);

char *EffEspDeb(char *S);
  /*-Return a string with leading white space removed*/


char *EffEspFin(char *S);
  /*-Return a string with trailing white space removed*/


char *EffBlancs(char *S);
  /*-Return a string with leading and trailing white space removed*/


static char NewName[80]; /* static variable for temporary strigs */

/*=========  Traitements sur les chaines de caractere =================*/
#ifndef __BORLANDC__
char *strupr (char *S)
{
	int i;
	for (i=0;S[i]!='\0';i++)
		if (islower(S[i])) S[i]=toupper(S[i]);
	return S;
}

char *strlwr (char *S)
{
	int i;
	for (i=0;S[i]!='\0';i++)
		if (isupper(S[i])) S[i]=tolower(S[i]);
	return S;
}
#endif

char *sdelete (char *s,int pos, int len)
/* Return s after deleting len chars beginning at pos */
{
	int i;

	if ((pos <0)||(len<=0)||((i=strlen(s))<=pos)) return s;
	if (i<pos+len)	s[pos]='\0';
	else memmove (s+pos,s+pos+len,i-pos-len+1);
	return s;
}
/*--------------------------------------------------------------------*/
char *NomSeul (char* NomFich)
{
	char *s;

	ChangeExt(NomFich,"");
	if (s=strrchr(NewName,dirchar)) return s+1;
	else return NewName;
}

char *ChangeExt (char* NomFich, char* Ext)
{
	char *s;

	strcpy (NewName, NomFich);
	if ((s=strrchr(NewName,'.'))&&!strchr(s,dirchar)) strcpy (s,Ext);
	else strcat (NewName, Ext);

	return NewName;
}

void Backup( char* NomFich )
{
	FILE* F;
	char Nom_Bak[256], *s;

	if(( F= fopen(NomFich, "r") )!=NULL)	/* si le fichier existe */
	{
		fclose( F );
		strcpy (Nom_Bak, NomFich);
		if ((s=strchr (NomFich,'.'))!=NULL) {
			if ((int)strlen(s)< 4) strcat (Nom_Bak,"~");
			else Nom_Bak[strlen(Nom_Bak)-1]='~';
		}
		else strcat (Nom_Bak, ".bak");

		unlink (Nom_Bak);
		rename( NomFich, Nom_Bak );
	}
}


