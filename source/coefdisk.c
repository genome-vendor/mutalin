/*------------------------------------------------------------------------
	Fichier COEFDISK.C

			Lecture des tables de comparaison de symoboles
------------------------------------------------------------------------*/
#include <ctype.h>
#include <string.h>
#include "disk.h"
#include "msgerr.h"
#include "parametr.h"





/*------------------------------------------------------------------------*/


static int LitMulFichCoeff( FILE* FICH, char* NomFich, t_pParamAlign Param )
{
	register int i,j;
	register int c;
	int cc;
	int DebutDesCoeff,OK;
	float lu;
	t_score s,ss;

	/*------------------------ passer les commentaires */
	DebutDesCoeff= faux;
	while( !DebutDesCoeff && !feof(FICH) )
	{
		/* passer les blancs */
		while( (c= fgetc(FICH))!=EOF && isspace(c) );

		/* passer les commentaires */
		if( c=='>')
			while( (c=fgetc(FICH))!='\n' && c!=EOF  ); /* passer a la ligne */

		else
			DebutDesCoeff= vrai;	/* on est au debut de la sequence */
	}

	/*------------------- lire la ligne des symboles */
	if (c==' ') while ((c= fgetc(FICH))==' ');
	for(i=Param->PremLettre-1; (i<DERLETTREMAX) && c!= EOF && c!='\n';)
	{
		if( Param->NumSymb[c]!= FreeSymb )
		{
			TraiteErr (12,NomFich);
			InitParam (Param);
			return faux;
		}
		i++;
		Param->Symb[i]= c;
		Param->NumSymb[c]=i;
		while ((c= fgetc(FICH))==' ');
	}

	Param->DerLettre= i;

	/*-------------------- lire les coefficients */
	for( i=Param->PremLettre,OK=1,ss=0;(OK==1) &&(i<=Param->DerLettre); i++ )
	{
		for( j=i;(OK==1) && (j<=Param->DerLettre); j++ )
		{
			while( ((OK=fscanf(FICH,"%g",&lu))!=1 ) && (OK!=EOF)) c=fgetc(FICH);
			Param->Coeff[i][j]= s=(t_score)lu;
			if (s > ss) ss=s;
		}
		while( (c=fgetc(FICH))!='\n' && c!=EOF  );
	}
	Param->CoeffMax = ss;

	if (i<=Param->DerLettre)
	{
		TraiteErr (12,NomFich);
		InitParam (Param);
		return faux;
	}
	/* rendre symetrique la matrice de coefficients */
	for( i=Param->PremLettre; i<Param->DerLettre; i++ )
		for( j=i+1; j<=Param->DerLettre; j++ )
			Param->Coeff[j][i]= Param->Coeff[i][j];

	while( (c=fgetc(FICH))!='\n' && c!=EOF  ); /* passer a la ligne */
	if (feof(FICH)) return vrai;

	/*------------------- lire la valeur de gap */
	if( fscanf(FICH,"%f",&lu)== 1 )
	{
		Param->Gap= lu;
		if( fscanf(FICH,"%f",&lu)== 1 )
			Param->Gap2= lu;
	}

	while( (c=fgetc(FICH))!='\n' && c!=EOF  ); /* passer a la ligne */



	/*------------------ initialisation des clusters de symboles */
	cc = Param->DerLettre;
	while( !feof(FICH) )
	{
		/* lire le symbole representant les symboles homologues dans le consensus */
		/* passer les blancs */
		while( ((c=fgetc(FICH))!=EOF) && isspace(c) );
		if (c==EOF) break;
		i=Param->NumSymb[c];
		if ((i==FreeSymb)&&(cc<DERLETTREMAX))
		{
			cc++;
			Param->NumSymb[c]=cc;
			i=cc;
			Param->Symb[i]=c;
		}
		if ((i>0)&&(i<=DERLETTREMAX))
			while ((c!=EOF) && (c!='\n'))
			{
				while ((c=fgetc(FICH))!=EOF && (c==' ') && !(c=='\n'));
				if ((c>0) && ((j = Param->NumSymb[c])>0)&&(j<=DERLETTREMAX))
					Param->Hom[j] = i;
			}
		else
			while( (c=fgetc(FICH))!='\n' && c!=EOF  ); /* passer a la ligne */
	}

	return vrai;
}
/*---------------------------------------------------------------------*/
static int OldLitGcgFichCoeff( FILE* FICH, char* NomFich, t_pParamAlign Param )
{
	register int i,j;
	register int c;
	float Flottant;
	int DebutDesCoeffs;
	int FinDesCommentaires;
	long PosDebLigne;
	t_score s, ss;

	/*------------------------ passer les commentaires */
	for( c=' ', FinDesCommentaires=faux; (c!=EOF) && !FinDesCommentaires;
				c=fgetc(FICH)  )
	{
		if( c=='.' )
		{
			if( (c=fgetc(FICH))=='.' )
				FinDesCommentaires= vrai;
			else
				ungetc( c, FICH );
		}
		/* memoriser les debuts de ligne */
		if( c=='\n' )
			PosDebLigne= ftell( FICH );

	}
	fseek( FICH, PosDebLigne, 0 );

	/*------------------- lire la ligne des symboles */
	for(c=' ',i=Param->PremLettre,DebutDesCoeffs=faux; (c!= EOF) && !DebutDesCoeffs;
					c= fgetc(FICH) )
	{
		if( c=='.' )
		{
			if( (c= fgetc(FICH))=='.')
				DebutDesCoeffs= vrai;
			else
				ungetc( c, FICH );
		}
		else if( !isspace(c) )
		{
			Param->Symb[i]= c;
			Param->NumSymb[c]=i;
			i++;
		}

	}
	Param->DerLettre= i-1;

	if(feof(FICH)||(i>=DERLETTREMAX))
	{
		TraiteErr(12, NomFich );
		InitParam (Param);
		return faux;
	}


	/*-------------------- lire les coefficients */
	for( i=Param->PremLettre,ss=0; i<=Param->DerLettre; i++ )
	{
		for( j=i; j<=Param->DerLettre; j++ )
			if( (fscanf(FICH,"%f",&Flottant))==EOF )
			{
				TraiteErr(12, NomFich );
				InitParam (Param);
				return faux;
			}
			else
			{
				Param->Coeff[i][j]= s = float_to_score(Flottant*10);
				if (s>ss) ss=s;
			}
		fscanf(FICH,"%*s \n");
	}
	Param->CoeffMax=ss;

	/*---------- rendre symetrique la matrice de coefficients */
	for( i=Param->PremLettre; i<Param->DerLettre; i++ )
		for( j=i+1; j<=Param->DerLettre; j++ )
			Param->Coeff[j][i]= Param->Coeff[i][j];

	return vrai;
}

static int LitGcgFichCoeff( FILE* FICH, char* NomFich, t_pParamAlign Param )
{
	char Line[256],*s;
	int i, j, OK, lu;
	t_score ss, sc;
	s=fgets(Line,256,FICH);
	if (s && strncmp(Line,"!!",2))
	{
		rewind(FICH);
		return OldLitGcgFichCoeff(FICH,NomFich,Param);
	}
	/*------------------------ passer les commentaires */
	while (s && !strstr(Line,"..\n")) s=fgets(Line,256,FICH);

	while (((s=fgets(Line,256,FICH))!=NULL) && (Line[0]=='\n'));
	if (s && (Line[0]=='{'))
	/*------------------- lire la valeur de gap */
	{
		while (fgets(Line,256,FICH) && (Line[0]!='}')) if (Line[0]!='!')
		{
			if (strstr(Line,"GAP_CREATE")) sscanf(Line,"%*s %d",&Param->Gap);
			else if (strstr(Line,"GAP_EXTEND")) sscanf(Line,"%*s %d",&Param->Gap2);
		}
		while (((s=fgets(Line,256,FICH))!=NULL) && (Line[0]=='\n'));
	}
	if (!s)
	{
		TraiteErr (12,NomFich);
		InitParam (Param);
		return faux;
	}

	/*------------------- lire la ligne des symboles */
	for(s=Line,i=Param->PremLettre; (i<DERLETTREMAX) && *s && *s!='\n';s++,i++)
	{
		while (isspace(*s)) s++;
		if (!*s) break;
		if( Param->NumSymb[*s]!= FreeSymb )
		{
			TraiteErr (12,NomFich);
			InitParam (Param);
			return faux;
		}
		Param->Symb[i]= *s;
		Param->NumSymb[*s]=i;
	}

	Param->DerLettre= i-1;
	/*-------------------- lire les coefficients */
	for( i=Param->PremLettre,OK=1,ss=0;OK &&(i<=Param->DerLettre);)
	{
		if (!fgets(Line,256,FICH)) break;
		s = strtok(Line," ");
		if (s && (Param->Symb[i]==s[0]))
		{
			for( j=Param->PremLettre;OK && (j<=i); j++ )
			{
				s = strtok(NULL," ");
				OK = s && (sscanf(s,"%d",&lu)==1);
				Param->Coeff[i][j]= sc=(t_score)lu;
				if (sc > ss) ss=sc;
			}
			i++;
		}
	}
	Param->CoeffMax = ss;

	if (i<=Param->DerLettre)
	{
		TraiteErr (12,NomFich);
		InitParam (Param);
		return faux;
	}
	/* rendre symetrique la matrice de coefficients */
	for( i=Param->PremLettre; i<Param->DerLettre; i++ )
	{
/*		fprintf (stderr,"\n%c ",Param->Symb[i]);*/
		for( j=Param->PremLettre; j<i; j++ )
		{
			Param->Coeff[j][i]= Param->Coeff[i][j];
/*			fprintf(stderr,"%4d",Param->Coeff[i][j]);*/
		}
	}
	return vrai;
}

int LitFichCoeff (FILE *FICH, char* FileName, t_pParamAlign Param)
{
	if (strstr (FileName, ".cmp") != NULL)
		return LitGcgFichCoeff( FICH, FileName, Param ) ;
	else
		return LitMulFichCoeff( FICH, FileName, Param ) ;
}
/*-----------------------------------------------------------------------*/




