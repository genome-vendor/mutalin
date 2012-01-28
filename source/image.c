/*------------------------------------------------------------------------
	Fichier IMAGE.C
	Enregistrement au format gif 

	Requiert un fichier de Configuration NomFich.cfg
	Cree un fichier Image NomFich.gif
------------------------------------------------------------------------*/
#include <string.h>
#include <ctype.h>

#include "util.h"
#include "commande.h"
#include "disk.h"
#include "ma.h"
#ifdef _gif
/* Modif 14/05/96 Tim Downs */

#include "gd.h"
#include "gdfontg.h"
#include "gdfonts.h"
#include "gdfontt.h"
#include "gdfontmb.h"
#include "gdfontl.h"
#ifdef _madomain
#	define seq_name_length 26
#	include "madomain.h"
#else
#	define seq_name_length 20
#endif
static gdImagePtr im_out;       /* declare image */
/* End Modif 14/05/96 Tim Downs */

typedef struct ParamImage_def {
	int BackGroundColour;
	int TextColour;
   int BoldColour;
	int HighColour;
	int LowColour;
	int NeutralColour;
	int LineSize;
	gdFontPtr FontSize;
	int FontHeight;
	int FontWidth;
	int GraduationStep;
	t_Style OutputStyle;
	char DocumentTitle[256];
		} ParamImage_t;

static ParamImage_t Image;
#ifndef _madomain
static const char IMSFformat[56] =
	" Name: %-*s Len: %5ld  Check: %4d  Weight:  %3.2f";
#endif
#define MinWidth   80

static int GetColour (char* Colour)
{
#ifdef HTML
	printf ("Looking for ->%s<-\n",Colour);
#endif
	if (!strncmp(Colour,"red",3))
		return(gdImageColorAllocate(im_out,255,0,0));
	else if (!strncmp(Colour,"blue",4))
		return(gdImageColorAllocate(im_out,0,0,255));
	else if (!strncmp(Colour,"green",5))
		return(gdImageColorAllocate(im_out,0,255,0));
	else if (!strncmp(Colour,"black",5))
		return(gdImageColorAllocate(im_out,0,0,0));
	else if (!strncmp(Colour,"grey",4))
		return(gdImageColorAllocate(im_out,47,79,79));
        else if (!strncmp(Colour,"white",5))
		return(gdImageColorAllocate(im_out,255,255,255));  
        else if (!strncmp(Colour,"grey",4))
		return(gdImageColorAllocate(im_out,190,190,190));  
        else if (!strncmp(Colour,"yellow",6))
		return(gdImageColorAllocate(im_out,255,255, 0));  
        else if (!strncmp(Colour,"pink",4))
		return(gdImageColorAllocate(im_out,255,192,203));  
        else if (!strncmp(Colour,"orange",6))
		return(gdImageColorAllocate(im_out,255,165,0));  
        else if (!strncmp(Colour,"violet",4))
		return(gdImageColorAllocate(im_out,238,130,238));  
		  else if (!strncmp(Colour,"brown",5))
		return(gdImageColorAllocate(im_out,165,42,42));  
        else if (!strncmp(Colour,"cyan",4))
		return(gdImageColorAllocate(im_out,0,255,255));  
	else
		{
#ifdef  HTML
		printf ("no colour corresponding\n");
#endif
		return(0);
		}
}

static void GetFontDimensions(char* Font)
{
	if (!strncmp(Font,"Tiny",4))
		{
#ifdef HTML
		printf("FontSize is :gdFontTiny\n");
#endif
		Image.FontSize = gdFontTiny;
					 Image.FontWidth = gdFontTiny->w;
					 Image.FontHeight = gdFontTiny->h;
		}
	else if (!strncmp(Font,"Small",5))
					 {
#ifdef HTML
		printf("FontSize is :gdFontSmall\n");
#endif
		Image.FontSize = gdFontSmall;
					 Image.FontWidth = gdFontSmall->w;
					 Image.FontHeight = gdFontSmall->h;
					 }
	else if (!strncmp(Font,"MediumBold",10))
					 {
#ifdef HTML
		printf("FontSize is :gdFontMediumBold\n");
#endif
		Image.FontSize = gdFontMediumBold;
					 Image.FontWidth = gdFontMediumBold->w;
					 Image.FontHeight = gdFontMediumBold->h;
		}
#ifndef __BORLANDC__
	else if (!strncmp(Font,"Large",5))
					 {
#ifdef HTML
		printf("FontSize is :gdFontLarge\n");
#endif
		Image.FontSize = gdFontLarge;
					 Image.FontWidth = gdFontLarge->w;
					 Image.FontHeight = gdFontLarge->h;
					 }
	else if (!strncmp(Font,"Giant",4))
					 {
#ifdef HTML
		printf("FontSize is :gdFontGiant\n");
#endif
		Image.FontSize = gdFontGiant;
					 Image.FontWidth = gdFontGiant->w;
					 Image.FontHeight = gdFontGiant->h;
					 }
#endif
	else
					 {
#ifdef HTML
		printf("FontSize is :gdFontMediumBold\n");
#endif
		Image.FontSize = gdFontMediumBold;
					 Image.FontWidth = gdFontMediumBold->w;
					 Image.FontHeight = gdFontMediumBold->h;
		}
}

static int OuvrirFichierConfig (FILE* *f, char *FileName)
{
	char Ligne[256], *l;
#ifndef __BORLANDC__
	if ((*f=fopen(FileName,"rt"))!=NULL) return 1;
#endif
#ifdef HTML
		printf("Unable to open \"%s\"\n",FileName);
		printf("Trying to open \"image.default\"\n");
#endif
	if ((*f=fopen("image.default","rt"))!=NULL)
	{
		strcpy (FileName,"image.default");
		return 1;
	}
	if ((l=getenv("MULTALIN"))!=NULL)
	{
		strcpy (Ligne,l);
		strcat (Ligne,"image.default");
		if((*f= fopen(Ligne,"rt"))!=NULL) strcpy(FileName,Ligne);
	}
	return *f!=NULL;
}

static void ChargerFontSize (char* FileName)
{
		  char Ligne[256];
        int c, UseDefaults = 0;
		  FILE* FICH;

	if (!OuvrirFichierConfig(&FICH,FileName))
	{
		UseDefaults = 1;
#ifdef HTML
		printf("Unable to open \"image.default\"\n");
		printf("Using internal defaults...\n");
#endif
	}
	if(!UseDefaults)
	{
#ifdef HTML
			printf("File Open\n");
#endif

			/* passer les espaces de debut */
			while( isspace(c= fgetc(FICH)) );
        	if( c==EOF )
                	UseDefaults = 1;
        	else
						ungetc( c, FICH );

			/* lire la ligne */
			if( !fgets( Ligne, 255, FICH ) )
						UseDefaults = 1;
			while (!feof(FICH))
			{

						/*--- FontSize  ---*/
						if( !strncmp( Ligne, "[FontSize]", 10 ) )
                	{
									fscanf(FICH,"%s",Ligne);
#ifdef HTML
									printf("FontSize is ->%s<-\n",Ligne);
#endif
									GetFontDimensions(Ligne);
                	}
						fscanf(FICH,"%s",Ligne);
			}
			fclose (FICH);
	}
	else GetFontDimensions("gdFontMediumBold");
}


static void ChargerParametres (char* FileName)
{
	char Ligne[256];
        int c,UseDefaults = 0;
	char Colour[256];
	FILE* FICH;
	printf("Creating a GIF image ");
	if (!OuvrirFichierConfig(&FICH,FileName))
	{
			UseDefaults = 1;
#ifdef HTML
			printf("Unable to open \"image.default\"\n");
#endif
			printf("Using internal defaults...\n");
	}
	if(!UseDefaults)
	{
		printf("Using %s\n",FileName);

				 /* passer les espaces de debut */
			while( isspace(c= fgetc(FICH)) );
			if( c==EOF )
			UseDefaults = 1;
			else
						ungetc( c, FICH );

			/* lire la ligne */
			if( !fgets( Ligne, 255, FICH ) )
			UseDefaults = 1;
		while (!feof(FICH))
		{

				/*--- BackGroundColour ---*/
			/* First colour defined is background color in image */
				if(! strncmp( Ligne, "[BackGroundColour]", 18 ) )
				{
               			fscanf(FICH,"%s",Ligne);
				strcpy(Colour,Ligne);
								Image.BackGroundColour= GetColour(Colour);
      		 	}
				/*--- TextColour  ---*/
				else if( !strncmp( Ligne, "[TextColour]", 12 ) )
				{
					fscanf(FICH,"%s",Ligne);
					strcpy(Colour,Ligne);
					Image.TextColour= GetColour(Colour);
				}
				/*--- BoldColour  ---*/
				else if( !strncmp( Ligne, "[BoldColour]", 12 ) )
				{
					fscanf(FICH,"%s",Ligne);
					strcpy(Colour,Ligne);
					Image.BoldColour= GetColour(Colour);
			}
				/*--- HighColour  ---*/
				else if( !strncmp( Ligne, "[HighColour]", 12 ) )
				{
               			fscanf(FICH,"%s",Ligne);
				strcpy(Colour,Ligne); 
							Image.HighColour= GetColour(Colour);
			}
				/*--- LowColour  ---*/
				else if( !strncmp( Ligne, "[LowColour]", 11 ) )
				{
								fscanf(FICH,"%s",Ligne);
				strcpy(Colour,Ligne);
							Image.LowColour= GetColour(Colour);
			}
				/*--- NeutralColour  ---*/
				else if( !strncmp( Ligne, "[NeutralColour]", 15 ) )
				{
								fscanf(FICH,"%s",Ligne);
				strcpy(Colour,Ligne);
							Image.NeutralColour= GetColour(Colour);
			}
						fscanf(FICH,"%s",Ligne);

		}
		fclose (FICH);
	}
	else
	{
		Image.BackGroundColour= GetColour("white");
		Image.TextColour= GetColour("black");
		Image.HighColour= GetColour("red");
		Image.LowColour= GetColour("blue");
		Image.NeutralColour= GetColour("grey");
		Image.GraduationStep=10;
		Image.OutputStyle= Normal;
	}
	if (Image.BoldColour == Image.BackGroundColour)
		Image.BoldColour= Image.HighColour;

}
static void IEcrit1LigneSeq(int x,int y,t_pSeq Seq,t_pSeq Cons,t_pSeq First_Seq,int IsCons,int IsFirst_Seq ,
	t_indL Start,t_indL Finish,t_pParamAlign Param,int BlockSize)
/* ecrit BlockSize residus de la sequence p_Seq */
{


	int c,c1,c2,h,h1,Colour;
	t_indL Pos;
	int  Count=0;
	t_Seq *S,*C,*F;
	char Gap[4]="  .";

	Gap[1]=Insert;
	x += (seq_name_length + 2) * Image.FontWidth ;


	for( Pos=Start, S=Seq+Pos,C=Cons+Pos,F=First_Seq+Pos; Pos<=Finish; Pos++,S++,C++,F++ )
	{
		c= S->Car;
		h= Param->Hom[Param->NumSymb[toupper(c)]];
		c1= C->Car;
		if (strchr(Gap,c1)) h1=c1;
		else h1= Param->Hom[Param->NumSymb[toupper(c1)]];
		c2= F->Car;
		if (!strchr(Gap,c)) /* le caractere a afficher n'est pas le '-' */
		{
			if (h==h1) /* residu = consensus */
				if (!islower(c1)) /* consensus fort */
					Colour=Image.HighColour;
				else /* consensus faible */
					Colour=Image.LowColour;
			else if (IsCons) Colour = Image.NeutralColour;
			else /* residu != consensus */
			{
				Colour = Image.NeutralColour;
				if (Image.OutputStyle==Case) c = tolower(c);
			}
		}
		else
		{
			Colour=Image.NeutralColour;
		}
		if ((Image.OutputStyle==Difference)&&(!IsFirst_Seq) && (!strchr(Gap,c)))
		{
			if (c == c2)
				c = '.';

		}
		gdImageChar(im_out,Image.FontSize,x,y,c,Colour);
		x += Image.FontWidth;
		Count++;
		if ( Count==BlockSize && Pos!=Finish)
		{
			x += Image.FontWidth;
			Count= 0;
		}
	}
}
#ifndef _madomain
static int IGCGchar (int c)
{
	if ((c==Insert)||(c==' ')) return '.';
	else return toupper(c);
}

static int IGCGCheckSum (t_pLigneSeq Seq)
{
	t_indL i;
	int Count;
	long Check;
	t_pSeq s;

	for (i=1, s= Seq->Seq+1, Count=1, Check=0; i<=Seq->Taille; i++,s++,Count++)
	{
		Check += Count * IGCGchar(s->Car);
		if (Count == 57) Count=0;
	}

	return (int)(Check % 10000);
}
#endif
/* ????????? */
void CreerImage(char* FileName ,t_pDescripteurSequence DS,t_pParamAlign Param,
	int complete,int Group2 )
/* sauvegarde au format MSF les sequences de DS dans OUTFILE */
{
	int BlockSize= 30000;     /* decoupage des residus par bloc de 10 */
	int x,y;                 /* x : colonne, y : ligne */
	int nc;						/* nb de symboles */
	int crn;		/* nb of copyright lines */
	int nh;	/* nb of header lines */
	char *s, name[seq_name_length+1];
	char chaine[80];           /* for conversion of fprintf to gdImageString */
	char FileCfgName[256],FileImageName[256];
	t_ExtensionFichier Extension;
		  /* color references
		  int white,black,blue,red,green;*/
		  /* declare image file */
		  FILE *image_file;

	t_indL Start, Finish, Limit ,j;
	int Spaces,Nb_Rows,First_Time=1;
	long MaxLen;
	t_indN i, *n, ind;
	t_pLigneSeq tmp;
#ifndef _madomain
	float Weight0 = ((float)DS->NbSeq) / DS->Sequence[0].Weight;
#endif
	t_Seq *First_Seq;
	long Width;
	/* initialisation de la fin de name */
	name[seq_name_length]=0;

	/* creation des noms de fichier */
	strcpy( Extension, ".gif" );
	strcpy(FileImageName, ChangeExt( FileName,Extension) );
	strcpy( Extension, ".cfg" );
	strcpy(FileCfgName, ChangeExt( FileName,Extension) );

	MaxLen= DS->Sequence[0].DerCol;

	/* create image */
	PEGetOutputStyleParam(&Image.OutputStyle,&Image.LineSize,&Image.GraduationStep);
	ChargerFontSize(FileCfgName); /* necessaire de charger avant creation */
	Nb_Rows = MaxLen / Image.LineSize + 1 ;
	crn = PrintCopyRight(NULL);
	nc = CalculeSymboles(Param);
	if (nc) nc++;
	if (Nb_Rows == 1)
		Image.LineSize = MaxLen;
	if (complete) nh = 9 +nc + crn
#ifndef _madomain
		+ DS->NbSeq
#endif
	;else nh = 0;
	Width = max (MinWidth,(seq_name_length+2+Image.LineSize+3));
		  im_out = gdImageCreate( Width*Image.FontWidth,
		(nh + Nb_Rows * (DS->NbSeq + 4)) *Image.FontHeight);
	/* end create image */

	/*	Initialise tous parametres de image */
	ChargerParametres(FileCfgName); /* Les couleurs ne peuvent etre definies qu'apres creation */

	x=0;y=0;
	if (complete)
	{
		  for (y=0,i=0;i<crn;i++,y+=Image.FontHeight)
			gdImageString(im_out,Image.FontSize,x,y,CopyRight[i],Image.TextColour);
		  sprintf(chaine,"Symbol comparison table: %s", Param->NomCoeff);
		  gdImageString(im_out,Image.FontSize,x,y,chaine,Image.TextColour);
		  y += Image.FontHeight;
		  sprintf(chaine,"Gap weight: %g", (float)Param->Gap);
		  gdImageString(im_out,Image.FontSize,x,y,chaine,Image.TextColour);
		  y += Image.FontHeight;
		  sprintf(chaine,"Gap length weight: %g", (float)Param->Gap2);
		  gdImageString(im_out,Image.FontSize,x,y,chaine,Image.TextColour);
		  y += Image.FontHeight;
		  sprintf(chaine,"Consensus levels: high=%d%% low=%d%%",(int)Param->ConsHlevel,
			(int)Param->ConsLlevel);
		  gdImageString(im_out,Image.FontSize,x,y,chaine,Image.TextColour);
		  y += Image.FontHeight;
		  if (nc)
		  {
			gdImageString(im_out,Image.FontSize,x,y,"Consensus symbols:",Image.TextColour);
			y += Image.FontHeight;
			for (nc=1;nc<=DERLETTREMAX;nc++) if ((s=ConsSymb[nc])[0])
			{
				sprintf(chaine," %c is anyone of %s",s[0],s+1);
				gdImageString(im_out,Image.FontSize,x,y,chaine,Image.TextColour);
				y += Image.FontHeight;
			}
		  }
		  sprintf(chaine," MSF: %6ld    Check:   0         ..",(long)MaxLen);
		  gdImageString(im_out,Image.FontSize,x,y,chaine,Image.TextColour);
		  y += Image.FontHeight;

#ifndef _madomain
	/* ecrire l'entete de fichier */
	for( i=1,n= DS->NoSeq; i<= DS->NbSeq; i++)
	{
		tmp=DS->Sequence+ (n ? *(++n) : i);
		sprintf(chaine, IMSFformat,LONGNOM,
			tmp->NomSeq,(long)MaxLen, IGCGCheckSum(tmp), tmp->Weight*Weight0);
		gdImageString(im_out,Image.FontSize,x,y,chaine,Image.TextColour);
			y += Image.FontHeight;

	}
	sprintf (chaine, IMSFformat,LONGNOM,
		DS->Sequence[0].NomSeq, (long)MaxLen, IGCGCheckSum(DS->Sequence), 0);
	gdImageString(im_out,Image.FontSize,x,y,chaine,Image.TextColour);
		  y += 2*Image.FontHeight;
	gdImageString(im_out,Image.FontSize,x,y,"//",Image.TextColour);
		  y += 2*Image.FontHeight;
#endif
	}
	/* ecrire les sequences */
	Start= 1;
	Limit= 10;
	do
	{

		while( Start>=Limit)
			Limit+= 10;
		First_Time=1;

		/* ecrire les numeros de debut de ligne */
		sprintf(chaine,"%*c%-5ld",seq_name_length+2,' ',(long)Start);
		gdImageString(im_out,Image.FontSize,x,y,chaine,Image.TextColour);
		/* was previously += */
		x = (seq_name_length + 2) * Image.FontWidth ;
		
		/* ecris le numero de fin de ligne */
		/* Tim Downs 15/05/96 */
		Finish = (Start+Image.LineSize)-1;
		if (Finish > MaxLen)
			Finish=MaxLen;
		Spaces=(Finish-Start)/BlockSize;
		for (j = 1; j < (Finish -Start + Spaces); j++)
		{
			if (!(j % Image.GraduationStep))
			{
				x -= 4*Image.FontWidth;
				sprintf(chaine,"%5ld",(long)Start+j-1);
				gdImageString(im_out,Image.FontSize,x,y,chaine,Image.TextColour);
				x += 4*Image.FontWidth;
			}
			x += Image.FontWidth;
		}
		x -= 3*Image.FontWidth;
		sprintf(chaine,"%5ld", (long)Finish);
		gdImageString(im_out,Image.FontSize,x,y,chaine,Image.TextColour);
		/* End Tim Downs 15/05/96 */

		/* display graduation */
		y += Image.FontHeight;
		x = (seq_name_length + 2) * Image.FontWidth ;
		gdImageChar(im_out,Image.FontSize,x,y,'|',Image.TextColour);
		for (j = 1; j <= (Finish -Start + Spaces); j++) 	
		{
			if (j != 1)
				if (j % Image.GraduationStep)  /* Step 10 between each + */
					gdImageChar(im_out,Image.FontSize,x,y,'-',Image.TextColour);
				else
					gdImageChar(im_out,Image.FontSize,x,y,'+',Image.TextColour);

			x += Image.FontWidth;
				
		}
		gdImageChar(im_out,Image.FontSize,x,y,'|',Image.TextColour);


		x =0;
		y += Image.FontHeight;

		/* ecrire les NbSeq lignes de sequences */
		for( i=1, n= DS->NoSeq; i<=DS->NbSeq; i++ )
		{
			ind = (n ? *(++n) : i);
			tmp= DS->Sequence+ind ;
			if (i==1) 	First_Seq = tmp->Seq;
			if( Start <= tmp->Taille )
			{
#ifdef _madomain
				strncpy(name,codearr[tmp->name],seq_name_length);
#else
				strncpy(name,tmp->NomSeq,seq_name_length);
#endif
				sprintf(chaine,"%*s  ",seq_name_length,name);
				if (ind < Group2)
					gdImageString(im_out,Image.FontSize,x,y,chaine,Image.BoldColour);
				else
					gdImageString(im_out,Image.FontSize,x,y,chaine,Image.TextColour);
				switch (Image.OutputStyle)
				{
				case 0: /* normal  */
					IEcrit1LigneSeq(x,y, tmp->Seq, DS->Sequence[0].Seq,First_Seq ,0
						,0, Start, Finish, Param,BlockSize);
					break;
				case 1: /* case */
					IEcrit1LigneSeq(x,y, tmp->Seq, DS->Sequence[0].Seq,First_Seq ,0
						,0, Start, Finish, Param,BlockSize);
					break;
				case 2: /* difference */
					if (First_Time)
					{
						IEcrit1LigneSeq(x,y, tmp->Seq, DS->Sequence[0].Seq,
						First_Seq ,0,1, Start, Finish, Param,BlockSize);
						First_Time=0;
					}
					else
					{
						IEcrit1LigneSeq(x,y, tmp->Seq, DS->Sequence[0].Seq,
						First_Seq ,0,0, Start, Finish, Param,BlockSize);
					}
					break;
			}
	
                       		y += Image.FontHeight;
			}
		}
		/* ecrire le consensus */
		sprintf(chaine,"%*s  ",seq_name_length,DS->Sequence[0].NomSeq);
		gdImageString(im_out,Image.FontSize,x,y,chaine,Image.TextColour);
		IEcrit1LigneSeq(x,y, DS->Sequence[0].Seq,DS->Sequence[0].Seq,First_Seq,1,0,
			Start, Finish, Param,BlockSize);
		y += Image.FontHeight;
		


	/*	fputc( '\n', OUTFILE ); */
		y += Image.FontHeight;
		Start+= Image.LineSize;
	} while( Finish < MaxLen);
	 /* save image */
        image_file = fopen (FileImageName,"wb");
        gdImageGif (im_out,image_file);
		  gdImageDestroy(im_out);
		  fclose (image_file);
}
#endif
