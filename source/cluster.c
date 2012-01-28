/*-------------------------------------------------------------------------
	CLUSTER.C
			Lecture et ecriture de fichiers d'arbres.

-------------------------------------------------------------------------*/
#include <string.h>
#include <ctype.h>

#include "ma.h"

#include "util.h"
#include "basproc1.h"
#include "commande.h"
#include "msgerr.h"


static int ReadArbre( t_ClassHier Arbre, t_DessinArb Dessin, t_pDescripteurSequence DS,
			t_CheminFichier FileName );
static int PrintArbre (FILE *FICH, t_pDescripteurSequence DS, t_ClassHier Arbre,
	t_DessinArb Dessin);
static int ReadScore( t_ClassHier Arbre, t_DessinArb Dessin, t_pDescripteurSequence DS,
			t_CheminFichier FileName );



/*====================  Lecture de l'arbre =======================*/

static t_indN Numero( t_pDescripteurSequence DS, t_StringNom Nom )
{
	t_indN Prot;
	t_pLigneSeq tmp;

	for (Prot=1,tmp=DS->Sequence+1;(Prot > 0) && (Prot <= DS->NbSeq)
			&& strcmp(tmp->NomSeq,Nom); Prot++,tmp++);

	if (Prot > DS->NbSeq)
		return 0;
	else
		return Prot;
}

static int Entre (t_indN pere, t_indN lui, t_indN fils, t_ClassHier Arbre)
{
	while( (pere != fils) && (pere != lui) )
		pere= Arbre[pere].SeqSuiv;
	return( pere == lui);
}

static void CalculLigne (t_ClassHier Arbre, t_DessinArb Dessin, t_indN NbSeq,
	t_indN Num, t_indN Largeur, char *Ligne)
{
	t_indN point,rang,papa,fils;
	int fin;
	char *l;
	t_DessinArb d;

	fin = min ((int)strlen(Ligne), DebutGraphe-2);
	sprintf (Ligne + fin, "%*s", Largeur-fin," ");

	point = DebutGraphe;
	l=Ligne+point-1;
	fin = 0;
	for (rang = NbSeq,d=Dessin+rang;(d->Col<point)&&(rang>0);rang--,d--);
	if (!rang) return;

	do {
		if (point==d->Col)
		{/* there is a node in this column */
			fils = d->Noeud;
			papa = Arbre[fils].Pere;
			if (Num==papa) *l= DepartBas;/* Num is the first line of the node */
			else if (Num==fils)/* Num is the second line of the node */
			{
				if ((*l==Verticale)||(*l==DepartGauche)||(*l==DepartBas))
					*l=DepartGauche;
				else *l=BasDroite;
				fin = 1;
			}
			else if (Entre(papa,Num,fils,Arbre))/*Num is between the 2 node lines*/
				if ((*l==' ')||(*l==Verticale)) *l=Verticale;
				else *l = DepartGauche;
			else if (!fin && (*l==' ')) *l=Trait;
			do {
				rang--; d--;
			} while ((d->Col<point)&&(rang>1));
		}
		else if (!fin) *l= Trait;/* no node in this column */
		if ((rang<2)||(point<d->Col)) {point++;l++;};
	} while (*l &&(rang>1));
}

/*====================  Sauvegarde de l'arbre =======================*/

static int ActEcrArbre( t_indN Groupe1, t_indN Groupe2, t_indN NbSeq, t_indN* Groupe,
				  t_ClassHier Arbre, void* V )
{
	t_indN Seq, i;
	t_indN *Profondeur = (t_indN *)V;

	Profondeur[Groupe1]++;
	for( i= 1,Seq=Groupe2; i<Groupe[Groupe2]; i++ )
		Seq= Arbre[Seq].SeqSuiv;
	Profondeur[Seq]--;
	Groupe[Groupe1]+= Groupe[Groupe2];
	Groupe[Groupe2]= 0;
	return vrai;
}


static int EcritArbre( t_ClassHier Arbre, t_pDescripteurSequence DS , FILE *FICH)
{
	t_indN *Profondeur, *Groupe, *i;
	t_indN N,Count,Seq1;
	int Err;
	t_pLigneSeq tmp;

	if ((Profondeur=calloc(DS->NbSeq+1,sizeof(t_indN)))==NULL)
	{
		TraiteErr(8);
		return faux;
	}
	if ((Groupe=calloc(DS->NbSeq+1,sizeof(t_indN)))==NULL)
	{
		TraiteErr(8);
		free (Profondeur);
		return faux;
	}
	memset( Profondeur, 0, (DS->NbSeq+1)*sizeof(t_indN) );
	for (N=0,i=Groupe;N<=DS->NbSeq;N++,i++) *i=1;
	TraiteBranche (Arbre, DS->NbSeq, Groupe, (void*)Profondeur,ActEcrArbre);
	for (Err=0,Seq1=Arbre[0].SeqSuiv,Count=0;Seq1 && !Err;)
	{
		N=Profondeur[Seq1];
		tmp = DS->Sequence+Seq1;
		Count+=abs(N)+ strlen(tmp->NomSeq);
		for (;N>0;N--)  if (Err = (fprintf(FICH,"(")==EOF)) break;
		if (Err) break;
		if (Err=(fprintf (FICH,"%s",tmp->NomSeq)==EOF)) break;
		for (;N<0;N++)  if (Err=(fprintf(FICH,")")==EOF)) break;
		if (Err) break;
		Seq1 = Arbre[Seq1].SeqSuiv;
		if (Seq1)
		{
			if (Err=(fprintf(FICH,",")==EOF)) break;
			if (Count > 72)
			{
				if (Err=(fprintf(FICH,"\n")==EOF)) break;
				Count=0;
			}
		}
		else if (Err=(fprintf (FICH,";\n")==EOF)) break;
	}
	free (Groupe);
	free(Profondeur);
	return !Err;
}

void SauvArbre (t_ClassHier Arbre, t_DessinArb Dessin,t_pDescripteurSequence DS,
	char* Extension)
{
	FILE *FICH;
	int OK;
	char *FileName = ChangeExt(PEGetNomFichSeq(),Extension);
	int format = PEGetFormatArbre();
	Backup (FileName);
	if ((FICH = fopen (FileName,"wt"))==NULL) TraiteErr (27,FileName);
	else
	{
		MessageAction (hcActWrite, FileName);
		if (format)  OK=PrintArbre(FICH,DS,Arbre,Dessin);
		else OK = EcritArbre (Arbre, DS, FICH);
		if (!OK) TraiteErr (16,FileName);
		fclose (FICH);
	}
}

static int ActDessin (t_indN Groupe1, t_indN Groupe2, t_indN Nbseq, t_indN *Groupe,
	t_ClassHier Arbre, void * V)
{
	t_indN* p1= Groupe+Groupe1;
	t_indN* p2= Groupe+Groupe2;

	if (*p1 > *p2) *p2 = *p1;
	*p1 = *p2 + 1;
	return vrai;
}

static int CalculeDessin (t_indN NbSeq, t_ClassHier Arbre, t_DessinArb Dessin)
{
	t_indN *Place, *p;
	t_DessinArb d;
	t_indN i,j,Mini,Fils;

	if ((Place = calloc (NbSeq+1,sizeof(t_indN)))==NULL) return faux;
	TraiteBranche (Arbre, NbSeq, Place, NULL, ActDessin);
	for (Fils=0,i=NbSeq,d=Dessin+i;i>1;i--,d--)
	{
		for (Mini=NbSeqMax+1,j=1,p=Place+1;j<=NbSeq;j++,p++)
			if (*p < Mini)
			{
				Mini=*p;
				Fils=j;
			};
		if (Mini > TailleGraphe) Mini = TailleGraphe;
		d->Col = Mini + DebutGraphe;
		d->Noeud =  Fils;
		Place[Fils] = NbSeqMax+1;
	}
	free (Place);
	return vrai;
}

static int FindCh( FILE* FICH, char *Ch )
{
	int lu = *Ch;
	if( lu != ')' )
		do
		{
			if( (lu= fgetc(FICH)) == EOF ) return faux;
		}
		while(isspace(lu)||(lu == '(')||(lu == ','));
	*Ch = lu;
	return vrai;
}

static int CompleteArb( t_ClassHier Arbre, t_indN *Seq, t_indN DPere, t_indN i )
{
	if( (i != 0) && (Arbre[i].Pere == 0) && (i != DPere) )
	{
		Arbre[*Seq].SeqSuiv= i;
		*Seq= i;
		Arbre[*Seq].Pere= DPere;
		return vrai;
	}
	else
		i= 0;
	return faux;
}


static int LitArbre( t_ClassHier Arbre, t_DessinArb Dessin, t_pDescripteurSequence DS, FILE* FICH )
{
	t_indN DPere, Seq, X, Count;
	char Ch, Chs[2];
	int OK;
	t_StringNom S;
	int Res= faux;

	memset( Arbre,0, (DS->NbSeq+1)*sizeof(t_Branche));
	memset(Dessin,0, (DS->NbSeq+1)*sizeof(t_Dessin));

	DPere= 0;
	Seq= 0;
	OK= vrai;
	Ch= '\0';
	Count= DS->NbSeq;


	while( OK && FindCh(FICH, &Ch) && (Ch != ';') )
	{
		if( Ch == ')' )
		{
			Ch= '\0';
			DPere= Arbre[DPere].Pere;
		}
		else
		{
			Count--;
			S[0]= '\0';
			do
			{
				sprintf (Chs, "%c",Ch);
				strcat (S,Chs);
				while ((Ch= fgetc(FICH))=='\n');
			}
			while( (Ch != ',') && (Ch != ')') && ((int)strlen(S) < LONGNOM) );
			X= Numero( DS, strupr(S) );
			if( OK= CompleteArb( Arbre, &Seq, DPere, X) )
				DPere= Seq;
		}
	}

	if( OK= OK && (DPere == Arbre[0].SeqSuiv) )
	{
		if( Count != 0 )
			for( Count= 1; Count<=DS->NbSeq; Count ++ )
				CompleteArb( Arbre, &Seq, DPere, Count );
		CalculeDessin (DS->NbSeq, Arbre, Dessin);
		Res= vrai;
	}
	else
	{
		memset( Arbre, 0, (DS->NbSeq+1)*sizeof(t_Branche) );
		memset(Dessin,0, (DS->NbSeq+1)*sizeof(t_Dessin));
	}
	return Res;
}

int GetArbre ( t_pDescripteurSequence DS, t_ClassHier Arbre, t_DessinArb Dessin)
{
	int FileFormat = PEGetKindOfTreeInput();
	char* TreeFile= PEGetTreeFile();
	if( FileFormat==1 ) 	/* on doit recuperer un arbre deja calcule */
		return ReadArbre( Arbre, Dessin, DS, TreeFile );
	else if (FileFormat==2)
	{
		if (!ReadScore (Arbre, Dessin, DS, TreeFile))
			return 0;
		else
		{
			SauvArbre( Arbre, Dessin,DS,".clu");
			return 1;
		}
	}
	else return -1;
}

static int LireScore2a2 (FILE *ScoreFile, t_score* Score)
{
	float lu;
	int Res;
	if (Res = (fscanf (ScoreFile,"%g", &lu)==1))
		*Score = (t_score)lu;
	return Res;
}

static int LitScore (FILE *ScoreFile, t_pDescripteurSequence DS,
	t_ClassHier Arbre, t_DessinArb Dessin)
{
	t_score ** Score, ** sz;
	t_indN i,Seq1, *Ordre, *o, *or;
	int OK;
	t_StringNom M;


	if ((Ordre=calloc (DS->NbSeq,sizeof (t_indN)))==NULL)
	{
		TraiteErr (8);
		return vrai;
	}
	if ((Score=calloc (DS->NbSeq-1,sizeof (t_score*)))==NULL)
	{
		TraiteErr (8);
		free (Ordre);
		return vrai;
	}
	for (OK=vrai,i=1,sz=Score;OK && (i<DS->NbSeq);i++,sz++)
		if ((*sz=calloc(i,sizeof(t_score)))==NULL) break;

	if (i==DS->NbSeq)
	{
		for (i=0,o=Ordre;i<DS->NbSeq;i++,o++)
		{
			if ((fscanf(ScoreFile,"%s",M)!=1)||((*o=Numero(DS,strupr(M)))==0))
			{
				OK=faux;
				break;
			}
			for (Seq1=0,or=Ordre;(Seq1<i)
				&&LireScore2a2(ScoreFile,Place(*o,*or,Score));Seq1++,or++);
			OK = (Seq1==i);
		}
		OK = OK && CalculArbre (Score, Arbre,Dessin, DS->NbSeq,0);
	}
	else { TraiteErr (8); OK=vrai;}
	for (i=1,sz=Score;(i<DS->NbSeq)&& *sz;i++,sz++) free (*sz);
	free (Score);
	free (Ordre);
	return OK;
}

static int PrintArbre (FILE *FICH, t_pDescripteurSequence DS, t_ClassHier Arbre,
	t_DessinArb Dessin)
{
	t_indN Seq1;
	char ligne[LargeurGraphe+1];


	for (Seq1=Arbre[0].SeqSuiv; Seq1; Seq1=Arbre[Seq1].SeqSuiv)
	{
		strncpy (ligne,DS->Sequence[Seq1].NomSeq,DebutGraphe-2);
		CalculLigne (Arbre,Dessin,DS->NbSeq, Seq1,LargeurGraphe,ligne);
		if (fprintf (FICH,"%s\n",ligne)==EOF) return faux;
	}
	return vrai;
}

static int ReadArbre( t_ClassHier Arbre, t_DessinArb Dessin, t_pDescripteurSequence DS,
			t_CheminFichier FileName )
{
	FILE* FICH;
	int Erreur= faux;

	FICH= fopen(FileName,"rt");

	MessageAction( hcActRead, FileName);

	Erreur= ( !FICH || feof(FICH) || !LitArbre( Arbre, Dessin, DS, FICH) );

	if (FICH) fclose( FICH );

	if( Erreur )
	{
		TraiteErr(28, FileName );
		return faux;
	}
		return vrai;
}

static int ReadScore( t_ClassHier Arbre, t_DessinArb Dessin, t_pDescripteurSequence DS,
			t_CheminFichier FileName )
{
	FILE* FICH;
	int Erreur= faux;

	FICH= fopen(FileName,"rt");

	MessageAction( hcActRead, FileName);

	Erreur= ( !FICH || feof(FICH) || !LitScore( FICH,DS,Arbre, Dessin) );

	if (FICH) fclose( FICH );

	if( Erreur )
	{
		TraiteErr(28, FileName );
		return faux;
	}
		return vrai;
}







