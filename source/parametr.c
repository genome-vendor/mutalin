/*------------------------------------------------------------------------
	PARAMETR.C
				Gestion des parametres d'alignement

------------------------------------------------------------------------*/
/*
14 juin 1994 : Gap penalty = Gap + length * Gap2
	Gap2 is 0 by default and can be user defined
*/

#include <string.h>
#include <ctype.h>

#include "parametr.h"
#include "msgerr.h" 


void AjouteSymb (char* Ens, int c)
{
	char Chs[2];
	if (strchr(Ens,c)) return;
	sprintf (Chs,"%c",c);
	strcat (Ens,Chs);
}

/*------------------------------------------------------------------------*/
static int SymbInterdits( int c )
{
	return !isprint(c)||islower(c)||isdigit(c)||isspace(c)
		||(c=='*')||(c==Insert)||(c==DebutLigne);
}

/*------------------------------------------------------------------------*/

void InitParam( t_pParamAlign Param )
/* initialisation de la structure de type t_paramAlign */
{
	int i;

	/* symboles */

	for( i= 0; i<= 127; i++ )     /* ! a controler: ou *(Param->Coeff) */
		if(  SymbInterdits (i ))
			Param->NumSymb[i]= NoSymb;		/* symbole interdit */
		else
			Param->NumSymb[i]= FreeSymb;		/* symbole libre */

	memset( Param->Symb, 0, DERLETTREMAX*sizeof(byte) );
	Param->Symb[0]= ' ';
	Param->NumSymb[' ']=0;
	Param->Symb[1]=Insert;
	Param->NumSymb[Insert]=1;
	Param->PremLettre = 2;

	/* table des coefficients */
	memset( Param->Coeff, 0, DERLETTREMAX*DERLETTREMAX*sizeof(t_score));

	/* table des homologies */
	for( i=0; i<=DERLETTREMAX; i++ )
		Param->Hom[i]= i;

	Param->DerLettre= 1;

	Param->Gap= 0;
	Param->Gap2=0;
	Param->GapExtT=0;

   Param->Weighted = vrai;
	Param->UNEITER= faux;
	Param->ScMethod= ABSOLU;
}

#ifndef clus2dom
static void SupprimDsSeq(t_pDescripteurSequence DS, char* Ens )
{
	t_indN i;
	t_indL pos;
	t_pLigneSeq tmp;
	t_Seq *s;

	for( i=0,tmp=DS->Sequence+1; i< DS->NbSeq; i++,tmp++ )
	{
		for( pos= tmp->PremCol,s=tmp->Seq+pos; pos<= tmp->DerCol;)
			if (strchr (Ens, s->Car))
			{
				if( tmp->Taille > pos )	memmove( s, s+1,tmp->Taille-pos );
				tmp->Taille--;
				tmp->DerCol--;
			}
			else
			{
				pos++;
				s++;
			}
	}
}

/*------------------------------------------------------------------------*/
void VerifiSymb(t_pDescripteurSequence DS, t_pParamAlign Param, char* EnsLet)
{
	char Ens[128];
	byte *c;
	int L;

	Ens[0]='\0';
	for( L= 1, c=Param->NumSymb+1; L<= 127; L++,c++ )
		if(((*c<0)||(*c > Param->DerLettre) )&& strchr (EnsLet, L ))
			{
				if( TraiteErr( 50,L) == cmOK ) AjouteSymb (Ens, L);
				else
				{
					InitParam( Param );
					Param->NomCoeff[0]= '\0';
					return;
				}
			}
	if( strcmp(Ens,"") )	SupprimDsSeq( DS, Ens );
}
#endif
#ifndef _domainer
int ExisteSequence( t_pDescripteurSequence DS, t_indN DerSeq )
{
	register int i;
	t_pLigneSeq tmp;
	char *Nom = DS->Sequence[DerSeq].NomSeq;

	for( i=1,tmp=DS->Sequence+1; i<DerSeq; i++,tmp++ )
		if(!strcmp(tmp->NomSeq,Nom))
			return vrai;

	return faux;
}
#endif

