#include <string.h>
#include <math.h>

#include "deftypes.h"
#include "util.h"
#include "parametr.h"
#include "disk.h"

#define TRUNC(X) (int)((X) + ((X>0) ? 0.5 : -0.5))
/*
float Frequency[26] = {	0.0759,0.0000,0.0197,0.0519,0.0616,0.0398,0.0727,0.0233,
								0.0525,0.0000,0.0587,0.0000,0.0917,0.0437,0.0000,0.0520,
								0.0408,0.0519,0.0710,0.0589,0.0000,0.0648,0.0138,0.0000,
								0.0326,0.0000};*/
 /* Frequency from Robinson and Robinson, 1991, PNAS, 88, 8880-8884*/
float Frequency[26] = {0.0781,0.0000,0.0192,0.0536,0.0630,0.0386,0.0738,0.0220,
							  0.0514,0.0000,0.0574,0.0902,0.0224,0.0449,0.0000,0.0520,
							  0.0426,0.0513,0.0712,0.0584,0.0000,0.0644,0.0133,0.0000,
							  0.0322,0.0000};

FILE *Infile, *Outfile;
char source[79], target[79], ope[5];
t_ParamAlign Param;

unsigned int TraiteErr(int NoErr, char *NomFich)
{
	printf ("%s contains no symbol comparison table.\n",NomFich);

	return !0;
}

int main (int argc, char* argv[])
{
	int Error = 1;
	int factor,i,j;
	byte *c, *ni, *nj;
	t_LignCoeff *d;
	t_score *k, mini, maxi;
	float *pi, *pj, s, s2, s1,q, ko, d1,d2,*Histo;

/* Frequency correction */
	for (i=0,pi=Frequency,q=0;i<26;i++,pi++) q+= *pi;
	for (i=0,pi=Frequency;i<26;i++,pi++) *pi /= q;

	strcpy (source, argv[2]);
	InitParam (&Param);
	if (argc<3) printf ("Usage : tabedit operation source [target]\n");
	else if ((Infile = fopen (source, "rt")) == NULL)
			printf ("Cannot read %s",source);
	else
	{
		Error = !LitFichCoeff(Infile, source, &Param);
		fclose(Infile);
		if (Error) return Error;
		Error = 1;
		strcpy (target,argv[argc>3 ? 3 : 2]);
		if (argc==3) strcpy (target, ChangeExt (target, ".tab"));
		if ((Outfile = fopen (target, "wt")) == NULL)
			printf ("Cannot open %s to save the new table\n",target);
		else
		{
			strcpy (ope, argv[1]);
			if (Error = (strspn (ope,"N+-x/") != 1)
				|| (sscanf (ope+1, "%d", &factor) != 1) || (factor==0))
			{
				printf ("first argument must be oX where \n");
				printf("o is one of N + - x/ and X is a non-zero integer value\n");
				printf ("exemple : +5  => 5 will be added to any entry\n");
				printf ("          N2  => the entries wille be normalized with mean 0 and s.d. 2\n");
			}
			else
			{
				fprintf (Outfile,">Modified table from %s with the operation %s\n",
					source, ope);
				printf ("Modified table from %s with the operation %s\n",
					source, ope);
				if (ope[0]=='N')
				{
					s  = 0.0;
					s2 = 0.0;
					for (i=0,pi=Frequency, ni=Param.NumSymb+'A';i<26;i++,pi++,ni++)
					if (*ni != -127)
					{
						d = Param.Coeff+ *ni;
						for (j=0,pj=Frequency,nj=Param.NumSymb+'A';j<26;j++,pj++,nj++)
						if (*nj != -127)
						{
							ko = (*d)[*nj];
							s += q = *pi * *pj * ko;
							s2+= q * ko;
						}
					}
					s2 = sqrt (s2- s*s) / factor;
					fprintf (Outfile,">Coefficient are normalized by (x%+.1g)/%.1g (mean 0, s.d. %d)\n",
						-s, s2, factor);
					printf ("Coefficient are normalized by (x%+.1g)/%.1g \n",
						-s, s2, factor);
				}
				for (i=Param.PremLettre, c=Param.Symb+i;i<=Param.DerLettre;i++,c++)
					fprintf (Outfile,"   %c",*c);
				mini = 500;
				maxi =-500;
				for (i=Param.PremLettre,d=Param.Coeff+i;i<=Param.DerLettre;i++,d++)
				{
					fprintf (Outfile,"\n");
					for (j=Param.PremLettre; j<i; j++) fprintf (Outfile,"    ");
					for (k=*d+j;j<=Param.DerLettre;j++,k++)
					{
						switch (ope[0])
						{
							 case 'N' : *k = TRUNC((*k - s) / s2); break;
							 case '+' : *k += factor; break;
							 case '-' : *k -= factor; break;
							 case 'x' : *k *= factor; break;
							 case '/' : *k /= factor; break;
						}
						fprintf (Outfile,"%4d",*k);
						if (mini >*k) mini = *k;
						if (maxi <*k) maxi = *k;
					}
				}
				Error = 0;
				s  = 0.0;
				s1 = 0.0;
				s2 = 0.0;
				d1=0.0;
				d2=0.0;
				for (i=0,pi=Frequency, ni=Param.NumSymb+'A';i<26;i++,pi++,ni++)
				if (*ni != -127)
				{
					d = Param.Coeff+ *ni;
					for (j=0,pj=Frequency,nj=Param.NumSymb+'A';j<26;j++,pj++,nj++)
					{
						if ((*nj == -127)||(*nj<*ni)) continue;
						if (*nj!=*ni)
						{
							q = 2.0 * *pi * *pj;
							ko=(*d)[*nj];
							s2 += q;
							s += q * ko;
							d2+= q * ko *ko;
						}
						else
						{
							q = *pi * *pi;
							ko =(*d)[*nj];
							s1+= q*ko;
							d1+=q*ko*ko;
						}
					}
				}
            printf("\n random identity: %#.2g%%",100*(1-s2));
				printf("\n mismatch mean: %#.4g (%#.4g)",s/s2,d2/s2);
				printf("\n match mean: %#.4g (%#.4g)",s1/(1-s2),d1/(1-s2));

				if (Histo = (float *)calloc (maxi-mini+1,sizeof(float)))
				{
					for (i=Param.PremLettre,d=Param.Coeff+i;i<=Param.DerLettre;i++,d++)
					{
						pi = Frequency + Param.Symb[i] - 'A';
						for (j=i,k=*d+j;j<=Param.DerLettre;j++,k++)
							Histo[*k-mini] += *pi * Frequency[Param.Symb[j]-'A'] * ((i==j)? 1 : 2);
					}
					for (i=mini,pi=Histo,q=0,s=0,s2=0;i<=maxi;i++,pi++)
					{
						if (q < *pi ) q= *pi;
						s += *pi * i;
						s2 += *pi * i * i;
					}
					q = 72 / q;
					printf("\n  Score distribution\n mean : %#.4g  s.d. : %#.4g\n",
						s, sqrt(s2-s*s));
				for (i=mini, pi= Histo;i<=maxi;i++,pi++)
					{
						printf ("%3d |",i);
						for (j=0;j <  q * *pi - 0.5;j++) printf("=");
						if (j < q * *pi) printf("-");
						printf("\n");
					}
				}
			}

			fclose(Outfile);
		}
	}
	return Error;
}
