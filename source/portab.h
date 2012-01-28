/*--------------------------------------------------------------------------
	Fichier PORTAB.H

		Parametres permettant d'assurer la portabilite
				du programme MA.
---------------------------------------------------------------------------*/
#ifndef PORTAB_H
#define PORTAB_H


/*--------- environnement de compilation -----------------------*/
/* ! Conserver une seule definition parmi les suivantes pour indiquer
l'environnement et le compilateur utilises.
la fonction MasqRech de disk.c fait intervenir la fonction system qui n'est pas
ANSI : vous pouvez avoir besion de la redefinir
la fonction InitMessages de msgerr.c fait intervenir la fonction isatty qui
n'est pas ANSI : vous pouvez avoir besoin de la redefinir */


/* for Borland C  __BORLANDC__ is defined,
	and _Windows is defined with Microsoft Windows */

/*#define _GCC_SUN_*/

/*--------- File name strings ------------*/
/* extension is the string beginning with the last point in a complete file name
	a file name is the string between the last dirchar and the extension */
#ifdef __BORLANDC__
#define dirchar '\\'
#else
#define dirchar '/'
#endif

/*---------  mul and Fasta output format -----------*/
/* you can choose mul as the output format : you will get all the sequences,
printed one after the other with '-' inserted in the gaps; by default, the
sequence position at the beginning of each line is written and the positions are
clustered by 10
>Sequence_name
	  1 AZERTYUIOP QSDF---KLM ...
	 51 WXC...

To get a true Fasta format, you can define NOCOMMENT; you will then get :
>Sequence_name
AZERTYUIOPQSDF---KLM...
WXC...
*/
#define NOCOMMENT


/*--------- shell command to list a directory ------------*/
#ifdef __BORLANDC__
#define CmdeListeFich "dir "
#else
#define CmdeListeFich "ls "
#endif

/*--------- parametres de fonctionnement ----------*/
/* Ces parametres doivent etre choisis apres tests d'occupation memoire
en utilisant des tailles de sequences croissantes.

ATTENTION: le programme ne detecte le manque de memoire qu'apres saturation
de la partition de swap du disque.
	Deux comportements peuvent alors se produire:
		- crash du systeme !!
		- mise em place d'un mode de swap propre au programme, le
repertoire de swap etant le repertoire temporaire par defaut, ou le repertoire
courant de l'utilisateur. Dans les deux cas, une saturation de disque est aussi
possible. */


/* Les valeurs maxi admissibles dependent de la taille maxi allouee par malloc*/

#ifdef __BORLANDC__
/* sizeof(size_t) = 2 */
#define NbSeqMax (0xFFF0 / sizeof(t_LigneSeq) - 1)
#define LongMax (0x7FF8/ sizeof(t_score) -1)
#else
/* sizeof(size_t) > 2 */
#define NbSeqMax (0x7FF0)
#define LongMax (0x7FF0)
#endif

/*----------  type 2 octets  ------------------------------------------*/
/* La plus grande partie des variables utilises dans MA necessitent un codage
sur 2 octets. Pour limiter l'occupation memoire, il est necessaire d'utiliser
le type C correspondant s'il existe:

		|	type_2_octets
	------------------------------------
	DOS/PC	|	short et int
	UNIX	|	short eventuellement 	(depend de la machine)  */

/* index for Sequence # */
typedef short t_indN;

/* index for Sequence position (length) */
typedef short t_indL;


/* caracteres graphiques */
/* if extended ASCII characters are printable, define _ASCII to use them */
#ifdef _ASCII
#	define Trait 'Ä'
#	define DepartBas 'Â'
#	define Verticale '³'
#	define DepartGauche  '´'
#	define BasDroite 'Ù'
#else
#	define Trait '-'
#	define DepartBas '+'
#	define Verticale '|'
#	define DepartGauche  '+'
#	define BasDroite '+'
#endif

/* if you have a gif library you can use the possibility to create gif images
 so define _gif */
/*#define _gif*/ 

/* _clus2dom is defined for a private use of MultAlin */
#ifdef clus2dom
#	undef _gif
#endif
#endif
