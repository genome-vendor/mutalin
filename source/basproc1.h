/*------------------------------------------------------------------------
	fichier BASPROC1.H

				fonctions utilisees pour l'alignement

-------------------------------------------------------------------------*/
/*
CORRECTIONS

1 juillet 1994 : MajTailleConsens has a third parameter "complet" that must be
	false when the sequences are not yet aligned

OPTIMISATIONS

7 juillet 1994 : utilisation de TabArbre[0] pour AncArbre


*/


#ifndef __BASPROC1_H

#define __BASPROC1_H
#include "deftypes.h"
#define LargeurGraphe 110
#define DebutGraphe 10
extern int TailleGraphe;

typedef int (*ProcAction) (t_indN, t_indN, t_indN, t_indN*, t_ClassHier, void*);

extern int ArbreAlignement(t_pParamAlign Param,t_pDescripteurSequence DS,
				t_ClassHier Arbre,t_DessinArb *Dessin,t_indN NbSeq);

extern int CalculArbre( t_score** Score, 		/* scores puis resultats */
			t_ClassHier Arbre,	/* nouvel arbre */
			t_DessinArb Dessin,
			t_indN NbSeq,
			t_score Maxi );

extern int EditSizeLine(t_pLigneSeq p, t_indL NlleTaille );

extern void Id2a2  (t_pParamAlign Param, t_pLigneSeq Sequence, t_score *x ,
			t_indN Prot, t_indN Prot2 );

extern void MajTailleConsens(t_pLigneSeq Sequence, t_indN NbSeq, char complet );

extern void MajWeightConsens (t_pLigneSeq Sequence, t_indN NbSeq);

extern t_score *Place(t_indN i,t_indN j,t_score**Score);

extern void Sc2a2( t_pParamAlign Param, t_pLigneSeq Sequence, t_score *x ,
			t_indN Prot, t_indN Prot2 );

extern void SupprimeGroupe( t_ClassHier Arbre, t_indN* Groupe, t_indN Groupe1,
			t_indN Groupe2 );

extern void TraducByteToLet( t_pDescripteurSequence DS, t_pParamAlign Param );

extern void TraducLetToByte( t_pDescripteurSequence DS, t_pParamAlign Param );

extern int TraiteBranche(t_ClassHier Arbre, t_indN NbSeq, t_indN *Groupe, void *V,
	 ProcAction Action);

#endif
