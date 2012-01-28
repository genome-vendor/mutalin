/*----------------------------------------------------------------------
				Alignement multiple des sequences
----------------------------------------------------------------------*/


#ifndef __MA_H

#define __MA_H
#include "deftypes.h"
/* dans copright.c */
extern char  CopyRight[5][57];
void ChargeCopyRight (void);
int PrintCopyRight (FILE *f);

/* dans maglob.c */
void InitSequenceDescription (t_pDescripteurSequence DS);
void LibereDescripteurSequence( t_pDescripteurSequence DS );
int MultipleAlignment (t_pDescripteurSequence DS, t_pParamAlign ParamAlign,
	int SeqDejaAligne);


/* dans aligne.c */
int Aligner (t_pDescripteurGroupe DG,t_pDescripteurSequence DS,
			t_ClassHier Arbre,t_pGrdTab GrdTab,t_pParamAlign ParamAlign);
int AligneAction( t_TabClassHier TabArbre,t_DessinArb *Dessin,
		t_pDescripteurSequence DS,t_pParamAlign Param,t_pDescripteurGroupe DG );
int Align2Groups (int N2,t_ClassHier Arbre, t_pDescripteurSequence DS,
	t_pDescripteurGroupe DG, t_pParamAlign Param);

/* dans fast.c */
int FastAction2( t_pDescripteurSequence DS, t_pParamAlign Param,
			t_ClassHier TabArbre, t_DessinArb Dessin  );

/* dans cluster.c */
int GetArbre (t_pDescripteurSequence DS, t_ClassHier Arbre, t_DessinArb Dessin);
/* try to get an already calculated tree, return 1 if success, 0 if failed
and -1 if not available */
void SauvArbre (t_ClassHier Arbre, t_DessinArb Dessin,t_pDescripteurSequence DS,
	char *Extension);

/* dans drivers.c */
extern char CtrlBreakHit;
void InitBreakControl (void);
void DoneBreakControl (void);


#endif
