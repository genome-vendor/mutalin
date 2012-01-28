/*------------------------------------------------------------------------
	PARAMETR.H
					Parametres d'alignement

-------------------------------------------------------------------------*/

#ifndef __PARAMETR_H

#define __PARAMETR_H
#include "deftypes.h"
void AjouteSymb (char* Ens, int c);

int ExisteSequence( t_pDescripteurSequence DS, t_indN DerSeq );
/* regarde si le nom de Sequence[DerSeq] est aussi celui d'une sequence precedente */

void InitParam( t_pParamAlign Param );
/* initialisation de la structure de type t_paramAlign */

void MajParam(t_pParamAlign Param );
/* maj de Param apres lecture de fichier et de la ligne de commande */

void VerifiSymb(t_pDescripteurSequence DS, t_pParamAlign Param,
	char* EnsLet);

#endif
