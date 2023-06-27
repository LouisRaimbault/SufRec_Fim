#ifndef SUFREC_FUNCTIONS_I_H
#define SUFREC_FUNCTIONS_I_H
#include "Sufrec_init.h"


void creabitfield_I ( int *ind,int* mod,uint64_t * ptul, uint64_t * nlgtabl,  int nbc);
void creabitfieldr_I ( int *ind,int * mod,uint64_t * ptul, uint64_t * nlgtabl,  int nbdases,  int reste);
void init_Rnodes_I (rnodes* outrnodes, int & ngroot ,int* vecsom,  int* tindices ,int * mod,uint64_t** Bitdata, int nbleft, int nbcases,  int reste);
void withwalk_I (pnodes * curfather, uint64_t ** Bitv, int * tab_ritem_ ,int nbc, int nb_desc);
void Bitsufrec_I (uint64_t** Bitdata, info_prog * ip, pnodes *** tab_pnodes);


#endif
