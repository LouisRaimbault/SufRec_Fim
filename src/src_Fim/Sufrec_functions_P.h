#ifndef SUFREC_FUNCTIONS_P_H
#define SUFREC_FUNCTIONS_P_H
#include "Sufrec_init.h"



void init_Rnodes_P (rnodes* outrnodes, int *  ngroot ,int* vecsom,  int* tindices ,int * mod,uint64_t** Bitdata, int nbleft, int nbcases,  int reste);
void withwalk_P (pnodes * Tree,pnodes ** curtree, uint64_t * curb, int cur_r ,rnodes * tab_r ,int nbc);
void Bitsufrec_P (uint64_t** Bitdata, int * sum_1freq, info_prog * ip, pnodes *& root_);

#endif
