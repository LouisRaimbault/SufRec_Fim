#ifndef SUFREC_FUNCTIONS_TO_H
#define SUFREC_FUNCTIONS_TO_H
#include <thread>
#include "Sufrec_init.h"


void init_Rnodes_TO (rnodes* outrnodes, int * ngroot  ,int* vecsom,  int* tindices ,int * mod,uint64_t** Bitdata, int nbleft, int nbcases,  int reste, int * nb_freq_);
void depthwalk_TO (pnodes * curfather, uint64_t ** Bitv, int * tab_ritem_ ,int nbc, int nb_desc, int * nb_freq_);
void thread_witdh2 ( uint64_t ** Bitdata, int * th_ind, int nvar_th, int maxul, int * nb_freq_, pnodes ** tab_subtree_);
void Bitsufrec_TO (uint64_t** Bitdata,int * sum_1freq, info_prog * ip, pnodes *** tab_subtree_);
#endif
