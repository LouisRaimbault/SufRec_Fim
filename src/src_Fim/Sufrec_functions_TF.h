#ifndef SUFREC_FUNCTIONS_TF_H
#define SUFREC_FUNCTIONS_TF_H
#include <thread>
#include <mutex>
#include "Sufrec_init.h"


void init_Rnodes_TF (rnodes* outrnodes, int * ngroot  ,int* vecsom,  int* tindices ,int * mod,uint64_t** Bitdata, int nbleft, int nbcases,  int reste, int  * nb_freq_);
void depthwalk_TF (pnodes * curfather, uint64_t ** Bitv, int * tab_ritem_ ,int nbc, int nb_desc, int * nb_freq_);
void safety_index (int * shared_int, int * local_int);
void thread_witdh ( uint64_t ** Bitdata, int maxul, int * nb_freq_, pnodes ** tab_subtree_, int *index , int nmax );
void Bitsufrec_TF (uint64_t** Bitdata,int * sum_1freq, info_prog * ip, pnodes *** tab_subtree_);

#endif
