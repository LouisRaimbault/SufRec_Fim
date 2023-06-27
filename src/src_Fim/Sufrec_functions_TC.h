#ifndef SUFREC_FUNCTIONS_TC_H
#define SUFREC_FUNCTIONS_TC_H
#include <thread>
#include <mutex>
#include "Sufrec_init.h"

void creabitfield_TC ( int *ind,int* mod,uint64_t * ptul, uint64_t * nlgtabl,  int nbc);
void creabitfieldr_TC ( int *ind,int * mod,uint64_t * ptul, uint64_t * nlgtabl,  int nbdases,  int reste);
void init_Rnodes_TC (rnodes* outrnodes, int & ngroot  ,int* vecsom,  int* tindices ,int * mod,uint64_t** Bitdata, int nbleft, int nbcases,  int reste, int * nb_freq_);
void widthwalk_TC (uint64_t ** Bitv, int * tab_ritem_ ,int nbc, int nb_desc, int * nb_freq_);
void safety_index_C (int * shared_int, int * local_int);
void thread_witdh ( uint64_t ** Bitdata, int maxul, int * nb_freq_, int *index , int nmax );
void Bitsufrec_TC (uint64_t** Bitdata,int * sum_1freq, info_prog * ip);

#endif
