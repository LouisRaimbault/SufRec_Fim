#ifndef SUFREC_FUNCTIONS_I_H
#define SUFREC_FUNCTIONS_I_H
#include <chrono>
#include "Sufrec_init.h"

void creabitfield_I ( int *ind,int* mod,uint64_t * ptul, uint64_t * nlgtabl,  int nbc);
void creabitfieldr_I ( int *ind,int * mod,uint64_t * ptul, uint64_t * nlgtabl,  int nbdases,  int &reste);
void init_Rnodes_I (pnodes * Tree , int * ngroot , int* tindices ,int * mod, int& nbcases,  int &reste);
void withwalk_I (pnodes * curfathers, uint64_t ** Bitv, int * tab_ritem_ ,int nbc, int nb_desc);
int Bitsufrec_I (uint64_t** Bitdata, std::vector<std::string> &varnames, int maxul, int nvar, char ordre, pnodes ** ret_alpha, int nvar_mod, int n_ajout, int minsup);



#endif
