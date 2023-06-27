#ifndef SUFREC_FUNCTIONS_P_H
#define SUFREC_FUNCTIONS_P_H
#include <chrono>
#include "Sufrec_init.h"


void creabitfield_P ( int *ind,int* mod,uint64_t * ptul, uint64_t * nlgtabl,  int nbc);
void creabitfieldr_P ( int *ind,int * mod,uint64_t * ptul, uint64_t * nlgtabl,  int nbdases,  int &reste);
void init_Rnodes_P (rnodes* outrnodes, int & ngroot ,int* vecsom, int * labels ,int* tindices,int * mod,uint64_t** Bitdata, int &nbleft, int& nbcases,  int &reste);
void withwalk_P (pnodes * Tree,pnodes *& curtree, uint64_t * curb,int nbc);
int Bitsufrec_P (uint64_t** Bitdata, std::vector<std::string> &varnames, int maxul, int nvar, char ordre, pnodes ** ret_alpha, int nb_var_moving, int nb_iter, int minsup);


#endif
