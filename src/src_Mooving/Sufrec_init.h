#ifndef SUFREC_INIT_H
#define SUFREC_INIT_H


#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <unordered_map>

extern const size_t SUL ;
extern const size_t SI  ;
extern const size_t SPN ;
extern const size_t SRN ;


struct rnodes;

struct pnodes 
{
int sup ;
rnodes * link_rn;
pnodes * brother;
pnodes * son ;
};

struct rnodes
{   int sup_data;
    int sup ;
    uint64_t * tab;
    rnodes * brother;
    pnodes * cor_bsb; // 
    uint64_t * tab_data;
    std::string *  var;
};

void Init_data (char * pathfile, int & nrows, int & maxul ,int & nvar, uint64_t **& Bitdata, std::vector<std::string> & varnames, char deli ) ;
void init_Rnodes (rnodes* outrnodes, int* vecsom,  int* tindices ,int * mod,uint64_t** Bitdata, int &nbleft, int& nbcases, int Mins,  int &reste);
void sortoutdata (uint64_t** Bdata, int * sumvec, int nb_elem, std::vector<std::string> &varnames, char ordre);
int init_prefixtree (uint64_t ** Bitdata, uint64_t **& Bdata,std::vector<std::string>& varnames, int support, int nvar, int maxul, char ordre );
void extract_and_erase_freq (pnodes*& Tree, int & nt);
void extract_and_erase_freq (pnodes*& Tree, std::string str, std::string *  listenom, std::string * namevalue, int * supvalue, int* sizevalue, int* ritemvalue , int size ,int& situeur );
void set_rnodes_in_order (rnodes * root_rnodes, rnodes * n_rnodes);
void erase_condition_father (pnodes * pn, rnodes  * n_rnodes, int * nbdel);
void erase_mobile (pnodes * pn, int * nbdel);
void erase_var (rnodes * rn, int * nbdel);

#endif
