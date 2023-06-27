#ifndef SUFREC_INIT_H
#define SUFREC_INIT_H


#include <fstream>
#include <string>
#include <iostream>
#include <unordered_map>
#include <chrono>

extern const size_t SUL ;
extern const size_t SI  ;
extern const size_t SPN ;
extern const size_t SRN ;

struct info_prog
{
  std::string * varnames;
  int * ind_mach;
  double rel_minsup;  
  double time;
  int nb_freq;
  int minsup;
  int nvar;
  int nrows; 
  int maxul;
  char type_algo;
  int nb_thread;
  int n_param;
  char delim;
  char ordre;
  int nvar_mach;
  int nb_mach;
  int num_mach;

};

struct pnodes 
{
int sup ;
int ritem ;
pnodes * brother;
pnodes * son ;
};

struct rnodes
{
    int sup ;
    uint64_t * tab;
};

void Init_data (char * pathfile, uint64_t *** Bitdata, int ** sum_1freq, info_prog * ip ) ;
void sortoutdata (uint64_t** Bdata, int * sumvec, int nb_elem, std::string * varnames, char ordre);
void init_sufixtree (uint64_t ** Bitdata,int * sum_1freq, std::string * varnames, int support, int nvar, int maxul, char ordre ) ;
void set_ind_machine (info_prog * ip);
void erase_freq (pnodes* Tree) ;
void erase_freq_spec (pnodes* Tree, int * nb_erase, int size);
void extract_and_erase_freq (pnodes* Tree, std::string str, std::string *  listenom, double nrs , int size ,std::ofstream  & flux_set );
void extract_and_erase_freq (pnodes* Tree, std::string str, std::string *  listenom, double nrs, int size, std::ofstream  & flux_set, std::ofstream  & flux_rules );
#endif