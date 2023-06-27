#include "Sufrec_init.h"
#include "Sufrec_functions_I.h"
#include "Sufrec_functions_TF.h"
#include "Sufrec_functions_TC.h"

#pragma optimization_level 3
#pragma GCC optimize("Ofast,no-stack-protector,unroll-loops,fast-math,O3")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx")
#pragma GCC optimize("Ofast")
#pragma GCC target("avx,avx2,fma")
#pragma GCC optimization ("unroll-loops")

int main (int argc , char ** argv )
{
  int nb_freq = 0;
  info_prog * ip = (info_prog*)malloc(sizeof(info_prog));
  ip->ordre = 'u';
  ip->n_param = 7;
  char  * transacpath = (char*)argv[1];
  std::string delimparam ((char*)argv[2]);
  std::string relsupparam ((char*)argv[3]);
  std::string typealgoparam ((char*)argv[4]);
  std::string nbmachparam ((char*)argv[5]);
  std::string nummachparam ((char*)argv[6]);
  std::string orderparam ;
  if (delimparam[0]=='d' && delimparam[1]=='=')
    {
      ip->delim = delimparam[2];
    }
  else { std::cout << "Wrong second argument, please set d=yourdelim for separate item transactions \n"; return 0;}

  if (relsupparam[0]=='s' && relsupparam[1]=='=')
    {relsupparam.erase(relsupparam.begin(),relsupparam.begin()+2);}
  else {std::cout << "Wrong third argument, please set s=yourrelativeminSup" << std::endl; return 0;}
  if (typealgoparam[0] == 't' && typealgoparam[1] == '=')
    {
      ip->type_algo = typealgoparam[2];
      if (ip->type_algo == 't') 
        {
          std::string t_a = typealgoparam;
          t_a.erase (t_a.begin(),t_a.begin()+4);
          ip->nb_thread = std::stoi(t_a);
         }
    }
  else {std::cout << "Wrong fourth argument, please set t=youretypealgo " << std::endl; return 0;}
  
  if (nbmachparam[0] == 'm' && nbmachparam[1] == '=')
     {
       nbmachparam.erase(nbmachparam.begin(),nbmachparam.begin()+2);
       ip->nb_mach = std::stoi(nbmachparam);
     }

  if (nummachparam[0] == 'x' && nummachparam[1] == '=')
     {
       nummachparam.erase(nummachparam.begin(),nummachparam.begin()+2);
       ip->num_mach = std::stoi(nummachparam);
     }


  if (argc > 7 ) 
    {orderparam = std::string ((char*)argv[7]);
     if (orderparam[0] == 'o' && orderparam[1] =='=') {ip->ordre = orderparam[2]; ip->n_param = ip->n_param+1; }
    } 
  ip->rel_minsup = std::stod(relsupparam);
  uint64_t ** Bitdata = NULL;
  int * sum_1freq = NULL;
  Init_data(transacpath,&Bitdata,&sum_1freq,ip);
  if (ip->nvar == 0) { std::cout << "Not a single frequent itemSet, please try with a lower relative minSup value. " << std::endl; return 0;} 
 ip->minsup = ip->nrows*ip->rel_minsup;
 init_sufixtree (Bitdata,sum_1freq,ip->varnames,ip->minsup,ip->nvar,ip->maxul,ip->ordre );

 pnodes ** tab_root_pnodes = NULL;
 set_ind_machine (ip);
 if (ip->type_algo == 'i') {Bitsufrec_I(Bitdata,ip,&tab_root_pnodes);}
 if (ip->type_algo == 't')
    {
     if (typealgoparam[3] == 'f') { Bitsufrec_TF(Bitdata,ip,&tab_root_pnodes);}
     if (typealgoparam[3] == 'c') { Bitsufrec_TC(Bitdata,ip);}
    } 
 free(sum_1freq);
 for (int s =0; s < ip->nvar;s++) {free(Bitdata[s]);} free(Bitdata);
  
  


 if (ip->type_algo == 't' && typealgoparam[3] == 'c') 
  {
    std::cout << "only count, no itemsets to erase \n"; 
    delete [] ip->varnames;
    free(ip);
    return 0;
  }

if (ip->n_param == argc)
  {
    std::cout << "Erase Frequent set .. \n " ;
    for (int o = 0; o < ip->nvar_mach;o++) {erase_freq (tab_root_pnodes[ip->ind_mach[o]]);}
    free (tab_root_pnodes); std::cout << "done" << std::endl; 
    delete [] ip->varnames;
    free(ip);
    return 0;
  }



  if (argc > ip->n_param)
    { 
      double nrs = (double)ip->nrows;
      char  * argoutpath = (char*)argv[ip->n_param];
      std::string outpath (argoutpath);
      outpath += ".txt";
      std::cout << "Filling Frequent sets informations in " << outpath << " and delete Tree ... \n";
      std::ofstream flux_set;
      flux_set.open(std::string(outpath));
      flux_set << "nameset\tSupport\trelative_support\tsize\n";
      int size = 1;
      std::string ststring = "";
      for (int o = 0; o < ip->nvar_mach; o++)
        {extract_and_erase_freq (tab_root_pnodes[ip->ind_mach[o]],ststring,ip->varnames,nrs,size,flux_set);}
      free(tab_root_pnodes);
      std::cout << "done \n";
      flux_set.close();
    }

  delete [] ip->varnames;
  free(ip);

  return 0;
}

