#include "Sufrec_init.h"
#include "Sufrec_functions_P.h"
#include "Sufrec_functions_I.h"
#include "Sufrec_functions_TF.h"
#include "Sufrec_functions_TO.h"
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
  ip->n_param = 5;
  char  * transacpath = (char*)argv[1];
  std::string delimparam ((char*)argv[2]);
  std::string relsupparam ((char*)argv[3]);
  std::string typealgoparam ((char*)argv[4]);
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
          t_a.erase(t_a.begin(),t_a.begin()+4);
          ip->nb_thread = std::stoi(t_a);
         }
    }
  else {std::cout << "Wrong fourth argument, please set t=youretypealgo " << std::endl; return 0;}
  if (argc > 5 ) 
    {orderparam = std::string ((char*)argv[5]);
     if (orderparam[0] == 'o' && orderparam[1] =='=') {ip->ordre = orderparam[2]; ip->n_param = ip->n_param+1; }
    } 
  ip->rel_minsup = std::stod(relsupparam);
  uint64_t ** Bitdata = NULL;
  int * sum_1freq = NULL;
  Init_data(transacpath,&Bitdata,&sum_1freq,ip);
  if (ip->nvar == 0) { std::cout << "Not a single frequent itemSet, please try with a lower relative minSup value. " << std::endl; return 0;}
 
 ip->minsup = ip->nrows*ip->rel_minsup;
 pnodes ** tab_root_pnodes = NULL;
 pnodes * root_pnodes = NULL;
 if (ip->type_algo == 'p') {Bitsufrec_P(Bitdata,sum_1freq,ip,root_pnodes);}
 if (ip->type_algo == 'i') {Bitsufrec_I(Bitdata,sum_1freq,ip,&tab_root_pnodes);}
 if (ip->type_algo == 't')
    {
     if (typealgoparam[3] == 'f') { Bitsufrec_TF(Bitdata,sum_1freq,ip,&tab_root_pnodes);}
     if (typealgoparam[3] == 'o') { Bitsufrec_TO(Bitdata,sum_1freq,ip,&tab_root_pnodes);}
     if (typealgoparam[3] == 'c') { Bitsufrec_TC(Bitdata,sum_1freq,ip);}
    } 
 free(sum_1freq);
 for (int s =0; s < ip->nvar;s++) {free(Bitdata[s]);} free(Bitdata);
  
  

 if (ip->n_param == argc)
  {

    std::cout << "Erase Frequent set .. \n " ;
    if (ip->type_algo == 't' && typealgoparam[3] == 'c') 
      {
        std::cout << "only count, no itemsets to erase or extract \n";
        delete [] ip->varnames;
        free(ip);
        return 0;

     }
    if (ip->type_algo != 'p') 
      { 
        for (int o = 0; o < ip->nvar;o++) {erase_freq(tab_root_pnodes[o]);}
        free (tab_root_pnodes); std::cout << "done" << std::endl; 
      }
    else {if (ip->type_algo == 'p') {erase_freq (root_pnodes); std::cout << "done \n";}}

    delete [] ip->varnames;
    free(ip);
    return 0;

  }

  
  
  
  if (argc == ip->n_param+1)
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
      for (int o = 0; o < ip->nvar; o++)
        {extract_and_erase_freq (tab_root_pnodes[o],ststring,ip->varnames,nrs,size,flux_set);}
      free(tab_root_pnodes);
      std::cout << "done \n";
      flux_set.close();
      delete [] ip->varnames;
      free(ip);
      return 0;
    }

  
  double nrs = (double)ip->nrows;  
  char  * argoutpath = (char*)argv[ip->n_param];
  ip->n_param = ip->n_param+1;
  char  * argoutpathrules = (char*)argv[ip->n_param];
  std::string outpath (argoutpath);
  std::string outpathrules (argoutpathrules);
  outpath += ".txt";
  outpathrules += "_variables.txt";
  std::ofstream flux_set;
  std::ofstream flux_rules;
  flux_rules.open(outpathrules);
  std::cout << "Writing frequent 1-item name in " << outpathrules << "\n";
  for(int j = 0; j < ip->nvar;j++){ flux_rules << ip->varnames[j] <<"\n";}
  flux_rules.close();

  outpathrules = std::string(argoutpathrules);
  outpathrules += ".txt";
  std::cout << "Filling Frequent sets informations in " << outpath <<  " informations for PrefRules in " << outpathrules <<" and delete Tree ... \n";
  flux_set.open(outpath);
  flux_rules.open(outpathrules);
  flux_set << "nameset\tSupport\trelative_support\tsize\n";
  int size = 1;
  std::string ststring = "";
  for (int o = 0; o < ip->nvar; o++)
    {extract_and_erase_freq (tab_root_pnodes[0],ststring,ip->varnames,nrs,size,flux_set,flux_rules);}
  free(tab_root_pnodes);
  std::cout << "done \n";
  flux_set.close();
  flux_rules.close();
  delete ip->varnames;
  free(ip);

  return 0;
}

