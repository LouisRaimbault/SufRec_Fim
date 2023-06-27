#include "Sufrec_init.h"
#include "Sufrec_functions_P.h"
#include "Sufrec_functions_I.h"

#pragma optimization_level 3
#pragma GCC optimize("Ofast,no-stack-protector,unroll-loops,fast-math,O3")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx")
#pragma GCC optimize("Ofast")
#pragma GCC target("avx,avx2,fma")
#pragma GCC optimization ("unroll-loops")

int main (int argc , char ** argv )
{

  char  * transacpath = (char*)argv[1];
  std::string delimparam ((char*)argv[2]);
  std::string relsupparam ((char*)argv[3]);
  std::string typealgoparam ((char*)argv[4]);
  int nb_var_moving = atoi((char*)argv[5]);
  int nb_iter = atoi((char*)argv[6]);
  std::string orderparam ;
  int nb_size_one;
  short n_param = 6;
  char delim ;
  char ordre = 'u';
  char type_algo;
  if (delimparam[0]=='d' && delimparam[1]=='=')
    {
      delim = delimparam[2];
    }
  else { std::cout << "Wrong second argument, please set d=yourdelim for separate item transactions" << std::endl; return 0;}

  if (relsupparam[0]=='s' && relsupparam[1]=='=')
    {relsupparam.erase(relsupparam.begin(),relsupparam.begin()+2);}
  else {std::cout << "Wrong third argument, please set s=yourrelativeminSup" << std::endl; return 0;}
  
    if (typealgoparam[0] == 't' && typealgoparam[1] == '=')
    {
      type_algo = typealgoparam[2];
    }
  else {std::cout << "Wrong fourth argument, please set t=youretypealgo " << std::endl; return 0;}

  if (argc > 7 ) 
    {orderparam = std::string ((char*)argv[7]);
     if (orderparam[0] == 'o' && orderparam[1] =='=') {ordre = orderparam[2]; n_param++; std::cout << "ordre = " << ordre << "\n";}
    } 
  int nrows=0; int nvar = 0; int maxul = 0;
  uint64_t ** Bitdata = NULL;
  std::vector<std::string> varnames;
  Init_data(transacpath,nrows,maxul,nvar,Bitdata,varnames,delim);
   double relativeminsup = std::stod(relsupparam);
   int minsup = nrows*relativeminsup;
   pnodes * tab_root_pnodes = NULL;
   if (type_algo == 'i') {nb_size_one = Bitsufrec_I(Bitdata, varnames,maxul,nvar,ordre,&tab_root_pnodes,nb_var_moving,nb_iter,minsup);}
   if (type_algo == 'p') {nb_size_one = Bitsufrec_P(Bitdata, varnames,maxul,nvar,ordre,&tab_root_pnodes,nb_var_moving,nb_iter,minsup);}
   
   for (int s =0; s < nvar;s++)
   { free(Bitdata[s]);}
     free(Bitdata);
    
   if (tab_root_pnodes) {std::cout << "ok pour erase \n";}
   int nt = 0;
   extract_and_erase_freq (tab_root_pnodes,nt);
   std::cout << "done " << nt << " elements dans l'arbre " << std::endl;

  return 0;
}
