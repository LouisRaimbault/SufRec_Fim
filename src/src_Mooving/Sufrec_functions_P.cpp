#include "Sufrec_functions_P.h"


#pragma optimization_level 3
#pragma GCC optimize("Ofast,no-stack-protector,unroll-loops,fast-math,O3")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx")
#pragma GCC optimize("Ofast")
#pragma GCC target("avx,avx2,fma")
#pragma GCC optimization ("unroll-loops")

int nb_freq_P ;
int minsup_P;
 static uint64_t Ultab_P [65];
 uint64_t * PT_Ultab_P = &Ultab_P[0];




void creabitfield_P ( int *ind,int* mod,uint64_t * ptul, uint64_t * nlgtabl,  int nbc)
{   auto * yb = &Ultab_P[0]; uint64_t u; int t = 0; auto * ybi = &Ultab_P[0];
    auto * yd = &Ultab_P[64];
    for (auto a =0; a < nbc ; a++)
    {  
      for (yb = ybi; yb != yd; yb++)
      { 
        u ^= (-((ptul[ind[t]]>>mod[t++]) & 1UL) ^ u) & *yb;
      }
      nlgtabl[a] = u;
    }
}

void creabitfieldr_P ( int *ind,int * mod,uint64_t * ptul, uint64_t * nlgtabl,  int nbdases,  int &reste)
{   
    auto * yb = &Ultab_P[0]; uint64_t u; int t = 0;
    auto * yd = &Ultab_P[64]; auto * ybi = &Ultab_P[0];
    for (auto a =0; a < nbdases ; a++)
    {  
      for (yb = ybi; yb != yd; yb++)
      { 
        u ^= (-((ptul[ind[t]]>>mod[t++]) & 1UL) ^ u) & *yb;
      }
      nlgtabl[a] = u;
    }
    u = 0;
    for ( yb = ybi; yb != &Ultab_P[reste]; yb++) 
    {
     u ^= (-((ptul[ind[t]]>>mod[t++]) & 1UL) ^ u) & *yb;}
    nlgtabl[nbdases] = u;

}



void init_Rnodes_P (pnodes * Tree , int * ngroot , int* tindices ,int * mod, int& nbcases,  int &reste)
{ 
  
  int sum = 0;
  Tree = Tree->brother; 
  if (reste ==0)
    { 
      while (Tree)
        { sum = Tree->link_rn->sup;         
          if (sum<minsup_P) {Tree->link_rn->sup = 0; Tree = Tree->brother;  continue;}
          (*ngroot)++;
          uint64_t * ptu = (uint64_t*) malloc(nbcases*SUL);
          creabitfield_P(tindices,mod,Tree->link_rn->tab_data,ptu,nbcases);
          Tree->link_rn->sup = sum;
          Tree->link_rn->tab = ptu;
          Tree = Tree->brother;
        }
     nb_freq_P += (*ngroot);   
    }

  else {
        
        while (Tree)
          { sum = Tree->link_rn->sup;
            if (sum<minsup_P) {Tree->link_rn->sup = 0; Tree = Tree->brother;  continue;}
            (*ngroot)++;
            uint64_t * ptu = (uint64_t*) malloc(nbcases*SUL);
            creabitfieldr_P(tindices,mod,Tree->link_rn->tab_data,ptu,nbcases-1,reste);
            Tree->link_rn->sup = sum;
            Tree->link_rn->tab = ptu;
            Tree = Tree->brother;
          }
        nb_freq_P += (*ngroot); 
      }
}



void withwalk_P (pnodes * Tree,pnodes *& curtree, uint64_t * curb,int nbc)
{ 
  
  if (!Tree->link_rn->sup)
    {  
       if (Tree->brother) {withwalk_P(Tree->brother,curtree,curb,nbc);}
       return;
    }
  
  auto * tib = &curb[0];
  auto * ed = &curb[nbc];
  auto * cb = &Tree->link_rn->tab[0];
  int sum = 0;


  for (tib; tib != ed; tib++, cb++ )
    {sum+= __builtin_popcountl(*tib&*cb);}
  if (sum < minsup_P) 
    {  
       Tree->link_rn->sup = 0;
       if (Tree->brother) {withwalk_P(Tree->brother,curtree,curb,nbc);}
       Tree->link_rn->sup = 1;
       return;        
    }
    nb_freq_P++;
    curtree = (pnodes*)malloc(SPN);
    curtree->sup = sum;
    curtree->link_rn = Tree->link_rn;
    curtree->son = NULL;
    curtree->brother = NULL;

    if (!Tree->son)
     {   
         if (Tree->brother) {withwalk_P(Tree->brother,curtree->brother,curb,nbc);}
         return;              
     }
     
    uint64_t * ncurb = (uint64_t*)malloc(nbc*SUL);
    auto * ad = &ncurb[0];
    cb = &Tree->link_rn->tab[0];
    for (tib = &curb[0]; tib != ed; tib++, cb++,ad++) {*ad = *tib&*cb;}

    withwalk_P(Tree->son,curtree->son,ncurb,nbc);
    free(ncurb);
    if (Tree->brother) {withwalk_P(Tree->brother,curtree->brother,curb,nbc);}
}



int Bitsufrec_P (uint64_t** Bitdata, std::vector<std::string> &varnames, int maxul, int nvar, char ordre, pnodes ** ret_alpha, int nb_var_moving, int nb_iter, int minsup)
{
 nb_freq_P = 0;
 minsup_P = minsup;
 for (int i = 0; i < 64; i++) {PT_Ultab_P[i] = (1UL <<i);}
 std::cout << "Support value is " << minsup_P << ". Checking unfrequent items ... \n ";
 int sum = 0;
 int nbc = 0;
 auto a = 0; int b= 0;
 auto e = 0;
 auto n = 0;
 uint64_t ** Bdata = NULL ;
 auto * pul = &Bitdata[0][0];
 auto * pulb = &Bitdata[0][0];
 int keep = init_prefixtree(Bitdata, Bdata ,varnames, minsup_P,nvar,maxul,ordre);
 pnodes * last_tree = NULL;
 std::cout <<"Done The Supportvalue is " << minsup_P << ". Start PrefRec with  " << keep << " frequents variables ..." <<  std::endl;
 if(keep ==0)
  {
    std::cout << "Not a single frequent itemSet, please try with a lower relative minsup_P value. " << std::endl;
    return 0;
  }
 
 (*ret_alpha) = (pnodes*)malloc(sizeof(pnodes));
 pnodes * rootalpha = (*ret_alpha);
 rootalpha->brother = NULL;
rootalpha->son = NULL;
 pnodes * lastree = rootalpha;
 pnodes * ttree = NULL;


 rnodes * rnodes_alpha = (rnodes*)malloc(sizeof(rnodes));
 rnodes_alpha->sup = 0;
 rnodes_alpha->sup_data = 0;
 rnodes_alpha->brother = (rnodes*)malloc(sizeof(rnodes));
 rnodes_alpha->brother->sup = maxul*64 +100;
 rnodes_alpha->brother->sup_data = maxul*64 +100;
 rnodes_alpha->brother->brother = NULL;

 rnodes * trn = NULL;

 int j, ngroot;
 uint64_t u; 

 for ( j=0; j<nb_var_moving;j++)
  { 
    nb_freq_P++;   
    pul = Bdata[j];
    sum = 0;
    e =0;
    for (a = 0; a < maxul;a++)
      {sum+= __builtin_popcountl(pul[a]);}
 
    int * indices = (int*)malloc(sum*SI);
    int * mod  = (int*)malloc(sum*SI);
    for (a=0; a < maxul; a++)
      { u = pul[a];
        for (n = 0; n < 64; n++)
          { if ((u >> n) & 1UL) {indices[e]= a; mod[e++]=n;} }
      }


    rnodes * root_ritem [j+1];    
    int nbc = sum/64;
    int reste = sum-64*nbc;
    if (reste) nbc ++;
    
    int s = 0;
    ttree = rootalpha->brother;
    while (ttree)
      {  
       s = 0; pulb = ttree->link_rn->tab_data;
       for (n = 0; n < maxul; n++) { s += __builtin_popcountl(pul[n]&pulb[n]);}
       ttree->link_rn->sup = s;
       ttree = ttree->brother;
      }  

    ngroot= 0;
    init_Rnodes_P(rootalpha,&ngroot,indices,mod,nbc,reste);
    free(indices);
    free(mod);
    rnodes * rnt = (rnodes*)malloc(sizeof(rnodes));
    rnt->sup_data = sum;
    rnt->tab = NULL;
    rnt->brother = NULL;
    rnt->tab_data = Bdata[j];
    rnt->var = &varnames[j];
    rnt->cor_bsb = lastree;

    pnodes * nouvnod = (pnodes*)malloc(SPN);
    nouvnod->brother=NULL;
    nouvnod->son=NULL;
    nouvnod->link_rn = rnt;
    nouvnod->sup = sum;

    set_rnodes_in_order (rnodes_alpha,rnt);

    pnodes ** trpnr = &nouvnod->son;
    ttree = rootalpha->brother;

    if (ngroot <1) 
    { 
      lastree->brother = nouvnod;
      lastree = nouvnod;
      continue;
    }

    if (ngroot ==1) 
      { 

        while (ttree)
        {
            if (ttree->link_rn->tab) 
            {
              *trpnr = (pnodes*)malloc(SPN);
             (*trpnr)->brother=NULL;
             (*trpnr)->son = NULL;
             (*trpnr)->link_rn = ttree->link_rn;
             (*trpnr)->sup=ttree->link_rn->sup;
             free(ttree->link_rn->tab);
             ttree->link_rn->tab = NULL;
             break;  
            }
          ttree = ttree->brother;
        }
       lastree->brother = nouvnod;
       lastree = nouvnod;
       continue;
      } 
    
    e = 0;

    while (ttree)
      {
        if (!ttree->link_rn->tab)
          {
            ttree = ttree->brother;
            continue;
          }

        *trpnr = (pnodes*)malloc(SPN);
        (*trpnr)->sup = ttree->link_rn->sup; 
        (*trpnr)->link_rn = ttree->link_rn ;
        (*trpnr)->son = NULL;
        (*trpnr)->brother=NULL;
        (*trpnr)->link_rn = ttree->link_rn;
        if (ttree->son)
          {
            withwalk_P(ttree->son,(*trpnr)->son,ttree->link_rn->tab,nbc);
          }        

        ttree = ttree->brother;
        trpnr = &(*trpnr)->brother; 
        
      }

    ttree = rootalpha->brother;  
    while (ttree)
     {if (ttree->link_rn->tab) {free(ttree->link_rn->tab); ttree->link_rn->tab = NULL;}
      ttree = ttree->brother;
      }
         
    lastree->brother = nouvnod;
    lastree = nouvnod;
   
  }  
 

  int nb_del = 0;
  int maxiter = j+nb_iter;
   std::cout << " End first step , number of items in Tree :  " << j << ". Number of frequent : " << nb_freq_P << "\n";
   nb_freq_P = 0;
  auto start = std::chrono::system_clock::now();
  for (j; j < maxiter;j++)
  { 
    nb_freq_P++;   
    pul = Bdata[j];
    sum = 0;
    e =0;
    for (a = 0; a < maxul;a++)
      {sum+= __builtin_popcountl(pul[a]);}
 
    int * indices = (int*)malloc(sum*SI);
    int * mod  = (int*)malloc(sum*SI);
    for (a=0; a < maxul; a++)
      { u = pul[a];
        for (n = 0; n < 64; n++)
          { if ((u >> n) & 1UL) {indices[e]= a; mod[e++]=n;} }
      }


    rnodes * root_ritem [j+1];    
    int nbc = sum/64;
    int reste = sum-64*nbc;
    if (reste) nbc ++;
    
    int s = 0;
    ttree = rootalpha->brother;
    while (ttree)
      {  
       s = 0; pulb = ttree->link_rn->tab_data;
       for (n = 0; n < maxul; n++) { s += __builtin_popcountl(pul[n]&pulb[n]);}
       ttree->link_rn->sup = s;
       ttree = ttree->brother;
      }  

    ngroot= 0;
    init_Rnodes_P(rootalpha,&ngroot,indices,mod,nbc,reste);
    free(indices);
    free(mod);
    rnodes * rnt = (rnodes*)malloc(sizeof(rnodes));
    rnt->sup_data = sum;
    rnt->tab = NULL;
    rnt->brother = NULL;
    rnt->tab_data = Bdata[j];
    rnt->var = &varnames[j];
    rnt->cor_bsb = lastree;

    pnodes * nouvnod = (pnodes*)malloc(SPN);
    nouvnod->brother=NULL;
    nouvnod->son=NULL;
    nouvnod->link_rn = rnt;
    nouvnod->sup = sum;

    set_rnodes_in_order (rnodes_alpha,rnt);

    pnodes ** trpnr = &nouvnod->son;
    ttree = rootalpha->brother;

    if (ngroot <1) 
    { 
      rnodes * todel = rnodes_alpha->brother;
      if (todel!= rnt)
        {
          lastree->brother = nouvnod;
          lastree = nouvnod;
          erase_var (todel, &nb_del);
          rnodes_alpha->brother = todel->brother;
          free (todel);
          continue;
        }
      free(nouvnod);
      nb_del++;
      rnodes_alpha->brother = todel->brother;
      free (todel);
      continue; 
    }

    if (ngroot ==1) 
      { 

        while (ttree)
        {
            if (ttree->link_rn->tab) 
            {
              *trpnr = (pnodes*)malloc(SPN);
             (*trpnr)->brother=NULL;
             (*trpnr)->son = NULL;
             (*trpnr)->link_rn = ttree->link_rn;
             (*trpnr)->sup=ttree->link_rn->sup;
             free(ttree->link_rn->tab);
             ttree->link_rn->tab = NULL;
             break;  
            }
          ttree = ttree->brother;
        }

        rnodes * todel = rnodes_alpha->brother;
        
        if (todel!= rnt)
          { lastree->brother = nouvnod;
            lastree = nouvnod;
            erase_var (todel, &nb_del);
            rnodes_alpha->brother = todel->brother;
            free (todel);
             continue;
          }

       if (nouvnod->son) 
        {erase_mobile (nouvnod->son,&nb_del);}
       nb_del++;
       free(nouvnod);
       rnodes_alpha->brother = todel->brother;
       
       free (todel);
       continue; 
      } 
    
    e = 0;

    while (ttree)
      {
        if (!ttree->link_rn->tab)
          {
            ttree = ttree->brother;
            continue;
          }

        *trpnr = (pnodes*)malloc(SPN);
        (*trpnr)->sup = ttree->link_rn->sup; 
        (*trpnr)->link_rn = ttree->link_rn ;
        (*trpnr)->son = NULL;
        (*trpnr)->brother=NULL;
        (*trpnr)->link_rn = ttree->link_rn;
        if (ttree->son)
          {
            withwalk_P(ttree->son,(*trpnr)->son,ttree->link_rn->tab,nbc);
          }        

        ttree = ttree->brother;
        trpnr = &(*trpnr)->brother; 
        
      }
      
    ttree = rootalpha->brother;  
    while (ttree)
     {
      if (ttree->link_rn->tab) {free(ttree->link_rn->tab); ttree->link_rn->tab = NULL;}
      ttree = ttree->brother;
     }
         
    rnodes * todel = rnodes_alpha->brother;
    if (todel!= rnt)
      { lastree->brother = nouvnod;
        lastree = nouvnod;
        erase_var (todel, &nb_del);
        rnodes_alpha->brother = todel->brother;
        free (todel);
         continue;
      }
    rnodes_alpha->brother = todel->brother;
    
    free (todel);
    if (nouvnod->son) {erase_mobile (nouvnod->son,&nb_del);}
    free(nouvnod);
    nb_del++;    
  }


  std::chrono::duration<double> diffe= std::chrono::system_clock::now() - start;
  std::cout << "End PrefRec,   " << nb_freq_P <<" has been added in " <<diffe.count() << " secs"<< " and " << nb_del << " has been deleted \n";
  free (Bdata);
  rnodes * rn_dell = rnodes_alpha;
  while (rn_dell)
  {rnodes_alpha = rn_dell->brother; free(rn_dell); rn_dell = rnodes_alpha;}
return keep; 
    
}