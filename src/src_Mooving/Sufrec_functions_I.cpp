#include "Sufrec_functions_I.h"


#pragma optimization_level 3
#pragma GCC optimize("Ofast,no-stack-protector,unroll-loops,fast-math,O3")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx")
#pragma GCC optimize("Ofast")
#pragma GCC target("avx,avx2,fma")
#pragma GCC optimization ("unroll-loops")

int nb_freq_I ;
int minsup_I;
 static uint64_t Ultab_I [65];
 uint64_t * PT_Ultab_I = &Ultab_I[0];




void creabitfield_I ( int *ind,int* mod,uint64_t * ptul, uint64_t * nlgtabl,  int nbc)
{   auto * yb = &Ultab_I[0]; uint64_t u; int t = 0; auto * ybi = &Ultab_I[0];
    auto * yd = &Ultab_I[64];
    for (auto a =0; a < nbc ; a++)
    {  
      for (yb = ybi; yb != yd; yb++)
      { 
        u ^= (-((ptul[ind[t]]>>mod[t++]) & 1UL) ^ u) & *yb;
      }
      nlgtabl[a] = u;
    }
}

void creabitfieldr_I ( int *ind,int * mod,uint64_t * ptul, uint64_t * nlgtabl,  int nbdases,  int &reste)
{   
    auto * yb = &Ultab_I[0]; uint64_t u; int t = 0;
    auto * yd = &Ultab_I[64]; auto * ybi = &Ultab_I[0];
    for (auto a =0; a < nbdases ; a++)
    {  
      for (yb = ybi; yb != yd; yb++)
      { 
        u ^= (-((ptul[ind[t]]>>mod[t++]) & 1UL) ^ u) & *yb;
      }
      nlgtabl[a] = u;
    }
    u = 0;
    for ( yb = ybi; yb != &Ultab_I[reste]; yb++) 
    {
     u ^= (-((ptul[ind[t]]>>mod[t++]) & 1UL) ^ u) & *yb;}
    nlgtabl[nbdases] = u;

}



void init_Rnodes_I (pnodes * Tree , int * ngroot , int* tindices ,int * mod, int& nbcases,  int &reste)
{ 
  int sit = 0;
  int sum = 0;
  Tree = Tree->brother; 
  if (reste ==0)
    { 
      while (Tree)
        { sum = Tree->link_rn->sup;         
          if (sum<minsup_I) {Tree = Tree->brother; continue;}
          (*ngroot)++;
          uint64_t * ptu = (uint64_t*) malloc(nbcases*SUL);
          creabitfield_I(tindices,mod,Tree->link_rn->tab_data,ptu,nbcases);
          Tree->link_rn->sup = sum;
          Tree->link_rn->tab = ptu;
          Tree = Tree->brother;
        }
     nb_freq_I += (*ngroot);   
    }

  else {
        
        while (Tree)
          { sum = Tree->link_rn->sup;
            if (sum<minsup_I) {Tree = Tree->brother; continue;}
            (*ngroot)++;
            uint64_t * ptu = (uint64_t*) malloc(nbcases*SUL);
            creabitfieldr_I(tindices,mod,Tree->link_rn->tab_data,ptu,nbcases-1,reste);
            Tree->link_rn->sup = sum;
            Tree->link_rn->tab = ptu;
            Tree = Tree->brother;
          }
        nb_freq_I += (*ngroot); 
      }
}



void withwalk_I (pnodes * curfathers, uint64_t ** Bitv, rnodes ** tab_ritem_ ,int nbc, int nb_desc)
{ 
  pnodes * trpn = curfathers->son->brother;
  
  pnodes ** move_trpn = NULL; int supcand = 0; rnodes* ritem_ ;
  int i=1; int j; int tr_desc = 0; 
  auto * tib = &Bitv[i][0];
  auto * ed = &Bitv[i][nbc];
  auto * eb = &Bitv[i][0];
  
  for (i; i < nb_desc; i++)
    { tib = &Bitv[i][0];
      ed = &Bitv[i][nbc];
      move_trpn = &trpn->son; tr_desc = 0;
      uint64_t ** trt_bns = (uint64_t **) malloc (i*sizeof(uint64_t*));
      rnodes ** tr_ritem_ = (rnodes **) malloc (i * sizeof(rnodes*));
      
      for (j = 0; j < i; j++)          
          { 
            supcand = 0;
            eb = &Bitv[j][0];
            for (auto * tob = tib; tob != ed; tob++, eb++) {supcand += __builtin_popcountl(*tob&*eb);}
            if (supcand < minsup_I) continue; 
            uint64_t * ptab = (uint64_t*) malloc (nbc*SUL);
            auto * pptab = &ptab[0];
            eb = &Bitv[j][0];
            for (auto * tob = tib; tob != ed; tob++, eb++, pptab++) {*pptab = *tob&*eb;}
            trt_bns[tr_desc] = ptab;
            ritem_ = tab_ritem_[j];
            tr_ritem_[tr_desc++] = ritem_;
            *move_trpn = (pnodes*)malloc(SPN);
            (*move_trpn)->sup=supcand;
            (*move_trpn)->link_rn = ritem_;
            (*move_trpn)->brother = NULL;
            (*move_trpn)->son=NULL;                
             move_trpn = &(*move_trpn)->brother;               
            }
        
      nb_freq_I += tr_desc;    
      if (tr_desc > 1 ) {withwalk_I (trpn, trt_bns, tr_ritem_ ,nbc, tr_desc);}
      for (j = 0; j < tr_desc; j++ ) { free(trt_bns[j]);} free(trt_bns); free(tr_ritem_);
        
      trpn = trpn->brother;  
      
    }
}



int Bitsufrec_I (uint64_t** Bitdata, std::vector<std::string> &varnames, int maxul, int nvar, char ordre, pnodes ** ret_alpha, int nvar_mod, int n_ajout, int minsup)
{
 nb_freq_I = 0;
 minsup_I = minsup;
 for (int i = 0; i < 64; i++) {PT_Ultab_I[i] = (1UL <<i);}
 std::cout << "Support value is " << minsup_I << ". Checking unfrequent items ... \n ";
 int sum = 0;
 int nbc = 0;
 auto a = 0; int b= 0;
 auto e = 0;
 auto n = 0;
 uint64_t ** Bdata = NULL ;
 auto * pul = &Bitdata[0][0];
 auto * pulb = &Bitdata[0][0];
 int keep = init_prefixtree(Bitdata, Bdata ,varnames, minsup_I,nvar,maxul,ordre);
 pnodes * last_tree = NULL;
 std::cout <<"Done The Supportvalue is " << minsup_I << ". Start PrefRec with  " << keep << " frequents variables ..." <<  std::endl;
 if(keep ==0)
  {
    std::cout << "Not a single frequent itemSet, please try with a lower relative minsup_I value. " << std::endl;
    return 0;
  }
 
 (*ret_alpha) = (pnodes*)malloc(sizeof(pnodes));
 pnodes * rootalpha = (*ret_alpha);
 rootalpha->brother = NULL;
 rootalpha->son = NULL;
 pnodes * lastree = rootalpha;
 pnodes * ttree = NULL;


 rnodes * rnodes_alpha = (rnodes*)malloc(sizeof(rnodes));
 rnodes_alpha->sup_data = 0;
 rnodes_alpha->brother = (rnodes*)malloc(sizeof(rnodes));
 rnodes_alpha->brother->sup_data = maxul*64 +100;
 rnodes_alpha->brother->brother = NULL;

 rnodes * trn = NULL;

 int t, ngroot;
 std::cout << "nvar moving = " << nvar_mod << "\n";
 uint64_t u;
 for ( t = 0; t < nvar_mod; t++)
  {
    
    nb_freq_I++;   
    pul = Bdata[t];
    sum = 0;
    e =0;
    for (a = 0; a < maxul;a++)
      {sum+= __builtin_popcountl(pul[a]);}
 
    int indices [sum];
    int mod [sum];
    for (a=0; a < maxul; a++)
      { u = pul[a];
        for (n = 0; n < 64; n++)
          { if ((u >> n) & 1UL) {indices[e]= a; mod[e++]=n;} }
      }


    rnodes * root_ritem [t+1];    
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
    init_Rnodes_I(rootalpha,&ngroot,indices,mod,nbc,reste);
    rnodes * rnt = (rnodes*)malloc(sizeof(rnodes));
    rnt->sup_data = sum;
    rnt->tab = NULL;
    rnt->brother = NULL;
    rnt->tab_data = Bdata[t];
    rnt->var = &varnames[t];
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
            {*trpnr = (pnodes*)malloc(SPN);
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
 
    uint64_t ** root_bns = (uint64_t**)malloc(ngroot*sizeof(uint64_t*));
    e = 0;
    while (ttree)
      {
        if (ttree->link_rn->tab) 
        {
          root_bns [e] = ttree->link_rn->tab;
          root_ritem[e++] = ttree->link_rn;
         *trpnr = (pnodes*)malloc(SPN);
         (*trpnr)->brother=NULL;
         (*trpnr)->son=NULL;
         (*trpnr)->sup = ttree->sup;
         (*trpnr)->link_rn = ttree->link_rn;
         trpnr = &(*trpnr)->brother ;          
        }
        ttree = ttree->brother;
      }
      
    withwalk_I (nouvnod,root_bns,root_ritem,nbc, ngroot); 
    free(root_bns); 
     
    ttree = rootalpha->brother;
    while (ttree)
     {if (ttree->link_rn->tab) {free(ttree->link_rn->tab); ttree->link_rn->tab = NULL;}
      ttree = ttree->brother;
      }
         
    lastree->brother = nouvnod;
    lastree = nouvnod;

    }

  
  int nb_del = 0;
  int maxiter = t+n_ajout;
  std::cout << " End first step , number of item in Tree : " << t << ". Number of frequent : " << nb_freq_I << "\n";
   nb_freq_I = 0;
  auto start = std::chrono::system_clock::now();

 for (t ; t < maxiter; t++)
  {

    nb_freq_I++;   
    pul = Bdata[t];
    sum = 0;
    e =0;
    for (a = 0; a < maxul;a++)
      {sum+= __builtin_popcountl(pul[a]);}
 
    int indices [sum];
    int mod [sum];
    for (a=0; a < maxul; a++)
      { u = pul[a];
        for (n = 0; n < 64; n++)
          { if ((u >> n) & 1UL) {indices[e]= a; mod[e++]=n;} }
      }


    rnodes * root_ritem [t+1];    
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
    init_Rnodes_I(rootalpha,&ngroot,indices,mod,nbc,reste);
    rnodes * rnt = (rnodes*)malloc(sizeof(rnodes));
    rnt->sup_data = sum;
    rnt->tab = NULL;
    rnt->brother = NULL;
    rnt->tab_data = Bdata[t];
    rnt->var = &varnames[t];
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
       { lastree->brother = nouvnod;
         lastree = nouvnod;
         erase_var (todel, &nb_del);
         rnodes_alpha->brother = todel->brother;
         free (todel);
         continue;
       }
      nb_del++; 
      free(nouvnod);
      rnodes_alpha->brother = todel->brother;
      free (todel);
      continue; 

    }

    if (ngroot ==1) 
      {
          while (ttree)
          {
            if (ttree->link_rn->tab) 
            {*trpnr = (pnodes*)malloc(SPN);
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
      free(nouvnod);
      nb_del++;
      rnodes_alpha->brother = todel->brother;
      free (todel);
      continue; 
      }
 
    uint64_t ** root_bns = (uint64_t**)malloc(ngroot*sizeof(uint64_t*));
    e = 0;
    while (ttree)
      {
        if (ttree->link_rn->tab) 
        {
          root_bns [e] = ttree->link_rn->tab;
          root_ritem[e++] = ttree->link_rn;
         *trpnr = (pnodes*)malloc(SPN);
         (*trpnr)->brother=NULL;
         (*trpnr)->son=NULL;
         (*trpnr)->sup = ttree->sup;
         (*trpnr)->link_rn = ttree->link_rn;
         trpnr = &(*trpnr)->brother ;          
        }
        ttree = ttree->brother;
      }

    withwalk_I (nouvnod,root_bns,root_ritem,nbc, ngroot);   
    free(root_bns);      
    
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

    nb_del++;
    rnodes_alpha->brother = todel->brother;
    free (todel);
    if (nouvnod->son) {erase_mobile (nouvnod->son,&nb_del);}
    free(nouvnod);         
  }   


  
 std::chrono::duration<double> diffe = std::chrono::system_clock::now() - start;
 std::cout << "End PrefRec,   " << nb_freq_I <<" has been added in " <<diffe.count() << " secs"<< " and " << nb_del << " has been deleted \n";
  free (Bdata);
 rnodes * rn_dell = rnodes_alpha;
 while (rn_dell)
  {rnodes_alpha = rn_dell->brother; free(rn_dell); rn_dell = rnodes_alpha;}

 return keep; 
    
}