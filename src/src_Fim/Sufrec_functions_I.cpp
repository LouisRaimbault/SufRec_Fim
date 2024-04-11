#include "Sufrec_functions_I.h"


#pragma optimization_level 3
#pragma GCC optimize("Ofast,no-stack-protector,unroll-loops,fast-math,O3")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx")
#pragma GCC optimize("Ofast")
#pragma GCC target("avx,avx2,fma")
#pragma GCC optimization ("unroll-loops")

int nb_freq_I;
int minsup_I;
static uint64_t Ultab_I [65];
uint64_t * PT_Ultab_I = &Ultab_I[0];



void init_Rnodes_I (rnodes* outrnodes, int * ngroot ,int* vecsom,  int* ind ,int * mod,uint64_t** Bitdata, int nbleft, int nbcases,  int reste)
{ 
  auto * yb = &Ultab_I[0]; uint64_t u; int t = 0; int j,nbc;
  auto * yd = &Ultab_I[64]; auto * ybi = &Ultab_I[0]; uint64_t * ptul;
  int sum = 0; 
  if (reste ==0)
    { nbc = nbcases;
      for (auto i=0; i < nbleft;i++)
        { sum = vecsom[i];
          if (sum<minsup_I) {outrnodes[i].tab=NULL; continue;}
          (*ngroot)++;
          uint64_t * ptu = (uint64_t*) malloc(nbc*SUL);
          ptul = Bitdata[i];
          t = 0;
          for (auto a =0; a < nbc ; a++)
          {  u = 0; j = 0;
            for (yb = ybi; yb != yd; yb++,j++)
            { 
              u |= (ptul[ind[t]] >> mod[t++] << j)&*yb;
            }
            ptu[a] = u;
          }
          outrnodes[i].sup = sum;
          outrnodes[i].tab=ptu;
        }
     nb_freq_I += (*ngroot);
     return; 
    }


  nbc = nbcases-1;
  
  for (auto i=0; i < nbleft;i++)
    { sum = vecsom[i];
      
      if (sum<minsup_I) {outrnodes[i].tab=NULL;  continue;}
      (*ngroot)++;
      uint64_t * ptu = (uint64_t*) malloc(nbcases*SUL);
      ptul = Bitdata[i];
      t = 0; 
      for (auto a =0; a < nbc ; a++)
      { j = 0; u = 0;
        for (yb = ybi; yb != yd; yb++,j++)
        { 
          u |= (ptul[ind[t]] >> mod[t++] << j)&*yb;
        }
        ptu[a] = u;
      }
      u = 0; j = 0;
      for ( yb = ybi; yb != &Ultab_I[reste]; yb++,j++) 
      {
       u |= (ptul[ind[t]] >> mod[t++] << j)&*yb;
      }
      ptu[nbc] = u;
      outrnodes[i].sup = sum;
      outrnodes[i].tab=ptu;     
    }
  nb_freq_I += (*ngroot);
      
}


void withwalk_I (pnodes * curfather, uint64_t ** Bitv, int * tab_ritem_ ,int nbc, int nb_desc)
{ 
  
  pnodes * trpn = curfather->son->brother;
  
  pnodes ** move_trpn = NULL; int supcand = 0; int ritem_ ;
  int i=1; int j; int tr_desc = 0; 
  auto * tib = &Bitv[i][0];
  auto * ed = &Bitv[i][nbc];
  auto * eb = &Bitv[i][0];
  
  for (i; i < nb_desc; i++)
    { tib = &Bitv[i][0];
      ed = &Bitv[i][nbc];
      move_trpn = &trpn->son;
      tr_desc = 0;
      uint64_t ** trt_bns = (uint64_t **) malloc (i*sizeof(uint64_t*));
      int * tr_ritem_ = (int *) malloc (i * sizeof(int));
      
      for (j = 0; j < i; j++)          
          { supcand = 0;
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
            (*move_trpn)->ritem = ritem_;
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



void Bitsufrec_I (uint64_t** Bitdata, int * sum_1freq, info_prog * ip, pnodes *** tab_pnodes)
{
 nb_freq_I = 0;
 minsup_I = ip->minsup;
 int maxul = ip->maxul;
 int nvar = ip->nvar;
 for (int i = 0; i < 64; i++) {PT_Ultab_I[i] = (1UL <<i);}

 std::cout << "The Support value is " << minsup_I << ". Checking order parameters and initialize data ... \n ";
 int sum = 0;
 int nbc = 0;
 auto a = 0; int b= 0;
 auto e = 0;
 auto n = 0;
 auto * pul = &Bitdata[0][0];
 auto * pulb = &Bitdata[0][0];
 init_sufixtree(Bitdata, sum_1freq, ip->varnames, minsup_I,nvar,maxul,ip->ordre);
 (*tab_pnodes) = (pnodes**)malloc(nvar*sizeof(pnodes*));
 std::cout <<"Done The Supportvalue is " << minsup_I << ". Start PrefRec with  " << nvar << " frequents variables ..." <<  std::endl;
 auto start = std::chrono::system_clock::now();
 pnodes rootalpha ;
 pnodes * lastree = NULL;
 lastree = &rootalpha; 
 int t, ngroot;
 uint64_t u;
 for ( t = 0; t < nvar; t++)
  {
    nb_freq_I++;   
    pul = Bitdata[t];
    sum = 0;
    e =0;
    for (a = 0; a < maxul;a++)
      {sum+= __builtin_popcountl(pul[a]);}
    rnodes tabrn [t+1]; 
    int  * indices = (int*) malloc (sum*SI);
    int * mod = (int*) malloc (sum*SI);
    for (a=0; a < maxul; a++)
      { u = pul[a];
        for (n = 0; n < 64; n++)
          { if ((u >> n) & 1UL) {indices[e]= a; mod[e++]=n;} }
      }

    
    int rootsum [t+1];
    int root_ritem [t+1];    
    int nbc = sum/64;
    int reste = sum-64*nbc;
    if (reste) nbc ++;
    
    int s = 0;
    for (a = 0; a < t; a++)
        {
          s = 0; pulb = Bitdata[a];
          for (n = 0; n < maxul; n++)
            { s += __builtin_popcountl(pul[n]&pulb[n]);}
          rootsum[a] = s;
        }
    pnodes * nouvnod = (pnodes*)malloc(SPN);
    nouvnod->brother=NULL;
    nouvnod->son=NULL;
    nouvnod->ritem = t;
    nouvnod->sup = sum;
    (*tab_pnodes)[t] = nouvnod ;

    ngroot= 0;
 
 
 init_Rnodes_I(tabrn,&ngroot,rootsum,indices,mod,Bitdata,t,nbc,reste);
 free(indices);
 free(mod);
  
  pnodes ** trpnr = &nouvnod->son;
  if (ngroot <1) {continue;}
  if (ngroot ==1) 
    {
      for (a =0; a < t; a++) 
        {if (tabrn[a].tab) 
          {*trpnr = (pnodes*)malloc(SPN);
           (*trpnr)->brother=NULL;
           (*trpnr)->son = NULL;
           (*trpnr)->ritem = a;
           (*trpnr)->sup=tabrn[a].sup;
           free(tabrn[a].tab); break;  
          }
        }
       continue;
    }
  uint64_t ** root_bns = (uint64_t**)malloc(ngroot*sizeof(uint64_t*));
  e = 0;
  for (a =0; a < t; a++)
    {if (tabrn[a].tab) 
      {root_bns [e] = tabrn[a].tab;
       root_ritem[e++] = a;
       *trpnr = (pnodes*)malloc(SPN);
       (*trpnr)->brother=NULL;
       (*trpnr)->son=NULL;
       (*trpnr)->sup = tabrn[a].sup;
       (*trpnr)->ritem = a;
       trpnr = &(*trpnr)->brother ;          
      }
    }  
  withwalk_I (nouvnod,root_bns,root_ritem,nbc, ngroot);

  free(root_bns);        

  for ( a = 0; a < t; a++)
   {if (tabrn[a].tab) free(tabrn[a].tab);}
  }  
  
 std::chrono::duration<double> diffe = std::chrono::system_clock::now() - start;
 std::cout << "End PrefRec, there is  " << nb_freq_I <<" sets, extracted int  " <<diffe.count() << " secs"<< std::endl;

 return; 
    
}