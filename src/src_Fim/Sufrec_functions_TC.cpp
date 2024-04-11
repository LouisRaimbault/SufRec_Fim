#include "Sufrec_functions_TC.h"
#pragma optimization_level 3
#pragma GCC optimize("Ofast,no-stack-protector,unroll-loops,fast-math,O3")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx")
#pragma GCC optimize("Ofast")
#pragma GCC target("avx,avx2,fma")
#pragma GCC optimization ("unroll-loops")

 int minsup_TC;
 static uint64_t Ultab_TC [65];
 uint64_t * PT_Ultab_TC = &Ultab_TC[0];
 std::mutex mtx_C;



void init_Rnodes_TC (rnodes* outrnodes, int * ngroot ,int* vecsom,  int* ind ,int * mod,uint64_t** Bitdata, int nbleft, int nbcases,  int reste, int * nb_freq_)
{ 
  auto * yb = &Ultab_TC[0]; uint64_t u; int t = 0; int j,nbc;
  auto * yd = &Ultab_TC[64]; auto * ybi = &Ultab_TC[0]; uint64_t * ptul;
  int sum = 0; 
  if (reste ==0)
    { nbc = nbcases;
      for (auto i=0; i < nbleft;i++)
        { sum = vecsom[i];
          if (sum<minsup_TC) {outrnodes[i].tab=NULL; continue;}
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
     (*nb_freq_) += (*ngroot);
     return; 
    }


  nbc = nbcases-1;
  
  for (auto i=0; i < nbleft;i++)
    { sum = vecsom[i];
      
      if (sum<minsup_TC) {outrnodes[i].tab=NULL;  continue;}
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
      for ( yb = ybi; yb != &Ultab_TC[reste]; yb++,j++) 
      {
       u |= (ptul[ind[t]] >> mod[t++] << j)&*yb;
      }
      ptu[nbc] = u;
      outrnodes[i].sup = sum;
      outrnodes[i].tab=ptu;     
    }
  (*nb_freq_) += (*ngroot); 
      
}



void widthwalk_TC (uint64_t ** Bitv, int * tab_ritem_ ,int nbc, int nb_desc, int * nb_freq_)
{ 

  int supcand = 0; int ritem_ ;
  int i=1; int j; int tr_desc = 0; 
  auto * tib = &Bitv[i][0];
  auto * ed = &Bitv[i][nbc];
  auto * eb = &Bitv[i][0];
  
  for (i; i < nb_desc; i++)
    { tib = &Bitv[i][0];
      ed = &Bitv[i][nbc];
      tr_desc = 0;
      uint64_t ** trt_bns = (uint64_t **) malloc (i*sizeof(uint64_t*));
      int * tr_ritem_ = (int *) malloc (i * sizeof(int));
      
      for (j = 0; j < i; j++)          
          { supcand = 0;
            eb = &Bitv[j][0];
            for (auto * tob = tib; tob != ed; tob++, eb++) {supcand += __builtin_popcountl(*tob&*eb);}
            if (supcand < minsup_TC) continue; 
            uint64_t * ptab = (uint64_t*) malloc (nbc*SUL);
            auto * pptab = &ptab[0];
            eb = &Bitv[j][0];
            for (auto * tob = tib; tob != ed; tob++, eb++, pptab++) {*pptab = *tob&*eb;}
            trt_bns[tr_desc] = ptab;
            ritem_ = tab_ritem_[j];
            tr_ritem_[tr_desc++] = ritem_;                
              
            }
        
      (*nb_freq_) += tr_desc;    
      if (tr_desc > 1 ) {widthwalk_TC (trt_bns, tr_ritem_ ,nbc, tr_desc,nb_freq_);}
      for (j = 0; j < tr_desc; j++ ) { free(trt_bns[j]);} free(trt_bns); free(tr_ritem_); 
      
    }
}

void safety_index_C (int * shared_int, int * local_int)
{ 
  *local_int = *shared_int;
  *shared_int = *shared_int+1;
}

void thread_witdh ( uint64_t ** Bitdata, int maxul, int * nb_freq_, int *index , int nmax )
{  
  int t, sum, e , ngroot;
  nmax = nmax -1;
  uint64_t * pul = NULL;
  uint64_t * pulb = NULL;
  auto a =0;
  uint64_t u;
  int n = 0;
  int cur_v;
  t = 1;
    while (t)
  { 
    mtx_C.lock();
    safety_index_C(index,&cur_v);
    mtx_C.unlock();
    if (cur_v > nmax) {t=0; continue;}
    
    (*nb_freq_)++;   
    pul = Bitdata[cur_v];
    sum = 0;
    e =0;
    for (a = 0; a < maxul;a++)
      {sum+= __builtin_popcountl(pul[a]);}
    rnodes tabrn [cur_v+1]; 
    int  * indices = (int*) malloc (sum*SI);
    int * mod = (int*) malloc (sum*SI);
    for (a=0; a < maxul; a++)
      { u = pul[a];
        for (n = 0; n < 64; n++)
          { if ((u >> n) & 1UL) {indices[e]= a; mod[e++]=n;} }
      }

    
    int rootsum [cur_v+1];
    int root_ritem [cur_v+1];    
    int nbc = sum/64;
    int reste = sum-64*nbc;
    if (reste) nbc ++;
    
    int s = 0;
    for (a = 0; a < cur_v; a++)
        {
          s = 0; pulb = Bitdata[a];
          for (n = 0; n < maxul; n++)
            { s += __builtin_popcountl(pul[n]&pulb[n]);}
          rootsum[a] = s;
        }

    ngroot= 0;
 
 
 init_Rnodes_TC(tabrn,&ngroot,rootsum,indices,mod,Bitdata,cur_v,nbc,reste,nb_freq_);
 free(indices);
 free(mod);

  if (ngroot <1) { continue;}
  if (ngroot ==1) 
    {for (a =0; a < cur_v; a++) 
      {if (tabrn[a].tab) 
        {  
          free(tabrn[a].tab); break;
        }
      }
      continue;
    }
  uint64_t ** root_bns = (uint64_t**)malloc(ngroot*sizeof(uint64_t*));
  e = 0;
  for (a =0; a < cur_v; a++)
    {if (tabrn[a].tab) 
      {root_bns [e] = tabrn[a].tab;
       root_ritem[e++] = a;        
      }
    }  
  widthwalk_TC (root_bns,root_ritem,nbc, ngroot, nb_freq_); 
  free(root_bns);         

  for ( a = 0; a < cur_v; a++)
   {if (tabrn[a].tab) free(tabrn[a].tab);}
  }  
}

void Bitsufrec_TC (uint64_t** Bitdata,int * sum_1freq, info_prog * ip)
{
 minsup_TC = ip->minsup;
 int nb_thread = ip->nb_thread;
 for (int i = 0; i < 64; i++) {PT_Ultab_TC[i] = (1UL <<i);}
  int nz = 0; 
  std::cout << "The Support value is " << minsup_TC << ". Checking order parameters and initialize data ... \n ";
  init_sufixtree(Bitdata, sum_1freq, ip->varnames, minsup_TC,ip->nvar,ip->maxul,ip->ordre);
  int max_thread = std::thread::hardware_concurrency();
  if (nb_thread > max_thread) {nb_thread = max_thread-2;}
  if (nb_thread < 1) {nb_thread = 1;}
  std::cout << "Start Sufreq with  " << ip->nvar << " frequents variables using   " << nb_thread << " threads \n";
  auto start = std::chrono::system_clock::now();
  
  int shared_int = 0;
  int th_nb_freq  [nb_thread]{0};
  std::thread vec_thread [nb_thread];
  for (nz = 0; nz < nb_thread; nz++)
    {
      vec_thread[nz] = std::thread (thread_witdh,Bitdata, ip->maxul, &th_nb_freq[nz],&shared_int,ip->nvar);
    }

  int tot_freq = 0;
  for (nz = 0; nz < nb_thread; nz++)
    {
      vec_thread[nz].join();
      tot_freq += th_nb_freq [nz];
    }



  std::chrono::duration<double> dur= std::chrono::system_clock::now() - start;
  std::cout << "End Sufreq, there is  " << tot_freq <<" sets, extracted int  " <<dur.count() << " secs"<< std::endl;
  ip->time = dur.count();
  ip->nb_freq = tot_freq;
    
}


