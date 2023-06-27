#include "Sufrec_functions_TO.h"
#pragma optimization_level 3
#pragma GCC optimize("Ofast,no-stack-protector,unroll-loops,fast-math,O3")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx")
#pragma GCC optimize("Ofast")
#pragma GCC target("avx,avx2,fma")
#pragma GCC optimization ("unroll-loops")

 int minsup_TO;
 static uint64_t Ultab_TO [65];
 uint64_t * PT_Ultab_TO = &Ultab_TO[0];




void creabitfield_TO ( int *ind,int* mod,uint64_t * ptul, uint64_t * nlgtabl,  int nbc)
{   auto * yb = &Ultab_TO[0]; uint64_t u; int t = 0; auto * ybi = &Ultab_TO[0];
    auto * yd = &Ultab_TO[64];
    for (auto a =0; a < nbc ; a++)
    {  
      for (yb = ybi; yb != yd; yb++)
      { 
        u ^= (-((ptul[ind[t]]>>mod[t++]) & 1UL) ^ u) & *yb;
      }
      nlgtabl[a] = u;
    }
}

void creabitfieldr_TO ( int *ind,int * mod,uint64_t * ptul, uint64_t * nlgtabl,  int nbdases,  int reste)
{   
    auto * yb = &Ultab_TO[0]; uint64_t u; int t = 0;
    auto * yd = &Ultab_TO[64]; auto * ybi = &Ultab_TO[0];
    for (auto a =0; a < nbdases ; a++)
    {  
      for (yb = ybi; yb != yd; yb++)
      { 
        u ^= (-((ptul[ind[t]]>>mod[t++]) & 1UL) ^ u) & *yb;
      }
      nlgtabl[a] = u;
    }
    u = 0;
    for ( yb = ybi; yb != &Ultab_TO[reste]; yb++) 
    {
     u ^= (-((ptul[ind[t]]>>mod[t++]) & 1UL) ^ u) & *yb;}
    nlgtabl[nbdases] = u;

}



void init_Rnodes_TO (rnodes* outrnodes, int & ngroot  ,int* vecsom,  int* tindices ,int * mod,uint64_t** Bitdata, int nbleft, int nbcases,  int reste, int * nb_freq_)
{ 
  
  int sum = 0; 
  if (reste ==0)
    {for (auto i=0; i < nbleft;i++)
        { sum = vecsom[i];
          if (sum<minsup_TO) {outrnodes[i].tab=NULL; continue;}
          ngroot++;
          uint64_t * ptu = (uint64_t*) malloc(nbcases*SUL);
          creabitfield_TO(tindices,mod,Bitdata[i],ptu,nbcases);
          outrnodes[i].sup = sum;
          outrnodes[i].tab=ptu;
        }
      (*nb_freq_) += ngroot; 

    }

  else {
        for (auto i=0; i < nbleft;i++)
          { sum = vecsom[i];
            if (sum<minsup_TO) {outrnodes[i].tab=NULL; continue;}
            ngroot++;
            uint64_t * ptu = (uint64_t*) malloc(nbcases*SUL);
            creabitfieldr_TO(tindices,mod,Bitdata[i],ptu,nbcases-1,reste);
            outrnodes[i].sup = sum;
            outrnodes[i].tab=ptu;
          }
        (*nb_freq_) += ngroot;
      }
}



void depthwalk_TO (pnodes * curfather, uint64_t ** Bitv, int * tab_ritem_ ,int nbc, int nb_desc, int * nb_freq_)
{ 
  pnodes * trpn = curfather->son->brother;
  
  pnodes ** move_trpn = NULL;
   int supcand = 0; int ritem_ ;
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
            if (supcand < minsup_TO) continue; 
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
        
      (*nb_freq_) += tr_desc;    
      if (tr_desc > 1 ) {depthwalk_TO (trpn, trt_bns, tr_ritem_ ,nbc, tr_desc,nb_freq_);}
      for (j = 0; j < tr_desc; j++ ) { free(trt_bns[j]);} free(trt_bns); free(tr_ritem_);
        
      trpn = trpn->brother;  
      
    }
}


void thread_witdh2 ( uint64_t ** Bitdata, int * th_ind, int nvar_th, int maxul, int * nb_freq_, pnodes ** tab_subtree_)
{  
  int t, sum, e , ngroot;
  uint64_t * pul = NULL;
  uint64_t * pulb = NULL;
  auto a =0;
  uint64_t u;
  int n = 0;
  int cur_v;
  t = 1;
  for (t = 0; t < nvar_th; t++)
  { 
    cur_v = th_ind[t];
   
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

    pnodes * nouvnod = (pnodes*)malloc(SPN);
    nouvnod->brother=NULL;
    nouvnod->son=NULL;
    nouvnod->ritem = cur_v;
    nouvnod->sup = sum;
    tab_subtree_[cur_v] = nouvnod ;

    ngroot= 0;
 
 
 init_Rnodes_TO(tabrn,ngroot,rootsum,indices,mod,Bitdata,cur_v,nbc,reste,nb_freq_);
 free(indices);
 free(mod);
  pnodes ** trpnr = &nouvnod->son;
  if (ngroot <1) { continue;}
  if (ngroot ==1) 
    {for (a =0; a < cur_v; a++) 
      {if (tabrn[a].tab) 
        {  *trpnr = (pnodes*)malloc(SPN);
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
  for (a =0; a < cur_v; a++)
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
  depthwalk_TO (nouvnod,root_bns,root_ritem,nbc, ngroot, nb_freq_); 
  free(root_bns);         

  for ( a = 0; a < cur_v; a++)
   {if (tabrn[a].tab) free(tabrn[a].tab);}
  }  
}


void set_ind_th (int * ind, int num, int nvar_thread, int nb_thread ,int a )
{

  int h = 0;
    if (a == 0)
      { 
        for (int j = 0; j <= (nvar_thread-2)/2; j++)
          {
            ind[h] = 2*(nb_thread)*j+num-1;
            h = h+2;
          }
        h = 1;
        for (int j = 1; j <= (nvar_thread)/2; j++)
          {
            ind[h] = 2*(nb_thread)*j-(num-1)-1;
            h = h+2;
          }       
      }
    if (a != 0)
      { 
        for (int j = 0; j <= (nvar_thread-1)/2; j++)
          {
            ind[h] = 2*(nb_thread)*j+num-1;
            h = h+2;
          }
        h = 1;
        for (int j = 1; j <= (nvar_thread-1)/2; j++)
          {
            ind[h] = 2*(nb_thread)*j-(num-1)-1;
            h = h+2;
          }       
      }


}

void Bitsufrec_TO (uint64_t** Bitdata,int * sum_1freq, info_prog * ip, pnodes *** tab_subtree_)
{
 minsup_TO = ip->minsup;
 int maxul = ip->maxul;
 int nvar = ip->nvar;
 int nb_thread = ip->nb_thread;
 for (int i = 0; i < 64; i++) {PT_Ultab_TO[i] = (1UL <<i);}
  int nz = 0; 
  std::cout << "The Support value is " << minsup_TO << ". Checking order parameters and initialize data ... \n ";
  init_sufixtree(Bitdata, sum_1freq, ip->varnames, minsup_TO,nvar,maxul,ip->ordre);
  int max_thread = std::thread::hardware_concurrency();
  if (nb_thread > max_thread) {nb_thread = max_thread-2;}
  if (nb_thread < 1) {nb_thread = 1;}

  std::cout << "Start Sufreq with  " << nvar << " frequents variables using   " << nb_thread << " threads \n";
  auto start = std::chrono::system_clock::now();
  (*tab_subtree_) = (pnodes**)malloc(nvar*sizeof(pnodes*));
  
  int * nvar_thread = (int*)calloc(nb_thread,sizeof(int));
  int ** th_ind = (int**)malloc(nb_thread*sizeof(int*));
  
  int th_nb_freq  [nb_thread]{0};
  std::thread vec_thread [nb_thread];
  int nvar_th = nvar/nb_thread;
  if (nvar % nb_thread != 0) {nvar_th = nvar_th +1;}
  for (int i = 0; i < nb_thread; i++) {th_ind[i] = (int*)malloc((nvar_th+2)*sizeof(int)); nvar_thread[i]=nvar_th;}  
  int a = nvar_th % 2;

  for (int i = 0; i < nb_thread; i++)
    {
      set_ind_th (th_ind[i],i+1,nvar_th,nb_thread ,a );
      if (th_ind[i][nvar_th-1] >= nvar) {nvar_thread[i] = nvar_th-1;}
    
    }
  
  for (nz = 0; nz < nb_thread; nz++)
    {
      vec_thread[nz] = std::thread (thread_witdh2,Bitdata, th_ind[nz], nvar_thread[nz], maxul, &th_nb_freq[nz], *tab_subtree_);
    }

  int tot_freq = 0;
  for (nz = 0; nz < nb_thread; nz++)
    {
      vec_thread[nz].join();
      tot_freq += th_nb_freq [nz];
    }



  std::chrono::duration<double> dur= std::chrono::system_clock::now() - start;
  std::cout << "End Sufreq, there is  " << tot_freq <<" sets, extracted int  " <<dur.count() << " secs"<< std::endl;
  ip->time = (double)dur.count();
  ip->nb_freq = tot_freq;
  for (int i = 0; i < nb_thread; i++) {free(th_ind[i]);} free(th_ind); free(nvar_thread);
  std::cout << "ok \n";
  
}


