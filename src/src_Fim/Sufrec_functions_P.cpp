#include "Sufrec_functions_P.h"


#pragma optimization_level 3
#pragma GCC optimize("Ofast,no-stack-protector,unroll-loops,fast-math,O3")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx")
#pragma GCC optimize("Ofast")
#pragma GCC target("avx,avx2,fma")
#pragma GCC optimization ("unroll-loops")

int nb_freq_P ; // global int to count number of frequent itemstes
int minsup_P; // global int support value
static uint64_t Ultab_P [65];
uint64_t * PT_Ultab_P = &Ultab_P[0];




void creabitfield_P ( int *ind,int* mod,uint64_t * ptul, uint64_t * nlgtabl,  int nbc) // function to reduce array : only index of the new item
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

void creabitfieldr_P ( int *ind,int * mod,uint64_t * ptul, uint64_t * nlgtabl,  int nbdases,  int reste) // function to reduce array : only index of the new item
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


// stock array of frequent 2-itemsets
// it will be use with labelmin info , to create new candidates
void init_Rnodes_P (rnodes* outrnodes, int * ngroot ,int* vecsom,  int* tindices ,int * mod,uint64_t** Bitdata, int nbleft, int nbcases,  int reste)
{ 
  
  int sum = 0; 
  if (reste ==0)
    {for (auto i=0; i < nbleft;i++)
        { sum = vecsom[i];
          if (sum<minsup_P) {outrnodes[i].tab=NULL; outrnodes[i].sup = 0; continue;}
          (*ngroot)++;
          
          uint64_t * ptu = (uint64_t*) malloc(nbcases*SUL);
          creabitfield_P(tindices,mod,Bitdata[i],ptu,nbcases);
          outrnodes[i].sup = sum;
          outrnodes[i].tab=ptu;
          ;
        }
       nb_freq_P += (*ngroot);   
    }

  else {
        for (auto i=0; i < nbleft;i++)
          { sum = vecsom[i];
            if (sum<minsup_P) {outrnodes[i].tab=NULL;  outrnodes[i].sup = 0; continue;}
            (*ngroot) ++;           
            uint64_t * ptu = (uint64_t*) malloc(nbcases*SUL);
            creabitfieldr_P(tindices,mod,Bitdata[i],ptu,nbcases-1,reste);
            outrnodes[i].sup = sum;
            outrnodes[i].tab=ptu;
            
          }
        nb_freq_P += (*ngroot); 
        
      }
}


// Principal function, check if itemsets is "to test" thanks to the pruning of nephew
// Check sup of the new candidate 
// if > minusp, add the new node in tree
// Continue visit of old tree and construction of the new tree
void withwalk_P (pnodes * Tree,pnodes ** curtree, uint64_t * curb, int cur_r ,rnodes * tab_r ,int nbc)
{
  
  if (!tab_r[cur_r].sup)
    {  
       if (Tree->brother) {withwalk_P(Tree->brother,curtree,curb,Tree->brother->ritem,tab_r,nbc);}
       return;
    }
  
  auto * tib = &curb[0];
  auto * ed = &curb[nbc];
  auto * cb = &tab_r[cur_r].tab[0];
  int sum = 0;


  for (tib; tib != ed; tib++, cb++ )
    {sum+= __builtin_popcountl(*tib&*cb);}
  if (sum < minsup_P) 
    {  tab_r[cur_r].sup = 0;
       if (Tree->brother) {withwalk_P(Tree->brother,curtree,curb,Tree->brother->ritem,tab_r,nbc);}
       tab_r[cur_r].sup = 1;
       return;        
    }

    nb_freq_P++;
    
    pnodes * ntree = (pnodes*)malloc(SPN);
    ntree->sup = sum;
    ntree->ritem = Tree->ritem;
    ntree->son = NULL;
    ntree->brother = NULL;
    (*curtree) = ntree;

    if (!Tree->son)
     {   
         if (Tree->brother) {withwalk_P(Tree->brother,&(ntree->brother),curb,Tree->brother->ritem,tab_r,nbc);}
         return;              
     }
     
    uint64_t * ncurb = (uint64_t*)malloc(nbc*SUL);
    auto * ad = &ncurb[0];
    cb = &tab_r[cur_r].tab[0];
    for (tib = &curb[0]; tib != ed; tib++, cb++,ad++) {*ad = *tib&*cb;}

    withwalk_P(Tree->son,&(ntree->son),ncurb,Tree->son->ritem,tab_r,nbc);
    free(ncurb);
    if (Tree->brother) {withwalk_P(Tree->brother,&(ntree->brother),curb,Tree->brother->ritem,tab_r,nbc);}

}



void Bitsufrec_P (uint64_t** Bitdata, int * sum_1freq, info_prog * ip, pnodes *& root_)
{
 nb_freq_P = 0;
 minsup_P = ip->minsup;
 int maxul = ip->maxul;
 int nvar = ip->nvar;
 char ordre = ip->ordre;
 for (int i = 0; i < 64; i++) {PT_Ultab_P[i] = (1UL <<i);}
 std::cout << "The Support value is " << minsup_P << ". Checking order parameters and initialize data ... \n ";
 int sum = 0;
 int nbc = 0;
 auto a = 0; int b= 0;
 auto e = 0;
 auto n = 0;
 uint64_t ** Bdata = Bitdata ;
 auto * pul = &Bitdata[0][0];
 auto * pulb = &Bitdata[0][0];
 auto * pulc = &Bitdata[0][0];
 int keep = nvar;
 init_sufixtree(Bitdata, sum_1freq, ip->varnames, minsup_P,nvar,maxul,ordre);
 pnodes * last_tree = NULL;

 std::cout <<" Start PrefRec with  " << keep << " frequents variables ..." <<  std::endl;
 if(keep ==0)
  {
    std::cout << "Not a single frequent itemSet, please try with a lower relative minsup_P value. " << std::endl;
    return;
  }
 root_ = (pnodes*)malloc(SPN);
 root_->brother = NULL;
 root_->son = NULL;
 auto start = std::chrono::system_clock::now();
 pnodes rootalpha ;
 pnodes * lastree = NULL;
 lastree = root_; 
 int j, ngroot;
 uint64_t u;
 for (int j=0; j<keep;j++)
  { 

    nb_freq_P++;   
    j = j;
    pul = Bdata[j];
    sum = 0;
    e =0;
    for (a = 0; a < maxul;a++)
      {sum+= __builtin_popcountl(pul[a]);}
    
    rnodes tabrn [j+1]; 
    int * indices = (int*)malloc(sum*sizeof(int));
    int * mod = (int*)malloc(sum*sizeof(int));
    for (a=0; a < maxul; a++)
      { u = pul[a];
        for (n = 0; n < 64; n++)
          { if ((u >> n) & 1UL) {indices[e]= a; mod[e++]=n;} }
      }
 
    
    int rootsum [j+1];
    int root_ritem [j+1];    
    int nbc = sum/64;
    int reste = sum-64*nbc;
    if (reste) nbc ++;
   
    int s = 0;
    for (a = 0; a < j; a++)
        {
          s = 0; pulb = Bdata[a];
          for (n = 0; n < maxul; n++)
            { s += __builtin_popcountl(pul[n]&pulb[n]);}
          rootsum[a] = s;
        }
    pnodes * nouvnod = (pnodes*) malloc (SPN);
    nouvnod->sup = sum;
    nouvnod->ritem = j;
    nouvnod->son = NULL;
    nouvnod->brother = NULL;
    lastree->brother = nouvnod;
    lastree  = nouvnod;  
    ngroot= 0;

    init_Rnodes_P(tabrn,&ngroot,rootsum,indices,mod,Bdata,j,nbc,reste);
    pnodes ** trpnr = &nouvnod->son;
    pnodes  * temp = root_->brother;
    free(indices);
    free(mod);
  
    
    if (ngroot <1) {continue;}
    if (ngroot ==1) 
      {for (a =0; a < j; a++)
        {if (tabrn[a].tab) 
            {*trpnr = (pnodes*)malloc(SPN); (*trpnr)->sup = tabrn[a].sup; (*trpnr)->ritem = a;(*trpnr)->son = NULL;(*trpnr)->brother=NULL;
              free(tabrn[a].tab); break;
            }

        }
        continue;
      } 

    e = 0;
    for (a =0; a < j; a++)
      { if (!tabrn[a].tab) 
          {   
            temp = temp->brother;
            continue;
          }
      
        *trpnr = (pnodes*)malloc(SPN); (*trpnr)->sup = tabrn[a].sup; (*trpnr)->ritem = a;(*trpnr)->son = NULL;(*trpnr)->brother=NULL;   
        if (temp->son) 
          { 
            
            withwalk_P(temp->son,&(*trpnr)->son,tabrn[a].tab,temp->son->ritem,tabrn,nbc);
          }  
        temp = temp->brother;
        trpnr = &(*trpnr)->brother; 
      }        
    for ( a = 0; a < j; a++)
      {if (tabrn[a].tab) free(tabrn[a].tab);}
  }  
  

  std::chrono::duration<double> diffe= std::chrono::system_clock::now() - start;
  std::cout << "End PrefRec, there is  " << nb_freq_P <<" sets, extracted int  " <<diffe.count() << " secs"<< std::endl;

return; 
    
}