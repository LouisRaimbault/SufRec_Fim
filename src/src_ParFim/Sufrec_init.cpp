#include "Sufrec_init.h"


const size_t SUL = sizeof(uint64_t);
const size_t SI = sizeof(int);
const size_t SPN = sizeof(pnodes);
const size_t SRN = sizeof(rnodes);



void Init_data (char * pathfile, uint64_t *** Bitdata, int ** sum_1freq, info_prog * ip ) 
{ 
  // init version , slower but safe for larger data
  auto start = std::chrono::system_clock::now(); // Start time init
  std::chrono::duration<double> diffe;
  std::cout << "import data from file as char array and count number of line ... \n";
  char * tabchar;
  long longueur = 0;long i = 0;
  char c;
  int w =0; int t = 0; int z = 0;
  int nvar = 0; int nrows = 0; int maxul = 0;
  char deli = ip->delim;
  std::string st ="";
  FILE *Fichier = fopen(pathfile,"r");  // open path to data
 
  if (Fichier != NULL)
    {                         // count number of lines and char
      while (!feof(Fichier))  
       {c = getc(Fichier);
        longueur ++;
        if (c == '\n') nrows++; }                      

    }
  fclose(Fichier); //close file
  tabchar = (char*)malloc(longueur * sizeof(char)); // create array of char
  Fichier = fopen(pathfile,"r"); // open again file
  if ( Fichier != NULL && tabchar  != NULL)
    {                           // Fill array
     while (i < longueur)  
      { c = getc(Fichier);
        tabchar[i++]=c;}
      }
  else perror ("\n\n problem in file ");
  fclose(Fichier); // close file
  
  diffe= std::chrono::system_clock::now() - start; // stop timer
  std::cout << "done in " << diffe.count() << " secondes \t";

  std::cout << "Assignment of a number to the different variables and create transation as integer array .. \n";
  start = std::chrono::system_clock::now();
  int ** TRS = (int**)malloc(nrows *sizeof(int*)); // transac items 
  int * n_TRS = (int*)malloc(nrows *sizeof(int)); // nb item per transac
  int sup_init = nrows * ip->rel_minsup;

  
  std::string a_st ="";
  int c_i = 1;
  int z_i;
  int c_nr = 0;
  int * P_TRS = NULL;
  i = 0;
  std::unordered_map<std::string,int> mappy; // give numer to id items

  // count sup for each item
  for (t = 0; t < nrows;t++)
    { 
      c_nr = 0;
      z = i;
      while (tabchar[i] != '\n')
        { 
          
          if (tabchar[i]==deli) 
            {
              c_nr++;
            }
         i++;
          
        }
      c_nr++;  
      TRS[t] = (int*) malloc(c_nr*sizeof(int));
      n_TRS[t] = c_nr;
      c_nr = 0;
      P_TRS = TRS[t];
      for (z; z < i; z++)
        {
          if (tabchar[z] == deli)
            {
              if (!mappy[a_st])
                {
                  mappy[a_st] = c_i;
                  P_TRS[c_nr] =c_i;
                  c_i++;
                  a_st = "";
                  c_nr++;          
                  continue;
                }

              P_TRS[c_nr] = mappy[a_st];
              a_st = "";
              c_nr++;
              continue;
            }
          a_st += tabchar[z];
        }
      
      if (!mappy[a_st])
        { 
          
          mappy[a_st] = c_i;
          P_TRS[c_nr] =c_i; 
          c_i++;
          a_st = "";
          i++;          
          continue;
        }

      P_TRS[c_nr] = mappy[a_st];
      a_st = "";
      i++;
    }

  free (tabchar);
  diffe= std::chrono::system_clock::now() - start;
  std::cout << "done in " << diffe.count() << " secondes \t";
  
  maxul = nrows/64;
  int rst = nrows-64*maxul;
  if (rst) maxul++;
  std::unordered_map<std::string,int>::iterator itmap;
  std::unordered_map<std::string,int>::iterator itendmap = mappy.end(); 
  int nvar_tot = mappy.size()+1;
  std::cout << "number of raws : " << nrows << " and total number of variables : " << nvar_tot-1 << "\n";

  
  std::cout << "Checking for frequent 1-item ... \n";
  start = std::chrono::system_clock::now();
  int * item_sup = (int*)calloc(nvar_tot,SI); // item support
  int * is_good_item = (int*) malloc((nvar_tot)*sizeof(int)); // if item sup > minsup
  uint64_t ** Bdata_tot = (uint64_t**)malloc ((nvar_tot+10) * sizeof(uint64_t*));


  for (t = 0; t < nrows; t++)
    {
      w = n_TRS[t];
      P_TRS = TRS[t];
      for (z = 0; z < w; z++)
        {
          item_sup[P_TRS[z]]++;
        }
    }

  diffe= std::chrono::system_clock::now() - start;
  std::cout << "Done in " << diffe.count() << " secondes, check for corresponding variables ...\n";
  
  for (t = 0; t < nvar_tot;t++)  {is_good_item[t] = 0;}
  t = 0;
  // Create first Bit data
  for (itmap = mappy.begin(); itmap != itendmap;itmap++) 
    { 

      if (item_sup[itmap->second] >= sup_init)
        {  
          Bdata_tot[itmap->second] = (uint64_t*)calloc(maxul,SUL);
          is_good_item[itmap->second] = 1;
          nvar++;
          continue;
        }    
    } 


  if (nvar == 0)
    {
      free (Bdata_tot);
      for (t = 0; t < nrows;t++) {free(TRS[t]);}
      free(TRS);
      free(n_TRS);
      free(is_good_item);
      free(item_sup);
      return;
    }
  t = 0;  
  std::unordered_map<int,std::string> map_names;

  for (itmap = mappy.begin(); itmap != itendmap;itmap++) 
    { 

      if (is_good_item[itmap->second])
        { 
          map_names[itmap->second] = itmap->first; 
        }    
    } 

  // Fill Bit data for frequent  item
  std::cout << "Create Bit Data store in uint64t aray with " << nvar << " frequent 1 -item ... \n";
  start = std::chrono::system_clock::now();
  z_i = 0; t = 0;
  int k = 0;
  uint64_t u =0;
  while(z_i < maxul-1)
  { 
    for (t =0; t < 64; t++,k++)
      { P_TRS = TRS[k];
        u = 1UL << t;
        for (z = 0; z < n_TRS[k]; z++)
          {
            if (is_good_item[P_TRS[z]])
              {
                Bdata_tot[P_TRS[z]][z_i]+= u ;
              }
          }
      }
    z_i++;
  }

  for (t = 0; t < rst; t++,k++)
      { 
        P_TRS = TRS[k];
        u = 1UL << t;
        for (z  = 0; z < n_TRS[k];z++)
          {
            if (is_good_item[P_TRS[z]])
              {
                Bdata_tot[P_TRS[z]][z_i] +=u;
              }
          }
      }     

  z= 0;    
  (*Bitdata) = (uint64_t**)malloc(nvar*sizeof(uint64_t*));
  std::string * varnames = new std::string [nvar];
  (*sum_1freq) = (int*)malloc(nvar*sizeof(int));
  for (t = 0; t < nvar_tot;t++)
    {
      if (is_good_item[t]) {(*sum_1freq)[z] = item_sup[t]; varnames[z] = map_names[t]; (*Bitdata)[z++] = Bdata_tot[t];}
    }

  // Reduce first Bitdata to only pointers of frequent 1 item
  free (Bdata_tot);
  for (t = 0; t < nrows;t++) {free(TRS[t]);}
  free(TRS);
  free(n_TRS);
  free(is_good_item);
  free(item_sup);
  diffe= std::chrono::system_clock::now() - start;
  std::cout << "Done in " << diffe.count() << " secondes\n";
  ip->minsup = sup_init;
  ip->nrows = nrows;
  ip->nvar = nvar;
  ip->maxul = maxul;
  ip->varnames = varnames;

 return ;
}

void set_ind_machine (info_prog * ip)
{
  ip->nvar_mach = ip->nvar/ip->nb_mach;
  if (ip->nvar%ip->nb_mach != 0) {ip->nvar_mach = ip->nvar_mach+1;}
  ip->ind_mach = (int*)malloc(ip->nvar_mach*sizeof(int));
  int h = 0; int num_m = ip->num_mach;

    if (ip->nvar_mach % 2 == 0)
    {
      for (int j = 0; j <= (ip->nvar_mach-2)/2; j++)
        {
          ip->ind_mach[h] = 2*(ip->nb_mach)*j+num_m-1;
          h = h+2;
        }
      h = 1;
      for (int j = 1; j <= (ip->nvar_mach)/2; j++)
        {
          ip->ind_mach[h] = 2*(ip->nb_mach)*j-(num_m-1)-1;
          h = h+2;
        }       
    }
  if (ip->nvar_mach % 2 != 0)
    {
      for (int j = 0; j <= (ip->nvar_mach-1)/2; j++)
        {
          ip->ind_mach[h] = 2*(ip->nb_mach)*j+num_m-1;
          h = h+2;
        }
      h = 1;
      for (int j = 1; j <= (ip->nvar_mach-1)/2; j++)
        {
          ip->ind_mach[h] = 2*(ip->nb_mach)*j-(num_m-1)-1;
          h = h+2;
        }       
    }  

  if (ip->ind_mach[ip->nvar_mach-1] >= ip->nvar ) {ip->nvar_mach = ip->nvar_mach-1;} 

}





void sortoutdata (uint64_t** Bdata, int * sumvec, int nb_elem, std::string * varnames, char ordre)
{
  // order either deacreasing or increasing sup 1-item.
  // Naive function, can be upgraded
  int i,j,treme,itreme,temp;
  uint64_t * tempt;
  std::string tempstr;
  if (ordre=='d')
  {  
    for (i = 0; i<nb_elem-1;i++)
    {itreme = i; treme=sumvec[i];
    for (j=i+1;j<nb_elem;j++) {if (sumvec[j]>treme) {treme = sumvec[j];itreme=j;}}
    temp=sumvec[itreme];sumvec[itreme]=sumvec[i];sumvec[i]=temp;
    tempt = Bdata[itreme]; Bdata[itreme] = Bdata[i];Bdata[i]=tempt;
    tempstr = varnames[itreme]; varnames[itreme] = varnames[i];varnames[i]=tempstr;
    }
  }
  if (ordre == 'i')
  {  
    for (i = 0; i<nb_elem-1;i++)
    {itreme = i; treme=sumvec[i];
    for (j=i+1;j<nb_elem;j++) {if (sumvec[j]<treme) {treme = sumvec[j];itreme=j;}}
    temp=sumvec[itreme];sumvec[itreme]=sumvec[i];sumvec[i]=temp;
    tempt = Bdata[itreme]; Bdata[itreme] = Bdata[i];Bdata[i]=tempt;
    tempstr = varnames[itreme]; varnames[itreme] = varnames[i];varnames[i]=tempstr;
    }
  }
}

void init_sufixtree (uint64_t ** Bitdata,int * sum_1freq, std::string * varnames, int support, int nvar, int maxul, char ordre ) // order if i or d
{ 
    if (ordre == 'i' || ordre == 'd') 
    {  
      sortoutdata(Bitdata, sum_1freq,nvar,varnames,ordre);
    }
    
}

void erase_freq_spec (pnodes* Tree, int * nb_erase, int size) // only erase all frequent itemsets
{ 


 if (Tree->son) {erase_freq_spec(Tree->son,nb_erase, size+1);}
 if (Tree->brother) {erase_freq_spec(Tree->brother,nb_erase,size);}
 (*nb_erase)++;
 free(Tree);

}


void erase_freq (pnodes* Tree) // only erase all frequent itemsets
{ 
 if (Tree->son) {erase_freq(Tree->son);}
 if (Tree->brother) {erase_freq(Tree->brother);}
 free(Tree);

}


void extract_and_erase_freq (pnodes* Tree, std::string str, std::string *  listenom, double nrs , int size ,std::ofstream  & flux_set ) // erase and extract fréquent itemsets
{ 
 
 std::string strtemp = str + listenom[Tree->ritem];
 flux_set << strtemp << "\t" << Tree->sup << "\t" <<(double)Tree->sup/nrs << "\t" << size <<"\n"; 

 if (Tree->son) {extract_and_erase_freq(Tree->son,strtemp,listenom,nrs,size+1,flux_set);}
 if (Tree->brother) {extract_and_erase_freq(Tree->brother,str,listenom,nrs,size,flux_set);}
 free(Tree);
}


void extract_and_erase_freq (pnodes* Tree, std::string str, std::string *  listenom, double nrs, int size, std::ofstream  & flux_set, std::ofstream  & flux_rules ) // erase and extract frequents itemsets with another file to launch Rules function
{
 std::string strtemp = str + listenom[Tree->ritem];
 double rel_sup = (double)Tree->sup/nrs;
 flux_set << strtemp << "\t" << Tree->sup << "\t" <<rel_sup << "\t" << size <<"\n"; 
 flux_rules << rel_sup << "\t" << size << "\t" << Tree->ritem << "\n";  

 if (Tree->son) {extract_and_erase_freq(Tree->son,strtemp,listenom,nrs,size+1,flux_set,flux_rules);}
 if (Tree->brother) {extract_and_erase_freq(Tree->brother,str,listenom,nrs,size,flux_set,flux_rules);}
 free(Tree);
}
