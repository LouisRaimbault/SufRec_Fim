# SufRec

## Description 

**n is always the number of thread to use for multi thread programms**

**src_Fim:**
**Frequent itemsets Mining. 5 algorithms.** 
* p : Based on Suffix Tree dependant construction, for each new item, a new subtree is rooted by traversing the existing tree
* i : Based on Suffix Tree independant construction, for each new item, a new subtree is rooted using only the set of items 
* tfn : Based on Suffix Tree independant construction, each new subtree is built by the first available thread.
* ton : Based on Suffix Tree independant construction, each thread now the subtrees it has to build
* tcn : Based on Suffix Tree independant construction, more specially tf. Only count the number of frequent itemsets. No storage of itemsets in RAM.

**src_ParFim:**
**Frequent itemSets Mining for clusters. 3 algorithms. 
The subtrees to be built are defined in advance by the machine number and the total number of machines.**
* i : Based on Suffix Tree independant construction, for each new item, a new subtree is rooted using only the set of items 
* tfn : Based on Suffix Tree independant construction, each new subtree is built byt the first available thread.
* tcn : Based on Suffix Tree independant construction, more specially tf. Only count the number of frequent itemsets. No storage in RAM.


**src_Mooving:**
**The tree always contains the frequent itemsets of a given number of items. At each iteration, the item wich generates the smallest number of frequent itemset is deleted (with all its itemsets) and a new item (with all its frequent itemsets) is added.**
* p : Based on Suffix Tree dependant construction, for each new item, a new subtree is rooted by traversing the existing tree
* i : Based on Suffix Tree independant construction, for each new item, a new subtree is rooted using only the set of items 


## Parameters for SufRec in src_Fim :
|param|required|note|
|--------------------|--------|--------|
|    transaction dataset    |    yes    | The Dataset transaction  |  
|    item delimitator "d="   |    yes    | the item separator you use with your dataset transaction | 
|    minimal relative support "s="   |    yes    | the minimal relativ support you wish     | 
|    Type of algorithm "t="   |    yes    | i, p, tfn, ton, tcn   | 
|    ordering frequent 1-itemset "o="   |    no    | u for unordered, i for ascendant order, d for decreasing order      | 
|    output file set infos    |    no    |  For file frequent set informations, put the directory file (with no extention type )    | 
|    output file for coefficient    |    no    |  If you want to use Prefrules, put the directory file (with no exention type)| 



## Parameters for SufRec in src_ParFim :
|param|required|note|
|--------------------|--------|--------|
|    transaction dataset    |    yes    | The Dataset transaction  |  
|    item delimitator "d="   |    yes    | the item separator you use with your dataset transaction | 
|    minimal relative support "s="   |    yes    | the minimal relativ support you wish     | 
|    Type of algorithm "t="   |    yes    | i, tfn, tcn   | 
|    Total number of machines "m="   |    yes    | integer, the total number of machines  | 
|    Identifiant Number of machine "x="   |    yes    | integer, the number of machine  | 
|    ordering frequent 1-itemset "o="   |    no    | u for unordered, i for ascendant order, d for decreasing order      | 
|    output file set infos    |    no    |  For file frequent set informations, put the directory file (with no extention type )    | 
 


## Parameters for SufRec in src_Mooving :
|param|required|note|
|--------------------|--------|--------|
|    transaction dataset    |    yes    | The Dataset transaction  |  
|    item delimitator "d="   |    yes    | the item separator you use with your dataset transaction | 
|    minimal relative support "s="   |    yes    | the minimal relativ support you wish     | 
|    Type of algorithm "t="   |    yes    | i, p | 
|    Constant number of item   |    yes    | integer  | 
|    Number of iterations   |    yes    | integer | 
|    ordering frequent 1-itemset "o="   |    no    | u for unordered, i for ascendant order, d for decreasing order      | 
 





