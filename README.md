# SufRec

## Description 

This page propose different versions of the **SufRec** C++ programs which can be used bu using GCC compiler and the command make.

These algorithms were built around the notion of *Mining Frequent ItemSets and associations rules*. This data analysis method was first introduced by *Agrawal et al. 1993* for mining transaction databases.


Its use and practices are close to those of PrefRec.
With SufRec, we propose mainly 2 uses, one of them allows an independent construction of trees and an efficient use of multi threads. 

Moreover, we give possibilities of recursive mooving applications, their details are given in the corresponding article.


**SufRec:**
* (recursiv) Mining frequent itemSets with a relative Support
* Get an output Tsv file for information about frequent itemSet
* Get 2 files necessary to launch Prefrules 


## Creating binaries and getting started
```
cd src/src_SufRecP && make
cd src/src_SufRecI && make
cd src/src_SufRecTh && make
cd src/src_Mooving_SufRecP && make
cd src/src_Mooving_SufRecI && make

SufRecP and SufRecI : ./SufRec <transaction dataset> < d=item_delimitator> <s=minimal_relative_support "s"> <output_file_set_infos> <output file for coefficient>
SufRecTH : ./SufRec <transaction dataset> < d=item_delimitator> <s=minimal_relative_support "s"> < n = number_of_threads> <output_file_set_infos> <output file for coefficient> 

Mooving_SufRecP and Mooving_SufRec_I :

```

## Dataset format 

For this version, the Dataset format accepted is a file of type transaction, such that the items are separated by a choosen separator
and the transaction by a new line.Please, now that the sep must be one of the ASCII Chart. 
The following tab is an example with sep=, item a b c and d, and five transactions:



|transactions|
|------------|
|a,b|
|c,d|
|a,c|
|a,b,d|
|c|


The folder sample contain a simple small dataset test. You can use it for the example below to see how the software works.


### Example
```
SufRecP et SufRecI

./SufRec ../sample/input/Fruits.txt d=, s=0.20  \\ Do the extraction on the bases "fruits" with a relativ minSup of 0.20, item are delimited by a ",". 
./SufRec ../sample/input/Fruits.txt d=, s=0.20 ../sample/output_Prefrec/infoset \\ Do the extraction and write frequent set informations in infoset   
./SufRec ../sample/input/Fruits.txt d=, s=0.20 ../sample/output_Prefrec/infoset ../sample/output_Prefrec/genrules \\ Do the extraction, write frequent set informations in info set,  and create 2 files , genrules.txt and genrules_item.txt, necessary to use Prefrules

SufRecTH

./SufRec ../sample/input/Fruits.txt d=, s=0.20 n=4  \\ Do the extraction on the bases "fruits" with a relativ minSup of 0.20, item are delimited by a ",", and the number of threads if egal to 4. 
./SufRec ../sample/input/Fruits.txt d=, s=0.20 n=4 ../sample/output_Prefrec/infoset \\ Do the extraction and write frequent set informations in infoset    
./SufRec ../sample/input/Fruits.txt d=, s=0.20 n=4 ../sample/output_Prefrec/infoset ../sample/output_Prefrec/genrules \\ Do the extraction, write frequent set informations in info set,  and create 2 files , genrules.txt and genrules_item.txt, necessary to use Prefrules
```

## Parameters for SufRec :
|param|required|note|
|--------------------|--------|--------|
|    transaction dataset    |    yes    | The Dataset transaction  |  
|    item delimitator "d="   |    yes    | the item separator you use with your dataset transaction | 
|    minimal relative support "s="   |    yes    | the minimal relativ support you wish     | 
|    ordering frequent 1-itemset "o="   |    no    | n for unordered, i for ascendant order, d for decreasing order      | 
|    output file set infos    |    no    |  For file frequent set informations, put the directory file (with no extention type )    | 
|    output file for coefficient    |    no    |  If you want to use Prefrules, put the directory file (with no exention type)| 



## Parameters for SufRecTH :
|param|required|note|
|--------------------|--------|--------|
|    transaction dataset    |    yes    | The Dataset transaction  |  
|    item delimitator "d="   |    yes    | the item separator you use with your dataset transaction | 
|    minimal relative support "s="   |    yes    | the minimal relativ support you wish     | 
|	 number of thread "n="	| yes | the number of thread you wan't to use |
|    ordering frequent 1-itemset "o="   |    no    | n for unordered, i for ascendant order, d for decreasing order      | 
|    output file set infos    |    no    |  For file frequent set informations, put the directory file (with no extention type )    | 
|    output file for coefficient    |    no    |  If you want to use Prefrules, put the directory file (with no exention type)| 





## Definitions of frequent itemSets :

Let us remind you the 2 mains definitions of this data analys method

Let I = {a1, . . . , an} be a finite set of items. A transaction database is a set of transactions T =
{t1, . . . , tN } where each transaction ti ⊂ I, 1 ≤ i ≤ N, represents a nonempty
subset of items. An itemset A is a subset of I; A is a k-itemset if it contains
k items. The support of an itemset A is denoted as supp(A) and is defined
as the number of transactions which contain A. The relative support of A is
freq(A) = supp(A)/N. A is frequent if freq(A) ≥ σ where σ is a user-specified minimum relative support threshold, called minSup.


**Definitions of confident associations rules  :**
An association rule is an implication A ⇒ B where A and B are two itemsets. The support of a rule A ⇒ B is defined as sup(A ⇒ B) = sup(A∪B).
The confidence of a rule A ⇒ B is defined as conf(A ⇒ B) = supp(A ⇒B)/supp(A) = freq(A∪B)/freq(A).
Considering a treshold minconf, a rule such that A ⇒ B is confident if conf (A ⇒ B) > minconf.






