# SufRec

## Description 

This page propose different versions of the **SufRec** C++ programs which can be used bu using GCC compiler and the command make.

These algorithms were built around the notion of *Mining Frequent ItemSets and associations rules*. This data analysis method was first introduced by *Agrawal et al. 1993* for mining transaction databases.


Its use and practices are close to those of PrefRec.
With SufRec, we propose mainly 3 uses, one of them allows a classic extraction of frequent itemsets. The second is a proposal to use the algorithm in parallel and simultaneously, using several machines. Finally, a last src folder is a proposal for mooving recursive applications. Some details are given in the readme of the src folder.

In this version, the init function is not the fastest (to retrieve the database), but is more efficient for high-volume databases.

The detail of the operation of these algorithms is specified in the corresponding scientific article (SufRec, an algorithm for mining association rules: Recursivity and task parallelism.  Expert    Systems with Applications, Volume 236, 121321  (2024).   https://doi.org/10.1016/j.eswa.2023.121321)

Available databases have been compressed. Information about these databases is given in the databases_test folder.



**src_Fim:**
* Mining frequent itemSets with a relative Support with 3 different recursive methods.
* Get an output Tsv file for information about frequent itemSets.
* Get 2 files necessary to launch Prefrules. 

**src_ParFim:**
* Mining frequent itemSets with a relative Support with 2 different methods.
* Intended for running with cluster.
* Get several output Tsv file for information about frequent itemSets.

**src_Mooving:**
* Mooving applications for extraction of frequent itemsets.


## Creating binaries and getting started
```
cd src/src_Fim && make
cd src/src_ParFim && make
cd src/src_Mooving && make

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
SufRec in src_Fim

./SufRec ../../sample/input/Fruits.txt d=, s=0.20 t=i o=d  \\ Do the extraction on the base "fruits" with a relativ minSup of 0.20 and the independant version of SufRec_Fim. Items are ordered by decreasing support.
./SufRec ../../sample/input/Fruits.txt d=, s=0.20 t=tf4 \\ Do the extraction on the base "fruits" with a relativ minSup of 0.20 and the multi thread version of SufRec_Fim. It uses fifo organisation and use 4 threads.
./SufRec ../../sample/input/Fruits.txt d=, s=0.20 t=tc4 \\ Only count frequent itemsets on the base "fruits" with a relativ minSup of 0.20 and the multi thread version of SufRec_Fim. It uses fifo organisation and use 4 threads.

SufRec in src_ParFim

./SufRec ../../sample/input/Fruits.txt d=, s=0.20 t=i m=5 x=1 o=d  \\ Do the extraction on the bases "fruits" with a relativ minSup of 0.20 and the independant version of SufRec_Fim. This extraction is supposed to be on the first machine, total of machine is 5. Items are ordered by decreasing support.

SufRec in src_Mooving

./SufRec ../../sample/input/Fruits.txt d=, s=0.20 t=i 5 3   \\ Do the extraction on the bases "fruits" with a relativ minSup of 0.20 and the independant version of SufRec_Fim. It starts by finding all the frequent itemsets with 5 items, then make 3 iterations of the mooving application.

```


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






