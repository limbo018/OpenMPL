/*  itemset search input/output common routines
            25/Nov/2007   by Takeaki Uno  e-mail:uno@nii.jp, 
    homepage:   http://research.nii.ac.jp/~uno/index.html  */
/* This program is available for only academic use, basically.
   Anyone can modify this program, but he/she has to write down 
    the change of the modification on the top of the source code.
   Neither contact nor appointment to Takeaki Uno is needed.
   If one wants to re-distribute this code, please
    refer the newest code, and show the link to homepage of 
    Takeaki Uno, to notify the news about the codes for the users. */

/* routines for itemset mining */

#ifndef _itemset_h_
#define _itemset_h_

#include"stdlib2.h"
#include"queue.h"
#define AHEAP_KEY_WEIGHT
#include"aheap.h"


typedef struct {
  int a;
  QUEUE itemset;   // current operating itemset
  QUEUE add;       // for equisupport (hypercube decomposition)
  int ub, lb;   // upper/lower bounds for the itemset size
  WEIGHT frq, pfrq, frq_ub, frq_lb;  // upper/lower bounds for the frequency
  WEIGHT rposi_lb, rposi_ub, posi_lb, posi_ub, nega_ub, nega_lb;  // upper/lower bounds for the sum of positive/negative weights
  WEIGHT setrule_lb;  // frequency lower bound for set rule
  double ratio, prob;   // confidence and independent probability of the current pattern
  double ratio_ub, ratio_lb, prob_ub, prob_lb;   // upper/lower bounds for confidence and independent probability
  QUEUE_INT target;  // target item for rule mining
  char *itemflag;       // 1 if it is include in the pattern (and 2 if included in add)
  WEIGHT *item_frq;    // frequency of each item
  WEIGHT total_weight;  // total weight of the input database
  int len_ub, len_lb;   // upper/lower bounds for the length of the pattern
  int gap_ub, gap_lb;   // upper/lower bounds for the gaps in the pattern
  LONG *sc;    // #itemsets classified by the sizes
  QUEUE_INT item_max, item_max_org;  // (original) maximum item
  AHEAP topk;  // heap for topk mining. valid if topk->h is not NULL
  int flag;    // flag for various functions
  PERM *perm;   // permutation array for output itemset: item => original item
  FILE *fp;    // file pointer to the output file
  char separator; // separator of items output
  int progress;
  LONG iters, iters2, iters3;  //iterations
  LONG solutions, solutions2;  // number of solutions output
  LONG outputs, outputs2;    // #calls of ITEMSET_output_itemset or ITEMSET_solusion
  LONG max_solutions; // maximum solutions to be output
  void *X;  // pointer to the original data
  int dir;  // direction flag for AGRAPH & SGRAPH

  int multi_core;  // number of processors
  LONG *multi_iters, *multi_iters2, *multi_iters3;  //iterations
  LONG *multi_solutions, *multi_solutions2;  // number of solutions output
  LONG *multi_outputs, *multi_outputs2;    // #calls of ITEMSET_output_itemset or ITEMSET_solusion
  FILE2 *multi_fp;  // output file2 pointer for multi-core mode
  WEIGHT *set_weight;  // the frequency of each prefix of current itemset
  QUEUE **set_occ;    // the occurrence of each prefix of current itemset

#ifdef MULTI_CORE
  pthread_spinlock_t lock_counter;   // couneter locker for jump counter
  pthread_spinlock_t lock_sc;   // couneter locker for score counter
  pthread_spinlock_t lock_output;   // couneter locker for #output 
#endif
} ITEMSET;

/* parameters for ITEMSET.flag */

#define ITEMSET_ITERS2 4  // output #iters2
#define ITEMSET_PRE_FREQ 8  // output frequency preceding to each itemset
#define ITEMSET_FREQ 16  // output frequency following to each itemset
#define ITEMSET_ALL 32 // concat all combinations of "add" to each itemset

#define ITEMSET_TRSACT_ID 64  // output transaction ID's in occurrences
#define ITEMSET_OUTPUT_EDGE 128  // output itemset as edge set (refer AGRAPH)
#define ITEMSET_IGNORE_BOUND 256 // ignore constraint for frequency
#define ITEMSET_RM_DUP_TRSACT 512 // remove duplicated transaction ID's
#define ITEMSET_MULTI_OCC_PRINT 1024 //print each component of occ
   // TRSACT_ID+MULTI_OCC_PRINT means print first two components of occ
#define ITEMSET_NOT_ITEMSET  2048 // do not print itemset to the output file
#define ITEMSET_RULE_SUPP  4096 // output confidence and item frquency by abusolute value
#define ITEMSET_OUTPUT_POSINEGA  8192 // output negative/positive frequencies
#define ITEMSET_MULTI_OUTPUT 16384 // for multi-core mode
#define ITEMSET_USE_ORG 32768 // use item_max_org to the size of use
#define ITEMSET_ITEMFRQ 65536 // allocate item_frq
#define ITEMSET_ADD 131072    // allocate add

#define ITEMSET_RULE_FRQ 262144
#define ITEMSET_RULE_INFRQ 524288
#define ITEMSET_RULE_RFRQ 1048576
#define ITEMSET_RULE_RINFRQ 2097152
#define ITEMSET_RFRQ 4194304
#define ITEMSET_RINFRQ 8388608
#define ITEMSET_POSI_RATIO 16777216
#define ITEMSET_SET_RULE 134217728

#define ITEMSET_APPEND 268435456   // append the output to the fiile
#define ITEMSET_RULE_ADD 536870912   // append items in add to the solution, for rule output

//#define ITEMSET_RULE (ITEMSET_RULE_FRQ + ITEMSET_RULE_INFRQ + ITEMSET_RULE_RFRQ + ITEMSET_RULE_RINFRQ + ITEMSET_RFRQ + ITEMSET_RINFRQ + ITEMSET_SET_RULE)  // for check any rule is true
#define ITEMSET_RULE (ITEMSET_RULE_FRQ + ITEMSET_RULE_INFRQ + ITEMSET_RULE_RFRQ + ITEMSET_RULE_RINFRQ + ITEMSET_SET_RULE)  // for check any rule is true

#ifndef ITEMSET_INTERVAL
#define ITEMSET_INTERVAL 500000
#endif

/* Output information about ITEMSET structure. flag&1: print frequency constraint */
void ITEMSET_print (ITEMSET *II, int flag);

/* topk.end>0 => initialize heap for topk mining */
/* all pointers will be set to 0, but not for */
/* if topK mining, set topk.end to "K" */
void ITEMSET_init (ITEMSET *I);
void ITEMSET_alloc (ITEMSET *I, char *fname, PERM *perm, QUEUE_INT item_max, size_t item_max_org);
void ITEMSET_end (ITEMSET *I);

/* sum the counters computed by each thread */
void ITEMSET_merge_counters (ITEMSET *I);

/*******************************************************************/
/* output at the termination of the algorithm */
/* print #of itemsets of size k, for each k */
/*******************************************************************/
void ITEMSET_last_output (ITEMSET *I);

/* output frequency, coverage */
void ITEMSET_output_frequency (ITEMSET *I, int core_id);

/* output an itemset to the output file */
void ITEMSET_output_itemset (ITEMSET *I, QUEUE *occ, int core_id);

/* output itemsets with adding all combination of "add"
   at the first call, i has to be "add->t" */
void ITEMSET_solution (ITEMSET *I, QUEUE *occ, int core_id);

/*************************************************************************/
/* ourput a rule */
/*************************************************************************/
void ITEMSET_output_rule (ITEMSET *I, QUEUE *occ, double p1, double p2, size_t item, int core_id);

/*************************************************************************/
/* check all rules for a pair of itemset and item */
/*************************************************************************/
void ITEMSET_check_rule (ITEMSET *I, WEIGHT *w, QUEUE *occ, size_t item, int core_id);

/*************************************************************************/
/* check all rules for an itemset and all items */
/*************************************************************************/
void ITEMSET_check_all_rule (ITEMSET *I, WEIGHT *w, QUEUE *occ, QUEUE *jump, WEIGHT total, int core_id);

#endif




