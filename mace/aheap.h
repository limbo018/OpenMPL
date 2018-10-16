/*
    array-based simple heap (fixed size)
            12/Apr/2001   by Takeaki Uno  e-mail:uno@nii.jp, 
    homepage:   http://research.nii.ac.jp/~uno/index.html  */
/* This program is available for only academic use, basically.
   Anyone can modify this program, but he/she has to write down 
    the change of the modification on the top of the source code.
   Neither contact nor appointment to Takeaki Uno is needed.
   If one wants to re-distribute this code, please
    refer the newest code, and show the link to homepage of 
    Takeaki Uno, to notify the news about the codes for the users. */

/* bench mark
  PentiumIII 500MHz, Memory 256MB Linux
  values 0-1,000,000 : set & del & get  1,000,000 times
  2.55sec

  # rotation  == 1/5 per 1 set/del

   *** simple array *** 
   value 0-1,000,000 set & set & set   1,000,000 times,
   0.88sec 
  */

#ifndef _aheap_h_
#define _aheap_h_

#include"stdlib2.h"

#ifndef AHEAP_KEY
 #ifdef AHEAP_KEY_DOUBLE
  #define AHEAP_KEY  double
  #define AHEAP_KEYHUGE DOUBLEHUGE
  #define AHEAP_KEYF "%f"
 #elif defined(AHEAP_KEY_WEIGHT)
  #define AHEAP_KEY  WEIGHT
  #define AHEAP_KEYHUGE WEIGHTHUGE
  #define AHEAP_KEYF WEIGHTF
 #else
#define AHEAP_KEY  int
  #define AHEAP_KEYHUGE INTHUGE
  #define AHEAP_KEYF "%d"
 #endif
#endif

#ifndef AHEAP_ID
 #define AHEAP_ID int
 #define AHEAP_ID_END INTHUGE
 #define AHEAP_IDF "%d"
#endif

#define AHEAP_IDX(H,i) (((i)+1-(H).base)%(H).end)
#define AHEAP_LEAF(H,i)   (((i)+(H).base)%(H).end+(H).end-1)
#define AHEAP_H(H,i)   (H).v[(((i)+(H).base)%(H).end+(H).end-1)]

typedef struct {
  unsigned char type;
  AHEAP_KEY *v;       /* array for heap key */
  int end;            /* the number of maximum elements */
  int base;           /* the constant for set 0 to the leftmost leaf */
} AHEAP;

QSORT_TYPE_HEADER (AHEAP_KEY, AHEAP_KEY)
QSORT_TYPE_HEADER (AHEAP_ID, AHEAP_ID)
extern AHEAP INIT_AHEAP;

/* initialization. allocate memory for H and fill it by +infinity */
void AHEAP_alloc (AHEAP *H, int num);
void AHEAP_end (AHEAP *H);

/* return the index of the leaf having the minimum key among the descendants
  of the given node i. If several leaves with the smallest key are there, 
  return the minimum index among them if f=0, maximum index if f=1, and
  random choice if f=2  */
AHEAP_ID AHEAP_findmin_node_ (AHEAP *H, AHEAP_ID i, int f);
AHEAP_ID AHEAP_findmin_node (AHEAP *H, AHEAP_ID i, int f);
AHEAP_ID AHEAP_findmin_head (AHEAP *H);
AHEAP_ID AHEAP_findmin_tail (AHEAP *H);
AHEAP_ID AHEAP_findmin_rnd (AHEAP *H);

/* return the index of the leaf having smaller value than a among the
  descendants of the given node i. If several leaves with the smallest key
  are there, return the minimum index among them if f=0, maximum index if f=1,
  and random choice if f=2  */
AHEAP_ID AHEAP_findlow_node (AHEAP *H, AHEAP_KEY a, AHEAP_ID i, int f);
AHEAP_ID AHEAP_findlow_head (AHEAP *H, AHEAP_KEY a);
AHEAP_ID AHEAP_findlow_tail (AHEAP *H, AHEAP_KEY a);
AHEAP_ID AHEAP_findlow_rnd (AHEAP *H, AHEAP_KEY a);

/* return the index of the leaf having smaller value than a next/previous to
  leaf i. return -1 if such a leaf does not exist  */
AHEAP_ID AHEAP_findlow_nxt (AHEAP *H, AHEAP_ID i, AHEAP_KEY a);
AHEAP_ID AHEAP_findlow_prv (AHEAP *H, AHEAP_ID i, AHEAP_KEY a);

/* change the key of node i to a /Add a to the key of node i, and update heap H */
void AHEAP_chg (AHEAP *H, AHEAP_ID i, AHEAP_KEY a);
void AHEAP_add (AHEAP *H, AHEAP_ID i, AHEAP_KEY a);

/* update the ancestor of node i */
void AHEAP_update (AHEAP *H, AHEAP_ID i);

/* find the leaf with the minimum key value among the leaves having index 
 smaller/larger than i, or between i and j */
AHEAP_ID AHEAP_upper_min (AHEAP *H, AHEAP_ID i);
AHEAP_ID AHEAP_lower_min (AHEAP *H, AHEAP_ID i);
AHEAP_ID AHEAP_interval_min (AHEAP *H, AHEAP_ID i, AHEAP_ID j);

/* print heap keys according to the structure of the heap */
void AHEAP_print (AHEAP *H);

#endif
