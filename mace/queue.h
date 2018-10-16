/* Library of queue: spped priority implementation 
            12/Apr/2001   by Takeaki Uno  e-mail:uno@nii.jp, 
    homepage:   http://research.nii.ac.jp/~uno/index.html  */
/* This program is available for only academic use, basically.
   Anyone can modify this program, but he/she has to write down 
    the change of the modification on the top of the source code.
   Neither contact nor appointment to Takeaki Uno is needed.
   If one wants to re-distribute this code, please
    refer the newest code, and show the link to homepage of 
    Takeaki Uno, to notify the news about the codes for the users. */

#ifndef _queue_h_
#define _queue_h_

#include"stdlib2.h"

#ifndef QUEUE_INT
 #ifdef QUEUE_INT_LONG
  #define QUEUE_INT LONG    // define the type before if change is needed 
  #define QUEUE_INTHUGE LONGHUGE    // comment out if QUEUE_INT is "short"
  #define QUEUE_INTF LONGF
 #else
  #define QUEUE_INT int    // define the type before if change is needed 
  #define QUEUE_INTHUGE INTHUGE    // comment out if QUEUE_INT is "short"
  #define QUEUE_INTF "%d"
 #endif
#endif

#ifndef QUEUE_ID
 #ifdef QUEUE_ID_LONG
  #define QUEUE_ID LONG    // define the type before if change is needed 
  #define QUEUE_IDHUGE LONGHUGE    // comment out if QUEUE_INT is "short"
  #define QUEUE_IDF LONGF
 #else
  #define QUEUE_ID int    // define the type before if change is needed 
  #define QUEUE_IDHUGE INTHUGE    // comment out if QUEUE_INT is "short"
  #define QUEUE_IDF "%d"
 #endif
#endif
#define SWAP_QUEUE_INT(a,b)  (common_QUEUE_INT=a,a=b,b=common_QUEUE_INT)
#define SWAP_QUEUE_ID(a,b)  (common_QUEUE_ID=a,a=b,b=common_QUEUE_ID)

typedef struct {
  unsigned char type;  // type of the structure
  QUEUE_INT *v;  // pointer to the array
  QUEUE_ID end;  // the length of the array
  QUEUE_ID t;  // end position+1
  QUEUE_ID s;  // start position
} QUEUE;

/* QUEUE stores at most end-1 elements. Overflow occurs after inserting end-1 elements */

#define QUEUE_INCREMENT(Q,i) ((i)=((i)>=(Q).end-1)?0:(i)+1)
#define QUEUE_DECREMENT(Q,i) ((i)=(i)==0?(Q).end-1:(i)-1)
#define QUEUE_LENGTH(Q) (((Q).t-(Q).s+(Q).end)%(Q).end)
#define QUEUE_LENGTH_(Q) ((Q).t-(Q).s)

/* macro for loop w.r.t., QUEUE */
#define QUEUE_F_LOOP(Q,i)  for((i)=(Q).s;(i)!=(Q).t;((i)=((i)>=(Q).end-1)?0:(i)+1))
#define QUEUE_F_LOOP_(Q,i)  for((i)=(Q).s;(i)<(Q).t;(i)++)
#define QUEUE_FE_LOOP(Q,i,x)  for((i)=(Q).s,x=(Q).v[i];(i)!=(Q).t;((i)=((i)>=(Q).end-1)?0:(i)+1),x=(Q).v[i])
#define QUEUE_FE_LOOP_(Q,i,x)  for((i)=(Q).s,x=(Q).v[i];(i)<(Q).t;(i)++,x=(Q).v[i])
#define QUEUE_B_LOOP(Q,i)  for((i)=(Q).t==0?(Q).end-1:(Q).t-1;(i)!=(Q).s;(i)=(i)==0?(Q).end-1:(i)-1)
#define QUEUE_B_LOOP_(Q,i)  for((i)=(Q).t-1;(i)>=(Q).s;(i)--)
#define QUEUE_BE_LOOP(Q,i,x)  for((i)=(Q).t==0?(Q).end-1:(Q).t-1,x=(Q).v[i];(i)!=(Q).s;(i)=(i)==0?(Q).end-1:(i)-1,x=(Q).v[i])
#define QUEUE_BE_LOOP_(Q,i,x)  for((i)=(Q).t-1;((i)>=(Q).s)?((x=(Q).v[i])||1):0;(i)--)

#define QUEUE_RMALL(Q) ((Q).t=(Q).s)
#define QUEUE_RMALL_(Q) ((Q).t=0)
#define QUEUE_HEAD(Q) ((Q).v[(Q).s])
#define QUEUE_TAIL_(Q) ((Q).v[(Q).t-1])

extern QUEUE INIT_QUEUE;
extern QUEUE_INT common_QUEUE_INT, *common_QUEUE_INTp;
QSORT_TYPE_HEADER(QUEUE_INT, QUEUE_INT)
QSORT_TYPE_HEADER(QUEUE_ID, QUEUE_ID)


/* initialization, not fill the memory by 0 */
void QUEUE_alloc (QUEUE *Q, QUEUE_ID siz);

/* termination processing */
void QUEUE_end (QUEUE *Q);

/* delivery: transpose that matrinx (transaction database) Q. Each row of the 
 transposed matrix is called occurrence.

variables to be set.
OQ:array for occurrences, c: for counting frequency, jump: list of items with non-empty OQ
if c!=NULL, count the frequency and set to c, and set occurrences to OQ, otherwise.
if jump==NULL, then the list of non-empty item will not be generated
Q:matrix, of an array of QUEUE, occ: list of rows of Q to be scaned, t; maximum ID of the
 row to be scaned; if occ==NULL, occ will be ignored, otherwise t will be ignored.
 M: end mark of each QUEUE. */
void QUEUE_delivery(QUEUE *OQ, VEC_ID *c, QUEUE *jump, QUEUE *Q, QUEUE *occ, VEC_ID t, QUEUE_INT M);
/* sort a QUEUE with WEIGHT, with already allocated memory (size have to no less than the size of QUEUE) */
void QUEUE_perm_WEIGHT (QUEUE *Q, WEIGHT *w, PERM *invperm, int flag);

/* remove (or unify) the consecutive same ID's in a QUEUE (duplication delete, if sorted) */
void QUEUE_rm_dup_WEIGHT (QUEUE *Q, WEIGHT *w);

/***********************************************************************/
/* duplicate occ's in jump, ( copy occ's to allocated QUEUE array) */
/* Q[i].end := original item, clear each original occ */
/* buffer size is multiplied by u */
/*******************************************************/


void QUEUE_occ_dup (QUEUE *jump, QUEUE **QQ, QUEUE *Q, WEIGHT **ww, WEIGHT *w, WEIGHT **ppw, WEIGHT *pw, int u);

/* return the position of the first element having value e. return -1 if no such element exists */
LONG QUEUE_ele (QUEUE *Q, QUEUE_INT e);

/* insert an element to the tail/head */
void QUEUE_ins_ (QUEUE *Q, QUEUE_INT e);
void QUEUE_ins (QUEUE *Q, QUEUE_INT e);
void QUEUE_ins_head_ (QUEUE *Q, QUEUE_INT e);
void QUEUE_ins_head (QUEUE *Q, QUEUE_INT e);

/* extract an element from the head/tail, without checking the underflow */
QUEUE_INT QUEUE_ext_ (QUEUE *Q);
QUEUE_INT QUEUE_ext (QUEUE *Q);
QUEUE_INT QUEUE_ext_tail_ (QUEUE *Q);
QUEUE_INT QUEUE_ext_tail (QUEUE *Q);

/* remove the j-th element and replace it by the tail/head or shift */
void QUEUE_rm_ (QUEUE *Q, QUEUE_ID j);
void QUEUE_rm (QUEUE *Q, QUEUE_ID j);
void QUEUE_rm_head_ (QUEUE *Q, QUEUE_ID j);
void QUEUE_rm_head (QUEUE *Q, QUEUE_ID j);
int QUEUE_rm_ele_ (QUEUE *Q, QUEUE_INT e);

/* Append Q2 to the tail of Q1. Q2 will (not) be deleted */
void QUEUE_append_ (QUEUE *Q1, QUEUE *Q2);
void QUEUE_append (QUEUE *Q1, QUEUE *Q2);
void QUEUE_concat_ (QUEUE *Q1, QUEUE *Q2);
void QUEUE_concat (QUEUE *Q1, QUEUE *Q2);

/* Append from j to jj th elements to the tail of Q1. Q2 will not be deleted */
void QUEUE_subconcat_ (QUEUE *Q1, QUEUE *Q2, QUEUE_ID j, QUEUE_ID jj);
void QUEUE_subconcat (QUEUE *Q1, QUEUE *Q2, QUEUE_ID j, QUEUE_ID jj);

/* initialize Q1 by length of Q2, and copy Q2 to Q1 */
void QUEUE_store_ (QUEUE *Q1, QUEUE *Q2);
void QUEUE_store (QUEUE *Q1, QUEUE *Q2);
/* copy Q2 to Q1 and delete Q2 */
void QUEUE_restore_ (QUEUE *Q1, QUEUE *Q2);
void QUEUE_restore (QUEUE *Q1, QUEUE *Q2);

/* copy Q2 to Q1 */
void QUEUE_cpy_ (QUEUE *Q1, QUEUE *Q2);
void QUEUE_cpy (QUEUE *Q1, QUEUE *Q2);
QUEUE QUEUE_dup_ (QUEUE *Q);
/* copy l elements of Q2 starting from s2 to the s1th position of Q1.
   size of Q1 is not increasing */
void QUEUE_subcpy_ (QUEUE *Q1, QUEUE_ID s1, QUEUE *Q2, QUEUE_ID s2, QUEUE_ID l);
void QUEUE_subcpy (QUEUE *Q1, QUEUE_ID s1, QUEUE *Q2, QUEUE_ID s2, QUEUE_ID l);

/* merge/minum/intersection of Q1 and Q2, and set Q1 to it.
 Both Q1 and Q2 have to be sorted in increasing order */
void QUEUE_merge_ (QUEUE *Q1, QUEUE *Q2);
void QUEUE_merge (QUEUE *Q1, QUEUE *Q2);
void QUEUE_minus_ (QUEUE *Q1, QUEUE *Q2);
void QUEUE_minus (QUEUE *Q1, QUEUE *Q2);
void QUEUE_and_ (QUEUE *Q1, QUEUE *Q2);
void QUEUE_and (QUEUE *Q1, QUEUE *Q2);
 
/* insertion sort */
void QUEUE_sort (QUEUE *Q);

  /* print */
void QUEUE_print (QUEUE *Q);
void QUEUE_print_ (QUEUE *Q);


#endif
