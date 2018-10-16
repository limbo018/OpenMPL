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

#ifndef _queue_c_
#define _queue_c_


#include"queue.h"
#include"stdlib2.c"

QSORT_TYPE(QUEUE_INT, QUEUE_INT)
QSORT_TYPE(QUEUE_ID, QUEUE_ID)
QUEUE INIT_QUEUE = {TYPE_QUEUE,NULL,0,0,0};
QUEUE_INT *common_QUEUE_INTp;

/* initialization, not fill the memory by 0 */
void QUEUE_alloc (QUEUE *Q, QUEUE_ID siz){
  *Q = INIT_QUEUE;
  Q->end = siz+1;
  malloc2 (Q->v, siz+1, EXIT);
}

/* termination processing */
void QUEUE_end (QUEUE *Q){
  free2 (Q->v);
  *Q = INIT_QUEUE;
}

/* tranpose the matrix ; counting/transpose/memory_allocate */
void QUEUE_delivery(QUEUE *OQ, VEC_ID *c, QUEUE *jump, QUEUE *Q, QUEUE *occ, VEC_ID t, QUEUE_INT M){  VEC_ID i, e;
  QUEUE_INT *x;
  FLOOP(i, 0, occ? occ->t: t){
    e = occ? occ->v[i]: i;
    if ( c ){
      if ( jump ){ MLOOP (x, Q[e].v, M){ if ( c[*x]==0 ) QUE_INS (*jump, *x); c[*x]++; }
      } else { MLOOP (x, Q[e].v, M){ c[*x]++; }}
    } else {
      if ( jump ){ MLOOP (x, Q[e].v, M){ if ( OQ[*x].t==0 ) QUE_INS (*jump, *x); QUE_INS (OQ[*x], e); }
      } else MLOOP (x, Q[e].v, M){ QUE_INS (OQ[*x], e); }
    }
  }
}

/* sort a QUEUE with WEIGHT, with already allocated memory */
void QUEUE_perm_WEIGHT (QUEUE *Q, WEIGHT *w, PERM *invperm, int flag){
  WEIGHT y;
  if ( w ){
    qsort_perm__QUEUE_INT (Q->v, Q->t, invperm, flag);
    ARY_INVPERMUTE_ (w, invperm, y, Q->t);
  }
  qsort_QUEUE_INT (Q->v, Q->t, flag);
}

/* remove (or unify) the consecutive same ID's in a QUEUE (duplication delete, if sorted) */
void QUEUE_rm_dup_WEIGHT (QUEUE *Q, WEIGHT *w){
  VEC_ID j, jj=0;
  if ( w ){
    FLOOP (j, 1, Q->t){
      if ( Q->v[j-1] != Q->v[j] ){
        Q->v[++jj] = Q->v[j];
        w[jj] = w[j];
      } else w[jj] += w[j];
    }
  } else FLOOP (j, 1, Q->t){
    if ( Q->v[j-1] != Q->v[j] ) Q->v[++jj] = Q->v[j];
  }
  if ( Q->t>0 ) Q->t = jj+1;
}

/***********************************************************************/
/* duplicate occ's in jump, ( copy occ's to allocated QUEUE array) */
/* Q[i].end := original item, clear each original occ */
/* buffer size is multiplied by u */
/*******************************************************/
void QUEUE_occ_dup (QUEUE *jump, QUEUE **QQ, QUEUE *Q, WEIGHT **ww, WEIGHT *w, WEIGHT **ppw, WEIGHT *pw, int u){
  QUEUE_ID i, l=QUEUE_LENGTH_(*jump);
  size_t cnt=0;
  QUEUE_INT e, *x;
  char *buf;
  int unit = sizeof(*Q) + (w?sizeof(*w):0) + (pw?sizeof(*pw):0);
 
  ENMAX (u, sizeof(*x));
  MQUE_FLOOP (*jump, x) cnt += Q[*x].t;
  if ( cnt == 0 ){ *QQ=NULL; return; }
  malloc2 (buf, l*unit + (cnt+l)*u, EXIT);
  *QQ = (QUEUE*)buf; buf += sizeof(*Q) *l;
  if ( w ){ *ww = (WEIGHT *)buf; buf += sizeof(*w)*l; }
  if ( pw ){ *ppw = (WEIGHT *)buf; buf += sizeof(*pw)*l; }
  for (i=0 ; i<jump->t ; i++){
    e = jump->v[i];
    (*QQ)[i].end = e;
    (*QQ)[i].v = (QUEUE_INT *)buf;
    (*QQ)[i].t = Q[e].t;
    memcpy (buf, Q[e].v, (Q[e].t+1)*u);
    buf += (Q[e].t+1) *u;
    if ( w ) (*ww)[i] = w[e];
    if ( pw ) (*ppw)[i] = pw[e];
  }
}


/* return the position of the first element having value e. return -1 if no such element exists */
LONG QUEUE_ele (QUEUE *Q, QUEUE_INT e){
  QUEUE_INT *x;
  MQUE_FLOOP (*Q, x)
    if ( *x == e ) return (x - Q->v);
  return (-1);
}

/* insert an element to the tail */
void QUEUE_ins_ (QUEUE *Q, QUEUE_INT e){
  Q->v[Q->t] = e;
  Q->t++;
}
void QUEUE_ins (QUEUE *Q, QUEUE_INT e){
  Q->v[Q->t] = e;
  QUEUE_INCREMENT (*Q, Q->t);
  if (Q->s == Q->t ) error_num ("QUEUE_ins: overflow", Q->s, EXIT);
}

/* insert an element to the head */
void QUEUE_ins_head_ (QUEUE *Q, QUEUE_INT e){
  Q->s--;
  Q->v[Q->s] = e;
}
void QUEUE_ins_head (QUEUE *Q, QUEUE_INT e){
  QUEUE_DECREMENT(*Q,Q->s);
  Q->v[Q->s] = e;
  if (Q->s == Q->t ) error_num ("QUEUE_ins_head: underflow", Q->s, EXIT);
}

/* extract an element from the head, without checking underflow */
QUEUE_INT QUEUE_ext_ (QUEUE *Q){
  (Q->s)++;
  return (Q->v[Q->s-1]);
}
QUEUE_INT QUEUE_ext (QUEUE *Q){
  QUEUE_INT e;
  if (Q->s == Q->t ) error_num ("QUEUE_ext: empty queue", Q->s, EXIT0);
  e = Q->v[Q->s];
  QUEUE_INCREMENT(*Q,Q->s);
  return ( e);
}

/* extract an element from the tail, without checking underflow */
QUEUE_INT QUEUE_ext_tail_ (QUEUE *Q){
  (Q->t)--;
  return (Q->v[Q->t]);
}
QUEUE_INT QUEUE_ext_tail (QUEUE *Q){
  if ( Q->s == Q->t ) error_num ("QUEUE_ext_tail: empty queue", Q->s, EXIT0);
  QUEUE_DECREMENT(*Q,Q->t);
  return (Q->v[Q->t]);
}

/* remove the j-th element and replace it by the tail */
void QUEUE_rm_ (QUEUE *Q, QUEUE_ID j){
  Q->t--;
  Q->v[j] = Q->v[Q->t];
}
void QUEUE_rm (QUEUE *Q, QUEUE_ID j){
  if ( Q->s <= Q->t ){
    if ( j < Q->s || j >= Q->t ) error ("QUEUE_rm: j is out of queue", EXIT);
  } else if ( j < Q->s && j >= Q->t ) error ("QUEUE_rm: j is out of queue", EXIT);
  QUEUE_DECREMENT(*Q,Q->t);
  Q->v[j] = Q->v[Q->t];
}

/* remove the j-th element and replace it by the head */
void QUEUE_rm_head_ (QUEUE *Q, QUEUE_ID j){
  Q->v[j] = Q->v[Q->s];
  Q->s++;
}
void QUEUE_rm_head (QUEUE *Q, QUEUE_ID j){
  if ( Q->s <= Q->t ){
    if ( j < Q->s || j >= Q->t ) error ("QUEUE_rm: j is out of queue", EXIT);
  } else if ( j < Q->s && j >= Q->t ) error ("QUEUE_rm: j is out of queue", EXIT);
  Q->v[j] = Q->v[Q->s];
  QUEUE_INCREMENT(*Q,Q->s);
}

/* remove the j-th element and shift the following elements to fill the gap */
int QUEUE_rm_ele_ (QUEUE *Q, QUEUE_INT e){
  QUEUE_ID i;
  QUEUE_F_LOOP (*Q, i){
    if ( Q->v[i] == e ){
      memcpy ( &(Q->v[i]), &(Q->v[i+1]), (Q->t-i-1)*sizeof(QUEUE_INT));
      Q->t--;
      return (1);
    }
  }
  return (0);
}  
/* insert e to the position determined by the increasing order of elements */
void QUEUE_ins_ele_ (QUEUE *Q, QUEUE_INT e){
  QUEUE_ID i;
  QUEUE_INT ee;
  QUEUE_BE_LOOP_ (*Q, i, ee){
    if ( ee<e ) break;
    Q->v[i+1] = ee;
  }
  Q->v[i+1] = e;
  Q->t++;
}

/* Append Q2 to the tail of Q1. Q2 will not be deleted */
void QUEUE_concat_ (QUEUE *Q1, QUEUE *Q2){
  memcpy ( &(Q1->v[Q1->t]), &(Q2->v[Q2->s]), (Q2->t-Q2->s)*sizeof(QUEUE_INT));
  Q1->t += Q2->t-Q2->s;
}
void QUEUE_concat (QUEUE *Q1, QUEUE *Q2){
  QUEUE_ID e = Q2->s;
  while ( e != Q2->t){
    QUEUE_ins (Q1, Q2->v[e]);
    QUEUE_INCREMENT(*Q2,e);
  }
}
/* Append Q2 to the tail of Q1. Q2 will be deleted */
void QUEUE_append_ (QUEUE *Q1, QUEUE *Q2){
  QUEUE_concat_ (Q1, Q2);
  QUEUE_RMALL (*Q2);
}
void QUEUE_append (QUEUE *Q1, QUEUE *Q2){ // more improvement can be
  while ( Q2->s != Q2->t )
      QUEUE_ins (Q1, QUEUE_ext(Q2));
}

/* Append from j to jj th elements to the tail of Q1. Q2 will not be deleted */
void QUEUE_subconcat_ (QUEUE *Q1, QUEUE *Q2, QUEUE_ID j, QUEUE_ID jj){
  for ( ; j<=jj ; j++){
    Q1->v[Q1->t] = Q2->v[j];
    Q1->t++;
  }
}
void QUEUE_subconcat (QUEUE *Q1, QUEUE *Q2, QUEUE_ID j, QUEUE_ID jj){
  while (1){
    QUEUE_ins (Q1, Q2->v[j]);
    if ( j == jj ) break;
    QUEUE_INCREMENT(*Q2,j);
  } 
}

/* initialize Q1 by length of Q2, and copy Q2 to Q1 */
void QUEUE_store_ (QUEUE *Q1, QUEUE *Q2){
  QUEUE_alloc (Q1, QUEUE_LENGTH(*Q2));
  QUEUE_concat_ (Q1, Q2);
}
void QUEUE_store (QUEUE *Q1, QUEUE *Q2){
  QUEUE_alloc (Q1, QUEUE_LENGTH(*Q2));
  QUEUE_concat (Q1, Q2);
}
/* copy Q2 to Q1 and delete Q2 */
void QUEUE_restore_ (QUEUE *Q1, QUEUE *Q2){
  QUEUE_RMALL (*Q1);
  QUEUE_concat_ (Q1, Q2);
  QUEUE_end (Q2);
}
void QUEUE_restore (QUEUE *Q1, QUEUE *Q2){
  QUEUE_RMALL (*Q1);
  QUEUE_concat (Q1, Q2);
  QUEUE_end (Q2);
}

/* copy Q2 to Q1 */
void QUEUE_cpy_ (QUEUE *Q1, QUEUE *Q2){
  Q1->s = Q1->t = 0;
  QUEUE_concat_ (Q1, Q2);
}
void QUEUE_cpy (QUEUE *Q1, QUEUE *Q2){
  QUEUE_RMALL (*Q1);
  QUEUE_concat (Q1, Q2);
}

/* compare two queues */
int QUEUE_cmp_ (QUEUE *Q1, QUEUE *Q2){
  QUEUE_INT *x, *y=Q2->v;
  MQUE_FLOOP (*Q1, x){
    if ( *x != *y ) return (0);
    y++;
  }
  return (1);
}

/* copy l elements of Q2 starting from s2 to the s1th position of Q1.
   size of Q1 is not increasing */
void QUEUE_subcpy_ (QUEUE *Q1, QUEUE_ID s1, QUEUE *Q2, QUEUE_ID s2, QUEUE_ID l){
  memcpy ( &(Q1->v[s1]), &(Q2->v[s2]), (l-s2)*sizeof(QUEUE_INT));
}
void QUEUE_subcpy (QUEUE *Q1, QUEUE_ID s1, QUEUE *Q2, QUEUE_ID s2, QUEUE_ID l){
  for ( ; s2!=l ; QUEUE_INCREMENT(*Q1,s1),QUEUE_INCREMENT(*Q2,s2) )
      Q1->v[s1] = Q2->v[s2];
  Q1->v[s1] = Q2->v[s2];
}

/* duplicate Q2 to Q1. The memory size will be the length of Q2 */
QUEUE QUEUE_dup_ (QUEUE *Q){
  QUEUE QQ;
  QUEUE_alloc (&QQ, MAX(Q->t+1, Q->end-1));
  QUEUE_cpy_ (&QQ, Q);
  return (QQ);
}

/* merge Q1 and Q2 by insert all elements of Q2 to Q1 with deleting duplications. Both Q1 and Q2 have to be sorted in increasing order */
void QUEUE_merge_ (QUEUE *Q1, QUEUE *Q2){
  QUEUE_ID i=Q1->t-1, j=Q2->t-1, t=i+j-Q2->s+1;
  QUEUE_INT ei, ej;
  if ( i+1 == Q1->s || j+1 == Q2->s ){
    QUEUE_concat_ (Q1, Q2);
    return;
  }
  Q1->t = t+1;
  ei = Q1->v[i];
  ej = Q2->v[j];
  while (1){
    if ( ei > ej ){
      Q1->v[t] = ei;
      if ( i == Q1->s ){
        QUEUE_subcpy_ (Q1, Q1->s, Q2, Q2->s, j);
        return;
      }
      i--;
      ei = Q1->v[i];
    } else {
      Q1->v[t] = ej;
      if ( j == Q2->s ) return;
      j--;
      ej = Q2->v[j];
    }
    t--;
  }
}
void QUEUE_merge (QUEUE *Q1, QUEUE *Q2){
  QUEUE_ID i=Q1->t, j=Q2->t;
  QUEUE_INT ei, ej;
  QUEUE_ID t = (Q1->t + QUEUE_LENGTH(*Q2)-1) % Q1->end;
  if ( QUEUE_LENGTH(*Q1) + QUEUE_LENGTH(*Q2) >= Q1->end ){
    print_err ("QUEUE_merge: overflow Q1->end="QUEUE_INTF", Q1length="QUEUE_INTF", Q2length="QUEUE_INTF"\n", Q1->end, QUEUE_LENGTH(*Q1), QUEUE_LENGTH(*Q2));
    exit (1);
  }
  if ( i == Q1->s || j == Q2->s ){
    QUEUE_concat (Q1, Q2);
    return;
  }

  Q1->t = t;
  QUEUE_DECREMENT(*Q1,i);
  QUEUE_DECREMENT(*Q2,j);
  ei = Q1->v[i];
  ej = Q2->v[j];
  while (1){
    if ( ei > ej ){
      Q1->v[t] = ei;
      if ( i == Q1->s ){
        QUEUE_subcpy (Q1, Q1->s, Q2, Q2->s, (j+Q2->end-Q2->s)%Q2->end);
        return;
      }
      QUEUE_DECREMENT(*Q1,i);
      ei = Q1->v[i];
    } else {  
      Q1->v[t] = ej;
      if ( j == Q2->s ) return;
      QUEUE_DECREMENT(*Q2,j);
      ej = Q2->v[j];
    }
    QUEUE_DECREMENT(*Q1,t);
  }
}

/* delete all elements of Q1 included in Q2.
 both Q1 and Q2 have to be sorted in increasing order */
void QUEUE_minus_ (QUEUE *Q1, QUEUE *Q2){
  QUEUE_ID i=Q1->s, i2 = Q2->s, ii=Q1->s;
  while ( i != Q1->t && i2 != Q2->t){
    if (Q1->v[i] > Q2->v[i2] ) i2++;
    else {
      if (Q1->v[i] < Q2->v[i2] ){
        Q1->v[ii] = Q1->v[i];
        ii++;
      }
      i++;
    }
  }
  while ( i != Q1->t ){
    Q1->v[ii] = Q1->v[i];
    i++;
    ii++;
  }
  Q1->t = ii;
}
void QUEUE_minus (QUEUE *Q1, QUEUE *Q2){
  QUEUE_ID i=Q1->s, i2 = Q2->s, ii=Q1->s;
  while ( i != Q1->t && i2 != Q2->t ){
    if ( Q1->v[i] > Q2->v[i2] ) QUEUE_INCREMENT (*Q2, i2);
    else {
      if ( Q1->v[i] < Q2->v[i2] ){
        Q1->v[ii] = Q1->v[i];
        QUEUE_INCREMENT (*Q1, ii);
      }
      QUEUE_INCREMENT (*Q1, i);
    }
  }
  while ( i != Q1->t ){
    Q1->v[ii] = Q1->v[i];
    QUEUE_INCREMENT (*Q1, i);
    QUEUE_INCREMENT (*Q1, ii);
  }
  Q1->t = ii;
}

/* Delete all elements of Q1 which are not included in Q2. 
 both Q1 and Q2 have to be sorted in increasing order */
QUEUE_ID QUEUE_intsec_ (QUEUE *Q1, QUEUE *Q2){
  QUEUE_ID i=Q1->s, i2 = Q2->s, c=0;
  while ( i != Q1->t ){
    if ( Q1->v[i] > Q2->v[i2] ){
      if ( ++i2 == Q2->t ) break;
    } else {
      if ( Q1->v[i] == Q2->v[i2] ) c++;
      i++;
    }
  }
  return (c);
}
void QUEUE_and_ (QUEUE *Q1, QUEUE *Q2){
  QUEUE_ID i=Q1->s, i2 = Q2->s, ii=Q1->s;
  while ( i != Q1->t ){
    if ( Q1->v[i] > Q2->v[i2] ){
      if ( ++i2 == Q2->t ) break;
    } else {
      if ( Q1->v[i] == Q2->v[i2] ) Q1->v[ii++] = Q1->v[i];
      i++;
    }
  }
  Q1->t = ii;
}
void QUEUE_and (QUEUE *Q1, QUEUE *Q2){
  QUEUE_ID i=Q1->s, i2 = Q2->s, ii=Q1->s;
  while ( i != Q1->t && i2 != Q2->t){
    if ( Q1->v[i] > Q2->v[i2] ) QUEUE_INCREMENT (*Q2, i2);
    else {
      if ( Q1->v[i] == Q2->v[i2] ){
        Q1->v[ii] = Q1->v[i];
        QUEUE_INCREMENT (*Q1, ii);
      }
      QUEUE_INCREMENT (*Q1, i);
    }
  }
  Q1->t = ii;
}

/* insertion sort */
void QUEUE_sort (QUEUE *Q){
  QUEUE_ID i = Q->s, j, jj;
  QUEUE_INT e;
  if ( i== Q->t ) return;
  QUEUE_INCREMENT(*Q,i);
  for ( ; i!=Q->t ; QUEUE_INCREMENT(*Q,i) ){
    e=Q->v[i]; 
    j=i; 
    while (1){
      jj = j;
      QUEUE_DECREMENT(*Q,j);
      if ( Q->v[j] <= e ) { Q->v[jj] = e; break; }
      Q->v[jj] = Q->v[j];
      if ( j == Q->s) { Q->v[j] = e; break; }
    }
  }
}


/* print a queue */
void QUEUE_print (QUEUE *Q){
  QUEUE_ID i;
  for ( i=Q->s ; i!=Q->t ; ){
    printf (QUEUE_INTF" ", Q->v[i]);
    QUEUE_INCREMENT(*Q,i);
  }
  printf ("\n");
}
/* permutation version */
void QUEUE_perm_print (QUEUE *Q, QUEUE_ID *q){
  QUEUE_ID i;
  for ( i=Q->s ; i!=Q->t ; ){
    printf (QUEUE_INTF" ", q[Q->v[i]]);
    QUEUE_INCREMENT(*Q,i);
  }
  printf ("\n");
}
void QUEUE_printn (QUEUE *Q){
  QUEUE_ID i;
  for ( i=Q->s ; i!=Q->t ; ){
    printf (QUEUE_INTF" ", Q->v[i]);
    QUEUE_INCREMENT(*Q,i);
  }
}
void QUEUE_perm_printn (QUEUE *Q, QUEUE_ID *q){
  QUEUE_ID i;
  for ( i=Q->s ; i!=Q->t ; ){
    printf (QUEUE_INTF" ",q[Q->v[i]]);
    QUEUE_INCREMENT(*Q,i);
  }
}
void QUEUE_print_ (QUEUE *Q){
  QUEUE_ID i;
  printf("s="QUEUE_IDF",t="QUEUE_INTF": ", Q->s, Q->t);
  for ( i=Q->s ; i!=Q->t ; ){
    printf (QUEUE_INTF" ",Q->v[i]);
    QUEUE_INCREMENT(*Q,i);
  }
  printf ("\n");
}

void QUEUE_print__ (QUEUE *Q){
  QUEUE_ID i;
  printf("s="QUEUE_IDF",t="QUEUE_IDF": ", Q->s, Q->t);
  for ( i=Q->s ; i!=Q->t ; i++ ) printf (QUEUE_INTF" ",Q->v[i]);
  printf ("\n");
}

#endif
