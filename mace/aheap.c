/*
    array-based simple heap (fixex size)
            12/Apr/2001   by Takeaki Uno  e-mail:uno@nii.jp, 
    homepage:   http://research.nii.ac.jp/~uno/index.html  */
/* This program is available for only academic use, basically.
   Anyone can modify this program, but he/she has to write down 
    the change of the modification on the top of the source code.
   Neither contact nor appointment to Takeaki Uno is needed.
   If one wants to re-distribute this code, please
    refer the newest code, and show the link to homepage of 
    Takeaki Uno, to notify the news about the codes for the users. */

#ifndef _aheap_c_
#define _aheap_c_

#include"aheap.h"

QSORT_TYPE (AHEAP_KEY, AHEAP_KEY)
QSORT_TYPE (AHEAP_ID, AHEAP_ID)
AHEAP INIT_AHEAP = {TYPE_AHEAP,NULL,0,0};

/* allocate memory */
void AHEAP_alloc (AHEAP *H, AHEAP_ID num){
  AHEAP_ID i;
#ifdef ERROR_CHECK
  if ( num<0 ) error_num ("size is out of range", num, EXIT);
#endif
  *H = INIT_AHEAP;
  if ( num>0 ) malloc2 (H->v, num*2, EXIT);
  H->end = num;
  ARY_FILL (H->v, 0, num*2, AHEAP_KEYHUGE);
  for (i=0 ; i<num-1 ; i=i*2+1);
  H->base = i - num + 1;
}

/* termination */
void AHEAP_end (AHEAP *H){
  free2 (H->v);
  *H = INIT_AHEAP;
}

/* return the index of the leaf having the minimum key among the descendants
  of the given node i. If several leaves with the smallest key are there, 
  return the minimum index among them if f=0, maximum index if f=1, and
  random choice if f=2  */
/* random choice version. choose one child to go down randomly for each node,
   thus it is not uniformly random */
/* node_ returns the ID of leaf */
AHEAP_ID AHEAP_findmin_node_ (AHEAP *H, AHEAP_ID i, int f){
  while ( i < H->end-1 ){
    if ( H->v[i*2+1] == H->v[i] )
      if ( H->v[i*2+2] == H->v[i] )
        if ( f == 2 ) i = i*2 + 1 + rand()%2;
        else i = i*2+1+f;
      else i = i*2+1;
    else i = i*2+2;
  }
  return (i);
}
AHEAP_ID AHEAP_findmin_node (AHEAP *H, AHEAP_ID i, int f){
  if ( H->end <= 0 ) return (-1);
  return (AHEAP_IDX(*H, AHEAP_findmin_node_ (H, i, f)));
}
AHEAP_ID AHEAP_findmin_head (AHEAP *H){ return (AHEAP_findmin_node (H, 0, 0) ); }
AHEAP_ID AHEAP_findmin_tail (AHEAP *H){ return (AHEAP_findmin_node (H, 0, 1) ); }
AHEAP_ID AHEAP_findmin_rnd (AHEAP *H){ return (AHEAP_findmin_node (H, 0, 2) ); }

/* return the index of the leaf having smaller value than a among the
  descendants of the given node i. If several leaves with the smallest key
  are there, return the minimum index among them if f=0, maximum index if f=1,
  and random choice if f=2  */
AHEAP_ID AHEAP_findlow_node (AHEAP *H, AHEAP_KEY a, AHEAP_ID i, int f){
  if ( H->end == 0 ) return (-1); 
  if ( H->v[0] > a ) return (-1);
  while ( i < H->end-1 ){
    if ( f == 2 ) {
      if ( H->v[i*2+1] <= a )
          if ( H->v[i*2+2] <= a ) i = i*2 + 1 + rand()%2;
          else i = i*2+1;
      else i = i*2+2;
    } else if ( H->v[i*2+1] <= a ) i = i*2+1+f; else i = i*2+2-f;
  }
  return (AHEAP_IDX(*H, i) );
}
AHEAP_ID AHEAP_findlow_head (AHEAP *H, AHEAP_KEY a){ return (AHEAP_findlow_node (H, a, 0, 0) ); }
AHEAP_ID AHEAP_findlow_tail (AHEAP *H, AHEAP_KEY a){ return (AHEAP_findlow_node (H, a, 0, 1) ); }
AHEAP_ID AHEAP_findlow_rnd (AHEAP *H, AHEAP_KEY a){ return (AHEAP_findlow_node (H, a, 0, 2) ); }

/* return the index of the leaf having smaller value than a next/previous to
  leaf i. return -1 if such a leaf does not exist  */
AHEAP_ID AHEAP_findlow_nxt (AHEAP *H, AHEAP_ID i, AHEAP_KEY a){
  if ( H->end == 0 ) return (-1);
  if ( i<0 || i>= H->end ) return ( AHEAP_findlow_head (H, a));
  for (i=AHEAP_LEAF(*H,i); i>0 ; i=(i-1)/2){
     /* i is the child of smaller index, and the key of the sibling of i is less than a */
    if ( i%2 == 1 && H->v[i+1] <= a ) return (AHEAP_findlow_node (H, a, i+1, 0) );
  }
  return (-1);
}
AHEAP_ID AHEAP_findlow_prv (AHEAP *H, AHEAP_ID i, AHEAP_KEY a){
  if ( H->end == 0 ) return (-1); 
  if ( i<0 || i>= H->end ) return ( AHEAP_findlow_head (H, a));
  for (i=AHEAP_LEAF(*H,i); i>0 ; i=(i-1)/2){
     /* i is the child of larger index, and the key of the sibling of i is less than a */
    if ( i%2 == 0 && H->v[i-1] <= a ) return (AHEAP_findlow_node (H, a, i-1, 1) );
  }
  return (-1);
}

/* change the key of node i to a /Add a to the key of node i, and update heap H */
void AHEAP_chg (AHEAP *H, AHEAP_ID i, AHEAP_KEY a){
  i = AHEAP_LEAF (*H, i);
  H->v[i] = a;
  AHEAP_update (H, i);
}
void AHEAP_add (AHEAP *H, AHEAP_ID i, AHEAP_KEY a){
  i = AHEAP_LEAF (*H, i);
  H->v[i] += a;
  AHEAP_update (H, i);
}

/* update the ancestor of node i */
void AHEAP_update (AHEAP *H, AHEAP_ID i){
  AHEAP_ID j;
  AHEAP_KEY a = H->v[i];
  while ( i>0 ){
    j = i - 1 + (i%2)*2;   /* j = the sibling of i */
    i = (i-1) / 2;
    if ( H->v[j] < a ) a = H->v[j];
    if ( a == H->v[i] ) break;
    H->v[i] = a;
  }
}

/* find the leaf with the minimum key value among the leaves having index 
 smaller/larger than i, or between i and j */
AHEAP_ID AHEAP_upper_min (AHEAP *H, AHEAP_ID i){
  AHEAP_ID fi=0, j = AHEAP_LEAF (*H, H->end - 1);
  AHEAP_KEY fm = AHEAP_KEYHUGE;
  if ( i == 0 ) return (AHEAP_findmin_head (H) );
  i = AHEAP_LEAF (*H, i-1);
  while ( i != j ){
    if ( i%2 ){ /* if i is the child with smaller index */
      if ( fm > H->v[i+1] ){
        fm = H->v[i+1];
        fi = i+1;
      }
    }
    i = (i-1)/2;
    if ( j == i ) break;  /* stop if the right pointer and the left pointer are the same */
    j = (j-1)/2;
  }
  while ( fi < H->end-1 ) fi = fi*2 + (H->v[fi*2+1]<=fm?1:2);
  return ( AHEAP_IDX(*H, fi) );
}
AHEAP_ID AHEAP_lower_min (AHEAP *H, AHEAP_ID i){
  AHEAP_ID fi=0, j = AHEAP_LEAF (*H, 0);
  AHEAP_KEY fm = AHEAP_KEYHUGE;
  if ( i == H->end-1 ) return (AHEAP_findmin_head (H) );
  i = AHEAP_LEAF (*H, i+1);
  while ( i != j ){
    if ( i%2 == 0 ){ /* if i is the child of larger index */
      if ( fm > H->v[i-1] ){
        fm = H->v[i-1];
        fi = i-1;
      }
    }
    j = (j-1)/2;
    if ( j == i ) break;  /* stop if the right pointer and the left pointer are the same */
    i = (i-1)/2;
  }
  while ( fi < H->end-1 ) fi = fi*2 + (H->v[fi*2+1]<=fm?1:2);
  return (AHEAP_IDX(*H, fi) );
}

/* find the index having the minimum among given two indices */
AHEAP_ID AHEAP_interval_min (AHEAP *H, AHEAP_ID i, AHEAP_ID j){
  AHEAP_ID fi=0;
  AHEAP_KEY fm = AHEAP_KEYHUGE;
  if ( i == 0 ) return (AHEAP_lower_min (H, j) );
  if ( j == H->end-1 ) return (AHEAP_upper_min (H, i) );
  i = AHEAP_LEAF (*H, i-1);
  j = AHEAP_LEAF (*H, j+1);
  while ( i != j && i != j-1 ){
    if ( i%2 ){ /* if i is the child of smaller index */
      if ( fm > H->v[i+1] ){
        fm = H->v[i+1];
        fi = i+1;
      }
    }
    i = (i-1)/2;
    if ( j == i || j == i+1 ) break;  /* stop if the right pointer and the left pointer are the same */
    if ( j%2 == 0 ){ /* if j is the child of larger index */
      if ( fm > H->v[j-1] ){
        fm = H->v[j-1];
        fi = j-1;
      }
    }
    j = (j-1)/2;
  }
  while ( fi < H->end-1 )
      fi = fi*2 + (H->v[fi*2+1] <= fm?1:2);
  return (AHEAP_IDX(*H, fi) );
}

/* print heap keys according to the structure of the heap */
void AHEAP_print (AHEAP *H){
  AHEAP_ID i, j=1;
  while ( j<=H->end*2-1 ){
    FLOOP (i, j-1, MIN(j, H->end)*2-1) printf (AHEAP_KEYF ",", H->v[i] );
    printf ("\n");
    j = j*2;
  }
}

#endif
