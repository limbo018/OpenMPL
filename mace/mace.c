/* MACE: MAximal Clique Enumerater */
/* ver 1.0 1/Sep/2005 Takeaki Uno   e-mail:uno@nii.jp, 
    homepage:   http://research.nii.ac.jp/~uno/index.html  */
/* This program is available for only academic use, basically.
   Anyone can modify this program, but he/she has to write down 
    the change of the modification on the top of the source code.
   Neither contact nor appointment to Takeaki Uno is needed.
   If one wants to re-distribute this code, do not forget to 
    refer the newest code, and show the link to homepage of 
    Takeaki Uno, to notify the news about this code for the users.
   For the commercial use, please make a contact to Takeaki Uno. */


#ifndef _mace_c_
#define _mace_c_

#define WEIGHT_DOUBLE

#include"sgraph.c"
#include"problem.c"

#define VBMMARK_MAX 16  /* MAXsize of BITMAP */ 
#define VBMINT unsigned long     /* variable type for BITMAP */
#define VBMINT_MAX 30  /* MAXsize of BITMAP */

typedef struct {
  VBMINT *edge;  /* BITMAP representation w.r.t. columns of vertices in the current clique */
  VBMINT *set, *reset;  /* array for BITMASKs */
  int *pos;   /* positions of vertices of the clique in the bitmap */
  QUEUE dellist;
  char *mark;
  int mark_max;
} MACEVBM;

//int **PP->shift;    /* pointers to the positions of the current processing items in each transaction */


void MACE_error (){
  ERROR_MES = "command explanation";
  print_err ("mace MCqVe [options] input-filename [output-filename]\n\
%%:show progress, _:no message, +:write solutions in append mode\n\
C:enumerate cliques, M:enumerate maximal cliques, e:edge_list format\n\
[options]\n\
-l [num]:output cliques with size at least [num]\n\
-u [num]:output cliques with size at most [num]\n\
-S [num]:stop after outputting [num] solutions\n\
-, [char]:give the separator of the numbers in the output\n\
-Q [filename]:replace the output numbers according to the permutation table given by [filename]\n\
if the 1st letter of input-filename is '-', be considered as 'parameter list'\n");
  EXIT;
}


/***************************************************/
/*  read parameters from command line              */
/***************************************************/
void MACE_read_param (int argc, char *argv[], PROBLEM *PP){
  ITEMSET *II=&PP->II;
  int c=1;
  if ( argc < c+2 ){ MACE_error (); return; }

  if ( !strchr (argv[c], '_') ){ II->flag |= SHOW_MESSAGE; PP->SG.flag |= SHOW_MESSAGE; }
  if ( strchr (argv[c], '%') ) II->flag |= SHOW_PROGRESS;
  if ( strchr (argv[c], '+') ) II->flag |= ITEMSET_APPEND;
  if ( strchr (argv[c], 'M') ) PP->problem = PROBLEM_MAXIMAL;
  else if ( strchr (argv[c], 'C') ) PP->problem = PROBLEM_FREQSET;
  else error ("M or C command has to be specified", EXIT);
  if ( strchr (argv[c], 'e') ) PP->SG.flag |= LOAD_ELE;
  c++;

  while ( argv[c][0] == '-' ){
    switch (argv[c][1]){
      case 'l': II->lb = atoi (argv[c+1]);
      break; default: goto NEXT;
    }
    c += 2;
    if ( argc < c+1 ){ MACE_error (); return; }
  }

  NEXT:;
  PP->SG.fname = argv[c];
  if ( argc>c+1 ) PP->output_fname = argv[c+1];
}



/***************************************************/
/*  initialization                                 */
/***************************************************/
void MACE_init (PROBLEM *PP, MACEVBM *VV){
  SGRAPH *G = &PP->SG;
  QUEUE_INT i;
  VBMINT p;

  PP->II.flag |= ITEMSET_ADD;
  G->flag |= LOAD_INCSORT + LOAD_RM_DUP + LOAD_EDGE;
  PROBLEM_load (PP);  if (ERROR_MES) return;
  PROBLEM_alloc (PP, G->edge.t, G->edge.t, G->edge.eles, NULL, PROBLEM_ITEMJUMP + PROBLEM_ITEMCAND + PROBLEM_SHIFT + PROBLEM_OCC_T);
  SGRAPH_rm_selfloop (G);
  FLOOP (i, 0, G->edge.t) G->edge.v[i].v[G->edge.v[i].t] = G->edge.t;

// delivery
  QUEUE_delivery (NULL, PP->occ_t, NULL, G->edge.v, NULL, G->edge.t, G->edge.t);
  MQUE_ALLOC (PP->OQ, G->edge.t, PP->occ_t, 0, 2, EXIT);

  if ( PP->problem & PROBLEM_CLOSED ){
    VV->edge = VV->set = VV->reset = NULL; VV->pos = NULL; VV->dellist.v = NULL;
    malloc2 (VV->edge, G->edge.t, goto ERR);
    malloc2 (VV->pos, G->edge.t, goto ERR);
    malloc2 (VV->set, VBMINT_MAX, goto ERR);
    malloc2 (VV->reset, VBMINT_MAX, goto ERR);
    QUEUE_alloc (&VV->dellist, VBMINT_MAX+2);
if ( ERROR_MES ) goto ERR;
    VV->dellist.t = VBMINT_MAX;
    ARY_FILL (VV->edge, 0, G->edge.t, 0);
    for (i=0,p=1 ; i<VBMINT_MAX ; i++,p*=2){
      VV->set[i] = p;
//        VV->reset[i] = 0xffffffff-p;
      VV->reset[i] = -1-p;
      VV->dellist.v[i] = i;
    }
//      for (i=1,MACEVBM_mark_max=1 ; i<VBMMARK_MAX ; i++ ) MACEVBM_mark_max*=2;
//      malloc2 (MACEVBM_mark, char, MACEVBM_mark_max, "MACE_init:MACEVBM_mark");
//      for (i=0 ; i<MACEVBM_mark_max ; i++) MACEVBM_mark[i] = 0;
  }
  return;
  ERR:;
  QUEUE_end (&VV->dellist);
  free2 (VV->edge);
  free2 (VV->pos);
  free2 (VV->set);
  free2 (VV->reset);
}


/******************************************************************/
/* iteration of clique enumeration   */
/******************************************************************/
void MACEclq_iter (PROBLEM *PP, QUEUE_INT v, QUEUE *occ){
  SGRAPH *G=&PP->SG;
  ITEMSET *II=&PP->II;
  QUEUE_INT *x, *xx, *y;
  QUEUE *Q = PP->OQ;

  QUE_INS (II->itemset, v);
  ITEMSET_output_itemset (II, NULL, 0);
  if ( II->itemset.t >= II->ub ) goto END;  // upper bound of clique
  
  MLOOP (y, occ->v, v){
    xx = G->edge.v[*y].v;
    MLOOP (x, occ->v, G->edge.t){
      while ( *x > *xx ) xx++;
      if ( *x == *xx ) QUE_INS (Q[*y], *x);
    }
    QUE_INS (Q[*y], G->edge.t);
    MACEclq_iter (PP, *y, &Q[*y]);
    Q[*y].t = 0;
  }
  END:;
  II->itemset.t--;
}


/******************************************************************/
/******************************************************************/
/******************************************************************/

/******************************************************************/
/* bitmap routines   */
/******************************************************************/
void MACEVBM_set_vertex (SGRAPH *G, QUEUE_INT v, MACEVBM *VV){
  QUEUE_INT *x;
  VBMINT p;
  VV->pos[v] =QUEUE_ext_tail_ (&VV->dellist);
  p = VV->set[VV->pos[v]];
  MLOOP (x, G->edge.v[v].v, G->edge.t) VV->edge[*x] |= p;
}
void MACEVBM_reset_vertex (SGRAPH *G, QUEUE_INT v, MACEVBM *VV){
  QUEUE_INT *x;
  VBMINT p;
  QUE_INS (VV->dellist, VV->pos[v]);
  p = VV->reset[VV->pos[v]];
  MLOOP (x, G->edge.v[v].v, G->edge.t) VV->edge[*x] &= p;
}
void MACEVBM_set_diff_vertexes (SGRAPH *G, QUEUE *K1, QUEUE *K2, MACEVBM *VV ){
  QUEUE_INT *x, *y = K2->v;
  MQUE_FLOOP (*K1, x){
    if ( *x == *y ) y++;
    else MACEVBM_set_vertex (G, *x, VV);
  }
}
void MACEVBM_reset_diff_vertexes (SGRAPH *G, QUEUE *K1, QUEUE *K2, MACEVBM *VV){
  QUEUE_INT *x, *y = K2->v;
  MQUE_FLOOP (*K1, x){
    if ( *x == *y ) y = y-K2->v<K2->t-1? y+1: y; 
    else MACEVBM_reset_vertex (G, *x, VV);
  }
}


/******************************************************************/
/* add a vertex v to clique K */
/******************************************************************/
void MACE_add_vertex (SGRAPH *G, PROBLEM *PP, QUEUE *K, QUEUE_INT v, MACEVBM *VV){
  QUE_INS (*K, v);
  if ( PP->problem & PROBLEM_CLOSED ){
    if ( K->t > VBMINT_MAX ) PP->problem -= PROBLEM_CLOSED;
    else MACEVBM_set_vertex (G, v, VV);
  }
}

/******************************************************************/
/* add a vertex v to clique K */
/******************************************************************/
void MACE_scan_vertex_list (SGRAPH *G, PROBLEM *PP, QUEUE_INT v, QUEUE_INT w){
  QUEUE_INT *xx;
  MQUE_MLOOP (G->edge.v[v], xx, w){
    if ( PP->OQ[*xx].t == 0 ) QUE_INS (PP->itemcand, *xx);
    QUE_INS (PP->OQ[*xx], v); 
  }
}

/* K := lex. muximum maximal clique including K (w.r.t. vertices <w ) */
/* MACE_occ[v] := N(v) \cap K */
void MACE_extend (SGRAPH *G, PROBLEM *PP, QUEUE *K, QUEUE_INT w, MACEVBM *VV){
  QUEUE_INT *x, v;

  MQUE_FLOOP (*K, x) MACE_scan_vertex_list (G, PP, *x, w);
  v = K->v[0];
  MQUE_MLOOP (G->edge.v[v], x, w);
         // x := position of vertex w in the list Q[v(= head of K)] 
  for (x-- ; x>=G->edge.v[v].v ; x--){
    if ( PP->OQ[*x].t == K->t ){
       MACE_scan_vertex_list (G, PP, *x, *x);
       MACE_add_vertex (G, PP, K, *x, VV);
    }
  }
}

/****************************************************************/
/* check the maximality of K\cap N(w) (=MACE_occ[w]),
   and whether the parent of C(K\cap N(w)) = K or not. */
/****************************************************************/
LONG MACE_parent_check (SGRAPH *G, PROBLEM *PP, QUEUE *K, QUEUE *ad, QUEUE *Q, QUEUE_INT w){
  QUEUE_INT j=0, e, i, flag =1;
  QUEUE_INT v=Q[w].v[0], *y = G->edge.v[w].v + G->edge.v[w].t-1, *x, *zz=Q[w].v, *Z;
  ad->t = 0;
  K->v[K->t] = -1;  // loop stopper

  FLOOP (i, 0, Q[w].t){
    e = Q[w].v[i];
     // pointers to the positions of the current processing items in each transaction
    PP->shift[i] = G->edge.v[e].v+G->edge.v[e].t-1;
  }
  for (x=G->edge.v[v].v + G->edge.v[v].t-1 ; *x>w ; x--){
    if ( *x <= (e=K->v[j]) ){ // skip if *x \in K (or w<*x)
      if ( *zz == e ) zz++; 
      else {  // insert *x to Qad if *x is not in K\cap N(w)
        PP->shift[Q[w].t + ad->t] = G->edge.v[e].v + G->edge.v[e].t-1-j;
        QUE_INS (*ad, e);
      }
      if ( *x < e ) x++;
      j++;
      continue;
    }
    i = 0;
    while (1){   // check *x is adjacent to all vertices in K\cap N(w) or not, one-by-one. if not, then break the loop
      while ( *PP->shift[i]>*x ) PP->shift[i]--;
      if ( *PP->shift[i] < *x ) goto LOOP_END;
      i++;
      if ( i== Q[w].t ){ // if *x is adjacent to all 
        if ( y<G->edge.v[w].v ) goto NEXT;
        while ( *y>*x ){
          if ( --y < G->edge.v[w].v ) goto NEXT;
        }
        if ( *y==*x ) return (*x);  //if *x is adjacent to w, then not maximal
        break;
      }
    }
    NEXT:;
    while (flag){   // check *x is adjacent to all vertices in K_{\le w}. If not, then break the loop
      if ( i== Q[w].t + ad->t ) return (*x); // if *x is adjacent to all, MACE_occ[w] is not a child
      Z = G->edge.v[ad->v[i-Q[w].t]].v;
      while ( *PP->shift[i]>*x ){
        PP->shift[i]--;
        if ( PP->shift[i] < Z ){ flag = 0; goto LOOP_END; } // reached to the end of the adjacency list of the i-th added vertex, thus no further vertex can be pass this check, and set flag to 0 not to come here again.
      }
      if ( *PP->shift[i] < *x ) goto LOOP_END;
      i++;
    }
    LOOP_END:;
  }
  return (-1);
}

/****************************************************************/
/* check the maximality of K\cap N(w) (=MACE_occ[w]),
   ad whether the parent of C(K\cap N(w)) = K or not. */
/*  BITMAP version */
/****************************************************************/
LONG MACEVBM_parent_check (SGRAPH *G, QUEUE *K, QUEUE *Q, QUEUE_INT w, MACEVBM *VV){
  QUEUE_INT v=Q[w].v[0];
  QUEUE_ID i;
  VBMINT p=0, pp;
  QUEUE_INT *y = G->edge.v[w].v + G->edge.v[w].t-1, *x, *z=K->v;
  K->v[K->t] = -1;  // loop stopper
  FLOOP (i, 0, Q[w].t) p |= VV->set[VV->pos[Q[w].v[i]]];
  pp = p;
  for (x=G->edge.v[v].v + G->edge.v[v].t-1 ; *x>w ; x--){
    while ( *x < *z ){ pp |= VV->set[VV->pos[*z]]; z++; }
//    if ( *z>=0 ) pp |= VV->set[VV->pos[*z]];

    if ( *x == *z ){ pp |= VV->set[VV->pos[*z]]; z++; continue; }
    if ( pp==(pp&VV->edge[*x]) ) return (*x);  // parentness

    if ( p == (p & VV->edge[*x]) ){ // maximality w.r.t P\cap N(w) (=occ[w]) 
      if ( y<G->edge.v[w].v ) goto NEXT;
      while ( *x < *y ){  // check *x is incident to w?  
        y--;
        if ( y<G->edge.v[w].v ) goto NEXT;
      }
      if ( *x==*y ) return (*x); // if *x is incident to w, parent is different.
    }
    NEXT:;
  }
  return (-1);
}

/*************************************************************************/
/*  simple routine for checking the maximality of a clique,
     and parent-child relation, for debugging */
/*************************************************************************/
LONG MACE_parent_check_max (SGRAPH *G, QUEUE *K, QUEUE *ad, QUEUE *Q, QUEUE_INT w){
  QUEUE_INT *x;
  QUEUE_cpy (ad, &G->edge.v[w]);
  MQUE_FLOOP (Q[w], x) QUEUE_and_ (ad, &G->edge.v[*x]);
  if ( ad->t==0 ) return (-1);
  if ( ad->v[ad->t-1] > w ) return (ad->v[ad->t-1]);
  return (-1);
}
LONG MACE_parent_check_parent (SGRAPH *G, QUEUE *K, QUEUE *ad, QUEUE *Q, QUEUE_INT w){
  QUEUE_INT t=0, *x, i;
  K->v[K->t] = -1; // loop stopper;
  QUEUE_cpy (ad, &G->edge.v[Q[w].v[0]]);
  MQUE_FLOOP (Q[w], x) QUEUE_and_ (ad, &G->edge.v[*x]);
  while ( ad->t > 0 ){
    i = QUEUE_ext_tail_ (ad);
    if ( i<w ) return (-1);
    while ( i<K->v[t] ) t++;
    if ( i > K->v[t] ) return (i);
    QUEUE_and_ (ad, &G->edge.v[i]);
  }
  return (-1);
}


/***************************************************************/
/*  under construction for future improvements */
/***************************************************************/
void MACE_make_same_list (SGRAPH *G, PROBLEM *PP, QUEUE_INT v, MACEVBM *VV){
  QUEUE *Q = PP->OQ;
  VBMINT p;
  QUEUE_ID i;
  QUEUE_INT u, *x;
  if ( Q[v].t > VBMMARK_MAX ) return;
  ARY_FILL (VV->mark, 0, VV->mark_max, 0);

  FLOOP (i, 0, Q[v].t) G->edge.v[Q[v].v[i]].s = i;
  FLOOP (i, PP->itemcand.s, PP->itemcand.t){
    u = PP->itemcand.v[i];
    p = 0;
    MQUE_FLOOP (Q[u], x) p |= VV->set[G->edge.v[*x].s];
//    if (MACEVBM_mark[p]) c1++;
//    else {MACEVBM_mark[p]++; c2++;}
  }
  MQUE_FLOOP (Q[v], x) G->edge.v[*x].s = 0;
}


/******************************************************************/
/******************************************************************/
/******************************************************************/


/*************************************************************************/
/* MACE main iteration */
/*************************************************************************/
void MACE_iter (PROBLEM *PP, int v, MACEVBM *VV){
  LONG ii;
  QUEUE_INT u;
  QUEUE_ID js = PP->itemcand.s;
  QUEUE *Q = PP->OQ;
  SGRAPH *G = &PP->SG;
  ITEMSET *II = &PP->II;

//printf ("%d:   ", II->iters);
//printf ("%d:  ", v); QUEUE_print__ (&Q[v]);

  II->iters++;
  PP->itemcand.s = PP->itemcand.t;
  MACE_add_vertex (G, PP, &Q[v], v, VV);
  MACE_extend (G, PP, &Q[v], v, VV);
//  MACE_make_same_list (v);

  II->itemset.t = 0;
  memcpy (II->itemset.v, Q[v].v, sizeof(QUEUE_INT)*Q[v].t);
  II->itemset.t = Q[v].t;
  ITEMSET_output_itemset (II, NULL, 0);
  qsort_QUEUE_INT (&PP->itemcand.v[PP->itemcand.s], PP->itemcand.t-PP->itemcand.s, -1);

  while ( PP->itemcand.t > PP->itemcand.s ){
    u = QUEUE_ext_tail_(&PP->itemcand);
    if ( u == QUEUE_TAIL_(Q[v]) ){
      Q[v].t--;
      if ( PP->problem & PROBLEM_CLOSED ) MACEVBM_reset_vertex (G, u, VV);
    } else {
      if ( PP->problem & PROBLEM_CLOSED ) ii = MACEVBM_parent_check (G, &Q[v], PP->OQ, u, VV);
      else ii = MACE_parent_check (G, PP, &Q[v], &II->add, PP->OQ, u);
       if ( ii==-1 ){
        if (PP->problem & PROBLEM_CLOSED) MACEVBM_reset_diff_vertexes (G, &Q[v], &Q[u], VV);
        MACE_iter (PP, u, VV); // recursive call for a child
        if (PP->problem & PROBLEM_CLOSED) MACEVBM_set_diff_vertexes (G, &Q[v], &Q[u], VV);
      }
    }
    Q[u].t = 0;
  }
  PP->itemcand.s = js;
  if ( PP->problem & PROBLEM_CLOSED ) MACEVBM_reset_vertex (G, v, VV);
}

/* MACE main */
void MACE (PROBLEM *PP, MACEVBM *VV){
  QUEUE *E = PP->SG.edge.v;
  ITEMSET *II=&PP->II;
  QUEUE_INT v;

  FLOOP (v, 0, PP->SG.edge.t){
    if ( E[v].t==0 ){
      II->itemset.t = 0;
      QUE_INS (II->itemset, v);
      ITEMSET_output_itemset (II, NULL, 0);
    } else if ( E[v].v[E[v].t-1] <= v ){
      MACE_iter (PP, v, VV);
    }
    PP->OQ[v].t = 0;
  }
}


// main routine
int MACE_main (int argc, char *argv[]){
  QUEUE_INT v, flag = 0;
  PROBLEM PP;
  SGRAPH *G=&PP.SG;
  MACEVBM VV;

  PROBLEM_init (&PP);
  MACE_read_param (argc, argv, &PP);
if ( ERROR_MES ) return (1);
//  if ( weightflag ) MACE_perm = SGRAPH_sort_node_w (&G, G->node_w, 1);
  if ( PP.problem & PROBLEM_MAXIMAL ){
    PP.problem |= PROBLEM_CLOSED;  // flag for the use of VBM
    flag = 1;
  }
  MACE_init (&PP, &VV);
//  if ( weightflag ) ENUWMCLQ ();
  if ( !ERROR_MES && G->edge.eles > 0 ){
    if ( PP.problem & PROBLEM_FREQSET ){
      FLOOP (v, 0, G->edge.t) MACEclq_iter (&PP, v, &G->edge.v[v]);
    } else MACE (&PP, &VV);
    ITEMSET_last_output (&PP.II);
  }
  
  PROBLEM_end (&PP);
  if ( flag ){
    free2 (VV.edge);
    free2 (VV.pos);
    free2 (VV.set);
    free2 (VV.reset);
    QUEUE_end (&VV.dellist);
  }

  return (ERROR_MES?1:0);
}

/*******************************************************************************/
#ifndef _NO_MAIN_
#define _NO_MAIN_
int main (int argc, char *argv[]){
  return (MACE_main (argc, argv) );
}
#endif
/*******************************************************************************/


#endif

