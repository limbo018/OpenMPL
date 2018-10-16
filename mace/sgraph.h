/*  graph library by array list
            12/Feb/2002    by Takeaki Uno
    homepage:   http://research.nii.ac.jp/~uno/index.html  */
/* This program is available for only academic use, basically.
   Anyone can modify this program, but he/she has to write down 
    the change of the modification on the top of the source code.
   Neither contact nor appointment to Takeaki Uno is needed.
   If one wants to re-distribute this code, please
    refer the newest code, and show the link to homepage of 
    Takeaki Uno, to notify the news about the codes for the users. */

/****************************************************************************/
/* d = degree of node i := G->edge/in/out.v[i].t
   d = max/min (in/out) degree := VEC_MAXT(d,G->edge.v/in.v/out.v,0,...->t)  (VEC_MINT, resp.)
   #nodes :=  SGRAPH_NODE_NUM(G)
   #edges :=  G->edge.eles/2
   #arcs :=   G->in.eles or G->out.eles
   load_node_weight := ARY_LOAD_WEIGHT(G->node_w,WEIGHT,filename,counter,"errormes", EXIT) 
   load_node_weight := ARY_LOAD_WEIGHT(G->node_w,WEIGHT,filename,counter,"errormes", EXIT) 

   sort_node by size := SGRAPH_sort_node_iter (G, qsort_perm_VECt ((VEC *)Q, G->node_end, flag)
   sort_node by weight := SGRAPH_sort_node_iter (G, qsort_perm_WEIGHT (w, G->node_end, flag)
*/
/****************************************************************************/

#ifndef _sgraph_h_
#define _sgraph_h_

#include"stdlib2.h"
#include"vec.h"


/*  structure for graph  */
typedef struct {
  unsigned char type;   // structure type flag
  char *fname;      // input file name
  int flag;         // flag for load routine

  SETFAMILY edge, in, out;  // setfamily for edge, in-arc, out-arc
  QUEUE_INT node1_num;   // the size of vertex set 1, bipartite case. otherwise 0
  WEIGHT *node_w, *wbuf;    // pointer to the node weight array(int)
  PERM *perm;       // node permutation (nodes->original)
  char *wfname;     // weight file name
} SGRAPH;
extern SGRAPH INIT_SGRAPH;

#define SGRAPH_NODE_NUM(G)    MAX((G).edge.t,(G).in.t)

/*************** initialization/termination ***********************/

/*  initialization, termination, allocate arrays for weights, copy and duplication */
void SGRAPH_alloc (SGRAPH *G, int node_num, size_t edge_num, size_t arc_num);
void SGRAPH_cpy (SGRAPH *G2, SGRAPH *G);
void SGRAPH_end (SGRAPH *G);


/**************** addition/deletion **********************************/

/*  make/take/remove edge e as connecting vertices u and v,
 and  edge (u,v). 
  do nothing if when make already existing edge, or delete non-existing edge.
  with range check of parameters */

void SGRAPH_edge_mk (SGRAPH *G, QUEUE_INT u, QUEUE_INT v, WEIGHT w);
void SGRAPH_edge_rm (SGRAPH *G, QUEUE_INT u, QUEUE_INT v);
void SGRAPH_arc_mk (SGRAPH *G, QUEUE_INT u, QUEUE_INT v, WEIGHT w);
void SGRAPH_arc_rm (SGRAPH *G, QUEUE_INT u, QUEUE_INT v);


/******************* node/edge sort, and duplication reduction *********************/

/* subroutine of sort_edge_list */
void SGRAPH_sort_edge_list_iter (QUEUE *Q, WEIGHT **w, PERM *invperm, VEC_ID i, int flag);

/* sort each array list, increasing if flag=1, and decreasing if flag=-1 */
void SGRAPH_sort_edge_list (SGRAPH *G, int flag);

/* replace node i by perm and invperm */
void SGRAPH_replace_index (SGRAPH *G, PERM *perm, PERM *invperm);

/* sort the nodes by permutation given by tmp */
PERM *SGRAPH_sort_node_iter (SGRAPH *G, PERM *tmp);

/* sort the nodes by degrees, increasing if flag=1, decreasing if flag=-1 */
PERM *SGRAPH_sort_node_t (SGRAPH *G, QUEUE *Q, int flag);

/* sort the nodes by node_weight, increasing if flag=1, decreasing if flag=-1 */
PERM *SGRAPH_sort_node_w (SGRAPH *G, WEIGHT *w, int flag);

/* remove multiple edges/arcs and self loops 
   it works only when after applying sort_incident_edges */
void SGRAPH_simple (SGRAPH *G, int flag);




/******************* print routines *************************************/


/*  print graph by numbers  */
void SGRAPH_print (SGRAPH *G);

/* Write the graph to file. Edges, arcs, and nodes from 0 to node_num/edge_num/arc_num are written to file. Parameters are
  (graph) (file name) (not write edge weight => 0) (not write node weight => 0) */
void SGRAPH_save (SGRAPH *G, char *fname);

/* graph load routine. Allocate memory as much as the size of input file.
   parameters are, 
   (graph) (file name) (read edges?) (read arcs?) (read node weights?) (read edge weight?) (bipartite?) */
/* In the row of each vertex, write only vertices larger than it connected by an edge */
void SGRAPH_load (SGRAPH *G);
void SGRAPH_load_node_weight (SGRAPH *G, char *filename);

void SGRAPH_rm_selfloop (SGRAPH *G);
/*
  format of file:(including notifications to make input file)
   
  the ith row corresponds to node i-1, and
    ID list of nodes adjacent to i, and having ID > i, for undirected graph
    ID list of nodes adjacent to i by out-going arc of i, for directed graph
   Separator is ",", but graph load routine accepts any letter for 
    separator but not a number.
   If the graph has both edges and arcs, write them in two lines separately,
    so a node then uses two lines, and #nodes = #lines/2.
  
    ==  Notifications to make input file ==
   Notice that if 0th line has node 2, and the 2nd line has 0, then there
    will be multiple edge (0,2) and (2,0).
   The read routine does not make error with multiple edges, it is allowed.

   The ID of nodes begin from 0. After reading graph, node_num is set to
    node_end.

   Input file example, without weights, E={(0,1),(0,2),(1,1),(1,3),(2,3)}
===========
   1,2
   1 3
   3
   
   [EOF]
=========
   Nodes are 0,1, and 2, both edges and arcs exist, with node/edge/arc weights)
   5000,1,30
   0,50,1,20,
   100,1,3
   2,20
   200
   
   [EOF]
=======
   where node weights are 5000, 100, 200, and edges and their weights are
    (0,1),30,   (1,1),3
    arcs and their weights are (0,0),50,   (0,1), 20,   (1,2), 20

    In the case of bipartite graph, write the adjacent-node lists only for 
     the node in node set one.
     
    
*/
  

/*************************************************************************/
   

#endif 








