/* library for sparse vector */
/* Takeaki Uno    27/Dec/2008 */

#ifndef _vec_h_
#define _vec_h_

//#define USE_MATH

#include"math.h"
#include"queue.h"

#ifndef SVEC_VAL
 #ifdef SVEC_VAL_INT
  #define SVEC_VAL int
  #define SVEC_VAL2 LONG
  #define SVEC_VAL_END INTHUGE
  #define SVEC_VAL2_END LONGHUGE
  #define SVEC_VALF "%d"
 #else
  #define SVEC_VAL double
  #define SVEC_VAL2 double
  #define SVEC_VAL_END DOUBLEHUGE
  #define SVEC_VAL2_END DOUBLEHUGE
  #define SVEC_VALF "%f"
 #endif
#endif

#define VEC_LOAD_BIN 16777216   // read binary file 
#define VEC_LOAD_BIN2 33554432   // read binary file with 2byte for each number 
#define VEC_LOAD_BIN4 67108864   // read binary file with 4byte for each number
#define VEC_LOAD_CENTERIZE 134217728  // read binary file, and minus the half(128) from each number
#define VEC_NORMALIZE 268435456  // read binary file, and minus the half(128) from each number

/* matrix */
typedef struct {
  unsigned char type;  // mark to identify type of the structure
  char *fname;      // input file name
  int flag;         // flag

  VEC *v;
  VEC_ID end;
  VEC_ID t;
  VEC_VAL *buf, *buf2;
  VEC_ID clms;
  size_t eles;
  VEC_VAL *buf_org, *buf2_org;
} MAT;

/* sparse vector, element */
typedef struct {
  QUEUE_ID i;
  SVEC_VAL a;
} SVEC_ELE;

/* sparse vector, vector */
typedef struct {
  unsigned char type;  // mark to identify type of the structure
  SVEC_ELE *v;
  VEC_ID end;
  VEC_ID t;
} SVEC;

/* sparse vector, matrix */
typedef struct {
  unsigned char type;  // mark to identify type of the structure
  char *fname;      // input file name
  int flag;         // flag

  SVEC *v;
  VEC_ID end;
  VEC_ID t;
  SVEC_ELE *buf, *buf2;
  VEC_ID clms;
  size_t eles, ele_end;
} SMAT;

/* set family */
typedef struct {
  unsigned char type;  // mark to identify type of the structure
  char *fname;      // input file name
  int flag;         // flag

  QUEUE *v;
  VEC_ID end;
  VEC_ID t;
  QUEUE_INT *buf, *buf2;
  VEC_ID clms;
  size_t eles, ele_end;
  WEIGHT *cw, *rw, **w, *wbuf;
  int unit;
  char *wfname, *cwfname, *rwfname;     // weight file name
} SETFAMILY;

#define INIT_SETFAMILY_ {TYPE_SETFAMILY,NULL,0,NULL,0,0,NULL,NULL,0,0,0,NULL,NULL,NULL,NULL,sizeof(QUEUE_INT),NULL,NULL,NULL}

extern MAT INIT_MAT;
extern SVEC INIT_SVEC;
extern SMAT INIT_SMAT;
extern SETFAMILY INIT_SETFAMILY;

QSORT_TYPE_HEADER (SVEC_VAL, SVEC_VAL)
QSORT_TYPE_HEADER (SVEC_VAL2, SVEC_VAL2)

#define   ARY_QUEUE_INPRO(f,U,V)  do{(f)=0;FLOOP(common_QUEUE_ID, 0, (QUEUE_ID)(U).t)(f)+=(V)[(U).v[common_QUEUE_ID]];}while(0)
#define   ARY_SVEC_INPRO(f,U,V)  do{(f)=0;FLOOP(common_VEC_ID, 0, (VEC_ID)(U).t)(f)+=((double)(U).v[common_VEC_ID].a)*(V)[(U).v[common_VEC_ID].i];}while(0)

/* terminate routine for VEC */
void VEC_end (VEC *V);
void MAT_end (MAT *M);
void SVEC_end (SVEC *V);
void SMAT_end (SMAT *M);
void SETFAMILY_end (SETFAMILY *M);

/* allocate memory according to rows and rowt */
void VEC_alloc (VEC *V, VEC_ID clms);
void MAT_alloc (MAT *M, VEC_ID rows, VEC_ID clms);
void SVEC_alloc (SVEC *V, VEC_ID end);
void SMAT_alloc (SMAT *M, VEC_ID rows, VEC_ID *rowt, VEC_ID clms, size_t eles);
void SETFAMILY_alloc (SETFAMILY *M, VEC_ID rows, VEC_ID *rowt, VEC_ID clms, size_t eles);
void SETFAMILY_alloc_weight (SETFAMILY *M);

/* count/read the number in file for MAT */
/* if *rows>0, only read count the numbers in a row, for the first scan. */
void MAT_load_bin (MAT *M, FILE2 *fp, int unit);
void MAT_file_load (MAT *M, FILE2 *fp);
void MAT_load (MAT *M);
void SMAT_load (SMAT *M);
void SETFAMILY_load (SETFAMILY *M);
void SETFAMILY_load_weight (SETFAMILY *M);
void SETFAMILY_load_row_weight (SETFAMILY *M);
void SETFAMILY_load_column_weight (SETFAMILY *M);

void MAT_print (FILE *fp, MAT *M);
void SVEC_print (FILE *fp, SVEC *M);
void SMAT_print (FILE *fp, SMAT *M);
void SETFAMILY_print (FILE *fp, SETFAMILY *M);
void SETFAMILY_print_weight (FILE *fp, SETFAMILY *M);
  

/* norm, normalization **************************/
double SVEC_norm (SVEC *V);
void SVEC_normalize (SVEC *V);

/* inner product **************************/
SVEC_VAL2 SVEC_inpro (SVEC *V1, SVEC *V2);

/** Euclidean distance routines *********************************/
double VEC_eucdist (VEC *V1, VEC *V2);
double SVEC_eucdist (SVEC *V1, SVEC *V2);
double VEC_SVEC_eucdist (VEC *V1, SVEC *V2);
double QUEUE_eucdist (QUEUE *Q1, QUEUE *Q2);
double VEC_QUEUE_eucdist (VEC *V, QUEUE *Q);

void VEC_rand_gaussian (VEC *V);

double VEC_linfdist (VEC *V1, VEC *V2);

/* compute the inner product, Euclidean distance for multi vector */
double MVEC_norm (void *V);
double MVEC_inpro (void *V, void *U);
double MVEC_double_inpro (void *V, double *p);
double MVEC_eucdist (void *V, void *U);

/* compute the inner product, euclidean distance for i,jth vector */
double MMAT_inpro_ij (void *M, int i, int j);
double MMAT_double_inpro_i (void *M, int i, double *p);
double MMAT_eucdist_ij (void *M, int i, int j);
double MMAT_norm_i (void *M, int i);


#endif
