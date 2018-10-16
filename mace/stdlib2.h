/* library for standard macros and functions 
 by Takeaki Uno 2/22/2002   e-mail: uno@nii.jp
    homepage:   http://research.nii.ac.jp/~uno/index.html  */
/* This program is available for only academic use, basically.
   Anyone can modify this program, but he/she has to write down 
    the change of the modification on the top of the source code.
   Neither contact nor appointment to Takeaki Uno is needed.
   If one wants to re-distribute this code, please 
    refer the newest code, and show the link to homepage of 
    Takeaki Uno, to notify the news about the codes for the users. */

#ifndef _stdlib2_h_
#define _stdlib2_h_

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<stdarg.h>

#ifdef USE_SIMD
#include<emmintrin.h>   // use MMX-SSE2;
#endif


#if defined(__cplusplus) && defined(__GNUC__)
 #define _cplusplus_
#endif

#ifdef MULTI_CORE
#include <sys/types.h>
#include <pthread.h>
#endif
#define CORE_MAX 16

/* comment out the following line if no error check is needed */
//#define ERROR_CHECK
/* comment out the following if exit is not needed after each error routine */
//#define ERROR_RET

#ifdef ERROR_RET   // definition of the process for errors
  #define EXIT return
  #define EXIT0 return(0)
#else 
  #define EXIT exit(1)
  #define EXIT0 exit(1)
#endif

// for dealing with files more than 2GB
#define _LARGEFILE_SOURCE
#define _FILE_OFFSET_BITS 64

#ifndef NULL
  #define NULL 0
#endif

#ifdef MTWISTER
 #define RANDOM ((long)(dsfmt_gv_genrand_close_open()*2147483648LL))
 #define RAND1 dsfmt_gv_genrand_close_open()
 #define RAND_INIT dsfmt_gv_init_gen_rand(514346237)
#elif defined(__GNUC__)
 #define RANDOM xor128()
 #define RAND1 ((double)xor128())/4294967296.0
 #define RAND_INIT xor128()
#else
 #define RANDOM rand()
 #define RAND1 ((double)rand())/2147483648.0
 #define RAND_INIT srand(0)
#endif


// 64bit integer
#ifdef LONG_32
 #define LONG int
 #define LONGHUGE INTHUGE
 #define LONGF "%d"
#elif !defined(LONG)
 #define LONG long long
 #define LONGHUGE 9000000000000000000LL
 #define LONGF "%lld"
#endif

// actual int (most proper sized integer, for the processor)
#ifdef INT_64
 #define INT LONG
 #define INTF LONGF
#else
 #define INT int
 #define INTF "%d"
#endif

#ifndef FILE_LONG 
 #define FILE_LONG LONG
 #define FILE_LONGHUGE LONGHUGE
 #define FILE_LONGF LONGF
#endif

#define UINTHUGE 4000000000U
#define INTHUGE 2000000000
#define USHORTHUGE 32767
#define SHORTHUGE 65535
#define DOUBLEHUGE 999999999999999999999999999999.9
#define ISEQUAL_VALUE   0.0000001
#define ISEQUAL_VALUE2  0.00001
#define PI 3.1415926535897932384647950288
#define PI_INT 31416
#define NPE 2.718281828459045235360287471352
#define NPE_INT 27183

#ifndef WEIGHT
 #ifdef WEIGHT_DOUBLE
  #define WEIGHT double
  #define WEIGHTHUGE DOUBLEHUGE
  #define WEIGHTF "%f"
 #else   //  define WEIGHT by int if it's undefined
  #define WEIGHT int
  #define WEIGHTHUGE INTHUGE
  #define WEIGHTF "%d"
 #endif
#endif

#ifndef PERM
 #ifdef PERM_LONG
  #define PERM LONG
  #define PERMHUGE LONGHUGE
  #define PERMF LONGF
 #else
  #define PERM int
  #define PERMHUGE INTHUGE
  #define PERMF "%d"
 #endif
#endif


extern size_t common_size_t;
extern INT common_INT, common_INT2;
extern char  *common_pnt, *common_charp;
extern FILE *common_FILE;
extern WEIGHT common_WEIGHT, *common_WEIGHTp;
extern char *ERROR_MES;
extern int print_time_flag;
extern char common_comm[], common_comm2[], *common_argv[];
typedef struct {
  int i1, i2, i3, i4, i5, i6, i7, i8, i9;
  LONG l1, l2, l3, l4, l5, l6, l7, l8, l9;
  double d1, d2, d3, d4, d5, d6, d7, d8, d9;
  char *s1, *s2, *s3, *s4, *s5, *s6, *s7, *s8, *s9;
  void *p1, *p2, *p3, *p4, *p5, *p6, *p7, *p8, *p9;
} PARAMS;
extern PARAMS internal_params;

/* lock&unlock for multi-core mode */
#ifdef MULTI_CORE
  extern int SPIN_LOCK_dummy;
  #define SPIN_LOCK(b,a)    (SPIN_LOCK_dummy=(((b)>1)&&pthread_spin_lock(&(a))))
  #define SPIN_UNLOCK(b,a)    (SPIN_LOCK_dummy=(((b)>1)&&pthread_spin_unlock(&(a))))
#else
  #define SPIN_LOCK(b,a)
  #define SPIN_UNLOCK(b,a)
#endif

#define TYPE_VEC 1
#define TYPE_MAT 2
#define TYPE_SVEC 3
#define TYPE_SMAT 4
#define TYPE_QUEUE 5
#define TYPE_SETFAMILY 6
#define TYPE_TRSACT 7
#define TYPE_ALIST 8
#define TYPE_MALIST 9
#define TYPE_AGRAPH 10
#define TYPE_SGRAPH 11
#define TYPE_AHEAP 12
#define TYPE_BASE 13
#define TYPE_FSTAR 14
#define TYPE_SEQ 15
#define TYPE_BARRAY 16

#define TYPE_FILE2 32

#define ADDR_FLOOR16(x)  do{common_charp=((char *)x)+15;x=(typeof(x))(common_charp-(((size_t)common_charp)&15));}while(0)
double SQRT(double x);

/*  random  */
#define rnd(a)    (random()%(a))
#define prob(a)   ((random()%65536)<(int)((a)*65536.0))

#define MARK      1
#define UNMARK    0
#define TRUE      1
#define FALSE     0

/*  equal/inequal with allowing numerical error for double  */
#define ISEQUAL(a,b)	((a)-(b)<ISEQUAL_VALUE&&(b)-(a)<ISEQUAL_VALUE)
#define ISGREAT(a,b)    ((a)-(b)>ISEQUAL_VALUE)
#define ISLESS(a,b) 	((b)-(a)>ISEQUAL_VALUE)
#define RANGE(a,b,c)  (((a)<=(b))&&((b)<=(c)))
#define BITRM(a,b)    ((a)-=((a)&(b)));

/*  macro for getting maximum/minimum of two values  */
#define MAX(a,b)      (((a)>(b))?(a):(b))
#define ENMAX(a,b)      ((a)=(((a)>(b))?(a):(b)))
#define MIN(a,b)      (((a)<(b))?(a):(b))
#define ENMIN(a,b)      ((a)=(((a)<(b))?(a):(b)))

/*  error routine  */
#define error(mes,x)        do{ERROR_MES=mes;fprintf(stderr,"%s\n",mes);x;}while(0)
#define error_num(mes,n,x)  do{ERROR_MES=mes;fprintf(stderr,"%s: %g\n",mes,(double)(n));x;}while(0)
#define error_range(f,x,y,mes)  do{if(!RANGE((x),(f),(y))){ERROR_MES=mes;fprintf(stderr,"%s: %g\n",mes,(double)(f));EXIT;}}while(0)
#define error_str(mes,s,x)  do{ERROR_MES=mes;fprintf(stderr,"%s: %s\n",mes,s);x;}while(0)
#define print_err(...)      fprintf(stderr,__VA_ARGS__)
#define print_mes(S,...)  do{if(((S)->flag)&1)fprintf(stderr,__VA_ARGS__);}while(0)
//#define print_fname(s,fname,...)  do{if(fname)fprintf(stderr,s,fname,__VA_ARGS__);}while(0)
#define print_fname(s,fname,...)  do{if(fname)fprintf(stderr,s,fname);}while(0)
#define mfree(...)          mfree_(NULL, __VA_ARGS__, (void *)1)
#define mfree2(...)          mfree2_(NULL, __VA_ARGS__, (void *)1)

  
/* basic array operations and loops */
#define   ARY_FILL(f,start,end,c) do{for(common_size_t=(size_t)(start);common_size_t<(size_t)(end);common_size_t++)(f)[common_size_t]=(c);}while(0)
//#define   ARY_FLOOP(V,i,x)  for( (i)=0,x=(V).v[0] ; (i)<(V).t ; (i)++,x=(V).v[i])
//#define   ARY_BLOOP(V,i,x)  for( (i)=(V).t-1,x=i>0?(V).v[i]:0 ; (i)>=0 ; (i)--,(x)=(i)>0?(V).v[i]:0)
#define   FLOOP(i,x,y)  for ((i)=(x) ; (i)<(y) ; (i)++)
#define   BLOOP(i,x,y)  for ((i)=(x) ; ((i)--)>(y) ; )
#define   MLOOP(z,x,M)  for ((z)=(x) ; *(z)<(M) ; (z)++)
/* binary search: return maximum index no larger than c */
#define BIN_SRCH(x,a,s,t,c) \
 do {\
  x=s; common_size_t=t; while ( x<common_size_t-1 ){\
   if ( a[(x+common_size_t)/2] <= c ) x = (x+common_size_t)/2; else common_size_t = (x+common_size_t)/2;\
 } while (0)\

/* allocate memory, and exit with error message if fault */
#ifdef _cplusplus_
#define   malloc2(f,b,x)     do{if(!((f)=(typeof(f))malloc(((size_t)sizeof((f)[0]))*(b)))){fprintf(stderr,"memory allocation error: line %d" #f " (" LONGF " byte)\n",__LINE__,(LONG)((LONG)sizeof((f)[0])*(b)));ERROR_MES="out of memory";x;}}while(0)
#define   calloc2(f,b,x)     do{if(!((f)=(typeof(f))calloc(sizeof((f)[0]),b))){fprintf(stderr,"memory allocation error: line %d" #f " (" LONGF " byte)\n",__LINE__,(LONG)((LONG)sizeof((f)[0])*(b)));ERROR_MES="out of memory";x;}}while(0)
#define   realloc2(f,b,x)    do{if(!(f=(typeof(f))realloc(f,((size_t)sizeof((f)[0]))*(b)))){fprintf(stderr,"memory allocation error: line %d" #f " (" LONGF " byte)\n",__LINE__,(LONG)((LONG)sizeof((f)[0])*(b)));ERROR_MES="out of memory";x;}}while(0)
#else
#define   malloc2(f,b,x)     do{if(!((f)=malloc(((size_t)sizeof((f)[0]))*(b)))){fprintf(stderr,"memory allocation error: line %d" #f " (" LONGF " byte)\n",__LINE__,(LONG)((LONG)sizeof((f)[0])*(b)));ERROR_MES="out of memory";x;}}while(0)
#define   calloc2(f,b,x)     do{if(!((f)=calloc(sizeof((f)[0]),b))){fprintf(stderr,"memory allocation error: line %d: " #f " (" LONGF " byte)\n",__LINE__,(LONG)((LONG)sizeof((f)[0])*(b)));ERROR_MES="out of memory";x;}}while(0)
#define   realloc2(f,b,x)    do{if(!(f=realloc(f,((size_t)sizeof((f)[0]))*(b)))){fprintf(stderr,"memory allocation error: line %d: " #f " (" LONGF " byte)\n",__LINE__, (LONG)((LONG)sizeof((f)[0])*(b)));ERROR_MES="out of memory";x;}}while(0)
#endif

#define   malloc2d(f,b,d,x)  do{malloc2(f,b,x);malloc2((f)[0],(b)*(d)+2,{free2(f);x;});FLOOP(common_size_t,0,(size_t)b)(f)[common_size_t]=&((f)[0][common_size_t*(d)]);}while(0)
#define   calloc2d(f,b,d,x)  do{malloc2(f,b,x);calloc2((f)[0],((size_t)(b))*((size_t)(d))+2,{free2(f);x;});FLOOP(common_size_t,0,(size_t)b)(f)[common_size_t]=&((f)[0][common_size_t*(d)]);}while(0)

/* reallocate memory and expand the memory size */
#define   reallocx_(f,end,end2,e,x)  do{realloc2(f,end2,x);FLOOP(common_size_t,(size_t)end,(size_t)end2)(f)[common_size_t]=(e);}while(0)
#define   reallocx(f,end,i,e,x)  do{if((size_t)(i)>=(size_t)(end)){reallocx_(f,end,MAX((end)*2+16,(i)+1),e,x);end=MAX((end)*2,(i)+1);}}while(0)
#define   reallocz(f,end,i,x)  do{if((size_t)(i)>=(size_t)(end)){end=MAX((end)*2,(i)+1);realloc2(f,end,x);}}while(0)
#define   realloci(f,i,x)  do{if(!((i)&((i)-1)))realloc2(f,(i)*2+1,x);}while(0)

/* basic array operations */
#define   ARY_REALLOCX(f,i,e,x)  reallocz((f).v,(f).end,i,e,x)
#define   ARY_REALLOCZ(f,i,x)  reallocz((f).v,(f).end,i,x)
#define   ARY_DUP(f,p,end,x)  do{malloc2(f,end,x);memcpy(f,p,sizeof(*(f))*(end));}while(0)
#define   ARY_MAX(m,i,f,x,y)   do{(m)=(f)[x];(i)=(x);FLOOP(common_INT,(x)+1,(y))if((m)<(f)[common_INT]){(i)=common_INT;(m)=(f)[i];}}while(0)
#define   ARY_MIN(m,i,f,x,y)   do{(m)=(f)[x];(i)=(x);FLOOP(common_INT,(x)+1,y)if((m)>(f)[common_INT]){(i)=common_INT;(m)=(f)[i];}}while(0)
#define   ARY_SUM(f,v,x,y)       do{(f)=0;FLOOP(common_INT,x,y)(f)+=(v)[common_INT];}while(0)
#define   ARY_NORM(f,v,b)      do{(f)=0;FLOOP(common_INT,0,b)(f)+=(v)[common_INT]*(v)[common_INT];(f)=sqrt(f);}while(0)
#define   ARY_NORMALIZE(v,b)   do{ARY_NORM(common_double,v,b);FLOOP(common_INT,0,b)(v)[common_INT]/=common_double;}while(0)
#define   ARY_INPRO(f,u,v,b)   do{(f)=0;for (common_INT=0 ; common_INT<(b)-3 ; common_INT+=4) (f)+=(u)[common_INT]*(v)[common_INT] + (u)[common_INT+1]*(v)[common_INT+1] + (u)[common_INT+2]*(v)[common_INT+2] + (u)[common_INT+3]*(v)[common_INT+3]; if (common_INT+1<(b)){(f)+=(u)[common_INT]*v[common_INT]+(u)[common_INT+1]*(v)[common_INT+1]; if (common_INT+2<(b)) (f)+=(u)[common_INT+2]*(v)[common_INT+2];} else if (common_INT<(b)) (f)+=(u)[common_INT]*(v)[common_INT];}while(0)
#define   ARY_DIST(f,u,v,b)   do{(f)=0;for (common_INT=0 ; common_INT<(b) ; common_INT++) (f)+=((u)[common_INT]-(v)[common_INT])*((u)[common_INT]-(v)[common_INT]); (f)=sqrt(f);}while(0)


/* macros for permutation arrays */
#define   ARY_INIT_PERM(f,end)   do{FLOOP(common_INT,0,(INT)end)(f)[common_INT]=common_INT;}while(0)
#define   ARY_INV_PERM_(f,p,end)      do{ARY_FILL(f,0,end,-1);FLOOP(common_INT,0,end)if((p)[common_INT]>=0&&(p)[common_INT]<(end))(f)[(p)[common_INT]]=common_INT;}while(0)
#define   ARY_INV_PERM(f,p,end,x) do{malloc2(f,end,x);ARY_INV_PERM_(f,p,end);}while(0)
#define   ARY_RND_PERM_(f,end)        do{(f)[0]=0;FLOOP(common_INT,1,end){common_INT2=rnd(common_INT+1);(f)[common_INT]=(f)[common_INT2];(f)[common_INT2]=common_INT;}}while(0)
#define   ARY_RND_PERM(f,end,x)   do{malloc2(f,end,x);ARY_RND_PERM_(f,end);}while(0)
  /* permute f so that f[i]=f[p[i]] (inverse perm). p will be destroyed (filled by end) */
#define   ARY_INVPERMUTE_(f,p,s,end)  do{ FLOOP(common_INT,0,end){ if ( (p)[common_INT]<(end) ){ (s)=(f)[common_INT]; do { common_INT2=common_INT; common_INT=(p)[common_INT]; (f)[common_INT2]=(f)[common_INT]; (p)[common_INT2]=end; }while ( (p)[common_INT]<(end) ); (f)[common_INT2] = (s);}}}while(0)
  /* permute f so that f[i]=f[p[i]] (inverse perm). not destroy p by allocating tmp memory */
#define   ARY_INVPERMUTE(f,p,s,end,x) do{ calloc2(common_pnt,end,x);FLOOP(common_INT,0,end){ if ( common_pnt[common_INT]==0 ){ (s)=(f)[common_INT]; do{ common_INT2=common_INT; common_INT=(p)[common_INT]; (f)[common_INT2]=(f)[common_INT]; common_pnt[common_INT2]=1; }while( common_pnt[common_INT]==0 ); (f)[common_INT2] = (s);}} free(common_pnt); }while(0)
//#define   ARY_PERM(f,p,s,mark,end)         do{FLOOP(common_size_t,0,end){  }}while(0)

/* macros for printing (writing to file) arrays */
#define   ARY_PRINT(f,x,y,a) do{FLOOP(common_size_t,x,y)printf(a,(f)[common_size_t]);printf("\n");}while(0)
#define   ARY_FPRINT(fp,f,x,y,a) do{FLOOP(common_size_t,(size_t)x,(size_t)y)fprintf((FILE *)fp,a,(f)[common_size_t]);fputc('\n',(FILE *)fp);}while(0)

#define   ST_MAX(m,i,S,a,x,y)   do{(m)=(S)[x].a;(i)=(x);FLOOP(common_INT,(x)+1,y)if((m)<(S)[common_INT].a){(i)=common_INT;(m)=(S)[i].a;}}while(0)
#define   ST_MIN(m,i,S,a,x,y)   do{(m)=(S)[x].a;(i)=(x);FLOOP(common_INT,(x)+1,y)if((m)>(S)[common_INT].a){(i)=common_INT;(m)=(S)[i].a;}}while(0)
#define   ST_SUM(k,S,a,x,y)       do{(k)=0;FLOOP(common_INT,x,y)(k)+=(S)[common_INT].a;}while(0)
#define   ST_FILL(S,a,start,end,c) do{for(common_INT=(start);common_INT<(end);common_INT++)(S)[common_INT].a = (c);}while(0)
#define   ST_PRINT(S,a,x,y,f) do{FLOOP(common_size_t,x,y)printf(f,(S)[common_size_t].a );printf("\n");}while(0)


/* macros for QUE type structure (have .s, .t, .v) */
#define   QUE_EXP(f,a,x)   do{reallocx((f).v,a,(f).end,(f).t,e,x);}while(0)
#define   QUE_INS(f,b)   do{(f).v[(f).t++]=(b);}while(0)
#define   QUE_INSZ(f,a,x)    do{reallocz((f).v,(f).end,(f).t+1,x);ARY_INS(f,a);}while(0)
#define   QUE_INIT_PERM(f,end)   do{ARY_INIT_PERM((f).v,(end));(f).t=(end);(f).s=0;}while(0)

/* a macro for open files with exiting if an error occurs */
#ifdef _MSC_
 #define   fopen2(f,a,b,x)     do{fopen_s(&f,a,b);if(!f){ERROR_MES="file open error";fprintf(stderr,"file open error: file name %s, open mode %s\n",a,b);x;}}while(0)
#else
 #define   fopen2(f,a,b,x)     do{if(!((f)=fopen(a,b))){ERROR_MES="file open error";fprintf(stderr,"file open error: file name %s, open mode %s\n",a,b);x;}}while(0)
#endif
#define   FILE2_open(f,a,b,x) do{if(a)fopen2((f).fp,a,b,x);else(f).fp=NULL;malloc2((f).buf_org,FILE2_BUFSIZ+1,x);(f).buf=(f).buf_org;(f).buf_end=(f).buf_org-1;(f).bit=0;*(f).buf=0;}while(0)
#define   FILE2_open_(f,a,x)  do{(f).fp=a;malloc2((f).buf_org,FILE2_BUFSIZ+1,x);(f).buf=(f).buf_org;(f).buf_end=(f).buf_org-1;(f).bit=0;*(f).buf=0;}while(0)

/* macros for allocating memory with exiting if an error occurs */
#define free2(a)   do{if(a){free(a);(a)=NULL;}}while(0)
#define free2d(a)  do{if(a){free2((a)[0]);free(a);(a)=NULL;}}while(0)
#define fclose2(a) do{if(a){fclose(a);(a)=NULL;}}while(0)

/* macros for reading integers from file, d=0 read one-line, d=1 read all file */
//#define   ARY_SCAN(k,a,fp,d)  do{(k)=0;do{do{FILE2_read_##a(&(fp));}while((FILE_err&6)==8-(d)*4);if(FILE_err&(4-2*(d)))break;(k)++;}while((FILE_err&(3-(d)))==0);}while(0)
#define   ARY_SCAN(num,type,fp,d)  do{(num)=0;do{do{FILE2_read_##type(&(fp));}while((FILE_err&((d)*5))==5);if(RANGE(5+(int)(d),FILE_err,6))break;(num)++;}while((FILE_err&(3-(int)(d)))==0);}while(0)
#define   ARY_READ_LINE(f,type,num,fp)  do{(num)=0;do{(f)[num]=FILE2_read_##type(&(fp));(num)+=(FILE_err&4)?0:1;}while((FILE_err&3)==0);}while(0)
#define   ARY_READ(f,type,num,fp)  do{FLOOP(common_size_t,0,(size_t)num){do{(f)[common_size_t]=FILE2_read_##type(&(fp));}while((FILE_err&6)==4);if(FILE_err&2)break;}}while(0)
#define   ARY_LOAD(f,type,num,fname,d,x)  do{FILE2_open(common_FILE2,fname,"r",x);ARY_SCAN(num,type,common_FILE2,d);malloc2(f,(num)+1,x);FILE2_reset(&common_FILE2);ARY_READ(f,type,num,common_FILE2);FILE2_close(&common_FILE2);}while(0)
#define   ARY_WRITE(n,f,num,q,x)  do{fopen2(common_FILE,n,"w",x);ARY_FPRINT(common_FILE,f,0,num,q);fclose(common_FILE);}while(0)

#ifndef MQUE_ONEMORE
 #define MQUE_ONEMORE 1
#endif

/* macros for generalized queue; end mark is necessary for INTSEC */
#define   MQUE_FLOOP(V,z)    for((z)=(V).v;(z)<(V).v+(V).t ; (z)++)
#define   MQUE_FLOOP2(V,z,zz)    for((z)=(V).v,(zz)=(z)+(V).t;(z)<(zz) ; (z)++)
#define   MQUE_BLOOP(V,z)    for((z)=(V).v+(V).t-1;(z)>=(V).v ; (z)--)
#define   MQUE_SLOOP(V,z)    for((z)=(V).v+(V).s;(z)<(V).v+(V).t ; (z)++)
#define   MQUE_SLOOP2(V,z,zz)    for((z)=(V).v+(V).s,(zz)=(z)+(V).t;(z)<(zz) ; (z)++)
#define   MQUE_SBLOOP(V,z)    for((z)=(V).v+(V).t-1;(z)>=(V).v+(V).s ; (z)--)

#ifdef _cplusplus_
 #define   MQUE_FLOOP_(V,z,s)    for((z)=(V).v ; (char *)(z)<((char *)(V).v)+(V).t*(s) ; (z)=(typeof(z))(((char *)(z))+(s)))
#else
 #define   MQUE_FLOOP_(V,z,s)    for((z)=(V).v ; (char *)(z)<((char *)(V).v)+(V).t*(s) ; (z)=(void *)(((char *)(z))+(s)))
#endif

#define   MQUE_MLOOP(V,z,M)  for((z)=(V).v; *((QUEUE_INT *)z)<(M) ; (z)++)

///  !!! errr MQUE_INTSEC !!!!!
#define   MQUE_INTSEC(f,U,V)  do{\
common_INT=0;(f)=0;\
FLOOP(common_INT2,0,(U).t){\
 while(*((QUEUE_INT *)(&((V).v[common_INT])))<*((QUEUE_INT *)(&((U).v[common_INT2])))&&common_INT<(V).t){ \
  if (++common_INT >= (V).t) break;\
 }if(*((QUEUE_INT *)(&((V).v[common_INT])))==*((QUEUE_INT *)(&((U).v[common_INT2]))))(f)++;\
}}while(0)
#define   MQUE_UNION(f,U,V)   do{MQUE_INTSEC(f,U,V);(f)=(U).t+(V).t-(f);}while(0)
#define   MQUE_DIF(f,U,V)     do{MQUE_INTSEC(f,U,V);(f)=(U).t+(V).t-(f)-(f);}while(0)
#define   MQUE_RM_DUP(V)      do{\
if((V).t>1){\
 common_INT=1;\
 FLOOP(common_INT2,1,(V).t){\
  if ( *((QUEUE_INT *)(&((V).v[common_INT2-1]))) != *((QUEUE_INT *)(&((V).v[common_INT2]))) ) (V).v[common_INT++]=(V).v[common_INT2];\
  } (V).t=common_INT;\
 }\
}while(0)

#define   MQUE_UNIFY(V,a)      do{\
if((V).t>1){\
 common_INT=0;\
 FLOOP(common_INT2,1,(V).t){\
  if ( *((QUEUE_INT *)(&((V).v[common_INT2-1]))) != *((QUEUE_INT *)(&((V).v[common_INT2]))) ) (V).v[++common_INT]=(V).v[common_INT2];\
  else *((a*)(((QUEUE_INT *)(&((V).v[common_INT2])))+1)) += *((a*)(((QUEUE_INT *)(&((V).v[common_INT2])))+1));\
 } (V).t=common_INT+1;\
}}while(0)




#ifndef VEC_VAL
 #ifdef VEC_VAL_CHAR
  #define VEC_VAL char
  #define VEC_VAL2 LONG
  #define VEC_VAL_END 128
  #define VEC_VAL2_END LONGHUGE
  #define VEC_VALF "%hhd"
 #elif defined(VEC_VAL_UCHAR)
  #define VEC_VAL unsigned char
  #define VEC_VAL2 LONG
  #define VEC_VAL_END 256
  #define VEC_VAL2_END LONGHUGE
  #define VEC_VALF "%hhu"
 #elif defined(VEC_VAL_INT)
  #define VEC_VAL int
  #define VEC_VAL2 LONG
  #define VEC_VAL_END INTHUGE
  #define VEC_VAL2_END LONGHUGE
  #define VEC_VALF "%d"
 #else
  #define VEC_VAL double
  #define VEC_VAL2 double
  #define VEC_VAL_END DOUBLEHUGE
  #define VEC_VAL2_END DOUBLEHUGE
  #define VEC_VALF "%f"
 #endif
#endif

#ifndef VEC_ID
 #ifdef VEC_ID_LONG
  #define VEC_ID LONG
  #define VEC_ID_END LONGHUGE
  #define VEC_IDF LONGF
 #else
  #define VEC_ID int
  #define VEC_ID_END INTHUGE
  #define VEC_IDF "%d"
 #endif
#endif

/* vector */
typedef struct {
  unsigned char type;  // mark to identify type of the structure
  VEC_VAL *v;
  VEC_ID end;
  VEC_ID t;
} VEC;

extern VEC INIT_VEC;
extern PERM common_PERM, *common_PERMp;
extern VEC_VAL common_VEC_VAL, *common_VEC_VALp;
extern VEC_ID common_VEC_ID;

/* tranpose the matrix ; counting/transpose/memory_allocate */
#define MQUE_DELIVERY_CNT(c,jump,f,y,M) do{ \
FLOOP(common_VEC_ID, 0, (f).t){ \
 MQUE_MLOOP( (f).v[common_VEC_ID], y, M){ \
  if( (c)[*((QUEUE_INT *)y)] == 0 ) ARY_INS(jump, *((QUEUE_INT *)y)); \
  (c)[*((QUEUE_INT *)y)]++; \
 } \
}}while(0)
#define MQUE_DELIVERY(occ,jump,f,y,M)  do{ \
FLOOP (common_VEC_ID, 0, (f).t){ \
 MQUE_MLOOP ((f).v[common_VEC_ID], y, M){ \
  if( (occ)[*((QUEUE_INT *)y)].t == 0 ) ARY_INS( jump, *((QUEUE_INT *)y)); \
  ARY_INS( (occ)[*((QUEUE_INT *)y)], common_VEC_ID); \
 } \
}}while(0)

  // allocate QUEUE of QUEUE, according to rows (cells of "ext" will be add to each row)
#ifdef _cplusplus_
#define MQUE_ALLOC(Q,rows,rowt,unit,ext,x)   do{ \
common_size_t=0; \
FLOOP (common_VEC_ID, 0, rows) common_size_t += rowt[common_VEC_ID]; \
calloc2 (Q, (rows)+1, x); \
malloc2 (common_pnt, (common_size_t+(rows)*(ext)+2)*((unit)<sizeof(QUEUE_INT)?sizeof(QUEUE_INT):(unit)), {free(Q);x;}); \
FLOOP (common_VEC_ID, 0, rows){ \
 (Q)[common_VEC_ID].end = rowt[common_VEC_ID]; \
 (Q)[common_VEC_ID].v = (typeof((Q)[common_VEC_ID].v))common_pnt; \
 common_pnt += ((unit)<sizeof(QUEUE_INT)?sizeof(QUEUE_INT):(unit))*(rowt[common_VEC_ID]+(ext));\
}}while(0)
#else
#define MQUE_ALLOC(Q,rows,rowt,unit,ext,x)   do{ \
common_size_t=0; \
FLOOP (common_VEC_ID, 0, rows) common_size_t += rowt[common_VEC_ID]; \
calloc2 (Q, (rows)+1, x); \
malloc2 (common_pnt, (common_size_t+(rows)*(ext)+2)*((unit)<sizeof(QUEUE_INT)?sizeof(QUEUE_INT):(unit)), {free(Q);x;}); \
FLOOP (common_VEC_ID, 0, rows){ \
 (Q)[common_VEC_ID].end = rowt[common_VEC_ID]; \
 (Q)[common_VEC_ID].v = (void *)common_pnt; \
 common_pnt += ((unit)<sizeof(QUEUE_INT)?sizeof(QUEUE_INT):(unit))*(rowt[common_VEC_ID]+(ext));\
}}while(0)
#endif

/**********************************************************/
/* macro for execute the programs as subroutine           */
/* do not nest the macro; <= variable will be used doubly */
/**********************************************************/
    // execute the main routine of the program as a subroutine
#define EXECSUB(func,flag,err,fmt,...) do{\
 sprintf (common_comm, fmt, __VA_ARGS__);\
 if ( flag&SHOW_MESSAGE ) print_err ("%s\n", common_comm);\
 common_INT = string_decompose (common_argv, common_comm, ' ', 100);\
 if ( func (common_INT, common_argv) ) err;\
}while(0)
    // remove a file on the specified directory
#define REMOV(dir,fn) do{ sprintf(common_comm,"%s%s",dir,fn); remove (common_comm); }while(0)
#define RENAM(dir,fn,fn2) do{ sprintf(common_comm,"%s%s",dir,fn); sprintf(common_comm2,"%s%s",dir,fn2); rename (common_comm, common_comm2); }while(0)
#define MREMOV(dir,...) mremove_(dir, __VA_ARGS__, NULL, NULL)


/*********************************************************/


#define SHOW_MESSAGE 1   // not print messages
#define SHOW_PROGRESS 2  // show progress of the computation
#define APPEND_OUTPUT 4  // append the output to existing file
#define NO_REMOVE 8      // do not remove the temporal files
#define LOAD_PERM 64      // permute the nodes/items by something
#define LOAD_RM_DUP 128   // duplicate items in each row, for loading data
#define LOAD_INCSORT 256  // sort rows in increasing order, for loading data
#define LOAD_DECSORT 512  // sort rows in decreasing order, for loading data
#define LOAD_ELE 1024     // load tuple-list file 
#define LOAD_TPOSE 2048     // transpose the file when load 
#define LOAD_DBLBUF 4096     // transpose the file when load 
#define LOAD_WSORT 8192     // sort rows by their weights 
#define LOAD_SIZSORT 16384     // sort rows by their sizes 
#define LOAD_DECROWSORT 32768   // sort rows in decreasing order
#define LOAD_NUM 65536   // read #columns, #rows and #elements from the 1st line of the file

#define LOAD_EDGEW 131072     // read edge weight
#define LOAD_ARCW 262144     // read arc weight
#define LOAD_NODEW 524288     // read node weight
#define LOAD_BIPARTITE 1048576     // read bipartite graph
#define LOAD_EDGE 2097152     // read edge 
#define LOAD_ARC 4194304     // read arc 
#define LOAD_GRAPHNUM 8388608   // read #vertices and #edges from the 1st line of the file
#define LOAD_ID1 16777216       // node ID begins from 1
#define LOAD_COMP 33554432   // read as a complement  

#define FILE_COUNT_ROWT 32   // count size of each row
#define FILE_COUNT_CLMT 64   // count size of each column
#define FILE_COUNT_NUM LOAD_NUM   // read #columns, #rows and #elements
#define FILE_COUNT_GRAPHNUM LOAD_GRAPHNUM   // read #vertices and #edges

#define FILE2_BUFSIZ 16384
typedef struct {   // structure for fast file reader routines
  unsigned char type;  // type definition
  FILE *fp;
  char *buf_org, *buf, *buf_end;   // head/current/tail of buffer
  char bit;
} FILE2;
extern FILE2 INIT_FILE2, common_FILE2;

void FILE2_flush (FILE2 *fp);
void FILE2_close (FILE2 *fp);
void FILE2_closew (FILE2 *fp);
void FILE2_reset (FILE2 *fp);

int FILE2_getc (FILE2 *fp);
void FILE2_puts (FILE2 *fp, char *s);
void FILE2_putc (FILE2 *fp, char c);
int FILE2_getbit (FILE2 *fp);
void FILE2_putbit (FILE2 *fp, int a);

/* message output */
//void print_mes (int flag, char *mes, ...);
//void print_err (char *mes, ...);

/* free many pointers*/
void mfree_(void *x, ...);
void mfree2_(void *x, ...);

/* remove many files */
void mremove_ (char *x, ...);

/* compute the minimum prime no less than n */
LONG min_prime (LONG n);

/* decompose the string by separator, and set v[i] to each resulted string.
  consecutive separators are regarded as one separator. */
int string_decompose (char **v, char *s, char sep, int max);

/* make two random numbers under normal distribution N(0,1) */
void rand_mk_2normal (double *a, double *b);
/* make a random point on a supersphere of d-dim., and set to double array already allocated */
void rand_d_gaussian (double *p, int d);
void rand_sphere (double *p, int d);



/* Read an integer from the file and return it,
    with read through the non-numeric letters.
   If it reaches to the end-of-file, then set FILE_err=2, if it reads a 
   newline, then set FILE_err=1.
   If read either the end-of-file or newline before reading an integer,
   return -1 */
FILE_LONG FILE2_read_int (FILE2 *fp);
double FILE2_read_double (FILE2 *fp);
WEIGHT FILE2_read_WEIGHT (FILE2 *fp);
LONG FILE2_read_mpint (FILE2 *fp, int max, int skip);  // multi-precise
LONG FILE2_read_intbit (FILE2 *fp, int k);  // read k bits as integer


/* fast file routine, print number, c is the char to be printed preceding to the number
   if c==0, nothing will be printed preceding the number
   if len<0 then the #digits following '.' does not change (filed by '0')  */
void FILE2_print_int (FILE2 *fp, LONG n, char c);
void FILE2_print_real (FILE2 *fp, double n, int len, char c);
void FILE2_print_WEIGHT (FILE2 *fp, WEIGHT w, int len, char c);
void FILE2_printf (FILE2 *fp, char *mes, ...);
void FILE2_print_mpint (FILE2 *fp, LONG a, int max, int skip);  // multi-precise
void FILE2_print_intbit (FILE2 *fp, LONG a, int k);  // print lowest k-bits

/* print a real number in a good style */
void fprint_real (FILE *fp, double f);
void print_real (double f);
void fprint_WEIGHT (FILE *fp, WEIGHT f);
void print_WEIGHT (WEIGHT f);

#define FILE_COUNT_INT VEC_ID
#define FILE_COUNT_INTF VEC_IDF
typedef struct {
  int flag;
  FILE_COUNT_INT clms, rows, eles, clm_end, row_end, row_btm, clm_btm; // #rows, #column, #elements, minimum elements
  FILE_COUNT_INT row_min, row_max, clm_min, clm_max;  // maximum/minimum size of column
  FILE_COUNT_INT *rowt, *clmt;   // size of each row/clmn
  WEIGHT total_rw, total_cw, *rw, *cw;  // WEIGHTs for rows/columns ... reserved.
  FILE_COUNT_INT rw_end, cw_end;
  PERM *rperm, *cperm;   // permutation (original->internal) of rows and columns
} FILE_COUNT;

extern FILE_COUNT INIT_FILE_COUNT;

/* count the clms, rows, items, each row size, each column size  */
/* file types can be, array list and element list*/
/* support transpose */
FILE_COUNT FILE2_count (FILE2 *fp, int flag, int skip_rows, int int_rows, int skip_clms, int int_clms, FILE_COUNT_INT row_limit);

/******************* integer array routines **************************/

/******************************* permutation routines ****************/
/* permutation is given by an integer array  */

/* sort an array of size "siz", composed of a structure of size "unit" byte
 in the order of perm */
/* use temporary memory of siz*unit byte */
//void perm_struct (void *a, int unit, int *perm, size_t siz);

/* SLIST:very simple one-sided list */
void SLIST_init (int *l, int num, int siz);
void SLIST_end (int *l);
#define SLIST_INS(l,m,e) (l[e]=l[m],l[m]=e);

#define QQSORT_ELE(a,x) ((a *)(&(common_pnt[(*((PERM *)(x)))*common_int])))
#define QQSORT_ELEt(a,x) (((a *)&(common_pnt[(*((PERM *)(x)))*common_int]))->t)
/* quick sort macros */
#define QSORT_TYPE(a,b) \
b common_##a; \
int qsort_cmp_##a (const void *x, const void *y){ \
 if ( *((b *)x) < *((b *)y) ) return (-1); else return ( *((b *)x) > *((b *)y) ); \
} \
int qsort_cmp__##a (const void *x, const void *y){ \
 if ( *((b *)x) > *((b *)y) ) return (-1); else return ( *((b *)x) < *((b *)y) ); \
} \
int qqsort_cmp_##a (const void *x, const void *y){ \
 b *xx=QQSORT_ELE(b,x), *yy=QQSORT_ELE(b,y); \
 if ( *xx < *yy ) return (-1); \
  else return ( *xx > *yy ); \
} \
int qqsort_cmp__##a (const void *x, const void *y){ \
 b *xx=QQSORT_ELE(b,x), *yy=QQSORT_ELE(b,y); \
 if ( *xx > *yy ) return (-1); \
 else return ( *xx < *yy ); \
} \
void qsort_##a (b *v, size_t siz, int unit){ \
 if ( unit == 1 || unit==-1 ) unit *= sizeof (b); \
 if (unit<0) qsort (v, siz, -unit, qsort_cmp__##a); else qsort (v, siz, unit, qsort_cmp_##a); \
} \
void qsort_perm__##a (b *v, size_t siz, PERM *perm, int unit){ \
 ARY_INIT_PERM(perm,(INT)siz); \
 if ( unit == 1 || unit==-1 ) unit *= sizeof (b); \
 common_int=MAX(unit,-unit); common_pnt=(char *)v; \
 if (unit<0) qsort (perm, siz, sizeof(PERM), qqsort_cmp__##a); \
 else qsort (perm, siz, sizeof(PERM), qqsort_cmp_##a); \
} \
PERM *qsort_perm_##a (b *v, size_t siz, int unit){ \
 PERM *perm; malloc2(perm, siz, EXIT0); \
 qsort_perm__##a (v, siz, perm, unit); return(perm); \
} \
size_t bin_search_##a (b *v, b u, size_t siz, int unit){ \
 size_t s=0, t; \
 b n; \
 if ( unit == 1 ) unit *= sizeof (b); \
 while (1){ \
   if ( siz-s <= 2 ){ \
     if ( u <= *((b *)(((char *)v)+unit*s)) ) return (s);\
     if ( siz == s+1 || u <= *((b *)(((char *)v)+unit*(s+1))) ) return (s+1);\
     return (s+2);\
   }\
   t = (s+siz) /2; \
   n = *((b *)(((char *)v)+unit*t));\
   if ( u < n ) siz = t; else s = t;\
 }\
}
// bin search returns the position that is equal to u, or the smallest in larger's

#define QSORT_TYPE_HEADER(a,b) \
extern b common_##a; \
int qsort_cmp_##a (const void *x, const void *y); \
int qsort_cmp__##a (const void *x, const void *y); \
int qqsort_cmp_##a (const void *x, const void *y); \
int qqsort_cmp__##a (const void *x, const void *y); \
void qsort_##a(b *v, size_t siz, int unit); \
void qsort_perm__##a (b *v, size_t siz, PERM *perm, int unit); \
PERM *qsort_perm_##a (b *v, size_t siz, int unit); \
size_t bin_search_##a (b *v, b u, size_t siz, int unit);

QSORT_TYPE_HEADER(int, int)
QSORT_TYPE_HEADER(uint, unsigned int)
QSORT_TYPE_HEADER(double, double)
QSORT_TYPE_HEADER(char, char)
QSORT_TYPE_HEADER(uchar, unsigned char)
QSORT_TYPE_HEADER(short, short)
QSORT_TYPE_HEADER(ushort, unsigned short)
QSORT_TYPE_HEADER(WEIGHT, WEIGHT)
QSORT_TYPE_HEADER(LONG, LONG)
QSORT_TYPE_HEADER(VEC_ID, VEC_ID)
QSORT_TYPE_HEADER(VEC_VAL, VEC_VAL)
QSORT_TYPE_HEADER(VEC_VAL2, VEC_VAL2)
QSORT_TYPE_HEADER(FILE_COUNT_INT, FILE_COUNT_INT)

int qsort_cmp_VECt (const void *x, const void *y);
int qsort_cmp__VECt (const void *x, const void *y);
void qsort_VECt (VEC *v, size_t siz, int unit);
int qqsort_cmp_VECt (const void *x, const void *y);
int qqsort_cmp__VECt (const void *x, const void *y);
void qsort_perm__VECt (VEC *v, size_t siz, PERM *perm, int unit);
PERM *qsort_perm_VECt (VEC *v, size_t siz, int unit);

/* swap macro for integer, double, char, and pointer */
#define SWAP_INT(a,b)  (common_LONG=a,a=b,b=common_LONG)
#define SWAP_DOUBLE(a,b)  (common_double=a,a=b,b=common_double)
#define SWAP_size_t(a,b)  (common_size_t=a,a=b,b=common_size_t)
#ifdef _cplusplus_
 #define SWAP_PNT(a,b)  (common_pnt=(typeof(common_pnt))a,a=(typeof(a))b,b=(typeof(b))common_pnt)
#else
 #define SWAP_PNT(a,b)  (common_pnt=(void *)a,a=(void *)b,b=(void *)common_pnt)
#endif

#define SWAP_VEC_ID(a,b)  (common_VEC_ID=a,a=b,b=common_VEC_ID)
#define SWAP_FILE_COUNT_INT(a,b)  (common_FILE_COUNT_INT=a,a=b,b=common_FILE_COUNT_INT)

/* sort index(int)/WEIGHT array and return the indices of the result */
/* perm[i*2] := rank of ith cell */
/* perm[i*2+1] := index of ith smallest cell */
/*  flag!=NULL  => opposite direction sort  */

//int *int_index_sort (int *w, size_t siz, int flag);
//int *WEIGHT_index_sort (WEIGHT *w, size_t siz, int flag);

/* radix sort for array of structures, by their integer members 
  ranging from mm to m */
/* sort array "perm" according to (int/void*) array "a".
   if perm==NULL, allocate memory and for perm */
/* return the permutation array of the result of the sorting
   in the decreasing order if unit<0 (unimplemented!!!) */
int *radix_sort (void *a, size_t siz, int unit, int mm, int m, int *perm);

/* permutate structure array *tt of unit size unit of size num, according to perm array *pp */
/* num has to be <INTHUGE/2 */
/* unit<0 means decreasing order (un-implemented!!!) */
void structure_permute (void *tt, int unit, int num, void *pp, int weight_siz);

#ifdef STDLIB2_RADIX_SORT     // radix sort with 1M byte static memory
void intarray_sort_iter (unsigned int *a, size_t siz, int unit);
void intarray_sort_iter_ (unsigned int *a, size_t siz, int unit);
void intarray_sort (unsigned *a, size_t siz, int unit);
#endif

  /* bitmasks, used for bit operations */
extern int BITMASK_UPPER1[32];
extern int BITMASK_UPPER1_[32];
extern int BITMASK_LOWER1[32];
extern int BITMASK_LOWER1_[32];
extern int BITMASK_1[32];
extern int BITMASK_31[32];
extern int BITMASK_16[8];
extern int BITMASK_UPPER16[8];
extern int BITMASK_LOWER16[8];
extern int BITMASK_FACT16[8];

/* bit operations */
#define BIT_SET(v,x)  (v[x/32]|=BITMASK_1[x%32])
#define BIT_DEL(v,x)  (v[x/32]&=BITMASK_31[x%32])
#define ISBIT(v,x)    (v[x/32]&BITMASK_1[x%32])

#endif


/******************************************/
/*  ==== terminology for comments  ====
  range check: to check the input parameter is valid, or in the valid range.
  If a function does not have this, its comment has "no range check"
  */
   
/*  ====  rules for the name of functions and routines  ====
  init: initialization for an object, structure, etc. memory is allocated 
   if needed.
  end: compared to init, termination of structures, etc.
   free allocated memory if it exists, but not free itself.
   different from ordinary new, create, del, free.

  cpy: copy an object without allocating memory
  dup: make a duplication of an object with allocating new memory

  new: new. allocate memory for new object. also used for re-allocation from 
   the list of deleted objects
  del: delete. free memory, or insert it to the list of deleted objects

  ins : insert. insert an element (active, not deleted) to an object, possible
     at the given position.
  out : extract. extract an element, possibly specified, from an object.
        it will be not deleted.
  rm  : extract, and delete
  rmall: delete all (specified) elements of an object
  mk  : make. new+insB
  mv  : move. move the elements from an object to another,
         or change the position.

  update : update an object, possibly of specified position, to the exact,
           or to the current one.
  chg    : change the status of an object to the specified one.
  
  prv: point the previous element
  nxt: point the next element
  end: the maximum size (allocated size) of an array, etc.
  num: the (current) number of elements in an array, etc.
  kth: for the use of "first k object"
  tkth: for the use of "first k object from the end". to kth.
  rnd: random
  print: print structure and etc.
  find: search or find an specified element from the set of structures, etc.
*/

/*  ==== rules for the name of variables ====
  - use i or j for the counter in loops
  - use e for elements
  - use k for the specified rank
  - use f or flag for flagment
  - use m for maximum value or minimum value
  - use c for counters
*/


