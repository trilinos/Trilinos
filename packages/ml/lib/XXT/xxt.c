/*************************************xxt.c************************************
Module Name: xxt
Module Info:

author:  Henry M. Tufo III
e-mail:  hmt@asci.uchicago.edu
contact:
+--------------------------------+--------------------------------+
|MCS Division - Building 221     |Department of Computer Science  |
|Argonne National Laboratory     |Ryerson 152                     |
|9700 S. Cass Avenue             |The University of Chicago       |
|Argonne, IL  60439              |Chicago, IL  60637              |
|(630) 252-5354/5986 ph/fx       |(773) 702-6019/8487 ph/fx       |
+--------------------------------+--------------------------------+

Last Modification: 3.27.00
**************************************xxt.c***********************************/


/*************************************xxt.c************************************
NOTES ON USAGE: 

**************************************xxt.c***********************************/


/*************************************xxt.c************************************
FILE FORMAT: 
------------------------------ Begin File -------------------------------------
------------------------------ End   File -------------------------------------
Note: 
**************************************xxt.c***********************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#if   defined NXSRC
#include <nx.h>
#elif defined MPISRC
#include <mpi.h>
#endif

#include "const.h"
#include "types.h"
#include "comm.h"
#include "error.h"
#include "ivec.h"
#include "bss_malloc.h"
#include "queue.h"
#include "gs.h"
#include "ml_include.h"
#include "xxt.h"

#define MAX_REA_NAME   80
#define MAX_READ_BUF   1000
#define STD_READ_BUF   80
#define MAX_SEPARATORS 4096

#define LEFT  0
#define RIGHT 0


typedef struct xxt_solver_info {
  int n,m;
  int num_sep_total;
  int *separators;
  int *col_sz;
  int *col_indices;
  REAL **col_vals;
  REAL *x;
  REAL *solve_uu, *solve_w;
  gs_ADT gs_handle;
  vfp A_matvec;
  void *grid_data;
} xxt_info;


struct xxt_CDT{
  int id;
  xxt_info *info;
};

static int n_xxt_handles=0;
static int fn_xxt_handles=0;
static xxt_ADT fhandles[10];

/* slated for extinction */
/* rea file handle */
char dir_name[MAX_REA_NAME+5];
char rea_name[MAX_REA_NAME+5];
char map_name[MAX_REA_NAME+5];
char sep_name[MAX_REA_NAME+5];

/* number of spectral elements (SE) over all P and # I own */
static int nel_global, nel_local;

/* number of dof w/n_global = nvu-nvo */
/* number of unique vertices I own = n - # outflow I own */
/* after separators read in equal to E coarse dimension n */
static int n_global, n_local;

/* depth of separator tree and max number of processors we can run on */
static int depth, max_proc;

/* number of corners a SE has ... 4 or 8 only! */
static int nc;

/* number of vertices = nel_global*nc */
static int nv;

/* number of unique vertices = nv - {copies over all SE} */
static int nvu;

/* number of outflow vertices */
static int nvo;

/* holds element vertex global numbering - size=n_global*nv */
static int *vertex;

/* holds element to processor map */
static int *map;

/* is vcor in effect? */
static int vcor;

/* hold E coarse matrix - note it's fortran indexed */
static int *a_local_to_global;
static int a_n, a_m;

/* X info ... it's nxm */
static int num_sep_total;
static int *separators;
static int n, m;
static int *col_sz;
static int *col_indices;
static REAL **col_vals;
static REAL *x;
static REAL *solve_uu, *solve_w;

/* number of things in each separator group */
static int *segments;

/* for xxt staged comm routine */
static int *stages;

/* xxt info */
static int level;
static int xxt_nnz;
static int xxt_max_nnz;

/* new xxt gen info */
static int *seps, *nsep, *lnsep, *fo, nfo=0;

static gs_ADT gs_handle;
static vfp A_matvec;
static void *grid_data;


#ifdef r8
double ddot(int n, double *x, int incx, double *y, int incy);
void daxpy(long int n, double da, double *dx, long int incx, 
           double *dy, long int incy);
#else
float sdot(int n, float *x, int incx, float *y, int incy);
void  saxpy(long int n, float da, float *dx, long int incx, 
           float *dy, long int incy);
#endif

/* prototypes */
double sqrt(double);
static void fix_a_map(int ub);
static void set_file_names(void);
static void xxt_set_local_map(int *, int, int, void *);
static void xxt_setup(int n);
static void xxt_load_firing_order(void);
static void xxt_solve (REAL *);
static void xxt_generate(void);
static void xxt_elm_to_proc();

/*new protos */
static void check_init(void);
static void det_separators(xxt_ADT xxt_handle,int *local2global, int n, int m,
			   void *mylocmatvec, void *grid_data);

void dump_vec(REAL *v, int n, int tag, int tag2, char * s);
void dump_ivec(int *v, int n, int tag, int tag2, char * s);

/*************************************xxt.c************************************
Function:  ML_Do_CoarseDirect

Input : 
Output: 
Return: 
Description:  
**************************************xxt.c***********************************/
int
ML_Do_CoarseDirect(xxt_ADT xxt_handle, double *in, double *out)
{
int i,index,col;
int *col_map;


#ifdef DEBUG
  error_msg_warning("ML_Do_CoarseDirect() :: begin\n");
#endif

  /*
  if (id!=xxt) 
    {error_msg_warning("ML_Do_CoarseDirect ... hey ML lost xxt handle!\n");}
  */

  /*
  rvec_zero(out,a_n);
  col_map = separators;
  for (i=0;i<a_n;i++)
    {
      col = col_map[i];
      index=ivec_linear_search(col, a_local_to_global, a_m);
      if (index==-1)
	{error_msg_fatal("didn't find %d\n",col);}
      out[i] = in[index];
    }

  rvec_copy(in,out,a_n);

  printf("in : ");
  for (i=0;i<a_n;i++)
    {printf(" %.3f ",in[i]);}
  printf("\n");
  */

  xxt_solve (in);

#ifdef DEBUG
  printf("%2d :: RES2 ::",my_id);
  for (i=0; i<n; i++)
    {printf("%5.3f ",1.0-in[i]);}
    /* {printf("%5.3f ",1.0*i-in[i]);} */
  printf("\n");
#endif

  /*
  rvec_zero(out,a_n);
  col_map = separators;
  for (i=0;i<a_n;i++)
    {
      col = col_map[i];
      index=ivec_linear_search(col, a_local_to_global, a_m);
      if (index==-1)
	{error_msg_fatal("didn't find %d\n",col);}
      out[index] = in[i];
    }
  
  printf("out : ");
  for (i=0;i<a_n;i++)
    {printf(" %.3f ",out[i]);}
  printf("\n");
  */

  rvec_copy(out,in,a_n);

#ifdef DEBUG
  error_msg_warning("ML_Do_CoarseDirect() :: end\n");
#endif

  return(TRUE);
}


/*************************************xxt.c************************************
Function: ML_Gen_CoarseDirect

Input : 
Output: 
Return: 
Description: get A_local, local portion of global coarse matrix which 
is a row dist. nxm matrix w/ n<m.
   o my_ml holds address of ML struct associated w/A_local and coarse grid
   o local2global holds global number of column i (i=0,...,m-1)
   o local2global holds global number of row    i (i=0,...,n-1)
   o mylocmatvec performs A_local . vec_local (note that gs is performed using 
   gs_init/gop).

mylocmatvec = my_ml->Amat[grid_tag].matvec->external;
mylocmatvec (void :: void *data, double *in, double *out)
**************************************xxt.c***********************************/
int
ML_Gen_CoarseDirect(xxt_ADT xxt_handle, int *local2global, void *g_data, 
		    int n, int m, void *mylocmatvec)
{
int i,alpha,w;
int op[] = {GL_ADD,0};
REAL tmp_uu[100];


#ifdef DEBUG
  error_msg_warning("ML_Gen_CoarseDirect() :: begin\n");
#endif

  /* hold grid date for matvec routine */
  grid_data = g_data;

  /* copy over local to global map */
  A_matvec = (vfp) mylocmatvec;
  xxt_set_local_map(local2global, n, m, mylocmatvec);

  /*
  for (i=0; i<100; i++)
    {gs_gop(gs_handle,tmp_uu,"+\0");}
  */

#ifdef DEBUG
  error_msg_warning("done gs_gop test 1\n");
#endif

  /* ok do it */
  xxt_setup(n);

  /*
  for (i=0; i<100; i++)
    {gs_gop(gs_handle,tmp_uu,"+\0");}
  */

  /* free communication handle */
  /* gs_free(gs_handle); */

#ifdef DEBUG
  error_msg_warning("done gs_gop test 2\n");
#endif

  bss_stats(); 
  perm_stats(); 

  alpha=1;
  for (i=0; i<100; i++)
    {
      gs_gop(gs_handle,tmp_uu,"+\0"); 
      giop(&alpha,&w,1,op); 
      fflush(NULL);
    }

#ifdef DEBUG
  error_msg_warning("ML_Gen_CoarseDirect() :: end\n");
#endif

  return(TRUE);
}



/*************************************xxt.c************************************
Function: 

Input : 
Output: 
Return: 
Description:  
**************************************xxt.c***********************************/
static 
void
xxt_setup(int N)
{
  int i,j,k;
  int max, rct, flag;
  REAL *ptr;
  

#ifdef DEBUG
  error_msg_warning("xxt_setup() :: begin\n");
  error_msg_warning("xxt_setup() :: hardwired vcor to FALSE\n");
#endif

  vcor = FALSE;
  n = N;
  
  /* load up firing order */
  xxt_generate();

  max = xxt_nnz;
  j=GL_MAX;
  giop(&max,&i,1,&j);
  if (!my_id) {printf("%d :: max   xxt_nnz=%d\n",my_id,max);}

  max = xxt_nnz;
  j=GL_ADD;
  giop(&max,&i,1,&j);
  if (!my_id) {printf("%d :: total xxt_nnz=%d\n",my_id,max);}

  if (!my_id)
    {
      printf("%d :: xxt   ratio  =%g\n",my_id,max/(pow(1.0*n_global,1.5)));
      printf("%d :: xxt   ratio  =%g\n",my_id,max/(pow(1.0*n_global,1.6667)));
      printf("%d :: xxt   max_buf=%d\n",my_id,xxt_max_nnz);
    }

#ifdef DEBUG
  error_msg_warning("xxt_setup() :: end\n");
#endif
}



/*************************************xxt.c************************************
Function: 

Input : 
Output: 
Return: 
Description:  

X is nxm
E is nxn

nxm mxn nx1 = nx1


**************************************xxt.c***********************************/
static
void
xxt_solve (register REAL *uc)
{
  register int off, len, *iptr;
  register REAL *x_ptr, *uu_ptr;
  REAL zero=0.0;


#ifdef DEBUG
  error_msg_warning("xxt_solve() :: begin\n");
#endif

#ifdef NOT
  if (n!=*n_crs)
    {error_msg_fatal("n=%d, n_crs=%d\n",n,*n_crs);}
#endif

  uu_ptr=solve_uu;
#if   BLAS&&r8
  dcopy(m,&zero,0,uu_ptr,1);
#elif BLAS
  scopy(m,&zero,0,uu_ptr,1);
#else
  rvec_zero(uu_ptr,m);
#endif

  for (x_ptr=x,iptr=col_indices; *iptr!=-1; x_ptr+=len)
    {
      off=*iptr++; len=*iptr++;
#if   BLAS&&r8
      *uu_ptr++ = ddot(len,uc+off,1,x_ptr,1);
#elif BLAS
      *uu_ptr++ = sdot(len,uc+off,1,x_ptr,1);
#else
      *uu_ptr++ = rvec_dot(uc+off,x_ptr,len);
#endif
    }

  uu_ptr=solve_uu;
  if (level) {ssgl_radd(uu_ptr, solve_w, level, stages);}

#if   BLAS&&r8
  dcopy(n,&zero,0,uc,1);
#elif BLAS
  scopy(n,&zero,0,uc,1);
#else
  rvec_zero(uc,n);
#endif

  for (x_ptr=x,iptr=col_indices; *iptr!=-1; x_ptr+=len)
    {
      off=*iptr++; len=*iptr++;
#if   BLAS&&r8
      daxpy(len,*uu_ptr++,x_ptr,1,uc+off,1);
#elif BLAS&&!DELTA
      saxpy(len,*uu_ptr++,x_ptr,1,uc+off,1);
#else
      rvec_axpy(uc+off,x_ptr,*uu_ptr++,len);
#endif
    }

#ifdef DEBUG
  error_msg_warning("xxt_solve() :: end\n");
#endif
}


/*************************************xxt.c************************************
Function: 

Input : 
Output: 
Return: 
Description:  

X is nxm
E is nxn

nxm mxn nx1 = nx1


**************************************xxt.c***********************************/
static
void
xxt_generate()
{
  int i,j,k,l,index;
  int dim, col;
  REAL *u, *uu, *v, *z, *w, alpha, alpha_w;
  int *col_map, *segs;
  int op[] = {GL_ADD,0};
  int max=0;
  int off, len;
  REAL *x_ptr;
  int *iptr, flag;
  int start=0, end=lnsep[0], work;
  int op2[] = {GL_MIN,0};
  int op3[] = {GL_MAX,0};
  int cmin;
  int id, mask;
#ifdef DEBUG
  int lseps[MAX_SEPARATORS];
  int nls=0;
#endif

  
#ifdef DEBUG
  error_msg_warning("xxt_generate() :: begin\n");
#endif

  /* is there a null space? */
  /* LATER add in ability to detect null space by checking alpha */
  for (i=0, j=0; i<=level; i++)
    {j+=nsep[i];}

  m = (vcor) ? (j-1) : (j);

  if (m!=j)
    {printf("xxt_generate() :: MJ m!=j\n");}

  error_msg_warning("xxt_generate() :: X(%d,%d)\n",n,m);    

  /* separator list */
  col_map = separators;

  /* get and initialize storage for x local         */
  /* note that x local is nxm and stored by columns */
  col_sz = (int *) perm_malloc(m*INT_LEN);
  col_indices = (int *) perm_malloc((2*m+1)*sizeof(int));
  col_vals = (REAL **) perm_malloc(m*sizeof(REAL *));
  for (i=j=0; i<m; i++, j+=2)
    {
      col_indices[j]=col_indices[j+1]=col_sz[i]=-1;
      col_vals[i] = NULL;
    }
  col_indices[j]=-1;

  /* size of separators for each sub-hc working from bottom of tree to top */
  /* this looks like nsep[]=segments */
  stages = (int *) perm_malloc((level+1)*INT_LEN);
  segs   = (int *) perm_malloc((level+1)*INT_LEN);
  ivec_zero(stages,level+1);
  ivec_copy(segs,nsep,level+1);
  for (i=0; i<level; i++)
    {segs[i+1] += segs[i];}
  stages[0] = segs[0];

  /* temporary vectors  */
  u  = (REAL *) bss_malloc(n*sizeof(REAL));
  z  = (REAL *) bss_malloc(n*sizeof(REAL));
  v  = (REAL *) bss_malloc(a_m*sizeof(REAL));
  uu = (REAL *) bss_malloc(m*sizeof(REAL));
  w  = (REAL *) bss_malloc(m*sizeof(REAL));

  /* extra nnz due to replication of vertices across separators */
  for (i=1, j=0; i<=level; i++)
    {j+=nsep[i];}

  /* storage for sparse x values */
  n_global = n;
  giop(&n_global,&work,1,op);
  xxt_max_nnz = (2.5*pow(1.0*n_global,1.6667) + j*n/2)/num_nodes;
  x = (REAL *) perm_malloc(xxt_max_nnz*sizeof(REAL));
  xxt_nnz = 0;

  /* LATER - can embed next sep to fire in gs */
  /* time to make the donuts - generate X factor */
#ifdef DEBUG	  
  printf("%2d XXT START :: %d, %d, %d, %d\n",my_id,level,m,n,xxt_max_nnz);
#endif
  for (id=dim=i=j=0;i<m;i++)
    {
      /* time to move to the next level? */
      while (i==segs[dim])
	{
#ifdef SAFE	  
	  if (dim==level)
	    {error_msg_fatal("dim about to exceed level\n"); break;}
#endif
	  if (!id)
	    {
	      cmin = col;
#ifdef DEBUG
	      printf("%2d SEPS %2d :: ",dim,my_id);
	      for (col=0;col<nls;col++)
		{printf("%2d ",lseps[col]);}
	      printf("\n");
	      col=cmin;
	      nls=0;
#endif
	    }

	  stages[dim++]=i;
	  end+=lnsep[dim];

	  cmin = num_nodes>>1;
	  id = my_id;
	  for (col=0;col<(level-dim);col++)
	    {
	      if (id>=cmin)
		{id-=cmin;}
	      cmin>>=1;
	    }
	}
      stages[dim]=i;

      printf("%2d %2d SXXT %2d %2d %2d %2d\n",dim,my_id,start,end,i,segs[dim]);

      /* which column are we firing? */
      /* i.e. set v_l */
      /* use new seps and do global min across hc to determine which one to fire */
      /* col = col_map[i]; */
      (start<end) ? (col=fo[start]) : (col=INT_MAX);
      giop_hc(&col,&work,1,op2,dim); 

      if (col==fo[start])
	{cmin=1;}
      else
	{cmin=0;}
      giop_hc(&cmin,&work,1,op,dim); 
      if (cmin!=1)
	{error_msg_fatal("more then one about to fire %d\n",cmin);}

      error_msg_warning("NXXT %d,%d,%d,%d\n",dim,start,end,col);

#ifdef DEBUG
      if (!id) {lseps[nls++]=col;}
#endif

      /* do I own it? I should */
      rvec_zero(v ,a_m);
      if (col==fo[start])
	{
	  start++;
	  index=ivec_linear_search(col, a_local_to_global, n);
	  if (index!=-1)
	    {
	      v[index] = 1.0; 
	      j++;
	      error_msg_warning("found %d\n",index);    
	    }
	  else
	    {error_msg_fatal("NOT FOUND!\n");}
	}
      else
	{
	  index=ivec_linear_search(col, a_local_to_global, a_m);
	  if (index!=-1)
	    {
	      v[index] = 1.0; 
	    }
	}

      /* distribute that info to others in sub-hc - no longer needed as we're searching cols */
      /* gs_gop_hc(gs_handle,v,"+\0",dim); */

#ifdef DEBUG
      dump_vec(v, a_m, dim, col, "V1"); 
#endif

      /* perform u = A.v_l */
      rvec_zero(u,n);
      A_matvec(grid_data,v,u,num_nodes-1);

#ifdef DEBUG
      dump_vec(u, n, dim, col, "U1"); 
#endif

      /* uu =  X^T.u_l (local portion) */
      /* technically only need to zero out first i entries */
      /* later turn this into an XXT_solve call ? */
      rvec_zero(uu,m);
      x_ptr=x;
      iptr = col_indices;
      for (k=0; k<i; k++)
	{
	  if (!col_vals[k])
	    {error_msg_fatal("x column %d is empty!\n",k);}
	  off = *iptr++;
	  len = *iptr++;

#if   BLAS&&r8
	  uu[k] = ddot(len,u+off,1,x_ptr,1);
#elif BLAS
	  uu[k] = sdot(len,u+off,1,x_ptr,1);
#else
	  uu[k] = rvec_dot(u+off,x_ptr,len);
#endif
	  x_ptr+=len;
	}

#ifdef DEBUG
      dump_vec(uu, 10, dim, col, "UU1"); 
#endif

      /* uu = X^T.u_l (comm portion) */
      /* not needed! rvec_zero(w,m); */
      ssgl_radd  (uu, w, dim, stages);

#ifdef DEBUG
      dump_ivec(stages, 2, dim, col, "UU2"); 
      dump_vec(uu, 10, dim, col, "UU2"); 
#endif

      /* z = X.uu */
      rvec_zero(z,n);
      x_ptr=x;
      iptr = col_indices;
      for (k=0; k<i; k++)
	{
	  if (!col_vals[k])
	    {error_msg_fatal("x column %d is empty!\n",k);}
	  off = *iptr++;
	  len = *iptr++;

#ifdef DEBUG
	  printf("%d:OFFSET %d,%d,%d\n",my_id,k,off,len);
#endif
#if   BLAS&r8
	  daxpy(len,uu[k],x_ptr,1,z+off,1);
#elif BLAS&&!DELTA
	  saxpy(len,uu[k],x_ptr,1,z+off,1);
#else
	  rvec_axpy(z+off,x_ptr,uu[k],len);
#endif
	  x_ptr+=len;
	}

#ifdef DEBUG
      dump_vec(z, n, dim, col, "Z1"); 
#endif

      /* compute v_l = v_l - z */
      rvec_zero(v+n,a_m-n);
#if   BLAS&&r8
      daxpy(n,-1.0,z,1,v,1);
#elif BLAS&&!DELTA
      saxpy(n,-1.0,z,1,v,1);
#else
      rvec_axpy(v,z,-1.0,n);
#endif

#ifdef DEBUG
      dump_vec(v, a_m, dim, col, "V2"); 
#endif

      /* compute u_l = A.v_l */
      gs_gop_hc(gs_handle,v,"+\0",dim);
      rvec_zero(u,n);
      A_matvec(grid_data,v,u,num_nodes-1);

      printf("%2d %2d %2d COL\n",dim,col,my_id);
      fflush(stdout);
#ifdef DEBUG
      dump_vec(u, n, dim, col, "U2"); 
#endif

      /* compute sqrt(alpha) = sqrt(v_l^T.u_l) - local portion */
#if   BLAS&&r8
      alpha = ddot(n,u,1,v,1);
#elif BLAS
      alpha = sdot(n,u,1,v,1);
#else
      alpha = rvec_dot(u,v,n);
#endif
      /* compute sqrt(alpha) = sqrt(v_l^T.u_l) - comm portion */
      grop_hc(&alpha, &alpha_w, 1, op, dim);
#ifdef r8
      alpha = sqrt(alpha);
#else
      alpha = (REAL) sqrt((double)alpha);
#endif

      /* check for small alpha                             */
      /* LATER use this to detect and determine null space */
#ifdef r8
      if (fabs(alpha)<1.0e-14)
	{error_msg_fatal("bad alpha! %g\n",alpha);}
#else
      if (fabs((double) alpha) < 1.0e-6)
	{error_msg_fatal("bad alpha! %g\n",alpha);}
#endif

      printf("%2d %2d %2d 1.0/alpha = %f\n",dim,col,my_id,1.0/alpha);

      /* compute v_l = v_l/sqrt(alpha) */
      rvec_scale(v,1.0/alpha,n);

      /* add newly generated column, v_l, to X */
      flag = 1;
      off=len=0;
      for (k=0; k<n; k++)
	{
	  if (v[k]!=0.0)
	    {
	      len=k;
	      if (flag)
		{off=k; flag=0;}
	    }
	}
      
      printf("%2dOFF %2d %2d %2d %7.4f\n",dim,my_id,i,off,v[off]);

      len -= (off-1);

      if (len>0)
	{
	  if ((xxt_nnz+len)>xxt_max_nnz)
	    {
#ifdef DEBUG	  
	      printf("increasing space for X by 2x!\n");
#endif
	      xxt_max_nnz *= 2;
	      x_ptr = (REAL *) perm_malloc(xxt_max_nnz*sizeof(REAL));
	      rvec_copy(x_ptr,x,xxt_nnz);
	      perm_free(x);
	      x = x_ptr;
	      x_ptr+=xxt_nnz;
	    }
	  xxt_nnz += len;      
	  rvec_copy(x_ptr,v+off,len);
#ifdef DEBUG	  
	  printf("XXT %d %d %d %d\n",col,i,j,len);
	  for (k=0; k<len; k++)
	    {printf("%2dXXTA %2d %2d %2d %2d %7.4f\n",dim,col,my_id,i,k,x_ptr[k]);}
#endif
	  col_indices[2*i] = off;
	  col_sz[i] = col_indices[2*i+1] = len;
	  col_vals[i] = x_ptr;
	}
      else
	{
	  col_indices[2*i] = 0;
	  col_sz[i] = col_indices[2*i+1] = 0;
	  col_vals[i] = x_ptr;
	}
    }

  /* close off stages for execution phase */
  stages[dim]=i;

  if (!my_id)
    {
      cmin = col;
#ifdef DEBUG
      printf("%2d SEPS %2d :: ",dim,my_id);
      for (col=0;col<nls;col++)
	{printf("%2d ",lseps[col]);}
      printf("\n");
#endif
      col=cmin;
    }

#ifdef DEBUG
  if (i!=m)
    {error_msg_fatal("didn't participate in all firings i=%d, m=%d, nfo=%d\n",i,m,nfo);}

  if ((j!=n)&&(!vcor))
    {error_msg_fatal("didn't fire all of mine j=%d, n=%d\n",j,n);}
#endif

  /* local nx1 vectors for fast execution phase */
  solve_uu = (REAL *) perm_malloc(m*sizeof(REAL));
  solve_w  = (REAL *) perm_malloc(m*sizeof(REAL));

  /*
  bss_free(u);
  bss_free(v);
  bss_free(uu);
  bss_free(z);
  bss_free(w);
  */

#ifdef DEBUG
  error_msg_warning("xxt_generate() :: end\n");
#endif
}


/*************************************xxt.c************************************
Function: 

Input : 
Output: 
Return: 
Description:  
**************************************xxt.c***********************************/
static
void
xxt_set_local_map(int *map, int n_loc, int m_loc, void *mylocmatvec)
{
  int i;


#ifdef DEBUG
  error_msg_warning("xxt_local_map() :: begin\n");
#endif

  a_n = n = n_loc;
  a_m = m_loc;
  a_local_to_global = (int *) perm_malloc((a_m+1)*INT_LEN);
  ivec_copy(a_local_to_global,map,a_m);
  a_local_to_global[a_m] = INT_MAX;
  A_matvec = (vfp) mylocmatvec;

#ifdef DEBUG
  for (i=0; i<a_m-1; i++)
    {
      if (a_local_to_global[i] >= a_local_to_global[i+1])
	{error_msg_warning("a_local_to_global isn't sorted?\n");}
      printf("%d :: %d, %d\n",my_id,i,a_local_to_global[i]);
    }
#endif

#ifdef DEBUG
  printf("%d :: a_map original:: ",my_id);
  for (i=0; i<a_m+1; i++)
    {printf("%d, ",a_local_to_global[i]);}
  printf("\n");
  fflush(stdout);
#endif

#ifdef DEBUG
  error_msg_warning("xxt_local_map() :: end\n");
#endif
}


/*************************************xxt.c************************************
Function: XXT_new()

Input :
Output:
Return:
Description:
**************************************xxt.c***********************************/
xxt_ADT 
XXT_new(void)
{
  xxt_ADT xxt_handle;

  n_xxt_handles++;
  check_init();
  xxt_handle=bss_malloc(sizeof(struct xxt_CDT));
  xxt_handle->id=n_xxt_handles;
  xxt_handle->info=bss_malloc(sizeof(xxt_info));
  return(xxt_handle);
};


/*************************************xxt.c************************************
Function: XXT_free()

Input :
Output:
Return:
Description:
**************************************xxt.c***********************************/
void 
XXT_free(xxt_ADT xxt_handle)
{
  if (!my_id) 
    {printf("%d :: XXT_free() :: id = %d\n",my_id,xxt_handle->id);}

  if (xxt_handle->id<=0)
    {error_msg_fatal("XXT_freee() :: bad xxt id = %d\n",xxt_handle->id);}

  n_xxt_handles--;
  bss_free(xxt_handle->info);
  bss_free(xxt_handle);
}


/*************************************xxt.c************************************
Function: XXT_generate()

Input :
Output:
Return:
Description:
**************************************xxt.c***********************************/
int 
XXT_generate(xxt_ADT xxt_handle,   /* prev. allocated xxt  handle */
	     int *local2global,    /* global column mapping       */
	     int n,                /* local num rows              */
	     int m,                /* local num cols              */
	     void *mylocmatvec,    /* b_loc=A_local.x_loc         */
	     void *grid_data       /* grid data for matvec        */
	     )
{
  /* initialize communication, malloc, etc. packages */
  check_init();

#ifdef DEBUG
  if (!my_id) 
    {printf("%d :: XXT_generate() :: id = %d\n",my_id,xxt_handle->id);}
#endif

  /* not really a robust check */
  if (xxt_handle->id<=0)
    {error_msg_fatal("ML_XXT_generate() :: bad xxt id = %d\n",xxt_handle->id);}

  /* LATER - data encapsulation */
  /* only 2^k for now */
  level = i_log2_num_nodes;
  if ((1<<level)!=num_nodes)
    {error_msg_fatal("only 2^k for now!!! %d != %d\n",1<<level,num_nodes);}

  /* set xxt communication handle to perform restricted matvec */
  gs_handle = gs_init(local2global, m, num_nodes);
  
  /* determine separators and generate firing order */
  det_separators(xxt_handle, local2global, n, m, mylocmatvec, grid_data);

  /* change name later */
  return(ML_Gen_CoarseDirect(xxt_handle,local2global, grid_data, n, m, 
			     mylocmatvec));
}


/*************************************xxt.c************************************
Function: ML_XXT_generate()

Input :
Output:
Return:
Description:
**************************************xxt.c***********************************/
int 
ML_XXT_generate(xxt_ADT xxt_handle,   /* prev. allocated xxt  handle */
		int *local2global,    /* global column mapping       */
		int n,                /* local num rows              */
		int m,                /* local num cols              */
		void *mylocmatvec,    /* b_loc=A_local.x_loc         */
		void *grid_data,      /* grid data for matvec        */
		int grid_tag,         /* grid tag for ML_Set_CSolve  */
		ML *my_ml             /* ML handle                   */
		)
{
#ifdef DEBUG
  if (!my_id) 
    {printf("%d :: ML_XXT_generate() :: id = %d\n",my_id,xxt_handle->id);}
#endif

  /* all on MPI_COMM_WORLD page? */
  if (my_ml->comm->ML_mypid!=my_id)
    {error_msg_fatal("ML_XXT_generate bad my_id %d\t%d\n",
		     my_ml->comm->ML_mypid,my_id);}
  if (my_ml->comm->ML_nprocs!=num_nodes)
    {error_msg_fatal("ML_XXT_generate bad np %d\t%d\n",
		     my_ml->comm->ML_nprocs,num_nodes);}

  /* check in solver w/ML and hold A matrix data struct ptr */
  /* ML_Set_CSolve(my_ml, grid_tag, grid_data, XXT_solve); */

  my_ml->SingleLevel[grid_tag].csolve->func->external = XXT_solve;
  my_ml->SingleLevel[grid_tag].csolve->func->ML_id = ML_EXTERNAL;
  my_ml->SingleLevel[grid_tag].csolve->data = xxt_handle;

  /* done ML specific stuff ... back to reg sched pgm */
  return(XXT_generate(xxt_handle, local2global, n, m, mylocmatvec, grid_data));
}


/*************************************xxt.c************************************
Function: XXT_solve

Input :
Output:
Return:
Description:
**************************************xxt.c***********************************/
void 

XXT_solve(xxt_ADT xxt_handle, int in_leng, double *x, int out_leng, double *b)
{
  int i, flag;

#ifdef DEBUG
  if (!my_id) 
    {printf("%d :: XXT_solve() :: id = %d\n",my_id,xxt_handle->id);}
#endif

  if (xxt_handle->id<=0)
    {error_msg_fatal("XXT_solve() :: bad xxt id = %d\n",xxt_handle->id);}

  /* change name later */
  /*
  if (!ML_Do_CoarseDirect(xxt_handle, b, x))
    {error_msg_fatal("XXT_solve() :: ML_Do failed\n");}
  */

#ifdef DEBUG
  printf("%2d bin :: ",my_id);
  for (i=0;i<n;i++)
    {printf("%f ",b[i]);}
  printf("\n");
#endif

  flag =ML_Do_CoarseDirect(xxt_handle, b, x);

#ifdef DEBUG
  printf("%2d xout :: ",my_id);
  for (i=0;i<n;i++)
    {printf("%f ",x[i]);}
  printf("\n");
#endif

  if (!flag)
    {error_msg_fatal("XXT_solve() :: ML_Do failed\n");}

}


/*************************************xxt.c************************************
Function: check_init

Input :
Output:
Return:
Description:
**************************************xxt.c***********************************/
static
void
check_init(void)
{
  /* initialize communication package */
  comm_init();

  /* initialize malloc packages - no!*/
  perm_init(); 
  bss_init();
}



/*************************************xxt.c************************************
Function: det_separators

Input :
Output:
Return:
Description:


  det_separators(xxt_handle, local2global, n, m, mylocmatvec, grid_data);

**************************************xxt.c***********************************/
static 
void 
det_separators(xxt_ADT xxt_handle,int *local2global, int n, int m,
	       void *mylocmatvec, void *grid_data)
{
  int i, ct, id;
  int mask, edge, *iptr; 
  int *dir, *used;
  int sum[4], w[4];
  int op[] = {GL_ADD,0};
  REAL *lhs, *rhs;


#ifdef DEBUG
  error_msg_warning("det_separators() :: start %d %d %d\n",level,n,m);
#endif

  dir  = bss_malloc(INT_LEN*(level+1));
  nsep = bss_malloc(INT_LEN*(level+1));
  lnsep= bss_malloc(INT_LEN*(level+1));
  fo   = bss_malloc(INT_LEN*(n+1));
  used = bss_malloc(INT_LEN*n);

  ivec_zero(dir  ,level+1);
  ivec_zero(nsep ,level+1);
  ivec_zero(lnsep,level+1);
  ivec_set (fo   ,-1,n+1);
  ivec_zero(used,n);

  lhs  = bss_malloc(REAL_LEN*m);
  rhs  = bss_malloc(REAL_LEN*m);

  /* determine separator sets top down */
  for (iptr=fo+n,id=my_id,mask=num_nodes>>1,edge=level;edge>0;edge--,mask>>=1)
    {
      /* set rsh of hc, fire, and collect lhs responses */
      (id<mask) ? rvec_zero(lhs,m) : rvec_set(lhs,1.0,m);
      gs_gop_hc(gs_handle,lhs,"+\0",edge);

      /* set lsh of hc, fire, and collect rhs responses */
      (id<mask) ? rvec_set(rhs,1.0,m) : rvec_zero(rhs,m);
      gs_gop_hc(gs_handle,rhs,"+\0",edge);

      /* count number of dofs I own that have signal and not in sep set */
      for (ivec_zero(sum,4),ct=i=0;i<n;i++)
	{
	  if (!used[i]) 
	    {
	      /* number of unmarked dofs on node */
	      ct++;
	      /* number of dofs to be marked on lhs hc */
	      if ((id< mask)&&(lhs[i]==1.0)) {sum[0]++;}
	      /* number of dofs to be marked on rhs hc */
	      if ((id>=mask)&&(rhs[i]==1.0)) {sum[1]++;}
#ifdef DEBUG
	      if ((id< mask)&&((lhs[i]!=1.0)&&(lhs[i]!=0.0)))
		{error_msg_fatal("pos bad lhs %d %d %f\n",edge,my_id,lhs[i]);}
	      if ((id>=mask)&&((rhs[i]!=1.0)&&(rhs[i]!=0.0)))
		{error_msg_fatal("pos bad rhs %d %d %f\n",edge,my_id,rhs[i]);}
#endif
	    }
	}

      /* go for load balance - choose half with most unmarked dofs, bias LHS */
      (id<mask) ? (sum[2]=ct) : (sum[3]=ct);
      giop_hc(sum,w,4,op,edge);

      /* lhs hc wins */
      if (sum[2]>=sum[3])
	{
	  if (id<mask)
	    {
	      /* mark dofs I own that have signal and not in sep set */
	      for (ct=i=0;i<n;i++)
		{
		  if ((!used[i])&&(lhs[i]==1.0))
		    {
		      ct++; nfo++;
#ifdef DEBUG
		      if (nfo>n)
			{error_msg_fatal("nfo about to exceed n\n");}
#endif
		      *--iptr = local2global[i];
		      used[i]=edge;
		    }
		}
	      if (ct>1) {ivec_sort(iptr,ct);}
	      lnsep[edge]=ct;
	    }
	  nsep[edge]=sum[0];
	  dir [edge]=LEFT;
	}
      /* rhs hc wins */
      else
	{
	  if (id>=mask)
	    {
	      /* mark dofs I own that have signal and not in sep set */
	      for (ct=i=0;i<n;i++)
		{
		  if ((!used[i])&&(rhs[i]==1.0))
		    {
		      ct++; nfo++;
#ifdef DEBUG
		      if (nfo>n)
			{error_msg_fatal("nfo about to exceed n\n");}
#endif
		      *--iptr = local2global[i];
		      used[i]=edge;
		    }
		}
	      if (ct>1) {ivec_sort(iptr,ct);}
	      lnsep[edge]=ct;
	    }
	  nsep[edge]=sum[1];
	  dir [edge]=RIGHT;
	}
      /* LATER or we can recur on these to order seps at this level */
      /* do we need full set of separators for this?                */

      /* fold rhs hc into lower */
      if (id>=mask)
        {id-=mask;}
    }

  /* level 0 is on processor case - so mark the remainder */
  for (ct=i=0;i<n;i++)
    {
      if (!used[i]) 
	{
	  ct++; nfo++;
#ifdef DEBUG
	  if (nfo>n)
	    {error_msg_fatal("nfo about to exceed n\n");}
#endif
	  *--iptr = local2global[i];
	  used[i]=edge;
	}
    }
  if (ct>1) {ivec_sort(iptr,ct);}
  lnsep[edge]=ct;
  nsep [edge]=ct;
  dir  [edge]=LEFT;

  /* LATER must free ? 
  bss_free(dir);
  bss_free(nsep);
  bss_free(lnsep);
  bss_free(lhs);
  bss_free(rhs);
  bss_free(fo);
  bss_free(used);
  */

#ifdef DEBUG  
  if (nfo!=n)
    {error_msg_fatal("hey nfo!= n?\n");}

  if (iptr!=fo)
    {error_msg_fatal("iptr and fo misaligned?\n");}

  /* dump out info */
  for (edge=0; edge<=level; edge++)
    {
      /* do I own it? I should */
      /* LATER expand search to columns of a ==> no gs_gop_hc before mat-vec */
      printf("%2d %2d SEPG %2d %2d :: ",edge,my_id,dir[edge],nsep[edge]);
      for (i=0;i<lnsep[edge];i++)
	{printf("%2d ",*iptr++);}
      printf("\n");
    }
#endif

#ifdef DEBUG  
  error_msg_warning("det_separators() :: end\n");
#endif
}

void
dump_vec(REAL *v, int n, int tag, int tag2, char * s)
{
  int i;
  printf("%2d %2d %s %2d :: ",tag,tag2,s,my_id);
  for (i=0;i<n;i++)
    {printf("%f ",v[i]);}
  printf("\n");
  fflush(stdout);
}

void
dump_ivec(int *v, int n, int tag, int tag2, char * s)
{
  int i;
  printf("%2d %2d %s %2d :: ",tag,tag2,s,my_id);
  for (i=0;i<n;i++)
    {printf("%2d ",v[i]);}
  printf("\n");
  fflush(stdout);
}
