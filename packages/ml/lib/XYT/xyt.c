/*************************************xyt.c************************************
Module Name: xyt
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
**************************************xyt.c***********************************/


/*************************************xyt.c************************************
NOTES ON USAGE: 

**************************************xyt.c***********************************/


#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>

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
#ifdef MLSRC
#include "ml_include.h"
#endif
#include "xyt.h"

#ifdef DEBUG
#define MAX_SEPARATORS 4096
#endif

#define LEFT  -1
#define RIGHT  1
#define BOTH   0
#define MAX_FORTRAN_HANDLES  10

typedef struct xyt_solver_info {
  int n, m, n_global, m_global;
  int nnz, max_nnz, msg_buf_sz;
  int *nsep, *lnsep, *fo, nfo, *stages;
  int *xcol_sz, *xcol_indices; 
  REAL **xcol_vals, *x, *solve_uu, *solve_w;
  int *ycol_sz, *ycol_indices; 
  REAL **ycol_vals, *y;
} xyt_info;

typedef struct matvec_info {
  int n, m, n_global, m_global;
  int *local2global;
  gs_ADT gs_handle;
  vfp matvec;
  void *grid_data;
} mv_info;

struct xyt_CDT{
  int id;
  int ns;
  int level;
  xyt_info *info;
  mv_info  *mvi;
};

static int n_xyt=0;
static int n_xyt_handles=0;
static int fn_xyt_handles=0;
static xyt_ADT fhandles[MAX_FORTRAN_HANDLES+1];

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
static void do_xyt_solve(xyt_ADT xyt_handle, REAL *rhs);
static int do_xyt_factor(xyt_ADT xyt_handle);
static int xyt_generate(xyt_ADT xyt_handle);
static void check_init(void);
static void check_handle(xyt_ADT xyt_handle);
static void det_separators(xyt_ADT xyt_handle);
static mv_info *set_mvi(int *local2global, int n, int m, void *matvec, void *grid_data);
static void do_matvec(mv_info *A, REAL *v, REAL *u);
#ifdef MLSRC
void ML_XYT_solve(xyt_ADT xyt_handle, int lx, double *x, int lb, double *b);
int  ML_XYT_factor(xyt_ADT xyt_handle, int *local2global, int n, int m,
		   void *matvec, void *grid_data, int grid_tag, ML *my_ml);
#endif

#ifdef NOT
/*************************************xyt.c************************************
Function: 

Input : 
Output: 
Return: 
Description:  
**************************************xyt.c***********************************/
#if defined UPCASE
void
XYT_MATVEC (void *matrix_data, double *in, double *out)
#else
xyt_matvec_(void *matrix_data, double *in, double *out)
#endif
{
;
}
#if defined UPCASE
void XYT_MATVEC(void *matrix_data, double *in, double *out);
else
void xyt_matvec_(void *matrix_data, double *in, double *out);
#endif

#endif




/*************************************xyt.c************************************
Function: do_xyt_factor

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
**************************************xyt.c***********************************/
static
int
do_xyt_factor(xyt_ADT xyt_handle)
{
  int flag;


#ifdef DEBUG
  error_msg_warning("do_xyt_factor() :: begin\n");
#endif

  flag=xyt_generate(xyt_handle);

#ifdef INFO
  xyt_stats(xyt_handle);
  bss_stats(); 
  perm_stats(); 
#endif

#ifdef DEBUG
  error_msg_warning("do_xyt_factor() :: end\n");
#endif

  return(flag);
}



/*************************************xyt.c************************************
Function: 

Input : 
Output: 
Return: 
Description:  
**************************************xyt.c***********************************/
static
int
xyt_generate(xyt_ADT xyt_handle)
{
  int i,j,k,l,index;
  int dim, col;
  REAL *u, *uu, *v, *z, *w, alpha, alpha_w;
  int *col_map, *segs;
  int op[] = {GL_ADD,0};
  int max=0;
  int off, len;
  REAL *x_ptr, *y_ptr;
  int *iptr, flag;
  int start=0, end, work;
  int op2[] = {GL_MIN,0};
  int op3[] = {GL_MAX,0};
  int cmin;
  int id, mask;
  gs_ADT gs_handle;
  int *nsep, *lnsep, *fo, nfo;
  int a_n=xyt_handle->mvi->n;
  int a_m=xyt_handle->mvi->m;
  int *a_local2global=xyt_handle->mvi->local2global;
  int level;
  int n, m;
  int *xcol_sz, *xcol_indices, *stages; 
  REAL **xcol_vals, *x;
  int *ycol_sz, *ycol_indices;
  REAL **ycol_vals, *y;
  int n_global;
  int xt_nnz=0, xt_max_nnz=0;
  int yt_nnz=0, yt_max_nnz=0;
  int xt_zero_nnz  =0;
  int xt_zero_nnz_0=0;
  int yt_zero_nnz  =0;
  int yt_zero_nnz_0=0;
#ifdef TIMINGS
  int op4[] = {NON_UNIFORM,GL_MIN,GL_MAX,GL_ADD};
#if   defined NXSRC
  double update_time=0.0, xt_time=0.0, yt_time=0.0, comm_time=0.0, mxm_time=0.0;
  double dclock(), time, vals[3], dwork[3];
#elif defined MPISRC
  double update_time=0.0, xt_time=0.0, yt_time=0.0, comm_time=0.0, mxm_time=0.0;
  double MPI_Wtime(), time, vals[3], dwork[3];
#endif
#endif


#ifdef DEBUG
  int lseps[MAX_SEPARATORS];
  int nls=0;
#endif

  
#ifdef DEBUG
  error_msg_warning("xyt_generate() :: begin\n");
#endif

  n=xyt_handle->mvi->n; 
  nsep=xyt_handle->info->nsep; 
  lnsep=xyt_handle->info->lnsep;
  fo=xyt_handle->info->fo;
  nfo=xyt_handle->info->nfo;
  end=lnsep[0];
  level=xyt_handle->level;
  gs_handle=xyt_handle->mvi->gs_handle;

  /* is there a null space? */
  /* LATER add in ability to detect null space by checking alpha */
  for (i=0, j=0; i<=level; i++)
    {j+=nsep[i];}

  m = j-xyt_handle->ns;
  if (m!=j)
    {printf("xyt_generate() :: null space exists %d %d %d\n",m,j,xyt_handle->ns);}

  error_msg_warning("xyt_generate() :: X(%d,%d)\n",n,m);    

  /* get and initialize storage for x local         */
  /* note that x local is nxm and stored by columns */
  xcol_sz = (int *) bss_malloc(m*INT_LEN);
  xcol_indices = (int *) bss_malloc((2*m+1)*sizeof(int));
  xcol_vals = (REAL **) bss_malloc(m*sizeof(REAL *));
  for (i=j=0; i<m; i++, j+=2)
    {
      xcol_indices[j]=xcol_indices[j+1]=xcol_sz[i]=-1;
      xcol_vals[i] = NULL;
    }
  xcol_indices[j]=-1;

  /* get and initialize storage for y local         */
  /* note that y local is nxm and stored by columns */
  ycol_sz = (int *) bss_malloc(m*INT_LEN);
  ycol_indices = (int *) bss_malloc((2*m+1)*sizeof(int));
  ycol_vals = (REAL **) bss_malloc(m*sizeof(REAL *));
  for (i=j=0; i<m; i++, j+=2)
    {
      ycol_indices[j]=ycol_indices[j+1]=ycol_sz[i]=-1;
      ycol_vals[i] = NULL;
    }
  ycol_indices[j]=-1;

  /* size of separators for each sub-hc working from bottom of tree to top */
  /* this looks like nsep[]=segments */
  stages = (int *) bss_malloc((level+1)*INT_LEN);
  segs   = (int *) bss_malloc((level+1)*INT_LEN);
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
  n_global = xyt_handle->info->n_global;
  xt_max_nnz = yt_max_nnz = (12*pow(1.0*n_global,1.6667) + j*n/2)/num_nodes;
  x = (REAL *) bss_malloc((xt_max_nnz+yt_max_nnz)*sizeof(REAL));
  y = x + xt_max_nnz;

  /* LATER - can embed next sep to fire in gs */
  /* time to make the donuts - generate X factor */
  for (id=dim=i=j=0;i<m;i++)
    {
      /* time to move to the next level? */
      while (i==segs[dim])
	{
	  stages[dim++]=i;
	  end+=lnsep[dim];
	}
      stages[dim]=i;

      /* which column are we firing? */
      /* i.e. set v_l */
      /* use new seps and do global min across hc to determine which one to fire */
      (start<end) ? (col=fo[start]) : (col=INT_MAX);
      giop_hc(&col,&work,1,op2,dim); 

      /* do I own it? I should */
      rvec_zero(v ,a_m);
      if (col==fo[start])
	{
	  start++;
	  index=ivec_linear_search(col, a_local2global, a_n);
	  v[index] = 1.0; 
	  j++;
	}
      else
	{
	  index=ivec_linear_search(col, a_local2global, a_m);
	  if (index!=-1)
	    {
	      v[index] = 1.0; 
	    }
	}

      /* perform u = A.v_l */
#ifdef TIMINGS
#if   defined NXSRC
      time = dclock();
#elif defined MPISRC
      time = MPI_Wtime();
#endif
#endif
      rvec_zero(u,n);
      do_matvec(xyt_handle->mvi,v,u);
#ifdef TIMINGS
#if   defined NXSRC
      mxm_time += dclock() - time;
#elif defined MPISRC
      mxm_time += MPI_Wtime() - time;
#endif
#endif

      /* uu =  X^T.u_l (local portion) */
      /* technically only need to zero out first i entries */
      /* later turn this into an XYT_solve call ? */
#ifdef TIMINGS
#if   defined NXSRC
      time = dclock();
#elif defined MPISRC
      time = MPI_Wtime();
#endif
#endif
      rvec_zero(uu,m);
      y_ptr=y;
      iptr = ycol_indices;
      for (k=0; k<i; k++)
	{
	  off = *iptr++;
	  len = *iptr++;

#if   BLAS&&r8
	  uu[k] = cblas_ddot(len,u+off,1,y_ptr,1);
#elif BLAS
	  uu[k] = cblas_sdot(len,u+off,1,y_ptr,1);
#else
	  uu[k] = rvec_dot(u+off,y_ptr,len);
#endif
	  y_ptr+=len;
	}
#ifdef TIMINGS
#if   defined NXSRC
      yt_time += dclock() - time;
#elif defined MPISRC
      yt_time += MPI_Wtime() - time;
#endif
#endif
      /* uu = X^T.u_l (comm portion) */
      /* not needed! rvec_zero(w,m); */
#ifdef TIMINGS
#if   defined NXSRC
      time = dclock();
#elif defined MPISRC
      time = MPI_Wtime();
#endif
#endif
      ssgl_radd  (uu, w, dim, stages);
#ifdef TIMINGS
#if   defined NXSRC
      comm_time += dclock() - time;
#elif defined MPISRC
      comm_time += MPI_Wtime() - time;
#endif
#endif

      /* z = X.uu */
#ifdef TIMINGS
#if   defined NXSRC
      time = dclock();
#elif defined MPISRC
      time = MPI_Wtime();
#endif
#endif
      rvec_zero(z,n);
      x_ptr=x;
      iptr = xcol_indices;
      for (k=0; k<i; k++)
	{
	  off = *iptr++;
	  len = *iptr++;

#if   BLAS&r8
	  cblas_daxpy(len,uu[k],x_ptr,1,z+off,1);
#elif BLAS
	  cblas_saxpy(len,uu[k],x_ptr,1,z+off,1);
#else
	  rvec_axpy(z+off,x_ptr,uu[k],len);
#endif
	  x_ptr+=len;
	}
#ifdef TIMINGS
#if   defined NXSRC
      xt_time += dclock() - time;
#elif defined MPISRC
      xt_time += MPI_Wtime() - time;
#endif
#endif

      /* compute v_l = v_l - z */
#ifdef TIMINGS
#if   defined NXSRC
      time = dclock();
#elif defined MPISRC
      time = MPI_Wtime();
#endif
#endif
      rvec_zero(v+a_n,a_m-a_n);
#if   BLAS&&r8
      cblas_daxpy(n,-1.0,z,1,v,1);
#elif BLAS
      cblas_saxpy(n,-1.0,z,1,v,1);
#else
      rvec_axpy(v,z,-1.0,n);
#endif

      /* compute u_l = A.v_l */
      if (a_n!=a_m)
	{gs_gop_hc(gs_handle,v,"+\0",dim);}
      rvec_zero(u,n);
      do_matvec(xyt_handle->mvi,v,u);
#ifdef TIMINGS
#if   defined NXSRC
      mxm_time += dclock() - time;
#elif defined MPISRC
      mxm_time += MPI_Wtime() - time;
#endif
#endif

      /* compute sqrt(alpha) = sqrt(u_l^T.u_l) - local portion */
#ifdef TIMINGS
#if   defined NXSRC
      time = dclock();
#elif defined MPISRC
      time = MPI_Wtime();
#endif
#endif
#if   BLAS&&r8
      alpha = cblas_ddot(n,u,1,u,1);
#elif BLAS
      alpha = cblas_sdot(n,u,1,u,1);
#else
      alpha = rvec_dot(u,u,n);
#endif

      /* compute sqrt(alpha) = sqrt(u_l^T.u_l) - comm portion */
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

      /* compute v_l = v_l/sqrt(alpha) */
      rvec_scale(v,1.0/alpha,n);
      rvec_scale(u,1.0/alpha,n);

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

      len -= (off-1);

      if (len>0)
	{
	  if ((xt_nnz+len)>xt_max_nnz)
	    {
	      xt_max_nnz *= 2;
	      x_ptr = (REAL *) bss_malloc(xt_max_nnz*sizeof(REAL));
	      rvec_copy(x_ptr,x,xt_nnz);
	      bss_free(x);
	      x = x_ptr;
	      x_ptr+=xt_nnz;
	    }
	  xt_nnz += len;      
	  rvec_copy(x_ptr,v+off,len);

          /* keep track of number of zeros */

#ifdef INFO
	  if (dim)
	    {
	      for (k=0; k<len; k++)
		{
		  if (x_ptr[k]==0.0)
		    {xt_zero_nnz++;}
		}
	    }
	  else
	    {
	      for (k=0; k<len; k++)
		{
		  if (x_ptr[k]==0.0)
		    {xt_zero_nnz_0++;}
		}
	    }
#endif

	  xcol_indices[2*i] = off;
	  xcol_sz[i] = xcol_indices[2*i+1] = len;
	  xcol_vals[i] = x_ptr;
	}
      else
	{
	  xcol_indices[2*i] = 0;
	  xcol_sz[i] = xcol_indices[2*i+1] = 0;
	  xcol_vals[i] = x_ptr;
	}


      /* add newly generated column, u_l, to Y */
      flag = 1;
      off=len=0;
      for (k=0; k<n; k++)
	{
	  if (u[k]!=0.0)
	    {
	      len=k;
	      if (flag)
		{off=k; flag=0;}
	    }
	}

      len -= (off-1);

      if (len>0)
	{
	  if ((yt_nnz+len)>yt_max_nnz)
	    {
	      yt_max_nnz *= 2;
	      y_ptr = (REAL *) bss_malloc(yt_max_nnz*sizeof(REAL));
	      rvec_copy(y_ptr,y,yt_nnz);
	      bss_free(y);
	      y = y_ptr;
	      y_ptr+=yt_nnz;
	    }
	  yt_nnz += len;      
	  rvec_copy(y_ptr,u+off,len);

#ifdef INFO
          /* keep track of number of zeros */
	  if (dim)
	    {
	      for (k=0; k<len; k++)
		{
		  if (y_ptr[k]==0.0)
		    {yt_zero_nnz++;}
		}
	    }
	  else
	    {
	      for (k=0; k<len; k++)
		{
		  if (y_ptr[k]==0.0)
		    {yt_zero_nnz_0++;}
		}
	    }
#endif
	  ycol_indices[2*i] = off;
	  ycol_sz[i] = ycol_indices[2*i+1] = len;
	  ycol_vals[i] = y_ptr;
	}
      else
	{
	  ycol_indices[2*i] = 0;
	  ycol_sz[i] = ycol_indices[2*i+1] = 0;
	  ycol_vals[i] = y_ptr;
	}
#ifdef TIMINGS
#if   defined NXSRC
      update_time += dclock() - time;
#elif defined MPISRC
      update_time += MPI_Wtime() - time;
#endif
#endif
    }

  /* close off stages for execution phase */
  while (dim!=level)
    {
      stages[dim++]=i;
      error_msg_warning("disconnected!!! dim(%d)!=level(%d)\n",dim,level);
    }
  stages[dim]=i;

  cmin=xt_zero_nnz+yt_zero_nnz;;
  giop(&cmin,&work,1,op);

  col=xt_zero_nnz_0 + yt_zero_nnz_0;
  giop(&col,&work,1,op);

  alpha = (double) (xt_nnz + yt_nnz);
  grop(&alpha,&alpha_w,1,op);
  if (!my_id)
    {printf("%d:xyt nnz 0's: %d %d %d (%f)\n",my_id,cmin,col,cmin+col,
	                 (1.0*(cmin+col))/alpha);
    }

  xyt_handle->info->n=xyt_handle->mvi->n;
  xyt_handle->info->m=m;
  xyt_handle->info->nnz=xt_nnz + yt_nnz;
  xyt_handle->info->max_nnz=xt_max_nnz + yt_max_nnz;
  xyt_handle->info->msg_buf_sz=stages[level]-stages[0];
  xyt_handle->info->solve_uu = (REAL *) bss_malloc(m*sizeof(REAL));
  xyt_handle->info->solve_w  = (REAL *) bss_malloc(m*sizeof(REAL));
  xyt_handle->info->x=x;
  xyt_handle->info->xcol_vals=xcol_vals;
  xyt_handle->info->xcol_sz=xcol_sz;
  xyt_handle->info->xcol_indices=xcol_indices;  
  xyt_handle->info->stages=stages;
  xyt_handle->info->y=y;
  xyt_handle->info->ycol_vals=ycol_vals;
  xyt_handle->info->ycol_sz=ycol_sz;
  xyt_handle->info->ycol_indices=ycol_indices;  

  bss_free(segs);
  bss_free(u);
  bss_free(v);
  bss_free(uu);
  bss_free(z);
  bss_free(w);

#ifdef TIMINGS
  vals[0]=vals[1]=vals[2]=yt_time;
  grop(vals,dwork,sizeof(op4)/sizeof(op4[0])-1,op4);

  if (!my_id)
    {
      printf("%d :: min   xyt_ytt=%g\n",my_id,vals[0]);
      printf("%d :: max   xyt_ytt=%g\n",my_id,vals[1]);
      printf("%d :: avg   xyt_ytt=%g\n",my_id,vals[2]/num_nodes);
    }

  vals[0]=vals[1]=vals[2]=xt_time;
  grop(vals,dwork,sizeof(op4)/sizeof(op4[0])-1,op4);

  if (!my_id)
    {
      printf("%d :: min   xyt_xtt=%g\n",my_id,vals[0]);
      printf("%d :: max   xyt_xtt=%g\n",my_id,vals[1]);
      printf("%d :: avg   xyt_xtt=%g\n",my_id,vals[2]/num_nodes);
    }

  vals[0]=vals[1]=vals[2]=comm_time;
  grop(vals,dwork,sizeof(op4)/sizeof(op4[0])-1,op4);

  if (!my_id)
    {
      printf("%d :: min   xyt_com=%g\n",my_id,vals[0]);
      printf("%d :: max   xyt_com=%g\n",my_id,vals[1]);
      printf("%d :: avg   xyt_com=%g\n",my_id,vals[2]/num_nodes);
    }

  vals[0]=vals[1]=vals[2]=mxm_time;
  grop(vals,dwork,sizeof(op4)/sizeof(op4[0])-1,op4);

  if (!my_id)
    {
      printf("%d :: min   xyt_mxm=%g\n",my_id,vals[0]);
      printf("%d :: max   xyt_mxm=%g\n",my_id,vals[1]);
      printf("%d :: avg   xyt_mxm=%g\n",my_id,vals[2]/num_nodes);
    }

  vals[0]=vals[1]=vals[2]=update_time;
  grop(vals,dwork,sizeof(op4)/sizeof(op4[0])-1,op4);

  if (!my_id)
    {
      printf("%d :: min   xyt_upd=%g\n",my_id,vals[0]);
      printf("%d :: max   xyt_upd=%g\n",my_id,vals[1]);
      printf("%d :: avg   xyt_upd=%g\n",my_id,vals[2]/num_nodes);
    }

#endif

#ifdef DEBUG
  error_msg_warning("xyt_generate() :: end\n");
#endif

  return(TRUE);
}



/*************************************xyt.c************************************
Function: 

Input : 
Output: 
Return: 
Description:  
**************************************xyt.c***********************************/
void
XYT_stats(xyt_ADT xyt_handle)
{
  int edge, i, *iptr, *nsep;
  int vals[9], work[9], op[] = {NON_UNIFORM,GL_MIN,GL_MAX,GL_ADD,GL_MIN,GL_MAX,GL_ADD,GL_MIN,GL_MAX,GL_ADD};


#ifdef DEBUG
  error_msg_warning("xyt_stats() :: begin\n");
#endif

  vals[0]=vals[1]=vals[2]=xyt_handle->info->nnz;
  vals[3]=vals[4]=vals[5]=xyt_handle->mvi->n;
  vals[6]=vals[7]=vals[8]=xyt_handle->info->msg_buf_sz;
  giop(vals,work,sizeof(op)/sizeof(op[0])-1,op);

  /* assume that it's square and the dofs aren't shared */
  /*
  xyt_handle->info->n_global=xyt_handle->info->m_global=vals[5];
  xyt_handle->mvi->n_global=xyt_handle->mvi->m_global=vals[5];
  */

  if (!my_id) 
    {
      printf("%d :: min   xyt_nnz=%d\n",my_id,vals[0]);
      printf("%d :: max   xyt_nnz=%d\n",my_id,vals[1]);
      printf("%d :: avg   xyt_nnz=%g\n",my_id,1.0*vals[2]/num_nodes);
      printf("%d :: tot   xyt_nnz=%d\n",my_id,vals[2]);
      printf("%d :: xyt   C(2d)  =%g\n",my_id,vals[2]/(pow(1.0*vals[5],1.5)));
      printf("%d :: xyt   C(3d)  =%g\n",my_id,vals[2]/(pow(1.0*vals[5],1.6667)));
      printf("%d :: min   xyt_n  =%d\n",my_id,vals[3]);
      printf("%d :: max   xyt_n  =%d\n",my_id,vals[4]);
      printf("%d :: avg   xyt_n  =%g\n",my_id,1.0*vals[5]/num_nodes);
      printf("%d :: tot   xyt_n  =%d\n",my_id,vals[5]);
      printf("%d :: min   xyt_buf=%d\n",my_id,vals[6]);
      printf("%d :: max   xyt_buf=%d\n",my_id,vals[7]);
      printf("%d :: avg   xyt_buf=%g\n",my_id,1.0*vals[8]/num_nodes);
    }

  iptr = xyt_handle->info->fo;
  nsep = xyt_handle->info->nsep;
  printf("%2d NSEPG  :: ",my_id);
  for (edge=0; edge<=xyt_handle->level; edge++)
    {
      printf("%2d ",nsep[edge]);
    }
  printf("\n");


#ifdef DEBUG
  error_msg_warning("xyt_stats() :: end\n");
#endif
}



/*************************************xyt.c************************************
Function: 

Input : 
Output: 
Return: 
Description:  
**************************************xyt.c***********************************/
static
void
do_xyt_solve(xyt_ADT xyt_handle, register REAL *uc)
{
  register int off, len, *iptr;
  int level       =xyt_handle->level;
  int n           =xyt_handle->info->n;
  int m           =xyt_handle->info->m;
  int *stages     =xyt_handle->info->stages;
  int *xcol_indices=xyt_handle->info->xcol_indices;
  int *ycol_indices=xyt_handle->info->ycol_indices;
  register REAL *x_ptr, *y_ptr, *uu_ptr;
  REAL zero=0.0;
  REAL *solve_uu=xyt_handle->info->solve_uu;
  REAL *solve_w =xyt_handle->info->solve_w;
  REAL *x       =xyt_handle->info->x;
  REAL *y       =xyt_handle->info->y;

#ifdef DEBUG
  error_msg_warning("do_xyt_solve() :: begin\n");
#endif

  uu_ptr=solve_uu;
#if   BLAS&&r8
  cblas_dcopy(m,&zero,0,uu_ptr,1);
#elif BLAS
  cblas_scopy(m,&zero,0,uu_ptr,1);
#else
  rvec_zero(uu_ptr,m);
#endif

  /* x  = X.Y^T.b */
  /* uu = Y^T.b */
  for (y_ptr=y,iptr=ycol_indices; *iptr!=-1; y_ptr+=len)
    {
      off=*iptr++; len=*iptr++;
#if   BLAS&&r8
      *uu_ptr++ = cblas_ddot(len,uc+off,1,y_ptr,1);
#elif BLAS
      *uu_ptr++ = cblas_sdot(len,uc+off,1,y_ptr,1);
#else
      *uu_ptr++ = rvec_dot(uc+off,y_ptr,len);
#endif
    }

  /* comunication of beta */
  uu_ptr=solve_uu;
  if (level) {ssgl_radd(uu_ptr, solve_w, level, stages);}

#if   BLAS&&r8
  cblas_dcopy(n,&zero,0,uc,1);
#elif BLAS
  cblas_scopy(n,&zero,0,uc,1);
#else
  rvec_zero(uc,n);
#endif

  /* x = X.uu */
  for (x_ptr=x,iptr=xcol_indices; *iptr!=-1; x_ptr+=len)
    {
      off=*iptr++; len=*iptr++;
#if   BLAS&&r8
      cblas_daxpy(len,*uu_ptr++,x_ptr,1,uc+off,1);
#elif BLAS
      cblas_saxpy(len,*uu_ptr++,x_ptr,1,uc+off,1);
#else
      rvec_axpy(uc+off,x_ptr,*uu_ptr++,len);
#endif
    }

#ifdef DEBUG
  error_msg_warning("do_xyt_solve() :: end\n");
#endif
}



/*************************************xyt.c************************************
Function: xyt_new_()

Input :
Output:
Return:
Description:
**************************************xyt.c***********************************/
#if defined UPCASE
int 
XYT_NEW (void)
#else
int 
xyt_new_(void)
#endif
{
  int i;
  xyt_ADT xyt_handle;


  if (fn_xyt_handles==MAX_FORTRAN_HANDLES)
    {error_msg_fatal("xyt_new_() :: too many xyt handles %d\n",MAX_FORTRAN_HANDLES);}

  fn_xyt_handles++;
  xyt_handle = XYT_new();

  for (i=1;i<MAX_FORTRAN_HANDLES;i++)
    {
      if (!fhandles[i])
	{fhandles[i]=xyt_handle; return(i);}
    }
  return(-1);
}



/*************************************xyt.c************************************
Function: XYT_new()

Input :
Output:
Return:
Description:
**************************************xyt.c***********************************/
xyt_ADT 
XYT_new(void)
{
  xyt_ADT xyt_handle;


#ifdef DEBUG
  error_msg_warning("XYT_new() :: start %d\n",n_xyt_handles);
#endif

  n_xyt_handles++;
  xyt_handle       = bss_malloc(sizeof(struct xyt_CDT));
  xyt_handle->id   = ++n_xyt;
  xyt_handle->info = NULL;
  xyt_handle->mvi  = NULL;

#ifdef DEBUG
  error_msg_warning("XYT_new() :: end   %d\n",n_xyt_handles);
#endif

  return(xyt_handle);
}



/*************************************xyt.c************************************
Function: XYT_factor()

Input :
Output:
Return:
Description:
**************************************xyt.c***********************************/
#if defined UPCASE
int 
XYT_FACTOR (int *ixyt_handle,   /* prev. allocated xyt  handle */
	    int *local2global, /* global column mapping       */
	    int *n,           /* local num rows              */
	    int *m,           /* local num cols              */
	    void *matvec,      /* b_loc=A_local.x_loc         */
	    void *grid_data    /* grid data for matvec        */
	    )
#else
int 
xyt_factor_(int *ixyt_handle,   /* prev. allocated xyt  handle */
	    int *local2global, /* global column mapping       */
	    int *n,            /* local num rows              */
	    int *m,            /* local num cols              */
	    void *matvec,      /* b_loc=A_local.x_loc         */
	    void *grid_data    /* grid data for matvec        */
	    )
#endif
{

#ifdef NOT
#if defined UPCASE
  if (!matvec)
    {matvec=(void *)XYT_MATVEC;}
#else
  if (!matvec)
    {matvec=(void *)xyt_matvec_;}
#endif
#endif

  return(XYT_factor(fhandles[*ixyt_handle],local2global,*n,*m,matvec,grid_data));
}
	 




/*************************************xyt.c************************************
Function: XYT_factor()

Input :
Output:
Return:
Description:
**************************************xyt.c***********************************/
int 
XYT_factor(xyt_ADT xyt_handle, /* prev. allocated xyt  handle */
	   int *local2global,  /* global column mapping       */
	   int n,              /* local num rows              */
	   int m,              /* local num cols              */
	   void *matvec,       /* b_loc=A_local.x_loc         */
	   void *grid_data     /* grid data for matvec        */
	   )
{
  int flag;
#ifdef TIMINGS
  int op[] = {NON_UNIFORM,GL_MIN,GL_MAX,GL_ADD};
#if   defined NXSRC
  double dclock(), time, vals[3], work[3];
#elif defined MPISRC
  double MPI_Wtime(), time, vals[3], work[3];
#endif
#endif
	

#ifdef DEBUG
  error_msg_warning("XYT_factor() :: start %d\n",n_xyt_handles);
#endif

  check_init();

#ifdef SAFE
  check_handle(xyt_handle);
#endif

  /* only 2^k for now and all nodes participating */
  if ((1<<(xyt_handle->level=i_log2_num_nodes))!=num_nodes)
    {error_msg_fatal("only 2^k for now and MPI_COMM_WORLD!!! %d != %d\n",
		     1<<i_log2_num_nodes,num_nodes);}

  /* space for X info */
  xyt_handle->info = bss_malloc(sizeof(xyt_info));

  /* set up matvec handles */
  xyt_handle->mvi  = set_mvi(local2global, n, m, matvec, grid_data);

  /* matrix is assumed to be of full rank */
  xyt_handle->ns=0;

  /* determine separators and generate firing order - NB xyt info set here */
  /* i.e. generate structural info */
#ifdef TIMINGS
#if   defined NXSRC
  time = dclock();
#elif defined MPISRC
  time = MPI_Wtime();
#endif
#endif

  det_separators(xyt_handle);

#ifdef TIMINGS
#if   defined NXSRC
  time = dclock() - time;
#elif defined MPISRC
  time = MPI_Wtime() - time;
#endif

  vals[0]=vals[1]=vals[2]=time;
  grop(vals,work,sizeof(op)/sizeof(op[0])-1,op);
  if (!my_id)
    {
      printf("%d :: min   xyt_str=%g\n",my_id,vals[0]);
      printf("%d :: max   xyt_str=%g\n",my_id,vals[1]);
      printf("%d :: avg   xyt_str=%g\n",my_id,vals[2]/num_nodes);
    }
#endif

  /* do the numerical factorization */
#ifdef TIMINGS
#if   defined NXSRC
  time = dclock();
#elif defined MPISRC
  time = MPI_Wtime();
#endif
#endif

  flag = do_xyt_factor(xyt_handle);

#ifdef TIMINGS
#if   defined NXSRC
  time = dclock() - time;
#elif defined MPISRC
  time = MPI_Wtime() - time;
#endif

  vals[0]=vals[1]=vals[2]=time;
  grop(vals,work,sizeof(op)/sizeof(op[0])-1,op);

  if (!my_id)
    {
      printf("%d :: min   xyt_num=%g\n",my_id,vals[0]);
      printf("%d :: max   xyt_num=%g\n",my_id,vals[1]);
      printf("%d :: avg   xyt_num=%g\n",my_id,vals[2]/num_nodes);
    }
#endif

#ifdef DEBUG
  error_msg_warning("XYT_factor() :: end   %d (flag=%d)\n",n_xyt_handles,flag);
#endif

  return(flag);
}




/*************************************xyt.c************************************
Function: XYT_solve

Input :
Output:
Return:
Description:
**************************************xyt.c***********************************/
void 
new_XYT_solve(xyt_ADT xyt_handle, double *b)
{
#ifdef TIMINGS
  int op[] = {NON_UNIFORM,GL_MIN,GL_MAX,GL_ADD};
#if   defined NXSRC
  double dclock(), time, vals[3], work[3];
#elif defined MPISRC
  double MPI_Wtime(), time, vals[3], work[3];
#endif
#endif


#ifdef DEBUG
  error_msg_warning("new_XYT_solve() :: start %d\n",n_xyt_handles);
#endif

#ifdef SAFE
  check_init();
  check_handle(xyt_handle);
#endif

#ifdef TIMINGS
#if   defined NXSRC
  time = dclock();
#elif defined MPISRC
  time = MPI_Wtime();
#endif
#endif

  do_xyt_solve(xyt_handle,b);

#ifdef TIMINGS
#if   defined NXSRC
  time = dclock() - time;
#elif defined MPISRC
  time = MPI_Wtime() - time;
#endif

  vals[0]=vals[1]=vals[2]=time;
  grop(vals,work,sizeof(op)/sizeof(op[0])-1,op);

  if (!my_id)
    {
      printf("%d :: min   xyt_slv=%g\n",my_id,vals[0]);
      printf("%d :: max   xyt_slv=%g\n",my_id,vals[1]);
      printf("%d :: avg   xyt_slv=%g\n",my_id,vals[2]/num_nodes);
    }
#endif

#ifdef DEBUG
  error_msg_warning("new_XYT_solve() :: end   %d\n",n_xyt_handles);
#endif
}



/*************************************xyt.c************************************
Function: xyt_solve_()

Input :
Output:
Return:
Description:
**************************************xyt.c***********************************/
#if defined UPCASE
void
XYT_SOLVE (int *n, double *b)
#else
void
new_xyt_solve_(int *n, double *b)
#endif
{
  new_XYT_solve(fhandles[*n],b);
}



/*************************************xyt.c************************************
Function: XYT_solve

Input :
Output:
Return:
Description:
**************************************xyt.c***********************************/
void 
XYT_solve(xyt_ADT xyt_handle, double *x, double *b)
{
#ifdef TIMINGS
  int op[] = {NON_UNIFORM,GL_MIN,GL_MAX,GL_ADD};
#if   defined NXSRC
  double dclock(), time, vals[3], work[3];
#elif defined MPISRC
  double MPI_Wtime(), time, vals[3], work[3];
#endif
#endif


#ifdef DEBUG
  error_msg_warning("XYT_solve() :: start %d\n",n_xyt_handles);
#endif

#ifdef SAFE
  check_init();
  check_handle(xyt_handle);
#endif

#ifdef TIMINGS
#if   defined NXSRC
  time = dclock();
#elif defined MPISRC
  time = MPI_Wtime();
#endif
#endif

  /* don't overwrite b!!! */
  rvec_copy(x,b,xyt_handle->mvi->n);
  do_xyt_solve(xyt_handle,x);

#ifdef TIMINGS
#if   defined NXSRC
  time = dclock() - time;
#elif defined MPISRC
  time = MPI_Wtime() - time;
#endif

  vals[0]=vals[1]=vals[2]=time;
  grop(vals,work,sizeof(op)/sizeof(op[0])-1,op);

  if (!my_id)
    {
      printf("%d :: min   xyt_slv=%g\n",my_id,vals[0]);
      printf("%d :: max   xyt_slv=%g\n",my_id,vals[1]);
      printf("%d :: avg   xyt_slv=%g\n",my_id,vals[2]/num_nodes);
    }
#endif


#ifdef DEBUG
  error_msg_warning("XYT_solve() :: end   %d\n",n_xyt_handles);
#endif
}



/*************************************xyt.c************************************
Function: xyt_free_()

Input :
Output:
Return:
Description:
**************************************xyt.c***********************************/
#if defined UPCASE
void
XYT_FREE (int n)
#else
void
xyt_free_(int n)
#endif
{
  xyt_ADT xyt_handle;


  fn_xyt_handles--;
  xyt_handle=fhandles[n];
  fhandles[n]=NULL; 
  XYT_free(xyt_handle);
}



/*************************************xyt.c************************************
Function: XYT_free()

Input :
Output:
Return:
Description:
**************************************xyt.c***********************************/
void 
XYT_free(xyt_ADT xyt_handle)
{
#ifdef DEBUG
  error_msg_warning("XYT_free() :: start %d\n",n_xyt_handles);
#endif

  check_init();
  check_handle(xyt_handle);
  n_xyt_handles--;

  bss_free(xyt_handle->info->nsep);
  bss_free(xyt_handle->info->lnsep);
  bss_free(xyt_handle->info->fo);
  bss_free(xyt_handle->info->stages);
  bss_free(xyt_handle->info->solve_uu);
  bss_free(xyt_handle->info->solve_w);
  bss_free(xyt_handle->info->x);
  bss_free(xyt_handle->info->xcol_vals);
  bss_free(xyt_handle->info->xcol_sz);
  bss_free(xyt_handle->info->xcol_indices);
  /*bss_free(xyt_handle->info->y);*/
  bss_free(xyt_handle->info->ycol_vals);
  bss_free(xyt_handle->info->ycol_sz);
  bss_free(xyt_handle->info->ycol_indices);
  bss_free(xyt_handle->mvi->local2global);
  bss_free(xyt_handle->info);
  gs_free(xyt_handle->mvi->gs_handle);
  bss_free(xyt_handle->mvi);
  bss_free(xyt_handle);
 
#ifdef DEBUG
  error_msg_warning("perm frees = %d\n",perm_frees());
  error_msg_warning("perm calls = %d\n",perm_calls());
  error_msg_warning("bss frees  = %d\n",bss_frees());
  error_msg_warning("bss calls  = %d\n",bss_calls());
  error_msg_warning("XYT_free() :: end   %d\n",n_xyt_handles);
#endif
}


#ifdef MLSRC
/*************************************xyt.c************************************
Function: ML_XYT_factor()

Input :
Output:
Return:
Description:

ML requires that the solver call be checked in
**************************************xyt.c***********************************/
int 
ML_XYT_factor(xyt_ADT xyt_handle,  /* prev. allocated xyt  handle */
	      int *local2global,   /* global column mapping       */
	      int n,               /* local num rows              */
	      int m,               /* local num cols              */
	      void *matvec,        /* b_loc=A_local.x_loc         */
	      void *grid_data,     /* grid data for matvec        */
	      int grid_tag,        /* grid tag for ML_Set_CSolve  */
	      ML *my_ml            /* ML handle                   */
	      )
{
  int flag;
#ifdef ML_TIMING
   double t0;

   t0 = GetClock();
#endif


#ifdef DEBUG
  error_msg_warning("ML_XYT_factor() :: start %d\n",n_xyt_handles);
#endif

#ifdef SAFE
  check_init();
  check_handle(xyt_handle);
  if (my_ml->comm->ML_mypid!=my_id)
    {error_msg_fatal("ML_XYT_factor bad my_id %d\t%d\n",
		     my_ml->comm->ML_mypid,my_id);}
  if (my_ml->comm->ML_nprocs!=num_nodes)
    {error_msg_fatal("ML_XYT_factor bad np %d\t%d\n",
		     my_ml->comm->ML_nprocs,num_nodes);}
#endif

  my_ml->SingleLevel[grid_tag].csolve->func->external = ML_XYT_solve;
  my_ml->SingleLevel[grid_tag].csolve->func->ML_id = ML_EXTERNAL;
  my_ml->SingleLevel[grid_tag].csolve->data = xyt_handle;

  /* done ML specific stuff ... back to reg sched pgm */

  flag = XYT_factor(xyt_handle, local2global, n, m, matvec, grid_data);
#ifdef ML_TIMING
   t0 = GetClock() - t0;
   my_ml->SingleLevel[grid_tag].csolve->build_time += t0;
   my_ml->timing->total_build_time += t0;
#endif

#ifdef DEBUG
  error_msg_warning("ML_XYT_factor() :: end   %d (flag=%d)\n",n_xyt_handles,flag);
#endif
  return(flag);
}



/*************************************xyt.c************************************
Function: ML_XYT_solve

Input :
Output:
Return:
Description:
**************************************xyt.c***********************************/
void 
ML_XYT_solve(xyt_ADT xyt_handle, int lx, double *sol, int lb, double *rhs)
{
  XYT_solve(xyt_handle, sol, rhs);
}



/*************************************xyt.c************************************
Function: ML_XYT_free()

Input :
Output:
Return:
Description:
**************************************xyt.c***********************************/
void *
ML_XYT_free(xyt_ADT xyt_handle)
{
  void *temp;

  temp = xyt_handle->mvi->grid_data;
  XYT_free(xyt_handle);
  return(temp);
}
#endif



/*************************************Xyt.c************************************
Function: check_init

Input :
Output:
Return:
Description:
**************************************xyt.c***********************************/
static
void
check_init(void)
{
#ifdef DEBUG
  error_msg_warning("check_init() :: start %d\n",n_xyt_handles);
#endif

  comm_init();
  /*
  perm_init(); 
  bss_init();
  */

#ifdef DEBUG
  error_msg_warning("check_init() :: end   %d\n",n_xyt_handles);
#endif
}



/*************************************xyt.c************************************
Function: check_handle()

Input :
Output:
Return:
Description:
**************************************xyt.c***********************************/
static
void 
check_handle(xyt_ADT xyt_handle)
{
  int vals[2], work[2], op[] = {NON_UNIFORM,GL_MIN,GL_MAX};


#ifdef DEBUG
  error_msg_warning("check_handle() :: start %d\n",n_xyt_handles);
#endif

  if (xyt_handle==NULL)
    {error_msg_fatal("check_handle() :: bad handle :: NULL %d\n",xyt_handle);}

  vals[0]=vals[1]=xyt_handle->id;
  giop(vals,work,sizeof(op)/sizeof(op[0])-1,op);
  if ((vals[0]!=vals[1])||(xyt_handle->id<=0))
    {error_msg_fatal("check_handle() :: bad handle :: id mismatch min/max %d/%d %d\n",
		     vals[0],vals[1], xyt_handle->id);}

#ifdef DEBUG
  error_msg_warning("check_handle() :: end   %d\n",n_xyt_handles);
#endif
}



/*************************************xyt.c************************************
Function: det_separators

Input :
Output:
Return:
Description:


  det_separators(xyt_handle, local2global, n, m, mylocmatvec, grid_data);

**************************************xyt.c***********************************/
static 
void 
det_separators(xyt_ADT xyt_handle)
{
  int i, ct, id;
  int mask, edge, *iptr; 
  int *dir, *used;
  int sum[4], w[4];
  REAL rsum[4], rw[4];
  int op[] = {GL_ADD,0};
  REAL *lhs, *rhs;
  int *nsep, *lnsep, *fo, nfo=0;
  gs_ADT gs_handle=xyt_handle->mvi->gs_handle;
  int *local2global=xyt_handle->mvi->local2global;
  int  n=xyt_handle->mvi->n;
  int  m=xyt_handle->mvi->m;
  int level=xyt_handle->level;
  int shared=0; 

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

  /* determine the # of unique dof */
  /* will this work for the non-symmetric case!?! */
  rvec_zero(lhs,m);
  rvec_set(lhs,1.0,n);
  gs_gop_hc(gs_handle,lhs,"+\0",level);
  rvec_zero(rsum,2);
  for (ct=i=0;i<n;i++)
    {
      if (lhs[i]!=0.0)
	{rsum[0]+=1.0/lhs[i]; rsum[1]+=lhs[i];}
    }
  grop_hc(rsum,rw,2,op,level);
  rsum[0]+=0.1;
  rsum[1]+=0.1;
  if (!my_id)
    {
      printf("xyt n unique = %d (%g)\n",(int) rsum[0], rsum[0]);
      printf("xyt n shared = %d (%g)\n",(int) rsum[1], rsum[1]);
    }

  if (fabs(rsum[0]-rsum[1])>EPS)
    {
      error_msg_warning("determine separators - SHARED!!!\n"); 
      shared=TRUE;
    }
  else
    {
      error_msg_warning("determine separators - NONSHARED!!!\n"); 
      shared=FALSE;
    }

  xyt_handle->info->n_global=xyt_handle->info->m_global=(int) rsum[0];
  xyt_handle->mvi->n_global =xyt_handle->mvi->m_global =(int) rsum[0];

  /* determine separator sets top down */
  if (shared)
    {
      /* solution is to do as in the symmetric shared case but then */
      /* pick the sub-hc with the most free dofs and do a mat-vec   */
      /* and pick up the responses on the other sub-hc from the     */
      /* initial separator set obtained from the symm. shared case  */
      error_msg_fatal("shared dof separator determination not ready ... see hmt!!!\n"); 
      for (iptr=fo+n,id=my_id,mask=num_nodes>>1,edge=level;edge>0;edge--,mask>>=1)
	{
	  /* set rsh of hc, fire, and collect lhs responses */
	  (id<mask) ? rvec_zero(lhs,m) : rvec_set(lhs,1.0,m);
	  gs_gop_hc(gs_handle,lhs,"+\0",edge);
	  
	  /* set lsh of hc, fire, and collect rhs responses */
	  (id<mask) ? rvec_set(rhs,1.0,m) : rvec_zero(rhs,m);
	  gs_gop_hc(gs_handle,rhs,"+\0",edge);
	  
	  for (i=0;i<n;i++)
	    {
	      if (id< mask)
		{		
		  if (lhs[i]!=0.0)
		    {lhs[i]=1.0;}
		}
	      if (id>=mask)
		{		
		  if (rhs[i]!=0.0)
		    {rhs[i]=1.0;}
		}
	    }

	  if (id< mask)
	    {gs_gop_hc(gs_handle,lhs,"+\0",edge-1);}
	  else
	    {gs_gop_hc(gs_handle,rhs,"+\0",edge-1);}

	  /* count number of dofs I own that have signal and not in sep set */
	  rvec_zero(rsum,4);
	  for (ivec_zero(sum,4),ct=i=0;i<n;i++)
	    {
	      if (!used[i]) 
		{
		  /* number of unmarked dofs on node */
		  ct++;
		  /* number of dofs to be marked on lhs hc */
		  if (id< mask)
		    {		
		      if (lhs[i]!=0.0)
			{sum[0]++; rsum[0]+=1.0/lhs[i];}
		    }
		  /* number of dofs to be marked on rhs hc */
		  if (id>=mask)
		    {		
		      if (rhs[i]!=0.0)
			{sum[1]++; rsum[1]+=1.0/rhs[i];}
		    }
		}
	    }

	  /* go for load balance - choose half with most unmarked dofs, bias LHS */
	  (id<mask) ? (sum[2]=ct) : (sum[3]=ct);
	  (id<mask) ? (rsum[2]=ct) : (rsum[3]=ct);
	  giop_hc(sum,w,4,op,edge);
	  grop_hc(rsum,rw,4,op,edge);
	  rsum[0]+=0.1; rsum[1]+=0.1; rsum[2]+=0.1; rsum[3]+=0.1;

	  if (id<mask)
	    {
	      /* mark dofs I own that have signal and not in sep set */
	      for (ct=i=0;i<n;i++)
		{
		  if ((!used[i])&&(lhs[i]!=0.0))
		    {
		      ct++; nfo++;

		      if (nfo>n)
			{error_msg_fatal("nfo about to exceed n\n");}

		      *--iptr = local2global[i];
		      used[i]=edge;
		    }
		}
	      if (ct>1) {ivec_sort(iptr,ct);}

	      lnsep[edge]=ct;
	      nsep[edge]=(int) rsum[0];
	      dir [edge]=LEFT;
	    }

	  if (id>=mask)
	    {
	      /* mark dofs I own that have signal and not in sep set */
	      for (ct=i=0;i<n;i++)
		{
		  if ((!used[i])&&(rhs[i]!=0.0))
		    {
		      ct++; nfo++;

		      if (nfo>n)
			{error_msg_fatal("nfo about to exceed n\n");}

		      *--iptr = local2global[i];
		      used[i]=edge;
		    }
		}
	      if (ct>1) {ivec_sort(iptr,ct);}

	      lnsep[edge]=ct;
	      nsep[edge]= (int) rsum[1];
	      dir [edge]=RIGHT;
	    }

	  /* LATER or we can recur on these to order seps at this level */
	  /* do we need full set of separators for this?                */

	  /* fold rhs hc into lower */
	  if (id>=mask)
	    {id-=mask;}
	}
    }
  else
    {
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
		  if ((id< mask)&&(lhs[i]!=0.0)) {sum[0]++;}
		  /* number of dofs to be marked on rhs hc */
		  if ((id>=mask)&&(rhs[i]!=0.0)) {sum[1]++;}
		}
	    }

	  /* for the non-symmetric case we need separators of width 2 */
	  /* so take both sides */
	  (id<mask) ? (sum[2]=ct) : (sum[3]=ct);
	  giop_hc(sum,w,4,op,edge);

	  ct=0;
	  if (id<mask)
	    {
	      /* mark dofs I own that have signal and not in sep set */
	      for (i=0;i<n;i++)
		{
		  if ((!used[i])&&(lhs[i]!=0.0))
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
	      /* LSH hc summation of ct should be sum[0] */
	    }
	  else
	    {
	      /* mark dofs I own that have signal and not in sep set */
	      for (i=0;i<n;i++)
		{
		  if ((!used[i])&&(rhs[i]!=0.0))
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
	      /* RSH hc summation of ct should be sum[1] */
	    }

	  if (ct>1) {ivec_sort(iptr,ct);}
	  lnsep[edge]=ct;
	  nsep[edge]=sum[0]+sum[1];
	  dir [edge]=BOTH;

	  /* LATER or we can recur on these to order seps at this level */
	  /* do we need full set of separators for this?                */

	  /* fold rhs hc into lower */
	  if (id>=mask)
	    {id-=mask;}
	}
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
  dir  [edge]=BOTH;

#ifdef DEBUG  
  if (nfo!=n)
    {error_msg_fatal("hey nfo!= n?\n");}

  if (iptr!=fo)
    {error_msg_fatal("iptr and fo misaligned?\n");}
#endif

#ifdef INFO
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

  xyt_handle->info->nsep=nsep;
  xyt_handle->info->lnsep=lnsep;
  xyt_handle->info->fo=fo;
  xyt_handle->info->nfo=nfo;

  bss_free(dir);
  bss_free(lhs);
  bss_free(rhs);
  bss_free(used);


  fflush(stdout);

#ifdef DEBUG  
  error_msg_warning("det_separators() :: end\n");
#endif
}



/*************************************xyt.c************************************
Function: set_mvi

Input :
Output:
Return:
Description:
**************************************xyt.c***********************************/
static
mv_info *set_mvi(int *local2global, int n, int m, void *matvec, void *grid_data)
{
  mv_info *mvi;


#ifdef DEBUG
  error_msg_warning("set_mvi() :: start\n");
#endif

  mvi = bss_malloc(sizeof(mv_info));
  mvi->n=n;
  mvi->m=m;
  mvi->n_global=-1;
  mvi->m_global=-1;
  mvi->local2global=bss_malloc((m+1)*INT_LEN);
  ivec_copy(mvi->local2global,local2global,m);
  mvi->local2global[m] = INT_MAX;
  mvi->matvec=matvec;
  mvi->grid_data=grid_data;

#ifdef NOT
  perm_init();
  bss_init();
  mvi->gs_handle = gs_init(local2global, m, num_nodes);
  gs_free(mvi->gs_handle);
  error_msg_warning("perm frees = %d\n",perm_frees());
  error_msg_warning("perm calls = %d\n",perm_calls());
  error_msg_warning("bss frees  = %d\n",bss_frees());
  error_msg_warning("bss calls  = %d\n",bss_calls());
#endif

  /* set xyt communication handle to perform restricted matvec */
  mvi->gs_handle = gs_init(local2global, m, num_nodes);

#ifdef DEBUG
  error_msg_warning("set_mvi() :: end   \n");
#endif
  
  return(mvi);
}



/*************************************xyt.c************************************
Function: set_mvi

Input :
Output:
Return:
Description:

      computes u = A.v 
      do_matvec(xyt_handle->mvi,v,u);
**************************************xyt.c***********************************/
static
void do_matvec(mv_info *A, REAL *v, REAL *u)
{
  A->matvec(A->grid_data,v,u);
}



