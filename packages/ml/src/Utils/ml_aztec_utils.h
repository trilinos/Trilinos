/********************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/********************************************************************* */
/*          Utilities for Aztec/ML users                               */
/********************************************************************* */

#ifndef __MLAZUTILS__
#define __MLAZUTILS__

#ifdef __cplusplus
extern "C" {
#endif

#ifdef AZTEC
#ifdef ML_MPI
#define AZ_MPI
#endif
#define AZ_ONLY_PRECONDITIONER  -7778

#include "az_aztec.h"
struct aztec_context {
   AZ_MATRIX *Amat; 
   AZ_PRECOND *Prec;
   int    *proc_config; 
   int    *options;
   double *params;
   double *status; 
   void   *getrowstuff;
   int    prec_or_solver;
#ifdef AZ_ver2_1_0_3
   struct AZ_SCALING *scaling;
#endif
};
#ifndef AZTEC2_0
extern void ML_Gen_SmootherAztec(ML *ml_handle, int level, int options[],
        double params[], int proc_config[], double status[], int N_iter,
	int pre_or_post, void (*prec_function)(
        double *, int *, int *, double *, struct AZ_MATRIX_STRUCT  *,
        struct AZ_PREC_STRUCT *));

extern int az_wrap_solvers(void *mydata, int in, double x[], int out, 
	double rhs[]);

extern void AZ_ML_SmootherClean(void *data);

#else
#ifdef ML_MPI
#include <mpi.h>
#define MPI_AZRequest MPI_Request
#define MPI_AZComm    MPI_Comm
#else
#define MPI_AZRequest int
#define MPI_AZComm    int
#endif

#define AZ_fixed_pt 33
extern struct AZ_MATRIX_STRUCT *AZ_matrix_create(int local);
extern struct AZ_PREC_STRUCT   *AZ_precond_create(struct AZ_MATRIX_STRUCT *Pmat,
                void (*prec_fun)( double *, int *, int *, double *,
                      struct AZ_MATRIX_STRUCT  *, struct AZ_PREC_STRUCT *) );

extern void AZ_matrix_destroy( struct AZ_MATRIX_STRUCT **Amat);
extern void AZ_precond_destroy(struct AZ_PREC_STRUCT **precond);
extern void AZ_set_MSR(AZ_MATRIX *Amat, int bindx[], double val[],
        int data_org[]);

extern void AZ_set_VBR(AZ_MATRIX *Amat, int rpntr[], int cpntr[], int bpntr[],
        int indx[], int bindx[], double val[], int data_org[]);


#define mdwrap_wait(a,b,c,x,y,z)   md_wrap_wait(a,b,c,(x),(y),(z))
#define mdwrap_iwrite(a,b,c,x,y,z) md_wrap_iwrite(a,b,c,(x),(y),(z))
#define mdwrap_iread(a,b,c,x,y)   md_wrap_iread((a),(b),(c),(x),(y))
#define mdwrap_write(a,b,c,x,y)   md_wrap_write((a),(b),(c),(x),(y))

#if   defined(caps)
#   define az_set_proc_config_        AZ_SET_PROC_CONFIG
#elif defined(matched)
#   define az_set_proc_config_        az_set_proc_config
#endif

extern void AZ_set_proc_config(int proc_config[], MPI_AZComm );


#endif
extern void AZ_mlcomm2data_org(ML_CommInfoOP *comm_info,int *data_org[]);

extern void AZ_set_ML_preconditioner(AZ_PRECOND **Precond, AZ_MATRIX *Amat,
				ML *ml_handle, int options[]);

extern int AZ_get_MSR_arrays(ML_Operator *, int **bindx, double **val);

extern int az_comm_wrapper(double vector[], void *data);

extern int az_msrgetrow_wrapper(void *data, int N_requested_rows, 
   int requested_rows[], int allocated_space, int columns[], double values[], 
   int row_lengths[]);

extern int az_vbrgetrow_wrapper(void *data, int N_requested_rows, 
   int requested_rows[], int allocated_space, int columns[], double values[], 
   int row_lengths[]);

extern int az_usergetrow_wrapper(void *data, int N_requested_rows, 
   int requested_rows[], int allocated_space, int columns[], double values[], 
   int row_lengths[]);

extern int az_matvec_wrapper(void *data, int, double p[], int, double ap[]);

extern int AZ_ML_Set_Amat(ML *ml_handle, int level, int isize, int osize,
	AZ_MATRIX *Amat, int *proc_config);

extern void AZ_ML_set_vbrdiagonal(ML *ml, int level, AZ_MATRIX *matrix);

extern void AZ_ML_set_userdiagonal(ML *ml, int level, AZ_MATRIX *matrix);

extern void ML_precondition(double ff[], int options[], int proc_config[],
                     double params[], AZ_MATRIX *mat, AZ_PRECOND *prec);

extern void AZ_ML_Clean(void *data);
#endif

#ifdef __cplusplus
}
#endif

#ifdef AZTEC_RAY_WILL_FIX

#include "ml_include.h"
#include "az_aztec.h"
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/
#endif
#endif

