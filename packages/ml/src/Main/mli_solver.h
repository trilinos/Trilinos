/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* ******************************************************************** */
/* Declaration of the MLI_Solver data structure                         */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : December, 1999                                       */
/* ******************************************************************** */

#ifndef _MLSOLVERIFACE_
#define _MLSOLVERIFACE_

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include "ml_common.h"
#include "ml_comm.h"
#include "ml_aggregate.h"
#include "ml_amg.h"

typedef struct
{
    int      Nrows;
    int      *rowptr;
    int      *colnum;
    int      *map;
    double   *values;
    int      sendProcCnt;
    int      *sendProc;
    int      *sendLeng;
    int      **sendList;
    int      recvProcCnt;
    int      *recvProc;
    int      *recvLeng;
}
MLI_CSRMatrix;

#ifdef ML_MPI
#include <mpi.h>
#define ML_MPI_Request MPI_Request
#define ML_MPI_Comm    MPI_Comm
#else
#define ML_MPI_Request int
#define ML_MPI_Comm    int
#endif

typedef struct
{
    MLI_CSRMatrix  *Amat;
    ML_MPI_Comm    comm;
    int            globalEqns;
    int            *partition;
}
MLI_Context;

typedef struct
{
    ML_MPI_Comm  comm;
    ML           *ml_ptr;
    int          nlevels;
    int          method;
    int          pre, post;
    int          pre_sweeps, post_sweeps;
    double       jacobi_wt;
    double       ag_threshold;
    int          ag_coarsen;
    int          ag_method;
    int          ndiag;
    double       *diag_scale;
    ML_Aggregate *ml_ag;
    ML_AMG       *ml_amg;
    MLI_Context  *contxt;
    int          nRows;
    int          *mat_ia;
    int          *mat_ja;
    double       *mat_a;
    int          nPDE;
    int          nNullVectors;
    double       *nullSpace;
    double       *rhs;
    double       *sol;
} MLI_Solver;

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

MLI_Solver *MLI_Solver_Create( ML_MPI_Comm );
int  MLI_Solver_Destroy( MLI_Solver * );
int  MLI_Solver_Setup(MLI_Solver *, double *x );
int  MLI_Solver_SetupDD(MLI_Solver *,int startRow, int Nrows, int *mat_ia,
                      int *mat_ja, double *mat_a, double *b, double *x );
int  MLI_Solver_Solve( MLI_Solver *solver );
int  MLI_Solver_Set_MLNumLevels( MLI_Solver *, int nlevels );
int  MLI_Solver_Set_KrylovMethod( MLI_Solver *, int method );
int  MLI_Solver_Set_StrongThreshold( MLI_Solver *, double strong_threshold );
int  MLI_Solver_Set_NumPreSmoothings( MLI_Solver *, int num_sweeps );
int  MLI_Solver_Set_NumPostSmoothings( MLI_Solver *, int num_sweeps );
int  MLI_Solver_Set_PreSmoother( MLI_Solver *, int smoother_type );
int  MLI_Solver_Set_PostSmoother( MLI_Solver *, int smoother_type );
int  MLI_Solver_Set_DampingFactor( MLI_Solver *, double factor );
int  MLI_Solver_Set_MGMethod( MLI_Solver *, int );
int  MLI_Solver_Set_CoarsenScheme( MLI_Solver *, int );
int  MLI_Solver_Get_IJAFromFile(MLI_Solver *, char *matfile, char *rhsfile);
int  MLI_Solver_Get_NullSpaceFromFile(MLI_Solver *, char *rbmfile);
int MLI_Irecv(void* buf, unsigned int count, int *src, int *mid,
             ML_MPI_Comm comm, ML_MPI_Request *request );
int MLI_SIrecv(void* buf, unsigned int count, int *src, int *mid,
             ML_MPI_Comm comm, ML_MPI_Request *request );
int MLI_Wait(void* buf, unsigned int count, int *src, int *mid,
             ML_MPI_Comm comm, ML_MPI_Request *request );
int MLI_SWait(void* buf, unsigned int count, int *src, int *mid,
             ML_MPI_Comm comm, ML_MPI_Request *request );
int MLI_Send(void* buf, unsigned int count, int dest, int mid, ML_MPI_Comm comm );
int MLI_SSend(void* buf, unsigned int count, int dest, int mid, ML_MPI_Comm comm );
int MLI_CSRExchBdry(double *vec, void *obj);
int MLI_CSRMatVec(ML_Operator *obj, int leng1, double p[], int leng2, double ap[]);
int MLI_CSRGetRow(ML_Operator *obj, int N_requested_rows, int requested_rows[],
    int allocated_space, int columns[], double values[], int row_lengths[]);
void MLI_Solver_Read_IJAFromFile(double **val, int **ia, int **ja, int *N,
                                double **rhs, char *matfile, char *rhsfile);
extern int MLI_Solver_Construct_CSRMatrix(int, int*, int*, double*,
                         MLI_CSRMatrix *, MLI_Solver*, int *,MLI_Context*);
extern int MLI_Solver_Construct_LocalCSRMatrix(int nrows, int *mat_ia,
                         int *mat_ja, double *mat_a, MLI_CSRMatrix *mli_mat,
                         MLI_Solver *solver, int *partition, MLI_Context *obj);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif
