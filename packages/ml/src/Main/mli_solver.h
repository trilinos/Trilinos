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

#include "ml_include.h"

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

#ifndef ML_MPI
#define MPI_Comm int
#endif

typedef struct
{
    MLI_CSRMatrix  *Amat;
    MPI_Comm       comm;
    int            globalEqns;
    int            *partition;
}
MLI_Context;

typedef struct 
{
    MPI_Comm     comm;
    ML           *ml_ptr;
    int          nlevels;
    int          method;
    int          pre, post;
    int          pre_sweeps, post_sweeps;
    double       jacobi_wt;
    double       ag_threshold;
    int          ag_coarsen;
    int          ndiag;
    double       *diag_scale;
    ML_Aggregate *ml_ag;
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

#ifdef __cplusplus
extern "C" {
#endif

MLI_Solver *MLI_Solver_Create( MPI_Comm );
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
int  MLI_Solver_Set_CoarsenScheme( MLI_Solver *, int );
int  MLI_Solver_Get_IJAFromFile(MLI_Solver *, char *matfile, char *rhsfile);
int  MLI_Solver_Get_NullSpaceFromFile(MLI_Solver *, char *rbmfile);

#ifdef __cplusplus
}
#endif

#endif

