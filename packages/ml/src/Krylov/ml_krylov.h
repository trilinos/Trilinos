/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Declaration of the ML_KrylovData structure                           */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : December, 1999                                       */
/* ******************************************************************** */

#ifndef __MLKRYLOVDATA_
#define __MLKRYLOVDATA_
#define ML_CG 0
#define ML_GMRES 1

/* ******************************************************************** */
/* local include files                                                  */
/* ******************************************************************** */

typedef struct ML_Krylov_Struct ML_Krylov;

/* #include "ml_struct.h" */
#include "ml_common.h"
#include "ml_comm.h"
#include "ml_operator.h"

/* -------------------------------------------------------------------- */
/* This data structure defines an enriched operator class for the       */
/* specification of the discretization matrix, the restriction and the  */
/* prolongation operator.                                               */
/* -------------------------------------------------------------------- */

struct ML_Krylov_Struct 
{
   int           ML_id;
   int           ML_method;
   int           ML_gmres_dim;
   int           ML_cgstabl_dim;
   int           ML_max_iterations;
   int           ML_print_freq;
   double        ML_tolerance;
   double        *diag_scale;
   ML_Operator   *ML_matrix;
   void          *ML_precon;
   ML_Comm       *ML_com;
   int           ML_eigen;
   int           ML_nonsym_eigen;
   double        ML_eigen_max;
   double        ML_eigen_min;
   int           (*ML_precfcn)(void*, int, double*, int, double*);
   int           ML_dont_scale_by_diag;
};

/* ******************************************************************** */
/* ******************************************************************** */
/*      User Interface Proto-types                                      */
/* ******************************************************************** */
/* ******************************************************************** */

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

extern ML_Krylov *ML_Krylov_Create(ML_Comm *);
extern int ML_Krylov_Destroy(ML_Krylov **);
extern ML_Comm *ML_Krylov_Get_Comm(ML_Krylov *data);
extern int ML_Krylov_Set_PrintFreq(ML_Krylov *, int n);
extern int ML_Krylov_Get_PrintFreq(ML_Krylov *);
extern int ML_Krylov_Set_Method(ML_Krylov *, int method);
extern int ML_Krylov_Set_GMRESSize(ML_Krylov *, int size);
extern int ML_Krylov_Get_GMRESSize(ML_Krylov *);
extern int ML_Krylov_Set_BICGSTABLSize(ML_Krylov *, int size);
extern int ML_Krylov_Get_BICGSTABLSize(ML_Krylov *);
extern int ML_Krylov_Set_Amatrix(ML_Krylov *, ML_Operator *mat);
extern ML_Operator* ML_Krylov_Get_Amatrix(ML_Krylov *);
extern int ML_Krylov_Set_MaxIterations(ML_Krylov *, int iter);
extern int ML_Krylov_Get_MaxIterations(ML_Krylov *);
extern int ML_Krylov_Set_Tolerance(ML_Krylov *, double tol);
extern int ML_Krylov_Set_Diagonal(ML_Krylov *, int leng, double *diag);
extern double ML_Krylov_Get_Tolerance(ML_Krylov *);
extern int ML_Krylov_Set_Precon(ML_Krylov *, void *);
extern int ML_Krylov_Set_PreconFunc(ML_Krylov*,
                                    int (*func)(void*,int,double*,int,double*));
extern void *ML_Krylov_Get_Precon(ML_Krylov *);
extern int ML_Krylov_Set_ComputeEigenvalues(ML_Krylov *);
extern int ML_Krylov_Set_ComputeNonSymEigenvalues(ML_Krylov *);
extern int ML_Krylov_Set_DiagScaling_Eig(ML_Krylov *data, int scale);

extern double ML_Krylov_Get_MaxEigenvalue(ML_Krylov *);

extern int ML_Krylov_Solve(ML_Krylov *, int, double *, double*);
extern int ML_MGVSolve_Wrapper(void *, int, double *, int, double*);
extern int ML_AMGVSolve_Wrapper(void *, int, double *, int, double*);
extern int ML_DiagScale_Wrapper(void *, int, double *, int, double*);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif


