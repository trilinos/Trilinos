/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
/* ******************************************************************** */
/* Declaration of the segregated structure                              */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : March, 1999                                          */
/* ******************************************************************** */

#ifndef __SEGSTRUCT__
#define __SEGSTRUCT__

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/* ******************************************************************** */
/* local include files                                                  */
/* ******************************************************************** */

#include "ml_common.h"
#include "ml_defs.h"
#include "ml_bdrypts.h"
#include "ml_mapper.h"
#include "ml_grid.h"
#include "ml_smoother.h"
#include "ml_comminfoop.h"
#include "ml_1level.h"
#include "ml_operator.h"
#include "ml_csolve.h"
#include "ml_operatoragx.h"
#include "ml_comm.h"
#include "ml_gridfunc.h"
#include "ml_vec.h"
#include "ml_rap.h"
#include "ml_utils.h"
#include "ml_mat_formats.h"
#include "ml_solver.h"
#include "az_aztec.h"
#include "az_aztec_defs.h"
#include <string.h>

#define ML_SEG_DIAGONAL 0
#define ML_SEG_UPPER_TRIANGULAR 1
#define ML_SEG_LOWER_TRIANGULAR 2
#define ML_SEG_DIAG_ELT 0
#define ML_SEG_OFFDIAG_ELT 1

/* ******************************************************************** */
/* ******************************************************************** */
/* data definition for the ML_SEG Class                                 */
/* ******************************************************************** */
/* -------------------------------------------------------------------- */
/* This is the primary data structure which users interact directly     */
/* with via the provided functions.                                     */
/* -------------------------------------------------------------------- */

struct ML_SEG_MATRIX_Struct
{
   AZ_MATRIX  *Amat;
   AZ_PRECOND *precond;
   int        *AZ_options;     /* Aztec's options array for preconditioning */
   double     *AZ_params;
   struct AZ_SCALING *AZ_scale;
   int        Amat_changed;
   int        precond_changed;
};

struct SEG_Struct
{
   int    SEG_nblocks;        /* number of block rows/cols in the matrix */
   int    SEG_noffdiags;      /* number of off-diagonal blocks */
   int    SEG_format;         /* SEG_DIAG, SEG_UPPER, or SEG_LOWER */
   int    **SEG_rowlists;     /* lists of rows in each block */
   int    *SEG_rowlist_lengs; /* lengths of each block row */
   int    SEG_total_nrows;    /* total number of rows in the matrix */
   struct ML_SEG_MATRIX_Struct **SEG_diag_list; /* matrices on the diagonal */
   struct ML_SEG_MATRIX_Struct **SEG_offdiag_list; /* matrices off the diagonal */
};

typedef   struct ML_SEG_Stuct           ML_SEG;
typedef   struct ML_SEG_MATRIX_Struct   ML_SEG_MATRIX;

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
extern int  ML_SEG_Create(struct ML_SEG_Struct **seg, AZ_MATRIX *Amat,
                          int nblocks, int **rowlist, int *rowlengs,
                          int format, int proc_config[]);
extern void ML_SEG_Destroy(struct ML_SEG_Struct **seg);
extern struct ML_SEG_MATRIX_Struct *ML_SEG_Matrix_Create(AZ_MATRIX *Amat,
                          int nrows, int *rowlist, int ncols, int *collist,
                          int format, int proc_config[]);
extern void ML_SEG_Matrix_Destroy(struct SEG_MATRIX_Struct **Smat);
extern void ML_SEG_Precondition(double ff[], int options[], int proc_config[],
                          double params[], AZ_MATRIX *mat, AZ_PRECOND *prec);
extern void AZ_Set_SEG_Preconditioner(AZ_PRECOND *Precond, struct SEG_Struct *seg,
                          int options[]);
extern void ML_SEG_Set_ML_Precond(struct ML_SEG_Struct *seg, ML *ml, int block_row,
                          double params[], int options[], int proc_config[]);
extern void ML_SEG_ML_Set_Amat(struct ML_SEG_Struct *seg, ML *ml, int level,
                          int block_row, int proc_config[]);
extern void ML_SEG_Replace_Submat(struct ML_SEG_Struct *seg, int block_row,
                          int block_col, AZ_MATRIX *newMat, int proc_config[]);
extern void ML_SEG_Set_AZ_Precond(struct ML_SEG_Struct *seg, int block_row,
                          double *params, int *options, int proc_config[]);

/*
extern int ML_Set_ResidualOutputFrequency(ML *ml, int output_freq);
extern int ML_Set_Tolerance(ML *ml, double tolerance);

extern int ML_Destroy(ML **ml);

extern int ML_Init_Comm(ML *ml);
extern int ML_Set_Comm_MyRank(ML *ml, int myrank);
extern int ML_Set_Comm_Nprocs(ML *ml, int nprocs);
extern int ML_Set_Comm_Communicator(ML *ml, USR_COMM com);
extern int ML_Set_Comm_Send(ML *ml, int (*send)());
extern int ML_Set_Comm_Recv(ML *ml, int (*recv)());
extern int ML_Set_Comm_Wait(ML *ml, int (*wait)());

extern int ML_Set_Comm(ML *ml, ML_Comm *comm);

extern int ML_Init_Grid(ML *, int nl, void *grid);
extern int ML_Set_Grid_GridFunc(ML *, int nl, ML_GridFunc *);
extern int ML_Set_Grid_MaxVertPerElmnt(ML *, int, int nvert);
extern int ML_Set_Grid_GetDimension(ML *, int nl, int (*func)());
extern int ML_Set_Grid_GetNVert(ML *, int nl, int (*func)());
extern int ML_Set_Grid_GetNElmnt(ML *, int nl, int (*func)());
extern int ML_Set_Grid_GetElmntNVert(ML *, int nl, int (*func)());
extern int ML_Set_Grid_GetElmntVertList(ML *, int nl, int (*func)());
extern int ML_Set_Grid_GetVertGlobalNum(ML *, int nl, int (*func)());
extern int ML_Set_Grid_GetElmntGlobalNum(ML *, int nl, ml_big_int (*func)());
extern int ML_Set_Grid_GetVertCoordinate(ML *, int nl, int (*func)());
extern int ML_Set_Grid_ComputeBasisCoef(ML *, int nl, int (*func)());
extern int ML_Set_Grid_GetElmntVolume(ML *, int nl, int (*func)());
extern int ML_Set_Grid_GetElmntMatrix(ML *, int nl, int (*func)());
extern int ML_Set_Grid_GetElmntNullSpace(ML *, int, int (*func)());

extern int ML_Gen_GridXsferUsingFEBasis(ML *, int L1, int L2, int stride);
extern int ML_Gen_MGHierarchyVanek(ML *, int start, int increment_or_decrement);

extern int ML_Set_Grid(ML *, int nl, void *grid, ML_GridFunc *);

extern int ML_Init_Amatrix(ML *,int level,int ilen,int olen,void *data);
extern int ML_Set_Amatrix_Matvec(ML*,int,
                int (*func)(void*,int,double*,int,double*));
		extern int ML_Set_Amatrix_Diag(ML*,int,int, double *);
		extern int ML_Set_Amatrix_Getrow(ML *ml, int level,
		int (*getrow)(void*,int,int*,int,int*,double*,int*),
		int (*comm  )(double *vec, void *data), int comm_vec_leng);

extern int ML_Set_Amatrix_GetrowNeighbors(ML*,int,int N_neigh,int *nlist);
extern int ML_Set_Amatrix_GetrowCommInfo(ML *, int level, int neighbor,
		int N_rcv, int *rcv_list, int N_send, int *send_list);
extern int ML_Set_Amatrix_NormalizationFactors(ML*,int,int n,double *fact);
extern int ML_Set_Amatrix_NullSpace(ML *, int, int, int, double *);

extern void ML_setup_grid_xsfer_op( void *f_grid, ML_GridFunc *fgrid_fcns,
		void *c_grid, ML_GridFunc *cgrid_fcns, void **xsfer,
		ML_Comm *comm );

extern int ML_Init_Restrictor(ML*,int L1,int L2,int,int,void *data);
extern int ML_Set_Restrictor_Matvec(ML*,int,
                int (*func)(void*,int,double*,int,double*));
		extern int ML_Set_Restrictor_Getrow(ML *ml, int level,
		int (*getrow)(void*,int,int*,int,int*,double*,int*),
		int (*comm  )(double *vec, void *data), int comm_vec_leng);
extern int ML_Set_Restrictor_GetrowNeighbors(ML *ml,int level,int N_neigh,
		int *neigh_list);
extern int ML_Set_Restrictor_GetrowCommInfo(ML *ml,int level,int neighbor,
		int N_rcv, int *rcv_list, int N_send, int *send_list);

extern int ML_Init_Prolongator(ML*,int L1,int L2,int,int,void *data);
extern int ML_Set_Prolongator_Matvec(ML *ml, int level,
		int (*func) (void *, int, double *, int, double *));
extern int ML_Set_Prolongator_Getrow(ML *ml, int level,
		int (*getrow)(void*,int,int*,int,int*,double*,int*),
		int (*comm  )(double *vec, void *data), int comm_vec_leng);
extern int ML_Set_Prolongator_GetrowNeighbors(ML *ml,int level,int N_neigh,
		int *neigh_list);
extern int ML_Set_Prolongator_GetrowCommInfo(ML *ml,int level,int neighbor,
		int N_rcv, int *rcv_list, int N_send, int *send_list);

extern int ML_Gen_CoarseSolverSuperLU(ML *ml_handle, int level);
extern int ML_Gen_SmootherJacobi( ML *, int nl, int pre_or_post,
		int ntimes, double omega );
extern int ML_Gen_SmootherGaussSeidel(ML*,int nl,int pre_post,int ntimes);
extern int ML_Gen_SmootherSymGaussSeidel(ML*,int nl,int pre_post,int ntimes,
		double omega);
extern int ML_Gen_BGSInverses(ML *ml,int nl,int blocksize,ML_Sm_BGS_Data **data);
extern int ML_Gen_SmootherBlockGaussSeidel(ML*,int nl,int pre_post,int ntimes,
		double omega, int blocksize);
extern int ML_Set_Smoother(ML *, int nl , int pre_post, void *obj,
		int (*func)(void *, int, double *, int, double *));
extern int ML_Set_CoarseSolver(ML *ml, int level, int leng, void *sol_obj,
                void (*solve)());

extern int ML_Gen_AmatrixRAP(ML *ml, int to_level, int from_level);

extern int ML_Set_EqnToridMapFunc(ML *ml, int level,int fleng,int tleng,
		void (*func)() );
extern int ML_Set_GridToqnMapFunc(ML *ml, int level,int fleng,int tleng,
		void (*func)() );
extern int ML_Set_BoundaryTypes(ML*,int level,int type,int n,int *data);

extern int ML_Gen_Solver(ML *ml, int method, int finest_level, int);
extern int ML_Iterate(ML *ml, double *sol, double *rhs);
extern int ML_Solve_MGV( ML *ml , double *din, double *dout);
extern double ML_Cycle_MGV(ML_1Level *curr, double *sol, double *rhs,
                int approx_all_zeros, ML_Comm *comm);
*/
#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif
#endif
