/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* ******************************************************************** */
/* Declaration of the ML structure                                      */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : March, 1999                                          */
/* ******************************************************************** */

#ifndef __MLSTRUCT__
#define __MLSTRUCT__

/* ******************************************************************** */
/* data structure type definition                                       */
/* ******************************************************************** */

typedef struct ML_Struct ML;

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
#include "ml_krylov.h"
#include "ml_mat_formats.h"
#include "ml_amg.h"
#include "ml_aggregate.h"
#include <string.h>

#ifdef WKC
/* WKC -- added header(s) for the new datastructures */
#include <Epetra_MultiVector.h>
#include <Epetra_LocalMap.h>
#endif

/* ******************************************************************** */
/* ******************************************************************** */
/* data definition for the ML Class                                     */
/* ******************************************************************** */
/* -------------------------------------------------------------------- */
/* This is the primary data structure which users interact directly     */
/* with via the provided functions.                                     */
/* -------------------------------------------------------------------- */

struct ML_Struct {
   int            id;
   int            ML_init_flag;    /* indicate initialization done */
   int            ML_scheme;       /* which ML scheme to pick      */
   int            ML_num_levels;   /* number of levels available   */
   int            ML_num_actual_levels;
                                   /* number of levels actually used */
                                   /* by the multigrid method.       */
   int            ML_num_transfers;/* number of transfers  */
   int            ML_finest_level, ML_coarsest_level;
   int            symmetrize_matrix;
   int            output_level;
   int            res_output_freq;
   double         tolerance;
   int            max_iterations;
   double         *spectral_radius;
   ML_Smoother    *pre_smoother;
   ML_Smoother    *post_smoother;
   ML_CSolve      *csolve;
   ML_Operator    *Amat, *Rmat, *Pmat;
   ML_Grid        *Grid;
   ML_BdryPts     *BCs;
   ML_Mapper      *eqn2grid;
   ML_Mapper      *grid2eqn;
   ML_1Level      *SingleLevel;
   ML_DVector     *Amat_Normalization;
   struct ML_Timing
                  *timing;       /* Used for timing information.    */
   ML_Comm        *comm;         /* communicator for ML             */
   int            *int_options;  /* optional integer parameters     */
   double         *dble_options; /* optional double parameters      */
   void           *void_options; /* optional other parameters       */
   int            (*func)(void);     /* optional function               */

};
struct ML_Timing {
   double precond_apply_time;
   double         total_build_time;
};

/* ******************************************************************** *
 * Control structure for the amount of information that ML prints.      *
 * ******************************************************************** */

typedef struct ML_PrintControl_Struct ML_PrintControl;

struct ML_PrintControl_Struct {
   int            output_level;
};

extern ML_PrintControl ML_PrintLevel;

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

extern int ML_Create(ML **ml, int Nlevels);
extern int ML_build_ggb(ML *ml, void *data);
extern int ML_Set_Symmetrize(ML *ml, int true_or_false);
extern int ML_Set_OutputLevel(ML *ml, int output_level);
extern int ML_Set_PrintLevel(int);
extern int ML_Get_PrintLevel(void);
extern int ML_Set_ResidualOutputFrequency(ML *ml, int output_freq);
extern int ML_Set_Tolerance(ML *ml, double tolerance);
extern int ML_Set_MaxIterations(ML *ml, int iterations);
extern int ML_Print_Timing(ML *ml);

extern int ML_Destroy(ML **ml);
#ifdef GREG
extern int ML_Destroy2(ML **ml);
#endif
extern void ML_Solve_SmootherDestroy(void *data);

extern int ML_Init_Comm(ML *ml);
extern int ML_Set_Comm_MyRank(ML *ml, int myrank);
extern int ML_Set_Comm_Nprocs(ML *ml, int nprocs);
extern int ML_Set_Comm_Communicator(ML *ml, USR_COMM com);
extern int ML_Set_Comm_Send(ML *ml, int (*send)(void*,unsigned int,int,int,USR_COMM));
extern int ML_Set_Comm_Recv(ML *ml, int (*recv)(void*,unsigned int,int*,int*,USR_COMM,USR_REQ*));
extern int ML_Set_Comm_Wait(ML *ml, int (*wait)(void*,unsigned int,int*,int*,USR_COMM,USR_REQ*));

extern int ML_Set_Comm(ML *ml, ML_Comm *comm);

extern int ML_Init_Grid(ML *, int nl, void *grid);
extern int ML_Set_Grid_GridFunc(ML *, int nl, ML_GridFunc *);
extern int ML_Set_Grid_MaxVertPerElmnt(ML *, int, int nvert);
extern int ML_Set_Grid_GetDimension(ML *, int nl, int (*func)(void *));
extern int ML_Set_Grid_GetNVert(ML *, int nl, int (*func)(void *));
extern int ML_Set_Grid_GetNElmnt(ML *, int nl, int (*func)(void *));
extern int ML_Set_Grid_GetElmntNVert(ML *, int nl, int (*func)(void *, int));
extern int ML_Set_Grid_GetElmntVertList(ML *, int nl, int (*func)(void *, int, int *));
extern int ML_Set_Grid_GetVertGlobalNum(ML *, int nl, int (*func)(void *, int));
extern int ML_Set_Grid_GetElmntGlobalNum(ML *, int nl, ml_big_int (*func)(void *, int));
extern int ML_Set_Grid_GetVertCoordinate(ML *, int nl, int (*func)(void *, int, double *));
extern int ML_Set_Grid_ComputeBasisCoef(ML *, int nl, int (*func)(void*,int,double*,int,double*,int*));
extern int ML_Set_Grid_GetElmntVolume(ML *, int nl, int (*func)(void*,int,int*,double*));
extern int ML_Set_Grid_GetElmntMatrix(ML *, int nl, int (*func)(void*,int,double**));
extern int ML_Set_Grid_GetElmntNullSpace(ML *, int, int (*func)(void*,int,double*));


extern int ML_Gen_GridXsferUsingFEBasis(ML *, int L1, int L2, int stride);
extern int ML_Gen_MGHierarchyVanek(ML *, int start, int increment_or_decrement);

extern int ML_Set_Grid(ML *, int nl, void *grid, ML_GridFunc *);

extern int ML_Init_Amatrix(ML *,int level,int ilen,int olen,void *data);
extern int ML_Get_Amatrix(ML *ml, int level, ML_Operator **matrix);
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

extern void ML_setup_grid_xsfer_op(void *, ML_GridFunc *, void *,
                     ML_GridFunc *, void **xsfer, ML_Comm *);

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

extern int ML_Gen_Blocks_Metis(ML *ml, int level, int *nblocks,int **block_list);

extern int ML_Gen_Smoother_Jacobi( ML *, int nl, int pre_or_post,
                     int ntimes, double omega );
extern int ML_Gen_Smoother_GaussSeidel(ML*,int nl,int pre_post,int ntimes,double);
extern int ML_Gen_Smoother_Hiptmair(ML*,int nl,int pre_post,int ntimes,
                     ML_Operator**, ML_Operator**, ML_Operator*,
				    void *, void **,
				    void *, void **, int);
extern int ML_Gen_Smoother_BlockHiptmair(ML*,int nl,int pre_post,int ntimes,
                     ML_Operator**, ML_Operator**, ML_Operator*,
				    void *, void **,
				    void *, void **, int);
extern int ML_Gen_Smoother_SymGaussSeidel(ML*,int nl,int pre_post,int ntimes,
		     double omega);
extern int ML_Gen_Smoother_SymGaussSeidelSequential(ML*,int nl,int pre_post,
                     int ntimes, double omega);
extern int ML_Gen_Smoother_MLS(ML*,int nl,int pre_post, double eig,
                     int degree);
extern int ML_Gen_Smoother_BlockDiagScaledCheby(ML *ml, int nl, 
					int pre_or_post,
					double eig_ratio, int deg,
					int nBlocks, int *blockIndices);

extern int ML_Gen_Smoother_OrderedSymGaussSeidel(ML *ml , int nl, int pre_or_post,
                     int ntimes, double omega);

extern int ML_Gen_Smoother_ParaSails(ML *ml, int nl, int pre_or_post, int ntimes,
   int sym, double thresh, int num_levels, double filter, int , int);

extern int ML_Gen_Smoother_BlockGaussSeidel(ML*,int nl,int pre_post,int ntimes,
		     double omega, int blocksize);
extern void BGS_Clean(void *data);
extern int ML_Gen_Smoother_VBlockJacobi(ML*,int nl,int pre_post, int ntimes,
                     double omeg, int Nblocks, int *blockList);
extern int ML_Gen_Smoother_VBlockSymGaussSeidel(ML*,int nl,int pre_post,
                     int ntimes, double omega, int Nblocks, int *blockList);
extern int ML_Gen_Smoother_VBlockSymGaussSeidelSequential(ML*,int nl, int,
                     int ntimes,double omega,int Nblocks,int *blockList);
extern int ML_Gen_Smoother_VBlockKrylovJacobi(ML*,int nl,int pre_post, int ntimes,
                     double omeg, int Nblocks, int *blockList);
extern int ML_Gen_Smoother_OverlappedDDILUT(ML*,int nl,int pre_post);
extern int ML_Gen_Smoother_VBlockAdditiveSchwarz(ML *,int nl,int pre_or_post,
                     int ntimes, int length, int *blkinfo);
extern int ML_Gen_Smoother_VBlockMultiplicativeSchwarz(ML *,int nl,
                     int pre_or_post, int ntimes, int length, int *blkinfo);

extern int ML_Gen_Smoother_GSextra( ML *ml , int nl, int pre_or_post,
		     int ntimes, double omega, int Nextra, int extra[]);
extern int ML_Set_Smoother(ML *, int nl , int pre_post, void *obj,
                     int (*func)(void *, int, double *, int, double *),
                     char *);

extern int ML_Gen_CoarseSolverSuperLU(ML *ml_handle, int level);
extern int ML_Gen_CoarseSolverAggregation(ML *ml_handle, int level,
                                          ML_Aggregate *ag);

extern int ML_Gen_AmatrixRAP(ML *ml, int to_level, int from_level);
extern int ML_Gen_Amatrix_Global(ML_Matrix_DCSR *inmat,
     ML_Matrix_DCSR *outmat, ML_Comm *comm, int *offset);

extern int ML_Set_EqnToGridMapFunc(ML *, int,int fleng,int tleng,
                                   int (*func)(void*,double*,double*));
extern int ML_Set_GridToEqnMapFunc(ML *, int,int fleng,int tleng,
                                   int (*func)(void*,double*,double*));
extern int ML_Set_BoundaryTypes(ML*,int level,int type,int n,int *data);

extern int ML_Setup(ML *ml, int method, int finest_level, int, void *);

extern int ML_Gen_Solver(ML *ml, int method, int finest_level, int);
extern int ML_Iterate(ML *ml, double *sol, double *rhs);
extern int ML_Solve(ML *ml, int inlen, double *sol, int outlen, double *rhs);

int ML_Solve_MGV( ML *ml , double *din, double *dout);

#ifdef WKC
/* WKC -- new prototype for V-cycle solve */
int ML_Solve_MGV( ML *ml , const Epetra_MultiVector &in ,
                         Epetra_MultiVector &out );
#endif

extern int ML_Solve_MGFull( ML *ml , double *din, double *dout);
extern int ML_Solve_Smoother(void *data, int isize, double *x, int osize,
			     double *rhs);

extern double ML_Cycle_MG(ML_1Level *curr, double *sol, double *rhs,
                     int approx_all_zeros, ML_Comm *comm, int, ML *ml);

#ifdef WKC
/* WKC -- new prototype for V-cycle solve */
extern double ML_Cycle_MG(ML_1Level *curr, Epetra_MultiVector &ep_sol,
                     Epetra_MultiVector &ep_rhs,
                     int approx_all_zeros, ML_Comm *comm, int, ML *ml);
#endif

extern double ML_Cycle_MGFull(ML_1Level *curr, double *sol, double *rhs,
                     int approx_all_zeros, ML_Comm *comm, int, ML *ml);
extern int ML_Solve_AMGV( ML *ml , double *din, double *dout);
extern double ML_Cycle_AMGV(ML_1Level *curr, double *sol, double *rhs,
                     int approx_all_zeros, ML_Comm *comm);
extern int ML_Solve_ProjectedAMGV( ML *ml , double *din, double *dout);
extern int ML_Gen_SmootherGSextra( ML *ml , int nl, int pre_or_post,
				   int ntimes, double omega, int Nextra,
				   int extra[]);
extern int ML_MLS_Setup_Coef(void *sm, int deg, int symmetrize);
extern int ML_Seg_Solve( ML *ml , double *din, double *dout);
extern int ML_Clean_CSolveSuperLU( void *vsolver, ML_CSolveFunc *func);
extern int ML_Solver_SetScheme(ML *ml, int scheme);
extern int ML_Smoother_Reinit(ML *ml);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif



