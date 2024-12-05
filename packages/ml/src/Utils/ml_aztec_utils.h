/********************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/********************************************************************* */
/*          Utilities for Aztec/ML users                               */
/********************************************************************* */

#ifndef __MLAZUTILS__
#define __MLAZUTILS__

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif


#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

#ifdef AZTEC
#if defined(ML_MPI) && !defined(AZ_MPI)
#define AZ_MPI
#endif

#define AZ_ONLY_PRECONDITIONER  -7778

#include "ml_common.h"
#include "az_aztec.h"

/******************************************************************************/
/* structure/function to wrap Aztec subblocks into a large Aztec matrix.      */
/******************************************************************************/
struct AZ_MAT_blockmat_data {
  int N, Nghost;
  AZ_MATRIX *Ke;
  AZ_MATRIX *M;
};



struct aztec_context {
   AZ_MATRIX *Amat;
   AZ_PRECOND *Prec;
   int    *proc_config;
   int    *options;
   double *params;
   double *status;
   void   *getrowstuff;
   int    prec_or_solver;
   ML_Comm *comm;
   int     offset;
   int    matrix_type;
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

extern int az_wrap_solvers(ML_Smoother *, int in, double x[], int out,
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

extern int AZ_block_MSR(int **param_bindx, double **param_val,
                 int N_update, int num_PDE_eqns, int *update);

#define mdwrap_wait(a,b,c,x,y,z)   md_wrap_wait(a,b,c,(x),(y),(z))
#define mdwrap_iwrite(a,b,c,x,y,z) md_wrap_iwrite(a,b,c,(x),(y),(z))
#define mdwrap_iread(a,b,c,x,y)   md_wrap_iread((a),(b),(c),(x),(y))
#define mdwrap_write(a,b,c,x,y)   md_wrap_write((a),(b),(c),(x),(y))

#if   defined(caps)
#   define az_set_proc_config_        AZ_SET_PROC_CONFIG
#else
#if defined(matched)
#   define az_set_proc_config_        az_set_proc_config
#endif
#endif

extern void AZ_set_proc_config(int proc_config[], MPI_AZComm );


#endif
extern void AZ_mlcomm2data_org(ML_CommInfoOP *comm_info,int *data_org[]);

extern void AZ_set_ML_preconditioner(AZ_PRECOND **Precond, AZ_MATRIX *Amat,
				ML *ml_handle, int options[]);


extern int az_comm_wrapper(double vector[], void *data);
extern int az_wrap_ml_comm(double vector[], AZ_MATRIX *Amat);

extern int az_msrgetrow_wrapper(ML_Operator *data, int N_requested_rows,
   int requested_rows[], int allocated_space, int columns[], double values[],
   int row_lengths[]);

extern int az_vbrgetrow_wrapper(ML_Operator *data, int N_requested_rows,
   int requested_rows[], int allocated_space, int columns[], double values[],
   int row_lengths[]);
extern int az_wrap_ml_getrow(int columns[], double values[], int row_lengths[],
		      struct AZ_MATRIX_STRUCT *Amat, int N_requested_rows,
			     int requested_rows[], int allocated_space);


extern int az_usergetrow_wrapper(ML_Operator *data, int N_requested_rows,
   int requested_rows[], int allocated_space, int columns[], double values[],
   int row_lengths[]);

extern int az_matvec_wrapper(ML_Operator *data, int, double p[], int, double ap[]);

extern void az_wrap_ml_matvec(double invec[], double outvec[], AZ_MATRIX *Amat,
			      int proc_config[]);
extern int AZ_convert_aztec_matrix_2ml_matrix(AZ_MATRIX *AZmat, ML_Operator *MLmat,
					      int *proc_config);
extern int AZ_ML_Set_Amat(ML *ml_handle, int level, int isize, int osize,
	AZ_MATRIX *Amat, int *proc_config);

extern void AZ_ML_set_vbrdiagonal(ML *ml, int level, AZ_MATRIX *matrix);

extern void AZ_ML_set_userdiagonal(ML *ml, int level, AZ_MATRIX *matrix);

extern void ML_precondition(double ff[], int options[], int proc_config[],
                     double params[], AZ_MATRIX *mat, AZ_PRECOND *prec);
extern void MLsmoother_precondition(double ff[], int options[],
				    int proc_config[], double params[],
				    AZ_MATRIX *mat, AZ_PRECOND *prec);

extern void AZ_ML_Clean(void *data);

extern int ML_MSR_sym_diagonal_scaling(AZ_MATRIX *Amat,
				       int proc_config[], double **);

extern int ML_MSR_scalesol(double *x, double *scale_vect,int length);
extern int ML_MSR_scalerhs(double *x, double *scale_vect,int length);
extern int AZ_block_MSR(int **param_bindx, double **param_val,
			int N_update, int num_PDE_eqns, int *update);


extern int  wrapper_DCSR_getrow(int columns[], double values[], int row_lengths[],
        	     struct AZ_MATRIX_STRUCT *Amat, int N_requested_rows,
                     int requested_rows[], int allocated_space);
extern void wrapper_DCSR_matvec(double *b, double *c,AZ_MATRIX *Amat,
                     int proc_config[]);
extern void AZ_transform_norowreordering(int proc_config[], int *external[],
    int bindx[], double val[], int update[], int *update_index[],
	int *extern_index[], int *data_org[], int N_update, int indx[], int bnptr[],
	int rnptr[], int *cnptr[], int mat_type);
extern void AZ_input_msr_matrix_nodiag(char datafile[], int update[],
	double **val, int **bindx,
	int N_update, int proc_config[]);
extern void AZ_add_new_row_nodiag(int therow, int *nz_ptr, int *current,
	double **val, int **bindx, char *input, FILE *dfp, int *msr_len,
	int *column0);
extern void ML_find_local_indices(int N_update, int bindx[], int update[],
	int *sorted_ext, int N_external, int map[], int start, int end);

extern void AZ_Tmat_transform2ml(int Nexterns, int global_node_externs[], int *reordered_node_externs,
			    int Tmat_bindx[], double Tmat_val[], int rowptr[], int Nlocal_nodes,
				 int global_node_inds[], ML_Comm *comm, int Nlocal_edges,
				 ML_Operator **Tmat);

extern void AZ_zeroDirichletcolumns(AZ_MATRIX *Amat, double rhs[],
				    int proc_config[] );
extern int ML_Tmat_applyDirichletBC(ML_Operator **Tmat, int *dirichlet_rows,
                             int num_dirichlet_rows);
extern void AZ_block_matvec(double *, double *, AZ_MATRIX *, int *);

extern int ML_Aggregate_AztecRead(ML_Aggregate *ag);
extern void new_norm(AZ_PRECOND *prec, double res[], double *result);
/* mgee */
extern void ML_AZ_Reader_ReadVariableBlocks(char *cmd_file_name, int *nblocks, int **blocks,
                                     int **block_pde, int *N_update, int **update,
                                     int proc_config[]);


/* The following definition and function declaration as used by
   MLAZ_iterate, which is supposed to replace AZ_iterate in code
   currently developed for Aztec (one-level prec), to test ML
   preconditioners (based on aggregation). Only a small part
   of ML's functionalities is currently supported. */

#define MLAZ_MAX_LEVELS 30
#define MLAZ_COARSE_LEVEL (MLAZ_MAX_LEVELS-1)
#define MLAZ_ALL -1
#define MLAZ_yes 1
#define MLAZ_no 0

/* options */
#define MLAZ_smoother                       1
#define MLAZ_max_levels                     2
#define MLAZ_num_smoother_steps             3
#define MLAZ_output                         7
#define MLAZ_coarsen_scheme                 8
#define MLAZ_metis_aggregation_property     9
#define MLAZ_metis_aggregation_value       10
#define MLAZ_max_coarse_size               11
#define MLAZ_is_problem_symmetric          12
#define MLAZ_pre_or_post_smoother          13
#define MLAZ_num_pde_eqns                  14
#define MLAZ_smoother_damping              15
#define MLAZ_amesos_solver                 16
#define MLAZ_max_procs                     17
#define MLAZ_req_aggre_per_proc            18
#define MLAZ_MLS_poly_order                19
#define MLAZ_MLS_alpha                     20
#define MLAZ_timing_detailed               21
#define MLAZ_prec_type                     22

/*  MLAZ_smoother */
#define MLAZ_Jacobi                 0 /* ML's Jacobi smoother */
#define MLAZ_GaussSeidel            1 /* ML's GS smoother */
#define MLAZ_SymGaussSeidel         2
#define MLAZ_MLS                    3
#define MLAZ_Aztec                  4
#define MLAZ_BlockGaussSeidel       5
#define MLAZ_IFPACK                 7
#define MLAZ_SuperLU                -1
#define MLAZ_Amesos                 -2

/* MLAZ_aggregates */
#define MLAZ_METIS                  0
#define MLAZ_ParMETIS               1
#define MLAZ_Uncoupled              2
#define MLAZ_MIS                    3
/* MLAZ_metis_aggregate_property */
#define MLAZ_NumLocalAggregates        0
#define MLAZ_NumGlobalAggregates       1
#define MLAZ_NumNodesPerAggregate 2

/* params */
#define MLAZ_smoothP_damping_factor 1
#define MLAZ_omega                  2
#define MLAZ_threshold              3
#define MLAZ_dumping_factor         4

extern void MLAZ_Defaults(void);
extern void MLAZ_Iterate( double delta_x[], double resid_vector[],
			  int options[], double params[],
			  double status[], int proc_config[],
			  AZ_MATRIX *Amat, struct AZ_SCALING *scaling );
extern void MLAZ_Set_Option(int, int);
extern void MLAZ_Set_Param( int,double);
extern void MLAZ_Set_LevelOption( int level,int option, int value);
extern void MLAZ_Set_LevelParam( int level, int option, double value);
extern void MLAZ_Set_LevelAztecSmoother( int level, int options[], double value[]);
extern int MLAZ_Setup_MLandAggregate( int N_update, int num_PDE_eqns,
			       int proc_config[AZ_PROC_SIZE],
				      ML *ml, ML_Aggregate *ag);
void MLAZ_Direct_Solve_Amesos( double delta_x[], double resid_vector[],
			       AZ_MATRIX * Amat, int proc_config[],
			       int choice, int max_procs );
void AZ_ML_Build_NodalCoordinates( int N, int N_update, int N_external,
				   int update[], int external[],
				   int update_index[], int extern_index[],
				   double x[], double y[], double z[],
				   int option );
#endif

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#ifdef AZTEC_RAY_WILL_FIX

#include "az_aztec.h"
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/
#endif
#endif
