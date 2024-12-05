/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* ******************************************************************** */
/* Declaration of the ML_Operator structure                             */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : March, 1999                                          */
/* ******************************************************************** */

#ifndef __MLOPERATOR__
#define __MLOPERATOR__

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/* ******************************************************************** */
/* data structure type definition                                       */
/* ******************************************************************** */

typedef struct ML_Operator_Subspace_Struct ML_Operator_Subspace;
typedef struct ML_Operator_Struct ML_Operator;
typedef struct ML_Function_Struct ML_Function;
typedef struct ML_GetrowFunc_Struct ML_GetrowFunc;

/* ******************************************************************** */
/* local include files                                                  */
/* ******************************************************************** */

#include "ml_common.h"
#include "ml_defs.h"
#include "ml_memory.h"
#include "ml_bdrypts.h"
#include "ml_1level.h"
#include "ml_operatoragx.h"
#include "ml_vec.h"
#include "ml_gridagx.h"
#include "ml_mls.h"
#include "ml_utils.h"
#include "ml_op_utils.h"

#ifdef WKC
#include <Epetra_Comm.h>
#include <Epetra_MultiVector.h>
#endif

/* -------------------------------------------------------------------- */
/*  data structure used to store pointers to functions such as matvec
    used by the operator class.                                         */
/* -------------------------------------------------------------------- */

struct ML_Function_Struct {
   int ML_id;
   int Nrows;
   int (*func_ptr)(ML_Operator *, int, double *, int, double *);
};

/* -------------------------------------------------------------------- */
/*  This data structure stores all information pertaining to performing
    the Getrow function on an operator object.                          */
/* -------------------------------------------------------------------- */

struct ML_GetrowFunc_Struct {
   int           ML_id;
   int           Nrows; /*used for all point rows even VBR point getrows*/
   int           N_block_rows; /*used to store number of block rows for VBR matrices*/
   ML_CommInfoOP *pre_comm; /*This is for all normal cases*/
   ML_CommInfoOP *post_comm; /*This is only used for the weird matvec/transpose*/
   int           (*func_ptr)(ML_Operator *,int,int*,int,int*,double*,int*);
   void          *data;
   int           use_loc_glob_map; /*If global when a point getrow is performed local indices are returned*/
   int           columns_loc_glob; /*Set to ML_LOCAL_INDICES if local or ML_GLOBAL_INDICES if global*/
   int           *loc_glob_map;
   int           *row_map;
};

/* -------------------------------------------------------------------- */
/*  This data structure stores all information necessary to be able to
    project out a subspace (e.g., a known nullspace).                   */
/* -------------------------------------------------------------------- */

struct ML_Operator_Subspace_Struct {
   double **basis_vectors;
   int    dimension;                /** number of basis vectors */
   int    vecleng;                  /** length of basis vectors */
   void   (*data_destroy)(void *);
   double *VAV;                     /** dimension by dimension system to solve */
   int    *pivots;                  /** pivots for VAV factorization */
   int    VAVdone;                  /** true if VAV is calculated already */
   double *res1,*res2,*vec1,*vec2;      /* work vectors */
};

typedef struct {
  double    threshold;
  int       (*aux_func_ptr)(ML_Operator *,int,int*,int,int*,double*,int*);
  int       enable;
  int       max_level;
  int**     filter;
  int       filter_size;
  double    m_threshold; /* material threshold */

  /* For tracking the filtered vs. unfiltered number of nonzeros */
  int filtered_nnz;
  int unfiltered_nnz;

} ML_Aux_Data;

void ML_Aux_Data_Create(ML_Aux_Data** ptr);

void ML_Aux_Data_Destroy(ML_Aux_Data** ptr);

ML_Aux_Data* ML_Aux_Data_Clone(ML_Aux_Data* original);

#define ML_TYPE_UNKNOWN 0
#define ML_TYPE_ROW_MATRIX 1
#define ML_TYPE_CRS_MATRIX 2
#define ML_TYPE_VBR_MATRIX 3

/* -------------------------------------------------------------------- */
/** This data structure defines an enriched operator class for the
    specification of the discretization matrix, the restriction and the
    prolongation operator.                                              */
/* -------------------------------------------------------------------- */

struct ML_Operator_Struct {
   int           ML_id; /*this and the getrow are not the pid that lives in comm
                           ML_id is used as the pid for communication in the
                           neighbor list however*/
   ML_Comm       *comm;
   ML_1Level     *to, *from;
   int           invec_leng, outvec_leng;
   void          *data;
   void          (*data_destroy)(void *);
   ML_Function   *matvec;
   ML_GetrowFunc *getrow;
   ML_DVector    *diagonal;      /** diagonal of matrix.     */
   int           N_nonzeros;
   int           max_nz_per_row; /* largest local stencil size */
   int           min_nz_per_row; /* minimum local stencil size */
   int           avg_nz_per_row; /* average local stencil size */
   int           blocks; /*only used for VBR matrices to say number of blocks*/
   int           from_an_ml_operator;
   ML_Operator   *sub_matrix;
   ML_BdryPts    *BCs;
   char          *DirichletRows; /* simple array of length outvec_leng
                                    to record Dirichlet rows */
   double        build_time, apply_time;
   double        apply_without_comm_time;
   int           ntimes, nflop;
   char          *label;
   int           num_PDEs, num_rigid;
   double        lambda_max, lambda_min, lambda_max_img;
   int           N_total_cols_est;
   int           halfclone;
   int           spectral_radius_scheme, spectral_radius_max_iters;
   int           NumZDir;
   int           Zorientation;      /* -1: not specified */
                                    /*  1: vertical      */
                                    /*  2: horizontal    */
   char          coarsencoord;    /* x,y,z for semicoarsening/line smoothing */
   int           sortColumnsAfterRAP; /* just for Paul Lin and a bit of an */
                                       /* ugly hack                        */

   ML_Operator_Subspace *subspace;
                /* This is just a hook into modes that we want to project out
                   before (after) invoking a MG cycle.  I couldn't think of
                   a more appropriate spot for these, especially as they need
                   to be available when ML is used as a preconditioner to a
                   Krylov method. */
   ML_Aux_Data   *aux_data;
                 /*!< General container for auxiliary matrix */
   int           type; /* simple ID that specifies the actual storage
                          used in data. It can be:
                          - ML_TYPE_UNKNOWN (default)
                          - ML_TYPE_ROW_MATRIX
                          - ML_TYPE_CRS_MATRIX
                          - ML_TYPE_VBR_MATRIX
                          By using this, we can same some wrapping, at least
                          for the finest-level operator.
                        */
};


/* -------------------------------------------------------------------- */
/*  This structure is used to implement both drop tolerances and matrix
    amalgamation (used in ML_aggregateCoarsenMIS()). The idea is to wrap
    the getrow() of the original matrix such that it handles the blocking
    and the dropping.                                                   */
/* -------------------------------------------------------------------- */

struct amalg_drop {
   void                 *original_data;
   struct ML_GetrowFunc_Struct *original_getrow;
   double               *scaled_diag;
   int                  block_size;
   double               drop_tolerance;
   ML_Operator          *Amat;
   int                  *blk_inds;
   /* used by ML_Operator_AmalgamateAndDropWeak_VBlocks */
   void                 *vblock_data;  /**< holds data structure aggr_vblock */
};

/* -------------------------------------------------------------------- */
/*  This structure is used to implicitly scale a matrix. The idea is to wrap
    the getrow() of the original matrix such that it handles the blocking
    and the dropping.                                                   */
/* -------------------------------------------------------------------- */

struct ml_matscale {
  ML_Operator *Amat;
  double      scalar;
  int         destroy_child;
};

/* -------------------------------------------------------------------- */
/*  This structure is used to implicitly scale a matrix with a vector.  *
 *  It extends ml_matscale.                                             */
/* -------------------------------------------------------------------- */

struct ml_matvscale {
  ML_Operator *Amat;
  double*     scale;
  int         destroy_child;
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

extern int ML_Operator_BlockPartition(ML_Operator *matrix, int n, int *nblks,
                          int *pnode_part, ML_Partitioner which_partitioner,
		          double *x_coord, double *y_coord, double *z_coord,int num_PDE_eqns);

extern ML_Operator *ML_Operator_Create(ML_Comm *comm);
extern int ML_Operator_Destroy(ML_Operator **);

extern ML_Operator *ML_Operator_halfClone( ML_Operator *original);
extern int ML_Operator_halfClone_Init(ML_Operator *mat,
					   ML_Operator *original);

extern int ML_Operator_halfClone_Clean( ML_Operator *mat);
extern int ML_Operator_halfClone_Destroy( ML_Operator **mat);


extern int ML_Operator_Init(ML_Operator *, ML_Comm *comm);
extern int ML_Operator_Clean(ML_Operator *);
extern int ML_Operator_Dump(ML_Operator *Ke, double *, double *,
			    char *str, int);

extern int ML_Operator_Set_Label(ML_Operator *, char *str);
extern int ML_Operator_Set_1Levels(ML_Operator *, ML_1Level*, ML_1Level*);
extern int ML_Operator_Set_BdryPts(ML_Operator *, ML_BdryPts *);
extern int ML_Operator_Set_ApplyFuncData(ML_Operator *, int, int, void*,
                      int,int (*func)(ML_Operator *,int,double*,int,double*),int);
extern int ML_Operator_Set_ApplyFunc(ML_Operator *,
                       int (*func)(ML_Operator *, int, double *, int, double *));
extern int ML_Operator_Set_Diag(ML_Operator *, int, double *);
extern int ML_Operator_Set_Getrow(ML_Operator *, int,
                       int (*func)(ML_Operator *,int,int*,int,int*,double*,int*));

extern int ML_Operator_Getrow(ML_Operator *, int, int *, int, int *,
                              double *, int*);
extern int ML_Operator_Get_Diag(ML_Operator *Amat, int length, double **diag);

extern int ML_Operator_Apply(ML_Operator *, int, double *, int, double *);

#ifdef WKC
/* WKC -- ADDED HEADER */
extern int ML_Operator_Apply(ML_Operator *, int, Epetra_MultiVector &,
                             int, Epetra_MultiVector & );
#endif


extern int ML_Operator_ApplyAndResetBdryPts(ML_Operator *, int, double *,
                                            int olen, double *);
extern int ML_Operator_Add(ML_Operator *A, ML_Operator *B, ML_Operator *C,
			   int matrix_type, double scalar);
#ifdef WKC
/* WKC -- ADDED HEADER */
extern int ML_Operator_ApplyAndResetBdryPts(ML_Operator *, int,
                     Epetra_MultiVector &, int olen, Epetra_MultiVector &);
#endif
extern int ML_Operator_MoveFromHierarchyAndClean(ML_Operator *newmat,
						 ML_Operator *hier);

extern int ML_Operator_Move2HierarchyAndDestroy(ML_Operator **newmat,
						ML_Operator *hier);

extern int ML_Operator_Transpose(ML_Operator *Amat, ML_Operator *Amat_trans );

extern int ML_Operator_Check_Getrow(ML_Operator *, int, char*);
extern double ML_Operator_MaxNorm(ML_Operator *matrix, int divide_diag);
extern double ML_Operator_FroNorm(ML_Operator *matrix, int divide_diag);
extern int ML_Operator_Print(ML_Operator *matrix, const char label[]);
extern int ML_Operator_ComputeNumNzs(ML_Operator *matrix);
/* Operator Scaling */
extern int ML_implicitscale_Getrow(ML_Operator *data, int N_requested_rows,
				   int requested_rows[], int allocated_space,
				   int columns[], double values[],
				   int row_lengths[]);
extern int ML_implicitscale_Matvec(ML_Operator *Amat_in, int ilen, double p[],
				   int olen, double ap[]);
extern ML_Operator *ML_Operator_ImplicitlyScale(ML_Operator *Amat,
						double scalar,
						int OnDestroy_FreeChild);
extern void ML_implicitscale_Destroy(void *data);
/* Operator Scaling with a vector */
extern int ML_implicitvscale_Getrow(ML_Operator *data, int N_requested_rows,
				   int requested_rows[], int allocated_space,
				   int columns[], double values[],
				   int row_lengths[]);
extern int ML_implicitvscale_Matvec(ML_Operator *Amat_in, int ilen, double p[],
				   int olen, double ap[]);
extern ML_Operator *ML_Operator_ImplicitlyVScale(ML_Operator *Amat,
                                                 double* scale,
                                                 int OnDestroy_FreeChild);
extern int ML_Operator_ExplicitDinvA(int BlockSize,
				      struct MLSthing *Dinv, ML_Operator *A);

extern ML_Operator *ML_Operator_ImplicitlyBlockDinvScale(ML_Operator *Amat);
extern void ML_implicitvscale_Destroy(void *data);
extern int ML_implicitvcscale_Getrow(ML_Operator *data, int N_requested_rows,
				   int requested_rows[], int allocated_space,
				   int columns[], double values[],
				   int row_lengths[]);
extern ML_Operator *ML_Operator_ImplicitlyVCScale(ML_Operator *Amat,
                                                 double* scale,
                                                 int OnDestroy_FreeChild);
extern int ML_CSR_DropSmall(ML_Operator *Pe, double AbsoluteDrop,
			    double RelativeRowDrop, double RelativeColDrop);
extern ML_Operator *ML_CSRmatrix_ColumnSubset(ML_Operator *Amat, int Nsubset,
					      int subset[]);
extern int ML_Operator_Set_SpectralNormScheme_Calc(       ML_Operator *mat);
extern int ML_Operator_Set_SpectralNormScheme_Anorm(      ML_Operator *mat);
extern int ML_Operator_Set_SpectralNormScheme_Anasazi(    ML_Operator *mat);
extern int ML_Operator_Set_SpectralNormScheme_PowerMethod(ML_Operator *mat);
extern int ML_Operator_Set_SpectralNorm_Iterations(ML_Operator *mat, int its);



/* amalagamation routines */
extern int ML_Operator_AmalgamateAndDropWeak(ML_Operator *Amat, int block_size,
               double drop_tolerance);

extern int ML_Operator_UnAmalgamateAndDropWeak(ML_Operator *Amat,
		int block_size, double drop_tolerance);

extern int ML_amalg_drop_getrow(ML_Operator *data, int N_requested_rows,
		int requested_rows[], int allocated_space, int columns[],
                double values[], int row_lengths[]);

extern int ML_Operator_GetDistributedDiagBlocks(ML_Operator *mat, int *blkinfo,
                                                int **new_ja, double **new_aa);

extern double ML_Operator_GetMaxEig(ML_Operator *Amat);

extern ML_Operator **ML_Operator_ArrayCreate( int length);
extern int ML_Operator_ArrayDestroy( ML_Operator **array, int length);
extern int ML_Operator_SetSubspace(ML *ml, double **vectors, int numvecs,
                                   int vecleng);
extern int ML_Operator_Amalgamate_Vec_Trans(ML_Operator *Amat, int *blocked,
                                            int **unblocked, int *size);
extern int AZ_get_MSR_arrays(ML_Operator *, int **bindx, double **val);
     /* define this here so we don't have to include ml_aztec_utils.h */
     /* in ml_struct.c and ml_smoother.c                              */
int ML_Operator_GetFlops(ML_Operator *mat);
void ML_Operator_GetGlobalDimensions(ML_Operator *A,int *nrows,int *ncols);

extern int ML_Operator_MisRootPts( ML_Operator *Amatrix,  int num_PDE_eqns,
				     int **);

#ifdef ML_WITH_EPETRA
extern int ML_Epetra_CRSinsert(ML_Operator *, int, int *, double *, int);
#endif

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif
