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

#ifdef WKC
#include <Epetra_MultiVector.h>
#endif

/* -------------------------------------------------------------------- */
/*! data structure used to store pointers to functions such as matvec   
    used by the operator class.                                         */
/* -------------------------------------------------------------------- */

struct ML_Function_Struct {
   int ML_id;
   int Nrows;
   int (*internal)(ML_Operator *, int, double *, int, double *);
   int (*external)(void *, int, double *, int, double *);
};

/* -------------------------------------------------------------------- */
/*! This data structure stores all information pertaining to performing  
    the Getrow function on an operator object.                          */
/* -------------------------------------------------------------------- */

struct ML_GetrowFunc_Struct {
   int           ML_id; 
   int           Nrows;
   ML_CommInfoOP *pre_comm;
   ML_CommInfoOP *post_comm;
   int           (*internal)(ML_Operator *,int,int*,int,int*,double*,int*);
   int           (*external)(void *,int,int*,int,int*,double*,int*);
   void          *data;
   int           use_loc_glob_map;
   int           *loc_glob_map;
   int           *row_map;
};

/* -------------------------------------------------------------------- */
/*! This data structure stores all information necessary to be able to   
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

/* -------------------------------------------------------------------- */
/*! This data structure defines an enriched operator class for the       
    specification of the discretization matrix, the restriction and the  
    prolongation operator.                                              */
/* -------------------------------------------------------------------- */

struct ML_Operator_Struct {
   int           ML_id;
   ML_Comm       *comm;
   ML_1Level     *to, *from;
   int           invec_leng, outvec_leng;
   void          *data;
   void          (*data_destroy)(void *);
   ML_Function   *matvec;
   ML_GetrowFunc *getrow;
   ML_DVector    *diagonal;      /** diagonal of matrix.     */
   int           N_nonzeros;
   int           max_nz_per_row;
   int           from_an_ml_operator;
   ML_Operator   *sub_matrix;
   ML_BdryPts    *bc;
   double        build_time, apply_time;
   int           ntimes, nflop;
   char          *label; 
   int           num_PDEs, num_rigid;
   double        lambda_max, lambda_min, lambda_max_img;
   int           N_total_cols_est;
   int           halfclone;
   ML_Operator_Subspace *subspace;
                /* This is just a hook into modes that we want to project out
                   before (after) invoking a MG cycle.  I couldn't think of
                   a more appropriate spot for these, especially as they need
                   to be available when ML is used as a preconditioner to a
                   Krylov method. */
};

/* -------------------------------------------------------------------- */
/*! This structure is used to implement both drop tolerances and matrix  
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

extern int ML_Operator_BlockPartition(ML_Operator *matrix, int nLocalNd, 
        int *nblk, int *pnode_part, int *ndwts /*=NULL*/, int *egwts/*=NULL*/, int nedges /*= 0*/ , int );
extern ML_Operator *ML_Operator_Create(ML_Comm *comm);
extern int ML_Operator_Destroy(ML_Operator **);

extern ML_Operator *ML_Operator_halfClone( ML_Operator *original);
extern int ML_Operator_halfClone_Init(ML_Operator *mat,
					   ML_Operator *original);

extern int ML_Operator_halfClone_Clean( ML_Operator *mat);
extern int ML_Operator_halfClone_Destroy( ML_Operator *mat);


extern int ML_Operator_Init(ML_Operator *, ML_Comm *comm);
extern int ML_Operator_Clean(ML_Operator *);
extern int ML_Operator_Dump(ML_Operator *Ke, double *, double *,
			    char *str, int);
extern int ML_hash_it( int new_val, int hash_list[], int hash_length, 
		       int *hash_used);

extern int ML_Operator_Set_Label(ML_Operator *, char *str);
extern int ML_Operator_Set_1Levels(ML_Operator *, ML_1Level*, ML_1Level*);
extern int ML_Operator_Set_BdryPts(ML_Operator *, ML_BdryPts *);
extern int ML_Operator_Set_ApplyFuncData(ML_Operator *, int, int, int, void*,
                      int,int (*func)(ML_Operator *,int,double*,int,double*),int);
extern int ML_Operator_Set_ApplyFunc(ML_Operator *, int, 
                       int (*func)(ML_Operator *, int, double *, int, double *));
extern int ML_Operator_Set_Diag(ML_Operator *, int, double *);
extern int ML_Operator_Set_Getrow(ML_Operator *, int, int, 
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
extern int ML_Operator_Print(ML_Operator *matrix, const char label[]);

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

#ifdef ML_WITH_EPETRA
extern int ML_Epetra_CRSinsert(ML_Operator *, int, int *, double *, int);
#endif

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif
