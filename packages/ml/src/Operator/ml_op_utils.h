/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* ******************************************************************** */
/* Declaration of the New stuff                                         */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : March, 1999                                          */
/* ******************************************************************** */

#ifndef __MLMATRIX__
#define __MLMATRIX__

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include "ml_common.h"

/************************************************************************/
/*       Structures used to implement block matrices where             */
/*       matrix[i+j*NBlockRows] is the (i,j)th block and NULL          */
/*       entries are interpreted as blocks of zeros.                   */
/***********************************************************************/

#define ML_DESTROY_EVERYTHING  2  /* destroy everything including      */
                                  /* entries in matrix                 */
                                  /* array  when destructor is invoked */
#define ML_DESTROY_SHALLOW     3  /* destroy everything but entries in */
                                  /* matrix                            */
#define ML_CLEAN_EVERYTHING    4  /* destroy everything but entries in */
                                  /* matrix which are instead cleaned. */

typedef struct ML_BlkMatData_Struct ML_BlkMatData;

/* Used to record local and global ids associated with each block matrix*/
struct ML_BlkMatData_Struct {
   int *GlobalId;
   int *LocalId;
};

struct MLBlkMat {
   int NBlockRows;        /* number of block rows in matrix               */
   int NBlockCols;        /* number of block cols in matrix               */
   int *RowStart;         /* RowStart[i] is row # of 1st row in blk row i */
   int *ColStart;         /* ColStart[i] is col # of 1st col in blk col i */
   int NGhost;            /* total number of ghosts associated with a     */
                          /* matvec on the entire matrix.                 */
   int invec;             /* total number of assigned columns             */
   int outvec;            /* total number of assigned rows                */
   ML_Operator **matrix;  /* matrix[i+j*NBlockRows] is (i,j)th block.     */
   ML_BlkMatData **matdata;/*matdata[i+j*NBlockRows] gives data for  */
                          /* (i,j)th block.                               */
   int destroy_level;     /* indicates whether *matrix should be          */
                          /* destroyed and/or matrix[k]'s should be       */
                          /* destroyed/clean when destructor invoked.     */
   int final_called;      /* indicates if ML_Operator_BlkMatFinalize()    */
                          /* has been invoked or not.                     */
};



#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif


extern int oldML_Mdfy_Prolongator_DirBdry(ML *, int , double *, double *);
extern int ML_Compute_Coarse_Bdry(ML *ml_handle, int level, int size,
           int fine_size);
extern int ML_Mdfy_Prolongator_DirBdry(ML *ml_handle, int level2, int size,
     int fine_size );

extern ML_Operator *ML_Operator_ExplicitlyScale(ML_Operator *matrix,
					 double scalar);

extern int ML_Operator_ChangeToSinglePrecision(ML_Operator *matrix);
extern int ML_Operator_ChangeToChar(ML_Operator *matrix);
extern int ML_Operator_ImplicitTranspose(ML_Operator *Rmat,
					 ML_Operator *Pmat,
					 int PostCommAlreadySet);
extern int ML_Gen_Restrictor_TransP(ML *, int, int, ML_Operator*);
extern int ML_Gen_Prolongator_Getrow(ML *, int , int , int , int ,
            int (*)(void* , int , int *, int , int *, double *, int *),
            int (*)(double *, void*), void *data, int);
  /*
extern int ML_Operator_Transpose(ML_Operator *Amat, ML_Operator *Amat_trans );
  */
extern int ML_Operator_ColPartition2RowPartition(ML_Operator *A, ML_Operator *Atrans);
extern int eye_getrows(ML_Operator *data, int N_requested_rows, int requested_rows[],
		       int allocated_space, int columns[], double values[],
		       int row_lengths[]);
extern	int eye_matvec(ML_Operator *Amat_in, int ilen, double p[], int olen, double ap[]);
extern int ML_Operator_Transpose_byrow(ML_Operator *A, ML_Operator *Atrans);
extern int ML_Operator_Getrow_Diag(ML_Operator *Amat, double **diagonal);
extern int ML_build_overlapped_pre_comm(ML_Operator *tempA, ML_CommInfoOP
					*old_comm, int max_per_proc,
					int *hash_list, int hash_length,
					int *hash_used, int old_Nrows,
					int *Nexternal, int *external[],
					int *Nexternal_allocated);
extern int ML_Operator_HashGlobalRcvList(ML_CommInfoOP *pre_comm, int Nrows,
					 int hash_list[], int hash_length,
					 int *hash_used, int Gid_assigned_to_proc[],
					 ML_Comm *comm,
					 int *Nexternal, int **external,
					 int *Nexternal_allocated);
extern int ML_overlap(ML_Operator *oldA, ML_Operator *newA, int overlap,
		      ML_CommInfoOP **nonOverlapped_2_Overlapped);
extern void ML_Operator_ReportStatistics(ML_Operator *mat, char *appendlabel,
                                         int perfAndCommStats);
extern void ML_Operator_Profile(ML_Operator *A, char *appendlabel);
extern void ML_Operator_Profile_SetIterations(int numits);
extern int ML_Operator_Profile_GetIterations();
extern int ML_Operator_Get_Nnz(ML_Operator *A);
extern char* ML_Operator_IdentifyDirichletRows(ML_Operator *A);

extern int  ML_Operator_BlkMatInit(ML_Operator *BlkMat, ML_Comm *comm,
       int NBlockRows, int NBlockCols, int destroy_level);

extern void  ML_Operator_BlkMatDestroy(void *data);

extern int ML_Operator_BlkMatInsert(ML_Operator *BlkMat, ML_Operator *Entry,
                        int Row, int Col);

extern ML_Operator *ML_Operator_BlkMatExtract(ML_Operator *BlkMat,
                        int Row, int Col);
extern int ML_Operator_BlkMatNumBlockRows(ML_Operator * BlkMat);
extern int ML_Operator_BlkMatNumBlockCols(ML_Operator * BlkMat);
extern int  ML_Operator_BlkMatFinalize(ML_Operator *BlkMat);

extern int ML_Operator_BlkMatMatvec(ML_Operator *BlkMat, int ilen,
        double p[], int olen, double ap[]);

extern int ML_Operator_BlkMatGetrow(ML_Operator *Amat, int N_requested_rows,
   int requested_rows[], int allocated_space, int columns[],
   double values[], int row_lengths[]);

extern int ML_Blkrap(ML_Operator *Rmat, ML_Operator *Amat, ML_Operator *Pmat,
              ML_Operator *Result, int matrix_type);

extern ML_Operator *ProjectMe(ML *mlptr, int BlockLocation, int FromLevel, ML_Operator *Matrix, int matrix_type);

extern void ML_Operator_Copy_Statistics(ML_Operator *source, ML_Operator *target);

#ifndef ML_CPP
#ifdef __cplusplus
  }
#endif
#endif

#endif
