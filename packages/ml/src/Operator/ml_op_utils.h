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

#include "ml_common.h"

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

extern int ML_Operator_ChangeToSinglePrecision(ML_Operator *matrix);
extern int ML_Operator_ChangeToChar(ML_Operator *matrix);
extern int ML_Operator_ImplicitTranspose(ML_Operator *Rmat, 
					 ML_Operator *Pmat,
					 int PostCommAlreadySet);
extern int ML_Gen_Restrictor_TransP(ML *, int, int);
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

#ifndef ML_CPP
#ifdef __cplusplus
  }
#endif
#endif

#endif
