#ifndef _ML_EPETRA_UTILS_H_
#define _ML_EPETRA_UTILS_H_

class Epetra_Comm;
class Epetra_BlockMap;
class Epetra_MultiVector;
class Epetra_RowMatrix;
class Epetra_Map;
class Epetra_Vector;
class Epetra_Import;
class Epetra_Object;
class Epetra_CrsMatrix;
class Epetra_RowMatrix;
class Epetra_LinearProblem;

#include "ml_common.h"

#ifdef ML_MPI
#ifndef EPETRA_MPI
#define EPETRA_MPI
#endif
#include "mpi.h"
#endif
#include "ml_include.h"

//! matvec function for Epetra_RowMatrix objects
int Epetra_ML_matvec(ML_Operator *data, int in, double *p, int out,
                 double *ap);

#ifdef WKC
int Epetra_ML_matvec_WKC(ML_Operator *data, int in, double *p, int out,
                 double *ap);
#endif 

//! Getrow function for Epetra_RowMatrix
/*! GetRow function for matrix of type Epetra_RowMatrix.  Supply local
  matrix (without ghost node columns) for rows given by requested_rows[0
  ... N_requested_rows-1].  Return this information in 'row_lengths,
  columns, values'.  If there is not enough space to complete this
  operation, return 0. Otherwise, return 1.

  \param \in data: pointer to user's data containing matrix
  values.

  \param \in N_requested_rows: number of rows for which nonzero are to
  be returned.

  \param \in requested_rows: requested_rows[0...N_requested_rows-1] give
  the row indices of the rows for which nonzero values are returned.

  \param \out row_lengths: row_lengths[i] is the number of nonzeros in
  the row 'requested_rows[i]'

  \param \out columns,values: columns[k] and values[k] contains the
  column number and value of a matrix nonzero where all nonzeros for
  requested_rows[i] appear before requested_rows[i+1]'s nonzeros.  NOTE:
  Arrays are of size 'allocated_space'.

  \param \in allocated_space: indicates the space available in 'columns'
  and 'values' for storing nonzeros. If more space is needed, return 0.
*/
int Epetra_ML_getrow(ML_Operator *data, int N_requested_rows,
                 int requested_rows[], int allocated_space, int columns[],
                 double values[], int row_lengths[]);

//! Updates ghost nodes.
/*! Update vec's ghost node via communication. 
  
    \param \inout vec: On input, vec contains data. On output, ghost values
    are updated.
 
    \param \inout data: On input, points to user's data containing
    matrix values and communication information.

*/
int Epetra_ML_comm_wrapper(double vec[], void *data);

//! Converts a (square) Epetra_RowMatrix into a ML operator.
int EpetraMatrix2MLMatrix(ML *ml_handle, int level,
                                Epetra_RowMatrix * Amat);

//! Wraps an Epetra_RowMatrix into a ML_Operator
int Epetra2MLMatrix(Epetra_RowMatrix * A, ML_Operator *Result);

//! Multiplies two Epetra_RowMatrix's.
Epetra_CrsMatrix *Epetra_MatrixMult(Epetra_RowMatrix *B, Epetra_RowMatrix *Bt);

//! Adds two Epetra_RowMatrix's.
Epetra_CrsMatrix *Epetra_MatrixAdd(Epetra_RowMatrix *B, Epetra_RowMatrix *Bt, double scalar);
int ML_Epetra_CRSinsert(ML_Operator *, int, int *, double *, int);

int ML_Operator2EpetraCrsMatrix(ML_Operator *Ke, Epetra_CrsMatrix * &
				CrsMatrix, int & MaxNumNonzeros,
				bool CheckNonzeroRow, double &);

/* This Proto-type is in ml_rap.h. This prevents ml_rap.c and ml_matmat_mult.c     */
/* from including ml_epetra_utils.h which would require the C++ compiler for these */
/* files.             
extern int  ML_back_to_epetraCrs(ML_Operator *Mat1Mat2,  ML_Operator *Result, 
				 ML_Operator *Mat1,  ML_Operator *Mat2); 
*/


#endif /* _ML_EPETRA_UTILS_H_ */
