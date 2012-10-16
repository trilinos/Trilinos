/*!
 *  \file ml_epetra_utils.h
 *
 *  \brief Interface to the Trilinos package Anasazi.
 *
 *  \date Last update to Doxygen: 22-Jul-04
 *
 */
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */


#ifndef _ML_EPETRA_UTILS_H_
#define _ML_EPETRA_UTILS_H_

class Epetra_Comm;
class Epetra_BlockMap;
class Epetra_MultiVector;
class Epetra_RowMatrix;
class Epetra_Map;
class Epetra_Vector;
class Epetra_IntVector;
class Epetra_Import;
class Epetra_Object;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
class Epetra_RowMatrix;
class Epetra_LinearProblem;
class Epetra_SerialDenseMatrix;
namespace Teuchos {
  class ParameterList;
}

#include "ml_common.h"

#ifdef ML_MPI
#ifndef EPETRA_MPI
#define EPETRA_MPI
#endif
#include "mpi.h"
#endif
#include "ml_include.h"
#include <iostream>

#ifdef HAVE_ML_EPETRAEXT
#include "EpetraExt_SolverMap_CrsMatrix.h"
#endif


#ifndef ML_CPP
extern "C" {
#endif

#include "Epetra_DataAccess.h"
  
// ====================================================================== 
//! Matrix-vector function for Epetra matrices.
/*! This is the ML matrix-vector wrap for Epetra matrices.
 */
int ML_Epetra_matvec(ML_Operator *data, int in, double *p,
                                        int out, double *ap);
int ML_Epetra_RowMatrix_matvec(ML_Operator *data, int in, double *p,
                                                  int out, double *ap);
int ML_Epetra_CrsMatrix_matvec(ML_Operator *data, int in, double *p,
                                                  int out, double *ap);
int ML_Epetra_VbrMatrix_matvec(ML_Operator *data, int in, double *p,
                                                  int out, double *ap);
int Epetra_ML_GetCrsDataptrs(ML_Operator *data, double **values, int **cols, int **rowptr);

#ifdef WKC
int ML_Epetra_matvec_WKC(ML_Operator *data, int in, double *p, int out,
                 double *ap);
#endif 

//! Getrow function for matrix of type Epetra_RowMatrix.
/*!
  Supply local matrix (without ghost node columns) for rows given by
  requested_rows[0 ... N_requested_rows-1].  Return this information in
  'row_lengths, columns, values'.  If there is not enough space to complete
  this operation, return 0. Otherwise, return 1.
 
  \param data (In)             
         Points to user's data containing matrix values. 
  \param N_requested_rows (In) Number of rows for which nonzero are to be
                   returned.
  \param requested_rows (In)  
         Requested_rows[0...N_requested_rows-1] give the
                   row indices of the rows for which nonzero values are
                   returned.
  \param row_lengths (Out)
         Row_lengths[i] is the number of nonzeros in the
         row 'requested_rows[i]'
  \param columns,values (Out)  
         Columns[k] and values[k] contains the column
         number and value of a matrix nonzero where all nonzeros for
         requested_rows[i] appear before requested_rows[i+1]'s
         nonzeros.  NOTE: Arrays are of size 'allocated_space'.
  \param allocated_space  (In)
         Indicates the space available in 'columns' and
         'values' for storing nonzeros. If more space is needed,
         return 0.
 */
int ML_Epetra_getrow(ML_Operator *data, int N_requested_rows,
                 int requested_rows[], int allocated_space, int columns[],
                 double values[], int row_lengths[]);

int ML_Epetra_RowMatrix_getrow(ML_Operator *data, int N_requested_rows, 
                               int requested_rows[], int allocated_space, 
                               int columns[], double values[],
                               int row_lengths[]);
  
int ML_Epetra_CrsMatrix_getrow(ML_Operator *data, int N_requested_rows,
            int requested_rows[], 
		    int allocated_space, int columns[], double values[],
                               int row_lengths[]);
int ML_Epetra_CrsMatrix_get_one_row(ML_Operator *data, int N_requested_rows,
                                    int requested_rows[], 
                                    int allocated_space, int columns[], double values[],
                                    int row_lengths[]);
  

int ML_Epetra_VbrMatrix_getrow(ML_Operator *data,
            int N_requested_rows, int requested_rows[], 
		    int allocated_space, int columns[], double values[],
		    int row_lengths[]);

void ML_Set_Filter(Teuchos::ParameterList& List);

int ML_Epetra_matvec_Filter(ML_Operator *data, int in, double *p, 
                            int out, double *ap);

int ML_Epetra_getrow_Filter(ML_Operator *data, int N_requested_rows,
                            int requested_rows[], int allocated_space, int columns[],
                            double values[], int row_lengths[]);
//! Update vec's ghost node via communication.
/*! Update vec's ghost node via communication. Note: the length of vec is
  given by N_local + N_ghost where Amat was created via
                  \c AZ_matrix_create(N_local);
  and a 'getrow' function was supplied via
                  \c AZ_set_MATFREE_getrow(Amat,,,,N_ghost,).
 
  \param vec Vector containing data. On output, ghost values
                    are updated.
 
  \param data  points to user's data containing matrix values.
                   and communication information.
 */ 
int ML_Epetra_comm_wrapper(double vec[], void *data);
int ML_Epetra_CrsMatrix_comm_wrapper(double vec[], void *data);
int ML_Epetra_VbrMatrix_comm_wrapper(double vec[], void *data);

// wrappers for Epetra_CrsGraph
int ML_Epetra_CrsGraph_comm_wrapper(double vec[], void *data);
int ML_Epetra_CrsGraph_matvec(ML_Operator *data, int in, double *p,
                              int out, double *ap);
int ML_Epetra_CrsGraph_getrow(ML_Operator *data, int N_requested_rows,
                              int requested_rows[], int allocated_space, 
                              int columns[], double values[],
                              int row_lengths[]);
int ML_Operator_WrapEpetraCrsGraph(Epetra_CrsGraph* Graph, ML_Operator *newMatrix);

#ifdef HAVE_ML_TEUCHOS
void ML_CreateSublists(const Teuchos::ParameterList &List,
                      Teuchos::ParameterList &newList,
                      int NumLevels);


#endif

#ifndef ML_CPP
}
#endif

//! Wraps an Epetra_RowMatrix into an ML_Operators.
/*! This function creates an ML_Operator that is based on the input 
 *  Epetra_RowMatrix. This is a "cheap" wrap in the sense that
 *  only function and pointers are created. Data is still coded as an
 *  Epetra_RowMatrix.
 *
 *  \note ML requires A->RowMatrixRowMap() == A->OperatorRangeMap()
 */
int EpetraMatrix2MLMatrix(ML *ml_handle, int level,
                                Epetra_RowMatrix * Amat);

//! Wraps an Epetra_RowMatrix into an ML_Operators, for the given level.
/*! This function creates an ML_Operator that is based on the input 
 *  Epetra_RowMatrix. This is a "cheap" wrap in the sense that
 *  only function and pointers are created. Data is still coded as an
 *  Epetra_RowMatrix. The ML_Operator is set in the specified level of the
 *  hierarchy.
 *
 *  \note ML requires A->RowMatrixRowMap() == A->OperatorRangeMap()
 */
int ML_Operator_WrapEpetraMatrix(Epetra_RowMatrix * A, ML_Operator *Result);


//! Wraps an Epetra_CrsMatrix into an ML_Operator, for the given level.
/*! This is an *ultra* cheap wrap in that I wrap the pointers that come out of
 *  Epetra's ExtractCrsDataPointers function.
 *
 *  You need to have remapped the Epetra Matrix to include all the columns
 *  before this routine gets called or else this won't work in parallel.
 *
 *  OptimizeStorage() *must* be called on the routines for the light wrap to be
 *  used.  Otherwise, the routine calls  ML_Operator_WrapEpetraMatrix as a fallback.
 *
 *  \note ML requires A->RowMatrixRowMap() == A->OperatorRangeMap()
 */
int ML_Operator_WrapEpetraCrsMatrix(Epetra_CrsMatrix * A, ML_Operator *newMatrix, bool verbose=false);

//! Wraps a ML_Operator into a Epetra_CrsMatrix.
/*! This is a somewhat cheap wrap in that the pointers get dumped into the Epetra_CrsMatrix.  
 */
void Epetra_CrsMatrix_Wrap_ML_Operator(ML_Operator * A, const Epetra_Comm &Comm, const Epetra_Map &RowMap,Epetra_CrsMatrix **Result,Epetra_DataAccess CV=View,int base=0);
//void Epetra_CrsMatrix_Wrap_ML_Operator(ML_Operator * A, Epetra_CrsMatrix *Result);


//! Multiplies two Epetra_RowMatrix's, returns the results as an Epetra_CrsMatrix.
Epetra_CrsMatrix *Epetra_MatrixMult(Epetra_RowMatrix *B, Epetra_RowMatrix *Bt);

//! Adds two Epetra_RowMatrix's, returns the result as an Epetra_CrsMatrix
Epetra_CrsMatrix *Epetra_MatrixAdd(Epetra_RowMatrix *B, Epetra_RowMatrix *Bt, double scalar);
int ML_Epetra_CRSinsert(ML_Operator *, int, int *, double *, int);

//! Converts an ML_Operator into an Epetra_CrsMatrix
/*! This function creates a new Epetra_CrsMatrix, and inserts all the nonzero
 * elements of the ML_Operator in it. This is an expensive conversion, in the
 * sense that the Epetra_RowMatrix is a \sl copy of the input ML_Operator.
 *
 * \note This function can be used with rectangular matrices.
 */
int ML_Operator2EpetraCrsMatrix(ML_Operator *Ke, Epetra_CrsMatrix * &
				CrsMatrix, int & MaxNumNonzeros,
				bool CheckNonzeroRow, double &, int base=0,bool verbose=false);

inline int ML_Operator2EpetraCrsMatrix(ML_Operator *Ke, Epetra_CrsMatrix * &
				CrsMatrix)
{
  int MaxNumNonzeros;
  double CPUTime;
  return(ML_Operator2EpetraCrsMatrix(Ke, CrsMatrix, MaxNumNonzeros, false, CPUTime));
}

Epetra_Map* Epetra_ML_readupdatevector(char* filename, Epetra_Comm& comm);
Epetra_CrsMatrix* Epetra_ML_readaztecmatrix(char* filename,Epetra_Map& map,
                                            Epetra_Comm& comm);



namespace ML_Epetra{

  //! Does an P^TAP for Epetra_CrsMatrices using ML's kernels.
  int ML_Epetra_PtAP(const Epetra_CrsMatrix & A, const Epetra_CrsMatrix & P, Epetra_CrsMatrix *&Result,bool verbose=false);

  //! Does an RAP for Epetra_CrsMatrices using ML's kernels.
  int ML_Epetra_RAP(const Epetra_CrsMatrix & A, const Epetra_CrsMatrix & P, const Epetra_CrsMatrix & R, Epetra_CrsMatrix *&Result,bool verbose);

  //! Finds the Dirichlet rows in a square matrix that got the one-and-zeros
  //treatment
  /*! Returns the local Dirichlet rows for a square matrix that go the
   *  ones-and-zeros treatment for BCs.
   */
  int* FindLocalDiricheltRowsFromOnesAndZeros(const Epetra_CrsMatrix & Matrix, int &numBCRows);

  //! Applies Dirichlet conditions to columns that rows already have.
  /*! Apply the Dirichlet BC's specified by the dirichletRows to remove all
   * columns (across all processors) that have entries zero'd by the BC's.  This
   * will symmetrize a square matrix, though the routine can be run on
   *  non-square matrices.
   */  
  void Apply_BCsToMatrixColumns(const int *dirichletRows, int numBCRows, const Epetra_CrsMatrix & Matrix);
  
  //! Applies Dirichlet conditions to matrix rows.
  void Apply_BCsToMatrixRows(const int *dirichletRows, int numBCRows, const Epetra_CrsMatrix & Matrix);

  //! Applies Dirichlet conditions to columns that rows already have.
  /*! Apply the Dirichlet BC's specified by BoundaryMatrix to remove all
   * columns (across all processors) that have entries zero'd by the BC's.  This
   * will symmetrize a square matrix, though the routine can be run on
   *  non-square matrices.
   */  
  void Apply_BCsToMatrixColumns(const Epetra_RowMatrix & iBoundaryMatrix, const Epetra_RowMatrix & iMatrix);
  void Apply_BCsToMatrixColumns(const Epetra_IntVector &dirichletColumns,const Epetra_CrsMatrix & Matrix);

  
  //! Applies boundary conditions to gradient matrix.  (Maxwell's equations)
  void Apply_BCsToGradient( const Epetra_RowMatrix & EdgeMatrix,
                            const Epetra_RowMatrix & T);


  //! Does Row/Column OAZ to a matrix.
  void Apply_OAZToMatrix(int *dirichletRows, int numBCRows, const Epetra_CrsMatrix & Matrix);

  //! Returns the local column numbers of the local rows passed in.
  Epetra_IntVector * LocalRowstoColumns(int *Rows, int numRows,const Epetra_CrsMatrix & Matrix);


    //! Finds Dirichlet the local Dirichlet columns, given the local Dirichlet rows
  Epetra_IntVector * FindLocalDirichletColumnsFromRows(const int *dirichletRows, int numBCRows,const Epetra_CrsMatrix & Matrix);

  //! Drops a 1 on the diagonal of zero'd our rows
  void Remove_Zeroed_Rows(const Epetra_CrsMatrix & Matrix, double tol=0.0);


  //! Transforms Epetra matrix column map (if necessary) to be compatible with
  /*! how ML handles column indices.  Any matrix that cannot be dynamically
      cast to an Epetra_CrsMatrix will not be changed.

      \param A - (In) Matrix that is to be transformed.
      \param transform - (In) EpetraExt widget that does the transformation.
      \param matrixName - (In) Optional label for the incoming matrix.
   */

#ifdef HAVE_ML_EPETRAEXT
  Epetra_RowMatrix* ModifyEpetraMatrixColMap( const Epetra_RowMatrix &A,
                                   EpetraExt::CrsMatrix_SolverMap &transform,
                                   const char* matrixName=0, bool verbose=false);
#endif

#ifdef HAVE_ML_TEUCHOS
int UpdateList(Teuchos::ParameterList &source, Teuchos::ParameterList &dest, bool OverWrite);
#endif
 
}


#ifdef FIXME
string ML_toString(const int& x);
string ML_toString(const double& x);
#endif


#ifdef __cplusplus
extern "C" {
#endif
int ML_Operator_Destroy_DiscreteLaplacian();
int ML_Operator_DiscreteLaplacian(ML_Operator* Op, int SymmetricPattern,
				  double* x_coord, double* y_coord,
				  double* z_coord, double theta,
				  ML_Operator** NewOp);
bool Epetra_ML_readaztecvector(char* filename, Epetra_MultiVector& Vector, 
                               Epetra_Map& map,Epetra_Comm& comm, int ivec);
bool Epetra_ML_readvariableblocks(char* filename, Epetra_Map& map,
                                  Epetra_Comm& comm, 
                                  int**blocks, int** block_pde);
bool Epetra_ML_writegidviz(char* filename, int label, 
                           Epetra_MultiVector& vector, int ivec, 
                           Epetra_Map& map, Epetra_Comm& comm);
                                  
#ifdef __cplusplus
}
#endif
  
//! Stops the code, waiting for a debugger to attach
/*! BreakForDebugger() is useful when the user wants to attach to the running
* process(es). This is a very easy task for serial runs -- just run gdb.
* Parallel runs may result more problematic. In this case, one can proceed as
* follows:
* - Either define the enviromental variable ML_BREAK_FOR_DEBUGGER
*   or create an empty file called ML_Debug_Now
* - Run the parallel code on a terminal (example, \c mpirun \c -np \c 4 \c
*   ml_example.exe )
* - The code will stop and print each process ID.
* - In another terminal, attach to the desired process.
* - Type one character to let the code continue, and debug as required.
*/
void ML_BreakForDebugger(const Epetra_Comm& comm);

#endif /* _ML_EPETRA_UTILS_H_ */
