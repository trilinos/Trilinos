// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

//
//  Coding tasks:
//  1)  Create the dense matrices in Solve()   DONE 
//  2)  Make the call to dgsvx() in Factor() dONE 
//  3)  Make the call to dgsvx() in Solve()  DONE
//
//  Later coding tasks:
//  0)  Factor called twice
//  1)  Refactor()
//  2)  Parameter list
//  3)  Transpose  - DONE 
//  4)  Destructor - In particular, need to call the SuperLU_FREE routines.
//  5)  Coments - especially in Amesos_Superlu.h
//
//  SymbolicFactorization() performs no action other than making sure that Factorization 
//    s performed
//  NumericFactorization() performs the factorization but no solve (right hand side is
//    is set to 0 vectors) reuses factors only if ReuseFactorization_ is set. 
//    If FactorizationOK_() && ReuseSymbolic_ 
//       ReFactor()
//    else
//       Factor() 
//
//  Solve() 
//
//  Factor() does everything from scratch:
//    Redistributes the data if necessary
//      Deletes any data structures left over from the previous call to Factor()
//    Copies the data into the format that SuperLU wants it
//    Calls dgssvx to factor the matrix with factor set to true
//  ReFactor()
//    Redistributes the data if necessary
//      - Attempting to check to make sure that the non-zero structure is unchanged
//    Copies the data into the format that SuperLU already has it
//       FIRST PASS - assert( false ) 
//    

#ifndef _AMESOS_SUPERLU_H_
#define _AMESOS_SUPERLU_H_

#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
#include "Teuchos_RefCountPtr.hpp"

class Epetra_Time;
class SLUData;
class Epetra_Comm;
class Epetra_CrsMatrix;
class Epetra_LinearProblem;

//! Amesos_Superlu:  Amesos interface to Xioye Li's SuperLU 3.0 serial code.
/*! 
 * Class Amesos_Superlu solves the linear systems of equations <TT>A X = B</TT>,
 * where A is defined as an Epetra_RowMatrix, and X and B are two 
 * Epetra_MultiVector's.
 *
 * \date Last updated on 28-Apr-05.
*/

class Amesos_Superlu: public Amesos_BaseSolver { 

public: 

  //@{ \name Constructor methods
  //! Amesos_Superlu Constructor.
  /*! Creates an Amesos_Superlu instance, using an Epetra_LinearProblem,
      passing in an already-defined Epetra_LinearProblem object. 

      Note: The operator in LinearProblem must be an
      Epetra_RowMatrix.

  */
  Amesos_Superlu(const Epetra_LinearProblem& LinearProblem );

  //! Amesos_Superlu Destructor.
  ~Amesos_Superlu();

  //@}

  //@{ \name Mathematical functions.

  //! Performs the symbolic factorization on the matrix A (do-nothing for this solver).
  int SymbolicFactorization();

  //! Performs the numeric factorization on the matrix A.
  /*! In addition to performing numeric factorization (and symbolic
      factorization if necessary) on the matrix A, the call to
      NumericFactorization() implies that no change will be made to
      the underlying matrix without a subsequent call to
      NumericFactorization().  

     \return Integer error code, set to 0 if successful.
  */
  int NumericFactorization();

  //! Solves A X = B (or A<SUP>T</SUP> X = B) 
  /*! 
   * Solves the linear system, after calling NumericFactorization()
   * if not yet done by the user.
   * \return Integer error code, set to 0 if successful.
   */
  int Solve();

  //@}
  
  //@{ \name Additional methods

  //! Get a pointer to the Problem.
  const Epetra_LinearProblem *GetProblem() const { return(Problem_); };

  //! Returns true if the matrix is square.
  bool MatrixShapeOK() const ;

  //! Specifies to solve the problem with A or its transpose.
  int SetUseTranspose(bool UseTranspose) {
    UseTranspose_ = UseTranspose; return(0);
  }

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(UseTranspose_);};

  //! Returns a pointer to the Epetra_Comm communicator associated with this matrix.
  const Epetra_Comm& Comm() const {return(GetProblem()->GetOperator()->Comm());};

  //! Sets the parameters as specified by the input list.
  int SetParameters( Teuchos::ParameterList &ParameterList )  ;

  //! Prints timing information.
  void PrintTiming() const;

  //! Prints status information.
  void PrintStatus() const;

  //@}

private:  

  //! Factors the matrix, no previous factorization available.
  int Factor();
  //! Re-factors the matrix.
  int ReFactor();

  //! Sets up the matrix on processor 0.
  int ConvertToSerial();

  //!  PerformNumericFactorization - Call Superlu to perform numeric factorization
  // Note:  All action is performed on process 0
  int PerformNumericFactorization(); 

  //! Returns a reference to the serial map.
  // Note: this method is delicate!
  const Epetra_Map& SerialMap() const
  {
    return(*(SerialMap_.get()));
  }

  //! Returns a reference to the importer.
  // Note: this method is delicate!
  const Epetra_Import& ImportToSerial() const
  {
    return(*(ImportToSerial_.get()));
  }

  void PrintLine() const
  {
    cout << "----------------------------------------------------------------------------" << endl;
  }

  //! Main structure for SuperLU.
  SLUData* data_;
  vector<double> berr_;
  vector<double> ferr_;
  vector<int> perm_r_;
  vector<int> perm_c_;
  vector<int> etree_;
  vector<double> R_;
  vector<double> C_;
  char equed_;
  // no idea of the following.
  double* DummyArray;

  //!< stores the matrix in SuperLU format.
  vector <int> Ap_;
  //!< stores the matrix in SuperLU format.
  vector <int> Ai_;
  //!< stores the matrix in SuperLU format.
  vector <double> Aval_; 
  //! Global size of the matrix.
  int NumGlobalRows_; 
  //! Global number of nonzeros in the matrix.
  int NumGlobalNonzeros_; 
  //! If \c true, solve the linear system with the transpose of the matrix.
  bool UseTranspose_;      
  //! If \c true, the factorization has been successfully computed.
  bool FactorizationDone_; 
  bool FactorizationOK_; 
  bool ReuseSymbolic_;
  //! Process number (i.e. Comm().MyPID() 
  int iam_;
  //! Contains a map with all elements assigned to processor 0.
  Teuchos::RefCountPtr<Epetra_Map> SerialMap_;
  //! Contains a matrix with all rows assigned to processor 0.
  Teuchos::RefCountPtr<Epetra_CrsMatrix> SerialCrsMatrixA_;
  //! Importer from distributed to SerialMap_.
  Teuchos::RefCountPtr<Epetra_Import> ImportToSerial_;
  //! For parallel runs, stores the matrix defined on SerialMap_.
  Epetra_RowMatrix* SerialMatrix_ ;
  //! Pointer to the user's defined linear problem.
  const Epetra_LinearProblem* Problem_;
  //! Pointer to the linear system matrix.
  Epetra_RowMatrix* RowMatrixA_;
  //! If \c true, the destructor prints out some status information.
  bool PrintStatus_;
  //! If \c true, the destructor prints out some timing information.
  bool PrintTiming_;
  //! Time spent in all calls to NumericFactorization().
  double NumTime_;
  //! Time spent in all calls to Solve().
  double SolTime_;
  Teuchos::RefCountPtr<Epetra_Time> Time_;
  //! Number of calls to NumericFactorization().
  int NumNumericFact_;
  //! Number of calls to Solve().
  int NumSolve_;
  //! If \c true, computes the norm of the true residual after solution.
  bool ComputeTrueResidual_;
  //! If \c true, computes the norm of the right-hand side and solution.
  bool ComputeVectorNorms_;

};  // End of  class Amesos_Superlu  
#endif /* _AMESOS_SUPERLU_H_ */
