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

#ifndef AMESOS_LAPACK_H
#define AMESOS_LAPACK_H

#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
#include "Epetra_Comm.h"
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Import.h"
class Epetra_RowMatrix;
class Epetra_LinearProblem;

//! Amesos_Lapack: an interface to LAPACK.
/*!
Class Amesos_Lapack enables the solution of the distributed linear
system, defined by an Epetra_LinearProblem, using LAPACK.

Amesos_Lapack stores the lineaar system matrix as an
Epetra_SerialDensMatrix. The linear problem is an Epetra_SerialDenseProblem.
Amesos_Lapack factorizes the matrix using DGETRF().

\date Sep-04

\author Marzio Sala, 9214.

*/
class Amesos_Lapack: public Amesos_BaseSolver { 

public: 

  //@{ \name Constructor methods
  //! Amesos_Lapack Constructor.
  /*! Creates an Amesos_Lapack instance, using an Epetra_LinearProblem,
      passing in an already-defined Epetra_LinearProblem object. 

      Note: The operator in LinearProblem must be an
      Epetra_RowMatrix.

  */
  Amesos_Lapack(const Epetra_LinearProblem& LinearProblem );

  //! Amesos_Lapack Destructor.
  /*! Completely deletes an Amesos_Lapack object.  
  */
  ~Amesos_Lapack(void);
  //@}

  //@{ \name Mathematical functions.

    //! Performs SymbolicFactorization on the matrix A.
    /*! 
      In addition to performing symbolic factorization on the matrix A, 
      the call to SymbolicFactorization() implies that no change will
      be made to the non-zero structure of the underlying matrix without 
      a subsequent call to SymbolicFactorization().
      
    \return Integer error code, set to 0 if successful.
  */
    int SymbolicFactorization() ;

    //! Performs NumericFactorization on the matrix A.
    /*!  In addition to performing numeric factorization (and symbolic
      factorization if necessary) on the matrix A, the call to
      NumericFactorization() implies that no change will be made to
      the underlying matrix without a subsequent call to
      NumericFactorization().  

     \return Integer error code, set to 0 if successful.
  */
    int NumericFactorization() ;

    //! Solves A X = B (or A<SUP>T</SUP> X = B) 
    /*! 

     \return Integer error code, set to 0 if successful.
     */
    int Solve();

  //@}
  
  //@{ \name Additional methods required to support the Epetra_Operator interface.

#ifdef FIXME
  //! Returns a character string describing the operator
  char * Label() const {return(Epetra_Object::Label());};
#endif
    
  //! Get a pointer to the Problem.
  const Epetra_LinearProblem *GetProblem() const { return(Problem_); };

  //! Returns true if Amesos_Lapack can handle this matrix shape 
  /*! Returns true if the matrix shape is one that KLU can
    handle. Amesos_Lapack only works with square matrices.  
  */
  bool MatrixShapeOK() const ;

  //! SetUseTranpose(true).
  int SetUseTranspose(bool UseTranspose) {
    UseTranspose_ = UseTranspose; 
    return(0);
  }

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(UseTranspose_);};

  //! Returns a pointer to the Epetra_Comm communicator associated with this matrix.
  const Epetra_Comm & Comm() const {
    return(GetProblem()->GetOperator()->Comm());
  }

  //!  Updates internal variables. 
  /*!  
    \return Integer error code, set to 0 if successful. 
   */
  int SetParameters( Teuchos::ParameterList &ParameterList )  ;

  //! Print timing information
  void PrintTiming();
  
  //! Print information about the factorization and solution phases.
  void PrintStatus();
  
  //@}

protected:

  //! Returns true if SymbolicFactorization() has been successfully called.
  bool IsSymbolicFactorizationOK()
  {
    return(IsSymbolicFactorizationOK_);
  }

  //! Returns true if SymbolicFactorization() has been successfully called.
  bool IsNumericFactorizationOK()
  {
    return(IsNumericFactorizationOK_);
  }
  
  //! Returns a pointer to the linear system matrix.
  const Epetra_RowMatrix* Matrix() const
  {
    return(dynamic_cast<Epetra_RowMatrix*>(Problem_->GetOperator()));
  }

  //! Returns the number of global rows, or -1 if Matrix() returns 0.
  int NumGlobalRows() const
  {
    if (Matrix())
      return(Matrix()->NumGlobalRows());
    else
      return(-1);
  }

  //! Returns the number of local rows, or -1 if Matrix() returns 0.
  int NumMyRows() const
  {
    if (Matrix())
      return(Matrix()->NumMyRows());
    else
      return(-1);
  }

  //! Returns the ID of calling process.
  int MyPID() const
  {
    return(Matrix()->Comm().MyPID());
  }

  //! Returns the number of processes in communicator.
  int NumProc() const
  {
    return(Matrix()->Comm().NumProc());
  }

  //! Returns a reference to serial map (that with all elements on process 0). Builds SerialMap_ if necessary or required.
  const Epetra_Map& SerialMap()
  {
    if (SerialMap_ == 0) {
      int NumElements = 0;
      if (MyPID() == 0)
	NumElements = NumGlobalRows();

      SerialMap_ = new Epetra_Map(-1,NumElements,0,Matrix()->Comm());
      assert (SerialMap_ != 0);
    }

    return(*SerialMap_);
  }

  //! Returns a reference to serial matrix (that with all rows on process 0). Builds SerialMap_ if necessary or required.
  Epetra_CrsMatrix& SerialMatrix()
  {
    if (SerialMatrix_ == 0) {
      SerialMatrix_ = new Epetra_CrsMatrix(Copy,SerialMap(),0);
      assert (SerialMatrix_ != 0);
    }

    return(*SerialMatrix_);
  }

  //! Returns a reference to the importer map. Builds SerialMap_ if necessary or required.
  const Epetra_Import& RowImporter()
  {
    if (RowImporter_ == 0) {
      RowImporter_ = new Epetra_Import(SerialMap(),Matrix()->RowMatrixRowMap());
      assert (RowImporter_ != 0);
    }

    return(*RowImporter_);
  }

  //! Solves the linear system, when only one process is used.
  int SolveSerial(Epetra_MultiVector& X,
		  const Epetra_MultiVector& B);
  
  //! Solves the linear system, when more than one process is used.
  int SolveDistributed(Epetra_MultiVector& X,
		       const Epetra_MultiVector& B);
  
  //! Convert a serial matrix to dense format. Only for Comm().NumProc() == 1.
  int SerialToDense();

  //! Convert a distributed matrix to dense format. Only for Comm().NumProc() > 1
  int DistributedToDense();

  //! Points to the Serial matrix (defined on process 0 only).
  Epetra_CrsMatrix* SerialMatrix_;
  //! Points to a Serial Map (to import LHS/RHS).
  Epetra_Map* SerialMap_;
  //! Importer from distributed map to SerialMap_.
  Epetra_Import* RowImporter_;
  //! Dense matrix.
  Epetra_SerialDenseMatrix DenseMatrix_;
  //! Dense LHS.
  Epetra_SerialDenseMatrix DenseLHS_;
  //! Dense RHS.
  Epetra_SerialDenseMatrix DenseRHS_;
  //! Linear problem for dense matrix and vectors.
  Epetra_SerialDenseSolver DenseSolver_;

  //! If \c true, the linear system with the transpose will be solved.
  bool UseTranspose_;
  //! Pointer to the linear problem.
  const Epetra_LinearProblem* Problem_;

  //! If \c true, some timing are printed in destructor.
  bool PrintTiming_;
  //! If \c true, some additional information are printed in destructor.
  bool PrintStatus_;
  //! If \c true, solution and rhs norms are computed in Solve().
  bool ComputeVectorNorms_;
  //! If \c true, the norm of the residual in computed in Solve().
  bool ComputeTrueResidual_;
  //! If greater than 0, prints out additional information.
  int verbose_;
  //! Ignore all elements of the matrix whose absolute value if less than Threshold_.
  double Threshold_;
  //! Value to be added to the diagonal.
  double AddToDiag_;
  //! If \c true, SymbolicFactorization() has been successfully called.
  bool IsSymbolicFactorizationOK_;
  //! If \c true, NumericFactorization() has been successfully called.
  bool IsNumericFactorizationOK_;

  //! Time to convert to LAPACK format.
  double ConTime_;
  //! Time for symbolic factorization.
  double SymTime_;
  //! Time for numeric factorization.
  double NumTime_;
  //! Time for solution.
  double SolTime_;
  //! Time to redistribute vectors.
  double VecTime_;
  //! Time to redistribute matrix.
  double MatTime_;

  //! Number of calls to SymbolicFactorization().
  int NumSymbolicFact_;
  //! Number of calls to NumericFactorization().
  int NumNumericFact_;
  //! Number of calls to Solver().
  int NumSolve_;

  //! Object used to track times.
  Epetra_Time * Time_;
  
};  // End of  class Amesos_Lapack
#endif /* AMESOS_LAPACK_H */
