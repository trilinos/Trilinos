
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

#ifndef AMESOS_UMFPACK_H
#define AMESOS_UMFPACK_H

#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "Epetra_Import.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_Comm.h"
#endif

//! Class Amesos_Umfpack:  An object-oriented wrapper for UMFPACK.
/*!  Amesos_Umfpack will solve a linear systems of equations: <TT>A X = B</TT>
   using Epetra objects and the UMFPACK solver library, where
  <TT>A</TT> is an Epetra_RowMatrix and <TT>X</TT> and <TT>B</TT> are 
  Epetra_MultiVector objects.

 
*/
class Amesos_Umfpack: public Amesos_BaseSolver { 

public: 

  //@{ \name Constructor methods
  //! Amesos_Umfpack Constructor.
  /*! Creates an Amesos_Umfpack instance, using an Epetra_LinearProblem,
      passing in an already-defined Epetra_LinearProblem object. 

      Note: The operator in LinearProblem must be an
      Epetra_RowMatrix.

  */
  Amesos_Umfpack( const Epetra_LinearProblem& LinearProblem );

  //! Amesos_Umfpack Destructor.
  /*! Completely deletes an Amesos_Umfpack object.  
  */
  ~Amesos_Umfpack(void);
  //@}

  //@{ \name Mathematical functions.

    //! Performs SymbolicFactorization on the matrix A.
    /*! 
      In addition to performing symbolic factorization on the matrix A, 
      the call to SymbolicFactorization() implies that no change will
      be made to the non-zero structure of the underlying matrix without 
      a subsequent call to SymbolicFactorization().
      
      preconditions:<ul>
      <li>GetProblem().GetOperator() != 0 (return -1)
      <li>MatrixShapeOk(GetProblem().GetOperator()) == true (return -6)
      </ul>

      postconditions:<ul>
      <li>Symbolic Factorization will be performed (or marked to be performed) 
      allowing NumericFactorization() and Solve() to be called.
      </ul>

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

    //! Solves A X = B (or A<SUP>T</SUP> x = B) 
    /*!
     \return Integer error code, set to 0 if successful.
  */
    int Solve();

  //@}
  
  //@{ \name Additional methods required to support the Epetra_Operator interface.

  //! Get a pointer to the Problem.
  const Epetra_LinearProblem *GetProblem() const { return(Problem_); };

  //! Returns true if UMFPACK can handle this matrix shape 
  /*! Returns true if the matrix shape is one that UMFPACK can
    handle. UMFPACK only works with square matrices.  
  */
  bool MatrixShapeOK() const ;

  int SetUseTranspose(bool UseTranspose) {UseTranspose_ = UseTranspose; return(0);};

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(UseTranspose_);};

  //! Returns a pointer to the Epetra_Comm communicator associated with this matrix.
  const Epetra_Comm & Comm() const {return(GetProblem()->GetOperator()->Comm());};

  //! Returns an estimate of the reciprocal of the condition number 
  /*  Rcond is an estimate of the reciprocal of the condition number of the 
      matrix at the time of the most recent call to NumericFactorization()
      Rcond = min(abs(diag))/max(abs(diag)) see Umfpack documentatoin
      for details.  
   */
  double GetRcond() const {return(Rcond_);}; 

  //! Sets parameters from the parameters list, returns 0 if successful.
  int SetParameters( Teuchos::ParameterList &ParameterList ) ;

  //! Prints timing information
  void PrintTiming();
  
  //! Prints information about the factorization and solution phases.
  void PrintStatus();

  //@}

private:  
  Epetra_RowMatrix* Matrix()
  {
    return(dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator()));
  }

  //! Converts matrix to a serial Epetra_CrsMatrix
  int ConvertToSerial();

  /*
    ConvertToUmfpackCRS - Convert matirx to form expected by Umfpack: Ai, Ap, Aval
    Preconditions:
      numentries_, NumGloalElements_ and SerialMatrix_ must be set.
    Postconditions:
      Ai, Ap, and Aval are resized and populated with a compresses row storage 
      version of the input matrix A.
  */
  int ConvertToUmfpackCRS();     

  /*
    PerformSymbolicFactorization - Call Umfpack to perform symbolic factorization
    Preconditions:
      IsLocal must be set to 1 if the input matrix is entirely stored on process 0
      Ap, Ai and Aval are a compressed row storage version of the input matrix A.
    Postconditions:
      Symbolic points to an UMFPACK internal opaque object containing the
        symbolic factorization and accompanying information.  
      SymbolicFactorizationOK_ = true; 
    Note:  All action is performed on process 0
  */
      
  int PerformSymbolicFactorization(); 

  /*
    PerformNumericFactorization - Call Umfpack to perform numeric factorization
    Preconditions:
      IsLocal must be set 
      Ap, Ai and Aval are a compressed row storage version of the input matrix A.
      Symbolic must be set
    Postconditions:
      Numeric points to an UMFPACK internal opaque object containing the
        numeric factorization and accompanying information.  
      NumericFactorizationOK_ = true; 
    Note:  All action is performed on process 0
  */
  int PerformNumericFactorization(); 

 protected:

  //! True if SymbolicFactorization has been done
  bool IsSymbolicFactorizationOK_;
  //! True if NumericFactorization has been done
  bool IsNumericFactorizationOK_;

  //! Umfpack internal opaque object
  void *Symbolic;
  //! Umfpack internal opaque object
  void *Numeric;

  //!  Ap, Ai, Aval form the compressed row storage used by Umfpack
  vector <int> Ap;
  vector <int> Ai;
  vector <double> Aval;

  //! Process number (i.e. Comm().MyPID() )
  int iam;
  //! 1 if Problem_->GetOperator() is stored entirely on process 0
  int IsLocal_;
  //! Number of non-zero entries in Problem_->GetOperator()
  int numentries_;
  //! Number of rows and columns in the Problem_->GetOperator()
  int NumGlobalElements_;

  //! Points to a Serial Map (unused if IsLocal == 1 ) 
  Epetra_Map *SerialMap_;
  //! Points to a Serial Copy of A
  /* If IsLocal==1 - Points to the original matrix 
   * If  IsLocal==0 - Points to SerialCrsMatrixA
   */
  Epetra_RowMatrix *SerialMatrix_;

  //! If \c true, solve the problem with the transpose.
  bool UseTranspose_;
  //! Pointer to the linear problem to solve.
  const Epetra_LinearProblem * Problem_;

  //! Reciprocal condition number estimate
  double Rcond_;

  //! If \c true, prints timing in the destructor.
  bool PrintTiming_;
  //! If \c true, prints additional information in the destructor.
  bool PrintStatus_;
  //! If \c true, prints the norm of LHS and RHS in Solve().
  bool ComputeVectorNorms_;
  //! If \c true, prints the norm of the computed residual in Solve().
  bool ComputeTrueResidual_;
  
  //! Toggles the output level.
  int verbose_;

  //! time to convert to MUMPS format
  double ConTime_;
  //! time for symbolic factorization
  double SymTime_;
  //! time for numeric factorization
  double NumTime_;
  //! time for solution
  double SolTime_;
  //! time to redistribute vectors
  double VecTime_;
  //! time to redistribute matrix
  double MatTime_;
  
  //! Number of symbolic factorizations.
  int NumSymbolicFact_;
  //! Number of numeric factorizations.
  int NumNumericFact_;
  //! Number of solves.
  int NumSolve_;  

  //! Used to track timing.
  Epetra_Time * Time_;
  //! Importer from distributed to serial (all rows on process 0).
  Epetra_Import * ImportToSerial_;
  
};  // class Amesos_Umfpack  
#endif /* AMESOS_UMFPACK_H */
