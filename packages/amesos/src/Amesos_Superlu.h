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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
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

#ifndef AMESOS_SUPERLU_H
#define AMESOS_SUPERLU_H

#if defined(Amesos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Amesos package is deprecated"
#endif
#endif

#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
#include "Amesos_NoCopiable.h"
#include "Amesos_Utils.h"
#include "Amesos_Time.h"
#include "Amesos_Status.h"
#include "Amesos_Control.h"
#include "Teuchos_RCP.hpp"

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

class Amesos_Superlu: public Amesos_BaseSolver,
                      private Amesos_Time,
                      private Amesos_NoCopiable,
                      private Amesos_Utils,
                      private Amesos_Control,
                      private Amesos_Status {

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

  int SymbolicFactorization();

  int NumericFactorization();

  int Solve();

  //@}
  //@{ \name Additional methods

  const Epetra_LinearProblem *GetProblem() const { return(Problem_); };

  bool MatrixShapeOK() const ;

  int SetUseTranspose (bool useTheTranspose) {
    UseTranspose_ = useTheTranspose; return(0);
  }

  bool UseTranspose() const {return(UseTranspose_);};

  const Epetra_Comm& Comm() const {return(GetProblem()->GetOperator()->Comm());};

  int SetParameters( Teuchos::ParameterList &ParameterList )  ;

  //! Returns the number of symbolic factorizations performed by this object.
  int NumSymbolicFact() const { return( Amesos_Status::NumSymbolicFact_ ); }

  //! Returns the number of numeric factorizations performed by this object.
  int NumNumericFact() const { return( Amesos_Status::NumNumericFact_ ); }

  //! Returns the number of solves performed by this object.
  int NumSolve() const { return( Amesos_Status::NumSolve_ ); }

  //! Prints timing information.
  void PrintTiming() const;

  //! Prints status information.
  void PrintStatus() const;

  //! Extracts timing information from the current solver and places it in the parameter list.
  void GetTiming( Teuchos::ParameterList &TimingParameterList ) const { Amesos_Time::GetTiming(TimingParameterList); }

private:

  //@}
  //@{ Utility methods

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

  //! Factors the matrix, no previous factorization available.
  int Factor();
  //! Re-factors the matrix.
  int ReFactor();

  //! Sets up the matrix on processor 0.
  int ConvertToSerial();

  //!  PerformNumericFactorization - Call Superlu to perform numeric factorization
  // Note:  All action is performed on process 0
  int PerformNumericFactorization();

  //@}

  //! Main structure for SuperLU.
  SLUData* data_;
  std::vector<double> berr_;
  std::vector<double> ferr_;
  std::vector<int> perm_r_;
  std::vector<int> perm_c_;
  std::vector<int> etree_;
  std::vector<double> R_;
  std::vector<double> C_;
  char equed_;
  // no idea of the following.
  double* DummyArray;

  //!< stores the matrix in SuperLU format.
  std::vector <int> Ap_;
  //!< stores the matrix in SuperLU format.
  std::vector <int> Ai_;
  //!< stores the matrix in SuperLU format.
  std::vector <double> Aval_;
  //! Global size of the matrix.
  long long NumGlobalRows_;
  //! Global number of nonzeros in the matrix.
  long long NumGlobalNonzeros_;
  //! If \c true, solve the linear system with the transpose of the matrix.
  bool UseTranspose_;
  //! If \c true, the factorization has been successfully computed.
  bool FactorizationOK_;
  bool FactorizationDone_;
  bool ReuseSymbolic_;
  //! Process number (i.e. Comm().MyPID()
  int iam_;
  //! Quick access pointer to internal timing data.
  int MtxConvTime_, MtxRedistTime_, VecRedistTime_;
  int NumFactTime_, SolveTime_, OverheadTime_;
  //! Contains a map with all elements assigned to processor 0.
  Teuchos::RCP<Epetra_Map> SerialMap_;
  //! Contains a matrix with all rows assigned to processor 0.
  Teuchos::RCP<Epetra_CrsMatrix> SerialCrsMatrixA_;
  //! Importer from distributed to SerialMap_.
  Teuchos::RCP<Epetra_Import> ImportToSerial_;
  //! For parallel runs, stores the matrix defined on SerialMap_.
  Epetra_RowMatrix* SerialMatrix_ ;
  //! Pointer to the user's defined linear problem.
  const Epetra_LinearProblem* Problem_;
  //! Pointer to the linear system matrix.
  Epetra_RowMatrix* RowMatrixA_;

};  // End of  class Amesos_Superlu
#endif /* AMESOS_SUPERLU_H */
