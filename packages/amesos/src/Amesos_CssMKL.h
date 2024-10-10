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

#ifndef AMESOS_CSSMKL_H
#define AMESOS_CSSMKL_H

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
#include "Amesos_Support.h"
#include "Amesos_Control.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Comm.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

//! Amesos_CssMKL: Interface to the CSSMKL package.

/*!
  \author Marzio Sala, SNL 9214

  \date Last updated on June 2005
*/

class Amesos_CssMKL:  public Amesos_BaseSolver, 
                      private Amesos_Time, 
                      private Amesos_NoCopiable, 
                      private Amesos_Utils, 
                      private Amesos_Control, 
                      private Amesos_Status { 

public: 

  //@{ \name Constructor methods
  //! Constructor.
  Amesos_CssMKL(const Epetra_LinearProblem& LinearProblem );

  //! Destructor.
  ~Amesos_CssMKL();
  //@}

  //@{ \name Mathematical functions.

  //! Performs SymbolicFactorization on the matrix A.
  int SymbolicFactorization() ;

  //! Performs NumericFactorization on the matrix A.
  int NumericFactorization() ;

  //! Solves A X = B (or A<SUP>T</SUP> X = B) 
  int Solve();
  //@}
  
  //@{ \name Additional methods required to support the Epetra_Operator interface.

  //! Get a pointer to the Problem.
  const Epetra_LinearProblem* GetProblem() const { return(Problem_); }

  //! Returns true if CSSMKL can handle this matrix shape 
  /*! Returns true if the matrix shape is one that CSSMKL can
    handle. CSSMKL only works with square matrices.  
  */
  bool MatrixShapeOK() const;

  //! SetUseTranpose()
  /*! 
    If SetUseTranspose() is set to true, 
    \f$A^T X = B\f$ is computed.
  */  
  int SetUseTranspose(bool UseTranspose) { UseTranspose_ = UseTranspose; return(0); }

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const { return(UseTranspose_); }

  //! Returns a pointer to the Epetra_Comm communicator associated with this matrix.
  const Epetra_Comm& Comm() const { return(GetProblem()->GetOperator()->Comm()); }

  //! Set parameters from the input parameters list, returns 0 if successful.
  int SetParameters( Teuchos::ParameterList &ParameterList );

  //! Returns the number of symbolic factorizations performed by this object.
  int NumSymbolicFact() const { return( Amesos_Status::NumSymbolicFact_ ); }

  //! Returns the number of numeric factorizations performed by this object.
  int NumNumericFact() const { return( Amesos_Status::NumNumericFact_ ); }

  //! Returns the number of solves performed by this object.
  int NumSolve() const { return( Amesos_Status::NumSolve_ ); }

  //! Prints timing information
  void PrintTiming() const;
  
  //! Prints information about the factorization and solution phases.
  void PrintStatus() const;

  //! Extracts timing information from the current solver and places it in the parameter list.
  void GetTiming( Teuchos::ParameterList &TimingParameterList ) const { Amesos_Time::GetTiming(TimingParameterList); }
 
  //@}

private:  
  int CheckError(const int error) const;

  inline const Epetra_Map& Map() const
  {
    return(Matrix_->RowMatrixRowMap());
  }
  
  inline const Epetra_RowMatrix& Matrix() const
  {
    return(*Matrix_);
  }

  int ConvertToCssMKL();
  int PerformSymbolicFactorization();
  int PerformNumericFactorization(); 


  Teuchos::RCP<Amesos_StandardIndex> StdIndex_;
  const Epetra_RowMatrix* OriginalMatrix_;
  const Epetra_RowMatrix* Matrix_;

  //! If \c true, the transpose of A is used.
  bool UseTranspose_;
  //! Pointer to the linear system problem.
  const Epetra_LinearProblem* Problem_;

  //! Quick access pointers to the internal timing data.
  int MtxConvTime_, MtxRedistTime_, VecRedistTime_;
  int SymFactTime_, NumFactTime_, SolveTime_;

  // Data for CSSMKL
  std::vector<double> aa_;
  std::vector<int>    ia_;
  std::vector<int>    ja_;
  std::vector<int>    perm_;

  void* pt_[64];

  int iparm_[64];
  double dparm_[64];
  int mtype_;
  int maxfct_; // Maximal number of factors with idential nonzero pattern (always 1)
  int mnum_; //! Actual matrix for solution phase (always 1)
  int msglvl_; //! Output level
  int nrhs_; // Number of RHS

  Teuchos::ParameterList param_;

  MPI_Fint CssComm_;
};  // class Amesos_CssMKL
#endif
