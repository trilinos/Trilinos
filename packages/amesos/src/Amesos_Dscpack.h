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

#ifndef AMESOS_DSCPACK_H
#define AMESOS_DSCPACK_H

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
#include "Epetra_LinearProblem.h"
#include "Epetra_Comm.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"

// Amesos_Dscpack_Pimpl contains a pointer to the structure defined in 
// dscmain.h.  This prevents Amesos_Dscpack.h from having to include dscmain.h.
//
//  Doxygen does not handle forward class references well.
#ifndef DOXYGEN_SHOULD_SKIP_THIS
class Amesos_Dscpack_Pimpl ; 
#endif


//! Amesos_Dscpack:  An object-oriented wrapper for Dscpack.
/*!  Amesos_Dscpack will solve a linear systems of equations: <TT>A X = B</TT>
   using Epetra objects and the Dscpack solver library, where
  <TT>A</TT> is an Epetra_RowMatrix and <TT>X</TT> and <TT>B</TT> are 
  Epetra_MultiVector objects.

*/
class Amesos_Dscpack: public Amesos_BaseSolver, 
                      private Amesos_Time, 
                      private Amesos_NoCopiable, 
                      private Amesos_Utils, 
                      private Amesos_Control, 
                      private Amesos_Status { 

public: 

  //@{ \name Constructor methods
  //! Amesos_Dscpack Constructor.
  /*! Creates an Amesos_Dscpack instance, using an Epetra_LinearProblem,
      passing in an already-defined Epetra_LinearProblem object. 

      Note: The operator in LinearProblem must be an
      Epetra_RowMatrix.

  */
  Amesos_Dscpack(const Epetra_LinearProblem& LinearProblem );

  //! Amesos_Dscpack Destructor.
  /*! Completely deletes an Amesos_Dscpack object.  
  */
  ~Amesos_Dscpack(void);
  //@}

  //@{ \name Mathematical functions.

  int SymbolicFactorization() ;

  int NumericFactorization() ;

  int Solve();

  //@}
  //@{ \name Additional methods required to support the Epetra_Operator interface.

  const Epetra_LinearProblem *GetProblem() const { return(Problem_); };

  //! Returns true if DSCPACK can handle this matrix shape 
  /*! Returns true if the matrix shape is one that DSCPACK can
    handle. DSCPACK only works with symetric matrices.  
  */
  bool MatrixShapeOK() const ;

  int SetUseTranspose(bool UseTranspose) 
  {
    return(0);
  }

  bool UseTranspose() const 
  {
    return(false);
  }

  const Epetra_Comm& Comm() const {return(GetProblem()->GetOperator()->Comm());};

  int SetParameters( Teuchos::ParameterList &ParameterList )  ;

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
  
  const Epetra_Import& Importer() const
  {
    return(*Importer_.get());
  }

  const Epetra_Map& DscRowMap() const
  {
    return(*DscRowMap_.get());
  }

  const Epetra_Map& DscColMap() const
  {
    return(*DscColMap_.get());
  }

  //! Performs the symbolic factorization.
  int PerformSymbolicFactorization();

  //! Performs the numeric factorization.
  int PerformNumericFactorization();

  //! Pointer to the linear problem.
  const Epetra_LinearProblem * Problem_;

  Teuchos::RCP<Amesos_Dscpack_Pimpl> PrivateDscpackData_; 

  bool FirstCallToSolve_;
  //! Tells us whether to free them
  bool A_and_LU_built;
  int *GlobalStructNewColNum; 
  int *GlobalStructNewNum; 
  int *GlobalStructOwner;
  int *LocalStructOldNum;

  int MyDscRank;
  int DscNumProcs;
  int NumLocalCols;
  int NumGlobalCols;
  int NumLocalStructs;
  int NumLocalNonz ; 

  RCP<Epetra_Import> Importer_;
  RCP<Epetra_Map>    DscColMap_;
  RCP<Epetra_Map>    DscRowMap_;

  int MaxProcs_;
 
  int MtxRedistTime_, MtxConvTime_, VecRedistTime_;
  int SymFactTime_, NumFactTime_, SolveTime_, OverheadTime_;
  
  // track memory (as reported by DSCPACK routines)
  int TotalMemory_;                       // estimates of the total memory requirements
                                          // for the factorization step as a
					  // whole number of Mbytes. As
					  // reported in the manual, this is a
					  // "fair" estimation, but not
					  // accurate at the last byte.

  
};  // class Amesos_Dscpack  
#endif /* AMESOS_DSCPACK_H */
