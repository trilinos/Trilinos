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

#ifndef AMESOS_DSCPACK_H
#define AMESOS_DSCPACK_H

#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_MpiComm.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Import.h"

//
//  dscmain.h does not check to make sure that it is not called twice,
//  hence the following check:
//
#ifndef DSC_LBLAS1
#define DBL_R_NUM
extern "C" {
#include "dscmain.h"
}
#endif

//! Amesos_Dscpack:  An object-oriented wrapper for Dscpack.
/*!  Amesos_Dscpack will solve a linear systems of equations: <TT>A X = B</TT>
   using Epetra objects and the Dscpack solver library, where
  <TT>A</TT> is an Epetra_RowMatrix and <TT>X</TT> and <TT>B</TT> are 
  Epetra_MultiVector objects.

*/
class Amesos_Dscpack: public Amesos_BaseSolver { 

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

    //! Solves A X = B (or A<SUP>T</SUP> x = B) 
    /*!
     \return Integer error code, set to 0 if successful.
  */
    int Solve();

  //@}
  
  //@{ \name Additional methods required to support the Epetra_Operator interface.

  //! Get a pointer to the Problem.
  const Epetra_LinearProblem *GetProblem() const { return(Problem_); };

  //! Returns true if DSCPACK can handle this matrix shape 
  /*! Returns true if the matrix shape is one that DSCPACK can
    handle. DSCPACK only works with symetric matrices.  
  */
  bool MatrixShapeOK() const ;

  int SetUseTranspose(bool UseTranspose) {UseTranspose_ = UseTranspose; return(0);};

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(UseTranspose_);};

  //! Returns a pointer to the Epetra_Comm communicator associated with this matrix.
  const Epetra_Comm & Comm() const {return(GetProblem()->GetOperator()->Comm());};

  //! Sets parameters as specified in the list, returns 0 if successful.
  int SetParameters( Teuchos::ParameterList &ParameterList )  ;

  //! Prints timing information
  void PrintTiming();
  
  //! Prints information about the factorization and solution phases.
  void PrintStatus();  

  //@}

protected:  
  
  //! Performs the symbolic factorization.
  int PerformSymbolicFactorization();
  //! Performs the numeric factorization.
  int PerformNumericFactorization();

  //! If \c true, SymbolicFactorization() has been successfully called.
  bool IsSymbolicFactorizationOK_; 
  //! If \c true, NumericFactorization() has been successfully called.
  bool IsNumericFactorizationOK_; 

  //! Distribution specified by DscOrder
  Epetra_CrsGraph * DscGraph_;

  //! Is \c true, the transpose of the matrix is used.
  bool UseTranspose_;
  //! Pointer to the linear problem.
  const Epetra_LinearProblem * Problem_;

  DSC_Solver	MyDSCObject;
  MPI_Comm MPIC ; 

  bool FirstCallToSolve_;
  //! Tells us whether to free them
  bool A_and_LU_built;
  int *GlobalStructNewColNum; 
  int *GlobalStructNewNum; 
  int *GlobalStructOwner;
  int *LocalStructOldNum;

  int MyDscRank ; 
  int DscNumProcs ; 
  int NumLocalCols; 
  int NumGlobalCols;
  int NumLocalStructs;
  int NumLocalNonz ; 

  //! If \c true, prints timing information in the destructor.
  bool PrintTiming_;
  //! If \c true, prints additinal information in the destructor.
  bool PrintStatus_;
  //! If \c true, computes the norm of rhs and solution.
  bool ComputeVectorNorms_;
  //! If \c true, compute the norm of the real residual.
  bool ComputeTrueResidual_;
  //! Toggles the output level.
  int verbose_;
  //! time to convert to DSCPACK format
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
  //! number of symbolic factorizations
  int NumSymbolicFact_;
  //! number of numeric factorizations
  int NumNumericFact_;
  //! number of solves
  int NumSolve_;
  //! used to track timing
  Epetra_Time * Time_;

  Epetra_Import * ImportToSerial_;

  Epetra_Map * DscMap_;

  int MaxProcs_;
  
  // track memory (as reported by DSCPACK routines)
  int TotalMemory_;                       // estimates of the total memory requirements
                                          // for the factorization step as a
					  // whole number of Mbytes. As
					  // reported in the manual, this is a
					  // "fair" estimation, but not
					  // accurate at the last byte.

  
};  // class Amesos_Dscpack  
#endif /* AMESOS_DSCPACK_H */
