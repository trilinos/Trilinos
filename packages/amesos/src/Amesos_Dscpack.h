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

#ifndef _AMESOS_DSCPACK_H_
#define _AMESOS_DSCPACK_H_

#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_MpiComm.h"
#include "Epetra_CrsGraph.h"

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

  Dscpack execution can be tuned through a variety of parameters.
  Amesos_Dscpack.h allows control of these parameters through the
  following named parameters, ignoring parameters with names that it
  does not recognize.  Where possible, the parameters are common to
  all direct solvers (although some may ignore them).  However, some
  parameters, in particular tuning parameters, are unique to each
  solver.
    
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
  Amesos_Dscpack(const Epetra_LinearProblem& LinearProblem, const AMESOS::Parameter::List &ParameterList );

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

      preconditions:<ul>
      <li>GetProblem().GetOperator() != 0 (return -1)
      <li>MatrixShapeOk(GetProblem().GetOperator()) == true (return -6)
      <li>The non-zero structure of the matrix should not have changed
          since the last call to SymbolicFactorization().  
      <li>The distribution of the matrix should not have changed 
          since the last call to SymbolicFactorization()
      </ul>

      postconditions:<ul>
      <li>Numeric Factorization will be performed (or marked to be performed) 
      allowing Solve() to be performed correctly despite a potential change in 
      in the matrix values (though not in the non-zero structure).
      </ul>

     \return Integer error code, set to 0 if successful.
  */
    int NumericFactorization() ;

    //! Solves A X = B (or A<SUP>T</SUP> x = B) 
    /*! 

      preconditions:<ul>
      <li>GetProblem().GetOperator() != 0 (return -1)
      <li>MatrixShapeOk(GetProblem().GetOperator()) == true (return -6)
      <li>GetProblem()->CheckInput (see Epetra_LinearProblem::CheckInput() for return values)
      <li>The non-zero structure of the matrix should not have changed
          since the last call to SymbolicFactorization().
      <li>The distribution of the matrix should not have changed 
          since the last call to SymbolicFactorization()
      <li>The matrix should not have changed
          since the last call to NumericFactorization().
      </ul>

      postconditions:<ul> 
      <li>X will be set such that A X = B (or
      A<SUP>T</SUP> X = B), within the limits of the accuracy of the
      underlying solver.  
      </ul>

     \return Integer error code, set to 0 if successful.
  */
    int Solve();

  //@}
  
  //@{ \name Additional methods required to support the Epetra_Operator interface.

#if 0
  //! Returns a character string describing the operator
  char * Label() const {return(Epetra_Object::Label());};
#endif
    
  //! Get a pointer to the Problem.
  const Epetra_LinearProblem *GetProblem() const { return(Problem_); };

  //! Get a pointer to the ParameterList.
  const AMESOS::Parameter::List *GetParameterList() const { return(ParameterList_); };

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

  //! Reads the parameter list and updates internal variables. 
  /*!
    ReadParameterList is called by SymbolicFactorization.  Hence, few codes 
    will need to make an explicit call to ReadParameterList.
   */
  int ReadParameterList() ;
  //@}

 private:  
  int PerformSymbolicFactorization();
  int PerformNumericFactorization();

 protected:

  bool SymbolicFactorizationOK_; 
  bool NumericFactorizationOK_; 

  Epetra_CrsGraph * DscGraph_ ; //  Distribution specified by DscOrder

  bool UseTranspose_;
  const Epetra_LinearProblem * Problem_;
  const AMESOS::Parameter::List * ParameterList_ ; 

  DSC_Solver	MyDSCObject;
  MPI_Comm MPIC ; 

  bool FirstCallToSolve_;
  bool A_and_LU_built ;            // Tells us whether to free them 
  int *GlobalStructNewColNum ; 
  int *GlobalStructNewNum ;  
  int *GlobalStructOwner ; 
  int *LocalStructOldNum ; 

  int MyDscRank ; 
  int DscNumProcs ; 
  int NumLocalCols; 
  int NumGlobalCols;
  int NumLocalStructs;
  int NumLocalNonz ; 



};  // End of  class Amesos_Dscpack  
#endif /* _AMESOS_DSCPACK_H_ */
