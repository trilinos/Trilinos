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

#ifndef _AMESOS_SCALAPACK_H_
#define _AMESOS_SCALAPACK_H_

#include "Amesos_SCALAPACK_wrappers.h"



#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
#include "Epetra_LinearProblem.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_Comm.h"
#endif
#include "Epetra_CrsGraph.h"

//! Amesos_Scalapck:  A parallel dense solver.  For now, we implement only the unsymmetric ScaLAPACK solver.
/*!  Amesos_ScaLAPACK, an object-oriented wrapper for ScaLAPACK, will solve a linear systems of equations: <TT>A X = B</TT>
   using Epetra objects and the ScaLAPACK library, where
  <TT>A</TT> is an Epetra_RowMatrix and <TT>X</TT> and <TT>B</TT> are 
  Epetra_MultiVector objects.

<br /><br /><p>AmesosScaLAPACK can be competitive for matrices 
that are not particularly sparse.  ScaLAPACK solves matrices
for which the fill-in is roughly 10% to 20% of the matrix size 
in time comparable to that achieve by other Amesos classes.

*/

class Amesos_Scalapack: public Amesos_BaseSolver { 

public: 

  //@{ \name Constructor methods
  //! Amesos_Scalapack Constructor.
  /*! Creates an Amesos_Scalapack instance, using an Epetra_LinearProblem,
      passing in an already-defined Epetra_LinearProblem object. 

      Note: The operator in LinearProblem must be an
      Epetra_RowMatrix.

  */
  Amesos_Scalapack(const Epetra_LinearProblem& LinearProblem, const AMESOS::Parameter::List &ParameterList );

  //! Amesos_Scalapack Destructor.
  /*! Completely deletes an Amesos_Scalapack object.  
  */
  ~Amesos_Scalapack(void);
  //@}

  //@{ \name Mathematical functions.

    //! Performs SymbolicFactorization on the matrix A.
    /*! 
      There is no symbolic factorization phase in ScaLAPACK, as it operates
      only on dense matrices.  Hence, Amesos_Scalapack::SymbolicFactorization()
      takes no action.

    \return Integer error code, set to 0 if successful.
  */
    int SymbolicFactorization() ;

    //! Performs NumericFactorization on the matrix A.
    /*!  In addition to performing numeric factorization 
      on the matrix A, the call to
      NumericFactorization() implies that no change will be made to
      the underlying matrix without a subsequent call to
      NumericFactorization().  

      preconditions:<ul>
      <li>GetProblem().GetOperator() != 0 (return -1)
      <li>MatrixShapeOk(GetProblem().GetOperator()) == true (return -6)
      <li>The non-zero structure of the matrix should not have changed
      since the last call to SymbolicFactorization().  Irrelevant for
      Amesos_Scalapack.
      <li>The distribution of the matrix should not have changed 
      since the last call to SymbolicFactorization(). Irrelevant for
      Amesos_Scalapack.
      </ul>

      postconditions:<ul>
      <li>nprow_, npcol_, DescA_ 
      <li>DenseA will be factored
      <li>Ipiv_ contains the pivots
      </ul>

     \return Integer error code, set to 0 if successful.
  */
    int NumericFactorization() ;

    //! Solves A X = B (or A<SUP>T</SUP> X = B) 
    /*! 

      preconditions:<ul>
      <li>GetProblem().GetOperator() != 0 (return -1)
      <li>MatrixShapeOk(GetProblem().GetOperator()) == true (return -6)
      <li>GetProblem()->CheckInput (see Epetra_LinearProblem::CheckInput() for return values)
      <li>The matrix should not have changed
          since the last call to NumericFactorization().
      </ul>

      postconditions:<ul> 
      <li>X will be set such that A X = B (or
      A<SUP>T</SUP> X = B), within the limits of the accuracy of the
      the Scalapack solver.  
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

  //! Returns true if SCALAPACK can handle this matrix shape 
  /*! Returns true if the matrix shape is one that SCALAPACK can
    handle. SCALAPACK only works with square matrices.  
  */
  bool MatrixShapeOK() const ;

  //! SetUseTranpose(true) is more efficient in Amesos_Scalapack
  /*! 
<ul>
  <li>If SetUseTranspose() is set to true, 
    <ul>
       <li><p class="code">A<sup>T</sup> X = B</p> is computed</li>
    </ul></li>
  <li>else
    <ul>
       <li><p class="code">A X = B</p> is computed</li>
    </ul></li>
</ul>
  */  
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
  /*
  RedistributeA - Convert matrix to a dense ScaLAPACK matrix
    Preconditions:
      Problem_ must be set 
      ReadParameterList() 
      ScaLAPACK1DMap and ScaLAPACK1DMatrix must either be 0 or be pointers to 
        appropriatly allocate objects.  If they are non-zero, those objects
	will be deleted and recreated.  
	
    Postconditions:
      nprow_, npcol_, DescA_ 
      RowMatrixA_ 
      ScaLAPACK1DMap_ 
      ScaLAPACK1DMatrix_ 
      ImportToScaLAPACK1D_
      ImportBackToOriginal_

   */
  int RedistributeA();

  /*
    ConvertToScalapack - Convert matirx to form expected by Scalapack: Ai, Ap, Aval
    Preconditions:
      Problem_ 
    Postconditions:
      nprow_, npcol_, 
  */
  int ConvertToScalapack();     

  /*
    PerformNumericFactorization - Call Scalapack to perform numeric factorization
    Preconditions:

    Postconditions:
      DenseA_, DescA_
  */
  int PerformNumericFactorization(); 

 protected:

  int MaxProcesses_;                     // default is -1 ; If positive, distribute 
                                         // problem over MaxProcesses

  int iam_;                              //  Process number (i.e. Comm().MyPID() 
  
  int NumGlobalElements_;                //  Number of rows and columns in the Problem_->GetOperator()

  //
  //  The following variables are required for the ScaLAPACK interface:
  //
  int nprow_ ;                           //  number of process rows: 1 for now
  int npcol_ ;                           //  number of process columns
  int ictxt_ ;                           //  BLACS context
  int m_per_p_;                          //  Number of columns per process
  int DescA_[10];                        //  ScaLAPACK array descriptor 

  Epetra_Map *ScaLAPACK1DMap_ ;          //  Points to a 1D Map which matches a ScaLAPACK 1D
                                         //  blocked (not block cyclic) distribution
  Epetra_CrsMatrix *ScaLAPACK1DMatrix_ ; //  Points to a  ScaLAPACK 1D
                                         //  blocked (not block cyclic) distribution
  Epetra_Map *SerialMap_ ;               //  Points to a 1D Map which matches a ScaLAPACK 1D
                                         //  blocked (not block cyclic) distribution
  vector<double> DenseA_;                //  The data in a ScaLAPACK 1D blocked format
  vector<int> Ipiv_ ;                    //  ScaLAPACK pivot information
 

  bool UseTranspose_;     
  const Epetra_LinearProblem * Problem_;
  const AMESOS::Parameter::List * ParameterList_ ; 

};  // End of  class Amesos_Scalapack  
#endif /* _AMESOS_SCALAPACK_H_ */
