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
#include "Epetra_LinearProblem.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_Comm.h"
#endif
#include "Epetra_CrsGraph.h"

class SLUData;

//! Amesos_Superlu:  Amesos interface to Xioye Li's SuperLU serial code.  
/*!  Amesos_Superlu, an object-oriented wrapper for Superlu, will solve a linear systems of equations: <TT>A X = B</TT>
   using Epetra objects and the Superlu solver library, where
  <TT>A</TT> is an Epetra_RowMatrix and <TT>X</TT> and <TT>B</TT> are 
  Epetra_MultiVector objects.
    
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
  /*! Completely deletes an Amesos_Superlu object.  
  */
  ~Amesos_Superlu(void);
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

    //! Solves A X = B (or A<SUP>T</SUP> X = B) 
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

  //! Returns true if SUPERLU can handle this matrix shape 
  /*! Returns true if the matrix shape is one that SUPERLU can
    handle. SUPERLU only works with square matrices.  
  */
  bool MatrixShapeOK() const ;

  //! SetUseTranpose(true) is more efficient in Amesos_Superlu
  /*! 
<ul>
  <li>If SetUseTranspose() is set to true, 
    <ul>
       <li><p class="code">A<sup>T</sup> X = B</p> is computed</li>
       <li>(This is the more efficient operation)</li>
    </ul></li>
  <li>else
    <ul>
       <li><p class="code">A X = B</p> is computed</li>
       <li>(This requires a matrix transpose)</li>
    </ul></li>
</ul>
  */  
  int SetUseTranspose(bool UseTranspose) {UseTranspose_ = UseTranspose; return(0);};

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(UseTranspose_);};

  //! Returns a pointer to the Epetra_Comm communicator associated with this matrix.
  const Epetra_Comm & Comm() const {return(GetProblem()->GetOperator()->Comm());};

  //!  Updates internal variables. 
  /*!  
      <br \>Preconditions:<ul>
      <li>None.</li>
      </ul>

      <br \>Postconditions:<ul> 
      <li>Internal variables controlling the factorization and solve will
      be updated and take effect on all subsequent calls to NumericFactorization() 
      and Solve().</li>
      <li>All parameters whose value are to differ from the default values must 
be included in ParameterList.  Parameters not specified in ParameterList 
revert to their default values.
      </ul>

    Amesos_Superlu accepts the following parameters:
    <ul>
      <li>"Verbose" - boolean:false - If true prints out some debug information
      <li>"Superlu" - list containing the following parameters:
      <ul>
        <li>"FactOption" - string:["SamePattern"] "SamePattern-SameRowPerm"
	<li>"ColPerm" - string:["COLAMD"] "..."
	<li>"Equil" - boolean:true 
        <li>"fill_fac" - int:-1
        <li>"panel_size" - int:-1
        <li>"relax" - int:-1
	<li>"pivot_thresh" - double:-1
	<li>"Iter_Refine" - string:["DOUBLE"] "NOREFINE" 
      </ul>
    </ul>

    
    \return Integer error code, set to 0 if successful. 

   */
  int SetParameters( Teuchos::ParameterList &ParameterList )  ;

  //@}

 private:  

  int Factor();
  int ReFactor();

  /*
  ConvertToSerial - Convert matrix to a serial Epetra_CrsMatrix
    Preconditions:
      Problem_ must be set 
      SerialMap and SerialCrsMatrix must either be 0 or be pointers to 
        appropriatly allocate objects.  If they are non-zero, those objects
	will be deleted (and possibly recreated).  
	
    Postconditions:
      IsLocal is set to 1 if the input matrix is entirely stored on process 0
      SerialMap points to a serial map if IsLocal==1
      SerialCrsMatrix contains a serial version of the matrix A if IsLocal==1
      SerialMatrix points to a serial copy of the matrix
      NumGlobalElements_   is set to the number of rows in the matrix
      numentries_ is set to the number of non-zeroes in the matrix 
   */
  int ConvertToSerial();

  /*
    PerformNumericFactorization - Call Superlu to perform numeric factorization
    Preconditions:
      IsLocal must be set 
      Ap, Ai and Aval are a compressed row storage version of the input matrix A.
      Symbolic must be set
    Postconditions:
      Numeric points to an SUPERLU internal opaque object containing the
        numeric factorization and accompanying information.  
      NumericFactorizationOK_ = true; 
    Note:  All action is performed on process 0
  */
  int PerformNumericFactorization(); 

 protected:

  SLUData * data_;
  double *DummyArray;


  //
  //  Ap, Ai, Aval form the compressed row storage used by Superlu
  //
  //  #define NOVEC
#ifdef NOVEC
  int* Ap_;
  int* Ai_;
  double* Aval_;
#else
  vector <int> Ap_;
  vector <int> Ai_;
  vector <double> Aval_;
#endif

  bool FactorizationDone_ ; 
  bool FactorizationOK_ ; 
  bool ReuseSymbolic_ ;
  bool UseTranspose_;      

  int iam_;                 //  Process number (i.e. Comm().MyPID() 
  
  int IsLocal_;            //  1 if Problem_->GetOperator() is stored entirely on process 0
                           //  Note:  Local Problems do not require redistribution of
                           //  the matrix A or vectors X and B.
  int numentries_;         //  Number of non-zero entries in Problem_->GetOperator()
  int NumGlobalElements_;  //  Number of rows and columns in the Problem_->GetOperator()

  Epetra_Map *SerialMap_ ;               //  Points to a Serial Map (unused if IsLocal == 1 ) 
  Epetra_CrsMatrix *SerialCrsMatrixA_ ;  //  Points to a Serial Copy of A (unused if IsLocal==1)
  Epetra_CrsMatrix *SerialMatrix_ ;      //  Points to a Serial Copy of A 
                                         //  IsLocal==1 - Points to the original matix 
                                         //  IsLocal==0 - Points to SerialCrsMatrixA
  const Epetra_LinearProblem * Problem_;

#ifdef NOVEC
  int* ColIndicesV_;
  double* RowValuesV_;
  //  int* Global_Columns_; 
  double* berr_;
  double* ferr_;

  int* perm_r_;
  int* perm_c_;
  int* etree_;
  double* R_;
  double* C_;
#else
  vector<int> ColIndicesV_;
  vector<double> RowValuesV_;
  //  vector<int> Global_Columns_; 
  vector<double> berr_;
  vector<double> ferr_;

  vector<int> perm_r_;
  vector<int> perm_c_;
  vector<int> etree_;
  vector<double> R_;
  vector<double> C_;
#endif
  char equed_;
  bool DestroyBandX_ ; 

};  // End of  class Amesos_Superlu  
#endif /* _AMESOS_SUPERLU_H_ */
