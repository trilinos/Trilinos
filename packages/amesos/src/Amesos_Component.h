/*
 I think that Amesos_Component should be an additional interface and
 hence functions which do not differ from the Amesos_BaseSolver class 
 are not included here.
 */
/*
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
*/

#ifndef _AMESOS_COMPONENTBASESOLVER_H_
#define _AMESOS_COMPONENTBASESOLVER_H_

#include "Teuchos_ParameterList.hpp"
#include "Epetra_LinearProblem.h"
class Epetra_LinearProblem;
class Epetra_MultiVector;
class Epetra_Map;
class Epetra_Comm;

//! Amesos_Component: A pure virtual class for direct solvers to be used within Amesos_Merikos to form a parallel direct solver.

/*! The Amesos_Component class is a pure virtual class specifying what Amesos_Merikos needs.   Every Amesos class named 
    Amesos_Comp_<i>SolverName</i> implements Amesos_Component.


    <H1>Differences betweeen Amesos_BaseSolver and  Amesos_Component.</H1> 

    <ul>
    <li>Amesos_Component performs minimal reordering.  
    <ul>
      <li>Amesos_Component expects row and coloumn ordering on input.
      <li>Amesos_Component performs LUx=b (not PLUQx = b).
    </ul>
    <li>Amesos_Component performs the forward solve (Lsolve) and 
the backward solve (Usolve) separately.
    <li>Amesos_Component performs a partial factorization.  
    <ul>
       <li>Amesos_Component performs factors at most the 
first SubMatrixSize_ rows and columns.  
       <li>Amesos_Component delays the factorization of any columns which generate zero pivots. 
       <li>Amesos_Component computes and returns the schur complement for the 
       <li>Lsolve and Usolve perform partial solves, treating the unfactored rows and columns as the identity.
    </ul>
    </ul>


    <H1>Usage Examples</H1> 

    <H2>Basic calling sequence</H2> 

<pre>
    Epetra_LinearProblem Problem(A,X,B);
    Amesos_SolverName Solver(Problem);

    Solver.PartialFactorization() ; 
      ... Ancestor factorization
    Solver.Lsolve() ; 
      ... Ancestor solves
    Solver.Usolve() ; 
</pre>

    <H2>Preconditions:
<ul>
<li>An ordering 
</ul>
    <H2>Postconditions:

    <H2>Constructor requirements</H2>
    Every Amesos_SolverName class should accept an
    Epetra_LinearProblem 


*/    

class Amesos_Component {
      
 public:

  //@{ \name Destructor.
    //! Destructor
    virtual ~Amesos_Component() {};
  //@}
  
  //@{ \name Mathematical functions.


    //! Performs partial factorization on the matrix A.
    /*! 
      Partial Factorization perfom

     \return Integer error code, set to 0 if successful.
  */
    virtual int PartialFactorization() = 0;

    //! Solves L X = B (or L<SUP>T</SUP> x = B) 
    /*! 


     \return Integer error code, set to 0 if successful.
  */
    virtual int Lsolve() = 0;

    //! Solves L X = B (or L<SUP>T</SUP> x = B) 
    /*! 


     \return Integer error code, set to 0 if successful.
  */
    virtual int Usolve() = 0;

  //@}
  
  //@{ \name Atribute access functions

    //! SetRowPermutation
    virtual int SetRowPermutation( int* RowPermutation ) = 0;

    //! SetColumnPermutation
    virtual int SetColumnPermutation( int* ColumnPermutation ) = 0;

    //! SetSubMatrixSize
    virtual int SetSubMatrixSize( int SubMatrixSize ) = 0;

    //! GetRowPermutation
    /*!
      RowPermutation reflects any row permutations performed by 
      PartialFactorization(). 
      Note:  It is not yet clear whether this row permutation 
      includes the RowPermuation upon input or whether it returns
      only the row permuations performed by the most recent 
      call to PartialFactorization().  In other words, in the 
      absence of pivoting, RowPerumation might be identical to 
      that given by SetRowPermutation() or it might be the 
      identity permutation.
     */
    virtual int GetRowPermutation( int** RowPermutation ) = 0;

    //! GetColumnPermutation
    /*!
      ColumnPermutation reflects any row permutations performed by 
      PartialFactorization(). 
      Note:  It is not yet clear whether this row permutation 
      includes the ColumnPermuation upon input or whether it returns
      only the row permuations performed by the most recent 
      call to PartialFactorization().  In other words, in the 
      absence of pivoting, ColumnPerumation might be identical to 
      that given by SetColumnPermutation() or it might be the 
      identity permutation.
     */
    virtual int GetColumnPermutation( int** ColumnPermutation ) = 0;

    //! GetSubMatrixSize
    /* 
       SubMatrixSize is the number of rows and columns in the matrix
       that was factored.  (i.e. the number of columns of L and the
       number of rows of U)
    */
    virtual int GetSubMatrixSize( int* SubMatrixSize ) = 0;

    //! GetSchurComplement
    /*
      SchurComplement is a square matrix with each side having size
      MatrixSize-SubMatrixSize which contains the Schur 
      complement based on the matrices L and U, i.e.
        U(1:SubMatrixSize,SubMatrixSize+1:MatrixSize)^T *
        L(SubMatrixSize+1:MatrixSize,1:SubMatrixSize)
     */ 
    virtual int GetSchurComplement( Epetra_CrsMatrix* SchurComplement ) = 0;

  //@}

};

#endif /* _AMESOS_COMPONENTBASESOLVER_H_ */
