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

/*! 
    <p>The Amesos_Component interface specifies what Amesos_Merikos needs.
Any Amesos class that implements Amesos_Component can be used by 
Amesos_Merikos to perform partial solves on subblocks of the matrix. 


    <H1>Member functions added by Amesos_Component.</H1> 

    <ul>
    <li>PartialFactorization()
    <ul>
       <li>PartialFactorization performs factors at most the 
first SubMatrixSize_ rows and columns.  
       <li>PartialFactorization delays the factorization of any columns which generate unstable (i.e. too small) pivots. 
       <li>PartialFactorization computes and returns the schur complement.
       <li>PartialFactorization does not need a symbolic factorization phase. 
It uses the permutation given by SetRowPermutation.
    </ul>
    <li>Lsolve performs a raw partial solve, treating the unfactored rows and 
columns as the identity without row or column permutation.  
    <li>Usolve performs a raw partial solve, treating the unfactored rows and 
columns as the identity without row or column permutation.  
    <li>SetRowPermutation - sets the row permutation
    <li>GetRowPermutation - gets the row permutation
    <li>SetColumnPermutation - sets the column permutation
    <li>GetColumnPermutation - gets the column permutation
    <li>SetSubMatrixSize - Sets the maximum number of rows (and columns)
to factor.
    <li>GetSubMatrixSize - Returns the number of rows (and columns) 
actually factored. 
    <li>SchurComplement - Returns the Schur complement, i.e. 
L21(SubMatrixSize+1:MatrixSize,1:SubMatrixSize) *
U12(1:SubMatrixSize,SubMatrixSize+1:MatrixSize)
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
      absence of pivoting, RowPermutation might be identical to 
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
      absence of pivoting, ColumnPermutation might be identical to 
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
        L(SubMatrixSize+1:MatrixSize,1:SubMatrixSize) *
        U(1:SubMatrixSize,SubMatrixSize+1:MatrixSize)
     */ 
    virtual int GetSchurComplement( Epetra_CrsMatrix* SchurComplement ) = 0;

  //@}

};

#endif /* _AMESOS_COMPONENTBASESOLVER_H_ */
