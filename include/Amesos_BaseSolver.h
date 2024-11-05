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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

/*!
 * \file Amesos_BaseSolver.h
 *
 * \class Amesos_BaseSolver
 *
 * \brief Pure virtual class for all Amesos concrete implementions
 *
 * \date Last updated on 24-May-05.
 */

#ifndef _AMESOS_BASESOLVER_H_
#define _AMESOS_BASESOLVER_H_

#if defined(Amesos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Amesos package is deprecated"
#endif
#endif

const int StructurallySingularMatrixError = -21;
const int NumericallySingularMatrixError = -22;

//#include "Amesos_ConfigDefs.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Epetra_LinearProblem.h"
class Epetra_LinearProblem;
class Epetra_MultiVector;
class Epetra_Map;
class Epetra_Comm;

//! Amesos_BaseSolver: A pure virtual class for direct solution of real-valued double-precision operators.
/*! 

The Amesos_BaseSolver class is a pure virtual class (that is, it specifies
interface only) that enables the use of real-valued
double-precision direct sparse solvers.  Every Amesos class named 
Amesos_<i>SolverName</i> derives from Amesos_BaseSolver.

<H1>Usage Examples</H1> 

<H2>Basic calling sequence</H2> The basic
calling sequence solves A x = b or A<SUP>T</SUP> x = b without
specifying how A has changed between each call to Solve().

\code
    Epetra_LinearProblem Problem(A,X,B);
    Amesos_SolverName Solver(Problem);
    Problem.SetUseTranspose(false); 
    while( ... ) { 

      <Here code which may change A or B>

      Solver.SymbolicFactorization(); 
      Solver.NumericFactorization(); 
      Solver.Solve();
    }
\endcode   
    

<H2>Re-using the symbolic factorization</H2> The following 
calling sequence performs multiple solves of A x = b or A<SUP>T</SUP> x = b 
in cases where the non-zero structure of A remains unchanged
between each call to Solve().

\code
    Epetra_LinearProblem Problem(A,X,B);
    Amesos_SolverName Solver(Problem);
    Problem.SetUseTranspose(false); 
    Solver.SymbolicFactorization(); 
    while( ... ) { 

      <Here code which may change B or the non-zero values
      (but not the non-zero structure) of A>

      Solver.NumericFactorization() ; 
      Solver.Solve() ; 
    }
\endcode
    


<H2>Re-using the numeric factorization</H2> The following 
calling sequence performs multiple solves of A x = b or A<SUP>T</SUP> x = b 
provided that A remains unchanged between each call to Solve().

\code
    Epetra_LinearProblem Problem(A,X,B);
    Amesos_SolverName Solver(Problem);
    Problem.SetUseTranspose( false ); 
    Solver.NumericFactorization(); 
    while( ... ) { 

      <Here code which may change B but not A>

      Solver.Solve();
    }
\endcode
    

<H2>Constructor requirements</H2>
 Every Amesos_SolverName class should accept an Epetra_LinearProblem 


<H2>Mathematical methods</H2> Four mathematical methods are
defined in the base class Amesos_BaseSolver:
SymbolicFactorization(), NumericFactorization(),
and Solve().  


<H2>Switching concrete classes</H2>
Different concrete
classes, each based on a different third party solver, will have
different performance characteristics and will accept different
parameters.


<H2>Changing the values of the underlying matrix operator.</H2> 

<p>Any changes to the values of a matrix must be accompanied by a call to
NumericFactorization() before the next call to Solve() or the behavior
of Solve() is undefined.  Any changes to the numerical structure of
the matrix must be followed by a call to SymbolicFactorization() and
NumericalFactorization() before the next call to Solve().

<p>Once SymbolicFactorization() has been called, classes implementing
this interface may assume that any change made to the non-zero
structure of the underlying matrix will be accompanied by a call
to SymbolicFactorization() prior to a subsequent call to
NumericFactorization or Solve().


<H2>Named Parameters</H2>
    
Parameters can be changed or added at any time by calling 
SetParameters(ParamList) with the new parameters specified in ParamList.


It is left to the user to be sure that changes made to the
parameters are appropriate for the concrete class that they are
using.

Examples of appropriate changes in parameters include:
<ul>
<li><b>Changing iterative refinement rules between calls to Solve()</b>
<li><b>Changing drop tolerance rules between calls to NumericFactorization()</b>
</ul>

Examples of inappropriate changes in parameters include: <ul>
<li><b>Changing drop tolerance rules between solve steps.</b>
<pre> Solver.NumericFactorization();
Solver.getList()->set("DropTolerance",.001);
Solver.Solve(); </pre> </ul> Results of making inappropriate
changes in parameters is unpredictable and could include an error
return, a bogus result or ignoring the parameter change.
    
<H2>Transpose solve</H2> Any class implementing Amesos_BaseSolver
should handle calls to SetUseTranspose() at any point.  However, the
result of a call to SetUseTranspose() which is not followed by a call
to SymbolicFactorization() and NumericFactorization() is
implementation dependent.  Some third party libraries are able to
solve A<SUP>T</SUP> x = b and Ax = b using the same factorization.
Others will require a new factorization anytime that a call to
SetUseTranspose() changes the intended solve from A<SUP>T</SUP> x = b
to Ax = b or vice-versa.


<H1>Performance expectations</H1>

The following is a list of performance guidelines that classes 
which implement the Amesos_BaseSolver class are expected to maintain.


<H1>Memory usage:</H1>
For serial codes, no more than one extra copy of the original matrix 
should be required.  Except that some codes require matrix transpostion 
which requires additional copies of the input matrix.  

For distributed memory codes, no serial copies of the original matrix 
should be required.  


<H1>Robustness requirements</H1>

<p>Failures should be caught by AMESOS_CHK_ERR().  
   The following error codes should be used:
     - 1: Singular matrix
     - 2: Non-symmetric matrix
     - 3: Matrix is not positive definite
     - 4: Insufficient memory

<p>Because we do not check to see if a matrix has changed
between the call to SymbolicFactorization() and the call to
NumericFactorization(), it is possible that a change to the 
matrix will cause a potentially catastrophic error.  
*/    

class Amesos_BaseSolver 
  : public Teuchos::ParameterListAcceptor
{

#if 0      
 private:
  Teuchos::RCP<Teuchos::ParameterList> paramList_ ; 
#endif

 public:

  /** \brief . */
  using Teuchos::ParameterListAcceptor::getParameterList;

  //@{ \name Destructor
  //! Destructor
  virtual ~Amesos_BaseSolver() {};
  
  //@}
  //@{ \name Mathematical functions

  //! Performs SymbolicFactorization on the matrix A.
  /*! 
      In addition to performing symbolic factorization on the matrix A, 
      the call to SymbolicFactorization() implies that no change will
      be made to the non-zero structure of the underlying matrix without 
      a subsequent call to SymbolicFactorization().
      
      <br \>Preconditions:<ul>
      <li>GetProblem().GetOperator() != 0 (return -1)
      <li>MatrixShapeOk(GetProblem().GetOperator()) == true (return -6)
      </ul>

      <br \>Postconditions:<ul> 
      <li>Symbolic Factorization will be performed (or marked to be performed) 
      allowing NumericFactorization() and Solve() to be called.
      </ul>

    \return Integer error code, set to 0 if successful.
  */
  virtual int SymbolicFactorization() = 0;

  //! Performs NumericFactorization on the matrix A.
  /*!  In addition to performing numeric factorization on the 
      matrix A, the call to
      NumericFactorization() implies that no change will be made to
      the underlying matrix without a subsequent call to
      NumericFactorization().  

      <br \>Preconditions:<ul>
      <li>GetProblem().GetOperator() != 0 (return -1)
      <li>MatrixShapeOk(GetProblem().GetOperator()) == true (return -6)
      <li>The non-zero structure of the matrix should not have changed
          since the last call to SymbolicFactorization().  
          (return -2 if the number of non-zeros changes) 
          Other changes can have arbitrary consequences.
      <li>The distribution of the matrix should not have changed 
          since the last call to SymbolicFactorization()
      <li>The matrix should be indexed from 0 to n-1, unless the parameter "Reindex" 
          was set to "true" prior to the call to SymbolicFactorization().  
          (return -3 - if caught)
      <li>The paremeter "Reindex" should not be set to "true" except on CrsMatrices.
          (return -4) 
      <li>The paremeter "Reindex" should not be set to "true" unless Amesos was built
          with EpetraExt, i.e. with --enable-epetraext on the configure line.
          (return -4) 
      <li>Internal errors retur -5.
      </ul>

      <br \>Postconditions:<ul> 
      <li>Numeric Factorization will be performed (or marked to be performed) 
      allowing Solve() to be performed correctly despite a potential change in 
      in the matrix values (though not in the non-zero structure).
      </ul>

     \return Integer error code, set to 0 if successful.
  */
  virtual int NumericFactorization() = 0;

    //! Solves A X = B (or A<SUP>T</SUP> x = B) 
    /*! 

      <br \>Preconditions:<ul>
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

      <br \>Postconditions:<ul> 
      <li>X will be set such that A X = B (or
      A<SUP>T</SUP> X = B), within the limits of the accuracy of the
      underlying solver.  
      </ul>

     \return Integer error code, set to 0 if successful.
  */
  virtual int Solve() = 0;

  //@}
  //@{ \name Attribute set methods

  //! If set true, X will be set to the solution of A<SUP>T</SUP> X = B (not A X = B) 
  /*! If the implementation of this interface does not support
           transpose use, this method should return a value of -1.
      
      <br \>Preconditions:<ul>
      <li>SetUseTranspose() should be called prior to the call to SymbolicFactorization() If NumericFactorization() or Solve() is called after SetUseTranspose() without an intervening call to SymbolicFactorization() the result is implementation dependent. </li>
      </ul>

      <br \>Postconditions:<ul> 
      <li>The next factorization and solve will be performed with
      the new value of UseTranspose. </li>
      </ul>

    \param UseTranspose -- (In) If true, solve A<SUP>T</SUP> X = B, otherwise
	   solve A X = B.  

    \return Integer error code, set to 0 if successful.  Set to -1 if
           this implementation does not support transpose.
  */
  virtual int SetUseTranspose(bool UseTranspose) = 0;
  
  //! Returns the current UseTranspose setting.
  virtual bool UseTranspose() const = 0;

  //!  Updates internal variables. 
  /*!  
      <br \>Preconditions:<ul>
      <li>None.</li>
      </ul>

      <br \>Postconditions:<ul> 
      <li>Internal variables controlling the factorization and solve will
      be updated and take effect on all subseuent calls to NumericFactorization() 
      and Solve().</li>
      <li>All parameters whose value are to differ from the default values must 
be included in ParameterList.  Parameters not specified in ParameterList 
revert to their default values.
      </ul>

    \return Integer error code, set to 0 if successful. 
  */
  virtual int SetParameters( Teuchos::ParameterList &ParameterList ) = 0 ;

  /** \brief Returns the Epetra_LinearProblem.
   *
   * <b>Warning!</b> Do not call <tt>return->SetOperator(...)</tt> to attempt
   * to change the <tt>Epetra_Operator</tt> object (even if the new matrix has
   * the same structure).  This new operator matrix will be ignored!
   */
  virtual const Epetra_LinearProblem* GetProblem() const = 0;

  //! Returns true if the solver can handle this matrix shape 
  /*! Returns true if the matrix shape is one that the underlying
    sparse direct solver can handle. Classes that work only on square
    matrices should return false for rectangular matrices.  Classes
    that work only on symmetric matrices whould return false for
    non-symmetric matrices.
  */
  virtual bool MatrixShapeOK() const = 0;

  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  virtual const Epetra_Comm & Comm() const = 0;

  //! Returns the number of symbolic factorizations performed by this object.
  virtual int NumSymbolicFact() const = 0;

  //! Returns the number of numeric factorizations performed by this object.
  virtual int NumNumericFact() const = 0;

  //! Returns the number of solves performed by this object.
  virtual int NumSolve() const = 0;

  //! Prints status information about the current solver.
  virtual void PrintStatus() const = 0;

  //! Prints timing information about the current solver. 
  virtual void PrintTiming() const = 0;

  //! Redefined from Teuchos::ParameterListAcceptor (Does Not Work)
  virtual void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList)
  {
    Teuchos::RCP<Teuchos::ParameterList> temp = paramList;
    //    paramList_ = paramlist ;
    //    this->SetParameters( *paramList_ );
  }

  //!  This is an empty stub 
  virtual Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList()
  {
    Teuchos::RCP<Teuchos::ParameterList> PL ;
    return PL ;
  }

  //!  This is an empty stub 
  virtual Teuchos::RCP<Teuchos::ParameterList> unsetParameterList()
  {
    Teuchos::RCP<Teuchos::ParameterList> PL ;
    //    this->SetParameters( *paramList_ );
    return PL ; 
  }

  //! Extracts timing information from the current solver and places it in the parameter list. (Does Not Work)
  virtual void GetTiming( Teuchos::ParameterList &TimingParameterList ) const 
  {
    Teuchos::ParameterList temp;
    TimingParameterList = temp.setName("NULL");
  }

  //@}



};

#endif /* _AMESOS_BASESOLVER_H_ */
