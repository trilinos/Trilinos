This file is out of date.  Is has not been refactored to use Amesos_Status. 

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

#ifndef _AMESOS_BTF_H_
#define _AMESOS_BTF_H_

#if defined(Amesos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Amesos package is deprecated"
#endif
#endif

#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_Comm.h"
#endif
#include "Epetra_CrsGraph.h"


//! Amesos_Btf:  Factors and solves a matrix after converting it to block triangular form.
/*!
    Amesos_Btf:
<ul>
<li>Compute an ordering which reduces the matrix to block upper 
triangular form. 
<li>Determine a task partitioning (static load balancing)
<li>Redistribute the data based on the owner computes rule
<ul>
<li>The process(es) to which a given diagonal block is assigned owns all rows in that diagonal block.   
</ul>
<li>Uses an Amesos solver on each of the diagonal blocks
<li>Uses Trilinos operations to handle the off-diagonal blocks.
</ul>
*/

class Amesos_Btf: public Amesos_BaseSolver { 

public: 

  //@{ \name Constructor methods
  //! Amesos_Btf Constructor.
  /*! Creates an Amesos_Btf instance, using an Epetra_LinearProblem,
      passing in an already-defined Epetra_LinearProblem object. 

      Note: The operator in LinearProblem must be an
      Epetra_RowMatrix.

  */
  Amesos_Btf(const Epetra_LinearProblem& LinearProblem );

  //! Amesos_Btf Destructor.
  /*! Completely deletes an Amesos_Btf object.  
  */
  ~Amesos_Btf(void);
  //@}

  //@{ \name Mathematical functions.

    //! Performs SymbolicFactorization on the matrix A.
    /*! 
<ul>
<li>Compute an ordering which reduces the matrix to block upper 
triangular form. 
<li>Determine a task partitioning (static load balancing)
<li>Redistribute the data based on the owner computes rule
<li>Instantiates an Amesos solver for each of the diagonal blocks
<li>Calls SymbolicFactorization() on each of the diagonal blocks
</ul>

    \return Integer error code, set to 0 if successful.
  */
    int SymbolicFactorization() ;

    //! Performs NumericFactorization on the matrix A.
    /*!  
<ul>
<li>Calls NumericFactorization() on each of the diagonal blocks
</ul>

     \return Integer error code, set to 0 if successful.
  */
    int NumericFactorization() ;

    //! Solves A X = B (or A<SUP>T</SUP> X = B) 
    /*! 

<pre>
Foreach block i:
   For each block j
     Compute x_i -= A_{i,j} x_j 
   Call Solve(x_i,b_i) 
   Broadcast x_i
</pre>

     \return Integer error code, set to 0 if successful.
  */
    int Solve();

  //@}
  
  //@{ \name Additional methods required to support the Epetra_Operator interface.

  //! Get a pointer to the Problem.
  const Epetra_LinearProblem *GetProblem() const { return(Problem_); };

  //! Returns true if BTF can handle this matrix shape 
  /*! Returns true if the matrix shape is one that BTF can
    handle. BTF only works with square matrices.  
  */
  bool MatrixShapeOK() const ;

  //! SetUseTranpose(true) causes Solve() To compute A^T X = B
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

    \return Integer error code, set to 0 if successful. 
   */
  int SetParameters( Teuchos::ParameterList &ParameterList )  ;

  //! Returns the number of symbolic factorizations performed by this object.
  int NumSymbolicFact() const { return( NumSymbolicFact_ ); }

  //! Returns the number of numeric factorizations performed by this object.
  int NumNumericFact() const { return( NumNumericFact_ ); }

  //! Returns the number of solves performed by this object.
  int NumSolve() const { return( NumSolve_ ); }

  //! Print timing information
  void PrintTiming();
  
  //! Print information about the factorization and solution phases.
  void PrintStatus();
  
  //@}

private:  
  bool UseTranspose_;
  const Epetra_LinearProblem * Problem_;

  bool PrintTiming_;
  bool PrintStatus_;
  bool ComputeVectorNorms_;
  bool ComputeTrueResidual_;
  
  int verbose_;
  int debug_;

  // some timing internal, copied from MUMPS
  double ConTime_;                        // time to convert to BTF format
  double SymTime_;                        // time for symbolic factorization
  double NumTime_;                        // time for numeric factorization
  double SolTime_;                        // time for solution
  double VecTime_;                        // time to redistribute vectors
  double MatTime_;                        // time to redistribute matrix
  
  int NumSymbolicFact_;
  int NumNumericFact_;
  int NumSolve_;  

  Epetra_Time * Time_;
  
};  // End of  class Amesos_Btf  
#endif /* _AMESOS_BTF_H_ */
