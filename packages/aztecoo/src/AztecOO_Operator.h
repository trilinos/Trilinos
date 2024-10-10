
/*
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

#ifndef _AZTECOO_OPERATOR_H_
#define _AZTECOO_OPERATOR_H_

#if defined(AztecOO_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The AztecOO package is deprecated"
#endif
#endif

class Epetra_MultiVector;
class Epetra_BlockMap;
class Epetra_Comm;
#include "Epetra_Operator.h"
#include "AztecOO.h"

//! AztecOO_Operator: An implementation of the Epetra_Operator class.
/*! The AztecOO_Operator class implements Epetra_Operator using a pre-constructed AztecOO solver object.
    Once constructed, an AztecOO_Operator can be used as a preconditioner within another AztecOO solver
    object.
*/    

class AztecOO_Operator: public virtual Epetra_Operator {
      
 public:

  //@{ \name Constructor.
    //! Uses an AztecOO instance to implement the Epetra_Operator interface.
  /*! Facilitates the use of an AztecOO solver instance as an operator.  This is particularly designed
      for using AztecOO as a preconditioner within another AztecOO instance.
    \param In - A fully-constructed AztecOO object.
    \param In - Number of iterations that should be performed.  \em Exactly this many iterations will be done if Tol = 0.0.
    \param In - Tolerance used for each application of AztecOO solver.
  */
  AztecOO_Operator(AztecOO * solver, int NumIters, double Tol = 0.0);
    //! Destructor
  ~AztecOO_Operator();
  //@}
  
  //@{ \name Attribute set methods.

    //! If set true, transpose of this operator will be applied.
    /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
        affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface 
	does not support transpose use, this method should return a value of -1.
      
    \param In
	   use_transpose - If true, multiply by the transpose of operator, otherwise just use operator.

    \warning - This returns -1 if use_transpose is true, because tranpse is not supported.
  */
  int SetUseTranspose(bool use_transpose)
    {
      if (use_transpose == true) return(-1);
      return(0);
    }
  //@}
  
  //@{ \name Mathematical functions.

    //! Returns the result of a AztecOO_Operator applied to a Epetra_MultiVector X in Y.
    /*! 
    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \warning - This method has no effect and returns -1 as error code.
  */
  int Apply(const Epetra_MultiVector& X,
            Epetra_MultiVector& Y) const
    {(void)X; (void)Y; return(-1);}

  //! Returns the result of a AztecOO_Operator inverse applied to an Epetra_MultiVector X in Y.
  /*! 
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectors containing result.
    
    \return Integer error code, set to 0 if successful.
  */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  
  //! Returns the infinity norm of the global matrix.
  /* Returns the quantity \f$ \| A \|_\infty\f$ such that
     \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].
     
     \warning This method must not be called unless HasNormInf() returns true.
  */ 
  double NormInf() const {return(0.0);};
  //@}
  
  //@{ \name Attribute access functions

  //! Returns a character string describing the operator
  const char * Label() const {return(Label_);};

  //! Returns a pointer to the AztecOO solver object that was used to create this AztecOO_Operator object.
  AztecOO * Solver() const {return(solver_);};

  //! Returns the number of iterations that will be performed with the AztecOO solver.
  int NumIters() const {return(NumIters_);};

  //! Returns the tolerance this will be used by the AztecOO solver.
  double Tol() const {return(Tol_);};
  
  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(false);};
  
  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  bool HasNormInf() const{return(false);};
  
  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  const Epetra_Comm & Comm() const{return(solver_->GetUserOperator()->Comm());};
  
  //! Returns the Epetra_BlockMap object associated with the domain of this matrix operator.
  const Epetra_Map & OperatorDomainMap() const
  {
    if (UseTranspose()) return(solver_->GetUserOperator()->OperatorRangeMap());
    else return(solver_->GetUserOperator()->OperatorDomainMap());
  }
  
  //! Returns the Epetra_BlockMap object associated with the range of this matrix operator.
  const Epetra_Map & OperatorRangeMap() const
  {
    if (UseTranspose()) return(solver_->GetUserOperator()->OperatorDomainMap());
    else return(solver_->GetUserOperator()->OperatorRangeMap());
  }
  //@}
  
 protected:

  AztecOO * solver_;
  int NumIters_;
  double Tol_;
  const char * Label_;
};

#endif /* _AZTECOO_OPERATOR_H_ */
