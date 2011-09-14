/*
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
//@HEADER
*/

#ifndef EPETRA_OPERATOR_H
#define EPETRA_OPERATOR_H

#include "Epetra_ConfigDefs.h"
class Epetra_MultiVector;
class Epetra_Map;
class Epetra_Comm;

//! Epetra_Operator: A pure virtual class for using real-valued double-precision operators.
/*! The Epetra_Operator class is a pure virtual class (specifies interface only) that 
    enable the use of real-valued double-precision operators. It is currently implemented by both the
    Epetra_CrsMatrix and Epetra_VbrMatrix classes and the Ifpack_CrsRiluk preconditioner class.

   
*/    

class EPETRA_LIB_DLL_EXPORT Epetra_Operator {
      
 public:

   //! @name Destructor
  //@{ 
    //! Destructor
    virtual ~Epetra_Operator() {};
  //@}
  
  //! @name Atribute set methods
  //@{ 

    //! If set true, transpose of this operator will be applied.
    /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
        affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface 
	does not support transpose use, this method should return a value of -1.
      
    \param In
	   UseTranspose -If true, multiply by the transpose of operator, otherwise just use operator.

    \return Integer error code, set to 0 if successful.  Set to -1 if this implementation does not support transpose.
  */
    virtual int SetUseTranspose(bool UseTranspose) = 0;
  //@}
  
  //! @name Mathematical functions
  //@{ 

    //! Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
    /*! 
    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
    virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;

    //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
    /*! 
    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning In order to work with AztecOO, any implementation of this method must 
              support the case where X and Y are the same object.
  */
    virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;

    //! Returns the infinity norm of the global matrix.
    /* Returns the quantity \f$ \| A \|_\infty\f$ such that
       \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].

       \warning This method must not be called unless HasNormInf() returns true.
    */ 
    virtual double NormInf() const = 0;
  //@}
  
  //! @name Atribute access functions
  //@{ 

    //! Returns a character string describing the operator
    virtual const char * Label() const = 0;

    //! Returns the current UseTranspose setting.
    virtual bool UseTranspose() const = 0;

    //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
    virtual bool HasNormInf() const = 0;

    //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
    virtual const Epetra_Comm & Comm() const = 0;

    //! Returns the Epetra_Map object associated with the domain of this operator.
    virtual const Epetra_Map & OperatorDomainMap() const = 0;

    //! Returns the Epetra_Map object associated with the range of this operator.
    virtual const Epetra_Map & OperatorRangeMap() const = 0;
  //@}

};

#endif /* EPETRA_OPERATOR_H */
