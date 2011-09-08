//@HEADER
// ***********************************************************************
// 
//                Komplex: Complex Linear Solver Package
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

#ifndef KOMPLEX_OPERATOR_H
#define KOMPLEX_OPERATOR_H

#include "Komplex_KForms.hpp"
#include "Epetra_DataAccess.h"

class Epetra_MultiVector;
class Epetra_Map;
class Epetra_Comm;
class Epetra_Operator;

//! Komplex_Operator: A class for using complex-valued double-precision operators stored in equivalent real formation.
/*! This class implements the pure virtual Epetra_Operator class.

    The Komplex_Operator class enables the construction and use of equivalent real formulations of complex-valued,
    double-precision operators in a distributed memory environment.  
    
    <b>Constructing Komplex_Operators</b>
  
    Komplex_Operator constructors have two data access modes:
    <ol>
    <li> Copy mode - Allocates memory and makes a copy of the user-provided data. In this case, the
    user data is not needed after construction.
    <li> View mode - Creates a "view" of the user data. In this case, the
    user data is required to remain intact for the life of the operator.
    </ol>
  
    \warning View mode is \e extremely dangerous from a data hiding perspective.
    Therefore, we strongly encourage users to develop code using Copy mode first and 
    only use the View mode in a secondary optimization phase.
  
    There are two different Komplex_Operator constructors:
    <ul>
    <li> Copy from or make view of a Epetra_Operator object, with the real and imaginary parts interleaved.
    <li> Copy from or make view of two Epetra_Operator objects, representing the real and imaginary parts.
    <\ul>

*/    

class Komplex_Operator : public virtual Epetra_Operator {
      
 public:

  //@{ \name Constructors/Destructor.

  //! Komplex_Operator constructor from one object.
  /*! Creates a Komplex_Operator object from one Epetra_Operator with interleaved real and
      imaginary values.

    \param DataAccess (In) Copy or view.
    \param Operator (In) Epetra_Operator containing the real and imaginary values interleaved.
    \param KForm (In) The Komplex_KForms to use for this operator; by default, it is set to K1.

    \return Pointer to a Komplex_Operator.
  */
  Komplex_Operator(Epetra_DataAccess CV, Epetra_Operator* Operator, Komplex_KForms KForm = K1);

  //! Komplex_Operator constructor from two objects.
  /*! Creates a Komplex_Operator object from two Epetra_Operator objects, one representing
      the real values and the other representing the imaginary values.

    \param DataAccess (In) Copy or view.
    \param Real (In) Epetra_Operator containing the real values.
    \param Imag (In) Epetra_Operator containing the imaginary values.
    \param KForm (In) The Komplex_KForms to use for this operator; by default, it is set to K1.

    \return Pointer to a Komplex_Operator.
  */
  Komplex_Operator(Epetra_DataAccess CV, Epetra_Operator* Real, Epetra_Operator* Imag,
                   Komplex_KForms KForm = K1);

  //! Destructor
  ~Komplex_Operator();
  //@}
  
  //@{ \name Atribute set methods.

  //! If set true, the transpose of this operator will be applied.
  /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
      affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface 
	does not support transpose use, this method should return a value of -1.      
    \param UseTranspose (In) If true, multiply by the transpose of the operator, otherwise just use the operator.

    \return Integer error code, set to 0 if successful.  Set to -1 if this implementation does not support transpose.
  */
  int SetUseTranspose(bool UseTranspose);
  //@}
  
  //@{ \name Mathematical functions.

  //! Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
  /*! 
    \param X (In) A Epetra_MultiVector of dimension NumVectors to multiply with the operator.
    \param Y (Out) A Epetra_MultiVector of dimension NumVectors containing the result.

    \return Integer error code, set to 0 if successful.
  */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
  /*! 
    \param X (In) A Epetra_MultiVector of dimension NumVectors to multiply with the inverse of the operator.
    \param Y (Out) A Epetra_MultiVector of dimension NumVectors containing the result.

    \return Integer error code, set to 0 if successful.

    \warning In order to work with AztecOO, any implementation of this method must 
             support the case where X and Y are the same object.
  */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the infinity norm of the global matrix.
  /* Returns the quantity \f$ \| A \|_\infty\f$ such that
     \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].

    \warning This method must not be called unless HasNormInf() returns true.
  */ 
  double NormInf() const;
  //@}
  
  //@{ \name Attribute access functions

  //! Returns a character string describing the operator
  const char * Label() const;

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const;

  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  bool HasNormInf() const;

  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  const Epetra_Comm & Comm() const;

  //! Returns the Epetra_Map object associated with the domain of this operator.
  const Epetra_Map & OperatorDomainMap() const;

  //! Returns the Epetra_Map object associated with the range of this operator.
  const Epetra_Map & OperatorRangeMap() const;
  //@}

  protected:

  private:
  Komplex_Ordering* Ordering_;
  Epetra_Operator* Real_;
  Epetra_Operator* Imag_;
  bool IsOneObject_;
};

#endif /* KOMPLEX_OPERATOR_H */
