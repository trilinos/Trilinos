/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
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

#ifndef IFPACK_REORDERING_H
#define IFPACK_REORDERING_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include "Ifpack_ConfigDefs.h"

namespace Teuchos {
  class ParameterList;
}
class Epetra_MultiVector;
class Ifpack_Graph;
class Epetra_RowMatrix;

//! Ifpack_Reordering: basic class for reordering for a Ifpack_Graph object.
/*!
Class Ifpack_Reordering is a pure virtual class that defines the
structure of all Ifpack reordering.

The Ifpack_Graph object is used \e only by method Compute().

A typical code reads as follows (using for instance RCM reordering):
\code
#include "Ifpack_Reordering.h"
#include "Ifpack_RCMReordering.h"
#include "Ifpack_Graph.h"
// A is an Epetra_RowMatrix pointer.
// Need to create a graph (which is a simple wrapper)
// This required include file Ifpack_Graph_Epetra_RowMatrix.h
Ifpack_Graph_Epetra_RowMatrix Graph(A);

// Construct the object
Ifpack_RCMReordering Reorder(Graph);
// Compute the reordering.
IFPACK_CHK_ERR(Reorder.Compute());
// Prints out some information
cout << Reorder;
\endcode

<P>An Ifpack_Reordering object is a tool used by class Ifpack_Preconditioner,
to reorder the localized matrix (with or without overlap). As its
basic usage is for localized matrices, this class takes care of
reordering the \e local rows only. It is also supposed that the
input graph contains no singletons. This is not a limitation, as
class Ifpack_AdditiveSchwarz will filter the graph using
Ifpack_SingletonFilter before using reordering.

<P>If IFPACK is configure with Teuchos support, method SetParameters()
should be adopted. Otherwise, users can set parameters (one at-a-time),
using methods SetParameter(), for integers and doubles.

<P>Ifpack_Preconditioner objects overload the << operator. Derived
classes should specify a Print() method, that will be used in
operator <<.

\author Marzio Sala, SNL 9214.

\date Last modified: Oct-04.
*/

class Ifpack_Reordering {

public:

  //! Destructor.
  virtual ~Ifpack_Reordering() {};

  //! Sets integer parameters `Name'.
  virtual int SetParameter(const std::string Name, const int Value) = 0;

  //! Sets double parameters `Name'.
  virtual int SetParameter(const std::string Name, const double Value) = 0;

  //! Sets all parameters.
  virtual int SetParameters(Teuchos::ParameterList& List) = 0;

  //! Computes all it is necessary to initialize the reordering object.
  virtual int Compute(const Ifpack_Graph& Graph) = 0;

  //! Computes all it is necessary to initialize the reordering object.
  virtual int Compute(const Epetra_RowMatrix& Matrix) = 0;

  //! Returns \c true is the reordering object has been successfully initialized, false otherwise.
  virtual bool IsComputed() const = 0;

  //! Returns the reordered index of row \c i.
  virtual int Reorder(const int i) const = 0;

  //! Returns the inverse reordered index of row \c i.
  virtual int InvReorder(const int i) const = 0;

  //! Applies reordering to multivector Xorig, whose local length equals the number of local rows, stores reordered vector in X.
  virtual int P(const Epetra_MultiVector& Xorig,
                Epetra_MultiVector& X) const = 0;

  //! Applies inverse reordering to multivector Xorig, whose local length equals the number of local rows, stores inverse reordered vector in X.
  virtual int Pinv(const Epetra_MultiVector& Xorig,
                   Epetra_MultiVector& X) const = 0;

  //! Prints basic information on iostream. This function is used by operator<<.
  virtual std::ostream& Print(std::ostream& os) const = 0;

};

inline std::ostream& operator<<(std::ostream& os, const Ifpack_Reordering& obj)
{
  return(obj.Print(os));
}

#endif
