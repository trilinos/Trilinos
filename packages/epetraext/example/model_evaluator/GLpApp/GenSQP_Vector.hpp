/*
//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// ***********************************************************************
//@HEADER
*/

#ifndef GENSQP_VECTOR_H
#define GENSQP_VECTOR_H

#include "Teuchos_RefCountPtr.hpp"

/** \class GenSQP::Vector
    \brief Provides the interface to generic abstract vector libraries.

    The interfaced functionality is very basic and includes routines for:\n
    \li linear combinations,
    \li inner products,
    \li scaling and copying operations,
    \li the cloning of a vector.
*/


namespace GenSQP {

class Vector {
public:

  virtual ~Vector() {}

  /** \brief Returns inner(*this,x).
  */
  virtual double innerProd( const Vector &x ) const = 0;

  /** \brief <tt>y = alpha*x + beta*y</tt> where <tt>y == *this</tt>.
  */
  virtual void linComb( const double &alpha, const Vector &x, const double &beta = 1.0 ) = 0;

  /** \brief <tt>y = alpha*y</tt> where <tt>y == *this</tt>.
  */
  virtual void Scale( const double &alpha ) = 0;

  /** \brief <tt>y = alpha</tt> where <tt>y == *this</tt>.
  */
  virtual void Set( const double &alpha ) = 0;

  /** \brief <tt>y = alpha*x</tt> where <tt>y == *this</tt>.
  */
  virtual void Set( const double &alpha, const Vector &x ) = 0;

  /** Clone to make a new (uninitialized) vector.
  */
  virtual Teuchos::RefCountPtr<Vector> createVector() const = 0;

}; // class Vector

} // namespace GenSQP

#endif
