// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_VECTOR_H
#define ROL_VECTOR_H

#include "Teuchos_RefCountPtr.hpp"

/** \class ROL::Vector
    \brief Provides the vector space interface.

    The basic interface to be supplied by the user includes:\n
    \li vector addition,
    \li scalar multiplication,
    \li dot (scalar) product of vectors,
    \li vector norm,
    \li cloning of vectors.

    The dot product can represent an inner product (in Hilbert space) or
    a duality pairing (in general Banach space).

    There are additional virtual member functions that the user
    may want to reimplement for added efficiency. 
*/


namespace ROL {

template <class Real>
class Vector {
public:

  virtual ~Vector() {}

  /** \brief Compute \f$y \leftarrow x + y\f$ where \f$y = \mbox{*this}\f$.
  */
  virtual void plus( const Vector &x ) = 0;

  /** \brief Compute \f$y \leftarrow \alpha y\f$ where \f$y = \mbox{*this}\f$.
  */
  virtual void scale( const Real alpha ) = 0;

  /** \brief Returns \f$ \langle y,x \rangle \f$ where \f$y = \mbox{*this}\f$.
  */
  virtual Real dot( const Vector &x ) const = 0;

  /** \brief Returns \f$ \| y \| \f$ where \f$y = \mbox{*this}\f$.
  */
  virtual Real norm() const = 0;

  /** \brief Clone to make a new (uninitialized) vector.
  */
  virtual Teuchos::RCP<Vector> clone() const = 0;


  /** \brief Compute \f$y \leftarrow \alpha x + y\f$ where \f$y = \mbox{*this}\f$.
  */
  virtual void axpy( const Real alpha, const Vector &x ) {
    Teuchos::RCP<Vector> ax = x.clone();
    ax->set(x);
    ax->scale(alpha);
    this->plus(*ax);
  }

  /**  \brief Set to zero vector.
  */
  virtual void zero() {
    this->scale( (Real)0 );
  }

  /** \brief Return i-th basis vector: define if finite-difference gradients
             and Hessians are used.
  */
  virtual Teuchos::RCP<Vector> basis( const int i ) const {return Teuchos::null;}

  /** \brief Return dimension of the vector space.
  */
  virtual int dimension() const {return 0;}

  /**  \brief Set \f$y \leftarrow x\f$ where \f$y = \mbox{*this}\f$.
  */
  virtual void set( const Vector &x ) {
    this->zero();
    this->plus(x);
  }

}; // class Vector

} // namespace ROL

#endif
