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

/** @ingroup la_group
    \class ROL::Vector
    \brief Defines the linear algebra or vector space interface.

    The basic linear algebra interface, to be implemented by the user, includes:\n
    \li #plus -- vector addition;
    \li #scale -- scalar multiplication;
    \li #dot -- dot (scalar) product of vectors;
    \li #norm -- vector norm;
    \li #clone -- cloning of existing vectors.

    The dot product can represent an inner product (in Hilbert space) or
    a duality pairing (in general Banach space).

    There are additional virtual member functions that can be
    overloaded for computational efficiency. 
*/

namespace ROL {

template <class Real>
class Vector {
public:

  virtual ~Vector() {}


  /** \brief Compute \f$y \leftarrow y + x\f$, where \f$y = \mathtt{*this}\f$.

             @param[in]      x  is the vector to be added to \f$\mathtt{*this}\f$.

             On return \f$\mathtt{*this} = \mathtt{*this} + x\f$.

             ---
  */
  virtual void plus( const Vector &x ) = 0;


  /** \brief Compute \f$y \leftarrow \alpha y\f$ where \f$y = \mathtt{*this}\f$.

             @param[in]      alpha is the scaling of \f$\mathtt{*this}\f$.

             On return \f$\mathtt{*this} = \alpha (\mathtt{*this}) \f$.

             ---
  */
  virtual void scale( const Real alpha ) = 0;


  /** \brief Compute \f$ \langle y,x \rangle \f$ where \f$y = \mathtt{*this}\f$.

             @param[in]      x  is the vector that forms the dot product with \f$\mathtt{*this}\f$.
             @return         The number equal to \f$\langle \mathtt{*this}, x \rangle\f$.

             ---
  */
  virtual Real dot( const Vector &x ) const = 0;


  /** \brief Returns \f$ \| y \| \f$ where \f$y = \mathtt{*this}\f$.

             @return         A nonnegative number equal to the norm of \f$\mathtt{*this}\f$.

             ---
  */
  virtual Real norm() const = 0;


  /** \brief Clone to make a new (uninitialized) vector.

             @return         A reference-counted pointer to the cloned vector.

             Provides the means of allocating temporary memory in ROL.

             ---             
  */
  virtual Teuchos::RCP<Vector> clone() const = 0;


  /** \brief Compute \f$y \leftarrow \alpha x + y\f$ where \f$y = \mathtt{*this}\f$.

             @param[in]      alpha is the scaling of @b x.
             @param[in]      x     is a vector.

             On return \f$\mathtt{*this} = \mathtt{*this} + \alpha x \f$.
             Uses #clone, #set, #scale and #plus for the computation.
             Please overload if a more efficient implementation is needed.

             ---
  */
  virtual void axpy( const Real alpha, const Vector &x ) {
    Teuchos::RCP<Vector> ax = x.clone();
    ax->set(x);
    ax->scale(alpha);
    this->plus(*ax);
  }

  /** \brief Set to zero vector.

             Uses #scale by zero for the computation.
             Please overload if a more efficient implementation is needed.

             ---
  */
  virtual void zero() {
    this->scale( (Real)0 );
  }


  /** \brief Return i-th basis vector.

             @param[in] i is the index of the basis function.
             @return A reference-counted pointer to the basis vector with index @b i.

             Overloading the basis is only required if the default gradient implementation
             is used, which computes a finite-difference approximation.

             ---
  */
  virtual Teuchos::RCP<Vector> basis( const int i ) const {return Teuchos::null;}


  /** \brief Return dimension of the vector space.

             @return The dimension of the vector space, i.e., the total number of basis vectors.

             Overload if the basis is overloaded.

             ---
  */
  virtual int dimension() const {return 0;}


  /** \brief Set \f$y \leftarrow x\f$ where \f$y = \mathtt{*this}\f$.

             @param[in]      x     is a vector.

             On return \f$\mathtt{*this} = x\f$.
             Uses #zero and #plus methods for the computation.
             Please overload if a more efficient implementation is needed.

             ---
  */
  virtual void set( const Vector &x ) {
    this->zero();
    this->plus(x);
  }


  /** \brief Return dual representation of \f$\mathtt{*this}\f$, for example,
             the result of applying a Riesz map, or change of basis, or
             change of memory layout.

             @return         A const reference to dual representation.

             By default, returns the current object.
             Please overload if you need a dual representation.

             ---
  */
  virtual const Vector & dual() const {
    return *this;
  }

}; // class Vector

} // namespace ROL

#endif
