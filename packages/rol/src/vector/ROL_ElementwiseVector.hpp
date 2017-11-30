
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

#ifndef ROL_ELEMENTWISE_VECTOR_H
#define ROL_ELEMENTWISE_VECTOR_H

#include "ROL_Vector.hpp"

/** @ingroup la_group
    \class ROL::ElementwiseVector
    \brief Intermediate abstract class which does not require users 
           implements plus, set, scale, axpy, norm, dot, or zero if 
           they implement the three elementwise functions: applyUnary,
           applyBinary, and reduce

           dot and norm are unweighted dot products and Euclidean norm 
           by default
*/

namespace ROL {

template< class Real>
class ElementwiseVector : public Vector<Real> {

public: 

  virtual ~ElementwiseVector() {}

  void plus( const Vector<Real> &x ) {
    this->applyBinary(Elementwise::Plus<Real>(),x);
  }

  void scale( const Real alpha ) {
    this->applyUnary(Elementwise::Scale<Real>(alpha));
  }

  virtual Real dot( const Vector<Real> &x ) const {
    ROL::Ptr<Vector<Real> > y = this->clone();
    y->set(*this);
    y->applyBinary(Elementwise::Multiply<Real>(),x);
    return y->reduce(Elementwise::ReductionSum<Real>());    
  }

  virtual Real norm() const {
      return std::sqrt(this->dot(*this));
  }

  void axpy( const Real alpha, const Vector<Real> &x ) {
    this->applyBinary(Elementwise::Axpy<Real>(alpha),x);
  }

  void zero() {
    this->applyUnary(Elementwise::Fill<Real>(Real(0)));
  }

  void set( const Vector<Real> &x ) {
    this->applyBinary(Elementwise::Set<Real>(),x);
  }

  // MUST overload these three functions
  virtual void applyUnary( const Elementwise::UnaryFunction<Real> &uf ) = 0;

  virtual void applyBinary( const Elementwise::BinaryFunction<Real> &bf,
                            const Vector<Real> &x ) = 0;

  virtual Real reduce( const Elementwise::ReductionOp<Real> &r ) const = 0;

}; // class ElementwiseVector


} // namespace ROL




#endif // ROL_ELEMENTWISE_VECTOR_H

