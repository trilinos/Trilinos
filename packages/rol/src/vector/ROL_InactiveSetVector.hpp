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

#ifndef ROL_INACTIVE_SET_VECTOR_HPP
#define ROL_INACTIVE_SET_VECTOR_HPP

#include "ROL_ScaledVector.hpp"
#include "ROL_BoundConstraint.hpp"

/** @ingroup la_group
    \class ROL::InactiveSet_PrimalVector
    \brief Defines the a Vector which has a diagonally scaled dot 
           product that neglects active set elements
           Used to simplify Semi-smooth Newton method implementation

    \class ROL::InactiveSet_DualVector
    \brief Defines the a Vector which has a diagonally scaled dot 
           product that neglects active set elements
           Used to simplify Semi-smooth Newton method implementation
*/

namespace ROL {

template<typename Real> class InactiveSet_PrimalVector;
template<typename Real> class InactiveSet_DualVector;


template<typename Real>
class InactiveSet_PrimalVector : public PrimalScaledVector<Real> {

  using V      = Vector<Real>;
  using Primal = InactiveSet_PrimalVector<Real>;
  using Dual   = InactiveSet_DualVector<Real>;
  using Bnd    = BoundConstraint<Real>;

private:

  mutable Ptr<V>  x_;            // Current optimization iterate
  Ptr<Bnd>        bnd_;          
  
public: 

  InactiveSet_PrimalVector( const Ptr<V>& vec,
                            const Ptr<V>& scaling_vec, 
                            const Ptr<V>& x,
                            const Ptr<Bnd>& bnd ) :
    PrimalScaledVector<Real>(vec,scaling_vec_), x_(x), bnd_(bnd) {}

  virtual ~InactiveSet_PrimalVector() {}


  Real dot( const V& x ) const override {

    auto& w = this->getWorkspace();  
    auto y  = w.copy(x);

    this->multiply_scaling( *y );

    // Set elements of y corresponsing the the active set of X to zero
    bnd_->pruneActive( *y, *x_ );

    return y->dot( *this->getVector() );    
  } 

  Ptr<V> clone() const override {
    return makePtr<Primal>( this->getVector()->clone(), 
                            this->getScalingVector(),
                            x_, bnd_ );
  }

  Ptr<V> basis( const int i ) const override { 
    return makePtr<Primal>( this->getVector()->basis(i),
                            this->getScalingVector(),
                            x_, bnd_ ); 
  }

  void const V& dual() const override {
    auto& w = this->getWorkspace();
    auto dual_vec = w.copy( this->getVector() );  
    this->multiply_scaling( dual_vec );
    return makePtr<Dual>( dual_vec, this->getScalingVector(), x_, bnd_ );
  } 

  void setIterateVector( const Ptr<V>& x ) const { x_->set(x); }
  const Ptr<V>& getIterateVector() { return x_; }
  const Ptr<const V>& getIterateVector() const { return x_; }


}; // class InactiveSet_PrimalVector 


template<typename Real>
class InactiveSet_DualVector : public DualScaledVector<Real> {

  using V      = Vector<Real>;
  using Primal = InactiveSet_PrimalVector<Real>;
  using Dual   = InactiveSet_DualVector<Real>;
  using Bnd    = BoundConstraint<Real>;

private:

  mutable Ptr<V>  x_;            // Current optimization iterate
  Ptr<Bnd>        bnd_;          
  
public: 

  InactiveSet_DualVector( const Ptr<V>& vec,
                          const Ptr<V>& scaling_vec, 
                          const Ptr<V>& x,
                          const Ptr<Bnd>& bnd ) :
    PrimalScaledVector<Real>(vec,scaling_vec_), x_(x), bnd_(bnd) {}

  virtual ~InactiveSet_PrimalVector() {}

  Real dot( const V& x ) const override {

    auto& w = this->getWorkspace();  
    auto y  = w.copy(x);
    this->divide_scaling( *y, this->getScalingVector() );

    // Set elements of y corresponsing the the active set of X to zero
    bnd_->pruneActive( *y, *x_ );

    return y->dot( *this->getVector() );    
  } 

  Ptr<V> clone() const override {
    return makePtr<Primal>( this->getVector()->clone(), 
                            this->getScalingVector(),
                            x_, bnd_ );
  }

  Ptr<V> basis( const int i ) const override { 
    return makePtr<Primal>( this->getVector()->basis(i),
                            this->getScalingVector(),
                            x_, bnd_ ); 
  }

  void const V& dual() const override {
    auto& w = this->getWorkspace();
    auto dual_vec = w.copy( this->getVector() );  
    this->multiply( dual_vec );
    return *( makePtr<Dual>( dual_vec, this->getScalingVector(), x_, bnd_ ) );
  } 

  void setIterateVector( const Ptr<V>& x ) const { x_->set(x); }
  const Ptr<V>& getIterateVector() { return x_; }
  const Ptr<const V>& getIterateVector() const { return x_; }

}; // class InactiveSet_PrimalVector 





} // namespace ROL



#endif // ROL_INACTIVE_SET_VECTOR_HPP
