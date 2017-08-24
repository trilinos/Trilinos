
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

#pragma once

#include "XROL_Vector.hpp"
#include "XROL_ElementwiseFunction.hpp"

namespace XROL {

// Forward declare class for traits
template<class> class ElementwiseVector;

// Traits specialization
template<class Derived> 
struct VectorTypeTraits<ElementwiseVector<Derived>> {
  using index_type   = typename VectorTypeTraits<Derived>::index_type;
  using element_type = typename VectorTypeTraits<Derived>::element_type;
}; 


template<class Derived>
class ElementwiseVector : public Vector<ElementwiseVector<Derived>> {
private:

  Derived* pdv_;
  const Derived* cpdv_;

public:

  using IndexType     = typename VectorTypeTraits<Derived>::index_type;
  using ElementType   = typename VectorTypeTraits<Derived>::element_type;
  using MagnitudeType = typename MagnitudeTypeTraits<ElementType>::type;

  ElementwiseVector() : pdv_(static_cast<Derived*>(this)),  
                        cpdv_(static_cast<const Derived*>(this)) {}

  virtual ~ElementwiseVector() {}  

  /********************************/
  /* Implemented in Derived class */
  /********************************/

  template<class Callable, class ...Args>
  void applyFunction( const Callable& f, const Args&... args ) { 
    pdv_->applyFunction(f,args...); 
  }

  auto basis( IndexType i )               const { return cpdv_->basis(i);    }
  auto clone()                            const { return cpdv_->clone();     }
  IndexType dimension()                   const { return cpdv_->dimension(); }
  auto dual()                             const { return getDual(*cpdv_);    }
  ElementType reduce( const ReduceOp &r ) const { return cpdv_->reduce(r);   }
  void print( std::ostream &os )          const {        cpdv_->print(os);   }

  /****************************************/
  /* Implemented in using derived methods */
  /****************************************/
  void axpy( const Real alpha, const Vector &x ) {
    Elementwise::Axpy<ElementType> f(alpha);
    pdv_->applyFunction(axpy,x);
  }

  ElementType dot( const Vector &x ) const {
        
  } 

  MagnitudeType norm() const {
  
  }

  void plus( const Vector &x ) { 
    pdv_->applyFunction( Elementwise::sum, x ); 
  } 

  void scale( const Real alpha ) {
    Elementwise::Scale<ElementType> f(alpha);
    pdv_->applyFunction( f );
  }

  void set( const Real &x ) {
    pdv_->applyFunction( Elementwise::set, x );      
  } 

  void zero() {
    pdv_->applyFunction( Elementwise::zero );
  }
  

};

} // namespace XROL 
