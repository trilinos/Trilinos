
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

#include <cstddef>
#include <ostream>
#include <type_traits>

#include "XROL_MagnitudeTypeTraits.hpp"
#include "XROL_VectorOptionalMethod.hpp"
#include "XROL_TypeDeduction.hpp"


namespace XROL {





// Forward delare type traits for vector - must be specialized for derived types
template<class> struct VectorTypeTraits;

template<class Derived>
class Vector { 

private: 

  Derived* pdv_;
  const Derived* cpdv_;

public:

  using IndexType     = typename VectorTypeTraits<Derived>::index_type;
  using ElementType   = typename VectorTypeTraits<Derived>::element_type;
  using MagnitudeType = typename MagnitudeTypeTraits<ElementType>::type;

  Vector() : pdv_(static_cast<Derived*>(this)),  
             cpdv_(static_cast<const Derived*>(this)) {}

  virtual ~Vector() {}

  /***************************/  
  /* Constant object methods */
  /***************************/  

  auto basis( IndexType i )               const { return cpdv_->basis(i);    }
  auto clone()                            const { return cpdv_->clone();     }
  IndexType dimension()                   const { return cpdv_->dimension(); }
  ElementType dot( const Vector &x )      const { return cpdv_->dot(x);      }
  MagnitudeType norm()                    const { return cpdv_->norm();      }

  /** \brief Reducer should be a Functor class where the function
             call operator() is defined and should also define a method 
             initialValue() which returns the seed value 
             ( e.g. 0 for sum, inf for min, etc ) 
  */

  template<Reducer>
  auto reduce( const Reducer &r ) const { return cpdv_->reduce(r);   }


  auto dual() const { 
    return VectorOptionalMethod::dual(*cpdv_);    
  }

  void print( std::ostream &os ) const {        
    VectorOptionalMethod::print(*cpdv_, os);   
  }


  /*******************************/
  /* Non-constant object methods */
  /*******************************/
 
  /** \brief Callable should be either Functor class where the function
             call operator() is defined or a lambda function
  */
  template<class Callable, class ...Args>
  void applyFunction( const Callable& f, const Args&... args ) { 
    UniformTypeCheck(args);
    pdv_->apply(f,args...); 
  }
  

  void axpy(  const Real alpha, const Vector &x ) { pdv_->axpy(alpha,x); }
  void plus(  const Vector &x )                   { pdv_->plus(x);       }
  void scale( const Real alpha )                  { pdv_->scale(alpha);  }
  void set(   const Real &x )                     { pdv_->set(x);        }

}; // class Vector


} // namespace XROL
