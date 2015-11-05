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

#ifndef ROL_UNARYFUNCTIONS_H
#define ROL_UNARYFUNCTIONS_H

#include <cstdlib>
#include <ctime>

#include "ROL_Elementwise_Function.hpp"

namespace ROL {
namespace Elementwise {

// Used to set every element in a vector to a specific value
template<class Real>
class Fill : public UnaryFunction<Real> {
public:
  Fill( const Real &value ) : value_(value) {}
  Real apply( const Real &x ) const {
    return value_;
  }  
private:  
  Real value_;
}; // class Fill


// Get the elementwise reciprocal of a vector
template<class Real> 
class Reciprocal : public UnaryFunction<Real> {
public:
  Real apply( const Real &x ) const {
    return static_cast<Real>(1)/x;
  }  
}; // class Reciprocal


// Generate a uniformly distributed random number
// between lower and upper
template<class Real> 
class UniformlyRandom : public UnaryFunction<Real> {
private:
  const Real lower_;
  const Real upper_;

public:
  UniformlyRandom( const Real &lower = 0.0, const Real &upper = 1.0) : 
    lower_(lower), upper_(upper) {
  }

  Real apply( const Real &x ) const {
    return (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX)) * (upper_-lower_) + lower_;
  }
}; // class UniformlyRandom


} // namespace Elementwise
} // namespace ROL

#endif // ROL_UNARYFUNCTIONS_H
