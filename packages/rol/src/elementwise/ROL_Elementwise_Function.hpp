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

#ifndef ROL_ELEMENTWISE_FUNCTION_H
#define ROL_ELEMENTWISE_FUNCTION_H

#include "Teuchos_CommHelpers.hpp"

namespace ROL {
namespace Elementwise {


template<class Real>
class ReductionOp {
public:
  virtual ~ReductionOp() {}
  virtual void reduce( const Real &input, Real &output ) = 0;
  virtual Real initialValue() = 0;
  virtual Teuchos::EReductionType reductionType() = 0;
};


// f(x)
template<class Real>
class UnaryFunction {
public:
  virtual ~UnaryFunction() {}
  virtual Real apply( const Real &x ) = 0;
};

template<class Real, class Func>
class UnaryFunctionWrapper : public UnaryFunction<Real> {
public:
  UnaryFunctionWrapper( Func f ) : f_(f) {}
  
  Real apply( const Real &x ) {
    return f_(x);
  }
private:
  Func f_;
};


// f(x,y)
template<class Real>
class BinaryFunction {
public:
  virtual ~BinaryFunction() {}
  virtual Real apply( const Real &x, const Real &y ) = 0;
};

template<class Real, class Func>
class BinaryFunctionWrapper : public BinaryFunction<Real> {
public:
  BinaryFunctionWrapper( Func f ) : f_(f) {}
  
  Real apply( const Real &x, const Real &y ) {
    return f_(x,y); 
  }
private:
  Func f_;   
};


} // namespace Elementwise 
} // namespace ROL

#include <ROL_Elementwise_Implementation.hpp>
#include <ROL_Elementwise_Reduce.hpp>

#endif
