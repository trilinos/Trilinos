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
#ifndef ROL_REDUCTIONOP_INTERFACE_H
#define ROL_REDUCTIONOP_INTERFACE_H

namespace ROL {
namespace Elementwise {

enum EReductionType {
  REDUCE_AND, ///< Logical AND
  REDUCE_EN2, ///< EuclideanNormSquared
  REDUCE_MAX, ///< Max
  REDUCE_MIN, ///< Min
  REDUCE_SUM ///< Sum
};

// Forward Declarations of ReductionOp derived types
template<typename> class EuclideanNormSquared;
template<typename> class ReductionAnd;
template<typename> class ReductionMax;
template<typename> class ReductionMin;
template<typename> class ReductionSum;

// Generic reduction operation
template<class Real>
class ReductionOp {
public:
  virtual ~ReductionOp() {}
  virtual void reduce( const Real& input, Real& output ) const = 0;
  virtual void reduce( const volatile Real& input, volatile Real& output ) const = 0; 
  virtual Real initialValue() const = 0;
  virtual EReductionType reductionType() const = 0;

  struct Visitor {
    virtual void visit( const EuclideanNormSquared<Real>& ) {}
    virtual void visit( const ReductionAnd<Real>&         ) {}
    virtual void visit( const ReductionMax<Real>& ) {}
    virtual void visit( const ReductionMin<Real>& ) {}
    virtual void visit( const ReductionSum<Real>& ) {}
  };

  virtual void accept( ReductionOp::Visitor& ) const {}
};

template<class Real>
class EuclideanNormSquared : public ReductionOp<Real> {
public:

  virtual ~EuclideanNormSquared() = default;

  void reduce( const Real& input, Real& output ) const override {
    output += input * input;
  }

  void reduce( const volatile Real& input, volatile Real& output ) const override {
    output += input * input;
  }

  Real initialValue() const override {
    return 0;
  }

  EReductionType reductionType() const override {
    return Elementwise::REDUCE_EN2;
  }

  void accept( typename ReductionOp<Real>::Visitor& visitor ) const override {
    visitor.visit(*this);
  }

};

template<typename Real>
class ReductionAnd : public ReductionOp<Real> {
public:
  virtual ~ReductionAnd() = default;

  void reduce( const Real& input, Real& output ) const override {
    output = (input*output)==0 ? 0.0 : 1.0;
  }

  void reduce( const volatile Real& input, volatile Real& output ) const override {
    output = (input*output)==0 ? 0.0 : 1.0;
  }

  Real initialValue() const override {
    return 1.0;
  }

  EReductionType reductionType() const override {
    return Elementwise::REDUCE_AND;
  }

  void accept( typename ReductionOp<Real>::Visitor& visitor ) const override {
    visitor.visit(*this);
  }
};


template<class Real>
class ReductionMax : public ReductionOp<Real> {
public:
  virtual ~ReductionMax() = default;

  ReductionMax() {
      ROL_TEST_FOR_EXCEPTION(!std::numeric_limits<Real>::is_specialized,
      std::logic_error,"ReductionMax() requires std::numeric_limits "
        "be specialized on supplied template parameter.\n" );
  }

  void reduce( const Real& input, Real& output ) const override {
    output = (input>output) ? input : output;
  }

  void reduce( const volatile Real& input, volatile Real& output ) const override {
    output = (input>output) ? input : output;
  }

  Real initialValue() const override {
    return std::numeric_limits<Real>::min();
  }

  EReductionType reductionType() const override {
    return Elementwise::REDUCE_MAX;
  }

  void accept( typename ReductionOp<Real>::Visitor& visitor ) const override {
    visitor.visit(*this);
  }
};

template<class Real>
class ReductionMin : public ReductionOp<Real> {
public:
  virtual ~ReductionMin() = default;

  ReductionMin() {
      ROL_TEST_FOR_EXCEPTION(!std::numeric_limits<Real>::is_specialized,
      std::logic_error,"ReductionMin() requires std::numeric_limits "
        "be specialized on supplied template parameter.\n" );
  }

  void reduce( const Real& input, Real& output ) const override {
    output = (input<output) ? input : output;
  }

  void reduce( const volatile Real& input, Real volatile& output ) const override {
    output = (input<output) ? input : output;
  }

  Real initialValue() const override {
    return std::numeric_limits<Real>::max();
  }

  EReductionType reductionType() const override {
    return Elementwise::REDUCE_MIN;
  }

  void accept( typename ReductionOp<Real>::Visitor& visitor ) const override {
    visitor.visit(*this);
  }
};



template<class Real>
class ReductionSum : public ReductionOp<Real> {
public:
  virtual ~ReductionSum() = default;

  void reduce( const Real& input, Real& output ) const override {
    output = output + input;
  }

  void reduce( const volatile Real& input, volatile Real& output ) const override {
    output = output + input;
  }

  Real initialValue() const override {
    return 0;
  }

  EReductionType reductionType() const override {
    return Elementwise::REDUCE_SUM;
  }

  void accept( typename ReductionOp<Real>::Visitor& visitor ) const override {
    visitor.visit(*this);
  }

};



} // namespace Elementwise 
} // namespace ROL

#endif  // ROL_REDUCTIONOP_INTERFACE_H
