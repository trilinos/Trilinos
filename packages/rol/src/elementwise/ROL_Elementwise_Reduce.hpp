// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_VECTOR_H
#include "ROL_Vector.hpp"
#else

#ifndef ROL_ELEMENTWISE_REDUCE_H
#define ROL_ELEMENTWISE_REDUCE_H

#include <ROL_Elementwise_Function.hpp>
#include <ROL_Types.hpp>

namespace ROL {
namespace Elementwise {

enum EReductionType {
  REDUCE_SUM, ///< Sum
  REDUCE_MIN, ///< Min
  REDUCE_MAX, ///< Max
  REDUCE_AND, ///< Logical AND
  REDUCE_BOR  ///< Bitwise OR
};

// Generic reduction operation
template<class Real>
class ReductionOp {
public:
  virtual ~ReductionOp() {}
  virtual void reduce( const Real &input, Real &output ) const = 0;
  virtual void reduce( const volatile Real &input, volatile Real &output ) const = 0; 
  virtual Real initialValue() const = 0;
  virtual EReductionType reductionType() const = 0;
};

template<class Real> 
class ReductionSum : public ReductionOp<Real> {
public: 
  void reduce( const Real &input, Real &output ) const {
    output = output + input;
  }

  void reduce( const volatile Real &input, volatile Real &output ) const {
    output = output + input;
  }

  Real initialValue() const {
    return 0;
  }

  EReductionType reductionType() const {
    return Elementwise::REDUCE_SUM;
  }

};

template<class Real> 
class ReductionAnd : public ReductionOp<Real> {
public:
  void reduce( const Real &input, Real &output ) const {
    output = (input*output)==0 ? 0.0 : 1.0;
  }

  void reduce( const volatile Real &input, volatile Real &output ) const {
    output = (input*output)==0 ? 0.0 : 1.0;
  }

  Real initialValue() const {
    return 1.0;
  }

  EReductionType reductionType() const {
    return Elementwise::REDUCE_AND;
  }
};

template<class Real>
class ReductionMin : public ReductionOp<Real> {
public:
  ReductionMin() {
      ROL_TEST_FOR_EXCEPTION(!std::numeric_limits<Real>::is_specialized, 
      std::logic_error,"ReductionMin() requires std::numeric_limits "  
        "be specialized on supplied template parameter.\n" ); 
  }

  void reduce( const Real &input, Real &output ) const {
    output = (input<output) ? input : output;
  }

  void reduce( const volatile Real &input, Real volatile &output ) const {
    output = (input<output) ? input : output;
  }
 
  Real initialValue() const {
    return std::numeric_limits<Real>::max();
  }

  EReductionType reductionType() const {
    return Elementwise::REDUCE_MIN;
  }

};

template<class Real>
class ReductionMax : public ReductionOp<Real> {
public:
  ReductionMax() {
      ROL_TEST_FOR_EXCEPTION(!std::numeric_limits<Real>::is_specialized, 
      std::logic_error,"ReductionMax() requires std::numeric_limits "  
        "be specialized on supplied template parameter.\n" ); 
  }

  void reduce( const Real &input, Real &output ) const {
    output = (input>output) ? input : output;
  }

  void reduce( const volatile Real &input, volatile Real &output ) const {
    output = (input>output) ? input : output;
  }
 
  Real initialValue() const {
    return std::numeric_limits<Real>::min();
  }

  EReductionType reductionType() const {
    return Elementwise::REDUCE_MAX;
  }
};


template<class Real> 
class EuclideanNormSquared : public ReductionOp<Real> {
public: 
  void reduce( const Real &input, Real &output ) const {
    output = output*output + input;
  }

  void reduce( const volatile Real &input, volatile Real &output ) const {
    output = output*output + input;
  }

  Real initialValue() const {
    return 0;
  }

  EReductionType reductionType() const {
    return Elementwise::REDUCE_SUM;
  }

};

} // namespace Elementwise 
} // namespace ROL

#endif 
#endif
