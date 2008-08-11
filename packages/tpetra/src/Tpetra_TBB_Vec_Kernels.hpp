
#ifndef TPETRA_TBB_KERNELS_HPP
#define TPETRA_TBB_KERNELS_HPP

// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Tpetra_ConfigDefs.hpp"

#ifdef HAVE_TPETRA_TBB

#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

namespace Tpetra {

//--------------------------------
template<typename OrdinalType, typename ScalarType>
struct update_op {
  ScalarType* y;
  ScalarType scalY;
  const ScalarType* x;
  ScalarType scalX;

  void operator()(const tbb::blocked_range<OrdinalType>& range) const {
    for(int i=range.begin(); i!=range.end(); ++i) {
      y[i] = scalY*y[i] + scalX*x[i];
    }
  }

};//struct update_op

//--------------------------------
template<typename OrdinalType, typename ScalarType>
struct dot_op {
  dot_op(const ScalarType* v1, const ScalarType* v2)
    : y(v1), x(v2), dot(0.0) {}
  dot_op(const dot_op& op, tbb::split)
    : y(op.y), x(op.x), dot(0.0) {}
  virtual ~dot_op(){}

  const ScalarType* y;
  const ScalarType* x;
  ScalarType dot;

  void operator()(const tbb::blocked_range<OrdinalType>& range) {
    for(int i=range.begin(); i!=range.end(); ++i) {
      dot += y[i] * x[i];
    }
  }

  void join(const dot_op& op) {dot += op.dot;}

};//struct update_op

//--------------------------------
//! perform axpy-like operation Y = scalarY*Y + scalarX*X
template<typename OrdinalType, typename ScalarType>
void threaded_vector_update(OrdinalType length,
                            ScalarType scalarY,
                            ScalarType* Y,
                            ScalarType scalarX,
                            const ScalarType* X)
{
  update_op<OrdinalType,ScalarType> op;
  op.y = Y;
  op.scalY = scalarY;
  op.x = X;
  op.scalX = scalarX;

  tbb::parallel_for(tbb::blocked_range<OrdinalType>(0, length, length/8), op);
}

//--------------------------------
//! perform dot-product operation result = sum(y[i]*x[i])
template<typename OrdinalType, typename ScalarType>
ScalarType threaded_vector_dot(OrdinalType length,
                         const ScalarType* Y,
                         const ScalarType* X)
{
  dot_op<OrdinalType,ScalarType> op(X, Y);

  tbb::parallel_reduce(tbb::blocked_range<OrdinalType>(0, length, length/8), op);

  return op.dot;
}

}//namespace Tpetra

#endif

#endif

