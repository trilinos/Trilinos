// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 *  Copyright 2008-2009 NVIDIA Corporation
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

/*! \file polynomial.h
 *  \brief polynomial relaxation.
 */

#pragma once

#include <cusp/detail/config.h>

#include <cusp/linear_operator.h>

namespace cusp
{
namespace precond
{
namespace aggregation
{
// forward definitions
template<typename MatrixType> struct sa_level;
} // end namespace aggregation
} // end namespace precond

namespace relaxation
{

template <typename ValueType, typename MemorySpace, typename Orientation>
class block_polynomial : public cusp::linear_operator<ValueType, MemorySpace>
{
public:

  typedef Orientation orientation;

  // note: default_coefficients lives on the host
  cusp::array1d<ValueType, cusp::host_memory> default_coefficients;
  cusp::array2d<ValueType, MemorySpace, Orientation> residual;
  cusp::array2d<ValueType, MemorySpace, Orientation> h;
  cusp::array2d<ValueType, MemorySpace, Orientation> y;

  block_polynomial();

  template <typename MatrixType, typename VectorType>
  block_polynomial(const MatrixType& A, const VectorType& coefficients);

  template <typename MemorySpace2>
  block_polynomial(const block_polynomial<ValueType,MemorySpace2,Orientation>& A);

  template <typename MatrixType>
  block_polynomial(const cusp::precond::aggregation::sa_level<MatrixType>& sa_level);

  // ignores initial x
  template<typename MatrixType, typename VectorType1, typename VectorType2>
  void presmooth(const MatrixType& A, const VectorType1& b, VectorType2& x);

  // smooths initial x
  template<typename MatrixType, typename VectorType1, typename VectorType2>
  void postsmooth(const MatrixType& A, const VectorType1& b, VectorType2& x);

  template <typename MatrixType, typename VectorType1, typename VectorType2>
  void operator()(const MatrixType& A, const VectorType1& b, VectorType2& x) const;

  template <typename MatrixType, typename VectorType1, typename VectorType2, typename VectorType3>
  void operator()(const MatrixType& A, const VectorType1& b, VectorType2& x, VectorType3& coeffients);
};

} // end namespace relaxation
} // end namespace cusp

#include <cusp/relaxation/block_polynomial.inl>
