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

/*! \file multilevel.h
 *  \brief Multilevel hierarchy
 *
 */

#pragma once

#include <cusp/detail/config.h>
#include <cusp/detail/lu.h>
#include <cusp/array2d.h>
#include <cusp/array1d.h>
#include <cusp/linear_operator.h>

namespace cusp
{

template <typename MatrixType, typename SmootherType, typename SolverType>
class block_multilevel :
  public cusp::linear_operator<typename MatrixType::value_type,
                               typename MatrixType::memory_space>
{
public:

  typedef typename MatrixType::index_type IndexType;
  typedef typename MatrixType::value_type ValueType;
  typedef typename MatrixType::memory_space MemorySpace;
  typedef typename SmootherType::orientation Orientation;

  struct level
  {
    MatrixType R;  // restriction operator
    MatrixType A;  // matrix
    MatrixType P;  // prolongation operator
    cusp::array2d<ValueType,MemorySpace,Orientation> x;   // per-level solution
    cusp::array2d<ValueType,MemorySpace,Orientation> b;   // per-level rhs
    cusp::array2d<ValueType,MemorySpace,Orientation> residual;  // per-level residual

    SmootherType smoother;

    level(){}

    template<typename Level_Type>
    level(const Level_Type& level) : R(level.R), A(level.A), P(level.P), x(level.x), b(level.b), residual(level.residual), smoother(level.smoother){}
  };

  SolverType solver;

  std::vector<level> levels;

  block_multilevel(){};

  template <typename MatrixType2, typename SmootherType2, typename SolverType2>
  block_multilevel(const block_multilevel<MatrixType2, SmootherType2, SolverType2>& M);

  template <typename Array1, typename Array2>
  void operator()(const Array1& x, Array2& y);

  template <typename Array1, typename Array2>
  void solve(const Array1& b, Array2& x);

  template <typename Array1, typename Array2, typename Monitor>
  void solve(const Array1& b, Array2& x, Monitor& monitor);

  void print( void );

  double operator_complexity( void );

  double grid_complexity( void );

protected:

  template <typename Array1, typename Array2>
  void _solve(const Array1& b, Array2& x, const size_t i);
};
/*! \}
 */

} // end namespace cusp

#include <cusp/block_multilevel.inl>
