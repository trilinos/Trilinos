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
#include <iostream>
#include <cusp/blas.h>
#include <cusp/MVmultiply.h>
#include <cusp/multiply.h>
#include <cusp/block_monitor.h>
#include <cusp/array1d.h>
#include <cusp/array2d.h>
#include <cusp/elementwise.h>

#include "Teuchos_TimeMonitor.hpp"

namespace cusp
{

template <typename MatrixType, typename SmootherType, typename SolverType>
template <typename MatrixType2, typename SmootherType2, typename SolverType2>
block_multilevel<MatrixType,SmootherType,SolverType>
::block_multilevel(const block_multilevel<MatrixType2,SmootherType2,SolverType2>& M)
{
  for( size_t lvl = 0; lvl < M.levels.size(); lvl++ )
    levels.push_back(M.levels[lvl]);
}

template <typename MatrixType, typename SmootherType, typename SolverType>
template <typename Array1, typename Array2>
void block_multilevel<MatrixType,SmootherType,SolverType>
::operator()(const Array1& b, Array2& x)
{
  CUSP_PROFILE_SCOPED();
  // perform 1 V-cycle
  _solve(b, x, 0);

}

template <typename MatrixType, typename SmootherType, typename SolverType>
template <typename Array1, typename Array2>
void block_multilevel<MatrixType,SmootherType,SolverType>
::solve(const Array1& b, Array2& x)
{
  CUSP_PROFILE_SCOPED();

  cusp::default_block_monitor<ValueType> monitor(b);

  solve(b, x, monitor);
}

template <typename MatrixType, typename SmootherType, typename SolverType>
template <typename Array1, typename Array2, typename Monitor>
void block_multilevel<MatrixType,SmootherType,SolverType>
::solve(const Array1& b, Array2& x, Monitor& monitor)
{
  CUSP_PROFILE_SCOPED();
  const MatrixType& A = levels[0].A;
  const size_t n = A.num_rows;
  const size_t numRHS = b.num_cols;

  // use simple iteration
  cusp::array2d<ValueType,MemorySpace,Orientation> update(n, numRHS);
  cusp::array2d<ValueType,MemorySpace,Orientation> residual(n, numRHS);

  // compute initial residual
  cusp::MVmultiply(A, x, residual);
  cusp::subtract(b, residual, residual);
  while(!monitor.finished(residual))
  {
    _solve(residual, update, 0);

    // x += M * r
    cusp::add(update, x, x);

    // update residual
    cusp::MVmultiply(A, x, residual);
    cusp::axpby_array(ValueType(1.0), b, ValueType(-1.0), residual, residual);
    ++monitor;
  }
}

template <typename MatrixType, typename SmootherType, typename SolverType>
template <typename Array1, typename Array2>
void block_multilevel<MatrixType,SmootherType,SolverType>
::_solve(const Array1& b, Array2& x, const size_t i)
{
  TEUCHOS_FUNC_TIME_MONITOR("CUSP Block Multilevel Solve");

  CUSP_PROFILE_SCOPED();
  const size_t numRHS = b.num_cols;

  if (i + 1 == levels.size())
  {
    TEUCHOS_FUNC_TIME_MONITOR_DIFF("CUSP Coarse-grid Solve", coarse_grid);

    // coarse grid solve
    // copy to host first
    cusp::array2d<ValueType,cusp::host_memory,Orientation> temp_b(b);
    cusp::array2d<ValueType,cusp::host_memory,Orientation> temp_x(x.num_rows, numRHS);
    solver(temp_b, temp_x);
    cusp::copy(temp_x, x);
  }
  else
  {
    const MatrixType& A = levels[i].A;
    // presmooth
    levels[i].smoother.presmooth(A, b, x);

    // compute residual <- b - A*x
    cusp::MVmultiply(A, x, levels[i].residual);
    cusp::axpby_array(ValueType(1.0), b, ValueType(-1.0), levels[i].residual, levels[i].residual);

    // restrict to coarse grid
    cusp::MVmultiply(levels[i].R, levels[i].residual, levels[i + 1].b);

    // compute coarse grid solution
    _solve(levels[i + 1].b, levels[i + 1].x, i + 1);

    // apply coarse grid correction
    cusp::MVmultiply(levels[i].P, levels[i + 1].x, levels[i].residual);
    cusp::axpby_array(ValueType(1.0), levels[i].residual, ValueType(1.0), x, x);

    // postsmooth
    levels[i].smoother.postsmooth(A, b, x);
  }
}

template <typename MatrixType, typename SmootherType, typename SolverType>
void block_multilevel<MatrixType,SmootherType,SolverType>
::print( void )
{
  size_t num_levels = levels.size();

  std::cout << "\tNumber of Levels:\t" << num_levels << std::endl;
  std::cout << "\tOperator Complexity:\t" << operator_complexity() << std::endl;
  std::cout << "\tGrid Complexity:\t" << grid_complexity() << std::endl;
  std::cout << "\tlevel\tunknowns\tnonzeros:\t" << std::endl;

  double nnz = 0;

  for(size_t index = 0; index < num_levels; index++)
    nnz += levels[index].A.num_entries;

  for(size_t index = 0; index < num_levels; index++)
  {
    double percent = levels[index].A.num_entries / nnz;
    std::cout << "\t" << index << "\t" << levels[index].A.num_cols << "\t\t" \
              << levels[index].A.num_entries << " \t[" << 100*percent << "%]" \
              << std::endl;
  }
}

template <typename MatrixType, typename SmootherType, typename SolverType>
double block_multilevel<MatrixType,SmootherType,SolverType>
::operator_complexity( void )
{
  size_t nnz = 0;

  for(size_t index = 0; index < levels.size(); index++)
    nnz += levels[index].A.num_entries;

  return (double) nnz / (double) levels[0].A.num_entries;
}

template <typename MatrixType, typename SmootherType, typename SolverType>
double block_multilevel<MatrixType,SmootherType,SolverType>
::grid_complexity( void )
{
  size_t unknowns = 0;
  for(size_t index = 0; index < levels.size(); index++)
    unknowns += levels[index].A.num_rows;

  return (double) unknowns / (double) levels[0].A.num_rows;
}

} // end namespace cusp

