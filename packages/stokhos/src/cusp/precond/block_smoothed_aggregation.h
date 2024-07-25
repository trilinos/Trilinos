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

/*! \file smoothed_aggregation.h
 *  \brief Algebraic multigrid preconditoner based on smoothed aggregation.
 *
 */

#pragma once

#include <cusp/detail/config.h>

#include <vector> // TODO replace with host_vector
#include <cusp/linear_operator.h>

#include <cusp/coo_matrix.h>
#include <cusp/csr_matrix.h>
#include <cusp/hyb_matrix.h>
#include <cusp/block_multilevel.h>
#include <cusp/relaxation/block_jacobi.h>
#include <cusp/relaxation/block_polynomial.h>

#include <cusp/detail/spectral_radius.h>
#include <cusp/detail/block_lu.h>

namespace cusp
{
namespace precond
{
namespace aggregation
{

/*! \addtogroup preconditioners Preconditioners
 *  \ingroup preconditioners
 *  \{
 */

template <typename IndexType, typename ValueType, typename MemorySpace>
struct amg_container {};

template <typename IndexType, typename ValueType>
struct amg_container<IndexType,ValueType,cusp::host_memory>
{
  // use CSR on host
  typedef typename cusp::csr_matrix<IndexType,ValueType,cusp::host_memory> setup_type;
  typedef typename cusp::csr_matrix<IndexType,ValueType,cusp::host_memory> solve_type;
};

template <typename IndexType, typename ValueType>
struct amg_container<IndexType,ValueType,cusp::device_memory>
{
  // use COO on device
//    typedef typename cusp::coo_matrix<IndexType,ValueType,cusp::device_memory> setup_type;
//    typedef typename cusp::hyb_matrix<IndexType,ValueType,cusp::device_memory> solve_type;
  typedef typename cusp::csr_matrix<IndexType,ValueType,cusp::device_memory> setup_type;
  typedef typename cusp::csr_matrix<IndexType,ValueType,cusp::device_memory> solve_type;

};

template<typename MatrixType>
struct sa_level
{
  typedef typename MatrixType::index_type IndexType;
  typedef typename MatrixType::value_type ValueType;
  typedef typename MatrixType::memory_space MemorySpace;

  MatrixType A_;                                        // matrix
  cusp::array1d<IndexType,MemorySpace> aggregates;      // aggregates
  cusp::array1d<ValueType,MemorySpace> B;               // near-nullspace candidates

  ValueType rho_DinvA;

  sa_level() {}

  template<typename SA_Level_Type>
  sa_level(const SA_Level_Type& sa_level) : A_(sa_level.A_), aggregates(sa_level.aggregates), B(sa_level.B), rho_DinvA(sa_level.rho_DinvA) {}
};


/*! \p smoothed_aggregation : algebraic multigrid preconditoner based on
 *  smoothed aggregation
 *
 */
//typename SmootherType = cusp::relaxation::block_polynomial<ValueType,MemorySpace>
template <typename IndexType, typename ValueType, typename MemorySpace,
          typename SmootherType,
          typename SolverType = cusp::detail::block_lu_solver<ValueType,cusp::host_memory> >
class block_smoothed_aggregation : public cusp::block_multilevel< typename amg_container<IndexType,ValueType,MemorySpace>::solve_type, SmootherType, SolverType>
{

  typedef typename amg_container<IndexType,ValueType,MemorySpace>::setup_type SetupMatrixType;
  typedef typename amg_container<IndexType,ValueType,MemorySpace>::solve_type SolveMatrixType;
  typedef typename cusp::block_multilevel<SolveMatrixType,SmootherType,SolverType> Parent;

public:

  ValueType theta;
  IndexType numRHS;
  std::vector< sa_level<SetupMatrixType> > sa_levels;

  template <typename MatrixType>
  block_smoothed_aggregation(const MatrixType& A, const IndexType numRHS, const ValueType theta=0);

  template <typename MatrixType, typename ArrayType>
  block_smoothed_aggregation(const MatrixType& A, const ArrayType& B, const IndexType numRHS, const ValueType theta=0);

  template <typename MemorySpace2,typename SmootherType2,typename SolverType2>
  block_smoothed_aggregation(const block_smoothed_aggregation<IndexType,ValueType,MemorySpace2,SmootherType2,SolverType2>& M);

protected:

  template <typename MatrixType, typename ArrayType>
  void init(const MatrixType& A, const ArrayType& B);

  void extend_hierarchy(void);
};
/*! \}
 */

} // end namespace aggregation
} // end namespace precond
} // end namespace cusp

#include <cusp/precond/block_smoothed_aggregation.inl>
