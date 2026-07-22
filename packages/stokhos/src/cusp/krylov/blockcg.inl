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

#include <cusp/array2d.h>
#include <cusp/array1d.h>
#include <cusp/MVmultiply.h>
#include <cusp/block_monitor.h>
#include <cusp/linear_operator.h>
#include <cusp/elementwise.h>



namespace cusp
{
namespace krylov
{

template <class LinearOperator,
          class Vector>
void blockcg(LinearOperator& A,
             Vector& x,
             Vector& b)
{
  typedef typename LinearOperator::value_type   ValueType;

  cusp::default_block_monitor<ValueType> monitor(b);

  cusp::krylov::blockcg(A, x, b, monitor);
}

template <class LinearOperator,
          class Vector,
          class Monitor>
void blockcg(LinearOperator& A,
             Vector& x,
             Vector& b,
             Monitor& monitor)
{
  typedef typename LinearOperator::value_type   ValueType;
  typedef typename LinearOperator::memory_space MemorySpace;

  cusp::identity_operator<ValueType,MemorySpace> M(A.num_rows, A.num_cols);

  cusp::krylov::blockcg(A, x, b, monitor, M);
}

template <class LinearOperator,
          class Vector,
          class Monitor,
          class Preconditioner>
void blockcg(LinearOperator& A,
             Vector& x,
             Vector& b,
             Monitor& monitor,
             Preconditioner& M)
{
  CUSP_PROFILE_SCOPED();

  typedef typename LinearOperator::value_type   ValueType;
  typedef typename LinearOperator::memory_space MemorySpace;
  typedef typename LinearOperator::index_type   IndexType;
  typedef typename Vector::orientation          Orientation;

  assert(A.num_rows == A.num_cols);        // sanity check

  const size_t N = A.num_rows;

  const size_t numRHS = b.num_cols;

  // allocate workspace
  cusp::array2d<ValueType,MemorySpace,Orientation> y(N, numRHS, 0);
  cusp::array2d<ValueType,MemorySpace,Orientation> r(N, numRHS, 0);
  cusp::array2d<ValueType,MemorySpace,Orientation> z(N, numRHS, 0);

  cusp::array1d<ValueType,MemorySpace> rTz(numRHS);
  cusp::array1d<ValueType,MemorySpace> rTz_old(numRHS);
  cusp::array1d<ValueType,MemorySpace> pAp(numRHS);

  cusp::array1d<ValueType,MemorySpace> alpha(numRHS);
  cusp::array1d<ValueType,MemorySpace> beta(numRHS);

  cusp::array2d<ValueType,MemorySpace,Orientation> p(N, numRHS,0);
  cusp::array2d<ValueType,MemorySpace,Orientation> tmp(N, numRHS,0);

  //Compute the initial alpha and (r^T)z
  // y <- Ax
  cusp::OVmultiply(A, x, y);
  // r <- b - A*x
  cusp::axpby_array(ValueType(1), b, ValueType(-1), y, r);

  cusp::axpby_array(ValueType(1.0), b, ValueType(-1.0), y, p);
  // z <- M*r
  cusp::MVmultiply(M, r, z);

  // p <- z
  cusp::copy(z, p);
  // rz = <r^H, z>
  //compute vector rHz where the components are the individual dot-products rHz[i] = r[i]^T * z[i] wherer[i] is the i-th column of r
  cusp::MVdot(r,z,rTz);

  //Iterate until moniter isfinished (if convg tol reached or max iters reached)
  while (!monitor.finished(r))
  {
    // y <- Ap
    cusp::OVmultiply(A, p, y);

    //pAp <- <p, Ap>
    cusp::MVdot(p, y, pAp);

    // Compute alpha <- <r,z>/<p,y> = (r^T)*z / (p^T)*A*p
    for (IndexType i = 0; i < numRHS; i++){
      alpha[i] = rTz[i] / pAp[i];
    }

    //Update the solution
    // x <- x + p *alpha

    cusp::MVmultiply(p, alpha, tmp);
    cusp::axpby_array(ValueType(1.0), x, ValueType(1.0), tmp, x);

    //Save the denomiator of beta before residual is updated to rHz_old
    cusp::copy(rTz, rTz_old);

    //Compute the new residual
    // r <- r - y * alpha = r - Ax * alpha
    cusp::MVmultiply(y, alpha, tmp);
    cusp::axpby_array(ValueType(1.0), r, ValueType(-1.0), tmp, r);


    //if precond
    // z <- M*r
    cusp::MVmultiply(M, r, z);
    //Compute beta <- new <r, z> / old <r, z>
    // rTz = <r, z>
    cusp::MVdot(r, z, rTz);
    for (IndexType i = 0; i < numRHS; i++){
      beta[i] = rTz[i] / rTz_old[i];
    }
    cusp::MVmultiply(p, beta, p);
    cusp::axpby_array(ValueType(1.0), z, ValueType(1.0), p, p);

    ++monitor;
  }
}

} // end namespace krylov
} // end namespace cusp

