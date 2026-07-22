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

/*! \file polynomial.inl
 *  \brief Inline file for polynomial.h
 */
#include <cusp/MVmultiply.h>
#include <cusp/multiply.h>

#include <cusp/detail/format_utils.h>
#include <cusp/detail/spectral_radius.h>

#include <math.h>

namespace cusp
{
namespace relaxation
{
namespace detail
{
template <typename ValueType>
void block_chebyshev_polynomial_coefficients(
  const ValueType rho,
  cusp::array1d<ValueType,cusp::host_memory>& coefficients,
  const ValueType lower_bound = 1.0/30.0,
  const ValueType upper_bound = 1.1)
{
  const size_t degree = 3;

  ValueType x0 = lower_bound * rho;
  ValueType x1 = upper_bound * rho;

  // Chebyshev roots for the interval [-1,1]
  cusp::array1d<ValueType,cusp::host_memory> std_roots(degree);

  for( size_t i=0; i<degree; i++ )
    std_roots[i] = std::cos( M_PI * (ValueType(i) + 0.5)/ degree );

  // Chebyshev roots for the interval [x0,x1]
  for( size_t i=0; i<degree; i++ )
    std_roots[i] = 0.5 * (x1-x0) * (1 + std_roots[i]) + x0;

  // Compute monic polynomial coefficients of polynomial with scaled roots
  // TODO: Implement convolution method for polynomial multiplication
  coefficients.resize(degree+1);
  ValueType a = std_roots[0];
  ValueType b = std_roots[1];
  ValueType c = std_roots[2];
  coefficients[0] = 1.0;
  coefficients[1] = -(a+b+c);
  coefficients[2] = (a*b) + (b*c) + (c*a);
  coefficients[3] = -(a*b*c);

  // Scale coefficients to enforce C(0) = 1.0
  ValueType scale_factor = 1.0/coefficients.back();
  cusp::blas::scal(coefficients, scale_factor);
}
}

// constructor
template <typename ValueType, typename MemorySpace, typename Orientation>
block_polynomial<ValueType,MemorySpace,Orientation>
::block_polynomial()
{
}

template <typename ValueType, typename MemorySpace, typename Orientation>
template<typename MatrixType, typename VectorType>
block_polynomial<ValueType,MemorySpace,Orientation>
::block_polynomial(const MatrixType& A, const VectorType& coefficients)
{
  size_t default_size = coefficients.size()-1;
  default_coefficients.resize( default_size );
  for( size_t index = 0; index < default_size; index++ )
    default_coefficients[index] = -coefficients[index];
}

template <typename ValueType, typename MemorySpace, typename Orientation>
template<typename MemorySpace2>
block_polynomial<ValueType,MemorySpace,Orientation>
::block_polynomial(const block_polynomial<ValueType,MemorySpace2,Orientation>& A) : default_coefficients(A.default_coefficients), residual(A.residual), h(A.h), y(A.y)
{
}

template <typename ValueType, typename MemorySpace, typename Orientation>
template<typename MatrixType>
block_polynomial<ValueType,MemorySpace,Orientation>
::block_polynomial(const cusp::precond::aggregation::sa_level<MatrixType>& sa_level)
{
  CUSP_PROFILE_SCOPED();

  ValueType rho = cusp::detail::ritz_spectral_radius_symmetric(sa_level.A_, 8);
  detail::block_chebyshev_polynomial_coefficients(rho, default_coefficients);
  default_coefficients.resize( default_coefficients.size() - 1 );

  for( size_t index = 0; index < default_coefficients.size(); index++ )
    default_coefficients[index] *= -1;

}

// linear_operator
template <typename ValueType, typename MemorySpace, typename Orientation>
template<typename MatrixType, typename VectorType1, typename VectorType2>
void block_polynomial<ValueType,MemorySpace,Orientation>
::operator()(const MatrixType& A, const VectorType1& b, VectorType2& x) const
{
  CUSP_PROFILE_SCOPED();

  block_polynomial<ValueType,MemorySpace,Orientation>::operator()(A,b,x,default_coefficients);
}

template <typename ValueType, typename MemorySpace, typename Orientation>
template<typename MatrixType, typename VectorType1, typename VectorType2>
void block_polynomial<ValueType,MemorySpace,Orientation>
::presmooth(const MatrixType& A, const VectorType1& b, VectorType2& x)
{
  CUSP_PROFILE_SCOPED();
  int N = A.num_rows;
  int numRHS = x.num_cols;
  residual.resize(N, numRHS);
  y.resize(N, numRHS);
  h.resize(N, numRHS);

  // Ignore the initial x and use b as the residual
  ValueType scale_factor = default_coefficients[0];
  //x <- scale_factor*b
  cusp::axpby_array(scale_factor, b, ValueType(0), x, x);
  for( size_t i = 1; i<default_coefficients.size(); i++ )
  {
    scale_factor = default_coefficients[i];
    cusp::MVmultiply(A, x, y);
    //x <- y + scale_factor*b
    cusp::axpby_array(scale_factor, b, ValueType(1.0), y, x);
  }
}

template <typename ValueType, typename MemorySpace, typename Orientation>
template<typename MatrixType, typename VectorType1, typename VectorType2>
void block_polynomial<ValueType,MemorySpace,Orientation>
::postsmooth(const MatrixType& A, const VectorType1& b, VectorType2& x)
{
  CUSP_PROFILE_SCOPED();

  // compute residual <- b - A*x
  cusp::MVmultiply(A, x, residual);
  cusp::axpby_array(ValueType(1), b, ValueType(-1), residual, residual);
  ValueType scale_factor = default_coefficients[0];
  cusp::axpby_array(scale_factor, residual, ValueType(0), h, h);
  for( size_t i = 1; i<default_coefficients.size(); i++ )
  {
    scale_factor = default_coefficients[i];

    cusp::MVmultiply(A, h, y);
    cusp::axpby_array(ValueType(1.0), y, scale_factor, residual, h);
  }

  cusp::axpby_array(ValueType(1.0), h, ValueType(1.0), x, x);
}

// override default coefficients
template <typename ValueType, typename MemorySpace, typename Orientation>
template<typename MatrixType, typename VectorType1, typename VectorType2, typename VectorType3>
void block_polynomial<ValueType,MemorySpace,Orientation>
::operator()(const MatrixType& A, const VectorType1& b, VectorType2& x, VectorType3& coefficients)
{
  if( cusp::blas::nrm2(x.column(0)) == 0.0 )
  {
    //residual = b;
    cusp::copy(b,residual);
  }
  else
  {
    // compute residual <- b - A*x
    cusp::MVmultiply(A, x, residual);
    cusp::blas::axpby(b.values, residual.values, residual.values, ValueType(1),
                      ValueType(-1));
  }
  
  ValueType scale_factor = coefficients[0];
  cusp::blas::axpby(residual.values, h.values, h.values, scale_factor,
                    ValueType(0));

  for( size_t i = 1; i<coefficients.size(); i++ )
  {
    scale_factor = coefficients[i];

    cusp::MVmultiply(A, h, y);
    cusp::blas::axpby(y.values, residual.values, h.values, ValueType(1.0),
                      scale_factor);
  }

  cusp::blas::axpy(h.values, x.values, ValueType(1.0));
}

} // end namespace relaxation
} // end namespace cusp

