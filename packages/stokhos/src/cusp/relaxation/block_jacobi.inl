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

/*! \file jacobi.inl
 *  \brief Inline file for jacobi.h
 */
#include <cusp/array2d.h>
#include <cusp/multiply.h>
#include <cusp/detail/format_utils.h>

#include <thrust/functional.h>
#include <thrust/transform.h>
#include <thrust/iterator/zip_iterator.h>

namespace cusp
{
namespace relaxation
{
namespace detail
{

template <typename ValueType>
struct block_jacobi_presmooth_functor
{
    ValueType omega;

    block_jacobi_presmooth_functor(ValueType omega) : omega(omega) {}

    __host__ __device__
    ValueType operator()(const ValueType& b, const ValueType& d) const
    {
        return omega * b / d;
    }
};

template <typename ValueType>
struct block_jacobi_postsmooth_functor
{
    ValueType omega;

    block_jacobi_postsmooth_functor(ValueType omega) : omega(omega) {}

    template <typename Tuple>
    __host__ __device__
    ValueType operator()(const Tuple& t)
    {
        const ValueType x = thrust::get<0>(t);
        const ValueType d = thrust::get<1>(t);
        const ValueType b = thrust::get<2>(t);
        const ValueType y = thrust::get<3>(t);

        return x + omega * (b - y) / d;
    }
};

} // end namespace detail


// constructor
template <typename ValueType, typename MemorySpace>
    block_jacobi<ValueType,MemorySpace>
    ::block_jacobi() : default_omega(0.0)
    {
    }

template <typename ValueType, typename MemorySpace>
template<typename MemorySpace2>
    block_jacobi<ValueType,MemorySpace>
    ::block_jacobi(const block_jacobi<ValueType,MemorySpace2>& A)
        : default_omega(A.default_omega), temp(A.temp), diagonal(A.diagonal)
    {
    }

template <typename ValueType, typename MemorySpace>
template<typename MatrixType>
    block_jacobi<ValueType,MemorySpace>
    ::block_jacobi(const MatrixType& A, ValueType omega)
        : default_omega(omega), temp(A.num_rows, 10)
    {
        CUSP_PROFILE_SCOPED();

        // extract the main diagonal
        cusp::detail::extract_diagonal(A, diagonal);
    }

template <typename ValueType, typename MemorySpace>
template<typename MatrixType>
    block_jacobi<ValueType,MemorySpace>
    ::block_jacobi(const cusp::precond::aggregation::sa_level<MatrixType>& sa_level, ValueType weight)
        : default_omega(weight/sa_level.rho_DinvA), temp(sa_level.A_.num_rows, 10)
    {
        CUSP_PROFILE_SCOPED();

        // extract the main diagonal
        cusp::detail::extract_diagonal(sa_level.A_, diagonal);
    }

// linear_operator
template <typename ValueType, typename MemorySpace>
template<typename MatrixType, typename VectorType1, typename VectorType2>
    void block_jacobi<ValueType,MemorySpace>
    ::operator()(const MatrixType& A, const VectorType1& b, VectorType2& x)
    {
        block_jacobi<ValueType,MemorySpace>::operator()(A,b,x,default_omega);
    }

// override default omega
template <typename ValueType, typename MemorySpace>
template<typename MatrixType, typename VectorType1, typename VectorType2>
    void block_jacobi<ValueType,MemorySpace>
    ::operator()(const MatrixType& A, const VectorType1& b, VectorType2& x, ValueType omega)
    {
        CUSP_PROFILE_SCOPED();

        // temp <- A*x
        cusp::multiply(A, x, temp);
	
	for (int i = 0; i < x.num_cols; i++){        
        // x <- x + omega * D^-1 * (b - temp)
        thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(x.column(i).begin(), diagonal.begin(), b.column(i).begin(), temp.column(i).begin())),
                          thrust::make_zip_iterator(thrust::make_tuple(x.column(i).end(),   diagonal.end(),   b.column(i).end(),   temp.column(i).end())),
                          x.column(i).begin(),
                          detail::block_jacobi_postsmooth_functor<ValueType>(omega));
    
	}
    }

template <typename ValueType, typename MemorySpace>
template<typename MatrixType, typename VectorType1, typename VectorType2>
    void block_jacobi<ValueType,MemorySpace>
    ::presmooth(const MatrixType&, const VectorType1& b, VectorType2& x)
    {
        CUSP_PROFILE_SCOPED();
        for (int i = 0; i < x.num_cols; i++){

        // x <- omega * D^-1 * b 
        thrust::transform(b.column(i).begin(), b.column(i).end(),
                          diagonal.begin(),
                          x.column(i).begin(),
                          detail::block_jacobi_presmooth_functor<ValueType>(default_omega));
	}

    }

template <typename ValueType, typename MemorySpace>
template<typename MatrixType, typename VectorType1, typename VectorType2>
    void block_jacobi<ValueType,MemorySpace>
    ::postsmooth(const MatrixType& A, const VectorType1& b, VectorType2& x)
    {
        CUSP_PROFILE_SCOPED();

        // y <- A*x
        cusp::multiply(A, x, temp);
        for (int i = 0; i < x.num_cols; i++){
        // x <- x + omega * D^-1 * (b - y)
        thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(x.column(i).begin(), diagonal.begin(), b.column(i).begin(), temp.column(i).begin())),
                          thrust::make_zip_iterator(thrust::make_tuple(x.column(i).end(),   diagonal.end(),   b.column(i).end(),   temp.column(i).end())),
                          x.column(i).begin(),
                          detail::block_jacobi_postsmooth_functor<ValueType>(default_omega));
	}
    

    }



} // end namespace relaxation
} // end namespace cusp

