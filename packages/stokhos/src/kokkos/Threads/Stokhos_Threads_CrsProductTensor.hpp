// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_THREADS_CRS_PRODUCT_TENSOR_HPP
#define STOKHOS_THREADS_CRS_PRODUCT_TENSOR_HPP

#include "Kokkos_Core.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_CrsProductTensor.hpp"
#include "Stokhos_BlockCrsMatrix.hpp"
#include "Stokhos_StochasticProductTensor.hpp"
#include "Stokhos_TinyVec.hpp"

namespace Stokhos {

template< typename ValueType , typename MatrixValue , typename VectorValue >
class Multiply<
  BlockCrsMatrix< StochasticProductTensor< ValueType, CrsProductTensor< ValueType , Kokkos::Threads > , Kokkos::Threads > , MatrixValue , Kokkos::Threads > ,
  Kokkos::View< VectorValue** , Kokkos::LayoutLeft , Kokkos::Threads > ,
  Kokkos::View< VectorValue** , Kokkos::LayoutLeft , Kokkos::Threads > >
{
public:
  typedef MultiplyImpl<ValueType, MatrixValue, VectorValue, Kokkos::Threads> impl_type;
  typedef typename impl_type::matrix_type matrix_type;
  typedef typename impl_type::block_vector_type block_vector_type;

  static void apply( const matrix_type & A ,
                     const block_vector_type & x ,
                     const block_vector_type & y ) {
    impl_type::apply(A,x,y);
  }
};

} // namespace Stokhos

#endif /* #ifndef STOKHOS_THREADS_CRS_PRODUCT_TENSOR_HPP */
