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

#ifndef STOKHOS_THREADS_COO_PRODUCT_TENSOR_HPP
#define STOKHOS_THREADS_COO_PRODUCT_TENSOR_HPP

#include "Kokkos_Threads.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_CooProductTensor.hpp"

namespace Stokhos {

template< typename ValueType, bool Pack >
class Multiply< CooProductTensor< ValueType, Kokkos::Threads, Pack >,
                void, void, DefaultSparseMatOps >
{
public:

  typedef Kokkos::Threads::size_type size_type ;
  typedef CooProductTensor< ValueType , Kokkos::Threads, Pack > tensor_type ;

  template< typename MatrixValue , typename VectorValue >
  static void apply( const tensor_type & tensor ,
                     const MatrixValue * const a ,
                     const VectorValue * const x ,
                           VectorValue * const y )
  {
    const size_type nEntry = tensor.entry_count();
    size_type i = 0, j = 0, k = 0, i_prev = -1;
    VectorValue val = 0.0, carry_val = 0.0;
    for ( size_type entry = 0 ; entry < nEntry ; ++entry ) {
      tensor.coord(entry, i, j, k);
      val = tensor.value(entry) * ( a[j] * x[k] + a[k] * x[j] );
      if (i == i_prev)
        carry_val += val;
      else {
        y[i_prev] += carry_val;
        carry_val = val;
      }
      i_prev = i;
    }
    y[i] += carry_val;
  }

  static size_type matrix_size( const tensor_type & tensor )
  { return tensor.dimension(); }

  static size_type vector_size( const tensor_type & tensor )
  { return tensor.dimension(); }
};

//----------------------------------------------------------------------------

} // namespace Stokhos

#endif /* #ifndef STOKHOS_THREADS_CRS_PRODUCT_TENSOR_HPP */
