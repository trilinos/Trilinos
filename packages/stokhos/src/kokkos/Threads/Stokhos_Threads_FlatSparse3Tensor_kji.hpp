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

#ifndef STOKHOS_THREADS_FLAT_SPARSE_3_TENSOR_KJI_HPP
#define STOKHOS_THREADS_FLAT_SPARSE_3_TENSOR_KJI_HPP

#include "Kokkos_Threads.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_FlatSparse3Tensor_kji.hpp"

namespace Stokhos {

template< typename ValueType >
class Multiply< FlatSparse3Tensor_kji< ValueType , Kokkos::Threads > , void , void , DefaultSparseMatOps >
{
public:

  typedef Kokkos::Threads::size_type size_type ;
  typedef FlatSparse3Tensor_kji< ValueType , Kokkos::Threads > tensor_type ;

  template< typename MatrixValue , typename VectorValue >
  static void apply( const tensor_type & tensor ,
                     const MatrixValue * const a ,
                     const VectorValue * const x ,
                           VectorValue * const y )
  {
    const size_type nk = tensor.num_k();

    // Loop over k
    for ( size_type k = 0; k < nk; ++k) {
      const MatrixValue ak = a[k];
      const VectorValue xk = x[k];

      // Loop over j for this k
      const size_type nj = tensor.num_j(k);
      const size_type jBeg = tensor.j_begin(k);
      const size_type jEnd = jBeg + nj;
      for (size_type jEntry = jBeg; jEntry < jEnd; ++jEntry) {
        const size_type j = tensor.j_coord(jEntry);
        VectorValue tmp = a[j] * xk + ak * x[j];

        // Loop over i for this k,j
        const size_type ni = tensor.num_i(jEntry);
        const size_type iBeg = tensor.i_begin(jEntry);
        const size_type iEnd = iBeg + ni;
        for (size_type iEntry = iBeg; iEntry < iEnd; ++iEntry) {
          const size_type i = tensor.i_coord(iEntry);
          y[i] += tensor.value(iEntry) * tmp;
        }
      }
    }
  }

  static size_type matrix_size( const tensor_type & tensor )
  { return tensor.dimension(); }

  static size_type vector_size( const tensor_type & tensor )
  { return tensor.dimension(); }
};

//----------------------------------------------------------------------------

} // namespace Stokhos

#endif /* #ifndef STOKHOS_THREADS_SPARSEPRODUCTTENSOR_KJI_HPP */
