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

#ifndef STOKHOS_HOST_LEXICOGRAPHIC_BLOCK_SPARSE_3_TENSOR_HPP
#define STOKHOS_HOST_LEXICOGRAPHIC_BLOCK_SPARSE_3_TENSOR_HPP

#include "KokkosArray_Host.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_LexicographicBlockSparse3Tensor.hpp"

namespace Stokhos {

template< typename ValueType >
class Multiply< LexicographicBlockSparse3Tensor< ValueType , KokkosArray::Host > , void , void , DefaultSparseMatOps >
{
public:

  typedef KokkosArray::Host::size_type size_type ;
  typedef LexicographicBlockSparse3Tensor< ValueType , KokkosArray::Host > tensor_type ;

  template< typename MatrixValue , typename VectorValue >
  static void apply( const tensor_type & tensor ,
                     const MatrixValue * const a ,
                     const VectorValue * const x ,
                           VectorValue * const y )
  {
    // const int max_size = 10;
    // MatrixValue ax[max_size][max_size];

    const size_type nBlock = tensor.num_coord();

    // Loop over coordinate blocks
    size_type value_entry = 0;
    for ( size_type block = 0; block < nBlock; ++block) {
      const size_type i_begin = tensor.get_i_begin(block);
      const size_type j_begin = tensor.get_j_begin(block);
      const size_type k_begin = tensor.get_k_begin(block);
      const size_type i_size = tensor.get_i_size(block);
      const size_type j_size = tensor.get_j_size(block);
      const size_type k_size = tensor.get_k_size(block);
      VectorValue * const y_block = y + i_begin;
      const MatrixValue * const a_block = a + j_begin;
      const VectorValue * const x_block = x + k_begin;

      // // Precompute a*x outer product
      // for (size_type j=0; j<j_size; ++j) {
      //   for (size_type k=0; k<k_size; ++k) {
      //     ax[j][k] = a_block[j]*x_block[k]; 
      //   }
      // }

      /*
      // Compute y_i = \sum_{j,k} c_{ijk} * a_j * x_k
      for (size_type i=0; i<i_size; ++i) {
        VectorValue ytmp = 0;
        for (size_type j=0; j<j_size; ++j) {
          const size_type imj = i-j;
          const size_type ipj = i+j+1;
          const size_type k_beg = 0      <= imj ? imj    : -imj;
          const size_type k_end = k_size <= ipj ? k_size :  ipj;
          const size_type k0 = k_beg % 2 == (i+j) % 2 ? k_beg : k_beg+1;
          for (size_type k=k0; k<k_end; ++k) {
            //ytmp += tensor.value(value_entry++) * ax[j][k];
            ytmp += tensor.value(value_entry++) * ( a_block[j] * x_block[k] );
          }
        }
        y_block[i] += ytmp ;
      }
      */

      // Compute y_i = \sum_{j,k} c_{ijk} * a_j * x_k
      for (size_type i=0; i<i_size; ++i) {
        VectorValue ytmp = 0;
        for (size_type j=0; j<j_size; ++j) {
          for (size_type k=((i+j)%2); k<k_size; k+=2) {
            ytmp += tensor.value(value_entry++) * ( a_block[j] * x_block[k] );
          }
        }
        y_block[i] += ytmp ;
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

#endif /* #ifndef STOKHOS_HOST_LEXICOGRAPHIC_BLOCK_SPARSE_3_TENSOR_HPP */
