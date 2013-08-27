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

#ifndef STOKHOS_THREADS_LEXICOGRAPHIC_BLOCK_SPARSE_3_TENSOR_HPP
#define STOKHOS_THREADS_LEXICOGRAPHIC_BLOCK_SPARSE_3_TENSOR_HPP

#include "Kokkos_Threads.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_LexicographicBlockSparse3Tensor.hpp"

namespace Stokhos {

template< typename ValueType >
class Multiply< LexicographicBlockSparse3Tensor< ValueType , Kokkos::Threads > , void , void , DefaultSparseMatOps >
{
public:

  typedef Kokkos::Threads::size_type size_type ;
  typedef LexicographicBlockSparse3Tensor< ValueType , Kokkos::Threads > tensor_type ;

  template< typename MatrixValue , typename VectorValue >
  static void apply( const tensor_type & tensor ,
                     const MatrixValue * const a ,
                     const VectorValue * const x ,
                           VectorValue * const y )
  {
    const size_type nBlock = tensor.num_coord();

    if (tensor.symmetric()) {

      // Loop over coordinate blocks
      size_type value_entry = 0;
      for ( size_type block = 0; block < nBlock; ++block) {
        const int i_begin = tensor.get_i_begin(block);
        const int j_begin = tensor.get_j_begin(block);
        const int k_begin = tensor.get_k_begin(block);
        const int p_i = tensor.get_p_i(block);
        const int p_j = tensor.get_p_j(block);
        const int p_k = tensor.get_p_k(block);
        VectorValue * const y_block = y + i_begin;
        const MatrixValue * const a_j_block = a + j_begin;
        const VectorValue * const x_k_block = x + k_begin;
        const MatrixValue * const a_k_block = a + k_begin;
        const VectorValue * const x_j_block = x + j_begin;

        // for (int i=0; i<=p_i; ++i) {
        //   VectorValue ytmp = 0;
        //   for (int j=0; j<=p_j; ++j) {
        //     int k0 = j_eq_k != 0 ? j : 0;
        //     if (symmetric && (k0 % 2 != (i+j) % 2)) ++k0;
        //     for (int k=k0; k<=p_k; k+=k_inc) {
        //       ytmp += tensor.value(value_entry++) *
        //         ( a_j_block[j] * x_k_block[k] + a_k_block[k] * x_j_block[j] );
        //     }
        //   }
        //   y_block[i] += ytmp ;
        // }

        if (tensor.get_j_eq_k(block) != 0) {
          for (int i=0; i<=p_i; ++i) {
            VectorValue ytmp = 0;
            for (int j=0; j<=p_j; ++j) {
              int k0 = j%2 != (i+j)%2 ? j+1 : j;
              for (int k=k0; k<=p_k; k+=2) {
                ytmp += tensor.value(value_entry++) *
                  ( a_j_block[j] * x_k_block[k] + a_k_block[k] * x_j_block[j] );
              }
            }
            y_block[i] += ytmp ;
          }
        }
        else {
          for (int i=0; i<=p_i; ++i) {
            VectorValue ytmp = 0;
            for (int j=0; j<=p_j; ++j) {
              for (int k=(i+j)%2; k<=p_k; k+=2) {
                ytmp += tensor.value(value_entry++) *
                  ( a_j_block[j] * x_k_block[k] + a_k_block[k] * x_j_block[j] );
              }
            }
            y_block[i] += ytmp ;
          }
        }
      }

    }

    else {

      // Loop over coordinate blocks
      size_type value_entry = 0;
      for ( size_type block = 0; block < nBlock; ++block) {
        const int i_begin = tensor.get_i_begin(block);
        const int j_begin = tensor.get_j_begin(block);
        const int k_begin = tensor.get_k_begin(block);
        const int p_i = tensor.get_p_i(block);
        const int p_j = tensor.get_p_j(block);
        const int p_k = tensor.get_p_k(block);
        VectorValue * const y_block = y + i_begin;
        const MatrixValue * const a_j_block = a + j_begin;
        const VectorValue * const x_k_block = x + k_begin;
        const MatrixValue * const a_k_block = a + k_begin;
        const VectorValue * const x_j_block = x + j_begin;

        // for (int i=0; i<=p_i; ++i) {
        //   VectorValue ytmp = 0;
        //   for (int j=0; j<=p_j; ++j) {
        //     int k0 = j_eq_k != 0 ? j : 0;
        //     if (symmetric && (k0 % 2 != (i+j) % 2)) ++k0;
        //     for (int k=k0; k<=p_k; k+=k_inc) {
        //       ytmp += tensor.value(value_entry++) *
        //         ( a_j_block[j] * x_k_block[k] + a_k_block[k] * x_j_block[j] );
        //     }
        //   }
        //   y_block[i] += ytmp ;
        // }

        if (tensor.get_j_eq_k(block) != 0) {
          for (int i=0; i<=p_i; ++i) {
            VectorValue ytmp = 0;
            for (int j=0; j<=p_j; ++j) {
              for (int k=j; k<=p_k; ++k) {
                ytmp += tensor.value(value_entry++) *
                  ( a_j_block[j] * x_k_block[k] + a_k_block[k] * x_j_block[j] );
              }
            }
            y_block[i] += ytmp ;
          }
        }
        else {
          for (int i=0; i<=p_i; ++i) {
            VectorValue ytmp = 0;
            for (int j=0; j<=p_j; ++j) {
              for (int k=0; k<=p_k; ++k) {
                ytmp += tensor.value(value_entry++) *
                  ( a_j_block[j] * x_k_block[k] + a_k_block[k] * x_j_block[j] );
              }
            }
            y_block[i] += ytmp ;
          }
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

#endif /* #ifndef STOKHOS_THREADS_LEXICOGRAPHIC_BLOCK_SPARSE_3_TENSOR_HPP */
