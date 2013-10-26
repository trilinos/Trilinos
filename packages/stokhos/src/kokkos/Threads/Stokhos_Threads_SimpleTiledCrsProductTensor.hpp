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

#ifndef STOKHOS_THREADS_SIMPLE_TILED_CRS_PRODUCT_TENSOR_HPP
#define STOKHOS_THREADS_SIMPLE_TILED_CRS_PRODUCT_TENSOR_HPP

#include "Kokkos_Threads.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_SimpleTiledCrsProductTensor.hpp"
#include "Stokhos_Threads_TinyVec.hpp"

namespace Stokhos {

template< typename ValueType >
class Multiply< SimpleTiledCrsProductTensor< ValueType , Kokkos::Threads > , void , void , DefaultSparseMatOps >
{
public:

  typedef Kokkos::Threads::size_type size_type ;
  typedef SimpleTiledCrsProductTensor< ValueType , Kokkos::Threads > tensor_type ;

  template< typename MatrixValue , typename VectorValue >
  static void apply( const tensor_type & tensor ,
                     const MatrixValue * const a ,
                     const VectorValue * const x ,
                           VectorValue * const y )
  {
    const size_type block_size = 2;
    typedef TinyVec<ValueType,block_size,false> TV;

    const size_type n_i_tile = tensor.num_i_tiles();
    for (size_type i_tile = 0; i_tile<n_i_tile; ++i_tile) {
      const size_type i_begin = tensor.i_begin(i_tile);
      const size_type i_size  = tensor.i_size(i_tile);

      const size_type n_j_tile = tensor.num_j_tiles(i_tile);
      for (size_type j_tile = 0; j_tile<n_j_tile; ++j_tile) {
        const size_type j_begin = tensor.j_begin(i_tile, j_tile);
        //const size_type j_size  = tensor.j_size(i_tile, j_tile);

        const size_type n_k_tile = tensor.num_k_tiles(i_tile, j_tile);
        for (size_type k_tile = 0; k_tile<n_k_tile; ++k_tile) {
          const size_type k_begin = tensor.k_begin(i_tile, j_tile, k_tile);
          //const size_type k_size  = tensor.k_size(i_tile, j_tile, k_tile);

          for (size_type i=0; i<i_size; ++i) {

            const size_type nEntry =
              tensor.num_entry(i_tile,j_tile,k_tile,i);
            const size_type iEntryBeg =
              tensor.entry_begin(i_tile,j_tile,k_tile,i);
            const size_type iEntryEnd = iEntryBeg + nEntry;
            size_type iEntry = iEntryBeg;

            VectorValue ytmp = 0 ;

            // Do entries with a blocked loop of size block_size
            if (block_size > 1) {
              const size_type nBlock = nEntry / block_size;
              const size_type nEntryB = nBlock * block_size;
              const size_type iEnd = iEntryBeg + nEntryB;

              TV vy;
              vy.zero();
              int j[block_size], k[block_size];

              for ( ; iEntry < iEnd ; iEntry += block_size ) {

                for (size_type ii=0; ii<block_size; ++ii) {
                  j[ii] = tensor.coord(iEntry+ii,0) + j_begin;
                  k[ii] = tensor.coord(iEntry+ii,1) + k_begin;
                }
                TV aj(a, j), ak(a, k), xj(x, j), xk(x, k),
                  c(&(tensor.value(iEntry)));

                // vy += c * ( aj * xk + ak * xj)
                aj.times_equal(xk);
                ak.times_equal(xj);
                aj.plus_equal(ak);
                c.times_equal(aj);
                vy.plus_equal(c);

              }

              ytmp += vy.sum();
            }

            // Do remaining entries with a scalar loop
            for ( ; iEntry<iEntryEnd; ++iEntry) {
              const size_type j = tensor.coord(iEntry,0) + j_begin;
              const size_type k = tensor.coord(iEntry,1) + k_begin;

              ytmp += tensor.value(iEntry) * ( a[j] * x[k] + a[k] * x[j] );
            }

            y[i+i_begin] += ytmp;

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

#endif /* #ifndef STOKHOS_THREADS_SIMPLE_TILED_CRS_PRODUCT_TENSOR_HPP */
