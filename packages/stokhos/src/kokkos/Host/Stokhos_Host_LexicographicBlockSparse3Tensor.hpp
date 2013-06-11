// @HEADER
// ***********************************************************************
//
//                     Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov)
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
