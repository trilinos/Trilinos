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

#ifndef STOKHOS_BLOCKCRSMATRIX_HPP
#define STOKHOS_BLOCKCRSMATRIX_HPP

#include "KokkosArray_CrsArray.hpp"
#include "KokkosArray_View.hpp"

#include "Stokhos_Multiply.hpp"

namespace Stokhos {

/** \brief  CRS matrix of dense blocks.
 *
 *  Matrix coefficients are stored by block and then by Crs entry.
 *    m_values( block.size() , m_graph.entry_count() )
 *
 *  Vectors are conformally stored as
 *    View( block.dimension() , m_graph.row_map.length() )
 */
template< class BlockSpec , typename ValueType , class Device >
class BlockCrsMatrix {
public:
 
  typedef Device                              device_type ;
  typedef typename device_type::size_type     size_type ;
  typedef ValueType                           value_type ;
  typedef BlockSpec                           block_spec ;
  typedef KokkosArray::CrsArray< size_type , device_type > graph_type ;
  typedef KokkosArray::View< value_type**, KokkosArray::LayoutLeft, device_type >  block_vector_type ;

  block_vector_type  values ;
  graph_type         graph ;
  block_spec         block ;
};

template< class BlockSpec ,
          typename MatrixValueType ,
          typename VectorValueType ,
          class Device >
void multiply(const BlockCrsMatrix<BlockSpec,MatrixValueType,Device> & A ,
	      const KokkosArray::View<VectorValueType**,KokkosArray::LayoutLeft,Device> & x ,
	      const KokkosArray::View<VectorValueType**,KokkosArray::LayoutLeft,Device> & y )
{
  typedef BlockCrsMatrix<BlockSpec,MatrixValueType,Device>  matrix_type ;
  typedef KokkosArray::View<VectorValueType**,KokkosArray::LayoutLeft,Device> block_vector_type ;

  Multiply<matrix_type,block_vector_type,block_vector_type,DefaultSparseMatOps>::apply( A , x , y );
}

template< class BlockSpec ,
          typename MatrixValueType ,
          typename VectorValueType ,
          class Device ,
	  class SparseMatOps >
void multiply(const BlockCrsMatrix<BlockSpec,MatrixValueType,Device> & A ,
	      const KokkosArray::View<VectorValueType**,KokkosArray::LayoutLeft,Device> & x ,
	      const KokkosArray::View<VectorValueType**,KokkosArray::LayoutLeft,Device> & y ,
	      const SparseMatOps& smo = SparseMatOps() )
{
  typedef BlockCrsMatrix<BlockSpec,MatrixValueType,Device>  matrix_type ;
  typedef KokkosArray::View<VectorValueType**,KokkosArray::LayoutLeft,Device> block_vector_type ;

  Multiply<matrix_type,block_vector_type,block_vector_type,SparseMatOps>::apply( A , x , y );
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Stokhos

#endif /* #ifndef STOKHOS_BLOCKCRSMATRIX_HPP */

