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

  Multiply<matrix_type,block_vector_type,block_vector_type>::apply( A , x , y );
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Stokhos

#endif /* #ifndef STOKHOS_BLOCKCRSMATRIX_HPP */

