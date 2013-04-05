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

#ifndef STOKHOS_MULTIPLY_HPP
#define STOKHOS_MULTIPLY_HPP

namespace Stokhos {

class DefaultSparseMatOps {};

template< class MatrixType ,
          class InputVectorType  = void ,
          class OutputVectorType = InputVectorType ,
	  class SparseMatOps = DefaultSparseMatOps > class Multiply ;

template< class MatrixType ,
          class InputVectorType  = void ,
          class OutputVectorType = InputVectorType ,
	  class SparseMatOps = DefaultSparseMatOps > class MMultiply ;

template < class ValueType, class Device > class MatrixMarketWriter ;

template< typename ValueType, typename VectorValue >
class Update
{
public:
  typedef VectorValue                               vector_type;
  typedef ValueType                                 value_type ;
  typedef typename vector_type::device_type         device_type ;
  typedef typename device_type::size_type           size_type ;
  

  const vector_type  m_x ;
  const vector_type  m_y ;
  const value_type   m_alpha ;
  const value_type   m_beta ;

  Update( const value_type&   alpha ,
	  const vector_type & x ,
	  const value_type &  beta ,
	  const vector_type & y )
  : m_x( x )
  , m_y( y )
  , m_alpha( alpha )
  , m_beta( beta )
  {}

  //--------------------------------------------------------------------------

  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const size_type iRow ) const
  {
    m_x(iRow) = m_alpha * m_x(iRow) + m_beta * m_y(iRow) ;
  }

  static void apply( const value_type&   alpha ,
		     const vector_type & x ,
		     const value_type &  beta ,
                     const vector_type & y )
  {
    const size_t row_count = x.dimension_0() ;
    KokkosArray::parallel_for( row_count , Update(alpha,x,beta,y) );
  }
};

} // namespace Stokhos

#endif

