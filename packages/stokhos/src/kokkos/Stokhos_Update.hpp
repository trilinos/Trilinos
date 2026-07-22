// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_UPDATE_HPP
#define STOKHOS_UPDATE_HPP

namespace Stokhos {

template< typename ValueType, typename VectorType >
class Update
{
public:
  typedef VectorType                        vector_type;
  typedef ValueType                         value_type;
  typedef typename vector_type::execution_space execution_space;
  typedef typename execution_space::size_type   size_type;


        vector_type  m_x;
  const vector_type  m_y;
  const value_type   m_alpha;
  const value_type   m_beta;

  Update(const value_type& alpha, vector_type& x,
         const value_type& beta,  const vector_type& y)
  : m_x( x )
  , m_y( y )
  , m_alpha( alpha )
  , m_beta( beta )
  {}

  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type iRow ) const
  {
    m_x(iRow) = m_alpha * m_x(iRow) + m_beta * m_y(iRow) ;
  }

  static void apply(const value_type& alpha, vector_type& x,
                    const value_type& beta,  const vector_type& y)
  {
    const size_t row_count = x.extent(0);
    Kokkos::parallel_for( row_count , Update(alpha,x,beta,y) );
  }
};

template <typename ValueType, typename VectorType>
void update(const ValueType& alpha, VectorType& x,
            const ValueType& beta,  const VectorType& y)
{
  Update<ValueType,VectorType>::apply( alpha , x , beta, y );
}

} // namespace Stokhos

#endif
