
#ifndef KOKKOS_IMPL_MULTIPLY_HPP
#define KOKKOS_IMPL_MULTIPLY_HPP

namespace Kokkos {
namespace Impl {

template< class MatrixType ,
          class InputVectorType  = void ,
          class OutputVectorType = InputVectorType > class Multiply ;

template< class MatrixType ,
          class InputVectorType  = void ,
          class OutputVectorType = InputVectorType > class MMultiply ;

template < class ValueType, class Device > class MatrixMarketWriter ;

template< typename ValueType, typename VectorValue >
class Update
{
public:
  typedef VectorValue                               vector_type;
  typedef ValueType                                 value_type ;
  typedef typename vector_type::execution_space         execution_space ;
  typedef typename execution_space::size_type           size_type ;
  

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

  KOKKOS_INLINE_FUNCTION
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
    parallel_for( row_count , Update(alpha,x,beta,y) );
  }
};

} // namespace Impl
} // namespace Kokkos

#endif

