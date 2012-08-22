
#ifndef KOKKOSARRAY_IMPL_MULTIPLY_HPP
#define KOKKOSARRAY_IMPL_MULTIPLY_HPP

namespace KokkosArray {
namespace Impl {

template< class MatrixType ,
          class InputVectorType  = void ,
          class OutputVectorType = InputVectorType > class Multiply ;

template< class MatrixType ,
          class InputVectorType  = void ,
          class OutputVectorType = InputVectorType > class MMultiply ;

template < class ValueType, class Device > class MatrixMarketWriter ;

template < class VectorType > class Update ;

} // namespace Impl
} // namespace KokkosArray

#endif

