
#include <iostream>
#include <KokkosArray_Macros.hpp>

template< class Device > int test();
template< class Device > int test_functor();

template<> int test< KokkosArray::Host >();
template<> int test< KokkosArray::Cuda >();
template<> int test_functor< KokkosArray::Host >();
template<> int test_functor< KokkosArray::Cuda >();

#if defined( HAVE_CUDA )
#include <KokkosArray_Cuda.hpp>

int test_cuda()
{
  KokkosArray::Cuda::initialize( KokkosArray::Cuda::SelectDevice(0) );

  std::cout << std::endl << "test< KokkosArray::Cuda >()" << std::endl ;
  test< KokkosArray::Cuda >();

  std::cout << std::endl << "test_functor< KokkosArray::Cuda >()" << std::endl ;
  test_functor< KokkosArray::Cuda >();

  KokkosArray::Cuda::finalize();

  return 0 ;
}

#else
int test_cuda() { return 0 ; }
#endif

int main()
{
  std::cout << std::endl << "test< KokkosArray::Host >()" << std::endl ;
  test< KokkosArray::Host >();

  std::cout << std::endl << "test_functor< KokkosArray::Host >()" << std::endl ;
  test_functor< KokkosArray::Host >();

  test_cuda();

  return 0 ;
}

