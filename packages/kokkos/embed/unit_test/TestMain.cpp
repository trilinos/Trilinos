
#include <iostream>
#include <KokkosArray_Macros.hpp>

template< class Device > int test();
template< class Device > int test_functor();

template<> int test< KokkosArray::HostSpace >();
template<> int test< KokkosArray::CudaSpace >();
template<> int test_functor< KokkosArray::HostSpace >();
template<> int test_functor< KokkosArray::CudaSpace >();

int test_cuda();

int main()
{
  std::cout << std::endl << "test< Host >()" << std::endl ;
  test< KokkosArray::HostSpace >();

  std::cout << std::endl << "test_functor< Host >()" << std::endl ;
  test_functor< KokkosArray::HostSpace >();

  std::cout << std::endl << "test< Cuda >()" << std::endl ;
  test< KokkosArray::CudaSpace >();

  std::cout << std::endl << "test_functor< Cuda >()" << std::endl ;
  test_functor< KokkosArray::CudaSpace >();

  return 0 ;
}

