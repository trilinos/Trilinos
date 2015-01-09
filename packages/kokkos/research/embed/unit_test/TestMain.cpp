
#include <iostream>
#include <Kokkos_Core.h>

template< class Space > int test();
template< class Space > int testdyn();
template< class Space > int test_functor();
template< class Space > int test_inner_product();

template<> int test< Kokkos::HostSpace >();
template<> int test< Kokkos::CudaSpace >();
template<> int testdyn< Kokkos::HostSpace >();
template<> int testdyn< Kokkos::CudaSpace >();
template<> int test_functor< Kokkos::HostSpace >();
template<> int test_functor< Kokkos::CudaSpace >();
template<> int test_inner_product< Kokkos::HostSpace >();
template<> int test_inner_product< Kokkos::CudaSpace >();

int test_cuda();

int main()
{
  std::cout << std::endl << "test< HostSpace >()" << std::endl ;
  test< Kokkos::HostSpace >();

  std::cout << std::endl << "testdyn< HostSpace >()" << std::endl ;
  testdyn< Kokkos::HostSpace >();

  std::cout << std::endl << "test_functor< HostSpace >()" << std::endl ;
  test_functor< Kokkos::HostSpace >();

  std::cout << std::endl << "test_inner_product< HostSpace >()" << std::endl ;
  test_inner_product< Kokkos::HostSpace >();

#if defined( KOKKOS_HAVE_CUDA )

  std::cout << std::endl << "test< CudaSpace >()" << std::endl ;
  test< Kokkos::CudaSpace >();

  std::cout << std::endl << "testdyn< CudaSpace >()" << std::endl ;
  testdyn< Kokkos::CudaSpace >();

  std::cout << std::endl << "test_functor< CudaSpace >()" << std::endl ;
  test_functor< Kokkos::CudaSpace >();

  std::cout << std::endl << "test_inner_product< CudaSpace >()" << std::endl ;
  test_inner_product< Kokkos::CudaSpace >();

#endif

  return 0 ;
}

