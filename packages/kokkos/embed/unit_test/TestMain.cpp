
#include <iostream>
#include <Kokkos_Macros.hpp>

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
  std::cout << std::endl << "test< Host >()" << std::endl ;
  test< Kokkos::HostSpace >();

  std::cout << std::endl << "testdyn< Host >()" << std::endl ;
  testdyn< Kokkos::HostSpace >();

  std::cout << std::endl << "test_functor< Host >()" << std::endl ;
  test_functor< Kokkos::HostSpace >();

  std::cout << std::endl << "test_inner_product< Host >()" << std::endl ;
  test_inner_product< Kokkos::HostSpace >();

  std::cout << std::endl << "test< Cuda >()" << std::endl ;
  test< Kokkos::CudaSpace >();

  std::cout << std::endl << "testdyn< Cuda >()" << std::endl ;
  testdyn< Kokkos::CudaSpace >();

  std::cout << std::endl << "test_functor< Cuda >()" << std::endl ;
  test_functor< Kokkos::CudaSpace >();

  std::cout << std::endl << "test_inner_product< Cuda >()" << std::endl ;
  test_inner_product< Kokkos::CudaSpace >();

  return 0 ;
}

