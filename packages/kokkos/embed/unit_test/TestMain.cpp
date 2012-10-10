
#include <iostream>
#include <KokkosArray_Macros.hpp>

template< class Device > int test();

template<> int test< KokkosArray::Host >();
template<> int test< KokkosArray::Cuda >();

int main()
{
  std::cout << std::endl << "test< KokkosArray::Host >()" << std::endl ;
  test< KokkosArray::Host >();

  std::cout << std::endl << "test< KokkosArray::Cuda >()" << std::endl ;
  test< KokkosArray::Cuda >();
  return 0 ;
}

