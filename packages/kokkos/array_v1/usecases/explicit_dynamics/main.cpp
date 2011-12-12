#include <iostream>
#include <cstdlib>

namespace Test{
  void test_Host(int beg, int end, int r, int threads);
  void test_Cuda(int beg, int end, int r);
}

int main(int argc, char ** argv)
{
  int beg = 4 ;
  int end = 12 ;
  int runs = 3 ;
  int host_threads = -1;

  if ( 1 < argc ) {
    host_threads = atoi(argv[1]);
  }
  if ( argc == 5) {
    beg = atoi(argv[2]);
    end = atoi(argv[3]);
    runs = atoi(argv[4]);
  }

std::cout << "\" " << argv[0]
                   << " host_threads begin end runs \"" << std::endl ;
std::cout << "\" " << argv[0]
          << " " << host_threads
          << " " << beg
          << " " << end
          << " " << runs
          << " \"" << std::endl ;


  if ( 0 <= host_threads ) {
    Test::test_Host(beg, end, runs, host_threads);
  }

#ifdef TEST_KOKKOS_CUDA
  Test::test_Cuda(beg , end, runs);
#endif

  return 0;
}
