#include <iostream>
#include <cstdlib>

namespace test{
  void test_Host(int beg, int end, int r);
  void test_Pthread (int beg, int end, int r, int t);
  void test_TPI (int beg, int end, int r, int t);
  void test_TBB(int beg, int end, int r, int t);
  void test_Cuda(int beg, int end, int r);
  void test_NUMA(int beg, int end, int r);
}

int main(int argc, char ** argv)
{
  int beg = 4 ;
  int end = 12 ;
  int runs = 3 ;
  int threads = 32;

  if ( argc == 5) {
    beg = atoi(argv[1]);
    end = atoi(argv[2]);
    runs = atoi(argv[3]);
    threads = atoi(argv[4]);
  }

#ifdef TEST_KOKKOS_HOST
  test::test_Host(beg, end, runs);
#endif
#ifdef TEST_KOKKOS_PTHREAD
  test::test_Pthread (beg, end, runs, threads);
#endif
#ifdef TEST_KOKKOS_NUMA
  test::test_NUMA (beg, end, runs);
#endif
#ifdef TEST_KOKKOS_TPI
  test::test_TPI (beg, end, runs, threads);
#endif
#ifdef TEST_KOKKOS_TBB
  test::test_TBB (beg, end, runs);
#endif
#ifdef TEST_KOKKOS_CUDA
  test::test_Cuda(beg , end, runs);
#endif

  return 0;
}
