/* Implementation for testing KokkosKernels on BCRS operations
   - block-tridiagonal factorization
   - block-tridiagonal solve
   - bcrs matvec

   StructuredBlock represents a 3D mesh having ni, nj, nk cells in each
   dimension. Variable ordering is such that the k index is the fastest and the
   i index is slowest. Smoothing lines are built in the k direction.
   BlockCrsMatrix is a simple block CRS data structure.
   BlockTridiagMatrices holds the block tridiagonal matrices.

   An example run is
   ./driver -ni 32 -nj 32 -nk 128 -bs 5 -c

   This runs a sequence of unit tests, then runs a problem having a 32x32x128
   structured block with the lines oriented along the third dimension (line
   length = 128). The block size is 5. -c adds a somewhat expensive check of the
   answer. It's good to run with -c once in a while, but the cheap unit tests
   that always run before the big problem already provide good coverage.
*/

#include "KokkosKernels_Test_BlockCrs.hpp"

using namespace KokkosKernels;

int main (int argc, char *argv[]) {
  Kokkos::initialize(argc, argv); 

  typedef Kokkos::DefaultHostExecutionSpace HostSpace;
  typedef Kokkos::DefaultHostExecutionSpace DeviceSpace;

  const bool detail = false;
  std::cout << "DeviceSpace::  "; DeviceSpace::print_configuration(std::cout, detail);
  std::cout << "HostSpace::    ";   HostSpace::print_configuration(std::cout, detail);
  std::cout << std::endl;

  // Unit tests
  // std::cout << " Unit Test::Flat :: Begin\n";
  // {
  //   Test::run<DeviceSpace>( 3,  4,  2, 25, 2);
  //   Test::run<DeviceSpace>(44, 63, 15,  4, 1);
  //   Test::run<DeviceSpace>( 2,  2, 15,  3, 3);
  //   Test::run<DeviceSpace>( 1,  1,  2, 63, 8);

  //   for (int nrhs=1;nrhs<=33;++nrhs)
  //     Test::run<DeviceSpace>(2, 2, 15, 3, nrhs);
  // }
  // std::cout << " Unit Test::Flat :: End\n";

  std::cout << " Unit Test::Vector :: Begin\n";
  {
#if defined(__AVX512F__)  
    //typedef VectorTag<SIMD<Test::scalar_type>,8> VectorType;
    typedef Vector<VectorTag<AVX<Test::scalar_type>,8> > VectorType;
#elif defined(__AVX2__) || defined(__AVX__)
    //typedef VectorTag<SIMD<Test::scalar_type>,4> VectorType;
    typedef Vector<VectorTag<AVX<Test::scalar_type>,4> > VectorType;
#endif    

    Test::run<DeviceSpace,VectorType>( 3,  4,  2, 25, 2);
    Test::run<DeviceSpace,VectorType>(44, 63, 15,  4, 1);
    Test::run<DeviceSpace,VectorType>( 2,  2, 15,  3, 3);
    Test::run<DeviceSpace,VectorType>( 1,  1,  2, 63, 8);

    for (int nrhs=1;nrhs<=33;++nrhs)
      Test::run<DeviceSpace,VectorType>(2, 2, 15, 3, nrhs);
  }
  std::cout << " Unit Test::Vector :: End\n";

#if defined(__KOKKOSKERNELS_INTEL_MKL_COMPACT_BATCHED__)
  std::cout << " Unit Test::MKLCmpt :: Begin\n";
  {
#if defined(__AVX512F__)  
    //typedef VectorTag<SIMD<Test::scalar_type>,8> VectorType;
    typedef Vector<VectorTag<AVX<Test::scalar_type>,8> > VectorType;
#elif defined(__AVX2__) || defined(__AVX__)
    //typedef VectorTag<SIMD<Test::scalar_type>,4> VectorType;
    typedef Vector<VectorTag<AVX<Test::scalar_type>,4> > VectorType;
#endif    

    const int test_mkl = 1;

    Test::run<DeviceSpace,VectorType>( 3,  4,  2, 25, 2, test_mkl);
    Test::run<DeviceSpace,VectorType>(44, 63, 15,  4, 1, test_mkl);
    Test::run<DeviceSpace,VectorType>( 2,  2, 15,  3, 3, test_mkl);
    Test::run<DeviceSpace,VectorType>( 1,  1,  2, 63, 8, test_mkl);

    for (int nrhs=1;nrhs<=33;++nrhs)
      Test::run<DeviceSpace,VectorType>(2, 2, 15, 3, nrhs, test_mkl);
  }
  std::cout << " Unit Test::MKLCmpt :: End\n";

  std::cout << " Unit Test::MKLBatch :: Begin\n";
  {
#if defined(__AVX512F__)  
    //typedef VectorTag<SIMD<Test::scalar_type>,8> VectorType;
    typedef Vector<VectorTag<AVX<Test::scalar_type>,8> > VectorType;
#elif defined(__AVX2__) || defined(__AVX__)
    //typedef VectorTag<SIMD<Test::scalar_type>,4> VectorType;
    typedef Vector<VectorTag<AVX<Test::scalar_type>,4> > VectorType;
#endif    

    const int test_mkl = 2;

    Test::run<DeviceSpace,VectorType>( 3,  4,  2, 25, 2, test_mkl);
    Test::run<DeviceSpace,VectorType>(44, 63, 15,  4, 1, test_mkl);
    Test::run<DeviceSpace,VectorType>( 2,  2, 15,  3, 3, test_mkl);
    Test::run<DeviceSpace,VectorType>( 1,  1,  2, 63, 8, test_mkl);

    for (int nrhs=1;nrhs<=33;++nrhs)
      Test::run<DeviceSpace,VectorType>(2, 2, 15, 3, nrhs, test_mkl);
  }
  std::cout << " Unit Test::MKLBatch:: End\n";
#endif

  // // Performance tests
  // std::cout << " Perf Test::Flat :: Begin\n";
  // {
  //   const Test::Input input(argc, argv);
  //   int r_val = Test::run<DeviceSpace>(input);
  //   r_val = 0;
  // } 
  // std::cout << " Perf Test::Flat :: End\n";

  std::cout << " Perf Test::Vector :: Begin\n";
  {
#if defined(__AVX512F__)  
    //typedef VectorTag<SIMD<Test::scalar_type>,8> VectorType;
    typedef Vector<VectorTag<AVX<Test::scalar_type>,8> > VectorType;
#elif defined(__AVX2__) || defined(__AVX__)
    //typedef VectorTag<SIMD<Test::scalar_type>,4> VectorType;
    typedef Vector<VectorTag<AVX<Test::scalar_type>,4> > VectorType;
#endif    

    const Test::Input input(argc, argv);
    int r_val = Test::run<DeviceSpace,VectorType>(input);
    r_val = 0;
  } 
  std::cout << " Perf Test::Vector :: End\n";

#if defined(__KOKKOSKERNELS_INTEL_MKL_COMPACT_BATCHED__)
  std::cout << " Perf Test::MKL Cmpt:: Begin\n";
  {
#if defined(__AVX512F__)  
    //typedef VectorTag<SIMD<Test::scalar_type>,8> VectorType;
    typedef Vector<VectorTag<AVX<Test::scalar_type>,8> > VectorType;
#elif defined(__AVX2__) || defined(__AVX__)
    //typedef VectorTag<SIMD<Test::scalar_type>,4> VectorType;
    typedef Vector<VectorTag<AVX<Test::scalar_type>,4> > VectorType;
#endif    

    const Test::Input input(argc, argv);
    const int test_mkl = 1;
    int r_val = Test::run<DeviceSpace,VectorType>(input, test_mkl);
    r_val = 0;
  } 
  std::cout << " Perf Test::MKL Cmpt:: End\n";

  std::cout << " Perf Test::MKL Batch:: Begin\n";
  {
#if defined(__AVX512F__)  
    //typedef VectorTag<SIMD<Test::scalar_type>,8> VectorType;
    typedef Vector<VectorTag<AVX<Test::scalar_type>,8> > VectorType;
#elif defined(__AVX2__) || defined(__AVX__)
    //typedef VectorTag<SIMD<Test::scalar_type>,4> VectorType;
    typedef Vector<VectorTag<AVX<Test::scalar_type>,4> > VectorType;
#endif    

    const Test::Input input(argc, argv);
    const int test_mkl = 2;
    int r_val = Test::run<DeviceSpace,VectorType>(input, test_mkl);
    r_val = 0;
  } 
  std::cout << " Perf Test::MKL Batch:: End\n";
#endif

  Kokkos::finalize();

  return 0;
}
