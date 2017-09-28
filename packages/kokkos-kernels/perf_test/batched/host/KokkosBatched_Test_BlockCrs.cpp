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

#include "Kokkos_Core.hpp"
#include "impl/Kokkos_Timer.hpp"

#include "KokkosBatched_Util.hpp"

#if (0)
typedef KokkosBatched::Experimental::Algo::LU::Unblocked   AlgoLU;
typedef KokkosBatched::Experimental::Algo::Trsm::Unblocked AlgoTrsm;
typedef KokkosBatched::Experimental::Algo::Gemm::Unblocked AlgoGemm;

typedef KokkosBatched::Experimental::Algo::Trsv::Unblocked AlgoTrsv;
typedef KokkosBatched::Experimental::Algo::Gemv::Unblocked AlgoGemv;
#else
typedef KokkosBatched::Experimental::Algo::LU::Blocked   AlgoLU;
typedef KokkosBatched::Experimental::Algo::Trsm::Blocked AlgoTrsm;
typedef KokkosBatched::Experimental::Algo::Gemm::Blocked AlgoGemm;

typedef KokkosBatched::Experimental::Algo::Trsv::Blocked AlgoTrsv;
typedef KokkosBatched::Experimental::Algo::Gemv::Blocked AlgoGemv;
#endif

#include "KokkosBatched_Test_BlockCrs.hpp"

using namespace KokkosBatched;

int main (int argc, char *argv[]) {
  Kokkos::initialize(argc, argv); 

  typedef Kokkos::DefaultHostExecutionSpace HostSpace;
  typedef Kokkos::DefaultExecutionSpace DeviceSpace;

  const bool detail = false;
  std::cout << "DeviceSpace::  "; DeviceSpace::print_configuration(std::cout, detail);
  std::cout << "HostSpace::    ";   HostSpace::print_configuration(std::cout, detail);
  std::cout << std::endl;

  // dummy vector length (this is explicitly used in cuda)
  constexpr int VectorLength = 0;
  constexpr int RangeTagOper = 0;

  // Unit tests
  bool profile = false;
  for (int i=1;i<argc;++i) {
    const std::string& token = argv[i];
    if (strncmp(token.c_str(), "-profile", 8) == 0) profile = true;
  }
  
  if (!profile) {
    // including compact layer, it is not possible to test 
    // scalar and vector in the same code with out templating
    // std::cout << " Unit Test::Range::Scalar :: Begin\n";
    // {
    //   Test::run<DeviceSpace,Test::scalar_type,VectorLength,RangeTagOper>( 3,  4,  2, 25, 2);
    //   Test::run<DeviceSpace,Test::scalar_type,VectorLength,RangeTagOper>(44, 63, 15,  4, 1);
    //   Test::run<DeviceSpace,Test::scalar_type,VectorLength,RangeTagOper>( 2,  2, 15,  3, 3);
    //   Test::run<DeviceSpace,Test::scalar_type,VectorLength,RangeTagOper>( 1,  1,  2, 63, 8);
      
    //   for (int nrhs=1;nrhs<=33;++nrhs)
    //     Test::run<DeviceSpace,Test::scalar_type,VectorLength,RangeTagOper>(2, 2, 15, 3, nrhs);
    // }
    // std::cout << " Unit Test::Range::Scalar :: End\n";

    std::cout << " Unit Test::Range::Vector :: Begin\n";
    {
#if defined(__AVX512F__)
      //typedef VectorTag<SIMD<Test::scalar_type>,8> VectorType;
      typedef Experimental::Vector<VectorTag<AVX<Test::scalar_type>,8> > VectorType;
#elif defined(__AVX2__) || defined(__AVX__)
      //typedef VectorTag<SIMD<Test::scalar_type>,4> VectorType;
      typedef Experimental::Vector<VectorTag<AVX<Test::scalar_type>,4> > VectorType;
#endif
      Test::run<DeviceSpace,VectorType,VectorLength,RangeTagOper>( 3,  4,  2, 25, 2);
      Test::run<DeviceSpace,VectorType,VectorLength,RangeTagOper>(44, 63, 15,  4, 1);
      Test::run<DeviceSpace,VectorType,VectorLength,RangeTagOper>( 2,  2, 15,  3, 3);
      Test::run<DeviceSpace,VectorType,VectorLength,RangeTagOper>( 1,  1,  2, 63, 8);
      
      for (int nrhs=1;nrhs<=33;++nrhs)
        Test::run<DeviceSpace,VectorType,VectorLength,RangeTagOper>(2, 2, 15, 3, nrhs);
    }
    std::cout << " Unit Test::Range::Vector :: End\n";
  }

  // MKL
#if defined(__KOKKOSBATCHED_INTEL_MKL_COMPACT_BATCHED__)
  std::cout << " Perf Test::CompactMKL Begin\n";
  {
#if defined(__AVX512F__)
    //typedef VectorTag<SIMD<Test::scalar_type>,8> VectorType;
    typedef Experimental::Vector<VectorTag<AVX<Test::scalar_type>,8> > VectorType;
#elif defined(__AVX2__) || defined(__AVX__)
    //typedef VectorTag<SIMD<Test::scalar_type>,4> VectorType;
    typedef Experimental::Vector<VectorTag<AVX<Test::scalar_type>,4> > VectorType;
#endif
    const bool test_mkl = true;
    const Test::Input input(argc, argv); 
    int r_val = Test::run<DeviceSpace,VectorType,VectorLength>(input, test_mkl);
    r_val = 0;
  } 
  std::cout << " Perf Test::CompactMKL End\n";  
#endif

  // Performance tests
  std::cout << " Perf Test::Vector Begin\n";
  {
#if defined(__AVX512F__)
    //typedef VectorTag<SIMD<Test::scalar_type>,8> VectorType;
    typedef Experimental::Vector<VectorTag<AVX<Test::scalar_type>,8> > VectorType;
#elif defined(__AVX2__) || defined(__AVX__)
    //typedef VectorTag<SIMD<Test::scalar_type>,4> VectorType;
    typedef Experimental::Vector<VectorTag<AVX<Test::scalar_type>,4> > VectorType;
#endif

    const Test::Input input(argc, argv); 
    int r_val = Test::run<DeviceSpace,VectorType,VectorLength>(input);
    r_val = 0;
  } 
  std::cout << " Perf Test::Vector End\n";

  Kokkos::finalize();

  return 0;
}
