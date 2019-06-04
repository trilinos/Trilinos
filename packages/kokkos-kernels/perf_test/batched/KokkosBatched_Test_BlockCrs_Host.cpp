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

//#define KOKKOSBATCHED_USE_UNBLOCKED_ALGO 1
#define KOKKOSBATCHED_USE_BLOCKED_ALGO 1

#if defined (KOKKOSBATCHED_USE_UNBLOCKED_ALGO)
typedef KokkosBatched::Experimental::Algo::LU::Unblocked   AlgoLU;
typedef KokkosBatched::Experimental::Algo::Trsm::Unblocked AlgoTrsm;
typedef KokkosBatched::Experimental::Algo::Gemm::Unblocked AlgoGemm;

typedef KokkosBatched::Experimental::Algo::Trsv::Unblocked AlgoTrsv;
typedef KokkosBatched::Experimental::Algo::Gemv::Unblocked AlgoGemv;
#endif
#if defined (KOKKOSBATCHED_USE_BLOCKED_ALGO)
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

#if !defined(__CUDA_ARCH__)
  typedef Kokkos::DefaultHostExecutionSpace HostSpaceType;
  const bool detail = false;

  Kokkos::print_configuration(std::cout, detail);

  enum : int { VectorLength = DefaultVectorLength<Test::scalar_type,typename HostSpaceType::memory_space>::value,
               RangeTagOper = 0 };

  // vector type
  typedef Experimental::Vector<SIMD<Test::scalar_type>,VectorLength> VectorType;  

  // Unit tests
  bool profile = false;
  for (int i=1;i<argc;++i) {
    const std::string& token = argv[i];
    if (strncmp(token.c_str(), "-profile", 8) == 0) profile = true;
  }


  if (!profile) {
    // including compact layer, it is not possible to test 
    // scalar and vector in the same code without templating
    std::cout << " Unit Test::Range::Vector :: Begin\n";
    {
      Test::run<HostSpaceType,VectorType,VectorLength,RangeTagOper>( 3,  4,  2, 25, 2);
      Test::run<HostSpaceType,VectorType,VectorLength,RangeTagOper>(44, 63, 15,  4, 1);
      Test::run<HostSpaceType,VectorType,VectorLength,RangeTagOper>( 2,  2, 15,  3, 3);
      
      for (int nrhs=1;nrhs<=33;++nrhs) 
        Test::run<HostSpaceType,VectorType,VectorLength,RangeTagOper>(2, 2, 15, 3, nrhs);
    }

    std::cout << " Unit Test::Range::Vector :: End\n";
  }
  
  // MKL
#if defined(__KOKKOSBATCHED_INTEL_MKL_COMPACT_BATCHED__)
  std::cout << " Perf Test::CompactMKL Begin\n";
  {
    const bool test_mkl = true;
    const Test::Input<HostSpaceType> input(argc, argv); 
    Test::run<HostSpaceType,VectorType,VectorLength>(input, test_mkl);
  } 
  std::cout << " Perf Test::CompactMKL End\n";  
#endif

  // Performance tests
  std::cout << " Perf Test::Vector Begin\n";
  {
    const Test::Input<HostSpaceType> input(argc, argv); 
    Test::run<HostSpaceType,VectorType,VectorLength>(input);
  } 
  std::cout << " Perf Test::Vector End\n";

#endif
  Kokkos::finalize();

  return 0;
}
