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

#if defined(KOKKOS_ENABLE_CUDA) 
#define __KOKKOSBATCHED_TEST_ENABLE_CUDA__

#include "KokkosBatched_Util.hpp"

#define KOKKOSBATCHED_USE_UNBLOCKED_ALGO 1
//#define KOKKOSBATCHED_USE_BLOCKED_ALGO 1

#if defined (KOKKOSBATCHED_USE_UNBLOCKED_ALGO)
typedef KokkosBatched::Algo::LU::Unblocked   AlgoLU;
typedef KokkosBatched::Algo::Trsm::Unblocked AlgoTrsm;
typedef KokkosBatched::Algo::Gemm::Unblocked AlgoGemm;

typedef KokkosBatched::Algo::Trsv::Unblocked AlgoTrsv;
typedef KokkosBatched::Algo::Gemv::Unblocked AlgoGemv;
#endif
#if defined (KOKKOSBATCHED_USE_BLOCKED_ALGO)
typedef KokkosBatched::Algo::LU::Blocked   AlgoLU;
typedef KokkosBatched::Algo::Trsm::Blocked AlgoTrsm;
typedef KokkosBatched::Algo::Gemm::Blocked AlgoGemm;

typedef KokkosBatched::Algo::Trsv::Blocked AlgoTrsv;
typedef KokkosBatched::Algo::Gemv::Blocked AlgoGemv;
#endif

#include "KokkosBatched_Test_BlockCrs.hpp"

using namespace KokkosBatched;

int main (int argc, char *argv[]) {
  Kokkos::initialize(argc, argv); 

  typedef Kokkos::DefaultExecutionSpace DeviceSpaceType;

  const bool detail = false;

  Kokkos::print_configuration(std::cout, detail);

  enum : int { VectorLength = DefaultVectorLength<Test::scalar_type,typename DeviceSpaceType::memory_space>::value,
               RangeTagOper = 0,
               TeamTagOper = 1 };
  
  // Unit tests
  bool profile = false;
  for (int i=1;i<argc;++i) {
    const std::string& token = argv[i];
    if (strncmp(token.c_str(), "-profile", 8) == 0) profile = true;
  }
  

  if (!profile) {
    // std::cout << " Unit Test::Range :: Begin\n";
    // {
    //   Test::run<DeviceSpaceType,Test::scalar_type,VectorLength,RangeTagOper>( 3,  4,  2, 25, 2);
    //   Test::run<DeviceSpaceType,Test::scalar_type,VectorLength,RangeTagOper>(44, 63, 15,  4, 1);
    //   Test::run<DeviceSpaceType,Test::scalar_type,VectorLength,RangeTagOper>( 2,  2, 15,  3, 3);
    //   Test::run<DeviceSpaceType,Test::scalar_type,VectorLength,RangeTagOper>( 1,  1,  2, 63, 8);
      
    //   for (int nrhs=1;nrhs<=33;++nrhs)
    //     Test::run<DeviceSpaceType,Test::scalar_type,VectorLength,RangeTagOper>(2, 2, 15, 3, nrhs);
    // }
    // std::cout << " Unit Test::Range :: End\n";
    
    std::cout << " Unit Test::Team :: Begin\n";
    {
      Test::run<DeviceSpaceType,Test::scalar_type,VectorLength,TeamTagOper>( 3,  4,  2, 25, 2);
      Test::run<DeviceSpaceType,Test::scalar_type,VectorLength,TeamTagOper>(44, 63, 15,  4, 1);
      Test::run<DeviceSpaceType,Test::scalar_type,VectorLength,TeamTagOper>( 2,  2, 15,  3, 3);
      Test::run<DeviceSpaceType,Test::scalar_type,VectorLength,TeamTagOper>( 1,  1,  2, 63, 8);
      
      for (int nrhs=1;nrhs<=33;++nrhs)
        Test::run<DeviceSpaceType,Test::scalar_type,VectorLength,TeamTagOper>(2, 2, 15, 3, nrhs);
    }
    std::cout << " Unit Test::Team :: End\n";
  }

  // Performance tests
  std::cout << " Perf Test:: Begin\n";
  {
    const Test::Input<DeviceSpaceType> input(argc, argv); 
    Test::run<DeviceSpaceType,Test::scalar_type,VectorLength>(input);
  } 
  std::cout << " Perf Test:: End\n";

  Kokkos::finalize();

  return 0;
}
#else

int main(int argc, char *argv[]) {
  std::cout << "Kokkos::Cuda is not enabled\n";
  return -1;
}

#endif
