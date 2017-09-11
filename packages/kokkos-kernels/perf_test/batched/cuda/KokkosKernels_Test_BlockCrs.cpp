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

#include "KokkosKernels_Util.hpp"

#if (1)
typedef KokkosKernels::Batched::Experimental::Algo::LU::Unblocked   AlgoLU;
typedef KokkosKernels::Batched::Experimental::Algo::Trsm::Unblocked AlgoTrsm;
typedef KokkosKernels::Batched::Experimental::Algo::Gemm::Unblocked AlgoGemm;

typedef KokkosKernels::Batched::Experimental::Algo::Trsv::Unblocked AlgoTrsv;
typedef KokkosKernels::Batched::Experimental::Algo::Gemv::Unblocked AlgoGemv;
#else
typedef KokkosKernels::Batched::Experimental::Algo::LU::Blocked   AlgoLU;
typedef KokkosKernels::Batched::Experimental::Algo::Trsm::Blocked AlgoTrsm;
typedef KokkosKernels::Batched::Experimental::Algo::Gemm::Blocked AlgoGemm;

typedef KokkosKernels::Batched::Experimental::Algo::Trsv::Blocked AlgoTrsv;
typedef KokkosKernels::Batched::Experimental::Algo::Gemv::Blocked AlgoGemv;
#endif

#include "KokkosKernels_Test_BlockCrs.hpp"

using namespace KokkosKernels;

int main (int argc, char *argv[]) {
  Kokkos::initialize(argc, argv); 

  typedef Kokkos::DefaultHostExecutionSpace HostSpace;
  typedef Kokkos::DefaultExecutionSpace DeviceSpace;

  const bool detail = false;
  std::cout << "DeviceSpace::  "; DeviceSpace::print_configuration(std::cout, detail);
  std::cout << "HostSpace::    ";   HostSpace::print_configuration(std::cout, detail);
  std::cout << std::endl;

  constexpr int VectorLength = 16;

  // Unit tests
  bool profile = false;
  for (int i=1;i<argc;++i) {
    const std::string& token = argv[i];
    if (strncmp(token.c_str(), "-profile", 8) == 0) profile = true;
  }
  
  if (!profile) {
    std::cout << " Unit Test::Range :: Begin\n";
    {
      Test::run<DeviceSpace,Test::scalar_type,VectorLength,0>( 3,  4,  2, 25, 2);
      Test::run<DeviceSpace,Test::scalar_type,VectorLength,0>(44, 63, 15,  4, 1);
      Test::run<DeviceSpace,Test::scalar_type,VectorLength,0>( 2,  2, 15,  3, 3);
      Test::run<DeviceSpace,Test::scalar_type,VectorLength,0>( 1,  1,  2, 63, 8);
      
      for (int nrhs=1;nrhs<=33;++nrhs)
        Test::run<DeviceSpace,Test::scalar_type,VectorLength,0>(2, 2, 15, 3, nrhs);
    }
    std::cout << " Unit Test::Range :: End\n";
    
    std::cout << " Unit Test::Team :: Begin\n";
    {
      Test::run<DeviceSpace,Test::scalar_type,VectorLength,1>( 3,  4,  2, 25, 2);
      Test::run<DeviceSpace,Test::scalar_type,VectorLength,1>(44, 63, 15,  4, 1);
      Test::run<DeviceSpace,Test::scalar_type,VectorLength,1>( 2,  2, 15,  3, 3);
      Test::run<DeviceSpace,Test::scalar_type,VectorLength,1>( 1,  1,  2, 63, 8);
      
      for (int nrhs=1;nrhs<=33;++nrhs)
        Test::run<DeviceSpace,Test::scalar_type,VectorLength,1>(2, 2, 15, 3, nrhs);
    }
    std::cout << " Unit Test::Team :: End\n";
    
    // std::cout << " Unit Test::TeamShmem :: Begin\n";
    // {
    //   //Test::run<DeviceSpace,Test::scalar_type,VectorLength,2>( 3,  4,  2, 25, 2);
    //   Test::run<DeviceSpace,Test::scalar_type,VectorLength,2>(44, 63, 15,  4, 1);
    //   Test::run<DeviceSpace,Test::scalar_type,VectorLength,2>( 2,  2, 15,  3, 3);
    //   //Test::run<DeviceSpace,Test::scalar_type,VectorLength,2>( 1,  1,  2, 63, 8);
      
    //   for (int nrhs=1;nrhs<=4;++nrhs)
    //   Test::run<DeviceSpace,Test::scalar_type,VectorLength,2>(2, 2, 15, 3, nrhs);
    // }
    // std::cout << " Unit Test::TeamShmem :: End\n";
  }

  // Performance tests
  std::cout << " Perf Test:: Begin\n";
  {
    const Test::Input input(argc, argv); 
    int r_val = Test::run<DeviceSpace,Test::scalar_type,VectorLength>(input);
    r_val = 0;
  } 
  std::cout << " Perf Test:: End\n";

  Kokkos::finalize();

  return 0;
}
