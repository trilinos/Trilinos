/// \author Kyungjoo Kim (kyukim@sandia.gov)

#define KokkosBatched_Test_Gemm_Host_Real 1
#include "KokkosBatched_Test_Gemm_Host.hpp"

using namespace KokkosBatched;

template<typename AlgoTagType>
void run(const int N) {
  typedef Kokkos::DefaultHostExecutionSpace HostSpaceType;

  Kokkos::print_configuration(std::cout);

  // Test::Gemm< 4, AlgoTagType>(N);
  // Test::Gemm< 8, AlgoTagType>(N);
  // Test::Gemm<16, AlgoTagType>(N);
  // Test::Gemm<20, AlgoTagType>(N);
  // Test::Gemm<32, AlgoTagType>(N);
  // Test::Gemm<64, AlgoTagType>(N);

  PerfTest::Gemm< 3, HostSpaceType, AlgoTagType>(N);
  PerfTest::Gemm< 5, HostSpaceType, AlgoTagType>(N);
  PerfTest::Gemm<10, HostSpaceType, AlgoTagType>(N);
  PerfTest::Gemm<15, HostSpaceType, AlgoTagType>(N);
}

int main(int argc, char *argv[]) {
  
  Kokkos::initialize(argc, argv);
#if !defined(__CUDA_ARCH__)
  const int ntest = 1;
  //const int N[6] = { 256, 512, 768, 1024, 1280, 1536 };
  int N[1] = { 128*128 };

  for (int i=1;i<argc;++i) {
    const std::string& token = argv[i];
    if (token == std::string("-N")) N[0] = std::atoi(argv[++i]);
  }

  {        
    for (int i=0;i<ntest;++i) {
      std::cout << " N = " << N[i] << std::endl;
      
      std::cout << "\n Testing Algo::Gemm::Unblocked\n";
      run<Algo::Gemm::Unblocked>(N[i]);
      
      std::cout << "\n Testing Algo::Gemm::Blocked\n";
      run<Algo::Gemm::Blocked>(N[i]);

#if defined(__KOKKOSBATCHED_INTEL_MKL_COMPACT_BATCHED__)
      std::cout << "\n Testing Algo::Gemm::CompactMKL\n";
      run<Algo::Gemm::CompactMKL>(N[i]);
#endif

    }
  }

#endif
  Kokkos::finalize();

  return 0;
}
