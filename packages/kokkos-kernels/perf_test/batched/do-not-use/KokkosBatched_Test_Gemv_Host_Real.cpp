/// \author Kyungjoo Kim (kyukim@sandia.gov)

#define KokkosBatched_Test_Gemv_Host_Real 1
#include "KokkosBatched_Test_Gemv_Host.hpp"

using namespace KokkosBatched;

template<typename AlgoTagType>
void run(const int N) {
  typedef Kokkos::DefaultHostExecutionSpace HostSpaceType;

  Kokkos::print_configuration(std::cout);

  // PerfTest::Gemv< 4, 1, ExecSpace,AlgoTagType>(N);
  // PerfTest::Gemv< 8, 1, ExecSpace,AlgoTagType>(N);
  // PerfTest::Gemv<16, 1, ExecSpace,AlgoTagType>(N);
  // PerfTest::Gemv<20, 1, ExecSpace,AlgoTagType>(N);
  // PerfTest::Gemv<32, 1, ExecSpace,AlgoTagType>(N);
  // PerfTest::Gemv<64, 1, ExecSpace,AlgoTagType>(N);

  PerfTest::Gemv< 3, 1, HostSpaceType,AlgoTagType>(N);
  PerfTest::Gemv< 5, 1, HostSpaceType,AlgoTagType>(N);
  PerfTest::Gemv<10, 1, HostSpaceType,AlgoTagType>(N);
  PerfTest::Gemv<15, 1, HostSpaceType,AlgoTagType>(N);
}

int main(int argc, char *argv[]) {
  
  Kokkos::initialize(argc, argv);
#if !defined(__CUDA_ARCH__)
  const int ntest = 1;
  //const int N[6] = { 256, 512, 768, 1024, 1280, 1536 };
  const int N[1] = { 128*128 };

  {        
    for (int i=0;i<ntest;++i) {
      std::cout << " N = " << N[i] << std::endl;
      
      std::cout << "\n Testing Algo::Gemv::Unblocked\n";
      run<Algo::Gemv::Unblocked>(N[i]);
      
      std::cout << "\n Testing Algo::Gemv::Blocked\n";
      run<Algo::Gemv::Blocked>(N[i]);
    }
  }
#endif
  Kokkos::finalize();
  
  return 0;
}
