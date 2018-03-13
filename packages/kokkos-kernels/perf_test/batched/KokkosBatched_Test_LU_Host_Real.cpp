/// \author Kyungjoo Kim (kyukim@sandia.gov)

#define KokkosBatched_Test_LU_Host_Real 1
#include "KokkosBatched_Test_LU_Host.hpp"

using namespace KokkosBatched::Experimental;

template<typename AlgoTagType>
void run(const int N) {
  typedef Kokkos::DefaultHostExecutionSpace HostSpaceType;

  Kokkos::print_configuration(std::cout, false);

  PerfTest::LU< 3, HostSpaceType,AlgoTagType>(N);
  PerfTest::LU< 5, HostSpaceType,AlgoTagType>(N);
  PerfTest::LU<10, HostSpaceType,AlgoTagType>(N);
  PerfTest::LU<15, HostSpaceType,AlgoTagType>(N);
}

int main(int argc, char *argv[]) {

  Kokkos::initialize(argc, argv);
  
  int N = 128*128;

  for (int i=1;i<argc;++i) {
    const std::string& token = argv[i];
    if (token == std::string("-N")) N = std::atoi(argv[++i]);
  }

  {
    std::cout << " N = " << N << std::endl;
    
    std::cout << "\n Testing Algo::LU::Unblocked\n";
    run<Algo::LU::Unblocked>(N);
    
    std::cout << "\n Testing Algo::LU::Blocked\n";
    run<Algo::LU::Blocked>(N);
  }

  Kokkos::finalize();

  return 0;
}
