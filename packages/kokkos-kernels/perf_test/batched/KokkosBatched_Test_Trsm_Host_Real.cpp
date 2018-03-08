/// \author Kyungjoo Kim (kyukim@sandia.gov)

#define KokkosBatched_Test_Trsm_Host_Real 1
#include "KokkosBatched_Test_Trsm_Host.hpp"

using namespace KokkosBatched::Experimental;

template<typename AlgoTagType>
void run(const int N) {
  typedef Kokkos::DefaultHostExecutionSpace HostSpaceType;

  Kokkos::print_configuration(std::cout, false);

  std::cout << "\n\n Used for Factorization \n\n";

  /// Left, Lower, NoTrans, UnitDiag (used in LU factorization and LU solve)

  PerfTest::Trsm<0, 3, 3, HostSpaceType,AlgoTagType>(N);
  PerfTest::Trsm<0, 5, 5, HostSpaceType,AlgoTagType>(N);
  PerfTest::Trsm<0,10,10, HostSpaceType,AlgoTagType>(N);
  PerfTest::Trsm<0,15,15, HostSpaceType,AlgoTagType>(N);

  /// Left, Lower, NoTrans, NonUnitDiag

  PerfTest::Trsm<1, 3, 3, HostSpaceType,AlgoTagType>(N);
  PerfTest::Trsm<1, 5, 5, HostSpaceType,AlgoTagType>(N);
  PerfTest::Trsm<1,10,10, HostSpaceType,AlgoTagType>(N);
  PerfTest::Trsm<1,15,15, HostSpaceType,AlgoTagType>(N);

  /// Right, Upper, NoTrans, UnitDiag

  PerfTest::Trsm<2, 3, 3, HostSpaceType,AlgoTagType>(N);
  PerfTest::Trsm<2, 5, 5, HostSpaceType,AlgoTagType>(N);
  PerfTest::Trsm<2,10,10, HostSpaceType,AlgoTagType>(N);
  PerfTest::Trsm<2,15,15, HostSpaceType,AlgoTagType>(N);

  /// Right, Upper, NoTrans, NonUnitDiag (used in LU factorization)

  PerfTest::Trsm<3, 3, 3, HostSpaceType,AlgoTagType>(N);
  PerfTest::Trsm<3, 5, 5, HostSpaceType,AlgoTagType>(N);
  PerfTest::Trsm<3,10,10, HostSpaceType,AlgoTagType>(N);
  PerfTest::Trsm<3,15,15, HostSpaceType,AlgoTagType>(N);

  std::cout << "\n\n Used for Solve \n\n";

  /// Left, Lower, NoTrans, UnitDiag (used in LU solve)

  PerfTest::Trsm<0, 3, 1, HostSpaceType,AlgoTagType>(N);
  PerfTest::Trsm<0, 5, 1, HostSpaceType,AlgoTagType>(N);
  PerfTest::Trsm<0,10, 1, HostSpaceType,AlgoTagType>(N);
  PerfTest::Trsm<0,15, 1, HostSpaceType,AlgoTagType>(N);

  /// Left, Upper, Notrans, NonUnitDiag (user in LU solve)

  PerfTest::Trsm<4, 3, 1, HostSpaceType,AlgoTagType>(N);
  PerfTest::Trsm<4, 5, 1, HostSpaceType,AlgoTagType>(N);
  PerfTest::Trsm<4,10, 1, HostSpaceType,AlgoTagType>(N);
  PerfTest::Trsm<4,15, 1, HostSpaceType,AlgoTagType>(N);
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

    std::cout << "\n Testing Algo::Trsm::Unblocked\n";
    run<Algo::Trsm::Unblocked>(N);

    std::cout << "\n Testing Algo::Trsm::Blocked\n";
    run<Algo::Trsm::Blocked>(N);
  }

  Kokkos::finalize();

  return 0;
}
