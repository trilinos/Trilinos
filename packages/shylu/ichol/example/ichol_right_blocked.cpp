#include <Kokkos_Core.hpp>

#include <Kokkos_TaskPolicy.hpp>
#include <impl/Kokkos_Serial_TaskPolicy.hpp>

#include <Kokkos_Qthread.hpp>
#include <Qthread/Kokkos_Qthread_TaskPolicy.hpp>

#include <Kokkos_Threads.hpp>
#include <Threads/Kokkos_Threads_TaskPolicy.hpp>

#include "util.hpp"
#include "graph_helper_scotch.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "crs_team_view.hpp"
#include "crs_task_view.hpp"

#include "sequential_for.hpp"
#include "parallel_for.hpp"

#include "task_policy_graphviz.hpp"
#include "team_factory.hpp"
#include "task_factory.hpp"
#include "task_team_factory.hpp"

#include "ichol.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

//#define USE_SEQUENTIAL_FOR
//typedef Kokkos::Serial space_type;

typedef Kokkos::Threads space_type;
//typedef Kokkos::Qthread space_type;

using namespace Example;

typedef space_type ExecSpace;

typedef CrsMatrixBase<value_type,ordinal_type,size_type,space_type> CrsMatrixBaseType;
typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;

#ifdef USE_SEQUENTIAL_FOR
typedef TaskTeamFactory<TeamPolicy,Future,TeamThreadLoopRegion> TaskFactoryType;
typedef SequentialFor ForType;
#else
typedef TaskTeamFactory<Kokkos::Experimental::TaskPolicy<space_type>,
                        Kokkos::Experimental::Future<int,space_type>,
                        Kokkos::Impl::TeamThreadRangeBoundariesStruct> TaskFactoryType;
typedef ParallelFor ForType;
#endif

typedef CrsTaskView<CrsMatrixBaseType,TaskFactoryType> CrsTaskViewType;
typedef GraphHelper_Scotch<CrsMatrixBaseType> GraphHelperType;

int main (int argc, char *argv[]) {
  if (argc < 4) {
    cout << "Usage: " << argv[0] << " filename blksize nthreads" << endl;
    return -1;
  }

  const int blocksize = atoi(argv[2]);
  IChol<Uplo::Upper,AlgoIChol::RightBlocked>::blocksize = blocksize;

  const int nthreads = atoi(argv[3]);
  ExecSpace::initialize(nthreads);
  cout << "Default execution space initialized = "
       << typeid(Kokkos::DefaultExecutionSpace).name()
       << endl;

  CrsMatrixBaseType AA("AA");

  ifstream in;
  in.open(argv[1]);
  if (!in.good()) {
    cout << "Error in open the file: " << argv[1] << endl;
    return -1;
  }
  AA.importMatrixMarket(in);
  cout << AA << endl;

  GraphHelperType S(AA);
  S.computeOrdering();

  CrsMatrixBaseType PA("Permuted AA");
  PA.copy(S.PermVector(), S.InvPermVector(), AA);

  CrsMatrixBaseType UU("Upper Triangular of AA");
  UU.copy(Uplo::Upper, PA);

  CrsTaskViewType U(UU);

  {
    int r_val = 0;
    typedef typename CrsTaskViewType::policy_type policy_type;

#ifdef USE_SEQUENTIAL_FOR
    IChol<Uplo::Upper,AlgoIChol::RightBlocked>
      ::TaskFunctor<CrsTaskViewType,ForType>(U).apply(policy_type::member_null(), r_val);
#else
    policy_type policy;
    auto future = policy.create_team(IChol<Uplo::Upper,AlgoIChol::RightBlocked>
                                     ::TaskFunctor<CrsTaskViewType,ForType>(U), 0);
    policy.spawn(future);
    Kokkos::Experimental::wait(future);
#endif

    if (r_val != 0) {
      cout << " Error = " << r_val << endl;
      return r_val;
    }
    cout << UU << endl;
  }

  ExecSpace::finalize();

  return 0;
}
