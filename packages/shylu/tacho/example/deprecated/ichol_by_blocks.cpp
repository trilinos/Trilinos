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

#include "crs_matrix_helper.hpp"

#include "team_view.hpp"
#include "task_view.hpp"

#include "sequential_for.hpp"
#include "parallel_for.hpp"

#include "task_policy_graphviz.hpp"
#include "task_factory.hpp"
#include "team_factory.hpp"
#include "task_team_factory.hpp"

#include "ichol.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

// space
//#define USE_SEQUENTIAL_FOR
//typedef Kokkos::Serial space_type;

typedef Kokkos::Threads space_type;
//typedef Kokkos::Qthread space_type;

using namespace Example;

// exec space
typedef space_type ExecSpace;

// flat matrix
typedef CrsMatrixBase<value_type,ordinal_type,size_type,space_type> CrsMatrixBaseType;
typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;

// scotch reordering
typedef GraphHelper_Scotch<CrsMatrixBaseType> GraphHelperType;

#ifdef USE_SEQUENTIAL_FOR
typedef TaskTeamFactory<TaskPolicy,Future,TeamThreadLoopRegion> TaskFactoryType;
typedef SequentialFor ForType;
#else
typedef TaskTeamFactory<Kokkos::Experimental::TaskPolicy<space_type>,
                        Kokkos::Experimental::Future<int,space_type>,
                        Kokkos::Impl::TeamThreadRangeBoundariesStruct> TaskFactoryType;
typedef ParallelFor ForType;
#endif

// block representation for CrsMatrix
typedef TaskView<CrsMatrixViewType,TaskFactoryType> CrsTaskViewType;

// hier matrix
typedef CrsMatrixBase<CrsTaskViewType,ordinal_type,size_type,space_type> CrsHierBaseType;
typedef CrsMatrixView<CrsHierBaseType> CrsHierViewType;

typedef TaskView<CrsHierViewType,TaskFactoryType> CrsHierTaskType;

int main (int argc, char *argv[]) {
  if (argc < 3) {
    cout << "Usage: " << argv[0] << " filename nthreads" << endl;
    return -1;
  }

  const int nthreads = atoi(argv[2]);
  ExecSpace::initialize(nthreads);

#ifdef USE_SEQUENTIAL_FOR
#else
  ExecSpace::print_configuration(std::cout, true);
#endif

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

  GraphHelperType S(AA);
  S.computeOrdering();

  CrsMatrixBaseType PA("Permuted AA");
  PA.copy(S.PermVector(), S.InvPermVector(), AA);

  CrsMatrixBaseType UU("UU");
  UU.copy(Uplo::Upper, PA);

  cout << UU << endl;

  CrsHierBaseType HH("HH");

  CrsMatrixHelper::flat2hier(Uplo::Upper, UU, HH,
                             S.NumBlocks(),
                             S.RangeVector(),
                             S.TreeVector());

  for (ordinal_type k=0;k<HH.NumNonZeros();++k)
    HH.Value(k).fillRowViewArray();

  cout << HH << endl;

  CrsHierTaskType H(&HH);

  {
    int r_val = 0;
    typedef typename CrsTaskViewType::policy_type policy_type;

    IChol<Uplo::Upper,AlgoIChol::ByBlocks>::
      TaskFunctor<ForType,CrsHierTaskType>(H).apply(policy_type::member_null(), r_val);

#ifdef USE_SEQUENTIAL_FOR
    // do nothing
#else
    policy_type policy;
    Kokkos::Experimental::wait(policy);
#endif

    cout << UU << endl;

#ifdef USE_SEQUENTIAL_FOR
    ofstream out;
    out.open("graph.gv");
    if (!out.good()) {
      cout << "Error in open the file: task_graph.gv" << endl;
      return -1;
    }

    policy_type policy;
    policy.graphviz(out);
#endif
  }

  ExecSpace::finalize();

  return 0;
}
