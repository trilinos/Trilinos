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

#include "dense_matrix_base.hpp"
#include "dense_matrix_view.hpp"

#include "team_view.hpp"
#include "task_view.hpp"

#include "sequential_for.hpp"
#include "parallel_for.hpp"

#include "task_policy_graphviz.hpp"
#include "team_factory.hpp"
#include "task_factory.hpp"
#include "task_team_factory.hpp"

#include "scale.hpp"
#include "tri_solve.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

#define USE_SEQUENTIAL_FOR
typedef Kokkos::Serial space_type; 

//typedef Kokkos::Threads space_type; 
//typedef Kokkos::Qthread space_type; 

using namespace Example;

typedef space_type ExecSpace;

typedef CrsMatrixBase<value_type,ordinal_type,size_type,space_type> CrsMatrixBaseType;
typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;

typedef DenseMatrixBase<value_type,ordinal_type,size_type,space_type> DenseMatrixBaseType;
typedef DenseMatrixView<DenseMatrixBaseType> DenseMatrixViewType;

#ifdef USE_SEQUENTIAL_FOR
typedef TaskTeamFactory<TeamPolicy,Future,TeamThreadLoopRegion> TaskFactoryType;
typedef SequentialFor ForType;
#else
typedef TaskTeamFactory<Kokkos::Experimental::TaskPolicy<space_type>,
                        Kokkos::Experimental::Future<int,space_type>,
                        Kokkos::Impl::TeamThreadRangeBoundariesStruct> TaskFactoryType;
typedef ParallelFor ForType;
#endif

typedef TaskView<CrsMatrixViewType,TaskFactoryType> CrsTaskViewType;
typedef TaskView<DenseMatrixViewType,TaskFactoryType> DenseTaskViewType;

typedef GraphHelper_Scotch<CrsMatrixBaseType> GraphHelperType;

int main (int argc, char *argv[]) {
  if (argc < 3) {
    cout << "Usage: " << argv[0] << " filename nthreads" << endl;
    return -1;
  }

  const int nthreads = atoi(argv[2]);
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
  
  CrsMatrixBaseType UU("UU");
  UU.copy(Uplo::Upper, PA);
  
  cout << UU << endl;
  
  CrsTaskViewType U(&UU);
  U.fillRowViewArray();

  DenseMatrixBaseType BB("BB", UU.NumRows(), 1);
  for (ordinal_type i=0;i<BB.NumRows();++i)
    BB.Value(i, 0) = 1.0;

  DenseTaskViewType B(&BB);

  {
    int r_val = 0;
    typedef typename CrsTaskViewType::policy_type policy_type;

#ifdef USE_SEQUENTIAL_FOR
    TriSolve<Uplo::Upper,Trans::ConjTranspose,AlgoTriSolve::Unblocked>
      ::TaskFunctor<ForType,CrsTaskViewType,DenseTaskViewType>
      (Diag::NonUnit, U, B).apply(policy_type::member_null(), r_val);
#else
    policy_type policy;
    auto future = policy.create_team(TriSolve<Uplo::Upper,Trans::ConjTranspose,AlgoTriSolve::Unblocked>
                                     ::TaskFunctor<ForType,CrsTaskViewType,DenseTaskViewType>
                                     (Diag::NonUnit, U, B), 0);
    
    policy.spawn(future);
    Kokkos::Experimental::wait(policy);
#endif

    if (r_val != 0)  {
      cout << " Error = " << r_val << endl;
      return r_val;
    }

    cout << BB << endl;
  }

  ExecSpace::finalize();

  return 0;
}
