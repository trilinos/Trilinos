#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

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

// task team factory
typedef TaskTeamFactory<Kokkos::Experimental::TaskPolicy<space_type>,
                        Kokkos::Experimental::Future<int,space_type>,
                        Kokkos::Impl::TeamThreadRangeBoundariesStruct> TaskFactoryType;
typedef ParallelFor ForType;

// block representation for CrsMatrix
typedef TaskView<CrsMatrixViewType,TaskFactoryType> CrsTaskViewType;

// hier matrix
typedef CrsMatrixBase<CrsTaskViewType,ordinal_type,size_type,space_type> CrsHierBaseType;
typedef CrsMatrixView<CrsHierBaseType> CrsHierViewType;

typedef TaskView<CrsHierViewType,TaskFactoryType> CrsHierTaskType;

#define DOTLINE "====================================="

// ---------------------------------------------------------------------------------
int main (int argc, char *argv[]) {
  if (argc < 4) {
    cout << "Usage: " << argv[0] << " filename nthreads teamsize" << endl;
    return -1;
  }

  Kokkos::Impl::Timer timer;

  double t[10];
  int cnt = 0;

  const int nthreads = stoi(argv[2]);
  const int teamsize = stoi(argv[3]);

  // initialization
  // --------------------------------------
  timer.reset();

  ExecSpace::initialize(nthreads);
  ExecSpace::print_configuration(std::cout, true);

  t[cnt] = timer.seconds();

  cout << "default execution space initialized = "
       << typeid(Kokkos::DefaultExecutionSpace).name()
       << endl
       << "# of threads = "
       << nthreads
       << endl;

  // import a matrix
  // --------------------------------------
  ifstream in;
  in.open(argv[1]);
  if (!in.good()) {
    cout << "Error in open the file: " << argv[1] << endl;
    return -1;
  }

  timer.reset();

  CrsMatrixBaseType AA("AA");
  AA.importMatrixMarket(in);

  t[++cnt] = timer.seconds();

  cout << "time for importing a sparse matrix = "
       << t[cnt]
       << endl;

  // reorder a matrix using Scotch
  // --------------------------------------
  timer.reset();

  GraphHelperType S(AA);
  S.computeOrdering();

  t[++cnt] = timer.seconds();

  cout << "time for computing ordering (scotch) = "
       << t[cnt]
       << endl
       << "# of rows = " << S.NumRows() << ", # of blocks " << S.NumBlocks() << " "
       << endl;

  // permute a matrix with a given ordering
  // --------------------------------------
  timer.reset();

  CrsMatrixBaseType PA("Permuted A");
  PA.copy(S.PermVector(), S.InvPermVector(), AA);

  t[++cnt] = timer.seconds();

  cout << "time for permuting matrix = "
       << t[cnt]
       << endl;

  {
    // run ichol by blocks
    // --------------------------------------
    CrsMatrixBaseType RR("RR");
    
    RR.copy(Uplo::Upper, PA);
    {
      typename TaskFactoryType::policy_type policy(20, teamsize);
      TaskFactoryType::setPolicy(&policy);

      timer.reset();
      
      CrsHierBaseType HH("HH");
      CrsMatrixHelper::flat2hier(Uplo::Upper, RR, HH,
                                 S.NumBlocks(),
                                 S.RangeVector(),
                                 S.TreeVector());
      
      for (ordinal_type k=0;k<HH.NumNonZeros();++k)
        HH.Value(k).fillRowViewArray();
      
      t[++cnt] = timer.seconds();
      
      cout << "time for creating partitioned matrices = "
           << t[cnt]
           << endl;
      
      timer.reset();
      
      CrsHierTaskType H(&HH);
      
      int r_val = 0;
      
      IChol<Uplo::Upper,AlgoIChol::ByBlocks>::
        TaskFunctor<ForType,CrsHierTaskType>(H).apply(r_val);
      
      t[++cnt] = timer.seconds();
      
      cout << "time for task gen in ichol_by_blocks = "
           << t[cnt]
           << endl;
      
      timer.reset();
      
      Kokkos::Experimental::wait(TaskFactoryType::Policy());
      
      t[++cnt] = timer.seconds();
      
      cout << "time for wait generated tasks = "
           << t[cnt]
           << endl;
      
      // cout << RR << endl;
    }
    
    // numeric factorization only
    double t_task_team_by_blocks = t[cnt] + t[cnt-1];
    
    // data parallel using team threads
    // --------------------------------------
    RR.copy(Uplo::Upper, PA);
    {
      typename TaskFactoryType::policy_type policy(20, nthreads);
      TaskFactoryType::setPolicy(&policy);

      CrsTaskViewType R(&RR);
      R.fillRowViewArray();
      
      timer.reset();
      
      auto future = TaskFactoryType::Policy().create_team(IChol<Uplo::Upper,AlgoIChol::UnblockedOpt1>
                                                          ::TaskFunctor<ForType,CrsTaskViewType>(R), 0);
      TaskFactoryType::Policy().spawn(future);
      Kokkos::Experimental::wait(TaskFactoryType::Policy());
      
      t[++cnt] = timer.seconds();
      
      cout << "time for ichol data parallel unblocked opt1 team parallel = "
           << t[cnt]
           << endl;
      
      // cout << RR << endl;
    }
    
    // numeric factorization only
    double t_team_parallel = t[cnt];

    // sequential
    // ----------
    RR.copy(Uplo::Upper, PA);
    {
      typename TaskFactoryType::policy_type policy(20, 1);
      TaskFactoryType::setPolicy(&policy);

      CrsTaskViewType R(&RR);
      R.fillRowViewArray();
      
      timer.reset();
      
      auto future = TaskFactoryType::Policy().create(IChol<Uplo::Upper,AlgoIChol::UnblockedOpt1>
                                                     ::TaskFunctor<ForType,CrsTaskViewType>(R), 0);
      TaskFactoryType::Policy().spawn(future);
      Kokkos::Experimental::wait(TaskFactoryType::Policy());
      
      t[++cnt] = timer.seconds();
      
      cout << "time for ichol sequential unblocked opt1 team parallel = "
           << t[cnt]
           << endl;
      
      // cout << RR << endl;
    }
    
    // numeric factorization only
    double t_sequential = t[cnt];
    
    // result
    // --------------------------------------
    cout << "task data scale [sequential/task-team-by-blocks] = "
         << t_sequential/t_task_team_by_blocks
         << endl;
    cout << "     data scale [sequential/team-parallel] = "
         << t_sequential/t_team_parallel
         << endl;
  }

  // finalization
  // --------------------------------------
  timer.reset();

  ExecSpace::finalize();

  t[0] += timer.seconds();

  cout << "time for Kokkos init/finalize = "
       << t[0]
       << endl;
  // --------------------------------------

  return 0;
}
