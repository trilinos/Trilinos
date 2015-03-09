#include <Kokkos_Core.hpp>
#include <Kokkos_TaskPolicy.hpp>
#include <impl/Kokkos_Serial_TaskPolicy.hpp>

#include <Kokkos_Qthread.hpp>
#include <Qthread/Kokkos_Qthread_TaskPolicy.hpp>
                                                           
#include "util.hpp"
#include "sequential_for.hpp"
#include "graph_helper_scotch.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "crs_matrix_helper.hpp"
#include "crs_task_view.hpp"

#include "task_policy_graphviz.hpp"
#include "task_factory.hpp"
#include "task_team_factory.hpp"

#include "ichol.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

// space 
typedef Kokkos::Serial space_type;
//typedef Kokkos::Qthread space_type;

using namespace Example;

// flat matrix 
typedef CrsMatrixBase<value_type,ordinal_type,size_type,space_type> CrsMatrixBaseType;
typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;

// scotch reordering
typedef GraphHelper_Scotch<CrsMatrixBaseType> GraphHelperType;

#define USE_GRAPHVIZ

#ifdef USE_GRAPHVIZ
typedef TaskTeamFactory<TaskPolicy,Future,TeamThreadLoopRegion> TaskFactoryType;
#else
typedef TaskTeamFactory<Kokkos::Experimental::TaskPolicy<space_type>,
                        Kokkos::Experimental::Future<int,space_type>,
                        TeamThreadLoopRegion> TaskFactoryType;
using Kokkos::Experimental::wait;
#endif

// block representation for CrsMatrix
typedef CrsTaskView<CrsMatrixBaseType,TaskFactoryType> CrsTaskViewType;

// hier matrix
typedef CrsMatrixBase<CrsTaskViewType,ordinal_type,size_type,space_type> CrsHierBaseType;
typedef CrsTaskView<CrsHierBaseType,TaskFactoryType> CrsHierViewType;

int main (int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " filename" << endl;
    return -1;
  }

  // Kokkos::initialize();
  const int threads_count = 16;
  Kokkos::Qthread::initialize( threads_count );
  Kokkos::Qthread::print_configuration( std::cout , true );

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

  cout << HH << endl;

  CrsHierViewType H(HH);

  {
    typedef typename CrsTaskViewType::policy_type::member_type member_type;

    IChol<Uplo::Upper,AlgoIChol::RightByBlocks>::invoke<CrsHierViewType,SequentialFor>(member_type(), H);
    
#ifdef USE_GRAPHVIZ
    // do nothing
#else
    for (ordinal_type k=0;k<HH.NumNonZeros();++k) 
      wait(HH.Value(k).Future());
#endif
    
    cout << UU << endl;
    
#ifdef USE_GRAPHVIZ
    ofstream out;
    out.open("graph_right.gv");
    if (!out.good()) {
      cout << "Error in open the file: task_graph.gv" << endl;
      return -1;
    }
    
    TaskFactoryType::policy_type policy;
    policy.graphviz(out);
#endif
  }
  
  Kokkos::Qthread::finalize();
  //Kokkos::finalize(); 
  
  return 0;
}
