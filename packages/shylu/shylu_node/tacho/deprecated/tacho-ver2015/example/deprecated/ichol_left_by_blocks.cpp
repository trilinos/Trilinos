#include <Kokkos_Core.hpp>
#include <Kokkos_TaskPolicy.hpp>
#include <impl/Kokkos_Serial_TaskPolicy.hpp>

#include <Kokkos_Qthread.hpp>
#include <Qthread/Kokkos_Qthread_TaskPolicy.hpp>
                                                           
#include "util.hpp"
#include "graph_helper_scotch.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "crs_matrix_helper.hpp"
#include "crs_task_view.hpp"

#include "task_policy_graphviz.hpp"
#include "task_factory.hpp"

#include "ichol.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

// space 
//typedef Kokkos::Serial space_type;
typedef Kokkos::Qthread space_type;

// flat matrix 
typedef Example::CrsMatrixBase<value_type,ordinal_type,size_type,space_type> CrsMatrixBase;
typedef Example::CrsMatrixView<CrsMatrixBase> CrsMatrixView;

// scotch reordering
typedef Example::GraphHelper_Scotch<CrsMatrixBase> GraphHelper;

//#define USE_GRAPHVIZ

#ifdef USE_GRAPHVIZ
typedef Example::TaskFactory<Example::TaskPolicy,
                             Example::Future> TaskFactory;
#else
typedef Example::TaskFactory<Kokkos::TaskPolicy<space_type>,
                             Kokkos::Future<int,space_type> > TaskFactory;
using Kokkos::wait;
#endif

// flat2hier
typedef Example::CrsMatrixHelper CrsMatrixHelper; 

// block representation for CrsMatrix
typedef Example::CrsTaskView<CrsMatrixBase,TaskFactory> CrsTaskView;

// hier matrix
typedef Example::CrsMatrixBase<CrsTaskView,ordinal_type,size_type,space_type> CrsHierBase;
typedef Example::CrsTaskView<CrsHierBase,TaskFactory> CrsHierView;

typedef Example::Uplo Uplo;
typedef Example::AlgoIChol AlgoIChol;

using Example::IChol;

int main (int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " filename" << endl;
    return -1;
  }

  // Kokkos::initialize();
  const int nthreads = 16;
  Kokkos::Qthread::initialize(nthreads);
  Kokkos::Qthread::print_configuration(std::cout, true);

  cout << "Default execution space initialized = "
       << typeid(Kokkos::DefaultExecutionSpace).name()
       << endl;

  CrsMatrixBase AA("AA");

  ifstream in;
  in.open(argv[1]);
  if (!in.good()) {
    cout << "Error in open the file: " << argv[1] << endl;
    return -1;
  }
  AA.importMatrixMarket(in);

  GraphHelper S(AA);
  S.computeOrdering();

  CrsMatrixBase PA("Permuted AA");
  PA.copy(S.PermVector(), S.InvPermVector(), AA);
  
  CrsMatrixBase LL("LL");
  LL.copy(Uplo::Lower, PA);

  cout << LL << endl;

  CrsHierBase HH("HH");

  CrsMatrixHelper::flat2hier(Uplo::Lower, LL, HH,
                             S.NumBlocks(),
                             S.RangeVector(),
                             S.TreeVector());

  cout << HH << endl;

  CrsHierView H(HH);

  IChol<Uplo::Lower,AlgoIChol::LeftByBlocks>::invoke(H);

#ifdef USE_GRAPHVIZ
  // do nothing
#else
  for (ordinal_type k=0;k<HH.NumNonZeros();++k) 
    wait(HH.Value(k).Future());
#endif

  cout << LL << endl;

#ifdef USE_GRAPHVIZ
  ofstream out;
  out.open("graph_left.gv");
  if (!out.good()) {
    cout << "Error in open the file: task_graph.gv" << endl;
    return -1;
  }

  TaskFactory::policy_type policy;
  policy.graphviz(out);
#endif

  Kokkos::Qthread::finalize();
  //Kokkos::finalize(); 

  return 0;
}
