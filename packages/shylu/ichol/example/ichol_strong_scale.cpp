#include <Kokkos_Core.hpp>

#include <Kokkos_TaskPolicy.hpp>
#include <impl/Kokkos_Serial_TaskPolicy.hpp>
#include <impl/Kokkos_Timer.hpp>

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

// helper to convert a flat matrix to a hier matrix
typedef Example::CrsMatrixHelper CrsMatrixHelper; 

// kokkos tasking with future
//#define USE_GRAPHVIZ
#ifdef USE_GRAPHVIZ
typedef Example::TaskFactory<Example::TaskPolicy,
                             Example::Future> TaskFactory;
#else
typedef Example::TaskFactory<Kokkos::Experimental::TaskPolicy<space_type>,
                             Kokkos::Experimental::Future<int,space_type> > TaskFactory;
using Kokkos::Experimental::wait;
#endif

// member submatrix with future inside
typedef Example::CrsTaskView<CrsMatrixBase,TaskFactory> CrsTaskView;

// hier matrix
typedef Example::CrsMatrixBase<CrsTaskView,ordinal_type,size_type,space_type> CrsHierBase;
typedef Example::CrsTaskView<CrsHierBase,TaskFactory> CrsHierView;

// constant space
typedef Example::Uplo Uplo;
typedef Example::AlgoIChol AlgoIChol;

// driver
using Example::IChol;


// ---------------------------------------------------------------------------------
int main (int argc, char *argv[]) {

  if (argc < 3) {
    cout << "Usage: " << argv[0] << " nthreads filename" << endl;
    return -1;
  }

  Kokkos::Impl::Timer timer;
  double t[10];
  int cnt = 0;

  const int nthreads = stoi(argv[1]);

  // --------------------------------------
  timer.reset();

  Kokkos::Qthread::initialize(nthreads);
  Kokkos::Qthread::print_configuration(cout, true);

  t[cnt] = timer.seconds();
  // --------------------------------------

  cout << "default execution space initialized = "
       << typeid(Kokkos::DefaultExecutionSpace).name()
       << endl
       << "# of threads = "
       << nthreads 
       << endl;

  ifstream in;
  in.open(argv[2]);
  if (!in.good()) {
    cout << "Error in open the file: " << argv[2] << endl;
    return -1;
  }

  // --------------------------------------
  timer.reset();

  CrsMatrixBase A("A");
  A.importMatrixMarket(in);

  t[++cnt] = timer.seconds();
  // --------------------------------------

  cout << "time for importing a sparse matrix = "
       << t[cnt] 
       << endl;

  // --------------------------------------
  timer.reset();

  GraphHelper S(A);
  S.computeOrdering();

  t[++cnt] = timer.seconds();
  // --------------------------------------
  
  cout << "time for computing ordering (scotch) = "
       << t[cnt]
       << endl
       << "# of rows = " << S.NumRows() << ", # of blocks " << S.NumBlocks() << " " 
       << endl;

  // --------------------------------------
  timer.reset();

  CrsMatrixBase PA("Permuted A");
  PA.copy(S.PermVector(), S.InvPermVector(), A);

  t[++cnt] = timer.seconds();
  // --------------------------------------
  
  cout << "time for permuting matrix = "
       << t[cnt]
       << endl;

  // --------------------------------------
  CrsMatrixBase R("R");

  R.copy(Uplo::Upper, PA);
  {
    timer.reset();
    
    CrsHierBase H("H");
    CrsMatrixHelper::flat2hier(Uplo::Upper, R, H,
                               S.NumBlocks(),
                               S.RangeVector(),
                               S.TreeVector());
    
    for (ordinal_type k=0;k<H.NumNonZeros();++k) 
      H.Value(k).fillRowViewArray();
    
    t[++cnt] = timer.seconds();
    // --------------------------------------
    
    cout << "time for creating partitioned matrices = "
         << t[cnt]
         << endl;
    
    // --------------------------------------
    timer.reset();
    
    CrsHierView HH(H);
    IChol<Uplo::Upper,AlgoIChol::RightByBlocks>::invoke(HH);
    
    t[++cnt] = timer.seconds();
    // --------------------------------------
    
    cout << "time for ichol qthreads task gen = "
         << t[cnt]
         << endl;

    // --------------------------------------
    timer.reset();
    
#ifdef USE_GRAPHVIZ
    ofstream out;
    out.open("graph.gv");
    if (!out.good()) {
      cout << "Error in open the file: graph.gv" << endl;
      return -1;
    }
    
    TaskFactory::policy_type policy;
    policy.graphviz(out);
#else
    for (ordinal_type k=0;k<H.NumNonZeros();++k) 
      wait(H.Value(k).Future());
#endif
    
    t[++cnt] = timer.seconds();
    
    // --------------------------------------
    
    cout << "time for ichol qthreads wait = "
         << t[cnt]
         << endl;
  }
  
  R.copy(Uplo::Upper, PA);
  {
    CrsMatrixView RR(R);
    RR.fillRowViewArray();
    
    // --------------------------------------
    timer.reset();
    
    IChol<Uplo::Upper,AlgoIChol::RightUnblockedOpt1>::invoke(RR);
    
    t[++cnt] = timer.seconds();
    // --------------------------------------
    
    cout << "time for ichol sequential unblocked opt1 = "
         << t[cnt] 
         << endl;
    
    cout << "scale [qthread/sequential] = "
         << t[cnt]/(t[cnt-1] + t[cnt-2])
         << endl;
  }
  
  timer.reset();

  Kokkos::Qthread::finalize();

  t[0] += timer.seconds();
  // --------------------------------------

  cout << "time for Kokkos init/finalize = "
       << t[0]
       << endl;
  
  return 0;
}
