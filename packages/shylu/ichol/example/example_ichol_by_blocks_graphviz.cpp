#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>     

#include <Kokkos_TaskPolicy.hpp>
#include <impl/Kokkos_Serial_TaskPolicy.hpp>

#include "Teuchos_CommandLineProcessor.hpp"

#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "graph_helper_scotch.hpp"
#include "crs_matrix_helper.hpp"

#include "team_view.hpp"
#include "task_view.hpp"

#include "sequential_for.hpp"
#include "task_policy_graphviz.hpp"

#include "team_factory.hpp"
#include "task_factory.hpp"
#include "task_team_factory.hpp"

#include "ichol.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

// space
typedef Kokkos::Serial space_type;

using namespace Example;

// exec space
typedef space_type ExecSpace;

// flat matrix
typedef CrsMatrixBase<value_type,ordinal_type,size_type,space_type> CrsMatrixBaseType;
typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;

// scotch reordering
typedef GraphHelper_Scotch<CrsMatrixBaseType> GraphHelperType;

typedef TaskTeamFactory<TaskPolicy,Future,TeamThreadLoopRegion> TaskFactoryType;
typedef SequentialFor ForType;

// block representation for CrsMatrix
typedef TaskView<CrsMatrixViewType,TaskFactoryType> CrsTaskViewType;

// hier matrix
typedef CrsMatrixBase<CrsTaskViewType,ordinal_type,size_type,space_type> CrsHierBaseType;
typedef CrsMatrixView<CrsHierBaseType> CrsHierViewType;

typedef TaskView<CrsHierViewType,TaskFactoryType> CrsHierTaskType;

// ---------------------------------------------------------------------------------
int main (int argc, char *argv[]) {

  Teuchos::CommandLineProcessor clp;
  clp.setDocString("This example program demonstrates ICholByBlocks task dag visualization");

  string file_input = "test.mtx";
  clp.setOption("file-input", &file_input, "Input file (MatrixMarket SPD matrix)");

  string file_output = "ichol.gv";
  clp.setOption("file-output", &file_output, "Output file (dot file)");

  clp.recogniseAllOptions(true);
  clp.throwExceptions(false);

  Teuchos::CommandLineProcessor::EParseCommandLineReturn r_parse= clp.parse( argc, argv );

  if (r_parse == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return 0;
  if (r_parse != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL  ) return -1;

  // initialization
  // --------------------------------------
  ExecSpace::initialize();

  Kokkos::Impl::Timer timer;
  double t = 0.0;

  // import a matrix
  // --------------------------------------
  cout << "ICholByBlocks::graphviz:: import input file = " << file_input << endl;
  CrsMatrixBaseType AA("AA");
  {
    timer.reset();

    ifstream in;
    in.open(file_input);
    if (!in.good()) {
      cout << "Error in open the file: " << file_input << endl;
      return -1;
    }
    AA.importMatrixMarket(in);

    t = timer.seconds();
  }
  cout << "ICholByBlocks::graphviz:: import input file::time = " << t << endl;

  // reorder a matrix using Scotch
  // --------------------------------------
  cout << "ICholByBlocks::graphviz:: reorder the matrix" << endl;
  CrsMatrixBaseType UU("UU");
  CrsHierBaseType HH("HH");
  {
    timer.reset();

    GraphHelperType S(AA);
    S.computeOrdering();

    cout << "ICholByBlocks::graphviz:: "
         << "# of rows = " << S.NumRows() << ", # of blocks " << S.NumBlocks()
         << endl;

    CrsMatrixBaseType PA("Permuted A");
    PA.copy(S.PermVector(), S.InvPermVector(), AA);

    UU.copy(Uplo::Upper, PA);

    CrsMatrixHelper::flat2hier(Uplo::Upper, UU, HH,
                               S.NumBlocks(),
                               S.RangeVector(),
                               S.TreeVector());

    t = timer.seconds();
  }
  cout << "ICholByBlocks::graphviz:: reorder the matrix::time = " << t << endl;

  // run ichol by blocks
  // --------------------------------------
  {
    typename TaskFactoryType::policy_type policy;
    TaskFactoryType::setPolicy(&policy);

    cout << "ICholByBlocks::graphviz:: factorize the matrix" << endl;
    timer.reset();

    int r_val = 0;
    CrsHierTaskType H(&HH);
    IChol<Uplo::Upper,AlgoIChol::ByBlocks>::
      TaskFunctor<ForType,CrsHierTaskType>(H).apply(r_val);

    ofstream out;
    out.open(file_output);
    if (!out.good()) {
      cout << "Error in open the file: ichol.gv" << endl;
      return -1;
    }

    TaskFactoryType::Policy().graphviz(out);
    TaskFactoryType::Policy().clear();

    t = timer.seconds();
    cout << "ICholByBlocks::graphviz:: factorize the matrix::time = " << t << endl;
  }

  // finalization
  // --------------------------------------
  ExecSpace::finalize();

  return 0;
}
