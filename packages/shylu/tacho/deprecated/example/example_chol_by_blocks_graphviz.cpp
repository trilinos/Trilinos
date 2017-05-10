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

#include "task_view.hpp"

#include "task_policy_graphviz.hpp"

#include "task_factory.hpp"

#include "chol.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

// space
typedef Kokkos::Serial space_type;

using namespace Tacho;

// exec space
typedef space_type ExecSpace;

// flat matrix
typedef CrsMatrixBase<value_type,ordinal_type,size_type,space_type> CrsMatrixBaseType;
typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;

// scotch reordering
typedef GraphHelper_Scotch<CrsMatrixBaseType> GraphHelperType;

typedef TaskFactory<TaskPolicy,Future> TaskFactoryType;

// block representation for CrsMatrix
typedef TaskView<CrsMatrixViewType,TaskFactoryType> CrsTaskViewType;

// hier matrix
typedef CrsMatrixBase<CrsTaskViewType,ordinal_type,size_type,space_type> CrsHierBaseType;
typedef CrsMatrixView<CrsHierBaseType> CrsHierViewType;

typedef TaskView<CrsHierViewType,TaskFactoryType> CrsHierTaskType;

// ---------------------------------------------------------------------------------
int main (int argc, char *argv[]) {

  Teuchos::CommandLineProcessor clp;
  clp.setDocString("This example program demonstrates CholByBlocks task dag visualization");

  string file_input = "test.mtx";
  clp.setOption("file-input", &file_input, "Input file (MatrixMarket SPD matrix)");

  string file_output = "chol.gv";
  clp.setOption("file-output", &file_output, "Output file (dot file)");

  int treecut = 10;
  clp.setOption("treecut", &treecut, "Level to cut tree from bottom");

  int minblksize = 0;
  clp.setOption("minblksize", &minblksize, "Minimum block size for internal reordering");

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
  cout << "CholByBlocks::graphviz:: import input file = " << file_input << endl;
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
  cout << "CholByBlocks::graphviz:: import input file::time = " << t << endl;

  // reorder a matrix using Scotch
  // --------------------------------------
  cout << "CholByBlocks::graphviz:: reorder the matrix" << endl;
  CrsMatrixBaseType UU("UU");
  CrsHierBaseType HH("HH");
  {
    timer.reset();

    GraphHelperType::size_type_array rptr(AA.Label()+"Graph::RowPtrArray", AA.NumRows() + 1);
    GraphHelperType::ordinal_type_array cidx(AA.Label()+"Graph::ColIndexArray", AA.NumNonZeros());

    AA.convertGraph(rptr, cidx);
    GraphHelperType S(AA.Label()+"ScotchHelper",
                      AA.NumRows(),
                      rptr,
                      cidx);
    S.computeOrdering(treecut, minblksize);

    cout << "CholByBlocks::graphviz:: "
         << "# of rows = " << S.NumRows() << ", # of blocks " << S.NumBlocks()
         << endl;

    CrsMatrixBaseType PA("Permuted A");
    PA.copy(S.PermVector(), S.InvPermVector(), AA);

    UU.copy(Uplo::Upper, PA);

    CrsMatrixHelper::flat2hier(Uplo::Upper, UU, HH,
                               S.NumBlocks(),
                               S.RangeVector(),
                               S.TreeVector());

    cout << "CholByBlocks::graphviz:: "
         << "# of nnz in Hier = " << HH.NumNonZeros()
         << endl;

    t = timer.seconds();
  }
  cout << "CholByBlocks::graphviz:: reorder the matrix::time = " << t << endl;

  // run chol by blocks
  // --------------------------------------
  {
    typename TaskFactoryType::policy_type policy;
    TaskFactoryType::setPolicy(&policy);

    cout << "CholByBlocks::graphviz:: factorize the matrix" << endl;
    timer.reset();

    int r_val = 0;
    CrsHierTaskType H(&HH);
    Chol<Uplo::Upper,AlgoChol::ByBlocks>::
      TaskFunctor<CrsHierTaskType>(H).apply(r_val);

    ofstream out;
    out.open(file_output);
    if (!out.good()) {
      cout << "Error in open the file: chol.gv" << endl;
      return -1;
    }

    TaskFactoryType::Policy().graphviz(out);
    cout << "CholByBlocks::graphviz:: size of queue = " << TaskFactoryType::Policy().size() << endl;
    TaskFactoryType::Policy().clear();

    t = timer.seconds();
    cout << "CholByBlocks::graphviz:: factorize the matrix::time = " << t << endl;
  }

  // finalization
  // --------------------------------------
  ExecSpace::finalize();

  return 0;
}
