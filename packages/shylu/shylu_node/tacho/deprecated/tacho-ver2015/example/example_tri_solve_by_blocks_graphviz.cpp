#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include <Kokkos_TaskPolicy.hpp>
#include <impl/Kokkos_Serial_TaskPolicy.hpp>

#include "Teuchos_CommandLineProcessor.hpp"

#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "dense_matrix_base.hpp"
#include "dense_matrix_view.hpp"

#include "graph_helper_scotch.hpp"
#include "crs_matrix_helper.hpp"
#include "dense_matrix_helper.hpp"

#include "task_view.hpp"

#include "task_policy_graphviz.hpp"

#include "task_factory.hpp"

#include "tri_solve.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

// space
typedef Kokkos::Serial space_type;

using namespace Tacho;

// exec space
typedef space_type ExecSpace;

// graphviz task policy
typedef TaskFactory<TaskPolicy,Future> TaskFactoryType;

// flat crs matrix
typedef CrsMatrixBase<value_type,ordinal_type,size_type,space_type> CrsMatrixBaseType;
typedef GraphHelper_Scotch<CrsMatrixBaseType> GraphHelperType;

// block representation of crs matrix
typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;
typedef TaskView<CrsMatrixViewType,TaskFactoryType> CrsTaskViewType;

// hier crs matrix
typedef CrsMatrixBase<CrsTaskViewType,ordinal_type,size_type,space_type> CrsHierMatrixBaseType;

// hier task view
typedef CrsMatrixView<CrsHierMatrixBaseType> CrsHierMatrixViewType;
typedef TaskView<CrsHierMatrixViewType,TaskFactoryType> CrsHierTaskViewType;

// flat dense matrix
typedef DenseMatrixBase<value_type,ordinal_type,size_type,space_type> DenseMatrixBaseType;

// block representation of dense matrix
typedef DenseMatrixView<DenseMatrixBaseType> DenseMatrixViewType;
typedef TaskView<DenseMatrixViewType,TaskFactoryType> DenseTaskViewType;

// hier dense matrix
typedef DenseMatrixBase<DenseTaskViewType,ordinal_type,size_type,space_type> DenseHierMatrixBaseType;

// hier task view
typedef DenseMatrixView<DenseHierMatrixBaseType> DenseHierMatrixViewType;
typedef TaskView<DenseHierMatrixViewType,TaskFactoryType> DenseHierTaskViewType;

// ---------------------------------------------------------------------------------
int main (int argc, char *argv[]) {

  Teuchos::CommandLineProcessor clp;
  clp.setDocString("This example program demonstrates TriSolveByBlocks task dag visualization");

  string file_input = "test.mtx";
  clp.setOption("file-input", &file_input, "Input file (MatrixMarket SPD matrix)");

  string file_output = "chol.gv";
  clp.setOption("file-output", &file_output, "Output file (dot file)");

  int nrhs = 1;
  clp.setOption("nrhs", &nrhs, "Number of right hand side");

  int nb = nrhs;
  clp.setOption("nb", &nb, "Blocksize of right hand side");

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
  cout << "TriSolveByBlocks::graphviz:: import input file = " << file_input << endl;
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
  cout << "TriSolveByBlocks::graphviz:: import input file::time = " << t << endl;

  // reorder a matrix using Scotch
  // --------------------------------------
  cout << "TriSolveByBlocks::graphviz:: reorder the matrix" << endl;
  CrsMatrixBaseType UU("UU");
  DenseMatrixBaseType BB("BB", AA.NumRows(), nrhs);

  CrsHierMatrixBaseType HU("HU");
  DenseHierMatrixBaseType HB("HB");
  {
    timer.reset();

    typename GraphHelperType::size_type_array rptr(AA.Label()+"Graph::RowPtrArray", AA.NumRows() + 1);
    typename GraphHelperType::ordinal_type_array cidx(AA.Label()+"Graph::ColIndexArray", AA.NumNonZeros());

    AA.convertGraph(rptr, cidx);
    GraphHelperType S(AA.Label()+"ScotchHelper",
                      AA.NumRows(),
                      rptr,
                      cidx);
    S.computeOrdering();

    cout << "TriSolveByBlocks::graphviz:: "
         << "# of rows = " << S.NumRows() << ", # of blocks " << S.NumBlocks()
         << endl;

    CrsMatrixBaseType PA("Permuted A");
    PA.copy(S.PermVector(), S.InvPermVector(), AA);

    UU.copy(Uplo::Upper, PA);

    CrsMatrixHelper::flat2hier(Uplo::Upper, UU, HU,
                               S.NumBlocks(),
                               S.RangeVector(),
                               S.TreeVector());

    DenseMatrixHelper::flat2hier(BB, HB,
                                 S.NumBlocks(),
                                 S.RangeVector(),
                                 nb);
    t = timer.seconds();
  }
  cout << "TriSolveByBlocks::graphviz:: reorder the matrix::time = " << t << endl;

  // run tri solve by blocks
  // --------------------------------------
  {
    typename TaskFactoryType::policy_type policy;
    TaskFactoryType::setPolicy(&policy);

    cout << "TriSolveByBlocks::graphviz:: forward/backward solve the matrix" << endl;
    timer.reset();

    int r_val = 0;
    CrsHierTaskViewType TU(&HU);
    DenseHierTaskViewType TB(&HB);

    TaskFactoryType::Policy().set_work_phase(2);
    TriSolve<Uplo::Upper,Trans::ConjTranspose,AlgoTriSolve::ByBlocks>
      ::TaskFunctor<CrsHierTaskViewType,DenseHierTaskViewType>
      (Diag::NonUnit, TU, TB).apply(r_val);

    TaskFactoryType::Policy().set_work_phase(3);
    TriSolve<Uplo::Upper,Trans::NoTranspose,AlgoTriSolve::ByBlocks>
      ::TaskFunctor<CrsHierTaskViewType,DenseHierTaskViewType>
      (Diag::NonUnit, TU, TB).apply(r_val);

    ofstream out;
    out.open(file_output);
    if (!out.good()) {
      cout << "Error in open the file: " << file_output << endl;
      return -1;
    }

    TaskFactoryType::Policy().graphviz(out);
    TaskFactoryType::Policy().clear();

    t = timer.seconds();
    cout << "TriSolveByBlocks::graphviz:: forward/backward solve the matrix::time = " << t << endl;
  }

  // finalization
  // --------------------------------------
  ExecSpace::finalize();

  return 0;
}
