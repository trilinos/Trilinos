#include "Teuchos_Assert.hpp"
#include <Teuchos_DefaultComm.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Operator.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>

#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosSolverManager.hpp>

#include <algorithm>

class TpetraManger
{
  public:
    // Singleton implementation
    static TpetraManger& theTpetraInitializer(int& argc, char**& argv)
    {
        static TpetraManger TpetraInit(argc, argv);
        return TpetraInit;
    }

  private:
    // Enforce Singleton pattern
    TpetraManger(int& argc, char**& argv)
    {
        Tpetra::initialize(&argc, &argv);
    }

    ~TpetraManger()
    {
        Tpetra::finalize();
    }

  private:
    TpetraManger(const TpetraManger&) = delete;
    TpetraManger& operator=(const TpetraManger&) = delete;
    TpetraManger(TpetraManger&&) = delete;
    TpetraManger& operator=(TpetraManger&&) = delete;
};

template <typename TeuchosCommType>
void runTest(TeuchosCommType teuchosComm, const bool callFinalize = false)
{
  using LO = Tpetra::Vector<>::local_ordinal_type;
  using GO = Tpetra::Vector<>::global_ordinal_type;
  using MV = Tpetra::MultiVector<double, LO, GO>;
  using OP = Tpetra::Operator<double, LO, GO>;

  auto fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  TEUCHOS_ASSERT_EQUALITY(teuchosComm->getSize(), 1);

  //  setup 4x4 linear system
  int numRows = 4;
  int indexBase = 0;

  std::vector<GO> myGlobalNodes;
  myGlobalNodes.push_back(0);
  myGlobalNodes.push_back(1);

  auto map = Teuchos::rcp( new Tpetra::Map<LO, GO>(numRows, Teuchos::arrayViewFromVector(myGlobalNodes), indexBase, teuchosComm));
  auto A = Teuchos::rcp( new Tpetra::CrsMatrix<double, LO, GO>(map, numRows) );
  auto x = Teuchos::rcp( new Tpetra::Vector<double, LO, GO>(map) );
  auto b = Teuchos::rcp( new Tpetra::Vector<double, LO, GO>(map) );

  std::vector<GO> cols;
  std::vector<double> vals;
  cols.clear();
  vals.clear();
  cols.push_back(0); cols.push_back(1);
  vals.push_back(2); vals.push_back(-1);
  A->insertGlobalValues(0, Teuchos::arrayViewFromVector(cols), Teuchos::arrayViewFromVector(vals));
  cols.clear();
  vals.clear();
  cols.push_back(0); cols.push_back(1);
  vals.push_back(-1); vals.push_back(4);
  A->insertGlobalValues(1, Teuchos::arrayViewFromVector(cols), Teuchos::arrayViewFromVector(vals));
  b->replaceGlobalValue(0, 1.);
  b->replaceGlobalValue(1, 2.);
  A->fillComplete();

  // A->describe(*fos, Teuchos::VERB_EXTREME);
  // b->describe(*fos, Teuchos::VERB_EXTREME);

  auto linearProblem = Teuchos::rcp( new Belos::LinearProblem<double, MV, OP>( A, x, b ) );

  // NOTE: Comment out everything below here and all variants of the test pass
  linearProblem->setProblem();

  Belos::SolverFactory<double, MV, OP> solverFactory;

  auto belosList = Teuchos::rcp( new Teuchos::ParameterList );
  belosList->set( "Maximum Iterations", numRows);
  belosList->set( "Convergence Tolerance", 1.e-12 );
  // belosList->set( "Verbosity", Belos::IterationDetails + Belos::OrthoDetails
  //                            + Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails);

  auto solverManager = solverFactory.create("CG", belosList);
  solverManager->setProblem(linearProblem);
  solverManager->solve();

  if(callFinalize){
    Tpetra::finalize();
  }
}

// Call Tpetra::finalize() right before end of main
void finalizeAtEnd(int argc, char *argv[]){
  Tpetra::initialize(&argc, &argv);
  auto teuchosComm = Teuchos::DefaultComm<int>::getComm();

  runTest(teuchosComm);

  Tpetra::finalize();
}

// Use singleton to call Tpetra::finalize() automatically
void singleton(int argc, char *argv[]){
  [[maybe_unused]] const TpetraManger& tpetraInit(TpetraManger::theTpetraInitializer(argc, argv));
  auto teuchosComm = Teuchos::DefaultComm<int>::getComm();

  runTest(teuchosComm);
}

// Call Tpetra::finalize() within scope of runTest
void finalizeInTestScope(int argc, char *argv[]){
  Tpetra::initialize(&argc, &argv);
  auto teuchosComm = Teuchos::DefaultComm<int>::getComm();

  runTest(teuchosComm, true);

  Tpetra::finalize();
}


int main(int argc, char *argv[])
{

  Teuchos::CommandLineProcessor clp(true);
  int test_case = 0;
  clp.setOption("case",&test_case);

  const Teuchos::CommandLineProcessor::EParseCommandLineReturn parseResult = clp.parse (argc, argv);
  switch (parseResult) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
  case Teuchos::CommandLineProcessor::PARSE_ERROR:
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
  }

  std::cout << "test case: " << test_case << std::endl;

  switch (test_case) {
  case 0:
    // This works
    finalizeAtEnd(argc, argv);
    break;

  case 1:
    /* This used to fail with the following error
       //   Kokkos::Impl::SharedAllocationRecord 'p6�' failed decrement count = 0
       //   Kokkos::Impl::SharedAllocationRecord failed decrement count
     */
    singleton(argc, argv);
    break;

  case 2:
    /*
      This fails. And that's the correct thing to do since the solver manager and its memory
      are deallocated after Kokkos::finalize.
     */
    finalizeInTestScope(argc, argv);
    break;

  default:
    return EXIT_FAILURE;
  }

  return 0;
}
