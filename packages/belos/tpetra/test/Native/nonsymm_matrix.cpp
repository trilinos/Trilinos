#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosTpetraAdapter.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Vector.hpp"
#include "KokkosBlas1_mult.hpp"
#include "Teuchos_ParameterList.hpp"
#include <iostream>

namespace { // (anonymous)

struct CommandLineOptions {
  std::string solverName {"TPETRA GMRES"};
  double offDiagDiff = 1.0 / 8.0;
  // mfh 14 Aug 2018: GMRES takes 20 iterations on this problem (with
  // offDiagDiff = 1/8).  We add 10 iterations to allow for rounding
  // error and differences in the algorithm.
  int maxAllowedNumIters {30};
  int maxNumIters {100};
  int restartLength {30};
  int stepSize {1};
  bool useCholQR {false};
  bool useCholQR2 {false};
  bool computeRitzValues {true};
  bool zeroInitialGuess {true};
  bool verbose {true};
};
CommandLineOptions commandLineOptions;

TEUCHOS_STATIC_SETUP()
{
  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions (true);
  clp.setOption ("solver", &commandLineOptions.solverName,
                 "Name of the solver to test.  Belos::SolverFactory::create "
                 "must accept this string.  Protect with double quotes if it "
                 "has spaces: e.g., \"TPETRA CG PIPELINE\".");
  clp.setOption ("offDiagDiff", &commandLineOptions.offDiagDiff,
                 "Value of the term that makes the matrix nonsymmetric");
  clp.setOption ("maxNumAllowedIters",
                 &commandLineOptions.maxAllowedNumIters,
                 "Maximum number of iterations that the solver is "
                 "allowed to take before converging, in order for "
                 "the test to pass.");
  clp.setOption ("maxNumIters", &commandLineOptions.maxNumIters,
                 "Maximum number of iterations that the solver is "
                 "allowed to take, over all restart cycles.  This has "
                 "nothing to do with the test passing.");
  clp.setOption ("restartLength", &commandLineOptions.restartLength,
                 "Maximum number of iterations per restart cycle.  "
                 "This corresponds to the standard Belos parameter "
                 "\"Num Blocks\".");
  clp.setOption ("stepSize", &commandLineOptions.stepSize,
                 "Step size; only applies to algorithms that take it.");
  clp.setOption ("useCholQR", "noCholQR", &commandLineOptions.useCholQR,
                 "Whether to use CholQR");
  clp.setOption ("useCholQR2", "noCholQR2", &commandLineOptions.useCholQR2,
                 "Whether to use CholQR2");
  clp.setOption ("computeRitzValues", "noRitzValues", &commandLineOptions.computeRitzValues,
                 "Whether to compute Ritz values");
  clp.setOption ("zeroInitialGuess", "nonzeroInitialGuess",
                 &commandLineOptions.zeroInitialGuess, "Whether to test "
                 "with a zero, or a nonzero, initial guess vector");
  clp.setOption ("verbose", "quiet", &commandLineOptions.verbose,
                 "Whether to print verbose output");
}

// Create a nonsymmetric tridiagonal matrix representing a
// discretization of a 1-D convection-diffusion operator.
// Stencil looks like this:
//
// [1/4 - offDiagDiff, 1, 1/4 + offDiagDiff]
//
// The point is to have the matrix be diagonally dominant, but still
// nonsymmetric.
template<class SC>
Teuchos::RCP<Tpetra::CrsMatrix<SC> >
createNonsymmTridiagMatrix (const Teuchos::RCP<const Tpetra::Map<> >& rowMap,
                            const SC offDiagDiff)
{
  using Teuchos::rcp;
  using crs_matrix_type = Tpetra::CrsMatrix<SC>;
  using map_type = Tpetra::Map<>;
  using LO = typename map_type::local_ordinal_type;
  using GO = typename map_type::global_ordinal_type;
  using crs_matrix_type = Tpetra::CrsMatrix<SC>;
  using STS = Teuchos::ScalarTraits<SC>;
  //using mag_type = typename Tpetra::CrsMatrix<SC>::mag_type;

  const LO lclNumRows = rowMap.is_null () ? LO (0) :
    LO (rowMap->getNodeNumElements ());
  const GO gblMinGblInd = rowMap->getMinAllGlobalIndex ();
  const GO gblMaxGblInd = rowMap->getMaxAllGlobalIndex ();
  auto A = rcp (new crs_matrix_type (rowMap, 3, Tpetra::StaticProfile));

  const SC ONE = STS::one ();
  const SC TWO = ONE + ONE;
  const SC FOUR = TWO + TWO;
  const SC baseOffDiagEnt = ONE / FOUR;
  const SC subDiagEnt = baseOffDiagEnt - offDiagDiff;
  const SC diagEnt = ONE;
  const SC superDiagEnt = baseOffDiagEnt + offDiagDiff;

  Teuchos::Array<GO> gblColIndsBuf (3);
  Teuchos::Array<SC> valsBuf (3);
  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    const GO gblRow = rowMap->getGlobalElement (lclRow);
    const GO gblCol = gblRow;
    LO numEnt = 0; // to be set below
    if (gblRow == gblMinGblInd && gblRow == gblMaxGblInd) {
      numEnt = 1;
      valsBuf[0] = diagEnt;
      gblColIndsBuf[0] = gblCol;
    }
    else if (gblRow == gblMinGblInd) {
      numEnt = 2;
      valsBuf[0] = diagEnt;
      valsBuf[1] = superDiagEnt;
      gblColIndsBuf[0] = gblCol;
      gblColIndsBuf[1] = gblCol + GO (1);
    }
    else if (gblRow == gblMaxGblInd) {
      numEnt = 2;
      valsBuf[0] = subDiagEnt;
      valsBuf[1] = diagEnt;
      gblColIndsBuf[0] = gblCol - GO (1);
      gblColIndsBuf[1] = gblCol;
    }
    else {
      numEnt = 3;
      valsBuf[0] = subDiagEnt;
      valsBuf[1] = diagEnt;
      valsBuf[2] = superDiagEnt;
      gblColIndsBuf[0] = gblCol - GO (1);
      gblColIndsBuf[1] = gblCol;
      gblColIndsBuf[2] = gblCol + GO (1);
    }
    Teuchos::ArrayView<GO> gblColInds = gblColIndsBuf.view (0, numEnt);
    Teuchos::ArrayView<SC> vals = valsBuf.view (0, numEnt);
    A->insertGlobalValues (gblRow, gblColInds, vals);
  }
  A->fillComplete ();
  return A;
}

void
testSolver (Teuchos::FancyOStream& out,
            bool& success,
            const std::string& solverName,
            const int maxAllowedNumIters,
            const bool verbose)
{
  using Teuchos::FancyOStream;
  using Teuchos::getFancyOStream;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using std::endl;
  using map_type = Tpetra::Map<>;
  using MV = Tpetra::MultiVector<>;
  using OP = Tpetra::Operator<>;
  using SC = MV::scalar_type;
  using GO = map_type::global_ordinal_type;
  using mag_type = MV::mag_type;
  using STS = Teuchos::ScalarTraits<SC>;
  using STM = Teuchos::ScalarTraits<mag_type>;
  // The Teuchos unit test framework likes to capture output to 'out',
  // and not print anything until the test is done.  This can hinder
  // debugging.  If the test crashes without useful output, try
  // setting this to 'true'.  That will change 'myOut' from an alias
  // to 'out', into a wrapper for std::cerr.
  constexpr bool debug = false;

  const SC ZERO = STS::zero ();
  const SC ONE = STS::one ();

  RCP<FancyOStream> myOutPtr =
    debug ? getFancyOStream (rcpFromRef (std::cerr)) : rcpFromRef (out);
  Teuchos::FancyOStream& myOut = *myOutPtr;

  myOut << "Test \"native\" Tpetra version of solver \"" << solverName << "\""
        << endl;
  myOut << (STS::isComplex ? "complex type" : "real type") << endl;
  Teuchos::OSTab tab1 (out);

  myOut << "Create the linear system to solve" << endl;

  auto comm = Tpetra::getDefaultComm ();
  const GO gblNumRows = 10000;
  const GO indexBase = 0;
  RCP<const map_type> map (new map_type (gblNumRows, indexBase, comm));
  auto A = createNonsymmTridiagMatrix (map, commandLineOptions.offDiagDiff);

  MV X_initial (A->getDomainMap (), 1);
  if (commandLineOptions.zeroInitialGuess) {
    X_initial.putScalar (ZERO); // (re)set initial guess to zero
  }
  else {
    X_initial.putScalar (ONE); // just something nonzero, to test
  }
  MV X (X_initial, Teuchos::Copy);
  MV B (A->getRangeMap (), 1);
  B.randomize ();

  myOut << "Create solver instance using Belos::SolverFactory" << endl;

  RCP<Belos::SolverManager<SC, MV, OP> > solver;
  try {
    Belos::SolverFactory<SC, MV, OP> factory;
    solver = factory.create (solverName, Teuchos::null);
  }
  catch (std::exception& e) {
    myOut << "*** FAILED: Belos::SolverFactory::create threw an exception: "
        << e.what () << endl;
    success = false;
    return;
  }

  TEST_ASSERT( solver.get () != nullptr );
  if (solver.get () == nullptr) {
    myOut << "Belos::SolverFactory returned a null solver." << endl;
    return;
  }

  myOut << "Set parameters" << endl;
  RCP<ParameterList> params = parameterList ("Belos");
  params->set ("Verbosity", verbose ? 1 : 0);
  params->set ("Maximum Iterations", commandLineOptions.maxNumIters);
  params->set ("Num Blocks", commandLineOptions.restartLength);
  if (solverName == "TPETRA GMRES S-STEP") {
    params->set ("Step Size", commandLineOptions.stepSize);
    params->set ("Compute Ritz Values", commandLineOptions.computeRitzValues);
    params->set ("CholeskyQR",  commandLineOptions.useCholQR);
    params->set ("CholeskyQR2", commandLineOptions.useCholQR2);
  }
  try {
    solver->setParameters (params);
  }
  catch (std::exception& e) {
    myOut << "*** FAILED: setParameters threw an exception: "
        << e.what () << endl;
    success = false;
    return;
  }
  catch (...) {
    myOut << "*** FAILED: setParameters threw an exception "
      "not a subclass of std::exception." << endl;
    success = false;
    return;
  }

  myOut << "Set up the linear system to solve" << endl;
  auto lp = rcp (new Belos::LinearProblem<SC, MV, OP> (A, rcpFromRef (X),
                                                       rcpFromRef (B)));
  lp->setProblem ();

  myOut << "Solve the linear system" << endl;
  solver->setProblem (lp);
  const Belos::ReturnType belosResult = solver->solve ();

  myOut << "Belos solver wrapper result: "
        << (belosResult == Belos::Converged ? "Converged" : "Unconverged")
        << endl
        << "Number of iterations: " << solver->getNumIters ()
        << endl;

  TEST_ASSERT( solver->getNumIters () <= maxAllowedNumIters );
  TEST_ASSERT( belosResult == Belos::Converged );

  myOut << "Check the explicit residual norm(s)" << endl;

  // Get the tolerance that the solver actually used.
  const mag_type tol = [&] () {
      const char tolParamName[] = "Convergence Tolerance";
      auto pl = solver->getCurrentParameters ();
      if (! pl->isType<mag_type> (tolParamName)) {
        pl = solver->getValidParameters ();
      }
      TEUCHOS_TEST_FOR_EXCEPTION
        (! pl->isType<mag_type> (tolParamName), std::logic_error,
       "Solver lacks \"" << tolParamName << "\" parameter, in either "
         "getCurrentParameters() or getValidParameters().");
      return pl->get<mag_type> (tolParamName);
    } ();

  MV R_initial (B.getMap (), B.getNumVectors ());
  A->apply (X_initial, R_initial);
  R_initial.update (ONE, B, -ONE);
  Teuchos::Array<mag_type> R_initial_norms (R_initial.getNumVectors ());
  R_initial.norm2 (R_initial_norms ());

  MV R_final (B.getMap (), B.getNumVectors ());
  A->apply (X, R_final);
  R_final.update (ONE, B, -ONE);
  Teuchos::Array<mag_type> R_final_norms (R_final.getNumVectors ());
  R_final.norm2 (R_final_norms ());

  for (size_t j = 0; j < R_final.getNumVectors (); ++j) {
    const mag_type relResNorm = (R_initial_norms[j] == STM::zero ()) ?
      R_final_norms[j] :
      R_final_norms[j] / R_initial_norms[j];
    myOut << "Column " << (j+1) << " of " << R_final.getNumVectors ()
          << ": Absolute residual norm: " << R_final_norms[j]
          << ", Relative residual norm: " << relResNorm
          << ", Tolerance: " << tol
          << endl;
    TEST_ASSERT( relResNorm <= tol );
  }
  myOut << endl;

  Teuchos::TimeMonitor::summarize();
}

TEUCHOS_UNIT_TEST( TpetraNativeSolvers, Diagonal )
{
  testSolver (out, success, commandLineOptions.solverName,
              commandLineOptions.maxAllowedNumIters,
              commandLineOptions.verbose);
}

} // namespace (anonymous)

int main (int argc, char* argv[])
{
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain (argc, argv);
}
