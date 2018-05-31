#include "Ifpack2_Factory.hpp"
#include "Ifpack2_Details_CanChangeMatrix.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosSolverFactory.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Core.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
//#include "Teuchos_ParameterXMLFileReader.hpp"

#include <algorithm> // std::transform
#include <cctype> // std::toupper
#include <sstream>
#include <functional>

namespace { // (anonymous)

// See example here:
//
// http://en.cppreference.com/w/cpp/string/byte/toupper
std::string stringToUpper (std::string s)
{
  std::transform (s.begin (), s.end (), s.begin (),
                  [] (unsigned char c) { return std::toupper (c); });
  return s;
}

std::vector<std::string>
splitIntoStrings (const std::string& s,
                  const char sep = ',')
{
  using size_type = std::string::size_type;

  size_type cur_pos;
  size_type last_pos = 0;
  size_type length = s.length ();

  std::vector<std::string> strings;
  while (last_pos < length + size_type (1)) {
    cur_pos = s.find_first_of(sep, last_pos);
    if (cur_pos == std::string::npos) {
      cur_pos = length;
    }
    if (cur_pos != last_pos) {
      auto token = std::string (s.data () + last_pos,
                                static_cast<size_type> (cur_pos - last_pos));
      strings.push_back (stringToUpper (token));
    }
    last_pos = cur_pos + size_type (1);
  }
  return strings;
}

template<class T>
std::vector<T>
splitIntoValues (const std::string& s,
                 const char sep = ',')
{
  using size_type = std::string::size_type;

  size_type cur_pos;
  size_type last_pos = 0;
  size_type length = s.length ();

  std::vector<T> values;
  while (last_pos < length + size_type (1)) {
    cur_pos = s.find_first_of(sep, last_pos);
    if (cur_pos == std::string::npos) {
      cur_pos = length;
    }
    if (cur_pos != last_pos) {
      auto token = std::string (s.data () + last_pos,
                                static_cast<size_type> (cur_pos - last_pos));
      T val {};
      std::istringstream is (token);
      is >> val;
      if (is) {
        values.push_back (val);
      }
    }
    last_pos = cur_pos + size_type (1);
  }
  return values;
}

// Values of command-line arguments.
struct CmdLineArgs {
  std::string matrixFilename;
  std::string rhsFilename;
  std::string solverTypes = "GMRES";
  std::string orthogonalizationMethod = "ICGS";
  std::string convergenceToleranceValues = "1.0e-2";
  std::string maxIterValues = "100";
  std::string restartLengthValues = "20";
  std::string preconditionerTypes = "RELAXATION";
};

// Read in values of command-line arguments.
bool
getCmdLineArgs (CmdLineArgs& args, int argc, char* argv[])
{
  Teuchos::CommandLineProcessor cmdp (false, true);
  cmdp.setOption ("matrixFilename", &args.matrixFilename, "Name of Matrix "
                  "Market file with the sparse matrix A");
  cmdp.setOption ("rhsFilename", &args.rhsFilename, "Name of Matrix Market "
                  "file with the right-hand side vector(s) B");
  cmdp.setOption ("solverTypes", &args.solverTypes,
                  "One or more Belos solver types, "
                  "separated by commas");
  cmdp.setOption ("convergenceTolerances", &args.convergenceToleranceValues,
                  "One or more doubles, separated by commas; each value "
                  "is a convergence tolerance to try");
  cmdp.setOption ("orthogonalizationMethod", &args.orthogonalizationMethod,
                  "Orthogonalization method (for GMRES solver only; "
                  "ignored otherwise)");
  cmdp.setOption ("maxIters", &args.maxIterValues,
                  "One or more integers, separated by commas; each value "
                  "is a maximum number of solver iterations to try");
  cmdp.setOption ("restartLengths", &args.restartLengthValues,
                  "One or more integers, separated by commas; each value "
                  "is a maximum restart length to try (for GMRES solver only; "
                  "ignored otherwise)");
  cmdp.setOption ("preconditionerTypes", &args.preconditionerTypes,
                  "One or more Ifpack2 preconditioner types, "
                  "separated by commas");

  auto result = cmdp.parse (argc, argv);
  return result == Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL;
}

template<class ScalarType>
struct BelosSolverResult {
  typename Teuchos::ScalarTraits<ScalarType>::magnitudeType achievedTolerance;
  int numIters;
  bool converged;
};

template<class CrsMatrixType>
class BelosIfpack2Solver {
private:
  using scalar_type = typename CrsMatrixType::scalar_type;
  using local_ordinal_type = typename CrsMatrixType::local_ordinal_type;
  using global_ordinal_type = typename CrsMatrixType::global_ordinal_type;
  using node_type = typename CrsMatrixType::node_type;
  using row_matrix_type =
    Tpetra::RowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;
  using MV = Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;
  using OP = Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;
  using problem_type = Belos::LinearProblem<scalar_type, MV, OP>;
  using solver_type = Belos::SolverManager<scalar_type, MV, OP>;
  using preconditioner_type =
    Ifpack2::Preconditioner<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;

  void createPreconditioner ()
  {
    rightPrec_ = Teuchos::null; // destroy old one first, to reduce peak memory use
    rightPrec_ =
      Ifpack2::Factory::create<row_matrix_type> (precType_, A_);
    rightPrec_->setParameters (* (precParams_));
  }

  void createSolver ()
  {
    Belos::SolverFactory<scalar_type, MV, OP> belosFactory;
    solver_ = Teuchos::null; // destroy old one first, to reduce peak memory use
    solver_ = belosFactory.create (solverType_, solverParams_);
  }

  void setPreconditionerMatrix (const Teuchos::RCP<const row_matrix_type>& A)
  {
    if (rightPrec_.get () != nullptr) {
      // Not all Ifpack2 preconditioners have a setMatrix method.  Use
      // it if it exists; otherwise, start over with a new instance.
      using can_change_matrix = Ifpack2::Details::CanChangeMatrix<row_matrix_type>;
      can_change_matrix* rightPrec = dynamic_cast<can_change_matrix*> (rightPrec_.get ());
      if (rightPrec != nullptr) {
        rightPrec_->setMatrix (A);
      }
      else {
        rightPrec_ = Teuchos::null; // blow it away; make a new one only on demand
      }
    }
  }

  void initializePreconditioner ()
  {
    if (rightPrec_.get () == nullptr) {
      createPreconditioner ();
    }
    rightPrec_->initialize ();
  }

  void computePreconditioner ()
  {
    if (rightPrec_.get () == nullptr) {
      createPreconditioner ();
      rightPrec_->initialize ();
    }
    rightPrec_->compute ();
  }

public:
  BelosIfpack2Solver () = default;

  BelosIfpack2Solver (const Teuchos::RCP<const CrsMatrixType>& A,
                      const std::string& solverType = "GMRES",
                      const std::string& precType = "RELAXATION") :
    A_ (A),
    solverType_ (solverType),
    precType_ (precType)
  {}

  void setMatrix (const Teuchos::RCP<const CrsMatrixType>& A)
  {
    if (A_.get () != A.get ()) {
      setPreconditionerMatrix (A);
      // Belos solvers don't deal well with a complete change of the matrix.
      solver_ = Teuchos::null;
    }
    A_ = A;
  }

  void
  setPreconditionerTypeAndParameters (const std::string& precType,
                                      const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    precType_ = precType;
    precParams_ = params;
    if (rightPrec_.get () != nullptr) {
      if (precType_ != precType) {
        rightPrec_ = Teuchos::null; // blow it away; make a new one only on demand
      }
      else {
        rightPrec_->setParameters (*params);
      }
    }
  }

  void
  setSolverTypeAndParameters (const std::string& solverType,
                              Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    solverType_ = solverType;
    solverParams_ = params;
    if (solver_.get () != nullptr) {
      if (solverType_ != solverType) {
        solver_ = Teuchos::null; // blow it away; make a new one only on demand
      }
      else {
        solver_->setParameters (params);
      }
    }
  }

  void initialize ()
  {
    TEUCHOS_TEST_FOR_EXCEPTION
      (A_.get () == nullptr, std::runtime_error, "Solver: You must call "
       "setMatrix with a nonnull matrix before you may call initialize.");
    // Calling this implies that the matrix's graph has changed.
    // Belos' solvers don't handle that very well, so best practice is
    // to recreate them in this case.
    solver_ = Teuchos::null;
    initializePreconditioner ();
  }

  void compute ()
  {
    TEUCHOS_TEST_FOR_EXCEPTION
      (A_.get () == nullptr, std::runtime_error, "Solver: You must call "
       "setMatrix with a nonnull matrix before you may call compute.");
    computePreconditioner ();
  }


  BelosSolverResult<scalar_type>
  solve (Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
         const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& B)
  {
    using Teuchos::RCP;
    using Teuchos::rcpFromRef;

    if (solver_.get () == nullptr) {
      createSolver ();
    }
    if (rightPrec_.get () == nullptr) {
      createPreconditioner ();
      initializePreconditioner ();
      computePreconditioner ();
    }

    RCP<problem_type> problem (new problem_type (A_, rcpFromRef (X), rcpFromRef (B)));
    if (rightPrec_.get () != nullptr) {
      problem->setRightPrec (rightPrec_);
    }
    problem->setProblem ();
    solver_->setProblem (problem);
    const Belos::ReturnType solveResult = solver_->solve ();

    typename Teuchos::ScalarTraits<scalar_type>::magnitudeType tol {-1.0};
    try {
      tol = solver_->achievedTol ();
    }
    catch (...) {}
    const int numIters = solver_->getNumIters ();
    const bool converged = (solveResult == Belos::Converged);
    return {tol, numIters, converged};
  }

private:
  Teuchos::RCP<const CrsMatrixType> A_;
  Teuchos::RCP<solver_type> solver_;
  Teuchos::RCP<preconditioner_type> rightPrec_;

  std::string solverType_;
  Teuchos::RCP<Teuchos::ParameterList> solverParams_;
  std::string precType_;
  Teuchos::RCP<Teuchos::ParameterList> precParams_;
};

class TpetraInstance {
public:
  TpetraInstance (int* argc, char*** argv) {
    Tpetra::initialize (argc, argv);
  }

  ~TpetraInstance () {
    Tpetra::finalize ();
  }
};

template<class CrsMatrixType, class MultiVectorType>
void
solveAndReport (BelosIfpack2Solver<CrsMatrixType>& solver,
                MultiVectorType& X,
                const MultiVectorType& B,
                const int myRank,
                const std::string& solverType,
                const std::string& precType,
                const std::string& orthogonalizationMethod,
                const typename MultiVectorType::mag_type convergenceTolerance,
                const int restartLength,
                const int maxIters)
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;

  RCP<ParameterList> solverParams (new ParameterList ("Belos"));
  RCP<ParameterList> precParams (new ParameterList ("Ifpack2"));

  solverParams->set ("Convergence Tolerance", convergenceTolerance);
  solverParams->set ("Maximum Iterations", maxIters);
  if (solverType == "GMRES") {
    solverParams->set ("Num Blocks", restartLength);
    solverParams->set ("Maximum Restarts", restartLength * maxIters);
    solverParams->set ("Orthogonalization", orthogonalizationMethod);
  }
  solver.setSolverTypeAndParameters (solverType, solverParams);
  solver.setPreconditionerTypeAndParameters (precType, precParams);
  solver.initialize ();
  solver.compute ();

  // Solve the linear system AX=B.
  auto result = solver.solve (X, B);

  if (myRank == 0) {
    using std::cout;
    using std::endl;
    cout << "Solver:" << endl
         << "  Solver type: " << solverType << endl
         << "  Preconditioner type: " << precType << endl
         << "  Convergence tolerance: " << convergenceTolerance << endl
         << "  Maximum number of iterations: " << maxIters << endl;
    if (solverType == "GMRES") {
      cout << "  Restart length: " << restartLength << endl
           << "  Orthogonalization method: " << orthogonalizationMethod << endl;
    }
    cout << "Results:" << endl
         << "  Converged: " << (result.converged ? "true" : "false") << endl
         << "  Number of iterations: " << result.numIters << endl
         << "  Achieved tolerance: " << result.achievedTolerance << endl
         << endl;
  }
}

} // namespace (anonymous)

int
main (int argc, char* argv[])
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using std::cerr;
  using std::endl;
  using crs_matrix_type = Tpetra::CrsMatrix<>;
  using MV = Tpetra::MultiVector<>;
  // using mag_type = MV::mag_type;
  using reader_type = Tpetra::MatrixMarket::Reader<crs_matrix_type>;

  TpetraInstance tpetraInstance (&argc, &argv);
  auto comm = Tpetra::getDefaultComm ();

  // Get command-line arguments.
  CmdLineArgs args;
  const bool gotCmdLineArgs = getCmdLineArgs (args, argc, argv);
  if (! gotCmdLineArgs) {
    if (comm->getRank () == 0) {
      cerr << "Failed to get command-line arguments!" << endl;
    }
    return EXIT_FAILURE;
  }

  if (args.matrixFilename == "") {
    if (comm->getRank () == 0) {
      cerr << "Must specify sparse matrix filename!" << endl;
    }
    return EXIT_FAILURE;
  }

  std::vector<std::string> solverTypes;
  if (args.solverTypes == "") {
    solverTypes = {"GMRES", "TFQMR", "BICGSTAB"};
  }
  else {
    solverTypes = splitIntoStrings (args.solverTypes);
  }

  std::vector<std::string> preconditionerTypes;
  if (args.preconditionerTypes == "") {
    preconditionerTypes = {"RELAXATION"};
  }
  else {
    preconditionerTypes = splitIntoStrings (args.preconditionerTypes);
  }

  std::vector<int> maxIterValues;
  if (args.maxIterValues == "") {
    maxIterValues = {100};
  }
  else {
    maxIterValues = splitIntoValues<int> (args.maxIterValues);
  }

  std::vector<int> restartLengthValues;
  if (args.restartLengthValues == "") {
    restartLengthValues = {20};
  }
  else {
    restartLengthValues = splitIntoValues<int> (args.restartLengthValues);
  }

  std::vector<double> convergenceToleranceValues;
  if (args.convergenceToleranceValues == "") {
    convergenceToleranceValues = {20};
  }
  else {
    convergenceToleranceValues =
      splitIntoValues<double> (args.convergenceToleranceValues);
  }

  // Read sparse matrix A from Matrix Market file.
  RCP<crs_matrix_type> A =
    reader_type::readSparseFile (args.matrixFilename, comm);
  if (A.get () == nullptr) {
    if (comm->getRank () == 0) {
      cerr << "Failed to load sparse matrix A from file "
        "\"" << args.matrixFilename << "\"!" << endl;
    }
    return EXIT_FAILURE;
  }

  // Read right-hand side vector(s) B from Matrix Market file, or
  // generate B if file not specified.
  RCP<MV> B;
  RCP<MV> X;
  if (args.rhsFilename == "") {
    B = Teuchos::rcp (new MV (A->getRangeMap (), 1));
    X = Teuchos::rcp (new MV (A->getDomainMap (), 1));
    X->randomize ();
    A->apply (*X, *B);
    X->putScalar (0.0);
  }
  else {
    auto map = A->getRangeMap ();
    B = reader_type::readDenseFile (args.rhsFilename, comm, map);
    if (B.get () == nullptr) {
      if (comm->getRank () == 0) {
        cerr << "Failed to load right-hand side vector(s) from file \""
             << args.rhsFilename << "\"!" << endl;
      }
      return EXIT_FAILURE;
    }
    X = Teuchos::rcp (new MV (A->getDomainMap (), B->getNumVectors ()));
  }

  // Create the solver.
  BelosIfpack2Solver<crs_matrix_type> solver (A);

  // Solve the linear system using various solvers and preconditioners.
  for (std::string solverType : solverTypes) {
    for (std::string precType : preconditionerTypes) {
      for (int maxIters : maxIterValues) {
        for (int restartLength : restartLengthValues) {
          for (double convTol : convergenceToleranceValues) {
            solveAndReport (solver, *X, *B, comm->getRank (),
                            solverType,
                            precType,
                            args.orthogonalizationMethod,
                            convTol,
                            restartLength,
                            maxIters);
          }
        }
      }
    }
  }

  return EXIT_SUCCESS;
}
