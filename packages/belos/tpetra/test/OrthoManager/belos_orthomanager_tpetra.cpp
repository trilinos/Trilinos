// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file belos_orthomanager_tpetra.cpp
/// \brief Test (Mat)OrthoManager subclass(es) with Tpetra
///
/// Test various subclasses of (Mat)OrthoManager, using
/// Tpetra::MultiVector as the multivector implementation,
/// and Tpetra::Operator as the operator implementation.
///
#include "belos_orthomanager_tpetra_util.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Tpetra_Core.hpp"

using std::endl;

//////////////////////////////////////////////////////////////////////
// Command-line arguments
//////////////////////////////////////////////////////////////////////

// The name of the (Mat)OrthoManager subclass to instantiate.
std::string orthoManName ("DGKS");

// For SimpleOrthoManager: the normalization method to use.
std::string normalization ("CGS");

// Name of the Harwell-Boeing sparse matrix file from which to read
// the inner product operator matrix.  If name is "" or not provided
// at the command line, use the standard Euclidean inner product.
std::string filename;

bool verbose = false;
bool debug = false;

// The OrthoManager is tested with three different multivectors: S
// (sizeS columns), X1 (sizeX1 columns), and X2 (sizeX2 columns).  The
// values below are defaults and may be changed by the corresponding
// command-line arguments.
int sizeS  = 5;
int sizeX1 = 11;
int sizeX2 = 13;

// Default _global_ number of rows.  The number of rows per MPI
// process must be no less than max(sizeS, sizeX1, sizeX2).  To ensure
// that the test always passes with default parameters, we must scale
// below by the number of processes.  The default value below may be
// changed by a command-line parameter with a corresponding name.
int numRows = 100;

// Set default values of command-line arguments,
// and get their actual values from the command line.
static void
getCmdLineArgs (const Teuchos::Comm<int>& comm, int argc, char* argv[])
{
  using Teuchos::CommandLineProcessor;
  typedef Tpetra::MultiVector<>::scalar_type ST;

  // Define a OrthoManagerFactory to use to get the names of valid
  // orthogonalization manager types.  We won't use this factory to
  // create them, so we should be able to pick the Scalar, MV, and
  // OP template parameters as we wish.
  typedef Belos::OrthoManagerFactory<double,
    Tpetra::MultiVector<ST>, Tpetra::Operator<ST> > factory_type;
  factory_type factory;

  ////////////////////////////////////////////////////////////////////
  // Set default values of command-line arguments.
  ////////////////////////////////////////////////////////////////////

  // The name of the (Mat)OrthoManager subclass to instantiate.
  orthoManName = factory.defaultName ();

  // For SimpleOrthoManager: the normalization method to use.  Valid
  // values: "MGS", "CGS".
  normalization = "CGS";

  // Name of the Harwell-Boeing sparse matrix file from which to read
  // the inner product operator matrix.  If name is "" or not provided
  // at the command line, use the standard Euclidean inner product.
  filename = "";

  verbose = false;
  debug = false;

  // The OrthoManager is tested with three different multivectors: S
  // (sizeS columns), X1 (sizeX1 columns), and X2 (sizeX2 columns).
  // The values below are defaults and may be changed by the
  // corresponding command-line arguments.
  sizeS  = 5;
  sizeX1 = 11;
  sizeX2 = 13;

  // Default _global_ number of rows.  The number of rows per MPI
  // process must be no less than max(sizeS, sizeX1, sizeX2).  To
  // ensure that the test always passes with default parameters, we
  // scale by the number of processes.  The default value below may be
  // changed by a command-line parameter with a corresponding name.
  numRows = 100 * comm.getSize ();

  ////////////////////////////////////////////////////////////////////
  // Define command-line arguments and parse them.
  ////////////////////////////////////////////////////////////////////

  CommandLineProcessor cmdp (false, true);
  cmdp.setOption ("verbose", "quiet", &verbose, "Print messages and results.");
  cmdp.setOption ("debug", "nodebug", &debug, "Print debugging information.");
  cmdp.setOption ("filename", &filename,
                  "Filename of a Harwell-Boeing sparse matrix, used as the "
                  "inner product operator by the orthogonalization manager."
                  "  If not provided, no matrix is read and the Euclidean "
                  "inner product is used.  This is only valid for Scalar="
                  "double.");
  {
    std::ostringstream os;
    const int numValid = factory.numOrthoManagers ();
    const bool plural = numValid > 1 || numValid == 0;

    os << "OrthoManager subclass to test.  There ";
    os << (plural ? "are " : "is ") << numValid << (plural ? "s: " : ": ");
    factory.printValidNames (os);
    os << ".";
    cmdp.setOption ("ortho", &orthoManName, os.str ().c_str ());
  }
  cmdp.setOption ("normalization", &normalization,
                  "For SimpleOrthoManager (--ortho=Simple): the normalization "
                  "method to use.  Valid values: \"MGS\", \"CGS\".");
  cmdp.setOption ("numRows", &numRows,
                  "Controls the number of rows of the test "
                  "multivectors.  If an input matrix is given, this "
                  "parameter\'s value is ignored, since the vectors must "
                  "be commensurate with the dimensions of the matrix.");
  cmdp.setOption ("sizeS", &sizeS, "Controls the number of columns of the "
                  "input multivector.");
  cmdp.setOption ("sizeX1", &sizeX1, "Controls the number of columns of the "
                  "first basis.");
  cmdp.setOption ("sizeX2", &sizeX2, "Controls the number of columns of the "
                  "second basis.  We require for simplicity of testing (the "
                  "routines do not require it) that sizeX1 >= sizeX2.");
  // Parse the command-line arguments.
  {
    const CommandLineProcessor::EParseCommandLineReturn parseResult =
      cmdp.parse (argc, argv);
    // If the caller asks us to print the documentation, or does not
    // provide the name of an OrthoManager subclass, we just keep
    // going.  Otherwise, we throw an exception.
    TEUCHOS_TEST_FOR_EXCEPTION(
      parseResult != CommandLineProcessor::PARSE_SUCCESSFUL,
      std::invalid_argument,
      "Failed to parse command-line arguments");
  }
  //
  // Validate command-line arguments
  //
  TEUCHOS_TEST_FOR_EXCEPTION(
    numRows <= 0, std::invalid_argument, "numRows <= 0 is not allowed");
  TEUCHOS_TEST_FOR_EXCEPTION(
    numRows <= sizeS + sizeX1 + sizeX2, std::invalid_argument,
    "numRows <= sizeS + sizeX1 + sizeX2 is not allowed");
}

// C++03 does not allow partial specialization of a template function,
// so we have to declare a class to load the sparse matrix and create
// the Map.

// We don't have a generic Harwell-Boeing sparse file reader yet.  If
// filename == "", then the user doesn't want to read in a sparse
// matrix, so we can just create the appropriate Map and be done with
// it.  Otherwise, raise an exception at run time.
template<class ScalarType, class LocalOrdinalType, class GlobalOrdinalType, class NodeType>
class SparseMatrixLoader {
public:
  typedef Tpetra::Map<LocalOrdinalType, GlobalOrdinalType, NodeType> map_type;
  typedef Tpetra::CrsMatrix<ScalarType, LocalOrdinalType, GlobalOrdinalType, NodeType> matrix_type;

  static void
  load (Teuchos::RCP<map_type>& map,
        Teuchos::RCP<matrix_type>& M,
        const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
        std::ostream& debugOut)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      filename != "", std::logic_error, "Sorry, reading in a Harwell-Boeing "
      "sparse matrix file for ScalarType="
      << Teuchos::TypeNameTraits<ScalarType>::name () << " is not yet "
      "implemented.  This currently only works for ScalarType=double.");
    const GlobalOrdinalType indexBase = 0;
    map = Teuchos::rcp (new map_type (numRows, indexBase, comm, Tpetra::GloballyDistributed));
    M = Teuchos::null;
  }
};

// We _do_ know how to read a Harwell-Boeing sparse matrix file with Scalar=double.
template<class LocalOrdinalType, class GlobalOrdinalType, class NodeType>
class SparseMatrixLoader<double, LocalOrdinalType, GlobalOrdinalType, NodeType> {
public:
  typedef Tpetra::Map<LocalOrdinalType, GlobalOrdinalType, NodeType> map_type;
  typedef Tpetra::CrsMatrix<double, LocalOrdinalType, GlobalOrdinalType, NodeType> matrix_type;

  static void
  load (Teuchos::RCP<map_type>& map,
        Teuchos::RCP<matrix_type>& M,
        const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
        std::ostream& debugOut)
  {
    // If the sparse matrix is loaded successfully, this call will
    // modify numRows to be the number of rows in the sparse matrix.
    // Otherwise, it will leave numRows alone.
    std::pair<Teuchos::RCP<map_type>, Teuchos::RCP<matrix_type> > results =
      Belos::Test::loadSparseMatrix<typename matrix_type::scalar_type,LocalOrdinalType,
      GlobalOrdinalType,
      NodeType> (comm, filename, numRows, debugOut);
    map = results.first;
    M = results.second;
  }
};


// Run the actual test.  The test has the same template parameters as
// MultiVector.  The most interesting template parameters for this
// test are ScalarType and NodeType.
//
// Return true if test passed, else return false.
template<class ScalarType, class LocalOrdinalType, class GlobalOrdinalType, class NodeType>
bool runTest (const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
{
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef ScalarType scalar_type;
  typedef LocalOrdinalType local_ordinal_type;
  typedef GlobalOrdinalType global_ordinal_type;
  typedef NodeType node_type;

  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> MV;
  typedef Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type, node_type> OP;
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;
  typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> crs_matrix_type;

  // Declare an output manager for handling local output.  Initialize,
  // using the caller's desired verbosity level.
  //
  // NOTE In Anasazi, this class is called BasicOutputManager.  In
  // Belos, this class is called OutputManager.  We should eventually
  // resolve the difference.
  RCP<Belos::OutputManager<scalar_type> > outMan =
    Belos::Test::makeOutputManager<scalar_type> (verbose, debug);

  // Stream for debug output.  If debug output is not enabled, then
  // this stream doesn't print anything sent to it (it's a "black
  // hole" stream).
  std::ostream& debugOut = outMan->stream (Belos::Debug);
  Belos::Test::printVersionInfo (debugOut);

  // Load the inner product operator matrix from the given filename.
  // If filename == "", use the identity matrix as the inner product
  // operator (the Euclidean inner product), and leave M as
  // Teuchos::null.  Also return an appropriate Map (which will always
  // be initialized, even if filename == ""; it should never be
  // Teuchos::null).
  RCP<map_type> map;
  RCP<crs_matrix_type> M;
  {
    typedef SparseMatrixLoader<scalar_type, local_ordinal_type,
      global_ordinal_type, node_type> loader_type;
    loader_type::load (map, M, comm, debugOut);
  }
  TEUCHOS_TEST_FOR_EXCEPTION(map.is_null (), std::runtime_error,
    "(Mat)OrthoManager test code failed to initialize the Map.");
  {
    // The maximum number of columns that will be passed to a
    // MatOrthoManager's normalize() routine.  Some MatOrthoManager
    // subclasses (e.g., Tsqr(Mat)OrthoManager) need to have the
    // number of columns no larger than the number of rows on any
    // process.  We check this _after_ attempting to load any sparse
    // matrix to be used as the inner product matrix, because if a
    // sparse matrix is successfully loaded, its number of rows will
    // override the number of rows specified on the command line (if
    // specified), and will also override the default number of rows.
    const size_t maxNormalizeNumCols = std::max (sizeS, std::max (sizeX1, sizeX2));
    // getLocalNumElements() returns a size_t, which is unsigned, and
    // you shouldn't compare signed and unsigned values.  This is why
    // we make maxNormalizeNumCols a size_t as well.
    if (map->getLocalNumElements () < maxNormalizeNumCols) {
      std::ostringstream os;
      os << "The number of elements on this process " << comm->getRank()
         << " is too small for the number of columns that you want to test."
         << "  There are " << map->getLocalNumElements() << " elements on "
        "this process, but the normalize() method of the MatOrthoManager "
        "subclass will need to process a multivector with "
         << maxNormalizeNumCols << " columns.  Not all MatOrthoManager "
        "subclasses can handle a local row block with fewer rows than "
        "columns.";
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str ());
    }
  }

  // This factory object knows how to make a (Mat)OrthoManager
  // subclass, given a name for the subclass and its parameters.
  Belos::OrthoManagerFactory<scalar_type, MV, OP> factory;

  // Using the factory object, instantiate the specified OrthoManager
  // subclass to be tested.  Specify "default" parameters (which
  // should favor accuracy over performance), but override the default
  // parameters to get the desired normalization method for
  // SimpleOrthoManaager.
  RCP<Belos::OrthoManager<scalar_type, MV> > orthoMan;
  {
    std::string label (orthoManName);
    RCP<ParameterList> params =
      parameterList (*(factory.getDefaultParameters (orthoManName)));
    if (orthoManName == "Simple") {
      params->set ("Normalization", normalization);
      label = label + " (" + normalization + " normalization)";
    } else if (orthoManName == "TSQR") {
      params->set ("randomizeNullSpace", false); // for testing the norm of zero-vector is zero
    }
    orthoMan = factory.makeOrthoManager (orthoManName, M, outMan, label, params);
  }

  // Whether the specific OrthoManager subclass promises to compute
  // rank-revealing orthogonalizations.  If yes, then test it on
  // rank-deficient multivectors, otherwise only test it on full-rank
  // multivectors.
  const bool isRankRevealing = factory.isRankRevealing (orthoManName);

  // "Prototype" multivector.  The test code will use this (via
  // Belos::MultiVecTraits) to clone other multivectors as
  // necessary.  (This means the test code doesn't need the Map, and
  // it also makes the test code independent of the idea of a Map.)
  RCP<MV> S = rcp (new MV (map, sizeS));

  // Test the OrthoManager subclass.  Return the number of tests
  // that failed.  None of the tests should fail (this function
  // should return zero).
  int numFailed = 0;
  {
    typedef Belos::Test::OrthoManagerTester<scalar_type, MV> tester_type;
    numFailed = tester_type::runTests (orthoMan, isRankRevealing, S,
                                       sizeX1, sizeX2, outMan);
  }

  if (numFailed != 0) {
    outMan->stream (Belos::Errors) << numFailed << " errors." << endl;
    return false;
  } else {
    return true;
  }
}

int
main (int argc, char *argv[])
{
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  bool success = false;

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  try {
    RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();

    // Get values of command-line arguments.
    getCmdLineArgs (*comm, argc, argv);

    typedef Tpetra::Map<>::local_ordinal_type local_ordinal_type;
    typedef Tpetra::Map<>::global_ordinal_type global_ordinal_type;
    typedef Tpetra::Map<>::node_type node_type;

    {
      typedef Tpetra::MultiVector<>::scalar_type scalar_type;
      success = runTest<scalar_type, local_ordinal_type,
              global_ordinal_type, node_type> (comm);
      if (success) {
        // The Trilinos test framework depends on seeing this message,
        // so don't rely on the OutputManager to report it correctly.
        if (comm->getRank () == 0) {
          std::cout << "End Result: TEST PASSED" << endl;
        }
      }
      else {
        if (comm->getRank () == 0) {
          std::cout << "End Result: TEST FAILED" << endl;
        }
      }
    }
#if defined(HAVE_TEUCHOS_COMPLEX) && defined(HAVE_TPETRA_COMPLEX)
    {
      typedef std::complex<Tpetra::MultiVector<>::scalar_type> scalar_type;
      success = runTest<scalar_type, local_ordinal_type,
              global_ordinal_type, node_type> (comm);
      if (success) {
        // The Trilinos test framework depends on seeing this message,
        // so don't rely on the OutputManager to report it correctly.
        if (comm->getRank () == 0) {
          std::cout << "End Result: TEST PASSED" << endl;
        }
      } else {
        if (comm->getRank () == 0) {
          std::cout << "End Result: TEST FAILED" << endl;
        }
      }
    }
#endif // defined(HAVE_TEUCHOS_COMPLEX) && defined(HAVE_TPETRA_COMPLEX)
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
