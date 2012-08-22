//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

/// \file belos_orthomanager_tpetra.cpp
/// \brief Test (Mat)OrthoManager subclass(es) with Tpetra
///
/// Test various subclasses of (Mat)OrthoManager, using
/// Tpetra::MultiVector as the multivector implementation,
/// and Tpetra::Operator as the operator implementation.
///
#include "belos_orthomanager_tpetra_util.hpp"
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>

using std::endl;

//
// These typedefs make main() as generic as possible.
//
typedef double scalar_type;
typedef int local_ordinal_type;
typedef int global_ordinal_type;

#ifdef HAVE_KOKKOSCLASSIC_TBB
typedef Kokkos::TBBNode node_type;
#else
typedef Kokkos::SerialNode node_type;
#endif // HAVE_KOKKOSCLASSIC_TBB

typedef Teuchos::ScalarTraits<scalar_type> SCT;
typedef SCT::magnitudeType magnitude_type;
typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> MV;
typedef Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type, node_type> OP;
typedef Belos::MultiVecTraits<scalar_type, MV> MVT;
typedef Belos::OperatorTraits<scalar_type, MV, OP> OPT;
typedef Teuchos::SerialDenseMatrix<int, scalar_type> serial_matrix_type;
typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;
typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> sparse_matrix_type;

/// \fn main
/// \brief Test driver for (Mat)OrthoManager subclasses
int 
main (int argc, char *argv[]) 
{
  using Belos::OrthoManager;
  using Belos::OrthoManagerFactory;
  using Belos::OutputManager;
  using Teuchos::CommandLineProcessor;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &std::cout);
  RCP<const Teuchos::Comm<int> > pComm = 
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  // This factory object knows how to make a (Mat)OrthoManager
  // subclass, given a name for the subclass.  The name is not the
  // same as the class' syntactic name: e.g., "TSQR" is the name of
  // TsqrOrthoManager.
  OrthoManagerFactory<scalar_type, MV, OP> factory;

  // The name of the (Mat)OrthoManager subclass to instantiate.
  std::string orthoManName (factory.defaultName());

  // For SimpleOrthoManager: the normalization method to use.  Valid
  // values: "MGS", "CGS".
  std::string normalization ("CGS");

  // Name of the Harwell-Boeing sparse matrix file from which to read
  // the inner product operator matrix.  If name is "" or not provided
  // at the command line, use the standard Euclidean inner product.
  std::string filename;

  bool verbose = false;
  bool debug = false;

  // The OrthoManager is tested with three different multivectors: S
  // (sizeS columns), X1 (sizeX1 columns), and X2 (sizeX2 columns).
  // The values below are defaults and may be changed by the
  // corresponding command-line arguments.
  int sizeS  = 5;
  int sizeX1 = 11;
  int sizeX2 = 13;

  // Default _global_ number of rows.  The number of rows per MPI
  // process must be no less than max(sizeS, sizeX1, sizeX2).  To
  // ensure that the test always passes with default parameters, we
  // scale by the number of processes.  The default value below may be
  // changed by a command-line parameter with a corresponding name.
  int numRows = 100 * pComm->getSize();

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption ("verbose", "quiet", &verbose,
		  "Print messages and results.");
  cmdp.setOption ("debug", "nodebug", &debug,
		  "Print debugging information.");
  cmdp.setOption ("filename", &filename,
		  "Filename of a Harwell-Boeing sparse matrix, used as the "
		  "inner product operator by the orthogonalization manager."
		  "  If not provided, no matrix is read and the Euclidean "
		  "inner product is used.");
  {
    std::ostringstream os;
    const int numValid = factory.numOrthoManagers();
    const bool plural = numValid > 1 || numValid == 0;

    os << "OrthoManager subclass to test.  There ";
    os << (plural ? "are " : "is ") << numValid << (plural ? "s: " : ": ");
    factory.printValidNames (os);
    os << ".";
    cmdp.setOption ("ortho", &orthoManName, os.str().c_str());
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
    const CommandLineProcessor::EParseCommandLineReturn parseResult = cmdp.parse (argc,argv);
    // If the caller asks us to print the documentation, or does not
    // provide the name of an OrthoManager subclass, we let the "test"
    // pass trivially.
    if (parseResult == CommandLineProcessor::PARSE_HELP_PRINTED)
      {
	if (Teuchos::rank(*pComm) == 0)
	  std::cout << "End Result: TEST PASSED" << endl;
	return EXIT_SUCCESS;
      }
    TEUCHOS_TEST_FOR_EXCEPTION(parseResult != CommandLineProcessor::PARSE_SUCCESSFUL, 
		       std::invalid_argument, 
		       "Failed to parse command-line arguments");
  }
  //
  // Validate command-line arguments
  //
  TEUCHOS_TEST_FOR_EXCEPTION(numRows <= 0, std::invalid_argument, "numRows <= 0 is not allowed");
  TEUCHOS_TEST_FOR_EXCEPTION(numRows <= sizeS + sizeX1 + sizeX2, std::invalid_argument, 
		     "numRows <= sizeS + sizeX1 + sizeX2 is not allowed");
    
  // Declare an output manager for handling local output.  Initialize,
  // using the caller's desired verbosity level.
  //
  // NOTE In Anasazi, this class is called BasicOutputManager.  In
  // Belos, this class is called OutputManager.  We should eventually
  // resolve the difference.
  RCP<OutputManager<scalar_type> > outMan = 
    Belos::Test::makeOutputManager<scalar_type> (verbose, debug);

  // Stream for debug output.  If debug output is not enabled, then
  // this stream doesn't print anything sent to it (it's a "black
  // hole" stream).
  std::ostream& debugOut = outMan->stream(Belos::Debug);
  Belos::Test::printVersionInfo (debugOut);

  // Load the inner product operator matrix from the given filename.
  // If filename == "", use the identity matrix as the inner product
  // operator (the Euclidean inner product), and leave M as
  // Teuchos::null.  Also return an appropriate Map (which will
  // always be initialized; it should never be Teuchos::null).
  RCP<map_type> map;
  RCP<sparse_matrix_type> M; 
  {
    using Belos::Test::loadSparseMatrix;
    // If the sparse matrix is loaded successfully, this call will
    // modify numRows to be the number of rows in the sparse matrix.
    // Otherwise, it will leave numRows alone.
    std::pair<RCP<map_type>, RCP<sparse_matrix_type> > results = 
      loadSparseMatrix<local_ordinal_type, global_ordinal_type, node_type> (pComm, filename, numRows, debugOut);
    map = results.first;
    M = results.second;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(map.is_null(), std::runtime_error,
		     "Error: (Mat)OrthoManager test code failed to "
		     "initialize the Map");
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
    // getNodeNumElements() returns a size_t, which is unsigned, and
    // you shouldn't compare signed and unsigned values.  This is why
    // we make maxNormalizeNumCols a size_t as well.
    if (map->getNodeNumElements() < maxNormalizeNumCols)
      {
	std::ostringstream os;
	os << "The number of elements on this process " << pComm->getRank() 
	   << " is too small for the number of columns that you want to test."
	   << "  There are " << map->getNodeNumElements() << " elements on "
	  "this process, but the normalize() method of the MatOrthoManager "
	  "subclass will need to process a multivector with " 
	   << maxNormalizeNumCols << " columns.  Not all MatOrthoManager "
	  "subclasses can handle a local row block with fewer rows than "
	  "columns.";
	throw std::invalid_argument(os.str());
      }
  }

  // Using the factory object, instantiate the specified OrthoManager
  // subclass to be tested.  Specify "default" parameters (which
  // should favor accuracy over performance), but override the default
  // parameters to get the desired normalization method for
  // SimpleOrthoManaager.
  RCP<OrthoManager<scalar_type, MV> > orthoMan;
  {
    std::string label (orthoManName);
    RCP<ParameterList> params = 
      parameterList (*(factory.getDefaultParameters (orthoManName)));
    if (orthoManName == "Simple") {
      params->set ("Normalization", normalization);
      label = label + " (" + normalization + " normalization)";
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

  // Only Rank 0 gets to write to cout.  The other processes dump
  // output to a black hole.
  //std::ostream& finalOut = (Teuchos::rank(*pComm) == 0) ? std::cout : Teuchos::oblackholestream;

  if (numFailed != 0)
    {
      outMan->stream(Belos::Errors) << numFailed << " errors." << endl;

      // The Trilinos test framework depends on seeing this message,
      // so don't rely on the OutputManager to report it correctly.
      if (Teuchos::rank(*pComm) == 0)
	std::cout << "End Result: TEST FAILED" << endl;	
      return EXIT_FAILURE;
    }
  else 
    {
      if (Teuchos::rank(*pComm) == 0)
	std::cout << "End Result: TEST PASSED" << endl;
      return EXIT_SUCCESS;
    }
}



