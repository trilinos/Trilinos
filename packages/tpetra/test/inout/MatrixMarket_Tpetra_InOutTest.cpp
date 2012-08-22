// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_DefaultPlatform.hpp>

#include <Kokkos_ConfigDefs.hpp>
#include <Kokkos_SerialNode.hpp>

#if defined(HAVE_KOKKOSCLASSIC_TBB)
#  include <Kokkos_TBBNode.hpp>
#endif // defined(HAVE_KOKKOSCLASSIC_TBB)

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <algorithm>

using std::endl;

namespace Tpetra {
  namespace MatrixMarket {
    namespace Test {

      /// \fn getNode
      /// \brief Return an RCP to a Kokkos Node
      ///
      template<class NodeType>
      Teuchos::RCP<NodeType>
      getNode() {
	throw std::runtime_error ("This Kokkos Node type not supported (compile-time error)");
      }

      template<>
      Teuchos::RCP<Kokkos::SerialNode>
      getNode() {
	Teuchos::ParameterList defaultParams;
	return Teuchos::rcp (new Kokkos::SerialNode (defaultParams));
      }

#if defined(HAVE_KOKKOSCLASSIC_TBB)
      template<>
      Teuchos::RCP<Kokkos::TBBNode>
      getNode() {
	// "Num Threads" specifies the number of threads.  Defaults to an
	// automatically chosen value.
	Teuchos::ParameterList defaultParams;
	return Teuchos::rcp (new Kokkos::TBBNode (defaultParams));
      }
#endif // defined(HAVE_KOKKOSCLASSIC_TBB)

      /// Test Tpetra::MatrixMarket::Reader::readSparseFile()
      ///
      /// \param inputFilename [in] Name of the Matrix Market format
      ///   sparse matrix file to read (on MPI Rank 0 only).
      /// \param outputFilename [in] Name of the Matrix Market format
      ///   sparse matrix file to write (on MPI Rank 0 only).  Ignored
      ///   if testWrite==false.
      /// \param pComm [in] Communicator, over whose MPI ranks to
      ///   distribute the returned Tpetra::CrsMatrix.
      /// \param echo [in] Whether or not to echo the resulting 
      ///   matrix to cout in Matrix Market format.
      /// \param tolerant [in] Whether or not to parse the file 
      ///   tolerantly.
      /// \param verbose [in] Whether to print verbose output.
      /// \param debug [in] Whether to print debugging output.
      void
      testReadAndWriteSparseFile (const std::string& inputFilename, 
				  const std::string& outputFilename, 
				  const Teuchos::RCP<const Teuchos::Comm<int> >& pComm,
				  const bool testWrite,
				  const bool echo,
				  const bool tolerant, 
				  const bool verbose,
				  const bool debug)
      {
	using Teuchos::RCP;
	using std::cerr;
	using std::cout;
	using std::endl;

	typedef double scalar_type;
	typedef int local_ordinal_type;
	typedef int global_ordinal_type;
	// typedef size_t global_ordinal_type;
	// #if defined(HAVE_KOKKOSCLASSIC_TBB)
	//       typedef Kokkos::TBBNode node_type;
	// #else
	typedef Kokkos::SerialNode node_type;
	// #endif // defined(HAVE_KOKKOSCLASSIC_TBB)
	typedef Teuchos::ScalarTraits<scalar_type> STS;

	const bool callFillComplete = true;

	// Get a Kokkos Node instance for the particular Node type.
	RCP<node_type> pNode = getNode<node_type>();
	const int myRank = Teuchos::rank (*pComm);

	if (verbose && myRank == 0)
	  cout << "About to read Matrix Market file \"" << inputFilename 
	       << "\""
	       << (callFillComplete ? " (calling fillComplete())" : "")
	       << ":" << endl;

	// Read the sparse matrix from the given Matrix Market file.
	// This routine acts like an MPI barrier.
	typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, 
	  global_ordinal_type, node_type> sparse_matrix_type;
	typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type> reader_type;
	RCP<sparse_matrix_type> pMatrix =
	  reader_type::readSparseFile (inputFilename, pComm, pNode, 
				       callFillComplete, tolerant, debug);
	TEUCHOS_TEST_FOR_EXCEPTION(pMatrix.is_null(), std::runtime_error,
			   "The Tpetra::CrsMatrix returned from "
			   "readSparseFile() is null.");
	TEUCHOS_TEST_FOR_EXCEPTION(callFillComplete && ! pMatrix->isFillComplete(), 
			   std::logic_error,
			   "We asked readSparseFile() to call fillComplete() "
			   "on the Tpetra::CrsMatrix before returning it, but"
			   " it did not.");

	if (! pMatrix.is_null() && verbose && myRank == 0)
	  cout << "Successfully read Matrix Market file \"" << inputFilename 
	       << "\"." << endl;

	typedef Tpetra::MatrixMarket::Writer<sparse_matrix_type> writer_type;
	if (testWrite && outputFilename != "")
	  writer_type::writeSparseFile (outputFilename, pMatrix, debug);
	if (echo)
	  writer_type::writeSparse (cout, pMatrix, debug);
      }
    } // namespace Test
  } // namespace MatrixMarket
} // namespace Tpetra

/// \fn main
/// \brief Benchmark driver
int 
main (int argc, char *argv[]) 
{
  using Teuchos::Comm;
  using Teuchos::CommandLineProcessor;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::cout;
  using std::endl;

  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &cout);
  RCP<const Comm<int> > pComm = 
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  std::string inputFilename;  // Matrix Market file to read
  std::string outputFilename; // Matrix Market file to write (if applicable)
  bool testWrite = false; // Test Matrix Market output?
  bool tolerant = false; // Parse the file tolerantly?
  bool echo = false;     // Echo the read-in matrix back?
  bool verbose = false;  // Verbosity of output
  bool debug = false;    // Print debugging info?

  CommandLineProcessor cmdp (false, true);
  cmdp.setOption ("inputFilename", &inputFilename,
		  "Name of the Matrix Market sparse matrix file to read.");
  cmdp.setOption ("outputFilename", &outputFilename,
		  "If --testWrite is true, then write the read-in matrix to "
		  "the given file in Matrix Market format on (MPI) Proc 0.  "
		  "Otherwise, this argument is ignored.  Note that symmetric"
		  " storage will have been expanded and comments will have "
		  "been stripped, so the output file may not be identical to"
		  " the input file.");
  cmdp.setOption ("testWrite", "noTestWrite", &testWrite,
		  "Whether to test Matrix Market file output.");
  cmdp.setOption ("tolerant", "strict", &tolerant, 
		  "Whether to parse the Matrix Market file tolerantly.");
  cmdp.setOption ("echo", "noecho", &echo,
		  "Whether to echo the read-in matrix back to stdout on Rank 0 "
		  "in Matrix Market format.  Symmetric storage will have been "
		  "expanded, so the result will not be identical to the input "
		  "file, though the matrix represented will be the same.");
  cmdp.setOption ("verbose", "quiet", &verbose,
		  "Print messages and results.");
  cmdp.setOption ("debug", "nodebug", &debug,
  		  "Print debugging information.");

  // Parse the command-line arguments.
  {
    const CommandLineProcessor::EParseCommandLineReturn parseResult = 
      cmdp.parse (argc,argv);
    // If the caller asks us to print the documentation, or does not
    // explicitly say to run the benchmark, we let this "test" pass
    // trivially.
    if (parseResult == CommandLineProcessor::PARSE_HELP_PRINTED)
      {
	if (Teuchos::rank(*pComm) == 0)
	  cout << "End Result: TEST PASSED" << endl;
	return EXIT_SUCCESS;
      }
    TEUCHOS_TEST_FOR_EXCEPTION(parseResult != CommandLineProcessor::PARSE_SUCCESSFUL, 
		       std::invalid_argument, 
		       "Failed to parse command-line arguments");
  }

  // Test reading in the sparse matrix.  If no filename or an empty
  // filename is specified, we don't invoke the test and report a
  // "TEST PASSED" message.
  if (inputFilename != "")
    {
      using Tpetra::MatrixMarket::Test::testReadAndWriteSparseFile;
      testReadAndWriteSparseFile (inputFilename, outputFilename, pComm, 
				  testWrite, echo, tolerant, verbose, debug);
    }

  // Only Rank 0 gets to write to cout.
  if (Teuchos::rank(*pComm) == 0)
    std::cout << "End Result: TEST PASSED" << endl;
  return EXIT_SUCCESS;
}



