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

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <algorithm>

namespace {
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

  // Test Tpetra::MatrixMarket::Reader::readSparseFile()
  //
  // \param inputFilename [in] Name of the Matrix Market format
  //   dense matrix file to read (on MPI Rank 0 only).
  // \param outputFilename [in] Name of the Matrix Market format
  //   dense matrix file to write (on MPI Rank 0 only).  Ignored
  //   if testWrite==false.
  // \param comm [in] Communicator, over whose MPI ranks to
  //   distribute the returned Tpetra::MultiVector.
  // \param echo [in] Whether or not to echo the resulting 
  //   matrix to cout in Matrix Market format.
  // \param tolerant [in] Whether or not to parse the file 
  //   tolerantly.
  // \param verbose [in] Whether to print verbose output.
  // \param debug [in] Whether to print debugging output.
  template<class ScalarType, class LO, class GO, class NodeType>
  void
  testReadAndWriteSparseFile (const std::string& inputFilename, 
			      const std::string& outputFilename, 
			      const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
			      const Teuchos::RCP<NodeType>& node,
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

    typedef ScalarType scalar_type;
    typedef LO local_ordinal_type;
    typedef GO global_ordinal_type;
    typedef NodeType node_type;
    typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> 
      map_type;
    typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, 
      global_ordinal_type, node_type> multivector_type;

    // The reader and writer classes are templated on the
    // Tpetra::CrsMatrix specialization, from which the
    // Tpetra::MultiVector specialization is derived.
    typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, 
      global_ordinal_type, node_type> sparse_matrix_type;

    const int myRank = comm->getRank ();

    if (verbose && myRank == 0) {
      cout << "About to read Matrix Market dense file \"" << inputFilename 
	   << "\":" << endl;
    }

    // Map describing multivector distribution; starts out null and is
    // an output argument of readDenseFile().
    RCP<map_type> map; 

    // Read the dense matrix from the given Matrix Market file.
    // This routine acts like an MPI barrier.
    typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type> reader_type;
    RCP<multivector_type> X =
      reader_type::readDenseFile (inputFilename, comm, node, map, tolerant, debug);

    TEUCHOS_TEST_FOR_EXCEPTION(X.is_null(), std::runtime_error,
      "The Tpetra::MultiVector returned from readDenseFile() is null.");
    if (! X.is_null() && verbose && myRank == 0) {
      cout << "Successfully read Matrix Market file \"" << inputFilename 
	   << "\"." << endl;
    }

    // If specified, write the read-in sparse matrix to a file and/or
    // echo it to stdout.
    typedef Tpetra::MatrixMarket::Writer<sparse_matrix_type> writer_type;
    if (testWrite && outputFilename != "") {
      writer_type::writeSparseFile (outputFilename, X, debug);
    }
    if (echo) {
      writer_type::writeSparse (cout, X, debug);
    }
  }
} // namespace (anonymous)


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

  typedef double scalar_type;
  typedef int local_ordinal_type;
  typedef int global_ordinal_type;
  typedef Kokkos::SerialNode node_type;

  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &cout);
  RCP<const Comm<int> > comm = 
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
		  "Name of the Matrix Market dense matrix file to read.");
  cmdp.setOption ("outputFilename", &outputFilename,
		  "If --testWrite is true, then write the read-in matrix to "
		  "the given file in Matrix Market format on (MPI) Proc 0.  "
		  "Otherwise, this argument is ignored.  Note that the output "
		  "file may not be identical to the input file.");
  cmdp.setOption ("testWrite", "noTestWrite", &testWrite,
		  "Whether to test Matrix Market file output.");
  cmdp.setOption ("tolerant", "strict", &tolerant, 
		  "Whether to parse the input Matrix Market file tolerantly.");
  cmdp.setOption ("echo", "noecho", &echo,
		  "Whether to echo the read-in matrix back to stdout on Rank 0 "
		  "in Matrix Market format.  Note that the echoed matrix may "
		  "not be identical to the input file.");
  cmdp.setOption ("verbose", "quiet", &verbose, "Print messages and results.");
  cmdp.setOption ("debug", "nodebug", &debug, "Print debugging information.");

  // Parse the command-line arguments.
  {
    const CommandLineProcessor::EParseCommandLineReturn parseResult = 
      cmdp.parse (argc,argv);
    // If the caller asks us to print the documentation, or does not
    // explicitly say to run the benchmark, we let this "test" pass
    // trivially.
    if (parseResult == CommandLineProcessor::PARSE_HELP_PRINTED) {
      if (Teuchos::rank(*pComm) == 0) {
	cout << "End Result: TEST PASSED" << endl;
      }
      return EXIT_SUCCESS;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(parseResult != CommandLineProcessor::PARSE_SUCCESSFUL, 
      std::invalid_argument, "Failed to parse command-line arguments.");
  }

  // Get a Kokkos Node instance for the particular Node type.
  RCP<node_type> node = getNode<node_type>();

  // Test reading in the sparse matrix.  If no filename or an empty
  // filename is specified, we don't invoke the test and report a
  // "TEST PASSED" message.
  if (inputFilename != "") {
    using Tpetra::MatrixMarket::Test::testReadAndWriteDenseFile;
    testReadAndWriteDenseFile<double, int, int, node_type> (inputFilename, outputFilename, comm, node, testWrite, echo, tolerant, verbose, debug);
  }

  // Only Rank 0 gets to write to cout.
  if (comm->getRank() == 0) {
    cout << "End Result: TEST PASSED" << endl;
  }
  return EXIT_SUCCESS;
}



