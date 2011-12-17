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

#include <MatrixMarket_raw.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
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

#if defined(HAVE_KOKKOS_TBB)
      template<>
      Teuchos::RCP<Kokkos::TBBNode>
      getNode() {
	// "Num Threads" specifies the number of threads.  Defaults to an
	// automatically chosen value.
	Teuchos::ParameterList defaultParams;
	return Teuchos::rcp (new Kokkos::TBBNode (defaultParams));
      }
#endif // defined(HAVE_KOKKOS_TBB)

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

  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &std::cout);
  RCP<const Comm<int> > pComm = 
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  // Name of the Matrix Market sparse matrix file to read.
  std::string filename;
  // Whether to echo the sparse matrix to stdout after reading it
  // successfully.
  bool echo = false;
  // Whether to parse the Matrix Market file tolerantly.
  bool tolerant = false;
  // Verbosity of output
  bool verbose = false; 
  // Whether to print debugging-level output
  bool debug = false;   

  CommandLineProcessor cmdp (false, true);
  cmdp.setOption ("filename", &filename,
		  "Name of the Matrix Market sparse matrix file to read.");
  cmdp.setOption ("echo", "noecho", &echo, 
		  "Whether to echo the sparse matrix contents to stdout "
		  "after reading it successfully.");
  cmdp.setOption ("tolerant", "strict", &tolerant, 
		  "Whether to parse the Matrix Market file tolerantly.");
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
	  std::cout << "End Result: TEST PASSED" << endl;
	return EXIT_SUCCESS;
      }
    TEUCHOS_TEST_FOR_EXCEPTION(parseResult != CommandLineProcessor::PARSE_SUCCESSFUL, 
		       std::invalid_argument, 
		       "Failed to parse command-line arguments");
  }

  // Test reading in the sparse matrix.  If no filename or an empty
  // filename is specified, the test passes trivially.
  int success;
  if (filename == "")
    success = 1;
  else
    {
      using Tpetra::MatrixMarket::Raw::Reader;
      typedef double scalar_type;
      typedef int ordinal_type;
      typedef Reader<scalar_type, ordinal_type> reader_type;
      const bool theSuccess = 
	reader_type::readFile (*pComm, filename, echo, tolerant, debug);
      success = theSuccess ? 1 : 0;
    }
  if (debug)
    {
      // Make sure that all the processes finish.  This should not be
      // necessary, since readFile is a collective for which all ranks
      // agree on the returned Boolean result.
      Teuchos::barrier (*pComm);
    }

  // Only Rank 0 gets to write to cout.
  if (Teuchos::rank(*pComm) == 0)
    {
      if (success)
	std::cout << "End Result: TEST PASSED" << endl;
      else
	std::cout << "End Result: TEST FAILED" << endl;
    }
  if (success)
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}



