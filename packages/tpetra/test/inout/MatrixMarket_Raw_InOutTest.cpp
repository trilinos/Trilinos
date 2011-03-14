// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
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

  // Whether to parse the Matrix Market file tolerantly.
  bool tolerant = false;

  bool verbose = false; // Verbosity of output
  bool debug = false;   // Whether to print debugging-level output

  CommandLineProcessor cmdp (false, true);
  cmdp.setOption ("verbose", "quiet", &verbose,
		  "Print messages and results.");
  cmdp.setOption ("debug", "nodebug", &debug,
		  "Print debugging information.");
  cmdp.setOption ("filename", &filename,
		  "Name of the Matrix Market sparse matrix file to read.");
  cmdp.setOption ("tolerant", "strict", &tolerant, 
		  "Whether to parse the Matrix Market file tolerantly.");
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
    TEST_FOR_EXCEPTION(parseResult != CommandLineProcessor::PARSE_SUCCESSFUL, 
		       std::invalid_argument, 
		       "Failed to parse command-line arguments");
  }

  // Test reading in the sparse matrix.  If no filename or an empty
  // filename is specified, the test passes trivially.
  if (filename != "")
    {
      using Tpetra::MatrixMarket::Raw::Reader;

      typedef double scalar_type;
      typedef int ordinal_type;
      typedef Reader<scalar_type, ordinal_type> reader_type;
      reader_type::readFile (filename, tolerant, debug);
    }

  // Only Rank 0 gets to write to cout.
  if (Teuchos::rank(*pComm) == 0)
    std::cout << "End Result: TEST PASSED" << endl;
  return EXIT_SUCCESS;
}



