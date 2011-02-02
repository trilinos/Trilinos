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

#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_DefaultPlatform.hpp>

#include <Kokkos_ConfigDefs.hpp>
#include <Kokkos_SerialNode.hpp>
#if defined(HAVE_KOKKOS_TBB)
#  include <Kokkos_TBBNode.hpp>
#endif // defined(HAVE_KOKKOS_TBB)

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <algorithm>

using std::endl;

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

    /// Test MatrixMarket::Tpetra::Reader::readFile()
    ///
    /// \param filename [in] Name of the Matrix Market format sparse
    ///   matrix file to read (on MPI Rank 0 only)
    /// \param pComm [in] Communicator, over whose MPI ranks to
    ///   distribute the returned Tpetra::CrsMatrix.
    /// \param tolerant [in] Whether or not to parse the file 
    ///   tolerantly
    ///
    /// \return Tpetra::CrsMatrix read in from the file 
    template<class SparseMatrixType>
    Teuchos::RCP<SparseMatrixType>
    readFile (const std::string& filename,
	      const Teuchos::RCP<const Teuchos::Comm<int> >& pComm, 
	      const bool tolerant,
	      const bool debug)
    {
      typedef typename SparseMatrixType::node_type node_type;
      Teuchos::RCP<node_type> pNode = getNode<node_type>();
      
      typedef MatrixMarket::Tpetra::Reader<SparseMatrixType> reader_type;
      return reader_type::readFile (filename, pComm, pNode, tolerant, debug);
    }

    void
    testReadFile (const std::string& filename, 
		  const Teuchos::RCP<const Teuchos::Comm<int> >& pComm,
		  const bool tolerant, 
		  const bool verbose,
		  const bool debug)
    {
      using Teuchos::RCP;
      using std::cout;
      using std::endl;

      typedef double scalar_type;
      typedef int local_ordinal_type;
      typedef size_t global_ordinal_type;
#if defined(HAVE_KOKKOS_TBB)
      typedef Kokkos::TBBNode node_type;
#else
      typedef Kokkos::SerialNode node_type;
#endif // defined(HAVE_KOKKOS_TBB)

      typedef ::Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> sparse_matrix_type;
      const int myRank = Teuchos::rank (*pComm);
      if (verbose && myRank == 0)
	cout << "About to read Matrix Market file \"" << filename << "\":" << endl;

      // Read the sparse matrix from the given Matrix Market file.
      // This routine acts like an MPI barrier.
      RCP<sparse_matrix_type> pMatrix = 
	readFile<sparse_matrix_type> (filename, pComm, tolerant, debug);
      if (! pMatrix.is_null())
	{
	  if (verbose && myRank == 0)
	    cout << "Successfully read Matrix Market file \"" << filename << "\"." << endl;
	}
      else 
	{
	  if (verbose && myRank == 0)
	    cout << "Failed to read Matrix Market file \"" << filename << "\"." << endl;
	}
    }
  } // namespace Test
} // namespace MatrixMarket


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

  std::string filename;  // Matrix Market file to read
  bool tolerant = false; // Parse the file tolerantly?
  bool verbose = false;  // Verbosity of output
  bool debug = false;    // Print debugging info?

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
  // filename is specified, we don't invoke the test and report a
  // "TEST PASSED" message.
  if (filename != "")
    MatrixMarket::Test::testReadFile (filename, pComm, tolerant, verbose, debug);

  // Only Rank 0 gets to write to cout.
  if (Teuchos::rank(*pComm) == 0)
    std::cout << "End Result: TEST PASSED" << endl;
  return EXIT_SUCCESS;
}



