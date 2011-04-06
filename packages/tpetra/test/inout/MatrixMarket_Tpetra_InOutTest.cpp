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

      /// Test Tpetra::MatrixMarket::Reader::readSparseFile()
      ///
      /// \param filename [in] Name of the Matrix Market format sparse
      ///   matrix file to read (on MPI Rank 0 only).
      /// \param pComm [in] Communicator, over whose MPI ranks to
      ///   distribute the returned Tpetra::CrsMatrix.
      /// \param echo [in] Whether or not to echo the resulting 
      ///   matrix to cout in Matrix Market format.
      /// \param tolerant [in] Whether or not to parse the file 
      ///   tolerantly.
      /// \param verbose [in] Whether to print verbose output.
      /// \param debug [in] Whether to print debugging output.
      ///
      void
      testReadSparseFile (const std::string& filename, 
			  const Teuchos::RCP<const Teuchos::Comm<int> >& pComm,
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
	// #if defined(HAVE_KOKKOS_TBB)
	//       typedef Kokkos::TBBNode node_type;
	// #else
	typedef Kokkos::SerialNode node_type;
	// #endif // defined(HAVE_KOKKOS_TBB)

	typedef Teuchos::ScalarTraits<scalar_type> STS;

	// Get a Kokkos Node instance for the particular Node type.
	RCP<node_type> pNode = getNode<node_type>();

	const int numProcs = Teuchos::size (*pComm);
	const int myRank = Teuchos::rank (*pComm);

	if (verbose && myRank == 0)
	  cout << "About to read Matrix Market file \"" << filename << "\":" << endl;

	// Read the sparse matrix from the given Matrix Market file.
	// This routine acts like an MPI barrier.
	const bool callFillComplete = true;
	typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, 
	  global_ordinal_type, node_type> sparse_matrix_type;
	typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type> reader_type;
	RCP<sparse_matrix_type> pMatrix =
	  reader_type::readSparseFile (filename, pComm, pNode, 
				       callFillComplete, tolerant, debug);
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

	if (echo)
	  {
	    // Number of errors encountered while printing out the matrix.
	    int numPrintErrors = 0;

	    // Rank 0: Print the banner line and the dimensions line to cout.
	    if (myRank == 0)
	      {
		cout << "%%MatrixMarket matrix coordinate ";
		if (STS::isComplex)
		  cout << "complex ";
		else
		  cout << "real ";
		cout << "general" << endl;

		// getGlobalNum{Rows,Cols}() does not return what you
		// think it should return.  Instead, ask the range
		// resp. domain map for the number of rows
		// resp. columns in the matrix.
		cout << pMatrix->getRangeMap()->getGlobalNumElements() 
		     << " "
		     << pMatrix->getDomainMap()->getGlobalNumElements() 
		     << " "
		     << pMatrix->getGlobalNumEntries()
		     << endl;
	      }
	    Teuchos::barrier (*pComm);

	    // Let each processor in turn print to cout its rows.  We
	    // assume that all processors can print to cout.  We do
	    // _not_ assume here that the row map is one-to-one;
	    // printing should work just fine, as long as nonzeros
	    // themselves are not stored redundantly.
	    for (int p = 0; p < numProcs; ++p)
	      {
		typedef ArrayView<global_ordinal_type>::size_type size_type;
		if (myRank == p)
		  {
		    // Storage for column indices and values in each
		    // row.  Will be resized as necessary.  (Proc p
		    // may not own any rows, in which case Proc p
		    // won't need to allocate these at all.
		    Array<global_ordinal_type> indices;
		    Array<scalar_type> values;

		    // List of the rows with storage on Proc p.
		    ArrayView<const global_ordinal_type> myRows = 
		      pMatrix->getRowMap()->getNodeElementList();
		    // Number of rows with storage on Proc p.
		    const size_type myNumRows = myRows.size();

		    // For each row that Proc p owns, print its
		    // entries to cout.
		    for (size_type k = 0; k < myNumRows; ++k)
		      {
			const global_ordinal_type curRow = myRows[k];
			size_t numEntries = 
			  pMatrix->getNumEntriesInGlobalRow (curRow);

			// Resize (if necessary) the arrays for
			// holding column indices and values for the
			// current row.
			//
			// Signed to unsigned integer conversion, for
			// integers of the same size, shouldn't
			// overflow.
			if (static_cast<size_t> (indices.size()) < numEntries)
			  indices.resize (numEntries);
			if (static_cast<size_t> (values.size()) < numEntries)
			  values.resize (numEntries);
			// This views are exactly the right length to
			// hold the data for the current row.  indices
			// and values may be longer than necessary;
			// that's an optimization, to avoid resizing
			// them with every row.
			ArrayView<global_ordinal_type> indicesView = 
			  indices.view (0, numEntries);
			ArrayView<scalar_type> valuesView = 
			  values.view (0, numEntries);

			// Make sure there were no surprises with
			// the number of entries.
			size_t newNumEntries = 0;
			pMatrix->getGlobalRowCopy (curRow, indicesView,
						   valuesView, newNumEntries);
			if (newNumEntries != numEntries)
			  numPrintErrors++;
			else
			  {
			    for (size_t j = 0; j < numEntries; ++j)
			      cout << curRow << " " 
				   << indicesView[j] << " "
				   << valuesView[j] << endl;
			  }
		      }
		  }
		Teuchos::barrier (*pComm);

		// If there were any errors on any processors, stop right away.
		int totalNumPrintErrors = 0;
		Teuchos::reduceAll (*pComm, Teuchos::REDUCE_SUM, numPrintErrors, 
				    Teuchos::Ptr<int> (&totalNumPrintErrors));
		TEST_FOR_EXCEPTION(totalNumPrintErrors > 0, std::runtime_error,
				   "Failed to print Tpetra::CrsMatrix.  Total "
				   "number of print errors thus far: " 
				   << numPrintErrors);
	      }
	  }
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

  std::string filename;  // Matrix Market file to read
  bool tolerant = false; // Parse the file tolerantly?
  bool echo = false;     // Echo the read-in matrix back?
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
  cmdp.setOption ("echo", "noecho", &echo,
		  "Whether to echo the read-in matrix back to stdout on Rank 0 "
		  "in Matrix Market format.  Symmetric storage will have been "
		  "expanded, so the result will not be identical to the input "
		  "file, though the matrix represented will be the same.");

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
    TEST_FOR_EXCEPTION(parseResult != CommandLineProcessor::PARSE_SUCCESSFUL, 
		       std::invalid_argument, 
		       "Failed to parse command-line arguments");
  }

  // Test reading in the sparse matrix.  If no filename or an empty
  // filename is specified, we don't invoke the test and report a
  // "TEST PASSED" message.
  if (filename != "")
    {
      using Tpetra::MatrixMarket::Test::testReadSparseFile;
      testReadSparseFile (filename, pComm, echo, tolerant, verbose, debug);
    }

  // Only Rank 0 gets to write to cout.
  if (Teuchos::rank(*pComm) == 0)
    std::cout << "End Result: TEST PASSED" << endl;
  return EXIT_SUCCESS;
}



