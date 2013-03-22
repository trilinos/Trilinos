/*
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
*/

#include <Teuchos_UnitTestHarness.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_ETIHelperMacros.h>

#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <MatrixMarket_Tpetra.hpp>

#include <Kokkos_DefaultNode.hpp>

#include <Teuchos_Array.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_TypeNameTraits.hpp>

//
// Test for Tpetra::CrsMatrix::sumIntoGlobalValues(), with nonowned
// rows.  This test is like CrsMatrix_NonlocalSumInto.cpp, except that
// it attempts to sum into remote entries that don't exist on the
// process that owns them.  Currently, CrsMatrix silently ignores
// these entries.  (This is how CrsMatrix implements Import and Export
// when the target matrix has a fixed column Map.  Data are
// redistributed between the two row Maps, and "filtered" by the
// target matrix's column Map.)  This unit test verifies that behavior
// by ensuring the following:
//
// 1. fillComplete() (actually globalAssemble()) does not throw an
//    exception when the incoming entries don't exist on the process
//    that owns their rows.
//
// 2. The ignored entries are actually ignored.  They must change
//    neither the structure nor the values of the matrix.
//
// mfh 16 Dec 2012: The one-template-argument version breaks explicit
// instantiation.  Ah well.
//
//TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( CrsMatrix, NonlocalSumInto_Ignore, CrsMatrixType )
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, NonlocalSumInto_Ignore, LocalOrdinalType, GlobalOrdinalType, ScalarType, NodeType )
{
  using Tpetra::createContigMapWithNode;
  using Tpetra::createNonContigMapWithNode;
  using Tpetra::global_size_t;
  using Tpetra::Map;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::av_const_cast;
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::OrdinalTraits;
  using Teuchos::outArg;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::reduceAll;
  using Teuchos::ScalarTraits;
  using Teuchos::tuple;
  using Teuchos::TypeNameTraits;
  using std::endl;

#if 0
  // Extract typedefs from the CrsMatrix specialization.
  typedef typename CrsMatrixType::scalar_type scalar_type;
  typedef typename CrsMatrixType::local_ordinal_type local_ordinal_type;
  typedef typename CrsMatrixType::global_ordinal_type global_ordinal_type;
  typedef typename CrsMatrixType::node_type node_type;
#endif // 0

  typedef ScalarType scalar_type;
  typedef LocalOrdinalType local_ordinal_type;
  typedef GlobalOrdinalType global_ordinal_type;
  typedef NodeType node_type;

  // Typedefs derived from the above canonical typedefs.
  typedef ScalarTraits<scalar_type> STS;
  typedef typename STS::magnitudeType magnitude_type;
  typedef ScalarTraits<magnitude_type> STM;
  typedef Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

  // Abbreviation typedefs.
  typedef scalar_type ST;
  typedef local_ordinal_type LO;
  typedef global_ordinal_type GO;
  typedef node_type NT;

  typedef Tpetra::CrsMatrix<ST, LO, GO, NT> CrsMatrixType;

  // CrsGraph specialization corresponding to CrsMatrixType (the
  // CrsMatrix specialization).
  typedef Tpetra::CrsGraph<LO, GO, NT, typename CrsMatrixType::mat_solve_type> crs_graph_type;

  ////////////////////////////////////////////////////////////////////  
  // HERE BEGINS THE TEST.
  ////////////////////////////////////////////////////////////////////

  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();

  // Get the default communicator.
  RCP<const Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
  const int numProcs = comm->getSize ();
  const int myRank = comm->getRank ();

  if (myRank == 0) {
    out << "Test with " << numProcs << " process" << (numProcs != 1 ? "es" : "") << endl;
  }

  // This test doesn't make much sense if there is only one MPI
  // process.  We let it pass trivially in that case.
  if (numProcs == 1) {
    out << "Number of processes in world is one; test passes trivially." << endl;
    return;
  }

  // Get a Kokkos Node instance.  It would be nice if we could pass in
  // parameters here, but threads don't matter for this test; it's a
  // test for distributed-memory capabilities.

  if (myRank == 0) {
    out << "Creating Kokkos Node of type " << TypeNameTraits<node_type>::name () << endl;
  }
  RCP<node_type> node;
  {
    ParameterList pl; // Kokkos Node types require a PL inout.
    node = rcp (new node_type (pl));
  }

  // Number of rows in the matrix owned by each process.
  const LO numLocalRows = 10;

  // Number of (global) rows and columns in the matrix.
  const GO numGlobalRows = numLocalRows * numProcs;
  const GO numGlobalCols = numGlobalRows;
  // Prevent compile warning for unused variable.
  // (It's not really "variable" if it's const, but oh well.)
  (void) numGlobalCols; 

  if (myRank == 0) {
    out << "Creating contiguous row Map" << endl;
  }

  // Create a contiguous row Map, with numLocalRows rows per process.
  RCP<const map_type> rowMap = createContigMapWithNode<LO, GO, NT> (INVALID, numLocalRows, comm, node);

  // For now, reuse the row Map for the domain and range Maps.  Later,
  // we might want to test using different domain or range Maps.
  RCP<const map_type> domainMap = rowMap;
  RCP<const map_type> rangeMap = rowMap;

  // Min and max row and column index of this process.  Use the row
  // Map for the row and column indices, since we're only inserting
  // indices into the graph for rows that the calling process owns.
  const GO globalMinRow = rowMap->getMinGlobalIndex ();
  const GO globalMaxRow = rowMap->getMaxGlobalIndex ();
  const GO globalMinCol = domainMap->getMinAllGlobalIndex ();
  const GO globalMaxCol = domainMap->getMaxAllGlobalIndex ();

  if (myRank == 0) {
    out << "Creating graph" << endl;
  }

  // Create a numGlobalRows by numGlobalCols graph and set its
  // structure.  Every process sets its diagonal entries (which it
  // owns).  Unlike in the CrsMatrix_NonlocalSumInto.cpp test, we
  // don't set any other entries.  As a result, the later calls to
  // sumIntoGlobalValues() for nonowned rows should fail.
  RCP<const crs_graph_type> graph;
  {
    // We have a good upper bound for the number of entries per row,
    // so use static profile.  Leave the upper bound as 2 (just as it
    // is in the CrsMatrix_NonlocalSumInto.cpp test) so that there
    // would actually be room for the incoming entries from remote
    // calls to sumIntoGlobalValues().
    RCP<crs_graph_type> nonconstGraph (new crs_graph_type (rowMap, 2, Tpetra::StaticProfile));

    TEUCHOS_TEST_FOR_EXCEPTION(globalMinRow >= globalMaxRow, std::logic_error, 
      "This test only works if globalMinRow < globalMaxRow.");

    // Insert all the diagonal entries, and only the diagonal entries
    // (unlike in the other test).
    for (GO globalRow = globalMinRow; globalRow <= globalMaxRow; ++globalRow) {
      nonconstGraph->insertGlobalIndices (globalRow, tuple (globalRow));
    }

    nonconstGraph->fillComplete (domainMap, rangeMap);
    graph = rcp_const_cast<const crs_graph_type> (nonconstGraph);
  }

  // Test whether the graph has the correct structure.
  bool localGraphSuccess = true;
  std::ostringstream graphFailMsg;
  {
    Array<GO> ind (2); // upper bound

    for (GO globalRow = globalMinRow; globalRow <= globalMaxRow; ++globalRow) {
      size_t numEntries = 0; // output argument of below line.
      graph->getGlobalRowCopy (globalRow, ind (), numEntries);

      // Revise view based on numEntries.
      ArrayView<GO> indView = ind.view (0, numEntries);

      // Sort the view.
      std::sort (indView.begin (), indView.end ());

      if (numEntries != as<size_t> (1)) {
	localGraphSuccess = false;
	graphFailMsg << "Proc " << myRank << ": globalRow = " << globalRow << ": numEntries = " << numEntries << " != 1" << endl;
      }
      if (numEntries > 0 && indView[0] != globalRow) {
	localGraphSuccess = false;
	graphFailMsg << "Proc " << myRank << ": globalRow = " << globalRow << ": indView[0] = " << indView[0] << " != globalRow = " << globalRow << endl;
      }
    }
  }
 
  // Make sure that all processes successfully created the graph.
  bool globalGraphSuccess = true;
  {
    int globalGraphSuccess_int = 1;
    reduceAll (*comm, Teuchos::REDUCE_MIN, localGraphSuccess ? 1 : 0, outArg (globalGraphSuccess_int));
    globalGraphSuccess = (globalGraphSuccess_int != 0);
  }
  if (! globalGraphSuccess) {
    if (myRank == 0) {
      out << "Graph structure not all correct:" << endl << endl;
    }
    // Print out the failure messages on all processes.
    for (int p = 0; p < numProcs; ++p) {
      if (p == myRank) {
	out << graphFailMsg.str () << endl;
	std::flush (out);
      }
      // Do some barriers to allow output to finish.
      comm->barrier ();
      comm->barrier ();
      comm->barrier ();
    }
  }
  TEUCHOS_TEST_FOR_EXCEPTION(! globalGraphSuccess, std::logic_error, "Graph structure test failed.");  

  if (myRank == 0) {
    out << "Creating matrix" << endl;
  }

  // Create the matrix, using the above graph.
  RCP<CrsMatrixType> matrix (new CrsMatrixType (graph));

  if (myRank == 0) {
    out << "Setting all matrix entries to 1" << endl;
  }

  // Set all the owned entries to one.  Later we'll set nonlocal
  // entries' values in a loop.
  matrix->setAllToScalar (STS::one ());

  // Attempt to sum into nonowned entries (which nevertheless exist in
  // the matrix, just not on this process) using this process' rank.
  // The sumIntoGlobalValues() calls will record the data, but the
  // globalAssemble() method (called by fillComplete()) will silently
  // ignore entries whose columns are not in the column Map.  The
  // comment at the top of this test explains why this behavior is
  // reasonable.
  //
  // mfh 15,16 Dec 2012: Silently ignoring columns not in the column
  // Map has implications for the implementation of
  // sumIntoGlobalValues() for nonowned rows.  In particular, a
  // version of Map's getRemoteIDList() that uses one-sided
  // communication could invoke MPI_Get to figure out what the remote
  // process owns, without asking it or otherwise requiring
  // synchronization.  Thus, sumIntoGlobalValues() could throw
  // immediately on the calling process, rather than deferring the
  // exception to the remote process in globalAssemble().  If we
  // switch to that implementation, this unit test must be changed
  // accordingly.
  if (globalMinRow > rowMap->getMinAllGlobalIndex ()) {
    // Attempt to write to the (numLocalRows-1,numLocalCols-1) local entry of the previous process.
    matrix->sumIntoGlobalValues (globalMinRow-1, tuple (globalMaxCol), tuple (as<ST> (myRank)));
  }
  if (globalMaxRow < rowMap->getMaxAllGlobalIndex ()) {
    // Attempt to write to the (0,0) local entry of the next process.
    matrix->sumIntoGlobalValues (globalMaxRow+1, tuple (globalMinCol), tuple (as<ST> (myRank)));
  }

  if (myRank == 0) {
    out << "Calling fillComplete on the matrix" << endl;
  }
  TEST_NOTHROW(matrix->fillComplete (domainMap, rangeMap)); // Tpetra::Details::InvalidGlobalIndex<GO>

  // mfh 15 Dec 2012: We currently don't make promises about the state
  // of the matrix if fillComplete() throws.  Later, we might like to
  // improve the exception guarantees of fillComplete().  In that
  // case, the commented-out code below should be allowed to run.

  if (myRank == 0) {
    out << "Testing the matrix values" << endl;
  }

  // Test whether the entries have their correct values.
  bool localSuccess = true;
  std::ostringstream failMsg;
  {
    Array<GO> ind (2); // upper bound
    Array<ST> val (2); // upper bound

    for (GO globalRow = globalMinRow; globalRow <= globalMaxRow; ++globalRow) {
      size_t numEntries = 0; // output argument of below line.
      matrix->getGlobalRowCopy (globalRow, ind (), val (), numEntries);

      // Revise views based on numEntries.
      ArrayView<GO> indView = ind.view (0, numEntries);
      ArrayView<ST> valView = val.view (0, numEntries);

      // Sort the views jointly by column index.
      Tpetra::sort2 (indView.begin (), indView.end (), valView.begin ());

      if (numEntries != as<size_t> (1)) {
	localSuccess = false;
	failMsg << "Proc " << myRank << ": globalRow = " << globalRow << ": numEntries = " << numEntries << " != 1" << endl;
      }
      if (numEntries > 0 && indView[0] != globalRow) {
	localSuccess = false;
	failMsg << "Proc " << myRank << ": globalRow = " << globalRow << ": indView[0] = " << indView[0] << " != globalRow = " << globalRow << endl;
      }
      if (numEntries > 0 && valView[0] != STS::one ()) {
	localSuccess = false;
	failMsg << "Proc " << myRank << ": globalRow = " << globalRow << ": valView[0] = " << valView[0] << " != 1" << endl;
      }
    }
  }

  bool globalSuccess = true;
  {
    int globalSuccess_int = 1;
    reduceAll (*comm, Teuchos::REDUCE_MIN, localSuccess ? 1 : 0, outArg (globalSuccess_int));
    globalSuccess = (globalSuccess_int != 0);
  }

  if (! globalSuccess) {
    // Print out the failure messages on all processes.
    for (int p = 0; p < numProcs; ++p) {
      if (p == myRank) {
	out << failMsg.str () << endl;
	out << "Proc " << myRank << ": localSuccess = " << localSuccess << ", globalSuccess = " << globalSuccess << endl;
	//      std::flush (out);
      }
      // Do some barriers to allow output to finish.
      comm->barrier ();
      comm->barrier ();
      comm->barrier ();
    }
  }

  TEST_EQUALITY_CONST(globalSuccess, true);
}

//////////////////////////////////////////////////////////////////////
// INSTANTIATE THE TEMPLATED UNIT TESTS
//////////////////////////////////////////////////////////////////////

// mfh 16 Dec 2012: The #if 0 .. #endif section only worked if
// explicit instantiation was turned off.  See note in the comment
// above the test.
#if 0
//
// Instantiations for default Kokkos::Node type.
//

typedef Tpetra::CrsMatrix<float, int, int> mat_float_int_int_type;
typedef Tpetra::CrsMatrix<float, int, long> mat_float_int_long_type;

typedef Tpetra::CrsMatrix<double, int, int> mat_double_int_int_type;
typedef Tpetra::CrsMatrix<double, int, long> mat_double_int_long_type;

// Some tests are commented out to save time.  I've aimed for an
// orthogonal test plan over the variables (Scalar, GlobalOrdinal).

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( CrsMatrix, NonlocalSumInto_Ignore, mat_float_int_int_type )
// TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( CrsMatrix, NonlocalSumInto_Ignore, mat_float_int_long_type )

// TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( CrsMatrix, NonlocalSumInto_Ignore, mat_double_int_int_type )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( CrsMatrix, NonlocalSumInto_Ignore, mat_double_int_long_type )

#ifdef HAVE_TEUCHOS_COMPLEX

typedef Tpetra::CrsMatrix<std::complex<float>, int, int> mat_complex_float_int_int_type;
typedef Tpetra::CrsMatrix<std::complex<float>, int, long> mat_complex_float_int_long_type;

typedef Tpetra::CrsMatrix<std::complex<double>, int, int> mat_complex_double_int_int_type;
typedef Tpetra::CrsMatrix<std::complex<double>, int, long> mat_complex_double_int_long_type;

// TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( CrsMatrix, NonlocalSumInto_Ignore, mat_complex_float_int_int_type )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( CrsMatrix, NonlocalSumInto_Ignore, mat_complex_float_int_long_type )

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( CrsMatrix, NonlocalSumInto_Ignore, mat_complex_double_int_int_type )
// TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( CrsMatrix, NonlocalSumInto_Ignore, mat_complex_double_int_long_type )

#endif // HAVE_TEUCHOS_COMPLEX

#endif // 0


#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, NonlocalSumInto_Ignore, LO, GO, SCALAR, NODE )

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )
