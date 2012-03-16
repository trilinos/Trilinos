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

// This test demonstrates an Export to a CrsMatrix with a static
// graph.  In that case, all combine modes but INSERT are valid.

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_VerboseObject.hpp>

// Tpetra includes
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>


namespace {

  template<class CrsMatrixType>
  class Tester : public Teuchos::VerboseObject<Tester<CrsMatrixType> > {
  public:
    typedef CrsMatrixType matrix_type;
    typedef typename CrsMatrixType::scalar_type scalar_type;
    typedef typename CrsMatrixType::local_ordinal_type local_ordinal_type;
    typedef typename CrsMatrixType::global_ordinal_type global_ordinal_type;
    typedef typename CrsMatrixType::node_type node_type;
    typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;
    typedef Tpetra::CrsGraph<local_ordinal_type, global_ordinal_type, node_type> graph_type;

    Tester (const Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_DEFAULT,
	    const Teuchos::RCP<Teuchos::FancyOStream>& outStream = Teuchos::null) :
      Teuchos::VerboseObject<Tester<CrsMatrixType> > (verbLevel, outStream)
    {
      using Teuchos::FancyOStream;
      using Teuchos::getFancyOStream;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::rcpFromRef;

      RCP<FancyOStream> out = this->getOStream();
      if (this->getVerbLevel() == Teuchos::VERB_DEFAULT) {
	this->setVerbLevel (Teuchos::VERB_NONE); // run silently by default
      }
      if (out.is_null()) {
	if (this->getVerbLevel() == Teuchos::VERB_NONE) {
	  this->setOStream (getFancyOStream (rcp (new Teuchos::oblackholestream)));
	} 
	else {
	  this->setOStream (getFancyOStream (rcpFromRef (std::cout)));
	}
      }
    }

  private:
    typedef scalar_type ST;
    typedef local_ordinal_type LO;
    typedef global_ordinal_type GO;
    typedef node_type NT;

    Teuchos::RCP<const map_type> 
    makeNonoverlappingRowMap (const Teuchos::RCP<const Teuchos::Comm<int> >& comm, 
			      const Teuchos::RCP<node_type>& node,
			      const size_t localNumElts) const
    {
      using Tpetra::createContigMapWithNode;
      using Teuchos::RCP;
      using Teuchos::rcp;
      
      Tpetra::global_size_t globalNumElts = 
	Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
      return createContigMapWithNode<LO, GO, NT> (globalNumElts, localNumElts, comm, node);
    }

    Teuchos::RCP<const map_type> 
    makeOverlappingRowMap (const Teuchos::RCP<const map_type>& nonoverlapRowMap) const
    {
      using Tpetra::createNonContigMapWithNode;
      using Teuchos::Array;
      using Teuchos::as;
      using Teuchos::RCP;
      using Teuchos::rcp;
      typedef typename Array<GO>::size_type size_type;

      RCP<const Teuchos::Comm<int> > comm = nonoverlapRowMap->getComm ();
      RCP<NT> node = nonoverlapRowMap->getNode ();

      GO myMinGID = nonoverlapRowMap->getMaxGlobalIndex();
      GO myMaxGID = nonoverlapRowMap->getMaxGlobalIndex();

      if (myMaxGID != nonoverlapRowMap->getMaxAllGlobalIndex()) {
	++myMaxGID; // My process gets an extra row.
      }

      const size_type localNumElts = as<size_type> (myMaxGID - myMinGID + as<GO>(1));
      Array<GO> myGIDs (localNumElts);
      for (size_type k = 0; k < localNumElts; ++k) {
	myGIDs[k] = myMinGID + as<GO>(k);
      }
      return createNonContigMapWithNode<LO, GO, NT> (myGIDs(), comm, node);
    }

    Teuchos::RCP<const graph_type>
    makeNonoverlappingTestGraph (const Teuchos::RCP<const map_type>& nonoverlapRowMap) const
    {
      using Teuchos::Array;
      using Teuchos::ArrayView;
      using Teuchos::as;
      using Teuchos::RCP;
      using Teuchos::rcp;
      typedef typename ArrayView<GO>::size_type size_type;
	
      const size_t maxIndicesPerRow = 3;
      RCP<graph_type> graph = rcp (new graph_type (nonoverlapRowMap, maxIndicesPerRow));
      Array<GO> curIndices (maxIndicesPerRow);
      
      for (LO localRow = nonoverlapRowMap->getMinLocalIndex(); localRow <= nonoverlapRowMap->getMaxLocalIndex(); ++localRow) {
	const GO globalRow = nonoverlapRowMap->getGlobalElement (localRow);

	size_type numEntriesInCurRow = 3;
	GO startOffset = as<GO> (-1);
	if (globalRow == nonoverlapRowMap->getMinAllGlobalIndex()) {
	  numEntriesInCurRow = 2;
	  startOffset = as<GO> (0);
	} 
	else if (globalRow == nonoverlapRowMap->getMaxAllGlobalIndex()) {
	  numEntriesInCurRow = 2;
	}
	ArrayView<GO> curIndView = curIndices.view (0, numEntriesInCurRow);
	for (GO j = 0; j < as<GO> (numEntriesInCurRow); ++j) {
	  curIndView[j] = globalRow + startOffset;
	}
	graph->insertGlobalIndices (globalRow, curIndView);
      }
      
      graph->fillComplete();
      return graph;
    }


    Teuchos::RCP<const matrix_type>
    makeOverlappingSourceMatrix (const Teuchos::RCP<const map_type>& overlapRowMap) const
    {
      using Teuchos::Array;
      using Teuchos::ArrayView;
      using Teuchos::as;
      using Teuchos::RCP;
      using Teuchos::rcp;
      typedef typename ArrayView<GO>::size_type size_type;

      const ST diagVal = as<ST> (8);
      const ST leftVal = as<ST> (+1);
      const ST rightVal = as<ST> (-1);

      const size_t maxIndicesPerRow = 3;
      RCP<matrix_type> A_src = rcp (new matrix_type (overlapRowMap, maxIndicesPerRow));
      Array<GO> curIndices (maxIndicesPerRow);
      Array<ST> curValues (maxIndicesPerRow);

      curValues[0] = leftVal;
      curValues[1] = diagVal;
      curValues[2] = rightVal;
      
      for (LO localRow = overlapRowMap->getMinLocalIndex(); localRow <= overlapRowMap->getMaxLocalIndex(); ++localRow) {
	const GO globalRow = overlapRowMap->getGlobalElement (localRow);

	size_type numEntriesInCurRow = 3;
	GO startOffset = as<GO> (-1);
	ArrayView<ST> curValView;

	// Handle overlap: Don't insert entries for the "middle max"
	// row index.  "Middle max" means it's the max global index
	// for this process, but not the max over all processes.  We
	// assume that the overlap map is contiguous and overlaps by
	// one global index, so we resolve the ambiguity by lettng
	// the last owning process insert.
	if (globalRow != overlapRowMap->getMaxAllGlobalIndex() &&
	    globalRow == overlapRowMap->getMaxGlobalIndex()) {
	  continue;
	}
	// Insert entries in the current row.  Choose the column
	// indices and values to insert based on the current global
	// row index.
	if (globalRow == overlapRowMap->getMinAllGlobalIndex()) {
	  numEntriesInCurRow = 2;
	  startOffset = as<GO> (0);
	  curValView = curValues.view (1, numEntriesInCurRow);
	} 
	else if (globalRow == overlapRowMap->getMaxAllGlobalIndex()) {
	  numEntriesInCurRow = 2;
	  startOffset = as<GO> (-1);
	  curValView = curValues.view (0, numEntriesInCurRow);
	}
	else {
	  numEntriesInCurRow = 3;
	  startOffset = as<GO> (-1);
	  curValView = curValues.view (0, numEntriesInCurRow);
	}
	ArrayView<GO> curIndView = curIndices.view (0, numEntriesInCurRow);
	for (GO j = 0; j < as<GO> (numEntriesInCurRow); ++j) {
	  curIndView[j] = globalRow + startOffset + j;
	}
	A_src->insertGlobalValues (globalRow, curIndView, curValView);
      }
      
      A_src->fillComplete();
      return A_src;
    }
    
    /// \brief Whether A_src == A_tgt, locally.
    ///
    /// \pre The row Maps of A_src and A_tgt share the same GIDs,
    ///   except that the row Map of A_src has overlap, and the row
    ///   Map of A_tgt does not.  This fits our Export test protocol.
    /// 
    /// \pre Both A_src and A_tgt are locally indexed.  Test for this
    ///   via \c isLocallyIndexed().
    ///
    /// \return whether A_src == A_tgt locally, and if not, the global index of the first nonequal row.
    std::pair<bool, global_ordinal_type>
    locallyEqual (const Teuchos::RCP<const matrix_type>& A_src,
		  const Teuchos::RCP<const matrix_type>& A_tgt) const
    {
      using Teuchos::as;
      using Teuchos::Array;
      using Teuchos::ArrayView;
      using Teuchos::OrdinalTraits;
      using Teuchos::RCP;
      
      TEUCHOS_TEST_FOR_EXCEPTION(! A_src->isLocallyIndexed() || ! A_tgt->isLocallyIndexed(),
        std::invalid_argument, "Both A_src and A_tgt must be locally indexed.");

      // const size_t srcMaxNumRowEntries = A_src->getNodeMaxNumRowEntries ();
      // const size_t tgtMaxNumRowEntries = A_tgt->getNodeMaxNumRowEntries ();
      ArrayView<const LO> srcIndView, tgtIndView;
      ArrayView<const ST> srcValView, tgtValView;

      // We assume that the row Maps of A_src and A_tgt share the same
      // GIDs, except that the row Map of A_src has overlap, and the
      // row Map of A_tgt does not.  This means we can iterate on each
      // process through the rows owned by this process in A_tgt's row
      // Map.
      RCP<const map_type> srcRowMap = A_src->getRowMap();
      RCP<const map_type> tgtRowMap = A_tgt->getRowMap();
      for (LO tgtLocalRow = tgtRowMap->getMinLocalIndex(); 
	   tgtLocalRow <= tgtRowMap->getMaxLocalIndex(); ++tgtLocalRow) {
	// There's no particular reason why A_src should have the same
	// local indices as A_tgt, so convert back and forth.
	const GO globalRow = tgtRowMap->getGlobalElement (tgtLocalRow);
	TEUCHOS_TEST_FOR_EXCEPTION(globalRow == OrdinalTraits<GO>::invalid(), 
          std::logic_error, "Local row index " << tgtLocalRow 
          << " of A_tgt has no global index.");
	const LO srcLocalRow = srcRowMap->getLocalElement (globalRow);
	TEUCHOS_TEST_FOR_EXCEPTION(srcLocalRow == OrdinalTraits<LO>::invalid(), 
          std::logic_error, "Global row index " << globalRow 
          << " has no local index of A_src.");

	const size_t srcNumEntries = A_src->getNumEntriesInLocalRow (srcLocalRow);
	const size_t tgtNumEntries = A_tgt->getNumEntriesInLocalRow (tgtLocalRow);
	if (srcNumEntries != tgtNumEntries) {
	  return std::make_pair (false, globalRow);
	}
	A_src->getLocalRowView (srcLocalRow, srcIndView, srcValView);
	A_tgt->getLocalRowView (tgtLocalRow, tgtIndView, tgtValView);

	// Assume for now that the entries are sorted by column index.
	if (! std::equal (srcIndView.begin(), srcIndView.end(), tgtIndView.begin())) {
	  return std::make_pair (false, globalRow);
	}
	// FIXME (mfh 15 Mar 2012) Should we include a small error
	// tolerance here for roundoff?
	if (! std::equal (srcValView.begin(), srcValView.end(), tgtValView.begin())) {
	  return std::make_pair (false, globalRow);
	}
      } // for each local row index of the target matrix on my process
      
      return std::make_pair (true, as<global_ordinal_type> (0));
    }

    /// \brief Whether A_src == A_tgt, globally.
    ///
    /// \pre The row Maps of A_src and A_tgt share the same GIDs,
    ///   except that the row Map of A_src has overlap, and the row
    ///   Map of A_tgt does not.  This fits our Export test protocol.
    /// 
    /// \pre Both A_src and A_tgt are locally indexed.  Test for this
    ///   via \c isLocallyIndexed().
    ///
    /// \return whether A_src == A_tgt globally, and if not, the global index of the first nonequal row.
    std::pair<bool, global_ordinal_type>
    globallyEqual (const Teuchos::RCP<const matrix_type>& A_src,
		   const Teuchos::RCP<const matrix_type>& A_tgt) const
    {
      using Teuchos::as;
      using Teuchos::Comm;
      using Teuchos::ptr;
      using Teuchos::RCP;
      using Teuchos::reduceAll;

      RCP<const Comm<int> > comm = A_src->getComm();
      
      std::pair<bool, GO> localResult = locallyEqual (A_src, A_tgt);
      // reduceAll doesn't work for bool, so convert to int first.
      const int equalLocally = localResult.first ? 1 : 0;
      GO firstNonequalIndex = as<GO> (0);

      int equalGlobally = 1;
      reduceAll (*comm, Teuchos::REDUCE_MIN, equalLocally, ptr (&equalGlobally));

      if (equalGlobally == 0) {
	reduceAll (*comm, Teuchos::REDUCE_MIN, localResult.second, ptr (&firstNonequalIndex));
      }
      return std::make_pair (equalGlobally == 1, firstNonequalIndex);
    }

  public:
    bool 
    testExportToCrsMatrixWithStaticGraph (const Teuchos::RCP<const Teuchos::Comm<int> >& comm, 
					  const Teuchos::RCP<node_type>& node,
					  const size_t localNumElts) const
    {
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::FancyOStream;
      using Teuchos::includesVerbLevel;
      using Teuchos::OSTab;
      using std::endl;
      typedef Tpetra::Export<LO, GO, NT> export_type;

      RCP<FancyOStream> out = this->getOStream();
      const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();

      if (includesVerbLevel (verbLevel, Teuchos::VERB_LOW)) {
	*out << "Test Export to Tpetra::CrsMatrix target with static graph" << endl;
      }
      OSTab tab = this->getOSTab ();

      if (includesVerbLevel (verbLevel, Teuchos::VERB_LOW)) {
	*out << "Making nonoverlapping and overlapping row Maps" << endl;
      }
      RCP<const map_type> nonoverlapRowMap = 
	makeNonoverlappingRowMap (comm, node, localNumElts);
      RCP<const map_type> overlapRowMap = 
	makeOverlappingRowMap (nonoverlapRowMap);

      if (includesVerbLevel (verbLevel, Teuchos::VERB_LOW)) {
	*out << "Making source matrix" << endl;
      }
      RCP<const matrix_type> A_src = makeOverlappingSourceMatrix (overlapRowMap);

      if (includesVerbLevel (verbLevel, Teuchos::VERB_LOW)) {
	*out << "Making target graph and matrix" << endl;
      }
      RCP<const graph_type> graph_tgt = makeNonoverlappingTestGraph (nonoverlapRowMap);
      RCP<matrix_type> A_tgt = rcp (new matrix_type (graph_tgt));

      if (includesVerbLevel (verbLevel, Teuchos::VERB_LOW)) {
	*out << "Making Export instance" << endl;
      }
      export_type exporter (A_src->getRowMap(), A_tgt->getRowMap());
      if (includesVerbLevel (verbLevel, Teuchos::VERB_LOW)) {
	*out << "Doing Export from source to target matrix" << endl;
      }
      A_tgt->doExport (*A_src, exporter, Tpetra::ADD);
      if (includesVerbLevel (verbLevel, Teuchos::VERB_LOW)) {
	*out << "Calling fillComplete() on target matrix" << endl;
      }
      A_tgt->fillComplete();

      if (includesVerbLevel (verbLevel, Teuchos::VERB_LOW)) {
	*out << "Comparing source and target matrix" << endl;
      }
      std::pair<bool, GO> result = globallyEqual (A_src, A_tgt);
      if (includesVerbLevel (verbLevel, Teuchos::VERB_LOW)) {
	if (result.first) {
	  *out << "Source and target matrix are the same." << endl;
	}
	else {
	  *out << "Source and target matrix differ at row(s) >= row " 
	       << result.second << "." << endl;
	} 
      }

      return result.first;
    }
  };
} // namespace (anonymous)


int 
main (int argc, char *argv[])
{
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::CommandLineProcessor;
  using Teuchos::EVerbosityLevel;
  using Teuchos::FancyOStream;
  using Teuchos::getFancyOStream;
  using Teuchos::includesVerbLevel;
  using Teuchos::OSTab;
  using Teuchos::ptr;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using std::cout;
  using std::endl;
  
  typedef double ST;
  typedef int LO;
  typedef long GO;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType NT;
  typedef Tpetra::CrsMatrix<ST, LO, GO, NT> matrix_type;

  bool success = true; // May be changed by tests

  RCP<Teuchos::oblackholestream> blackHole = rcp (new Teuchos::oblackholestream);
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, blackHole.getRawPtr());
  RCP<const Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  RCP<NT> node = Tpetra::DefaultPlatform::getDefaultPlatform().getNode();
  const int numProcs = comm->getSize();
  const int myRank = comm->getRank();

  RCP<FancyOStream> out = getFancyOStream (rcpFromRef (std::cout));
  RCP<FancyOStream> procZeroOut = (myRank == 0) ? 
    getFancyOStream (rcpFromRef (std::cout)) : 
    getFancyOStream (blackHole);

  //
  // Default values of command-line options.
  //
  bool verbose = false;
  CommandLineProcessor cmdp (false, true);

  //
  // Set command-line options.
  //
  cmdp.setOption ("verbose", "quiet", &verbose, "Print verbose output.");
  // Parse command-line options.
  if (cmdp.parse (argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
    *procZeroOut << "End Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }

  const EVerbosityLevel verbLevel = 
    verbose ? Teuchos::VERB_EXTREME : Teuchos::VERB_NONE;

  if (includesVerbLevel (verbLevel, Teuchos::VERB_LOW)) {
    *out << "Running test on " << numProcs << " process" 
	 << (numProcs != 1 ? "es" : "") << "." << endl;
  }

  typedef Tpetra::CrsMatrix<ST, LO, GO, NT> matrix_type;
  Tester<matrix_type> tester (verbLevel, out);
  tester.testExportToCrsMatrixWithStaticGraph (comm, node, 10);

  *procZeroOut << "End Result: TEST " << (success ? "PASSED" : "FAILED") << endl;
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
