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

#include <Tpetra_TestingUtilities.hpp>

#include <Tpetra_BlockCrsGraph.hpp>

namespace {

  using Tpetra::TestingUtilities::getNode;
  using Tpetra::TestingUtilities::getDefaultComm;

  using Teuchos::as;
  using Teuchos::null;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Teuchos::arcpClone;
  using Tpetra::BlockMap;
  using Tpetra::OptimizeOption;
  using Tpetra::DoOptimizeStorage;
  using Tpetra::DoNotOptimizeStorage;
  using Tpetra::DefaultPlatform;
  using Tpetra::global_size_t;
  using std::sort;
  using Teuchos::arrayView;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::Comm;
  using std::endl;
  using std::swap;
  using std::min;
  using std::max;
  using Teuchos::Array;
  using Teuchos::TypeTraits::is_same;
  using Teuchos::ArrayView;
  using Tpetra::BlockCrsGraph;
  using Tpetra::global_size_t;
  using Teuchos::arcp;
  using std::string;
  using std::unique;
  using Teuchos::tuple;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;
  using Tpetra::ProfileType;
  using Tpetra::StaticProfile;
  using Tpetra::DynamicProfile;
  using Tpetra::Array_size_type;

  double errorTolSlack = 1e+1;
  string filedir;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.setOption(
        "filedir",&filedir,"Directory of expected input files.");
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &Tpetra::TestingUtilities::testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }

  //
  // UNIT TESTS
  //

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( BlockCrsGraph, ColMap1, LO, GO, Node )
  {
    //This test fills a (block-tri-diagonal) block-crs-graph such that in parallel
    //the column-map should have an overlapping set of entries (i.e.,
    //different than the row-map), and verifies that the column-map is correct.

    RCP<Node> node = getNode<Node>();
    typedef BlockCrsGraph<LO,GO,Node> BGRAPH;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    GO indexBase = 0;
    // create a Map
    const size_t numLocalBlocks = 2;
    GO maxGlobalBlock = numLocalBlocks*comm->getSize();
    const LO blockSize = 2;
    const size_t maxEntriesPerRow = 3;
    RCP<BlockMap<LO,GO,Node> > rowmap = rcp( new BlockMap<LO,GO,Node>(INVALID,numLocalBlocks,blockSize,indexBase,comm,node) );

    //now set up the list of block-column-ids that we expect the
    //column-map to contain after fillComplete:
    size_t numLocalColBlocks = numLocalBlocks;
    if (comm->getRank() != 0) ++numLocalColBlocks;
    if (comm->getRank() != comm->getSize()-1) ++numLocalColBlocks;
    Array<GO> blockColIDs(numLocalColBlocks);
    typedef typename Array<GO>::size_type Tsize_t;
    Teuchos::ArrayView<const GO> blk_rows = rowmap->getNodeBlockIDs();
    GO first_row = blk_rows[0];
    Tsize_t offset = 0;
    if (comm->getRank() != 0) {
      blockColIDs[offset++] = first_row - 1;
    }
    GO last_row = 0;
    for(LO i=0; i<blk_rows.size(); ++i) {
      blockColIDs[offset++] = blk_rows[i];
      last_row = blk_rows[i];
    }
    if (offset < blockColIDs.size()) blockColIDs[offset++] = last_row + 1;

    // create the graph
    RCP<BGRAPH> bgrph = rcp( new BGRAPH(rowmap,maxEntriesPerRow,DynamicProfile) );
    for(int i=0; i<blk_rows.size(); ++i) {
      GO row = blk_rows[i];
      {
        GO col = row;
        bgrph->insertGlobalIndices(row, Teuchos::arrayView(&col, 1));
      }
      if (row > indexBase) {
        GO col = row - 1;
        bgrph->insertGlobalIndices(row, Teuchos::arrayView(&col, 1));
      }
      if (row < maxGlobalBlock-1) {
        GO col = row + 1;
        bgrph->insertGlobalIndices(row, Teuchos::arrayView(&col, 1));
      }
    }

    bgrph->fillComplete();
    RCP<const BlockMap<LO,GO,Node> > colmap = bgrph->getBlockColMap();
    ArrayView<const GO> blk_cols = colmap->getNodeBlockIDs();
    TEST_EQUALITY(blk_cols.size(), blockColIDs.size());
    TEST_COMPARE_ARRAYS(blk_cols, blockColIDs() );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( BlockCrsGraph, Queries, LO, GO, Node )
  {
    //This test fills a (block-tri-diagonal) block-crs-graph such that
    //in parallel the column-map should have an overlapping set of
    //entries (i.e., different than the row-map), and verifies that
    //some queries work for attributes such as row-lengths, etc.

    RCP<Node> node = getNode<Node>();
    typedef BlockCrsGraph<LO,GO,Node> BGRAPH;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    GO indexBase = 0;
    // create a Map
    const size_t numLocalBlocks = 2;
    GO maxGlobalBlock = numLocalBlocks*comm->getSize();
    const LO blockSize = 2;
    const size_t maxEntriesPerRow = 3;
    RCP<BlockMap<LO,GO,Node> > rowmap = rcp( new BlockMap<LO,GO,Node>(INVALID,numLocalBlocks,blockSize,indexBase,comm,node) );

    //now set up the list of block-column-ids that we expect the
    //column-map to contain after fillComplete:
    size_t numLocalColBlocks = numLocalBlocks;
    if (comm->getRank() != 0) ++numLocalColBlocks;
    if (comm->getRank() != comm->getSize()-1) ++numLocalColBlocks;
    Array<GO> blockColIDs(numLocalColBlocks);
    typedef typename Array<GO>::size_type Tsize_t;
    Teuchos::ArrayView<const GO> blk_rows = rowmap->getNodeBlockIDs();
    GO first_row = blk_rows[0];
    Tsize_t offset = 0;
    if (comm->getRank() != 0) {
      blockColIDs[offset++] = first_row - 1;
    }
    GO last_row = 0;
    for(LO i=0; i<blk_rows.size(); ++i) {
      blockColIDs[offset++] = blk_rows[i];
      last_row = blk_rows[i];
    }
    if (offset < blockColIDs.size()) blockColIDs[offset++] = last_row + 1;

    // create the graph
    RCP<BGRAPH> bgrph = rcp( new BGRAPH(rowmap,maxEntriesPerRow,DynamicProfile) );
    for(int i=0; i<blk_rows.size(); ++i) {
      GO row = blk_rows[i];
      {
        GO col = row;
        bgrph->insertGlobalIndices(row, Teuchos::arrayView(&col, 1));
      }
      if (row > indexBase) {
        GO col = row - 1;
        bgrph->insertGlobalIndices(row, Teuchos::arrayView(&col, 1));
      }
      if (row < maxGlobalBlock-1) {
        GO col = row + 1;
        bgrph->insertGlobalIndices(row, Teuchos::arrayView(&col, 1));
      }
    }

    bgrph->fillComplete();
    RCP<const BlockMap<LO,GO,Node> > colmap = bgrph->getBlockColMap();
    ArrayView<const GO> blk_cols = colmap->getNodeBlockIDs();
    TEST_EQUALITY(blk_cols.size(), blockColIDs.size());
    TEST_COMPARE_ARRAYS(blk_cols, blockColIDs() );

    size_t map_blk_elems = blk_rows.size();
    TEST_EQUALITY(bgrph->getNodeNumBlockRows(), map_blk_elems);
    TEST_EQUALITY(bgrph->getGlobalNumBlockRows(), rowmap->getGlobalNumBlocks());

    for(int i=0; i<blk_rows.size(); ++i) {
      GO row = blk_rows[i];
      if (row > indexBase && row < maxGlobalBlock-1) {
        size_t row_len = bgrph->getGlobalBlockRowLength(row);
        size_t expected_row_len = 3;
        TEST_EQUALITY(row_len, expected_row_len);
      }
      else {
        size_t row_len = bgrph->getGlobalBlockRowLength(row);
        size_t expected_row_len = 2;
        TEST_EQUALITY(row_len, expected_row_len);
      }
    }
  }

  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP( LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( BlockCrsGraph, ColMap1  , LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( BlockCrsGraph, Queries  , LO, GO, NODE )

    TPETRA_ETI_MANGLING_TYPEDEFS()

    TPETRA_INSTANTIATE_LGN_NOGPU( UNIT_TEST_GROUP )
}
