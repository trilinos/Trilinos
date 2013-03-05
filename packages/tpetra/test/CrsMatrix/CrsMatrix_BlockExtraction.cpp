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

#include <numeric>
#include <algorithm>

#include <Tpetra_TestingUtilities.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <TpetraExt_BlockExtraction.hpp>

namespace {

  using Tpetra::TestingUtilities::getNode;
  using Tpetra::TestingUtilities::getDefaultComm;

  using std::string;
  using Teuchos::RCP;
  using Teuchos::null;
  using Teuchos::rcp;
  using Teuchos::Comm;
  using Teuchos::ScalarTraits;
  using Teuchos::OrdinalTraits;
  using Teuchos::ArrayRCP;
  using Teuchos::tuple;
  using Tpetra::CrsMatrix;
  using Tpetra::RowMatrix;
  using Tpetra::global_size_t;

  // string filedir;
  // double errorTolSlack = 1e+1;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    // clp.setOption(
    //     "filedir",&filedir,"Directory of expected matrix files.");
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &Tpetra::TestingUtilities::testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
    // clp.setOption(
    //     "error-tol-slack", &errorTolSlack,
    //     "Slack off of machine epsilon used to check test results" );
  }


  //
  // UNIT TESTS
  //

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockDiagonalExtraction, RuntimeExceptions, LO, GO, Scalar, Node )
  {
    typedef Tpetra::Map<LO,GO,Node> Map;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    RCP<Node> node = getNode<Node>();
    // set the block sizes
    // note one block of zero size, to test capability
    Teuchos::Tuple<int,7> block_sizes = Teuchos::tuple<int>(1,3,5,0,5,3,1) ;
    const int maxBlockSize = *std::max_element( block_sizes.begin(), block_sizes.end() );
    // create a Map
    const size_t numLocal = std::accumulate( block_sizes.begin(), block_sizes.end(), (size_t)0 );
    RCP<const Map> map = Tpetra::createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm,node);
    RCP<const RowMatrix<Scalar,LO,GO,Node> > mat;
    {
      RCP<CrsMatrix<Scalar,LO,GO,Node> > mat_crs = Tpetra::createCrsMatrix<Scalar>( map );
      for (GO gid=map->getMinGlobalIndex(); gid <= map->getMaxGlobalIndex(); ++gid) {
        // add diagonal entries
        mat_crs->insertGlobalValues( gid, tuple<GO>(gid), tuple<Scalar>(1.0) );
        // add some entries outside of the diagonal block
        if (gid >= map->getMinGlobalIndex() + maxBlockSize) mat_crs->insertGlobalValues( gid, tuple<GO>(gid - maxBlockSize), tuple<Scalar>(1.0) );
        if (gid + maxBlockSize <= map->getMaxGlobalIndex()) mat_crs->insertGlobalValues( gid, tuple<GO>(gid + maxBlockSize), tuple<Scalar>(1.0) );
      }
      mat_crs->fillComplete();
      mat = mat_crs;
    }
    //
    // bad block offsets
    //
    // block sizes too small (corresponds to latter first_points too small)
    {
      Teuchos::Tuple<int,7> bad_bsizes = Teuchos::tuple<int>(1, /* BAD */ 1,5,0,5,3,1) ;
      Teuchos::Array<LO> bad_bfirsts( bad_bsizes.size()+1 );
      bad_bfirsts[0] = 0;
      for (int i=0; i < (int)bad_bsizes.size(); ++i) {
        bad_bfirsts[i+1] = bad_bfirsts[i] + bad_bsizes[i];
      }
      Teuchos::ArrayRCP<Scalar> out_diags;
      Teuchos::ArrayRCP<LO>     out_offsets;
      TEST_THROW( Tpetra::Ext::extractBlockDiagonals( *mat, bad_bfirsts().getConst(), out_diags, out_offsets ) , std::runtime_error );
    }
    // block sizes too large (corresponds to latter first_points too large)
    {
      Teuchos::Tuple<int,7> bad_bsizes = Teuchos::tuple<int>(1, /* BAD */ 5 ,5,0,5,3,1) ;
      Teuchos::Array<LO> bad_bfirsts( bad_bsizes.size()+1 );
      bad_bfirsts[0] = 0;
      for (int i=0; i < (int)bad_bsizes.size(); ++i) {
        bad_bfirsts[i+1] = bad_bfirsts[i] + bad_bsizes[i];
      }
      Teuchos::ArrayRCP<Scalar> out_diags;
      Teuchos::ArrayRCP<LO>     out_offsets;
      TEST_THROW( Tpetra::Ext::extractBlockDiagonals( *mat, bad_bfirsts().getConst(), out_diags, out_offsets ) , std::runtime_error );
    }
    // negative block size (corresponds to non-monotonically increasing first_points)
    {
      Teuchos::Tuple<int,7> bad_bsizes = Teuchos::tuple<int>(1,3,5, /* BAD */ -1 ,5,3,1) ;
      Teuchos::Array<LO> bad_bfirsts( bad_bsizes.size()+1 );
      bad_bfirsts[0] = 0;
      for (int i=0; i < (int)bad_bsizes.size(); ++i) {
        bad_bfirsts[i+1] = bad_bfirsts[i] + bad_bsizes[i];
      }
      Teuchos::ArrayRCP<Scalar> out_diags;
      Teuchos::ArrayRCP<LO>     out_offsets;
      TEST_THROW( Tpetra::Ext::extractBlockDiagonals( *mat, bad_bfirsts().getConst(), out_diags, out_offsets ) , std::runtime_error );
    }
    // first first_point required to be zero
    {
      Teuchos::Tuple<int,7> bad_bsizes = Teuchos::tuple<int>(1,3,5,0,5,3,1) ;
      Teuchos::Array<LO> bad_bfirsts( bad_bsizes.size()+1 );
      bad_bfirsts[0] = 1; /* BAD */
      for (int i=0; i < (int)bad_bsizes.size(); ++i) {
        bad_bfirsts[i+1] = bad_bfirsts[i] + bad_bsizes[i];
      }
      Teuchos::ArrayRCP<Scalar> out_diags;
      Teuchos::ArrayRCP<LO>     out_offsets;
      TEST_THROW( Tpetra::Ext::extractBlockDiagonals( *mat, bad_bfirsts().getConst(), out_diags, out_offsets ) , std::runtime_error );
    }
    // matrix is required to be fillComplete()
    {
      Teuchos::Tuple<int,7> bad_bsizes = Teuchos::tuple<int>(1,3,5,0,5,3,1) ;
      Teuchos::Array<LO> bfirsts( bad_bsizes.size()+1 );
      bfirsts[0] = 0;
      for (int i=0; i < (int)bad_bsizes.size(); ++i) {
        bfirsts[i+1] = bfirsts[i] + bad_bsizes[i];
      }
      Teuchos::ArrayRCP<Scalar> out_diags;
      Teuchos::ArrayRCP<LO>     out_offsets;
      RCP<CrsMatrix<Scalar,LO,GO,Node> > not_fill_complete = Tpetra::createCrsMatrix<Scalar>( map );
      TEST_THROW( Tpetra::Ext::extractBlockDiagonals( *not_fill_complete, bfirsts().getConst(), out_diags, out_offsets ) , std::runtime_error );
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockDiagonalExtraction, SimpleExtraction, LO, GO, Scalar, Node )
  {
    using Teuchos::as;
    typedef Tpetra::Map<LO,GO,Node>           Map;
    typedef Tpetra::BlockMap<LO,GO,Node> BlockMap;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    RCP<Node> node = getNode<Node>();

    //
    // set the block sizes
    // note one block of zero size, to test capability
    Teuchos::Array<int> block_sizes = Teuchos::tuple<int>(1,3,5,0,5,3,1);
    const int maxBlockSize = *std::max_element( block_sizes.begin(), block_sizes.end() );
    const size_t expected_alloc_size = std::inner_product( block_sizes.begin(), block_sizes.end(), block_sizes.begin(), 0 );
    //
    // create a point Map
    //
    const size_t numLocal = std::accumulate( block_sizes.begin(), block_sizes.end(), as<size_t>(0) );
    RCP<const Map> map = Tpetra::createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm,node);
    //
    // fill matrix for testing
    //
    RCP<const RowMatrix<Scalar,LO,GO,Node> > mat;
    {
      RCP<CrsMatrix<Scalar,LO,GO,Node> > mat_crs = Tpetra::createCrsMatrix<Scalar>( map );
      for (GO gid=map->getMinGlobalIndex(); gid <= map->getMaxGlobalIndex(); ++gid) {
        // add diagonal entries
        mat_crs->insertGlobalValues( gid, tuple<GO>(gid), tuple<Scalar>(1.0) );
        out << gid << " ";
        // add some entries outside of the diagonal block
        if (gid >= map->getMinGlobalIndex() + maxBlockSize) {mat_crs->insertGlobalValues( gid, tuple<GO>(gid - maxBlockSize), tuple<Scalar>(1.0) ); out << gid - maxBlockSize << " ";}
        if (gid + maxBlockSize <= map->getMaxGlobalIndex()) {mat_crs->insertGlobalValues( gid, tuple<GO>(gid + maxBlockSize), tuple<Scalar>(1.0) ); out << gid + maxBlockSize << " ";}
      }
      mat_crs->fillComplete();
      mat = mat_crs;
    }

    //
    // create block_firsts for first extraction
    //
    Teuchos::Array<LO> block_firsts( block_sizes.size()+1 );
    block_firsts[0] = 0;
    for (int i=0; i < (int)block_sizes.size(); ++i) {
      block_firsts[i+1] = block_firsts[i] + block_sizes[i];
    }
    //
    // perform first extraction
    //
    Teuchos::ArrayRCP<Scalar> block_diagonals1;
    Teuchos::ArrayRCP<LO>     block_offsets1;
    Tpetra::Ext::extractBlockDiagonals<Scalar,LO,GO,Node>( *mat, block_firsts(), block_diagonals1, block_offsets1 );
    //
    // independently test first extraction
    //
    {
      TEST_EQUALITY( (size_t)expected_alloc_size, (size_t)block_diagonals1.size() );
      TEST_EQUALITY( (size_t)block_sizes.size(), (size_t)block_offsets1.size() );
      const int num_zeros_extracted    = (int)std::count( block_diagonals1.begin(), block_diagonals1.end(), ScalarTraits<Scalar>::zero() );
      const int num_nonzeros_extracted = (int)block_diagonals1.size() - num_zeros_extracted;
      TEST_EQUALITY( num_nonzeros_extracted, (int)mat->getNodeNumDiags() );
      TEST_EQUALITY_CONST( num_nonzeros_extracted < (int)mat->getNodeNumEntries(), true );
    }

    //
    // create a BlockMap for use in second extraction
    //
    Teuchos::Tuple<GO,7> globalBlockIDs = Teuchos::tuple<GO>(1,2,3,4,5,6,7) ;
    RCP<const BlockMap> bmap = rcp(new BlockMap(map,globalBlockIDs,block_sizes,map->getNode()));
    //
    // perform second extraction
    //
    Teuchos::ArrayRCP<Scalar> block_diagonals2;
    Teuchos::ArrayRCP<LO>     block_offsets2;
    Tpetra::Ext::extractBlockDiagonals<Scalar,LO,GO,Node>( *mat, *bmap, block_diagonals2, block_offsets2 );
    //
    // independently test second extraction
    //
    {
      TEST_EQUALITY( (size_t)expected_alloc_size, (size_t)block_diagonals2.size() );
      TEST_EQUALITY( (size_t)block_sizes.size(), (size_t)block_offsets2.size() );
      const int num_zeros_extracted    = (int)std::count( block_diagonals2.begin(), block_diagonals2.end(), ScalarTraits<Scalar>::zero() );
      const int num_nonzeros_extracted = (int)block_diagonals2.size() - num_zeros_extracted;
      TEST_EQUALITY( num_nonzeros_extracted, (int)mat->getNodeNumDiags() );
      TEST_EQUALITY_CONST( num_nonzeros_extracted < (int)mat->getNodeNumEntries(), true );
    }

    //
    // compare first extraction against second
    //
    TEST_COMPARE_ARRAYS( block_diagonals1, block_diagonals2 );
    TEST_COMPARE_ARRAYS( block_offsets1,   block_offsets2 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockRowExtraction, DiagonalExtraction, LO, GO, Scalar, Node )
  {
    typedef Tpetra::Map<LO,GO,Node>           Map;
    typedef Tpetra::BlockMap<LO,GO,Node> BlockMap;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    //
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    RCP<Node> node = getNode<Node>();
    const int myImageID = comm->getRank();
    //
    // set the block sizes
    // try to hit the border cases: some zero block at the outside and inside,
    Teuchos::Tuple<int,6> block_sizes = Teuchos::tuple<int>( myImageID%2 , 2 , 0 , myImageID+1 , 3, (myImageID+1)%2 );
    const int numBlocks = (int)block_sizes.size();
    const int maxBlockSize = *std::max_element( block_sizes.begin(), block_sizes.end() );
    //
    // create a point Map
    //
    const size_t numLocal = std::accumulate( block_sizes.begin(), block_sizes.end(), 0 );
    RCP<const Map> map = Tpetra::createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm,node);
    //
    // fill matrix for testing
    //
    //
    RCP<const RowMatrix<Scalar,LO,GO,Node> > mat;
    {
      RCP<CrsMatrix<Scalar,LO,GO,Node> > mat_crs = Tpetra::createCrsMatrix<Scalar>( map );
      for (GO gid=map->getMinGlobalIndex(); gid <= map->getMaxGlobalIndex(); ++gid) {
        // add diagonal entries
        mat_crs->insertGlobalValues( gid, tuple<GO>(gid), tuple<Scalar>(1.0) );
        // add some entries outside of the diagonal block; max
        if (gid >= map->getMinGlobalIndex() + maxBlockSize) mat_crs->insertGlobalValues( gid, tuple<GO>(gid - maxBlockSize), tuple<Scalar>(1.0) );
        if (gid + maxBlockSize <= map->getMaxGlobalIndex()) mat_crs->insertGlobalValues( gid, tuple<GO>(gid + maxBlockSize), tuple<Scalar>(1.0) );
      }
      mat_crs->fillComplete();
      mat = mat_crs;
    }

    //
    // create a BlockMap for the row and column partitioning (since we will be testing against the block diagonls
    //
    Teuchos::Tuple<GO,6> globalBlockIDs = Teuchos::tuple<GO>(1,2,3,4,5,6);
    RCP<const BlockMap> bmap = rcp(new BlockMap(map,globalBlockIDs,block_sizes,map->getNode()));

    //
    // perform block diagonal extraction
    //
    Teuchos::ArrayRCP<Scalar> block_diagonals;
    Teuchos::ArrayRCP<LO>     block_offsets;
    Tpetra::Ext::extractBlockDiagonals<Scalar,LO,GO,Node>( *mat, *bmap, block_diagonals, block_offsets );

    //
    // perform block row extractions
    //
    int total_num_nonzeros_extracted = 0;
    for (int b=0; b < numBlocks; ++b)
    {
      Teuchos::ArrayRCP< Teuchos::ArrayRCP<Scalar> > block_entries;
      Teuchos::ArrayRCP< LO >                        block_indices;
      Tpetra::Ext::extractBlockRow(b, *mat, *bmap, *bmap, block_entries, block_indices);
      TEST_EQUALITY( block_entries.size(), block_indices.size() );
      if (block_sizes[b] == 0) {
        // trivial block, nothing to extract
        TEST_EQUALITY_CONST( block_entries.size() == 0, true );
      }
      else {
        TEST_EQUALITY_CONST( block_entries.size() > 0, true );
        // find the diagonal, compare it against block_diagonals[b]
        bool diagFound = false;
        for (Teuchos_Ordinal jj=0; jj < block_entries.size(); ++jj)
        {
          // block row partitioning is block column partitioning for this test
          // therefore, we don't need to compare global block IDs; local block IDs will suffice
          if ( block_indices[jj] == b )
          {
            // can't find the diagonal block twice
            TEST_EQUALITY_CONST( diagFound, false );
            TEST_COMPARE_ARRAYS( block_entries[jj], block_diagonals(block_offsets[b],block_sizes[b]*block_sizes[b]) );
            diagFound = true;
          }
          // count the number of non-zeros
          const int num_zeros_extracted = (int)std::count( block_entries[jj].begin(), block_entries[jj].end(), ScalarTraits<Scalar>::zero() );
          total_num_nonzeros_extracted += (int)block_entries[jj].size() - num_zeros_extracted;
        }
        TEST_EQUALITY_CONST( diagFound, true );
      }
    }
    // should have extracted all non-zeros, no more or less
    TEST_EQUALITY( total_num_nonzeros_extracted , (int)mat->getNodeNumEntries() );
  }


  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP(SCALAR, LO, GO, NODE) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockDiagonalExtraction, SimpleExtraction, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockDiagonalExtraction, RuntimeExceptions, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockRowExtraction,      DiagonalExtraction, LO, GO, SCALAR, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN_NOGPU( UNIT_TEST_GROUP )

}

