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

#include <Tpetra_BlockMultiVector.hpp>
#include <Tpetra_Import.hpp>
#include <Tpetra_BlockCrsGraph.hpp>
#include <Tpetra_VbrMatrix.hpp>

// TODO: add test where some nodes have zero rows
// TODO: add test where non-"zero" graph is used to build matrix; if no values are added to matrix, the operator effect should be zero. This tests that matrix values are initialized properly.
// TODO: add test where dynamic profile initially has no allocation, then entries are added. this will test new view functionality.

namespace Teuchos {
  template <>
    ScalarTraits<int>::magnitudeType
    relErr( const int &s1, const int &s2 )
    {
      typedef ScalarTraits<int> ST;
      return ST::magnitude(s1-s2);
    }

  template <>
    ScalarTraits<char>::magnitudeType
    relErr( const char &s1, const char &s2 )
    {
      typedef ScalarTraits<char> ST;
      return ST::magnitude(s1-s2);
    }
}

namespace {

  using Tpetra::TestingUtilities::getNode;
  using Tpetra::TestingUtilities::getDefaultComm;

  using std::endl;
  using std::swap;
  using std::string;

  using Teuchos::TypeTraits::is_same;
  using Teuchos::FancyOStream;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Teuchos::arcpClone;
  using Teuchos::arrayView;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::tuple;
  using Teuchos::null;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;
  using Teuchos::ETransp;
  using Teuchos::NO_TRANS;
  using Teuchos::TRANS;
  using Teuchos::CONJ_TRANS;
  using Teuchos::EDiag;
  using Teuchos::UNIT_DIAG;
  using Teuchos::NON_UNIT_DIAG;
  using Teuchos::EUplo;
  using Teuchos::UPPER_TRI;
  using Teuchos::LOWER_TRI;

  using Tpetra::BlockMap;
  using Tpetra::BlockMultiVector;
  using Tpetra::BlockCrsGraph;
  using Tpetra::Operator;
  using Tpetra::VbrMatrix;
  using Tpetra::Import;
  using Tpetra::global_size_t;
  using Tpetra::DefaultPlatform;
  using Tpetra::ProfileType;
  using Tpetra::StaticProfile;
  using Tpetra::DynamicProfile;
  using Tpetra::OptimizeOption;
  using Tpetra::DoOptimizeStorage;
  using Tpetra::DoNotOptimizeStorage;
  using Tpetra::GloballyDistributed;
  using Tpetra::INSERT;

  double errorTolSlack = 1e+1;
  string filedir;


  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.setOption(
        "filedir",&filedir,"Directory of expected matrix files.");
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &Tpetra::TestingUtilities::testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }

  template<class Scalar,class Ordinal>
  void zero_lower_triangle(Ordinal N, Teuchos::Array<Scalar>& A)
  {
    Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
    for(Ordinal r=0; r<N; ++r) {
      for(Ordinal c=0; c<r; ++c) {
        A[c*N+r] = zero;
      }
    }
  }


  //
  // UNIT TESTS
  //

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( VbrMatrix, Basic, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef BlockMultiVector<Scalar,LO,GO,Node> BMV;
    typedef VbrMatrix<Scalar,LO,GO,Node> MAT;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    GO indexBase = 0;
    // create a Map
    const size_t numLocalBlocks = 2;
    const LO blockSize = 5;
    const size_t maxEntriesPerRow = 3;
    RCP<BlockMap<LO,GO,Node> > rowmap = rcp( new BlockMap<LO,GO,Node>(INVALID,numLocalBlocks,blockSize,indexBase,comm,node) );
    // create the matrix
    {
      RCP<MAT> vbr = rcp( new MAT(rowmap,maxEntriesPerRow,DynamicProfile) );
      TEST_EQUALITY(vbr->hasTransposeApply(), true);
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( VbrMatrix, SetAndGetBlockEntry1, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef BlockMultiVector<Scalar,LO,GO,Node> BMV;
    typedef VbrMatrix<Scalar,LO,GO,Node> MAT;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    GO indexBase = 0;
    // create a Map
    const size_t numLocalBlocks = 2;
    const LO blockSize = 2;
    const size_t maxEntriesPerRow = 2;
    RCP<BlockMap<LO,GO,Node> > rowmap = rcp( new BlockMap<LO,GO,Node>(INVALID,numLocalBlocks,blockSize,indexBase,comm,node) );
    // create the matrix
    {
      RCP<MAT> vbr = rcp( new MAT(rowmap,maxEntriesPerRow,DynamicProfile) );
      Teuchos::ArrayView<const GO> blk_rows = rowmap->getNodeBlockIDs();
      for(int i=0; i<blk_rows.size(); ++i) {
        GO row = blk_rows[i];
        for(int j=0; j<blk_rows.size(); ++j) {
          GO col = blk_rows[j];
          Teuchos::SerialDenseMatrix<int,Scalar> blkEntry(blockSize, blockSize);
          blkEntry.putScalar(row+col+1);
          vbr->setGlobalBlockEntry(row, col, blkEntry);
        }
      }
      for(int i=0; i<blk_rows.size(); ++i) {
        GO row = blk_rows[i];
        for(int j=0; j<blk_rows.size(); ++j) {
          GO col = blk_rows[j];
          LO numPtRows, numPtCols;
          Teuchos::ArrayRCP<const Scalar> blockEntry;
          vbr->getGlobalBlockEntryView(row, col, numPtRows, numPtCols, blockEntry);

          Teuchos::SerialDenseMatrix<int,Scalar> blk(blockSize,blockSize);
          blk.putScalar(row+col+1);

          Teuchos::ArrayRCP<const Scalar> blk_values(blk.values(), 0, blockSize*blockSize, false);
          TEST_COMPARE_FLOATING_ARRAYS( blockEntry, blk_values, 2*Teuchos::ScalarTraits<Scalar>::eps());
        }
      }
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( VbrMatrix, SetAndGetBlockEntry2, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef BlockMultiVector<Scalar,LO,GO,Node> BMV;
    typedef VbrMatrix<Scalar,LO,GO,Node> MAT;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    GO indexBase = 0;
    // create a Map
    const size_t numLocalBlocks = 2;
    const LO blockSize = 2;
    const size_t maxEntriesPerRow = 2;
    RCP<BlockMap<LO,GO,Node> > rowmap = rcp( new BlockMap<LO,GO,Node>(INVALID,numLocalBlocks,blockSize,indexBase,comm,node) );
    // create the matrix
    {
      RCP<MAT> vbr = rcp( new MAT(rowmap,maxEntriesPerRow,DynamicProfile) );
      Teuchos::ArrayView<const GO> blk_rows = rowmap->getNodeBlockIDs();
      for(int i=0; i<blk_rows.size(); ++i) {
        GO row = blk_rows[i];
        for(int j=0; j<blk_rows.size(); ++j) {
          GO col = blk_rows[j];
          Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
          vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
        }
      }
      for(int i=0; i<blk_rows.size(); ++i) {
        GO row = blk_rows[i];
        for(int j=0; j<blk_rows.size(); ++j) {
          GO col = blk_rows[j];
          LO numPtRows, numPtCols;
          Teuchos::ArrayRCP<const Scalar> blockEntry;
          vbr->getGlobalBlockEntryView(row, col, numPtRows, numPtCols, blockEntry);

          Teuchos::SerialDenseMatrix<int,Scalar> blk(blockSize,blockSize);
          blk.putScalar(row+col+1);

          Teuchos::ArrayRCP<const Scalar> blk_values(blk.values(), 0, blockSize*blockSize, false);
          TEST_COMPARE_FLOATING_ARRAYS( blockEntry, blk_values, 2*Teuchos::ScalarTraits<Scalar>::eps());
        }
      }
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( VbrMatrix, SetAndGetBlockEntry3, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef BlockMultiVector<Scalar,LO,GO,Node> BMV;
    typedef VbrMatrix<Scalar,LO,GO,Node> MAT;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    GO indexBase = 0;
    // create a Map
    const size_t numLocalBlocks = 2;
    const LO blockSize = 2;
    const size_t maxEntriesPerRow = 2;
    RCP<BlockMap<LO,GO,Node> > rowmap = rcp( new BlockMap<LO,GO,Node>(INVALID,numLocalBlocks,blockSize,indexBase,comm,node) );
    // create the matrix
    {
      RCP<MAT> vbr = rcp( new MAT(rowmap,maxEntriesPerRow,DynamicProfile) );
      Teuchos::ArrayView<const GO> blk_rows = rowmap->getNodeBlockIDs();
      for(int i=0; i<blk_rows.size(); ++i) {
        GO row = blk_rows[i];
        for(int j=0; j<blk_rows.size(); ++j) {
          GO col = blk_rows[j];
          Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
          vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
          //setLocalBlockEntry should throw since fillComplete hasn't been called yet:
          TEST_THROW(vbr->setLocalBlockEntry(i, j, blockSize, blockSize, blockSize, blkEntry()), std::runtime_error);
        }
      }

      vbr->fillComplete();

      for(int i=0; i<blk_rows.size(); ++i) {
        GO row = blk_rows[i];
        for(int j=0; j<blk_rows.size(); ++j) {
          GO col = blk_rows[j];
          LO numPtRows, numPtCols;
          Teuchos::ArrayRCP<const Scalar> blockEntry;
          vbr->getGlobalBlockEntryView(row, col, numPtRows, numPtCols, blockEntry);

          Teuchos::ArrayRCP<Scalar> nonconstblockEntry;
          vbr->getGlobalBlockEntryViewNonConst(row, col, numPtRows, numPtCols, nonconstblockEntry);

          Teuchos::SerialDenseMatrix<int,Scalar> blk(blockSize,blockSize);
          blk.putScalar(row+col+1);

          Teuchos::ArrayRCP<const Scalar> blk_values(blk.values(), 0, blockSize*blockSize, false);
          TEST_COMPARE_FLOATING_ARRAYS( blockEntry, blk_values, 2*Teuchos::ScalarTraits<Scalar>::eps());
          TEST_COMPARE_FLOATING_ARRAYS( nonconstblockEntry, blk_values, 2*Teuchos::ScalarTraits<Scalar>::eps());
        }
      }
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( VbrMatrix, putScalarAndSumInto1, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef BlockMultiVector<Scalar,LO,GO,Node> BMV;
    typedef VbrMatrix<Scalar,LO,GO,Node> MAT;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    GO indexBase = 0;
    // create a Map
    const size_t numLocalBlocks = 2;
    const LO blockSize = 2;
    const size_t maxEntriesPerRow = 2;
    RCP<BlockMap<LO,GO,Node> > rowmap = rcp( new BlockMap<LO,GO,Node>(INVALID,numLocalBlocks,blockSize,indexBase,comm,node) );
    // create the matrix
    {
      RCP<MAT> vbr = rcp( new MAT(rowmap,maxEntriesPerRow,DynamicProfile) );
      Teuchos::ArrayView<const GO> blk_rows = rowmap->getNodeBlockIDs();
      for(int i=0; i<blk_rows.size(); ++i) {
        GO row = blk_rows[i];
        for(int j=0; j<blk_rows.size(); ++j) {
          GO col = blk_rows[j];
          Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
          vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
        }
      }

      vbr->fillComplete();
      vbr->putScalar(10.0);

      for(int i=0; i<blk_rows.size(); ++i) {
        GO row = blk_rows[i];
        for(int j=0; j<blk_rows.size(); ++j) {
          GO col = blk_rows[j];
          Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
          vbr->sumIntoGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
        }
      }

      for(int i=0; i<blk_rows.size(); ++i) {
        GO row = blk_rows[i];
        for(int j=0; j<blk_rows.size(); ++j) {
          GO col = blk_rows[j];
          Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
          vbr->sumIntoLocalBlockEntry(i, j, blockSize, blockSize, blockSize, blkEntry());
        }
      }

      for(int i=0; i<blk_rows.size(); ++i) {
        GO row = blk_rows[i];
        for(int j=0; j<blk_rows.size(); ++j) {
          GO col = blk_rows[j];
          LO numPtRows, numPtCols;
          Teuchos::ArrayRCP<const Scalar> blockEntry;
          vbr->getGlobalBlockEntryView(row, col, numPtRows, numPtCols, blockEntry);

          Teuchos::ArrayRCP<Scalar> nonconstblockEntry;
          vbr->getGlobalBlockEntryViewNonConst(row, col, numPtRows, numPtCols, nonconstblockEntry);

          Teuchos::SerialDenseMatrix<int,Scalar> blk(blockSize,blockSize);
          blk.putScalar(2*(row+col+1)+10.0);

          Teuchos::ArrayRCP<const Scalar> blk_values(blk.values(), 0, blockSize*blockSize, false);
          TEST_COMPARE_FLOATING_ARRAYS( blockEntry, blk_values, 2*Teuchos::ScalarTraits<Scalar>::eps());
          TEST_COMPARE_FLOATING_ARRAYS( nonconstblockEntry, blk_values, 2*Teuchos::ScalarTraits<Scalar>::eps());
        }
      }
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( VbrMatrix, GlobalLocalRowView1, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef BlockMultiVector<Scalar,LO,GO,Node> BMV;
    typedef VbrMatrix<Scalar,LO,GO,Node> MAT;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    GO indexBase = 0;
    // create a Map
    const size_t numLocalBlocks = 2;
    const LO blockSize = 2;
    const size_t maxEntriesPerRow = 2;
    RCP<BlockMap<LO,GO,Node> > rowmap = rcp( new BlockMap<LO,GO,Node>(INVALID,numLocalBlocks,blockSize,indexBase,comm,node) );

    Teuchos::Array<Teuchos::ArrayRCP<Scalar> > global_row_0;

    // create the matrix
    {
      RCP<MAT> vbr = rcp( new MAT(rowmap,maxEntriesPerRow,DynamicProfile) );
      Teuchos::ArrayView<const GO> blk_rows = rowmap->getNodeBlockIDs();
      for(int i=0; i<blk_rows.size(); ++i) {
        GO row = blk_rows[i];
        for(int j=0; j<blk_rows.size(); ++j) {
          GO col = blk_rows[j];
          if (i==0) {
            global_row_0.push_back(Teuchos::ArrayRCP<Scalar>(blockSize*blockSize, row+col+1));
          }
          Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
          vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
        }
      }

      LO numPtRows;
      Teuchos::ArrayView<const GO> blockCols;
      Teuchos::Array<LO> ptColsPerBlockCol;
      Teuchos::Array<Teuchos::ArrayRCP<const Scalar> > gl_row_0;
      typedef typename Teuchos::Array<Teuchos::ArrayRCP<const Scalar> >::size_type Tsize_t;

      vbr->getGlobalBlockRowView(blk_rows[0], numPtRows, blockCols, ptColsPerBlockCol, gl_row_0);

      TEST_EQUALITY(gl_row_0.size(), global_row_0.size());
      for(Tsize_t i=0; i<gl_row_0.size(); ++i) {
        TEST_EQUALITY(gl_row_0[i].size(), global_row_0[i].size());
        TEST_COMPARE_FLOATING_ARRAYS(gl_row_0[i], global_row_0[i], 2*Teuchos::ScalarTraits<Scalar>::eps());
      }

      vbr->fillComplete();

      Teuchos::ArrayView<const LO> lblockCols;
      Teuchos::Array<LO> lptColsPerBlockCol;
      Teuchos::ArrayRCP<const Scalar> lblockEntries;

      vbr->getLocalBlockRowView(rowmap->getLocalBlockID(blk_rows[0]),
                                numPtRows, lblockCols, lptColsPerBlockCol,
                                lblockEntries);

      TEST_COMPARE_ARRAYS(ptColsPerBlockCol, lptColsPerBlockCol);
      LO offset = 0;
      for(Tsize_t i=0; i<gl_row_0.size(); ++i) {
        LO len = gl_row_0[i].size();
        TEST_COMPARE_FLOATING_ARRAYS(gl_row_0[i], Teuchos::ArrayView<const Scalar>(lblockEntries.getRawPtr()+offset, len), 2*Teuchos::ScalarTraits<Scalar>::eps());
        offset += len;
      }

      vbr->putScalar(10.0);

      for(int i=0; i<blk_rows.size(); ++i) {
        GO row = blk_rows[i];
        for(int j=0; j<blk_rows.size(); ++j) {
          GO col = blk_rows[j];
          Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
          vbr->sumIntoGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
        }
      }

      for(int i=0; i<blk_rows.size(); ++i) {
        GO row = blk_rows[i];
        for(int j=0; j<blk_rows.size(); ++j) {
          GO col = blk_rows[j];
          Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
          vbr->sumIntoLocalBlockEntry(i, j, blockSize, blockSize, blockSize, blkEntry());
        }
      }

      for(int i=0; i<blk_rows.size(); ++i) {
        GO row = blk_rows[i];
        for(int j=0; j<blk_rows.size(); ++j) {
          GO col = blk_rows[j];
          LO numPtRows, numPtCols;
          Teuchos::ArrayRCP<const Scalar> blockEntry;
          vbr->getGlobalBlockEntryView(row, col, numPtRows, numPtCols, blockEntry);

          Teuchos::ArrayRCP<Scalar> nonconstblockEntry;
          vbr->getGlobalBlockEntryViewNonConst(row, col, numPtRows, numPtCols, nonconstblockEntry);

          Teuchos::SerialDenseMatrix<int,Scalar> blk(blockSize,blockSize);
          blk.putScalar(2*(row+col+1)+10.0);

          Teuchos::ArrayRCP<const Scalar> blk_values(blk.values(), 0, blockSize*blockSize, false);
          TEST_COMPARE_FLOATING_ARRAYS( blockEntry, blk_values, 2*Teuchos::ScalarTraits<Scalar>::eps());
          TEST_COMPARE_FLOATING_ARRAYS( nonconstblockEntry, blk_values, 2*Teuchos::ScalarTraits<Scalar>::eps());
        }
      }
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( VbrMatrix, applySimple, LO, GO, Scalar, Node )
  {
    //This test builds a block-diagonal matrix and tests the apply method.
    //Since the matrix is block-diagonal, the apply does not require any
    //communication. (No import of the X vector is required.)
    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef BlockMultiVector<Scalar,LO,GO,Node> BMV;
    typedef VbrMatrix<Scalar,LO,GO,Node> MAT;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    GO indexBase = 0;
    // create a Map
    const size_t numLocalBlocks = 2;
    const LO blockSize = 2;
    const size_t maxEntriesPerRow = 2;
    RCP<BlockMap<LO,GO,Node> > rowmap = rcp( new BlockMap<LO,GO,Node>(INVALID,numLocalBlocks,blockSize,indexBase,comm,node) );
    // create the matrix
    {
      RCP<BMV> bmv1 = rcp(new BMV(rowmap, 1));
      RCP<BMV> bmv2 = rcp(new BMV(rowmap, 1));
      RCP<BMV> bmv3 = rcp(new BMV(rowmap, 1));
      bmv1->putScalar(1.0);
      ArrayRCP<Scalar> v3 = bmv3->get1dViewNonConst();
      RCP<MAT> vbr = rcp( new MAT(rowmap,maxEntriesPerRow,DynamicProfile) );
      Teuchos::ArrayView<const GO> blk_rows = rowmap->getNodeBlockIDs();
      LO row_offset = 0;
      for(int i=0; i<blk_rows.size(); ++i) {
        Scalar val = 0;
        GO row = blk_rows[i];
        for(int j=0; j<blk_rows.size(); ++j) {
          GO col = blk_rows[j];
          Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
          val += (row+col+1)*blockSize;
          vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
        }
        for(int k=0; k<blockSize; ++k) v3[row_offset++] = val;
      }

      vbr->fillComplete();
      vbr->apply(*bmv1, *bmv2);
      ArrayRCP<Scalar> v2 = bmv2->get1dViewNonConst();
      TEST_COMPARE_FLOATING_ARRAYS( v2, v3, 2*Teuchos::ScalarTraits<Scalar>::eps());
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( VbrMatrix, solveUpperNonUnit, LO, GO, Scalar, Node )
  {
    //This test builds a block-diagonal, upper-triangular matrix and tests the
    //solve method.
    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef BlockMultiVector<Scalar,LO,GO,Node> BMV;
    typedef VbrMatrix<Scalar,LO,GO,Node> MAT;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    GO indexBase = 0;
    // create a Map
    const size_t numLocalBlocks = 3;
    const LO blockSize = 3;
    const size_t maxEntriesPerRow = 3;
    RCP<BlockMap<LO,GO,Node> > rowmap = rcp( new BlockMap<LO,GO,Node>(INVALID,numLocalBlocks,blockSize,indexBase,comm,node) );
    // create the matrix
    {
      RCP<MAT> vbr = rcp( new MAT(rowmap,maxEntriesPerRow,DynamicProfile) );
      Teuchos::ArrayView<const GO> blk_rows = rowmap->getNodeBlockIDs();
      for(int i=0; i<blk_rows.size(); ++i) {
        Scalar val = 0;
        GO row = blk_rows[i];
        for(int j=i; j<blk_rows.size(); ++j) {
          GO col = blk_rows[j];
          Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
          if (row==col) zero_lower_triangle(blockSize, blkEntry);
          val += (row+col+1)*blockSize;
          vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
        }
      }

      vbr->fillComplete();

      RCP<BMV> x = rcp(new BMV(rowmap, 1));
      x->putScalar(1.0);
      RCP<BMV> y = rcp(new BMV(rowmap, 1));
      RCP<BMV> x2 = rcp(new BMV(rowmap, 1));

      vbr->apply(*x, *y);
      vbr->applyInverse(*y, *x2, Teuchos::NO_TRANS);
      ArrayRCP<const Scalar> v_x = x->get1dView();
      ArrayRCP<const Scalar> v_x2 = x2->get1dView();
      TEST_COMPARE_FLOATING_ARRAYS( v_x, v_x2, 2*Teuchos::ScalarTraits<Scalar>::eps());
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( VbrMatrix, solveTransposeUpperNonUnit, LO, GO, Scalar, Node )
  {
    //This test builds a block-diagonal, upper-triangular matrix and tests the
    //transpose-solve method.
    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef BlockMultiVector<Scalar,LO,GO,Node> BMV;
    typedef VbrMatrix<Scalar,LO,GO,Node> MAT;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    GO indexBase = 0;
    // create a Map
    const size_t numLocalBlocks = 3;
    const LO blockSize = 3;
    const size_t maxEntriesPerRow = 3;
    RCP<BlockMap<LO,GO,Node> > rowmap = rcp( new BlockMap<LO,GO,Node>(INVALID,numLocalBlocks,blockSize,indexBase,comm,node) );
    // create the matrix
    {
      RCP<MAT> vbr = rcp( new MAT(rowmap,maxEntriesPerRow,DynamicProfile) );
      Teuchos::ArrayView<const GO> blk_rows = rowmap->getNodeBlockIDs();
      for(int i=0; i<blk_rows.size(); ++i) {
        Scalar val = 0;
        GO row = blk_rows[i];
        for(int j=i; j<blk_rows.size(); ++j) {
          GO col = blk_rows[j];
          Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
          if (row==col) zero_lower_triangle(blockSize, blkEntry);
          val += (row+col+1)*blockSize;
          vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
        }
      }

      vbr->fillComplete();

      RCP<BMV> x = rcp(new BMV(rowmap, 1));
      x->putScalar(1.0);
      RCP<BMV> y = rcp(new BMV(rowmap, 1));
      RCP<BMV> x2 = rcp(new BMV(rowmap, 1));

      vbr->apply(*x, *y, Teuchos::TRANS);
      vbr->applyInverse(*y, *x2, Teuchos::TRANS);
      ArrayRCP<const Scalar> v_x = x->get1dView();
      ArrayRCP<const Scalar> v_x2 = x2->get1dView();
      TEST_COMPARE_FLOATING_ARRAYS( v_x, v_x2, 2*Teuchos::ScalarTraits<Scalar>::eps());
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( VbrMatrix, ColMap1, LO, GO, Scalar, Node )
  {
    //This test fills a (block-tri-diagonal) matrix such that in parallel
    //the column-map should have an overlapping set of entries (i.e.,
    //different than the row-map), and verify that the column-map is correct.

    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef BlockMultiVector<Scalar,LO,GO,Node> BMV;
    typedef VbrMatrix<Scalar,LO,GO,Node> MAT;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
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

    // create the matrix
    RCP<MAT> vbr = rcp( new MAT(rowmap,maxEntriesPerRow,DynamicProfile) );
    for(int i=0; i<blk_rows.size(); ++i) {
      GO row = blk_rows[i];
      {
        GO col = row;
        Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
        vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
      }
      if (row > indexBase) {
        GO col = row - 1;
        Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
        vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
      }
      if (row < maxGlobalBlock-1) {
        GO col = row + 1;
        Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
        vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
      }
    }

    vbr->fillComplete();
    RCP<const BlockMap<LO,GO,Node> > colmap = vbr->getBlockColMap();
    ArrayView<const GO> blk_cols = colmap->getNodeBlockIDs();
    TEST_EQUALITY(blk_cols.size(), blockColIDs.size());
    TEST_COMPARE_ARRAYS(blk_cols, blockColIDs() );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( VbrMatrix, applyParallel, LO, GO, Scalar, Node )
  {
    //This test fills a (block-tri-diagonal) matrix such that in parallel the
    //column-map should have an overlapping set of entries (i.e., different than
    //the row-map), and verify that apply works correctly. If the column-map
    //has an overlapping set of entries, then apply must do an import of the x
    //vector in order to get a correct result.

    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef BlockMultiVector<Scalar,LO,GO,Node> BMV;
    typedef VbrMatrix<Scalar,LO,GO,Node> MAT;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
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
    LO row_offset = 0;
    for(LO i=0; i<blk_rows.size(); ++i) {
      blockColIDs[offset++] = blk_rows[i];
      last_row = blk_rows[i];
    }
    if (offset < blockColIDs.size()) blockColIDs[offset++] = last_row + 1;

    RCP<BMV> bmv1 = rcp(new BMV(rowmap, 1));
    RCP<BMV> bmv2 = rcp(new BMV(rowmap, 1));
    RCP<BMV> bmv3 = rcp(new BMV(rowmap, 1));
    bmv1->putScalar(1.0);
    ArrayRCP<Scalar> v3 = bmv3->get1dViewNonConst();

    // create the matrix
    RCP<MAT> vbr = rcp( new MAT(rowmap,maxEntriesPerRow,DynamicProfile) );
    for(int i=0; i<blk_rows.size(); ++i) {
      Scalar val = 0;
      GO row = blk_rows[i];
      {
        GO col = row;
        Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
        vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
        val += (row+col+1)*blockSize;
      }
      if (row > indexBase) {
        GO col = row - 1;
        Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
        vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
        val += (row+col+1)*blockSize;
      }
      if (row < maxGlobalBlock-1) {
        GO col = row + 1;
        Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
        vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
        val += (row+col+1)*blockSize;
      }
      for(int k=0; k<blockSize; ++k) v3[row_offset++] = val;
    }

    vbr->fillComplete();

    vbr->apply(*bmv1, *bmv2);

    ArrayRCP<Scalar> v2 = bmv2->get1dViewNonConst();
    TEST_COMPARE_FLOATING_ARRAYS( v2, v3, 2*Teuchos::ScalarTraits<Scalar>::eps());
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( VbrMatrix, applyTransParallel, LO, GO, Scalar, Node )
  {
    //This test fills a (block-tri-diagonal) matrix such that in parallel the
    //column-map should have an overlapping set of entries (i.e., different than
    //the row-map), and verify that apply transpose works correctly. If the column-map
    //has an overlapping set of entries, then apply must do an export of the y
    //vector in order to get a correct result.

    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef BlockMultiVector<Scalar,LO,GO,Node> BMV;
    typedef VbrMatrix<Scalar,LO,GO,Node> MAT;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    GO indexBase = 0;
    // create a Map
    const size_t numLocalBlocks = 2;
    GO numGlobalBlocks = numLocalBlocks*comm->getSize();
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
    LO row_offset = 0;
    for(LO i=0; i<blk_rows.size(); ++i) {
      blockColIDs[offset++] = blk_rows[i];
      last_row = blk_rows[i];
    }
    if (offset < blockColIDs.size()) blockColIDs[offset++] = last_row + 1;

    RCP<BMV> bmv1 = rcp(new BMV(rowmap, 1));
    RCP<BMV> bmv2 = rcp(new BMV(rowmap, 1));
    RCP<BMV> bmv3 = rcp(new BMV(rowmap, 1));
    bmv1->putScalar(1.0);
    ArrayRCP<Scalar> v3 = bmv3->get1dViewNonConst();

    // create the matrix
    RCP<MAT> vbr = rcp( new MAT(rowmap,maxEntriesPerRow,DynamicProfile) );
    for(int i=0; i<blk_rows.size(); ++i) {
      Scalar val = 0;
      GO row = blk_rows[i];
      {
        GO col = row;
        Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
        vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
        val += (row+col+1)*blockSize;
      }
      if (row > indexBase) {
        GO col = row - 1;
        Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
        vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
        val += (row+col+1)*blockSize;
      }
      if (row < numGlobalBlocks-1) {
        GO col = row + 1;
        Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
        vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
        val += (row+col+1)*blockSize;
      }
      for(int k=0; k<blockSize; ++k) v3[row_offset++] = val;
    }

    vbr->fillComplete();

    vbr->apply(*bmv1, *bmv2, Teuchos::TRANS);

    ArrayRCP<Scalar> v2 = bmv2->get1dViewNonConst();
    TEST_COMPARE_FLOATING_ARRAYS( v2, v3, 2*Teuchos::ScalarTraits<Scalar>::eps());
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( VbrMatrix, applyParallelMV, LO, GO, Scalar, Node )
  {
    //This test fills a (block-tri-diagonal) matrix such that in parallel the
    //column-map should have an overlapping set of entries (i.e., different than
    //the row-map), and verify that apply works correctly. If the column-map
    //has an overlapping set of entries, then apply must do an import of the x
    //vector in order to get a correct result.
    //This test uses MultiVectors that have more than 1 column.

    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef BlockMultiVector<Scalar,LO,GO,Node> BMV;
    typedef VbrMatrix<Scalar,LO,GO,Node> MAT;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
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
    LO row_offset = 0;
    for(LO i=0; i<blk_rows.size(); ++i) {
      blockColIDs[offset++] = blk_rows[i];
      last_row = blk_rows[i];
    }
    if (offset < blockColIDs.size()) blockColIDs[offset++] = last_row + 1;

    size_t numVecs = 3;
    size_t numPointRows = rowmap->getPointMap()->getNodeNumElements();
    RCP<BMV> bmv1 = rcp(new BMV(rowmap, numVecs));
    RCP<BMV> bmv2 = rcp(new BMV(rowmap, numVecs));
    RCP<BMV> bmv3 = rcp(new BMV(rowmap, numVecs));
    bmv1->putScalar(1.0);
    ArrayRCP<Scalar> v3 = bmv3->get1dViewNonConst();

    // create the matrix
    RCP<MAT> vbr = rcp( new MAT(rowmap,maxEntriesPerRow,DynamicProfile) );
    for(int i=0; i<blk_rows.size(); ++i) {
      Scalar val = 0;
      GO row = blk_rows[i];
      {
        GO col = row;
        Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
        vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
        val += (row+col+1)*blockSize;
      }
      if (row > indexBase) {
        GO col = row - 1;
        Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
        vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
        val += (row+col+1)*blockSize;
      }
      if (row < maxGlobalBlock-1) {
        GO col = row + 1;
        Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
        vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
        val += (row+col+1)*blockSize;
      }
      for(int k=0; k<blockSize; ++k) {
        for(size_t v=0; v<numVecs; ++v) {
          v3[v*numPointRows+row_offset] = val;
        }
        ++row_offset;
      }
    }

    vbr->fillComplete();

    vbr->apply(*bmv1, *bmv2);

    ArrayRCP<Scalar> v2 = bmv2->get1dViewNonConst();
    TEST_COMPARE_FLOATING_ARRAYS( v2, v3, 2*Teuchos::ScalarTraits<Scalar>::eps());
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( VbrMatrix, Import1, LO, GO, Scalar, Node )
  {
    //This test fills a (block-diagonal) matrix and then imports it to a
    //different matrix with a different row-map.

    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef VbrMatrix<Scalar,LO,GO,Node> MAT;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();

    //For now, this is just a 2-proc test:
    if (comm->getSize() != 2) return;

    GO indexBase = 0;
    // create a Map
    const size_t numLocalBlocks_src = comm->getRank()==0 ? 2 : 1;
    const size_t numLocalBlocks_dest= comm->getRank()==0 ? 1 : 2;

    const LO blockSize = 2;
    const size_t maxEntriesPerRow = 1;
    RCP<const BlockMap<LO,GO,Node> > rowmap_src = rcp( new BlockMap<LO,GO,Node>(INVALID,numLocalBlocks_src,blockSize,indexBase,comm,node) );
    RCP<const BlockMap<LO,GO,Node> > rowmap_dest = rcp( new BlockMap<LO,GO,Node>(INVALID,numLocalBlocks_dest,blockSize,indexBase,comm,node) );

    typedef typename Array<GO>::size_type Tsize_t;
    Teuchos::ArrayView<const GO> blk_rows = rowmap_src->getNodeBlockIDs();

    // create the matrix
    RCP<MAT> vbr = rcp( new MAT(rowmap_src,maxEntriesPerRow,DynamicProfile) );
    for(int i=0; i<blk_rows.size(); ++i) {
      GO row = blk_rows[i];
      GO col = row;
      Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
      vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
    }

//    vbr->fillComplete();

    Tpetra::Import<LO,GO,Node> importer(convertBlockMapToPointMap(*rowmap_src),
                                        convertBlockMapToPointMap(*rowmap_dest));

    //Construct the matrix that will be imported into:
    RCP<MAT> vbr_dest = rcp( new MAT(rowmap_dest, maxEntriesPerRow,DynamicProfile) );

    vbr_dest->doImport(*vbr, importer, Tpetra::REPLACE);

    Teuchos::ArrayView<const GO> blk_rows_dest = rowmap_dest->getNodeBlockIDs();
    for(int i=0; i<blk_rows_dest.size(); ++i) {
      GO row = blk_rows_dest[i];
      GO col = row;
      Teuchos::Array<Scalar> expected_blockEntry(blockSize*blockSize, row+col+1);
      Teuchos::ArrayRCP<const Scalar> blockEntry;
      LO tmp;
      vbr_dest->getGlobalBlockEntryView(row, col, tmp, tmp, blockEntry);
      TEST_COMPARE_FLOATING_ARRAYS( expected_blockEntry, blockEntry, 2*Teuchos::ScalarTraits<Scalar>::eps());
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( VbrMatrix, Overlap0, LO, GO, Scalar, Node )
  {
    //This test creates a distributed matrix with 1 block-row per proc.
    //Block-equations each have 2 point-equations.
    //Each block-row will have just 1 block-column, with block-index==0.
    //Each proc will put a block-entry into the matrix for its own
    //block-row, and also one for each neighboring proc (proc-1 and proc+1).
    //This will create a matrix with overlapping data, which gets
    //assembled to the owning processors by globalAssemble which is called
    //by fillComplete.

    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef VbrMatrix<Scalar,LO,GO,Node> MAT;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    GO indexBase = 0;
    // create a Map
    const size_t numLocalBlocks = 1;
    const LO blockSize = 2;
    RCP<const BlockMap<LO,GO,Node> > rowmap = rcp( new BlockMap<LO,GO,Node>(INVALID,numLocalBlocks,blockSize,indexBase,comm,node) );

    typedef typename Array<GO>::size_type Tsize_t;
    const size_t maxEntriesPerRow = 1;

    // create the matrix
    RCP<MAT> vbr = rcp( new MAT(rowmap,maxEntriesPerRow,DynamicProfile) );

    // fill the matrix
    GO row = comm->getRank();
    GO col = 0;
    Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, 1);
    vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());

    if (comm->getRank() > 0) {
      GO row = comm->getRank()-1;
      Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, 1);
      vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
    }

    if (comm->getRank() < comm->getSize()-1) {
      GO row = comm->getRank()+1;
      Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, 1);
      vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
    }

    vbr->fillComplete();

    double num_contributions = 1;//self
    if (comm->getRank() > 0) ++num_contributions;//proc-1 contributed to my row
    if (comm->getRank() < comm->getSize()-1) ++num_contributions;//proc+1 contributed to my row
    Teuchos::ArrayView<const GO> blk_rows = rowmap->getNodeBlockIDs();
    for(Tsize_t i=0; i<blk_rows.size(); ++i) {
      GO row = blk_rows[i];
      GO col = row;
      Teuchos::Array<Scalar> expected_blockEntry(blockSize*blockSize, num_contributions);
      Teuchos::ArrayRCP<const Scalar> blockEntry;
      LO tmp;
      vbr->getGlobalBlockEntryView(row, col, tmp, tmp, blockEntry);
      TEST_COMPARE_FLOATING_ARRAYS( expected_blockEntry, blockEntry, 2*Teuchos::ScalarTraits<Scalar>::eps());
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( VbrMatrix, Overlap1, LO, GO, Scalar, Node )
  {
    //This test creates a (block-diagonal) matrix that is distributed,
    //then fills it by loading all data on proc 0, then verifies that after
    //fillComplete (which calls globalAssemble) the data is spread out to
    //the other procs.

    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef VbrMatrix<Scalar,LO,GO,Node> MAT;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();

    //For now, this is just a 2-proc test:
    if (comm->getSize() != 2) return;

    GO indexBase = 0;
    // create a Map
    const size_t numLocalBlocks = 2;

    const LO blockSize = 2;
    const size_t maxEntriesPerRow = 1;
    RCP<const BlockMap<LO,GO,Node> > rowmap = rcp( new BlockMap<LO,GO,Node>(INVALID,numLocalBlocks,blockSize,indexBase,comm,node) );

    typedef typename Array<GO>::size_type Tsize_t;

    // create the matrix
    RCP<MAT> vbr = rcp( new MAT(rowmap,maxEntriesPerRow,DynamicProfile) );

    // fill the matrix
    if (comm->getRank() == 0) {
      for(int i=0; i<int(comm->getSize()*numLocalBlocks); ++i) {
        GO row = i;
        GO col = row;
        Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
        vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
      }
    }

    vbr->fillComplete();

    Teuchos::ArrayView<const GO> blk_rows = rowmap->getNodeBlockIDs();
    for(Tsize_t i=0; i<blk_rows.size(); ++i) {
      GO row = blk_rows[i];
      GO col = row;
      Teuchos::Array<Scalar> expected_blockEntry(blockSize*blockSize, row+col+1);
      Teuchos::ArrayRCP<const Scalar> blockEntry;
      LO tmp;
      vbr->getGlobalBlockEntryView(row, col, tmp, tmp, blockEntry);
      TEST_COMPARE_FLOATING_ARRAYS( expected_blockEntry, blockEntry, 2*Teuchos::ScalarTraits<Scalar>::eps());
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( VbrMatrix, CtorBCrsGraph, LO, GO, Scalar, Node )
  {
    //This test creates a (block-tri-diagonal) matrix using a pre-filled
    //BlockCrsGraph such that in parallel the column-map should have an
    //overlapping set of entries (i.e., different than the row-map),
    //and verify that apply works correctly. If the column-map
    //has an overlapping set of entries, then apply must do an import of the x
    //vector in order to get a correct result.

    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef BlockCrsGraph<LO,GO,Node> BGRAPH;
    typedef BlockMultiVector<Scalar,LO,GO,Node> BMV;
    typedef VbrMatrix<Scalar,LO,GO,Node> MAT;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    GO indexBase = 0;
    // create a Map
    const size_t numLocalBlocks = 2;
    GO maxGlobalBlock = numLocalBlocks*comm->getSize();
    const LO blockSize = 2;
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

    const size_t maxEntriesPerRow = 3;
    // create and fill the graph
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

    RCP<BMV> bmv3 = rcp(new BMV(rowmap, 1));
    ArrayRCP<Scalar> v3 = bmv3->get1dViewNonConst();

    // create and fill the matrix
    RCP<MAT> vbr = rcp( new MAT(bgrph) );
    LO row_offset = 0;
    for(int i=0; i<blk_rows.size(); ++i) {
      Scalar val = 0;
      GO row = blk_rows[i];
      {
        GO col = row;
        Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
        vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
        val += (row+col+1)*blockSize;
      }
      if (row > indexBase) {
        GO col = row - 1;
        Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
        vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
        val += (row+col+1)*blockSize;
      }
      if (row < maxGlobalBlock-1) {
        GO col = row + 1;
        Teuchos::Array<Scalar> blkEntry(blockSize * blockSize, row+col+1);
        vbr->setGlobalBlockEntry(row, col, blockSize, blockSize, blockSize, blkEntry());
        val += (row+col+1)*blockSize;
      }
      for(int k=0; k<blockSize; ++k) v3[row_offset++] = val;
    }

    vbr->fillComplete();

    RCP<BMV> bmv1 = rcp(new BMV(rowmap, 1));
    RCP<BMV> bmv2 = rcp(new BMV(rowmap, 1));
    bmv1->putScalar(1.0);

    vbr->apply(*bmv1, *bmv2);

    ArrayRCP<Scalar> v2 = bmv2->get1dViewNonConst();
    TEST_COMPARE_FLOATING_ARRAYS( v2, v3, 2*Teuchos::ScalarTraits<Scalar>::eps());
  }

//
// INSTANTIATIONS
//


#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( VbrMatrix, Basic, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( VbrMatrix, SetAndGetBlockEntry1, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( VbrMatrix, SetAndGetBlockEntry2, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( VbrMatrix, SetAndGetBlockEntry3, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( VbrMatrix, putScalarAndSumInto1, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( VbrMatrix, GlobalLocalRowView1, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( VbrMatrix, applySimple, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( VbrMatrix, solveUpperNonUnit, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( VbrMatrix, solveTransposeUpperNonUnit, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( VbrMatrix, ColMap1, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( VbrMatrix, applyParallel, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( VbrMatrix, applyTransParallel, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( VbrMatrix, applyParallelMV, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( VbrMatrix, CtorBCrsGraph, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( VbrMatrix, Import1, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( VbrMatrix, Overlap0, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( VbrMatrix, Overlap1, LO, GO, SCALAR, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV_NOGPU( UNIT_TEST_GROUP )

}
