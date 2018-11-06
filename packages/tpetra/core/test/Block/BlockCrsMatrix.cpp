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

#include "Tpetra_TestingUtilities.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_BlockCrsMatrix.hpp"
#include "Tpetra_BlockCrsMatrix_Helpers.hpp"
#include "Tpetra_BlockVector.hpp"

namespace {
  using Tpetra::TestingUtilities::getDefaultComm;
  using Tpetra::Details::gathervPrint;
  using Teuchos::Array;
  using Teuchos::Comm;
  using Teuchos::outArg;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::endl;
  typedef Tpetra::global_size_t GST;

  //
  // UNIT TESTS
  //

  // Test BlockCrsMatrix's constructors.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockCrsMatrix, ctor, Scalar, LO, GO, Node )
  {
    typedef Tpetra::BlockMultiVector<Scalar, LO, GO, Node> BMV;
    typedef Tpetra::BlockCrsMatrix<Scalar, LO, GO, Node> BCM;
    typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;

    out << "Testing Tpetra::BlockCrsMatrix" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const Comm<int> > comm = getDefaultComm ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    out << "Creating mesh row Map" << endl;

    const size_t numLocalMeshPoints = 12;
    const GO indexBase = 1;
    // mfh 16 May 2014: Tpetra::CrsGraph still needs the row Map as an
    // RCP.  Later interface changes will let us pass in the Map by
    // const reference and assume view semantics.
    RCP<const map_type> meshRowMapPtr =
      rcp (new map_type (INVALID, numLocalMeshPoints, indexBase, comm));

    const LO blockSize = 4;
    //const LO numVecs = 3;

    // mfh 16 May 2014: Make a graph.  This graph is empty; that's OK.
    // We just want to check that the constructors of BlockCrsMatrix
    // work.
    out << "Creating mesh graph" << endl;
    graph_type graph (meshRowMapPtr, 0, Tpetra::StaticProfile);
    graph.fillComplete ();

    // Test the default constructor.  We can't wrap this in a
    // TEST_NOTHROW, because the variable need to be in scope for
    // tests below.
    out << "Testing default constructor" << endl;
    BCM blockMat;

    // Test the two-argument constructor.
    out << "Testing two-argument constructor" << endl;
    TEST_NOTHROW( blockMat = BCM (graph, blockSize) );

    // Test that the point domain and range Maps are correct.
    map_type pointDomainMap = BMV::makePointMap (* (graph.getDomainMap ()), blockSize);
    TEST_ASSERT( ! blockMat.getDomainMap ().is_null () &&
                 pointDomainMap.isSameAs (* (blockMat.getDomainMap ())) );
    map_type pointRangeMap = BMV::makePointMap (* (graph.getRangeMap ()), blockSize);
    TEST_ASSERT( ! blockMat.getRangeMap ().is_null () &&
                 pointRangeMap.isSameAs (* (blockMat.getRangeMap ())) );

    // Test that the result of getCrsGraph() has the same Maps.
    {
      graph_type graph2 = blockMat.getCrsGraph ();
      TEST_ASSERT( ! graph.getDomainMap ().is_null () &&
                   ! graph2.getDomainMap ().is_null () &&
                   graph.getDomainMap ()->isSameAs (* (graph2.getDomainMap ())) );
      TEST_ASSERT( ! graph.getRangeMap ().is_null () &&
                   ! graph2.getRangeMap ().is_null () &&
                   graph.getRangeMap ()->isSameAs (* (graph2.getRangeMap ())) );
      TEST_ASSERT( ! graph.getRowMap ().is_null () &&
                   ! graph2.getRowMap ().is_null () &&
                   graph.getRowMap ()->isSameAs (* (graph2.getRowMap ())) );
      TEST_ASSERT( ! graph.getColMap ().is_null () &&
                   ! graph2.getColMap ().is_null () &&
                   graph.getColMap ()->isSameAs (* (graph2.getColMap ())) );
    }

    // Test the four-argument constructor.
    out << "Testing four-argument constructor" << endl;
    TEST_NOTHROW( blockMat = BCM (graph, pointDomainMap, pointRangeMap, blockSize ) );

    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
  }

  // Test for basic functionality of a BlockCrsMatrix.  We create a
  // BlockCrsMatrix with a nontrivial graph; exercise getLocalRowView,
  // getLocalRowCopy, and replaceLocalValues; and test applyBlock and apply.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockCrsMatrix, basic, Scalar, LO, GO, Node )
  {
    typedef Tpetra::BlockMultiVector<Scalar, LO, GO, Node> BMV;
    typedef Tpetra::BlockVector<Scalar, LO, GO, Node> BV;
    typedef Tpetra::BlockCrsMatrix<Scalar, LO, GO, Node> BCM;
    typedef Tpetra::MultiVector<Scalar, LO, GO, Node> mv_type;
    typedef Tpetra::Vector<Scalar, LO, GO, Node> vec_type;
    typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    // The typedef below is also a test.  BlockCrsMatrix must have
    // this typedef, or this test won't compile.
    typedef typename BCM::little_block_type little_block_type;
    typedef typename BV::little_vec_type little_vec_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType MT;

    const Scalar two = STS::one () + STS::one ();
    const Scalar three = STS::one () + STS::one () + STS::one ();

    out << "Testing Tpetra::BlockCrsMatrix basic "
      "functionality" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const Comm<int> > comm = getDefaultComm ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    out << "Creating mesh row Map" << endl;

    const size_t numLocalMeshPoints = 12;
    const GO indexBase = 1;
    // Use a block size that is not a power of 2, to test correctness
    // in case the matrix pads blocks for SIMD-ization.
    const LO blockSize = 3;
    const size_t entriesPerBlock = blockSize * blockSize;

    // mfh 20 May 2014: Tpetra::CrsGraph still needs the row Map as an
    // RCP.  Later interface changes will let us pass in the Map by
    // const reference and assume view semantics.
    RCP<const map_type> meshRowMapPtr =
      rcp (new map_type (INVALID, numLocalMeshPoints, indexBase, comm));
    const map_type& meshRowMap = *meshRowMapPtr;

    // Make a graph.  It will have two entries per global row i_gbl:
    // (i_gbl, (i_gbl+1) % N) and (i_gbl, (i_gbl+2) % N), where N is
    // the global number of rows and columns.  We don't include the
    // diagonal, to make the Maps more interesting.
    out << "Creating mesh graph" << endl;

    const size_t maxNumEntPerRow = 2;
    graph_type graph (meshRowMapPtr, maxNumEntPerRow, Tpetra::StaticProfile);

    // Fill the graph.
    Teuchos::Array<GO> gblColInds (maxNumEntPerRow);
    const GO globalNumRows = meshRowMap.getGlobalNumElements ();
    for (LO lclRowInd = meshRowMap.getMinLocalIndex ();
         lclRowInd <= meshRowMap.getMaxLocalIndex (); ++lclRowInd) {
      const GO gblRowInd = meshRowMap.getGlobalElement (lclRowInd);
      for (size_t k = 0; k < maxNumEntPerRow; ++k) {
        const GO gblColInd = indexBase +
          ((gblRowInd - indexBase) + static_cast<GO> (k + 1)) %
          static_cast<GO> (globalNumRows);
        gblColInds[k] = gblColInd;
      }
      graph.insertGlobalIndices (gblRowInd, gblColInds ());
    }
    graph.fillComplete ();

    // Get the graph's column Map (the "mesh column Map").
    TEST_ASSERT( ! graph.getColMap ().is_null () );
    map_type meshColMap = * (graph.getColMap ());

    out << "Creating BlockCrsMatrix" << endl;

    // Construct the BlockCrsMatrix.
    BCM blockMat;
    TEST_NOTHROW( blockMat = BCM (graph, blockSize) );

    // Test that the matrix's block size is correct.
    TEST_ASSERT( blockMat.getBlockSize () == blockSize );

    // Test that the point domain and range Maps are correct.
    map_type pointDomainMap =
      BMV::makePointMap (* (graph.getDomainMap ()), blockSize);
    TEST_ASSERT( ! blockMat.getDomainMap ().is_null () &&
                 pointDomainMap.isSameAs (* (blockMat.getDomainMap ())) );
    map_type pointRangeMap =
      BMV::makePointMap (* (graph.getRangeMap ()), blockSize);
    TEST_ASSERT( ! blockMat.getRangeMap ().is_null () &&
                 pointRangeMap.isSameAs (* (blockMat.getRangeMap ())) );

    // Test that the result of getCrsGraph() has the same Maps.
    {
      graph_type graph2 = blockMat.getCrsGraph ();
      TEST_ASSERT( ! graph.getDomainMap ().is_null () &&
                   ! graph2.getDomainMap ().is_null () &&
                   graph.getDomainMap ()->isSameAs (* (graph2.getDomainMap ())) );
      TEST_ASSERT( ! graph.getRangeMap ().is_null () &&
                   ! graph2.getRangeMap ().is_null () &&
                   graph.getRangeMap ()->isSameAs (* (graph2.getRangeMap ())) );
      TEST_ASSERT( ! graph.getRowMap ().is_null () &&
                   ! graph2.getRowMap ().is_null () &&
                   graph.getRowMap ()->isSameAs (* (graph2.getRowMap ())) );
      TEST_ASSERT( ! graph.getColMap ().is_null () &&
                   ! graph2.getColMap ().is_null () &&
                   graph.getColMap ()->isSameAs (* (graph2.getColMap ())) );
    }

    out << "Test getLocalRowView, getLocalRowCopy, and replaceLocalValues" << endl;

    Array<Scalar> tempBlockSpace (maxNumEntPerRow * entriesPerBlock);

    // Test that getLocalRowView returns the right column indices.
    Array<LO> lclColInds (maxNumEntPerRow);
    Array<LO> myLclColIndsCopy (maxNumEntPerRow);
    Array<Scalar> myValsCopy (maxNumEntPerRow*entriesPerBlock);
    Array<LO> myLclColIndsSorted (maxNumEntPerRow);
    for (LO lclRowInd = meshRowMap.getMinLocalIndex ();
         lclRowInd <= meshRowMap.getMaxLocalIndex (); ++lclRowInd) {
      const LO* myLclColInds = NULL;
      Scalar* myVals = NULL;
      LO numEnt = 0;
      LO err = blockMat.getLocalRowView (lclRowInd, myLclColInds, myVals, numEnt);
      TEST_ASSERT( err == 0 );
      TEST_ASSERT( numEnt == static_cast<LO> (maxNumEntPerRow) );
      TEST_ASSERT( myLclColInds != NULL );
      TEST_ASSERT( myVals != NULL );

      // Compute what the local column indices in this row _should_ be.
      const GO gblRowInd = meshRowMap.getGlobalElement (lclRowInd);
      for (LO k = 0; k < numEnt; ++k) {
        const GO gblColInd = indexBase +
          ((gblRowInd - indexBase) + static_cast<GO> (k + 1)) %
          static_cast<GO> (globalNumRows);
        lclColInds[k] = meshColMap.getLocalElement (gblColInd);
      }
      // CrsGraph doesn't technically need to promise to sort by local
      // column indices, so we sort both arrays before comparing.
      std::sort (lclColInds.begin (), lclColInds.end ());
      std::copy (myLclColInds, myLclColInds + 2, myLclColIndsSorted.begin ());
      std::sort (myLclColIndsSorted.begin (), myLclColIndsSorted.end ());
      TEST_COMPARE_ARRAYS( lclColInds, myLclColIndsSorted );

      // Test that getLocalRowCopy works.
      size_t numEntries;
      blockMat.getLocalRowCopy (lclRowInd, myLclColIndsCopy(), myValsCopy(), numEntries);
      numEnt = static_cast<LO>(numEntries);
      TEST_ASSERT( err == 0 );
      TEST_ASSERT( numEnt == static_cast<LO> (maxNumEntPerRow) );

      // CrsGraph doesn't technically need to promise to sort by local
      // column indices, so we sort both arrays before comparing.
      std::copy (myLclColIndsCopy.getRawPtr(), myLclColIndsCopy.getRawPtr() + 2, myLclColIndsSorted.begin ());
      std::sort (myLclColIndsSorted.begin (), myLclColIndsSorted.end ());
      TEST_COMPARE_ARRAYS( lclColInds, myLclColIndsSorted );

      // Fill the entries in the row with zeros.
      std::fill (tempBlockSpace.begin (), tempBlockSpace.end (), STS::zero ());
      err = blockMat.replaceLocalValues (lclRowInd, lclColInds.getRawPtr (),
                                         tempBlockSpace.getRawPtr (), numEnt);
      TEST_ASSERT( err == numEnt );
      // Make sure that the input Scalar values didn't change (are
      // still all zero).
      for (LO k = 0; k < numEnt; ++k) {
        Scalar* const tempBlockPtr = tempBlockSpace.getRawPtr () +
          k * blockSize * blockSize;
        little_block_type tempBlock ((typename little_block_type::value_type*) tempBlockPtr, blockSize, blockSize);
        for (LO j = 0; j < blockSize; ++j) {
          for (LO i = 0; i < blockSize; ++i) {
            TEST_ASSERT( static_cast<Scalar> (tempBlock(i,j)) == STS::zero () );
          }
        }
      } // for each entry in the row

      // Create a block pattern which verifies that the matrix stores
      // blocks in row-major order.
      for (LO k = 0; k < numEnt; ++k) {
        Scalar* const tempBlockPtr = tempBlockSpace.getRawPtr () +
          k * blockSize * blockSize;
        little_block_type tempBlock ((typename little_block_type::value_type*) tempBlockPtr, blockSize, blockSize);
        for (LO j = 0; j < blockSize; ++j) {
          for (LO i = 0; i < blockSize; ++i) {
            tempBlock(i,j) = static_cast<Scalar> (static_cast<MT> (j + i * blockSize));
          }
        }
      } // for each entry in the row
      err = blockMat.replaceLocalValues (lclRowInd, lclColInds.getRawPtr (),
                                         tempBlockSpace.getRawPtr (), numEnt);
      TEST_ASSERT( err == numEnt );

      // Get a view of the current row again, and test that the
      // entries were modified as expected.  This tests that the
      // method assumes that the input blocks are row major.
      err = blockMat.getLocalRowView (lclRowInd, myLclColInds, myVals, numEnt);
      TEST_ASSERT( err == 0 );
      TEST_ASSERT( numEnt == static_cast<LO> (maxNumEntPerRow) );
      TEST_ASSERT( myLclColInds != NULL );
      TEST_ASSERT( myVals != NULL );

      for (LO k = 0; k < numEnt; ++k) {
        Scalar* curBlkPtr = myVals + k * blockSize * blockSize;
        little_block_type curBlk ((typename little_block_type::value_type*) curBlkPtr, blockSize, blockSize);

        for (LO j = 0; j < blockSize; ++j) {
          for (LO i = 0; i < blockSize; ++i) {
            TEST_ASSERT( static_cast<Scalar> (curBlk(i,j)) ==
                         static_cast<Scalar> (static_cast<MT> (j + i * blockSize)) );
          }
        }
      } // for each entry in the row
    } // for each local row

    out << "Test applyBlock for a single vector" << endl;

    // Fill a BlockVector and test applyBlock() for a single vector.
    {
      // Each block of the BlockVector will have i-th entry = b - i,
      // where is the block size.  For example, with blockSize = 3, each
      // block will be [3 2 1]^T.  Given how we filled the BlockVector
      // above, this will make each block dense mat-vec give the
      // following results:
      //
      // [0 1 2]   [3]   [ 4]
      // [3 4 5] * [2] = [22]
      // [6 7 8]   [1]   [40]
      //
      // In general, A(i,j) = j + i*b, and X(i) = b - i, so
      //
      // Y(i) = sum_j A(i,j) * X(j)
      //      = sum_j( (j + ib)(b - j) )
      //      = sum_j( b(1 - i)j - j^2 + ib^2 )
      //
      // Since each row of the matrix has two block entries, each block
      // of the result of the global mat-vec will be twice the above
      // result.
      BV X (* (graph.getDomainMap ()), blockSize);
      BV Y (* (graph.getRangeMap ()), blockSize);
      Y.putScalar (STS::zero ());

      const map_type& meshDomainMap = * (graph.getDomainMap ());
      for (LO lclDomIdx = meshDomainMap.getMinLocalIndex ();
           lclDomIdx <= meshDomainMap.getMaxLocalIndex (); ++lclDomIdx) {
        little_vec_type X_lcl = X.getLocalBlock (lclDomIdx);
        TEST_ASSERT( X_lcl.data () != NULL );
        TEST_ASSERT( static_cast<size_t> (X_lcl.extent (0)) == static_cast<size_t> (blockSize) );
        for (LO i = 0; i < blockSize; ++i) {
          X_lcl(i) = static_cast<Scalar> (static_cast<MT> (blockSize - i));
        }
      }

      TEST_NOTHROW( blockMat.applyBlock (X, Y) );

      const map_type& meshRangeMap = * (graph.getRangeMap ());
      for (LO lclRanIdx = meshRangeMap.getMinLocalIndex ();
           lclRanIdx <= meshRangeMap.getMaxLocalIndex (); ++lclRanIdx) {
        little_vec_type Y_lcl = Y.getLocalBlock (lclRanIdx);
        TEST_ASSERT( Y_lcl.data () != NULL );
        TEST_ASSERT( static_cast<size_t> (Y_lcl.extent (0)) == static_cast<size_t> (blockSize) );

        // Test that each actual output value matches its expected value.
        for (LO i = 0; i < blockSize; ++i) {
          // Compute the expected value for the current output index i.  I
          // could have worked out the above formula, but why not let the
          // computer do it?
          LO expectedVal = 0;
          for (LO j = 0; j < blockSize; ++j) {
            expectedVal += blockSize*(1 - i)*j - j*j + i*blockSize*blockSize;
          }
          expectedVal *= static_cast<LO> (2);
          out << "Y_lcl(" << i << ") = " << Y_lcl(i)
              << "; expectedVal = " << expectedVal << std::endl;
          TEST_ASSERT( static_cast<Scalar> (Y_lcl(i)) ==
                       static_cast<Scalar> (static_cast<MT> (expectedVal)) );
        }
      }

      // Repeat this test for alpha = 2 and beta = -3, where we
      // initially fill Y with ones.
      Y.putScalar (STS::one ());
      Scalar alpha = STS::one () + STS::one ();
      Scalar beta = -(STS::one () + STS::one () + STS::one ());

      TEST_NOTHROW( blockMat.applyBlock (X, Y, Teuchos::NO_TRANS, alpha, beta) );

      for (LO lclRanIdx = meshRangeMap.getMinLocalIndex ();
           lclRanIdx <= meshRangeMap.getMaxLocalIndex (); ++lclRanIdx) {
        little_vec_type Y_lcl = Y.getLocalBlock (lclRanIdx);
        TEST_ASSERT( Y_lcl.data () != NULL );
        TEST_ASSERT( static_cast<size_t> (Y_lcl.extent (0)) == static_cast<size_t> (blockSize) );

        // Test that each actual output value matches its expected value.
        for (LO i = 0; i < blockSize; ++i) {
          // Compute the expected value for the current output index i.  I
          // could have worked out the above formula, but why not let the
          // computer do it?
          LO expectedVal = 0;
          for (LO j = 0; j < blockSize; ++j) {
            expectedVal += blockSize*(1 - i)*j - j*j + i*blockSize*blockSize;
          }
          // Taking the real part doesn't do anything here, since
          // alpha is real.  However, it returns MT, which is safe to
          // cast to LO, even if Scalar is complex.
          expectedVal *= (static_cast<LO> (2) * static_cast<LO> (STS::real (alpha)));
          // See above note about taking the real part.
          expectedVal += static_cast<LO> (STS::real (beta));
          out << "Y_lcl(" << i << ") = " << Y_lcl(i)
              << "; expectedVal = " << expectedVal << std::endl;
          TEST_ASSERT( static_cast<Scalar> (Y_lcl(i)) ==
                       static_cast<Scalar> (static_cast<MT> (expectedVal)) );
        }
      }
    } // done with single-vector applyBlock test

    out << "Test applyBlock for multiple vectors (right-hand sides)" << endl;

    // Fill a BlockMultiVector and test applyBlock() for multiple vectors.
    {
      const LO numVecs = 3;

      // Each block of the BlockMultiVector will have i-th entry = b -
      // i, where is the block size.  For example, with blockSize = 3,
      // each block will be [3 2 1]^T.  Given how we filled the
      // BlockVector above, this will make each block dense mat-vec
      // give the following results:
      //
      // [0 1 2]   [3]   [ 4]
      // [3 4 5] * [2] = [22]
      // [6 7 8]   [1]   [40]
      //
      // In general, A(i,j) = j + i*b, and X(i) = b - i, so
      //
      // Y(i) = sum_j A(i,j) * X(j)
      //      = sum_j( (j + ib)(b - j) )
      //      = sum_j( b(1 - i)j - j^2 + ib^2 )
      //
      // Since each row of the matrix has two block entries, each block
      // of the result of the global mat-vec will be twice the above
      // result.
      //
      // For the multiple-vector case, we revise the above test so
      // that column j of X is scaled by j+1.
      BMV X (* (graph.getDomainMap ()), blockSize, numVecs);
      BMV Y (* (graph.getRangeMap ()), blockSize, numVecs);
      Y.putScalar (STS::zero ());

      const map_type& meshDomainMap = * (graph.getDomainMap ());
      for (LO lclDomIdx = meshDomainMap.getMinLocalIndex ();
           lclDomIdx <= meshDomainMap.getMaxLocalIndex (); ++lclDomIdx) {
        for (LO j = 0; j < numVecs; ++j) {
          little_vec_type X_lcl = X.getLocalBlock (lclDomIdx, j);
          TEST_ASSERT( X_lcl.data () != NULL );
          TEST_ASSERT( static_cast<size_t> (X_lcl.extent (0)) == static_cast<size_t> (blockSize) );
          for (LO i = 0; i < blockSize; ++i) {
            X_lcl(i) = static_cast<Scalar> (static_cast<MT> ((blockSize - i) * (j + 1)));
          }
        }
      }

      TEST_NOTHROW( blockMat.applyBlock (X, Y) );

      const map_type& meshRangeMap = * (graph.getRangeMap ());
      for (LO lclRanIdx = meshRangeMap.getMinLocalIndex ();
           lclRanIdx <= meshRangeMap.getMaxLocalIndex (); ++lclRanIdx) {
        for (LO col = 0; col < numVecs; ++col) {
          little_vec_type Y_lcl = Y.getLocalBlock (lclRanIdx, col);
          TEST_ASSERT( Y_lcl.data () != NULL );
          TEST_ASSERT( static_cast<size_t> (Y_lcl.extent (0)) == static_cast<size_t> (blockSize) );

          // Test that each actual output value matches its expected value.
          for (LO i = 0; i < blockSize; ++i) {
            // Compute the expected value for the current output index i.  I
            // could have worked out the above formula, but why not let the
            // computer do it?
            LO expectedVal = 0;
            for (LO j = 0; j < blockSize; ++j) {
              expectedVal += blockSize*(1 - i)*j - j*j + i*blockSize*blockSize;
            }
            expectedVal *= static_cast<LO> (2);
            expectedVal *= static_cast<LO> (col + 1);
            out << "Y_lcl(" << i << ") = " << Y_lcl(i)
                << "; expectedVal = " << expectedVal << std::endl;
            TEST_ASSERT( static_cast<Scalar> (Y_lcl(i)) ==
                         static_cast<Scalar> (static_cast<MT> (expectedVal)) );
          }
        }
      }

      // Repeat this test for alpha = 2 and beta = -3, where we
      // initially fill Y with ones.
      Y.putScalar (STS::one ());
      Scalar alpha = two;
      Scalar beta = -three;

      TEST_NOTHROW( blockMat.applyBlock (X, Y, Teuchos::NO_TRANS, alpha, beta) );

      for (LO lclRanIdx = meshRangeMap.getMinLocalIndex ();
           lclRanIdx <= meshRangeMap.getMaxLocalIndex (); ++lclRanIdx) {
        for (LO col = 0; col < numVecs; ++col) {
          little_vec_type Y_lcl = Y.getLocalBlock (lclRanIdx, col);
          TEST_ASSERT( Y_lcl.data () != NULL );
          TEST_ASSERT( static_cast<size_t> (Y_lcl.extent (0)) == static_cast<size_t> (blockSize) );

          // Test that each actual output value matches its expected value.
          for (LO i = 0; i < blockSize; ++i) {
            // Compute the expected value for the current output index i.  I
            // could have worked out the above formula, but why not let the
            // computer do it?
            LO expectedVal = 0;
            for (LO j = 0; j < blockSize; ++j) {
              expectedVal += blockSize*(1 - i)*j - j*j + i*blockSize*blockSize;
            }
            // Taking the real part doesn't do anything here, since
            // alpha is real.  However, it returns MT, which is safe to
            // cast to LO, even if Scalar is complex.
            expectedVal *= (static_cast<LO> (2) * static_cast<LO> (STS::real (alpha)));
            expectedVal *= static_cast<LO> (col + 1);
            // See above note about taking the real part.
            expectedVal += static_cast<LO> (STS::real (beta));
            out << "Y_lcl(" << i << ") = " << Y_lcl(i)
                << "; expectedVal = " << expectedVal << std::endl;
            TEST_ASSERT( static_cast<Scalar> (Y_lcl(i)) ==
                         static_cast<Scalar> (static_cast<MT> (expectedVal)) );
          }
        } // for each column (vector) of the BlockMultiVector
      } // for each local (mesh) row of the output BlockMultiVector
    } // done with multiple-vector applyBlock test

    out << "Test apply for a single vector" << endl;

    // Fill a BlockVector and test apply() for a single vector.
    {
      // Each block of the BlockVector will have i-th entry = b - i,
      // where is the block size.  For example, with blockSize = 3,
      // each block will be [3 2 1]^T.  Given how we filled the
      // BlockVector above, this will make each block dense mat-vec
      // give the following results:
      //
      // [0 1 2]   [3]   [ 4]
      // [3 4 5] * [2] = [22]
      // [6 7 8]   [1]   [40]
      //
      // In general, A(i,j) = j + i*b, and X(i) = b - i, so
      //
      // Y(i) = sum_j A(i,j) * X(j)
      //      = sum_j( (j + ib)(b - j) )
      //      = sum_j( b(1 - i)j - j^2 + ib^2 )
      //
      // Since each row of the matrix has two block entries, each block
      // of the result of the global mat-vec will be twice the above
      // result.
      BV X (* (graph.getDomainMap ()), blockSize);
      BV Y (* (graph.getRangeMap ()), blockSize);
      Y.putScalar (STS::zero ());

      const map_type& meshDomainMap = * (graph.getDomainMap ());
      for (LO lclDomIdx = meshDomainMap.getMinLocalIndex ();
           lclDomIdx <= meshDomainMap.getMaxLocalIndex (); ++lclDomIdx) {
        little_vec_type X_lcl = X.getLocalBlock (lclDomIdx);
        TEST_ASSERT( X_lcl.data () != NULL );
        TEST_ASSERT( static_cast<size_t> (X_lcl.extent (0)) == static_cast<size_t> (blockSize) );
        for (LO i = 0; i < blockSize; ++i) {
          X_lcl(i) = static_cast<Scalar> (static_cast<MT> (blockSize - i));
        }
      }

      vec_type X_vec = X.getVectorView ();
      vec_type Y_vec = Y.getVectorView ();

      TEST_NOTHROW( blockMat.apply (X_vec, Y_vec) );

      // This test also exercises whether getVectorView really does
      // return a view, since we access and test results using the
      // BlockVector.

      const map_type& meshRangeMap = * (graph.getRangeMap ());
      for (LO lclRanIdx = meshRangeMap.getMinLocalIndex ();
           lclRanIdx <= meshRangeMap.getMaxLocalIndex (); ++lclRanIdx) {
        little_vec_type Y_lcl = Y.getLocalBlock (lclRanIdx);
        TEST_ASSERT( Y_lcl.data () != NULL );
        TEST_ASSERT( static_cast<size_t> (Y_lcl.extent (0)) == static_cast<size_t> (blockSize) );

        // Test that each actual output value matches its expected value.
        for (LO i = 0; i < blockSize; ++i) {
          // Compute the expected value for the current output index i.  IL
          // could have worked out the above formula, but why not let the
          // computer do it?
          LO expectedVal = 0;
          for (LO j = 0; j < blockSize; ++j) {
            expectedVal += blockSize*(1 - i)*j - j*j + i*blockSize*blockSize;
          }
          expectedVal *= static_cast<LO> (2);
          out << "Y_lcl(" << i << ") = " << Y_lcl(i)
              << "; expectedVal = " << expectedVal << std::endl;
          TEST_ASSERT( static_cast<Scalar> (Y_lcl(i)) ==
                       static_cast<Scalar> (static_cast<MT> (expectedVal)) );
        }
      }

      // Repeat this test for alpha = 2 and beta = -3, where we
      // initially fill Y with ones.
      Y.putScalar (STS::one ());
      Scalar alpha = STS::one () + STS::one ();
      Scalar beta = -(STS::one () + STS::one () + STS::one ());

      TEST_NOTHROW( blockMat.apply (X_vec, Y_vec, Teuchos::NO_TRANS, alpha, beta) );

      for (LO lclRanIdx = meshRangeMap.getMinLocalIndex ();
           lclRanIdx <= meshRangeMap.getMaxLocalIndex (); ++lclRanIdx) {
        little_vec_type Y_lcl = Y.getLocalBlock (lclRanIdx);
        TEST_ASSERT( Y_lcl.data () != NULL );
        TEST_ASSERT( static_cast<size_t> (Y_lcl.extent (0)) == static_cast<size_t> (blockSize) );

        // Test that each actual output value matches its expected value.
        for (LO i = 0; i < blockSize; ++i) {
          // Compute the expected value for the current output index i.  I
          // could have worked out the above formula, but why not let the
          // computer do it?
          LO expectedVal = 0;
          for (LO j = 0; j < blockSize; ++j) {
            expectedVal += blockSize*(1 - i)*j - j*j + i*blockSize*blockSize;
          }
          // Taking the real part doesn't do anything here, since
          // alpha is real.  However, it returns MT, which is safe to
          // cast to LO, even if Scalar is complex.
          expectedVal *= (static_cast<LO> (2) * static_cast<LO> (STS::real (alpha)));
          // See above note about taking the real part.
          expectedVal += static_cast<LO> (STS::real (beta));
          out << "Y_lcl(" << i << ") = " << Y_lcl(i)
              << "; expectedVal = " << expectedVal << std::endl;
          TEST_ASSERT( static_cast<Scalar> (Y_lcl(i)) ==
                       static_cast<Scalar> (static_cast<MT> (expectedVal)) );
        }
      }
    } // done with single-vector apply test

    out << "Test apply for multiple vectors (right-hand sides)" << endl;

    // Fill a BlockMultiVector and test apply() for multiple vectors.
    {
      const LO numVecs = 3;

      // Each block of the BlockMultiVector will have i-th entry = b -
      // i, where is the block size.  For example, with blockSize = 3,
      // each block will be [3 2 1]^T.  Given how we filled the
      // BlockVector above, this will make each block dense mat-vec
      // give the following results:
      //
      // [0 1 2]   [3]   [ 4]
      // [3 4 5] * [2] = [22]
      // [6 7 8]   [1]   [40]
      //
      // In general, A(i,j) = j + i*b, and X(i) = b - i, so
      //
      // Y(i) = sum_j A(i,j) * X(j)
      //      = sum_j( (j + ib)(b - j) )
      //      = sum_j( b(1 - i)j - j^2 + ib^2 )
      //
      // Since each row of the matrix has two block entries, each block
      // of the result of the global mat-vec will be twice the above
      // result.
      //
      // For the multiple-vector case, we revise the above test so
      // that column j of X is scaled by j+1.
      BMV X (* (graph.getDomainMap ()), blockSize, numVecs);
      BMV Y (* (graph.getRangeMap ()), blockSize, numVecs);
      Y.putScalar (STS::zero ());

      const map_type& meshDomainMap = * (graph.getDomainMap ());
      for (LO lclDomIdx = meshDomainMap.getMinLocalIndex ();
           lclDomIdx <= meshDomainMap.getMaxLocalIndex (); ++lclDomIdx) {
        for (LO j = 0; j < numVecs; ++j) {
          little_vec_type X_lcl = X.getLocalBlock (lclDomIdx, j);
          TEST_ASSERT( X_lcl.data () != NULL );
          TEST_ASSERT( static_cast<size_t> (X_lcl.extent (0)) == static_cast<size_t> (blockSize) );
          for (LO i = 0; i < blockSize; ++i) {
            X_lcl(i) = static_cast<Scalar> (static_cast<MT> ((blockSize - i) * (j + 1)));
          }
        }
      }

      mv_type X_mv = X.getMultiVectorView ();
      mv_type Y_mv = Y.getMultiVectorView ();

      TEST_NOTHROW( blockMat.apply (X_mv, Y_mv) );

      // This test also exercises whether getMultiVectorView really
      // does return a view, since we access and test results using
      // the BlockMultiVector.

      const map_type& meshRangeMap = * (graph.getRangeMap ());
      for (LO lclRanIdx = meshRangeMap.getMinLocalIndex ();
           lclRanIdx <= meshRangeMap.getMaxLocalIndex (); ++lclRanIdx) {
        for (LO col = 0; col < numVecs; ++col) {
          little_vec_type Y_lcl = Y.getLocalBlock (lclRanIdx, col);
          TEST_ASSERT( Y_lcl.data () != NULL );
          TEST_ASSERT( static_cast<size_t> (Y_lcl.extent (0)) == static_cast<size_t> (blockSize) );

          // Test that each actual output value matches its expected value.
          for (LO i = 0; i < blockSize; ++i) {
            // Compute the expected value for the current output index i.  I
            // could have worked out the above formula, but why not let the
            // computer do it?
            LO expectedVal = 0;
            for (LO j = 0; j < blockSize; ++j) {
              expectedVal += blockSize*(1 - i)*j - j*j + i*blockSize*blockSize;
            }
            expectedVal *= static_cast<LO> (2);
            expectedVal *= static_cast<LO> (col + 1);
            out << "Y_lcl(" << i << ") = " << Y_lcl(i)
                << "; expectedVal = " << expectedVal << std::endl;
            TEST_ASSERT( static_cast<Scalar> (Y_lcl(i)) ==
                         static_cast<Scalar> (static_cast<MT> (expectedVal)) );
          }
        }
      }

      // Repeat this test for alpha = 2 and beta = -3, where we
      // initially fill Y with ones.
      Y.putScalar (STS::one ());
      Scalar alpha = two;
      Scalar beta = -three;

      TEST_NOTHROW( blockMat.apply (X_mv, Y_mv, Teuchos::NO_TRANS, alpha, beta) );

      for (LO lclRanIdx = meshRangeMap.getMinLocalIndex ();
           lclRanIdx <= meshRangeMap.getMaxLocalIndex (); ++lclRanIdx) {
        for (LO col = 0; col < numVecs; ++col) {
          little_vec_type Y_lcl = Y.getLocalBlock (lclRanIdx, col);
          TEST_ASSERT( Y_lcl.data () != NULL );
          TEST_ASSERT( static_cast<size_t> (Y_lcl.extent (0)) == static_cast<size_t> (blockSize) );

          // Test that each actual output value matches its expected value.
          for (LO i = 0; i < blockSize; ++i) {
            // Compute the expected value for the current output index i.  I
            // could have worked out the above formula, but why not let the
            // computer do it?
            LO expectedVal = 0;
            for (LO j = 0; j < blockSize; ++j) {
              expectedVal += blockSize*(1 - i)*j - j*j + i*blockSize*blockSize;
            }
            // Taking the real part doesn't do anything here, since
            // alpha is real.  However, it returns MT, which is safe to
            // cast to LO, even if Scalar is complex.
            expectedVal *= (static_cast<LO> (2) * static_cast<LO> (STS::real (alpha)));
            expectedVal *= static_cast<LO> (col + 1);
            // See above note about taking the real part.
            expectedVal += static_cast<LO> (STS::real (beta));
            out << "Y_lcl(" << i << ") = " << Y_lcl(i)
                << "; expectedVal = " << expectedVal << std::endl;
            TEST_ASSERT( static_cast<Scalar> (Y_lcl(i)) ==
                         static_cast<Scalar> (static_cast<MT> (expectedVal)) );
          }
        } // for each column (vector) of the BlockMultiVector
      } // for each local (mesh) row of the output BlockMultiVector
    } // done with multiple-vector apply test

    // Finishing with an all-reduce ensures that the test won't run to
    // completion on Process 0 (and therefore report a "pass") if
    // there is deadlock somewhere.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    out << "Hooray, got to the end of the BlockCrsMatrix test!" << std::endl;
    TEST_EQUALITY_CONST( gblSuccess, 1 );
  }

  // Test writing of a BlockCrsMatrix.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockCrsMatrix, write, Scalar, LO, GO, Node )
  {
    typedef Tpetra::BlockCrsMatrix<Scalar, LO, GO, Node> BCM;
    typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    // The typedef below is also a test.  BlockCrsMatrix must have
    // this typedef, or this test won't compile.
    typedef typename BCM::little_block_type little_block_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType MT;

    out << "Testing output of a Tpetra::BlockCrsMatrix" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const Comm<int> > comm = getDefaultComm ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    out << "Creating mesh row Map" << endl;

    const size_t numLocalMeshPoints = 12;
    const GO indexBase = 1;
    // Use a block size that is not a power of 2, to test correctness
    // in case the matrix pads blocks for SIMD-ization.
    const LO blockSize = 3;
    const size_t entriesPerBlock = blockSize * blockSize;

    // mfh 20 May 2014: Tpetra::CrsGraph still needs the row Map as an
    // RCP.  Later interface changes will let us pass in the Map by
    // const reference and assume view semantics.
    RCP<const map_type> meshRowMapPtr =
      rcp (new map_type (INVALID, numLocalMeshPoints, indexBase, comm));
    const map_type& meshRowMap = *meshRowMapPtr;

    // Make a graph.  It will have two entries per global row i_gbl:
    // (i_gbl, (i_gbl+1) % N) and (i_gbl, (i_gbl+2) % N), where N is
    // the global number of rows and columns.  We don't include the
    // diagonal, to make the Maps more interesting.
    out << "Creating mesh graph" << endl;

    const size_t maxNumEntPerRow = 2;
    graph_type graph (meshRowMapPtr, maxNumEntPerRow, Tpetra::StaticProfile);

    // Fill the graph.
    Teuchos::Array<GO> gblColInds (maxNumEntPerRow);
    const GO globalNumRows = meshRowMap.getGlobalNumElements ();
    for (LO lclRowInd = meshRowMap.getMinLocalIndex ();
         lclRowInd <= meshRowMap.getMaxLocalIndex (); ++lclRowInd) {
      const GO gblRowInd = meshRowMap.getGlobalElement (lclRowInd);
      for (size_t k = 0; k < maxNumEntPerRow; ++k) {
        const GO gblColInd = indexBase +
          ((gblRowInd - indexBase) + static_cast<GO> (k + 1)) % static_cast<GO> (globalNumRows);
        gblColInds[k] = gblColInd;
      }
      graph.insertGlobalIndices (gblRowInd, gblColInds ());
    }
    graph.fillComplete ();

    // Get the graph's column Map (the "mesh column Map").
    map_type meshColMap = * (graph.getColMap ());

    out << "Creating BlockCrsMatrix" << endl;

    // Construct the BlockCrsMatrix.
    BCM blockMat = BCM (graph, blockSize);

    out << "Test getLocalRowView and replaceLocalValues" << endl;

    Array<Scalar> tempBlockSpace (maxNumEntPerRow * entriesPerBlock);

    // Test that getLocalRowView returns the right column indices.
    Array<LO> lclColInds (maxNumEntPerRow);
    Array<LO> myLclColIndsCopy (maxNumEntPerRow);
    for (LO lclRowInd = meshRowMap.getMinLocalIndex ();
         lclRowInd <= meshRowMap.getMaxLocalIndex (); ++lclRowInd) {
      const LO* myLclColInds = NULL;
      Scalar* myVals = NULL;
      LO numEnt = 0;
      blockMat.getLocalRowView (lclRowInd, myLclColInds, myVals, numEnt);

      // Compute what the local column indices in this row _should_ be.
      const GO gblRowInd = meshRowMap.getGlobalElement (lclRowInd);
      for (LO k = 0; k < numEnt; ++k) {
        const GO gblColInd = indexBase +
          ((gblRowInd - indexBase) + static_cast<GO> (k + 1)) %
          static_cast<GO> (globalNumRows);
        lclColInds[k] = meshColMap.getLocalElement (gblColInd);
      }
      // CrsGraph doesn't technically need to promise to sort by local
      // column indices, so we sort both arrays before comparing.
      std::sort (lclColInds.begin (), lclColInds.end ());
      std::copy (myLclColInds, myLclColInds + 2, myLclColIndsCopy.begin ());
      std::sort (myLclColIndsCopy.begin (), myLclColIndsCopy.end ());
      TEST_COMPARE_ARRAYS( lclColInds, myLclColIndsCopy );

      // Fill the entries in the row with zeros.
      std::fill (tempBlockSpace.begin (), tempBlockSpace.end (), STS::zero ());
      blockMat.replaceLocalValues (lclRowInd, lclColInds.getRawPtr (),
                                         tempBlockSpace.getRawPtr (), numEnt);

      // Create a block pattern which verifies that the matrix stores
      // blocks in row-major order.
      for (LO k = 0; k < numEnt; ++k) {
        Scalar* const tempBlockPtr = tempBlockSpace.getRawPtr () +
          k * blockSize * blockSize;
        little_block_type tempBlock ((typename little_block_type::value_type*) tempBlockPtr, blockSize, blockSize);
        for (LO j = 0; j < blockSize; ++j) {
          for (LO i = 0; i < blockSize; ++i) {
            tempBlock(i,j) = static_cast<Scalar> (static_cast<MT> (j + i * blockSize) + 0.0123);
          }
        }
      } // for each entry in the row
      blockMat.replaceLocalValues (lclRowInd, lclColInds.getRawPtr (),
                                         tempBlockSpace.getRawPtr (), numEnt);
    } // for each local row

    blockMat.describe(out,Teuchos::VERB_EXTREME);

    Teuchos::ParameterList pl;
    //pl.set("precision",3);
    //pl.set("zero-based indexing",true);
    //pl.set("always use parallel algorithm",true);
    //pl.set("print MatrixMarket header",false);
    Tpetra::Experimental::blockCrsMatrixWriter(blockMat, "savedBlockMatrix.m", pl);
  }


  // Test BlockCrsMatrix::getLocalDiagCopy.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockCrsMatrix, getLocalDiagCopy, Scalar, LO, GO, Node )
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    typedef Tpetra::BlockCrsMatrix<Scalar, LO, GO, Node> BCM;
    typedef typename BCM::device_type device_type;
    typedef typename BCM::impl_scalar_type IST;
    typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    // The typedef below is also a test.  BlockCrsMatrix must have
    // this typedef, or this test won't compile.
    typedef typename BCM::little_block_type little_block_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;

    int lclSuccess = 1;
    int gblSuccess = 1;

    out << "Test Tpetra::BlockCrsMatrix::getLocalDiagCopy" << endl;
    Teuchos::OSTab tab0 (out);

    out << "Create mesh row Map" << endl;
    RCP<const Comm<int> > comm = getDefaultComm ();
    const size_t numLclMeshPoints = 7;
    const GST numGblMeshPoints = comm->getSize () * numLclMeshPoints;
    const GO indexBase = 0;

    // mfh 12 Dec 2015: Tpetra::CrsGraph still needs the row Map as an
    // RCP.  Later interface changes will let us pass in the Map by
    // const reference and assume view semantics.
    RCP<const map_type> meshRowMapPtr =
      rcp (new map_type (numGblMeshPoints, numLclMeshPoints, indexBase, comm));
    const map_type& meshRowMap = *meshRowMapPtr;

    // Make a graph.  It will have two entries per global row i_gbl:
    // the diagonal entry (i_gbl, i_gbl), and the off-diagonal entry
    // (i_gbl, (i_gbl+1) % N), where N is the global number of
    // columns.  The mod ensures that the relative order of the
    // diagonal and off-diagonal entries differs in different rows (so
    // that getLocalDiagOffsets actually has to search).
    out << "Creating mesh graph" << endl;

    const size_t maxNumEntPerRow = 2;
    graph_type graph (meshRowMapPtr, maxNumEntPerRow, Tpetra::StaticProfile);

    // Fill the graph.
    Teuchos::Array<GO> gblColInds (maxNumEntPerRow);
    const GO globalNumRows = meshRowMap.getGlobalNumElements ();
    for (LO lclRowInd = meshRowMap.getMinLocalIndex ();
         lclRowInd <= meshRowMap.getMaxLocalIndex (); ++lclRowInd) {
      const GO gblRowInd = meshRowMap.getGlobalElement (lclRowInd);
      gblColInds[0] = gblRowInd; // diagonal entry
      const GO gblOffDiagColInd = indexBase +
        ((gblRowInd - indexBase) + static_cast<GO> (1)) % static_cast<GO> (globalNumRows);
      gblColInds[1] = gblOffDiagColInd; // off-diagonal entry
      graph.insertGlobalIndices (gblRowInd, gblColInds ());
    }
    graph.fillComplete ();

    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    // Get the graph's column Map (the "mesh column Map").
    map_type meshColMap = * (graph.getColMap ());

    // Use a block size that is not a power of 2, to test correctness
    // in case BlockCrsMatrix ever were to pad blocks for SIMD.
    const LO blockSize = 3;
    //const size_t entriesPerBlock = blockSize * blockSize;

    // Construct the BlockCrsMatrix.
    out << "Create BlockCrsMatrix" << endl;
    BCM blockMat = BCM (graph, blockSize);

    // Get the (mesh) offsets of the diagonal blocks; we'll need them
    // later for getLocalDiagCopy.
    typedef typename Node::device_type DT;
    Kokkos::View<size_t*, DT> diagMeshOffsets ("offsets", numLclMeshPoints);
    try {
      graph.getLocalDiagOffsets (diagMeshOffsets);
    } catch (std::exception& e) {
      success = false;
      std::ostringstream os;
      os << "Proc " << comm->getRank () << ": getLocalDiagOffsets "
        "threw an exception: " << e.what () << endl;
      std::cerr << os.str ();
    }

    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "Test FAILED before or at getLocalDiagOffsets call; "
        "doesn't make sense to continue." << endl;
      return;
    }

    // Test that the mesh offsets are correct.  Use the Kokkos local
    // graph to test this.  (Mesh offsets are "locally absolute" --
    // this is technically an implementation detail (the result of
    // getLocalDiagOffsets is supposed to be opaque) but I think it's
    // OK to use here.)
    bool meshOffsetsCorrect = true;
    if (numLclMeshPoints == 0) {
      TEST_EQUALITY( diagMeshOffsets.extent (0), 0 );
      if (diagMeshOffsets.extent (0) != 0) {
        meshOffsetsCorrect = false;
      }
    }
    else {
      TEST_ASSERT( diagMeshOffsets.extent (0) != 0 );
      auto localGraph = graph.getLocalGraph ();
      const auto& colMap = * (graph.getColMap ());

      TEST_EQUALITY( static_cast<size_t> (numLclMeshPoints + 1),
                     static_cast<size_t> (localGraph.row_map.extent (0)) );
      if (static_cast<size_t> (numLclMeshPoints + 1) ==
          static_cast<size_t> (localGraph.row_map.extent (0))) {
        for (LO lclRowInd = 0;
             lclRowInd < static_cast<LO> (numLclMeshPoints);
             ++lclRowInd) {
          const GO gblRowInd = meshRowMap.getGlobalElement (lclRowInd);
          const GO gblColInd = gblRowInd;
          bool diagOffsetCorrect = false;

          const LO* lclColInds = NULL;
          Scalar* lclVals = NULL;
          LO numEnt = 0;
          LO err = blockMat.getLocalRowView (lclRowInd, lclColInds, lclVals, numEnt);
          TEST_ASSERT( err == 0 );
          if (err == 0) {
            const size_t offset = diagMeshOffsets[lclRowInd];
            if (offset >= static_cast<size_t> (numEnt)) {
              diagOffsetCorrect = false;
            }
            else {
              const LO actualLclColInd = lclColInds[offset];
              const GO actualGblColInd = colMap.getGlobalElement (actualLclColInd);
              diagOffsetCorrect = (actualGblColInd == gblColInd);
            }
          }

          TEST_ASSERT( diagOffsetCorrect );
          if (! diagOffsetCorrect ) {
            meshOffsetsCorrect = false;
          }
        }
      }
    }
    TEST_ASSERT( meshOffsetsCorrect );

    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "Test FAILED early; doesn't make sense to continue." << endl;
      return;
    }

    // FIXME (mfh 13 Dec 2015) If the mesh offsets are wrong on some
    // process, then we should stop the test early.

    // Fill the BlockCrsMatrix in such a way that we can easily tell,
    // by inspecting the values,
    //
    //   1. which block is the diagonal block, and
    //   2. whether the layout of the returned blocks is correct.
    out << "Fill the BlockCrsMatrix" << endl;
    for (LO lclRowInd = meshRowMap.getMinLocalIndex ();
         lclRowInd <= meshRowMap.getMaxLocalIndex (); ++lclRowInd) {
      const GO gblRowInd = meshRowMap.getGlobalElement (lclRowInd);
      const LO* lclColInds = NULL;
      Scalar* myVals = NULL;
      LO numEnt = 0;
      blockMat.getLocalRowView (lclRowInd, lclColInds, myVals, numEnt);

      // Fill the diagonal block D such that D(i,j) = (lclRowInd+1) *
      // (1 + i + j*blockSize).  Fill the off-diagonal block with -1.
      // This ensures that we can tell we got the right blocks, and
      // that we copied them in the correct order.
      for (LO k = 0; k < numEnt; ++k) {
        const LO offset = blockSize * blockSize * k;
        little_block_type curBlock (reinterpret_cast<IST*> (myVals) + offset,
                                    blockSize, blockSize); // row major
        const GO gblColInd = meshColMap.getGlobalElement (lclColInds[k]);
        if (gblColInd == gblRowInd) { // the diagonal block
          IST curVal = STS::one ();
          for (LO j = 0; j < blockSize; ++j) {
            for (LO i = 0; i < blockSize; ++i) {
              curBlock(i,j) = curVal;
              curVal += static_cast<IST> (STS::one ());
            }
          }
        }
        else { // not the diagonal block
          Kokkos::deep_copy (curBlock, static_cast<IST> (-STS::one ()));
        }
      }
    } // for each local mesh row

    // Now that we've filled the BlockCrsMatrix, use the previously
    // computed offsets to get the local diagonal copy.
    typedef Kokkos::View<IST***, device_type> diag_blocks_type;
    diag_blocks_type diagBlocks ("diagBlocks", numLclMeshPoints,
                                 blockSize, blockSize);
    blockMat.getLocalDiagCopy (diagBlocks, diagMeshOffsets);

    bool allBlocksGood = true;
    for (LO lclRowInd = 0; lclRowInd < static_cast<LO> (numLclMeshPoints); ++lclRowInd) {
      const GO gblRowInd = meshRowMap.getGlobalElement (lclRowInd);
      const LO* lclColInds = NULL;
      Scalar* myVals = NULL;
      LO numEnt = 0;
      blockMat.getLocalRowView (lclRowInd, lclColInds, myVals, numEnt);

      // Make sure that the diagonal blocks from getLocalDiagCopy
      // match those in the matrix.
      for (LO k = 0; k < numEnt; ++k) {
        const LO offset = blockSize * blockSize * k;
        little_block_type curBlock (reinterpret_cast<IST*> (myVals) + offset,
                                    blockSize, blockSize); // row major
        const GO gblColInd = meshColMap.getGlobalElement (lclColInds[k]);
        if (gblColInd == gblRowInd) { // the diagonal block
          auto diagBlock = subview (diagBlocks, lclRowInd, ALL (), ALL ());
          for (LO j = 0; j < blockSize; ++j) {
            for (LO i = 0; i < blockSize; ++i) {
              if (curBlock(i,j) != diagBlock(i,j)) {
                allBlocksGood = false;
                break;
              }
            }
          }
        }
      }
    } // for each local mesh row

    TEST_ASSERT( allBlocksGood );
  }



  // Test BlockCrsMatrix::setAllToScalar.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockCrsMatrix, SetAllToScalar, Scalar, LO, GO, Node )
  {
    typedef Tpetra::BlockMultiVector<Scalar, LO, GO, Node> BMV;
    typedef Tpetra::BlockCrsMatrix<Scalar, LO, GO, Node> BCM;
    typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;

    out << "Testing Tpetra::BlockCrsMatrix::setAllToScalar" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const Comm<int> > comm = getDefaultComm ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    out << "Creating mesh row Map" << endl;

    const size_t numLocalMeshPoints = 12;
    const GO indexBase = 1;
    // mfh 16 May 2014: Tpetra::CrsGraph still needs the row Map as an
    // RCP.  Later interface changes will let us pass in the Map by
    // const reference and assume view semantics.
    RCP<const map_type> meshRowMapPtr =
      rcp (new map_type (INVALID, numLocalMeshPoints, indexBase, comm));
    const GO numGlobalMeshPoints = meshRowMapPtr->getGlobalNumElements ();
    const LO blockSize = 4;

    // Make a graph.  It happens to have two entries per row.
    out << "Creating mesh graph" << endl;
    graph_type graph (meshRowMapPtr, 2, Tpetra::StaticProfile);

    if (meshRowMapPtr->getNodeNumElements () > 0) {
      const GO myMinGblRow = meshRowMapPtr->getMinGlobalIndex ();
      const GO myMaxGblRow = meshRowMapPtr->getMaxGlobalIndex ();
      for (GO gblRow = myMinGblRow; gblRow <= myMaxGblRow; ++gblRow) {
        // Insert two entries, neither of which are on the diagonal.
        Teuchos::Array<GO> gblCols (2);
        gblCols[0] = indexBase + ((gblRow - indexBase) + 1) % numGlobalMeshPoints;
        gblCols[1] = indexBase + ((gblRow - indexBase) + 2) % numGlobalMeshPoints;
        graph.insertGlobalIndices (gblRow, gblCols ());
      }
    }
    graph.fillComplete ();

    out << "Calling BlockCrsMatrix two-argument constructor" << endl;
    BCM blockMat (graph, blockSize);

    // Test that the point domain and range Maps are correct.
    out << "Test the matrix's point domain and range Maps" << endl;
    map_type pointDomainMap = BMV::makePointMap (* (graph.getDomainMap ()), blockSize);
    TEST_ASSERT( ! blockMat.getDomainMap ().is_null () &&
                 pointDomainMap.isSameAs (* (blockMat.getDomainMap ())) );
    map_type pointRangeMap = BMV::makePointMap (* (graph.getRangeMap ()), blockSize);
    TEST_ASSERT( ! blockMat.getRangeMap ().is_null () &&
                 pointRangeMap.isSameAs (* (blockMat.getRangeMap ())) );

    // Test that the result of getGraph() has the same Maps.
    out << "Test that the result of getGraph() has the same Maps "
      "as the original graph" << endl;
    {
      graph_type graph2 = blockMat.getCrsGraph ();
      TEST_ASSERT( ! graph.getDomainMap ().is_null () &&
                   ! graph2.getDomainMap ().is_null () &&
                   graph.getDomainMap ()->isSameAs (* (graph2.getDomainMap ())) );
      TEST_ASSERT( ! graph.getRangeMap ().is_null () &&
                   ! graph2.getRangeMap ().is_null () &&
                   graph.getRangeMap ()->isSameAs (* (graph2.getRangeMap ())) );
      TEST_ASSERT( ! graph.getRowMap ().is_null () &&
                   ! graph2.getRowMap ().is_null () &&
                   graph.getRowMap ()->isSameAs (* (graph2.getRowMap ())) );
      TEST_ASSERT( ! graph.getColMap ().is_null () &&
                   ! graph2.getColMap ().is_null () &&
                   graph.getColMap ()->isSameAs (* (graph2.getColMap ())) );
    }

    // Fill all entries of the matrix with 3.
    out << "Fill all entries of the matrix with 3 (setAllToScalar)" << endl;
    const Scalar three = STS::one () + STS::one () + STS::one ();
    blockMat.setAllToScalar (three);

    // Y := A*X, where X is a block multivector (with one column) full
    // of 1s.  Since there are two block entries per row, each of
    // which is all 3s, we know that each entry of the result Y will
    // be 6*blockSize.
    out << "Test applyBlock" << endl;
    const Scalar requiredValue = static_cast<Scalar> (6 * blockSize);
    BMV X (* (graph.getDomainMap ()), pointDomainMap, blockSize, static_cast<LO> (1));
    X.putScalar (STS::one ());
    BMV Y (* (graph.getRangeMap ()), pointRangeMap, blockSize, static_cast<LO> (1));
    blockMat.applyBlock (X, Y, Teuchos::NO_TRANS, STS::one (), STS::zero ());

    out << "Make sure applyBlock got the right answer" << endl;
    const LO myMinLclMeshRow = Y.getMap ()->getMinLocalIndex ();
    const LO myMaxLclMeshRow = Y.getMap ()->getMaxLocalIndex ();
    for (LO lclMeshRow = myMinLclMeshRow; lclMeshRow <= myMaxLclMeshRow; ++lclMeshRow) {
      typename BMV::little_vec_type Y_lcl = Y.getLocalBlock (lclMeshRow, 0);
      for (LO i = 0; i < blockSize; ++i) {
        TEST_EQUALITY( static_cast<Scalar> (Y_lcl(i)), requiredValue );
      }
    }

    TEST_NOTHROW( blockMat.setAllToScalar (STS::zero ()) );
    blockMat.applyBlock (X, Y, Teuchos::NO_TRANS, STS::one (), STS::zero ());
    for (LO lclMeshRow = myMinLclMeshRow; lclMeshRow <= myMaxLclMeshRow; ++lclMeshRow) {
      typename BMV::little_vec_type Y_lcl = Y.getLocalBlock (lclMeshRow, 0);
      for (LO i = 0; i < blockSize; ++i) {
        TEST_EQUALITY( static_cast<Scalar> (Y_lcl(i)), STS::zero () );
      }
    }

    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
  }


  // Test BlockCrsMatrix Import for the same graphs.  This is really
  // just a test of the "copy" part of copyAndPermute.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockCrsMatrix, ImportCopy, Scalar, LO, GO, Node )
  {
    typedef Tpetra::BlockMultiVector<Scalar, LO, GO, Node> BMV;
    typedef Tpetra::BlockCrsMatrix<Scalar, LO, GO, Node> BCM;
    typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::Import<LO, GO, Node> import_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;

    out << "Testing Tpetra::BlockCrsMatrix Import with same "
      "graphs" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const Comm<int> > comm = getDefaultComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    out << "Creating mesh row Map" << endl;

    const size_t numLocalMeshPoints = 12;
    const GO indexBase = 1;
    // mfh 16 May 2014: Tpetra::CrsGraph still needs the row Map as an
    // RCP.  Later interface changes will let us pass in the Map by
    // const reference and assume view semantics.
    RCP<const map_type> meshRowMapPtr =
      rcp (new map_type (INVALID, numLocalMeshPoints, indexBase, comm));
    const GO numGlobalMeshPoints = meshRowMapPtr->getGlobalNumElements ();
    const LO blockSize = 4;

    // Make a graph.  It happens to have two entries per row.
    out << "Creating mesh graph" << endl;
    graph_type graph (meshRowMapPtr, 2, Tpetra::StaticProfile);

    if (meshRowMapPtr->getNodeNumElements () > 0) {
      const GO myMinGblRow = meshRowMapPtr->getMinGlobalIndex ();
      const GO myMaxGblRow = meshRowMapPtr->getMaxGlobalIndex ();
      for (GO gblRow = myMinGblRow; gblRow <= myMaxGblRow; ++gblRow) {
        // Insert two entries, neither of which are on the diagonal.
        Teuchos::Array<GO> gblCols (2);
        gblCols[0] = indexBase + ((gblRow - indexBase) + 1) % numGlobalMeshPoints;
        gblCols[1] = indexBase + ((gblRow - indexBase) + 2) % numGlobalMeshPoints;
        graph.insertGlobalIndices (gblRow, gblCols ());
      }
    }
    graph.fillComplete ();

    // Create the two matrices.  They happen to have the same graph.
    out << "Create the matrices" << endl;
    BCM A1 (graph, blockSize);
    // We don't have to create the domain and range Maps all over
    // again.  Just reuse those of A1.
    BCM A2 (graph, * (A1.getDomainMap ()), * (A1.getRangeMap ()), blockSize);

    // Fill all entries of the first matrix with 3.
    const Scalar three = STS::one () + STS::one () + STS::one ();
    A1.setAllToScalar (three);

#ifdef HAVE_TPETRA_DEBUG
    if (! std::is_same<typename Kokkos::HostSpace, typename BCM::device_type::memory_space>::value) {
      // The above setAllToScalar should have run on device.
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
      TEST_ASSERT( A1.template need_sync<Kokkos::HostSpace> () );
      TEST_ASSERT( ! A1.template need_sync<typename BCM::device_type> () );
#else
      TEST_ASSERT( ! A1.need_sync_host () );
      TEST_ASSERT( ! A1.need_sync_device () );
#endif
    }
#endif // HAVE_TPETRA_DEBUG

    out << "The matrix A1, after construction:" << endl;
    A1.describe (out, Teuchos::VERB_EXTREME);

    // Fill all entries of the second matrix with -2.
    const Scalar minusTwo = -STS::one () - STS::one ();
    out << "Fill all entries of A2 with " << minusTwo << endl;
    A2.setAllToScalar (minusTwo);

    out << "The matrix A2, after construction:" << endl;
    A2.describe (out, Teuchos::VERB_EXTREME);

    out << "Create the Import" << endl;
    import_type imp (graph.getMap (), graph.getMap ());

    out << "Import A1 into A2" << endl;
    bool importSuccess = true;
    try {
      // The CombineMode doesn't matter for this example, since it
      // amounts to a matrix copy.  We use ADD arbitrarily.
      A2.doImport (A1, imp, Tpetra::ADD);
    } catch (std::exception& e) {
      importSuccess = false;
      if (myRank == 0) {
        out << "Import FAILED by throwing an exception" << endl;
      }
      for (int p = 0; p < numProcs; ++p) {
        if (p == myRank) {
          std::ostringstream os;
          os << "Process " << myRank << ": error messages from A1: "
             << A1.errorMessages () << endl
             << "Process " << myRank << ": error messages from A2: "
             << A2.errorMessages () << endl;
          std::cerr << os.str ();
        }
        comm->barrier (); // give time for output to complete
        comm->barrier ();
        comm->barrier ();
      }
    }
    if (A1.localError () || A2.localError ()) {
      if (myRank == 0) {
        out << "Import FAILED by reporting local error" << endl;
      }
      importSuccess = false;
    }

    TEST_ASSERT( importSuccess );

    int lclImportSuccess = importSuccess ? 1 : 0;
    int gblImportSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclImportSuccess, outArg (gblImportSuccess));
    importSuccess = (gblImportSuccess == 1);

    if (! importSuccess) {
      for (int p = 0; p < numProcs; ++p) {
        if (p == myRank) {
          std::ostringstream os;
          os << "Process " << myRank << ": error messages from A1: "
             << A1.errorMessages () << endl
             << "Process " << myRank << ": error messages from A2: "
             << A2.errorMessages () << endl;
          std::cerr << os.str ();
        }
        comm->barrier (); // give time for output to complete
        comm->barrier ();
        comm->barrier ();
      }
    }
    else { // doImport claims that it succeeded
      out << "Import claims that it succeeded; test the matrix" << endl;

      // Y := A2*X, where X is a block multivector (with one column)
      // full of 1s.  Since there are two block entries per row, each
      // of which is all 3s, we know that each entry of the result Y
      // will be 3*2*blockSize = 6*blockSize.
      const Scalar requiredValue = static_cast<Scalar> (6 * blockSize);
      BMV X (* (graph.getDomainMap ()), * (A2.getDomainMap ()), blockSize, static_cast<LO> (1));
      X.putScalar (STS::one ());
      BMV Y (* (graph.getRangeMap ()), * (A2.getRangeMap ()), blockSize, static_cast<LO> (1));
      A2.applyBlock (X, Y, Teuchos::NO_TRANS, STS::one (), STS::zero ());

      const LO myMinLclMeshRow = Y.getMap ()->getMinLocalIndex ();
      const LO myMaxLclMeshRow = Y.getMap ()->getMaxLocalIndex ();
      bool valsMatch = true;
      for (LO lclMeshRow = myMinLclMeshRow; lclMeshRow <= myMaxLclMeshRow; ++lclMeshRow) {
        typename BMV::little_vec_type Y_lcl = Y.getLocalBlock (lclMeshRow, 0);
        for (LO i = 0; i < blockSize; ++i) {
          if (static_cast<Scalar> (Y_lcl(i)) != requiredValue) {
            valsMatch = false;
          }
        }
      }
      TEST_ASSERT( valsMatch );
    }

    out << "The matrix A2, after Import (should be same as A1):" << endl;
    A2.describe (out, Teuchos::VERB_EXTREME);

    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
  }

  // Test BlockCrsMatrix Export for different graphs with different
  // row Maps.  This tests packAndPrepare and unpackAndCombine.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockCrsMatrix, ExportDiffRowMaps, Scalar, LO, GO, Node )
  {
    // typedef Tpetra::BlockMultiVector<Scalar, LO, GO, Node> BMV;
    typedef Tpetra::BlockCrsMatrix<Scalar, LO, GO, Node> BCM;
    typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::Export<LO, GO, Node> export_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    int lclSuccess = 1;
    int gblSuccess = 1;

    out << "Testing Tpetra::BlockCrsMatrix Import "
      "with different graphs with different row Maps" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const Comm<int> > comm = getDefaultComm ();
    const int myRank = comm->getRank ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    out << "Create nonoverlapping mesh row Map" << endl;
    const size_t numLocalMeshPoints = 12;
    const GO indexBase = 1;
    // mfh 16 May 2014: Tpetra::CrsGraph still needs the row Map as an
    // RCP.  Later interface changes will let us pass in the Map by
    // const reference and assume view semantics.
    RCP<const map_type> meshRowMapPtr =
      rcp (new map_type (INVALID, numLocalMeshPoints, indexBase, comm));
    const GO numGlobalMeshPoints = meshRowMapPtr->getGlobalNumElements ();
    const LO blockSize = 4;

    RCP<const map_type> meshDomainMapPtr = meshRowMapPtr;
    RCP<const map_type> meshRangeMapPtr = meshRowMapPtr;

    out << "Create overlapping mesh row Map" << endl;
    const size_t numRemoteMeshPoints = 2;
    // Put the "remote" indices at the beginning.  This makes the test
    // exercise the "permute" part of copyAndPermute.  We've already
    // exercised the "copy" part in the previous test.  Also, make
    // sure that the remote indices aren't necessarily consecutive or
    // in order.
    Array<GO> overlapMeshRowIndices (numLocalMeshPoints + numRemoteMeshPoints);
    for (size_t k = 0; k < numRemoteMeshPoints; ++k) {
      overlapMeshRowIndices[k] =
        (meshRowMapPtr->getMinGlobalIndex () - static_cast<GO> (2*k)) %
        numGlobalMeshPoints;
    }
    const size_t numOverlapMeshPoints = numLocalMeshPoints + numRemoteMeshPoints;
    RCP<const map_type> overlapMeshRowMapPtr =
      rcp (new map_type (INVALID, numOverlapMeshPoints, indexBase, comm));

    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      out << "*** FAILED to create maps! ***" << endl;
      return;
    }

    out << "Create graph with overlapping mesh row Map" << endl;
    // Make a graph.  It happens to have two entries per row.
    graph_type overlapGraph (overlapMeshRowMapPtr, 2, Tpetra::StaticProfile);
    if (overlapMeshRowMapPtr->getNodeNumElements () > 0) {
      const GO myMinGblRow = overlapMeshRowMapPtr->getMinGlobalIndex ();
      const GO myMaxGblRow = overlapMeshRowMapPtr->getMaxGlobalIndex ();
      for (GO gblRow = myMinGblRow; gblRow <= myMaxGblRow; ++gblRow) {
        // Insert two entries, neither of which are on the diagonal.
        Teuchos::Array<GO> gblCols (2);
        gblCols[0] = indexBase + ((gblRow - indexBase) + 1) % numGlobalMeshPoints;
        gblCols[1] = indexBase + ((gblRow - indexBase) + 2) % numGlobalMeshPoints;
        overlapGraph.insertGlobalIndices (gblRow, gblCols ());
      }
    }
    overlapGraph.fillComplete (meshDomainMapPtr, meshRangeMapPtr);

    {
      const LO lclNumRowsOverlap = static_cast<LO> (overlapMeshRowMapPtr->getNodeNumElements ());
      for (LO lclRow = 0; lclRow < lclNumRowsOverlap; ++lclRow) {
        const LO numEnt = static_cast<LO> (overlapGraph.getNumEntriesInLocalRow (lclRow));
        TEST_EQUALITY( numEnt, static_cast<LO> (2) );
      }
      const LO maxNumRowEnt = static_cast<LO> (overlapGraph.getNodeMaxNumRowEntries ());
      TEST_EQUALITY( maxNumRowEnt, static_cast<LO> (2) );
    }

    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      out << "*** FAILED to create or fill-complete graph with overlapping mesh row Map! ***" << endl;
      return;
    }

    RCP<export_type> theExport;
    RCP<graph_type> graph;
    try {
      out << "Create empty graph with nonoverlapping mesh row Map" << endl;
      graph = rcp (new graph_type (meshRowMapPtr, 0, Tpetra::DynamicProfile));

      out << "Create Export from overlap to nonoverlap mesh row Map" << endl;
      theExport = rcp (new export_type (overlapMeshRowMapPtr, meshRowMapPtr));

      out << "Import overlap graph into nonoverlap graph" << endl;
      graph->doExport (overlapGraph, *theExport, Tpetra::INSERT);

      out << "Call fillComplete on nonoverlap graph" << endl;
      graph->fillComplete (meshDomainMapPtr, meshRangeMapPtr);
    }
    catch (std::exception&) {
      lclSuccess = 0;
    }

    // Export can only add entries to each row, not remove them.
    {
      const LO lclNumRowsNonoverlap =
        static_cast<LO> (meshRowMapPtr->getNodeNumElements ());
      for (LO lclRow = 0; lclRow < lclNumRowsNonoverlap; ++lclRow) {
        const LO numEnt = static_cast<LO> (graph->getNumEntriesInLocalRow (lclRow));
        TEST_ASSERT( numEnt >= static_cast<LO> (2) );
      }
      const LO maxNumRowEnt = static_cast<LO> (graph->getNodeMaxNumRowEntries ());
      TEST_ASSERT( maxNumRowEnt >= static_cast<LO> (2) );
    }

    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      out << "*** CrsGraph Export or fillComplete failed! ***" << endl;
      return;
    }

    // Create the two matrices.  A_overlap has the overlapping mesh
    // row Map, and A has the nonoverlapping mesh row Map.
    RCP<BCM> A_overlap, A;
    try {
      out << "Create matrix with overlapping mesh row Map" << endl;
      A_overlap = rcp (new BCM (overlapGraph, blockSize));

      out << "Create matrix with nonoverlapping mesh row Map" << endl;
      A = rcp (new BCM (*graph, blockSize));
    }
    catch (std::exception&) {
      lclSuccess = 0;
    }
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      out << "*** FAILED to create the two BlockCrsMatrix instances! ***" << endl;
      return;
    }

    // Fill all entries of the first matrix with 3.
    try {
      const Scalar three = STS::one () + STS::one () + STS::one ();
      A_overlap->setAllToScalar (three);
    }
    catch (std::exception&) {
      lclSuccess = 0;
    }
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      out << "*** A_overlap->setAllToScalar(3) FAILED! ***" << endl;
      return;
    }

    // out << "The matrix A_overlap, after construction:" << endl;
    // A_overlap->describe (out, Teuchos::VERB_EXTREME);

    // Fill all entries of the second matrix with -2.
    try {
      const Scalar minusTwo = -STS::one () - STS::one ();
      A->setAllToScalar (minusTwo);
    } catch (std::exception&) {
      lclSuccess = 0;
    }
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      out << "*** A->setAllToScalar(-1) FAILED! ***" << endl;
      return;
    }

    // out << "The matrix A, after construction:" << endl;
    // A->describe (out, Teuchos::VERB_EXTREME);

    out << "Export A_overlap into A" << endl;
    try {
      A->doExport (*A_overlap, *theExport, Tpetra::REPLACE);
    }
    catch (std::exception& e) {
      lclSuccess = 0;
      std::ostringstream os;
      os << "Proc " << myRank << ": A->doExport(...) threw an exception: "
         << e.what () << endl;
      std::cerr << os.str ();
    }

    {
      const int lclErr = A->localError () ? 1 : 0;
      int gblErr = 0;
      reduceAll<int, int> (*comm, Teuchos::REDUCE_MAX, lclErr, outArg (gblErr));
      TEST_EQUALITY( gblErr, 0 );
      if (gblErr != 0) {
        out << "A reports a local error on some process!" << endl;
        gathervPrint (out, A->errorMessages (), *comm);
      }
    }

    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      using ::Tpetra::Details::gathervPrint;

      out << "*** Export FAILED!" << endl
          << "Error messages from A_overlap:" << endl;
      gathervPrint (out, A_overlap->errorMessages (), *comm);
      out << "Error messages from A:" << endl;
      gathervPrint (out, A->errorMessages (), *comm);
    }
    else { // doExport claims that it succeeded
      out << "Export claims that it succeeded" << endl;

      // // Y := A2*X, where X is a block multivector (with one column)
      // // full of 1s.  Since there are two block entries per row, each
      // // of which is all 3s, we know that each entry of the result Y
      // // will be 3*2*blockSize = 6*blockSize.
      // const Scalar requiredValue = static_cast<Scalar> (6 * blockSize);
      // BMV X (* (graph.getDomainMap ()), * (A2.getDomainMap ()), blockSize, static_cast<LO> (1));
      // X.putScalar (STS::one ());
      // BMV Y (* (graph.getRangeMap ()), * (A2.getRangeMap ()), blockSize, static_cast<LO> (1));
      // A2.applyBlock (X, Y, Teuchos::NO_TRANS, STS::one (), STS::zero ());

      // const LO myMinLclMeshRow = Y.getMap ()->getMinLocalIndex ();
      // const LO myMaxLclMeshRow = Y.getMap ()->getMaxLocalIndex ();
      // bool valsMatch = true;
      // for (LO lclMeshRow = myMinLclMeshRow; lclMeshRow <= myMaxLclMeshRow; ++lclMeshRow) {
      //   typename BMV::little_vec_type Y_lcl = Y.getLocalBlock (lclMeshRow, 0);
      //   for (LO i = 0; i < blockSize; ++i) {
      //     if (Y_lcl(i) != requiredValue) {
      //       valsMatch = false;
      //     }
      //   }
      // }
      // TEST_ASSERT( valsMatch );
    }

    // out << "The matrix A2, after Export (should be same as A1):" << endl;
    // A2.describe (out, Teuchos::VERB_EXTREME);

    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
  }

  //
  // Test BlockCrsMatrix's localGaussSeidel with a block diagonal matrix.
  //
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockCrsMatrix, localGSDiagonalMatrix, Scalar, LO, GO, Node )
  {
    using Kokkos::ALL;
    typedef Scalar ST;
    typedef Tpetra::BlockVector<ST, LO, GO, Node> BV;
    typedef Tpetra::BlockCrsMatrix<ST, LO, GO, Node> BCM;
    typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef typename graph_type::device_type device_type;
    typedef typename BCM::impl_scalar_type IST;
    typedef Teuchos::ScalarTraits<ST> STS;

    const ST two = STS::one () + STS::one ();
    const ST three = two + STS::one ();
    const auto tol = 10.0 * STS::eps ();

    out << "Testing Tpetra::BlockCrsMatrix::localGaussSeidel "
      "with a matrix whose graph is diagonal" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const Comm<int> > comm = getDefaultComm ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    out << "Creating mesh row Map" << endl;

    const size_t numLocalMeshPoints = 12;
    const GO indexBase = 1;
    // Use a block size that is not a power of 2, to test correctness
    // in case the matrix pads blocks for SIMD-ization.
    const LO blockSize = 3;

    // mfh 20 May 2014: Tpetra::CrsGraph still needs the row Map as an
    // RCP.  Later interface changes will let us pass in the Map by
    // const reference and assume view semantics.
    RCP<const map_type> meshRowMapPtr =
      rcp (new map_type (INVALID, numLocalMeshPoints, indexBase, comm));
    const map_type& meshRowMap = *meshRowMapPtr;

    out << "Creating mesh graph" << endl;

    const size_t maxNumEntPerRow = 1;
    graph_type graph (meshRowMapPtr, maxNumEntPerRow, Tpetra::StaticProfile);

    // Fill the graph with only diagonal entries
    Teuchos::Array<GO> gblColInds (maxNumEntPerRow);
    for (LO lclRowInd = meshRowMap.getMinLocalIndex ();
         lclRowInd <= meshRowMap.getMaxLocalIndex (); ++lclRowInd) {
      const GO gblRowInd = meshRowMap.getGlobalElement (lclRowInd);
      gblColInds[0] = gblRowInd;
      graph.insertGlobalIndices (gblRowInd, gblColInds ());
    }
    graph.fillComplete ();

    // Get the graph's column Map (the "mesh column Map").
    TEST_ASSERT( ! graph.getColMap ().is_null () );
    map_type meshColMap = * (graph.getColMap ());

    out << "Creating BlockCrsMatrix" << endl;

    // Construct the BlockCrsMatrix.
    BCM blockMat;
    TEST_NOTHROW( blockMat = BCM (graph, blockSize) );
    BV residual;
    TEST_NOTHROW( residual = BV(meshRowMap, blockSize));
    BV solution;
    TEST_NOTHROW( solution = BV(meshRowMap, blockSize));

    Teuchos::Array<ST> basematrix (blockSize*blockSize, STS::zero ());
    basematrix[0] = two;
    basematrix[2] = three;
    basematrix[3] = three;
    basematrix[4] = two;
    basematrix[7] = three;
    basematrix[8] = two;

    Teuchos::Array<ST> baseResidual(blockSize, STS::zero ());
    baseResidual[0] = STS::one ();
    baseResidual[1] = three;
    baseResidual[2] = -two;

    // FIXME (mfh 25 Aug 2014) This will likely only work with Scalar
    // = float or double.  On the other hand, the author of these
    // tests understood that and only instantiated them for
    // Scalar=double (see the instantiations list below).
    Teuchos::Array<ST> exactSolution(blockSize, STS::zero());
    exactSolution[0] = 43.0/35.0;
    exactSolution[1] = -12.0/35.0;
    exactSolution[2] = -17.0/35.0;

    // NOTE (mfh 26 May 2016) We may start modifying the matrix on
    // host now, because we haven't yet done anything to it on device.

    Teuchos::Array<LO> lclColInds(1);
    for (LO lclRowInd = meshRowMap.getMinLocalIndex ();
         lclRowInd <= meshRowMap.getMaxLocalIndex (); ++lclRowInd) {
      for (LO k = 0; k < blockSize*blockSize; ++k) {
        basematrix[k] *= two;
      }
      for (LO k = 0; k < blockSize; ++k) {
        baseResidual[k] *= two;
      }
      lclColInds[0] = lclRowInd;
      blockMat.replaceLocalValues (lclRowInd, lclColInds.getRawPtr (),
                                   basematrix.getRawPtr (), 1);
      residual.replaceLocalValues (lclRowInd, baseResidual.getRawPtr ());
      solution.replaceLocalValues (lclRowInd, baseResidual.getRawPtr ());
    }

    Kokkos::View<size_t*, device_type> diagonalOffsets ("offsets", numLocalMeshPoints);
    graph.getLocalDiagOffsets (diagonalOffsets);

    // Sync the matrix to device, since getLocalDiagCopy runs there.
    blockMat.template sync<device_type> ();

    typedef Kokkos::View<IST***, device_type> block_diag_type;
    block_diag_type blockDiag ("blockDiag", numLocalMeshPoints,
                               blockSize, blockSize);
    blockMat.getLocalDiagCopy (blockDiag, diagonalOffsets);

    Kokkos::View<int**, device_type> pivots ("pivots", numLocalMeshPoints, blockSize);
    // That's how we found this test: the pivots array was filled with ones.
    Kokkos::deep_copy (pivots, 1);
    out << "pivots size = " << pivots.extent(0) << endl;

    for (LO lclMeshRow = 0; lclMeshRow < static_cast<LO> (numLocalMeshPoints); ++lclMeshRow) {
      auto diagBlock = Kokkos::subview (blockDiag, lclMeshRow, ALL (), ALL ());

      // Make sure that the diagonal block is correct.
      Scalar* blkVals = NULL;
      const LO* blkColInds = NULL;
      LO blkNumEnt = 0;
      blockMat.getLocalRowView (lclMeshRow, blkColInds, blkVals, blkNumEnt);

      TEST_EQUALITY( blkNumEnt, static_cast<LO> (1) );
      if (blkNumEnt == 1) {
        typename BCM::const_little_block_type diagBlock2 ((typename BCM::const_little_block_type::value_type*) blkVals, blockSize, blockSize);
        for (LO j = 0; j < blockSize; ++j) {
          for (LO i = 0; i < blockSize; ++i) {
            TEST_EQUALITY( diagBlock(i,j), diagBlock2(i,j) );
          }
        }
      }

      auto ipiv = Kokkos::subview (pivots, lclMeshRow, ALL ());
      int info = 0;
      Tpetra::Experimental::GETF2 (diagBlock, ipiv, info);
      TEST_EQUALITY( info, 0 );

      // GETRI needs workspace.  Use host space for now.
      Teuchos::Array<IST> workVec (blockSize);
      Kokkos::View<IST*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
        work (workVec.getRawPtr (), blockSize);
      Tpetra::Experimental::GETRI (diagBlock, ipiv, work, info);
      TEST_EQUALITY( info, 0 );
    }

    blockMat.localGaussSeidel (residual, solution, blockDiag,
                               STS::one(), Tpetra::Forward);

    for (LO lclRowInd = meshRowMap.getMinLocalIndex ();
         lclRowInd <= meshRowMap.getMaxLocalIndex (); ++lclRowInd) {
      typename BV::little_vec_type xlcl = solution.getLocalBlock (lclRowInd);
      ST* x = reinterpret_cast<ST*> (xlcl.data ());
      out << "row = " << lclRowInd << endl;
      for (LO k = 0; k < blockSize; ++k) {
        TEST_FLOATING_EQUALITY( x[k], exactSolution[k], tol );
        x[k] = -STS::one ();
      }
    }

    blockMat.localGaussSeidel (residual, solution, blockDiag,
                               STS::one (), Tpetra::Backward);
    for (LO lclRowInd = meshRowMap.getMinLocalIndex ();
         lclRowInd <= meshRowMap.getMaxLocalIndex (); ++lclRowInd) {
      typename BV::little_vec_type xlcl = solution.getLocalBlock (lclRowInd);
      ST* x = reinterpret_cast<ST*> (xlcl.data ());
      for (LO k = 0; k < blockSize; ++k) {
        TEST_FLOATING_EQUALITY( x[k], exactSolution[k], tol );
      }
    }

    // Final output
    int lclSuccess = success;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
  }

  //
  // Test BlockCrsMatrix's localGaussSeidel with a triangular matrix (???)
  //
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockCrsMatrix, localGSTriangularMatrices, ST, LO, GO, Node )
  {
    using Kokkos::ALL;
    typedef Tpetra::BlockVector<ST, LO, GO, Node> BV;
    typedef Tpetra::BlockCrsMatrix<ST, LO, GO, Node> BCM;
    typedef typename BCM::impl_scalar_type IST;
    typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
    typedef typename graph_type::device_type device_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Teuchos::ScalarTraits<ST> STS;

    const ST two = STS::one () + STS::one ();
    const ST three = STS::one () + STS::one () + STS::one ();
    const auto tol = 10.0 * STS::eps ();

    out << "Testing Tpetra::BlockCrsMatrix::localGaussSeidel "
      "with a triangular matrix ( ??? )" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const Comm<int> > comm = getDefaultComm ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    out << "Creating mesh row Map" << endl;

    const size_t numLocalMeshPoints = 3;
    const GO indexBase = 1;
    // Use a block size that is not a power of 2, to test correctness
    // in case the matrix pads blocks for SIMD-ization.
    const LO blockSize = 2;

    // mfh 20 May 2014: Tpetra::CrsGraph still needs the row Map as an
    // RCP.  Later interface changes will let us pass in the Map by
    // const reference and assume view semantics.
    RCP<const map_type> meshRowMapPtr =
      rcp (new map_type (INVALID, numLocalMeshPoints, indexBase, comm));
    const map_type& meshRowMap = *meshRowMapPtr;

    out << "Creating mesh graph" << endl;

    const size_t maxNumEntPerRow = 3;
    graph_type graph (meshRowMapPtr, maxNumEntPerRow, Tpetra::StaticProfile);

    Teuchos::Array<GO> gblColInds (maxNumEntPerRow);

    for (LO lclRowInd = meshRowMap.getMinLocalIndex ();
         lclRowInd <= meshRowMap.getMaxLocalIndex (); ++lclRowInd) {

      const GO gblRowInd = meshRowMap.getGlobalElement (lclRowInd);

      for (size_t k = 0; k < maxNumEntPerRow; ++k) {
        const LO lclColInd = meshRowMap.getMinLocalIndex () + k;
        const GO gblColInd = meshRowMap.getGlobalElement (lclColInd);
        gblColInds[k] = gblColInd;
      }
      graph.insertGlobalIndices (gblRowInd, gblColInds ());
    }

    graph.fillComplete ();

    // Get the graph's column Map (the "mesh column Map").
    TEST_ASSERT( ! graph.getColMap ().is_null () );
    map_type meshColMap = * (graph.getColMap ());

    out << "Creating BlockCrsMatrix" << endl;

    // Construct the BlockCrsMatrix.
    BCM blockMat;
    TEST_NOTHROW( blockMat = BCM (graph, blockSize) );
    BV residual;
    TEST_NOTHROW( residual = BV (meshRowMap, blockSize));
    BV solution;
    TEST_NOTHROW( solution = BV (meshRowMap, blockSize));

    Teuchos::Array<ST> basematrix (maxNumEntPerRow * maxNumEntPerRow, STS::zero ());
    basematrix[0] = two;
    basematrix[3] = three;
    basematrix[4] = two;
    basematrix[6] = STS::one ();
    basematrix[7] = three;
    basematrix[8] = two;

    Teuchos::Array<ST> baseResidual (maxNumEntPerRow, STS::zero ());
    baseResidual[0] = STS::one ();
    baseResidual[1] = three;
    baseResidual[2] = -two;

    // FIXME (mfh 25 Aug 2014) This will likely only work with Scalar
    // = float or double.  On the other hand, the author of these
    // tests understood that and only instantiated them for
    // Scalar=double (see the instantiations list below).
    Teuchos::Array<ST> exactSolution (maxNumEntPerRow, STS::zero ());
    exactSolution[0] = 0.5;
    exactSolution[1] = 0.75;
    exactSolution[2] = -19.0/8.0;

    Teuchos::Array<ST> assembleMatrix (blockSize * blockSize, STS::zero ());
    Teuchos::Array<ST> assembleResidual (blockSize, STS::zero ());

    Teuchos::Array<LO> lclColInds (1);

    // lower triangular matrix
    for (LO lclRowInd = meshRowMap.getMinLocalIndex ();
         lclRowInd <= meshRowMap.getMaxLocalIndex (); ++lclRowInd) {
      const LO rowOffset = lclRowInd - meshRowMap.getMinLocalIndex ();
      for (LO k = 0; k < blockSize; ++k) {
        assembleResidual[k] = baseResidual[rowOffset];
      }
      residual.replaceLocalValues (lclRowInd, &assembleResidual[0]);

      for (LO k = 0; k < blockSize; ++k) {
        assembleResidual[k] = STS::zero ();
      }
      solution.replaceLocalValues (lclRowInd, &assembleResidual[0]);

      for (size_t i = 0; i < maxNumEntPerRow; ++i) {
        for (LO j = 0; j < blockSize*blockSize; ++j) {
          assembleMatrix[j] = STS::zero ();
        }
        const size_t indexBaseMatrix = (maxNumEntPerRow)*rowOffset+i;
        for (LO j = 0; j < blockSize; ++j) {
          assembleMatrix[(blockSize+1)*j] = basematrix[indexBaseMatrix];
        }

        lclColInds[0] = meshRowMap.getMinLocalIndex () + i;
        blockMat.replaceLocalValues(lclRowInd, lclColInds.getRawPtr(), &assembleMatrix[0], 1);
      }
    }

    Kokkos::View<size_t*, device_type> diagonalOffsets ("offsets", numLocalMeshPoints);
    graph.getLocalDiagOffsets(diagonalOffsets);

    typedef Kokkos::View<IST***, device_type> block_diag_type;
    block_diag_type blockDiag ("blockDiag", numLocalMeshPoints,
                               blockSize, blockSize);
    blockMat.getLocalDiagCopy (blockDiag, diagonalOffsets);

    Kokkos::View<int**, device_type> pivots ("pivots", numLocalMeshPoints, blockSize);
    // That's how we found this test: the pivots array was filled with ones.
    Kokkos::deep_copy (pivots, 1);

    for (LO lclMeshRow = 0; lclMeshRow < static_cast<LO> (numLocalMeshPoints); ++lclMeshRow) {
      auto diagBlock = Kokkos::subview (blockDiag, lclMeshRow, ALL (), ALL ());
      auto ipiv = Kokkos::subview (pivots, lclMeshRow, ALL ());
      int info = 0;
      Tpetra::Experimental::GETF2 (diagBlock, ipiv, info);
      TEST_EQUALITY( info, 0 );

      // GETRI needs workspace.  Use host space for now.
      Teuchos::Array<IST> workVec (blockSize);
      Kokkos::View<IST*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
        work (workVec.getRawPtr (), blockSize);
      Tpetra::Experimental::GETRI (diagBlock, ipiv, work, info);
      TEST_EQUALITY( info, 0 );
    }

    blockMat.localGaussSeidel (residual, solution, blockDiag,
                               STS::one (), Tpetra::Forward);

    for (LO lclRowInd = meshRowMap.getMinLocalIndex ();
         lclRowInd <= meshRowMap.getMaxLocalIndex(); ++lclRowInd) {
      const LO rowOffset = lclRowInd - meshRowMap.getMinLocalIndex ();
      typename BV::little_vec_type xlcl = solution.getLocalBlock (lclRowInd);
      ST* x = reinterpret_cast<ST*> (xlcl.data ());
      for (LO k = 0; k < blockSize; ++k) {
        TEST_FLOATING_EQUALITY( x[k], exactSolution[rowOffset], tol );
        x[k] = -STS::one ();
      }
    }

    // upper triangular matrix
    //
    // FIXME (mfh 25 Aug 2014) This will likely only work with Scalar
    // = float or double.  On the other hand, the author of these
    // tests understood that and only instantiated them for
    // Scalar=double (see the instantiations list below).
    exactSolution[0] = -3.5;
    exactSolution[1] = 3.0;
    exactSolution[2] = -1.0;

    for (LO lclRowInd = meshRowMap.getMinLocalIndex ();
         lclRowInd <= meshRowMap.getMaxLocalIndex (); ++lclRowInd) {

      const LO rowOffset = lclRowInd - meshRowMap.getMinLocalIndex ();
      for (size_t i = 0; i < maxNumEntPerRow; ++i) {
        for (LO k = 0; k < blockSize; ++k) {
          assembleResidual[k] = STS::zero ();
        }
        solution.replaceLocalValues(lclRowInd, &assembleResidual[0]);

        for (LO j = 0; j < blockSize*blockSize; ++j) {
          assembleMatrix[j] = STS::zero ();
        }
        const size_t indexBaseMatrix = maxNumEntPerRow * i + rowOffset;
        for (LO j = 0; j < blockSize; ++j) {
          assembleMatrix[(blockSize+1)*j] = basematrix[indexBaseMatrix];
        }

        lclColInds[0] = meshRowMap.getMinLocalIndex () + i;
        blockMat.replaceLocalValues (lclRowInd, lclColInds.getRawPtr (),
                                     &assembleMatrix[0], 1);
      }
    }

    blockMat.getLocalDiagCopy (blockDiag, diagonalOffsets);

    for (LO lclMeshRow = 0; lclMeshRow < static_cast<LO> (numLocalMeshPoints); ++lclMeshRow) {
      auto diagBlock = Kokkos::subview (blockDiag, lclMeshRow, ALL (), ALL ());
      auto ipiv = Kokkos::subview (pivots, lclMeshRow, ALL ());
      int info = 0;
      Tpetra::Experimental::GETF2 (diagBlock, ipiv, info);
      TEST_EQUALITY( info, 0 );

      // GETRI needs workspace.  Use host space for now.
      Teuchos::Array<IST> workVec (blockSize);
      Kokkos::View<IST*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
        work (workVec.getRawPtr (), blockSize);
      Tpetra::Experimental::GETRI (diagBlock, ipiv, work, info);
      TEST_EQUALITY( info, 0 );
    }

    blockMat.localGaussSeidel (residual, solution, blockDiag,
                               STS::one (), Tpetra::Symmetric);

    for (LO lclRowInd = meshRowMap.getMinLocalIndex ();
         lclRowInd <= meshRowMap.getMaxLocalIndex(); ++lclRowInd) {
      const LO rowOffset = lclRowInd - meshRowMap.getMinLocalIndex ();
      typename BV::little_vec_type xlcl = solution.getLocalBlock (lclRowInd);
      ST* x = reinterpret_cast<ST*> (xlcl.data ());
      for (LO k = 0; k < blockSize; ++k) {
        TEST_FLOATING_EQUALITY( x[k], exactSolution[rowOffset], tol );
        x[k] = -STS::one ();
      }
    }

    // Final output
    int lclSuccess = success;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
  }

  //
  // Test conversion from CrsMatrix to BlockCrsMatrix.  This test is valid only on
  // 1, 2, or 4 processes.  (See note below.)
  //
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockCrsMatrix, point2block, ST, LO, GO, Node )
  {
    typedef Tpetra::BlockCrsMatrix<ST, LO, GO, Node> block_matrix_type;
    typedef Tpetra::CrsMatrix<ST, LO, GO, Node>                    crs_matrix_type;
    typedef Tpetra::Import<LO, GO, Node>                           import_type;
    typedef Tpetra::MultiVector<ST, LO, GO, Node>                  mv_type;
    typedef Tpetra::Map<LO, GO, Node>                              map_type;
    typedef Tpetra::MatrixMarket::Reader<crs_matrix_type>          reader_type;
    typedef Teuchos::ScalarTraits<ST>                              STS;
    typedef typename STS::magnitudeType                            magnitude_type;
    ST zero = STS::zero(), one = STS::one();

    constexpr bool printToCerr = true;

    Teuchos::OSTab tab0 (out);
    if (printToCerr) {
      std::cerr << "Test conversion from (point) CrsMatrix to BlockCrsMatrix" << endl;
    }
    else {
      out << "Test conversion from (point) CrsMatrix to BlockCrsMatrix" << endl;
    }
    Teuchos::OSTab tab1 (out);

    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    std::ostringstream errStrm;

    RCP<const Comm<int> > comm = getDefaultComm();
    std::string matrixFile;
    if (STS::isComplex) {
      matrixFile = "blockA-complex.mm";
    }
    else {
      matrixFile = "blockA.mm";
    }
    if (printToCerr) {
      std::cerr << "Read CrsMatrix from file \"" << matrixFile << "\"" << endl;
    }
    else {
      out << "Read CrsMatrix from file \"" << matrixFile << "\"" << endl;
    }
    RCP<crs_matrix_type> pointMatrix;
    try {
      pointMatrix = reader_type::readSparseFile(matrixFile, comm);
    }
    catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Proc " << comm->getRank () << ": readSparseFile threw an "
        "std::exception: " << e.what ();
    }
    catch (...) {
      lclSuccess = 0;
      errStrm << "Proc " << comm->getRank () << ": readSparseFile threw an "
        "exception not a subclass of std::exception";
    }
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
    if (gblSuccess != 1) {
      if (printToCerr) {
        gathervPrint (std::cerr, errStrm.str (), *comm);
      }
      else {
        gathervPrint (out, errStrm.str (), *comm);
      }
      success = false;
      return;
    }

    if (printToCerr) {
      std::cerr << "Migrate input CrsMatrix to final parallel distribution" << endl;
    }
    else {
      out << "Migrate input CrsMatrix to final parallel distribution" << endl;
    }

    // Migrate pointMatrix to final parallel distribution.
    // Note that the input matrix has 12 point rows, with block size 3.  Point rows associated with a mesh node
    // must stay together.  This means the serial matrix can only be migrated to 1,2 or 4 processes.  3 processes
    // would split up dofs associate with a mesh node.
    Tpetra::global_size_t numGlobElts = pointMatrix->getRowMap()->getGlobalNumElements();
    const GO indexBase = 0;
    RCP<const map_type> parPointMap =
      rcp (new map_type (numGlobElts, indexBase, comm));
    RCP<crs_matrix_type> parPointMatrix =
      rcp (new crs_matrix_type (parPointMap, pointMatrix->getGlobalMaxNumRowEntries ()));
    RCP<const import_type> importer =
      rcp (new import_type (pointMatrix->getRowMap(), parPointMap));

    try {
      parPointMatrix->doImport(*pointMatrix, *importer, Tpetra::INSERT);
    }
    catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Proc " << comm->getRank () << ": doImport (point matrix to "
        "point matrix) threw an std::exception: " << e.what ();
    }
    catch (...) {
      lclSuccess = 0;
      errStrm << "Proc " << comm->getRank () << ": doImport (point matrix to "
        "point matrix) threw an exception not a subclass of std::exception";
    }
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
    if (gblSuccess != 1) {
      if (printToCerr) {
        gathervPrint (std::cerr, errStrm.str (), *comm);
      }
      else {
        gathervPrint (out, errStrm.str (), *comm);
      }
      success = false;
      return;
    }

    parPointMatrix->fillComplete();
    pointMatrix.swap(parPointMatrix);

    if (printToCerr) {
      std::cerr << "Convert CrsMatrix to BlockCrsMatrix" << endl;
    }
    else {
      out << "Convert CrsMatrix to BlockCrsMatrix" << endl;
    }

    int blockSize = 3;
    RCP<block_matrix_type> blockMatrix;
    try {
      blockMatrix = Tpetra::Experimental::convertToBlockCrsMatrix(*pointMatrix,blockSize);
    }
    catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Proc " << comm->getRank () << ": convertToBlockCrsMatrix "
        "threw an std::exception: " << e.what ();
    }
    catch (...) {
      lclSuccess = 0;
      errStrm << "Proc " << comm->getRank () << ": convertToBlockCrsMatrix "
        "threw an exception not a subclass of std::exception";
    }
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
    if (gblSuccess != 1) {
      gathervPrint (out, errStrm.str (), *comm);
      success = false;
      return;
    }

    if (printToCerr) {
      std::cerr << "Test resulting BlockCrsMatrix by comparing mat-vec result "
        "against CrsMatrix mat-vec result" << endl;
    }
    else {
      out << "Test resulting BlockCrsMatrix by comparing mat-vec result against "
        "CrsMatrix mat-vec result" << endl;
    }

    //normalized pseudo-random vector
    RCP<mv_type> randVec = rcp(new mv_type(pointMatrix->getDomainMap(),1));
    randVec->randomize();
    Teuchos::Array<magnitude_type> normVec1(1);
    randVec->norm2(normVec1);
    randVec->scale(1.0/normVec1[0]);

    RCP<mv_type> resultVec1 = rcp(new mv_type(pointMatrix->getRangeMap(),1));
    if (printToCerr) {
      std::cerr << "CrsMatrix::apply" << endl;
    }
    else {
      out << "CrsMatrix::apply" << endl;
    }
    try {
      pointMatrix->apply(*randVec, *resultVec1, Teuchos::NO_TRANS, one, zero);
    }
    catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Proc " << comm->getRank () << ": CrsMatrix::apply "
        "threw an std::exception: " << e.what ();
    }
    catch (...) {
      lclSuccess = 0;
      errStrm << "Proc " << comm->getRank () << ": CrsMatrix::apply "
        "threw an exception not a subclass of std::exception";
    }
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
    if (gblSuccess != 1) {
      if (printToCerr) {
        gathervPrint (std::cerr, errStrm.str (), *comm);
      }
      else {
        gathervPrint (out, errStrm.str (), *comm);
      }
      success = false;
      return;
    }
    if (printToCerr) {
      std::cerr << "Compute norm of result" << endl;
    }
    else {
      out << "Compute norm of result" << endl;
    }
    resultVec1->norm2(normVec1);

    RCP<mv_type> resultVec2 = rcp(new mv_type(blockMatrix->getRangeMap(),1));
    if (printToCerr) {
      std::cerr << "BlockCrsMatrix::apply" << endl;
    }
    else {
      out << "BlockCrsMatrix::apply" << endl;
    }
    try {
      blockMatrix->apply(*randVec, *resultVec2, Teuchos::NO_TRANS, one, zero);
    }
    catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Proc " << comm->getRank () << ": BlockCrsMatrix::apply "
        "threw an std::exception: " << e.what ();
    }
    catch (...) {
      lclSuccess = 0;
      errStrm << "Proc " << comm->getRank () << ": BlockCrsMatrix::apply "
        "threw an exception not a subclass of std::exception";
    }
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
    if (gblSuccess != 1) {
      if (printToCerr) {
        gathervPrint (std::cerr, errStrm.str (), *comm);
      }
      else {
        gathervPrint (out, errStrm.str (), *comm);
      }
      success = false;
      return;
    }
    Teuchos::Array<magnitude_type> normVec2(1);
    resultVec2->norm2(normVec2);

    resultVec2->update(-1.0,*resultVec1,1.0);
    Teuchos::Array<magnitude_type> normDelta(1);
    resultVec2->norm2(normDelta);
    Teuchos::Array<magnitude_type> relativeError(1);
    relativeError[0] = STS::magnitude(normDelta[0] / normVec1[0]);

    std::ostringstream normStr;
    normStr << "||CSR*xrand|| = " << normVec1[0] << ", ||CSR*xrand - BCSR*xrand|| / ||CSR*xrand|| = " << relativeError[0];
    if (printToCerr) {
      std::cerr << normStr.str() << std::endl;
    }
    else {
      out << normStr.str() << std::endl;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(relativeError[0]>1e-8, std::runtime_error, "BlockCrsMatrix matvec does not produce same result as CrsMatrix matvec.");
  }


//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockCrsMatrix, ctor, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockCrsMatrix, basic, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockCrsMatrix, write, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockCrsMatrix, getLocalDiagCopy, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockCrsMatrix, SetAllToScalar, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockCrsMatrix, ImportCopy, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockCrsMatrix, ExportDiffRowMaps, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockCrsMatrix, localGSDiagonalMatrix, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockCrsMatrix, localGSTriangularMatrices, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockCrsMatrix, point2block, SCALAR, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  // NOTE (mfh 24 Sep 2015) It only makes sense to test over Scalar
  // types which have a Teuchos::LAPACK implementation.

  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP )

} // namespace (anonymous)


