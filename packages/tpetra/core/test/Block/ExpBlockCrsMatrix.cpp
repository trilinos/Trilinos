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

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_TestingUtilities.hpp>
#include <Tpetra_Experimental_BlockCrsMatrix.hpp>
#include <Tpetra_Experimental_BlockVector.hpp>

namespace {
  using Tpetra::TestingUtilities::getNode;
  using Tpetra::TestingUtilities::getDefaultComm;
  using Teuchos::Array;
  using Teuchos::Comm;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::endl;
  typedef Tpetra::global_size_t GST;

  //
  // UNIT TESTS
  //

  // Test BlockCrsMatrix's constructors.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( ExpBlockCrsMatrix, ctor, Scalar, LO, GO, Node )
  {
    typedef Tpetra::Experimental::BlockMultiVector<Scalar, LO, GO, Node> BMV;
    typedef Tpetra::Experimental::BlockCrsMatrix<Scalar, LO, GO, Node> BCM;
    typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;

    out << "Testing Tpetra::Experimental::BlockCrsMatrix" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const Comm<int> > comm = getDefaultComm ();
    RCP<Node> node = getNode<Node> ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    out << "Creating mesh row Map" << endl;

    const size_t numLocalMeshPoints = 12;
    const GO indexBase = 1;
    // mfh 16 May 2014: Tpetra::CrsGraph still needs the row Map as an
    // RCP.  Later interface changes will let us pass in the Map by
    // const reference and assume view semantics.
    RCP<const map_type> meshRowMapPtr =
      rcp (new map_type (INVALID, numLocalMeshPoints, indexBase, comm, node));

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
  // BlockCrsMatrix with a nontrivial graph, exercise getLocalRowView
  // and replaceLocalValues, and test applyBlock and apply.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( ExpBlockCrsMatrix, basic, Scalar, LO, GO, Node )
  {
    typedef Tpetra::Experimental::BlockMultiVector<Scalar, LO, GO, Node> BMV;
    typedef Tpetra::Experimental::BlockVector<Scalar, LO, GO, Node> BV;
    typedef Tpetra::Experimental::BlockCrsMatrix<Scalar, LO, GO, Node> BCM;
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

    out << "Testing Tpetra::Experimental::BlockCrsMatrix basic "
      "functionality" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const Comm<int> > comm = getDefaultComm ();
    RCP<Node> node = getNode<Node> ();
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
      rcp (new map_type (INVALID, numLocalMeshPoints, indexBase, comm, node));
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
      std::copy (myLclColInds, myLclColInds + 2, myLclColIndsCopy.begin ());
      std::sort (myLclColIndsCopy.begin (), myLclColIndsCopy.end ());
      TEST_COMPARE_ARRAYS( lclColInds, myLclColIndsCopy );

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
        little_block_type tempBlock (tempBlockPtr, blockSize, blockSize, 1);
        for (LO j = 0; j < blockSize; ++j) {
          for (LO i = 0; i < blockSize; ++i) {
            TEST_ASSERT( tempBlock(i,j) == STS::zero () );
          }
        }
      } // for each entry in the row

      // Create a block pattern which verifies that the matrix stores
      // blocks in row-major order.
      for (LO k = 0; k < numEnt; ++k) {
        Scalar* const tempBlockPtr = tempBlockSpace.getRawPtr () +
          k * blockSize * blockSize;
        little_block_type tempBlock (tempBlockPtr, blockSize, blockSize, 1);
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
        little_block_type curBlk (curBlkPtr, blockSize, blockSize, 1);

        for (LO j = 0; j < blockSize; ++j) {
          for (LO i = 0; i < blockSize; ++i) {
            TEST_ASSERT( curBlk(i,j) == static_cast<Scalar> (static_cast<MT> (j + i * blockSize)) );
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
        TEST_ASSERT( X_lcl.getRawPtr () != NULL );
        TEST_ASSERT( X_lcl.getBlockSize () == blockSize );
        for (LO i = 0; i < blockSize; ++i) {
          X_lcl(i) = static_cast<Scalar> (static_cast<MT> (blockSize - i));
        }
      }

      TEST_NOTHROW( blockMat.applyBlock (X, Y) );

      const map_type& meshRangeMap = * (graph.getRangeMap ());
      for (LO lclRanIdx = meshRangeMap.getMinLocalIndex ();
           lclRanIdx <= meshRangeMap.getMaxLocalIndex (); ++lclRanIdx) {
        little_vec_type Y_lcl = Y.getLocalBlock (lclRanIdx);
        TEST_ASSERT( Y_lcl.getRawPtr () != NULL );
        TEST_ASSERT( Y_lcl.getBlockSize () == blockSize );

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
          TEST_ASSERT( Y_lcl(i) == static_cast<Scalar> (static_cast<MT> (expectedVal)) );
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
        TEST_ASSERT( Y_lcl.getRawPtr () != NULL );
        TEST_ASSERT( Y_lcl.getBlockSize () == blockSize );

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
          TEST_ASSERT( Y_lcl(i) == static_cast<Scalar> (static_cast<MT> (expectedVal)) );
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
          TEST_ASSERT( X_lcl.getRawPtr () != NULL );
          TEST_ASSERT( X_lcl.getBlockSize () == blockSize );
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
          TEST_ASSERT( Y_lcl.getRawPtr () != NULL );
          TEST_ASSERT( Y_lcl.getBlockSize () == blockSize );

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
            TEST_ASSERT( Y_lcl(i) == static_cast<Scalar> (static_cast<MT> (expectedVal)) );
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
          TEST_ASSERT( Y_lcl.getRawPtr () != NULL );
          TEST_ASSERT( Y_lcl.getBlockSize () == blockSize );

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
            TEST_ASSERT( Y_lcl(i) == static_cast<Scalar> (static_cast<MT> (expectedVal)) );
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
        TEST_ASSERT( X_lcl.getRawPtr () != NULL );
        TEST_ASSERT( X_lcl.getBlockSize () == blockSize );
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
        TEST_ASSERT( Y_lcl.getRawPtr () != NULL );
        TEST_ASSERT( Y_lcl.getBlockSize () == blockSize );

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
          TEST_ASSERT( Y_lcl(i) == static_cast<Scalar> (static_cast<MT> (expectedVal)) );
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
        TEST_ASSERT( Y_lcl.getRawPtr () != NULL );
        TEST_ASSERT( Y_lcl.getBlockSize () == blockSize );

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
          TEST_ASSERT( Y_lcl(i) == static_cast<Scalar> (static_cast<MT> (expectedVal)) );
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
          TEST_ASSERT( X_lcl.getRawPtr () != NULL );
          TEST_ASSERT( X_lcl.getBlockSize () == blockSize );
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
          TEST_ASSERT( Y_lcl.getRawPtr () != NULL );
          TEST_ASSERT( Y_lcl.getBlockSize () == blockSize );

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
            TEST_ASSERT( Y_lcl(i) == static_cast<Scalar> (static_cast<MT> (expectedVal)) );
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
          TEST_ASSERT( Y_lcl.getRawPtr () != NULL );
          TEST_ASSERT( Y_lcl.getBlockSize () == blockSize );

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
            TEST_ASSERT( Y_lcl(i) == static_cast<Scalar> (static_cast<MT> (expectedVal)) );
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
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( ExpBlockCrsMatrix, write, Scalar, LO, GO, Node )
  {
    typedef Tpetra::Experimental::BlockCrsMatrix<Scalar, LO, GO, Node> BCM;
    typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    // The typedef below is also a test.  BlockCrsMatrix must have
    // this typedef, or this test won't compile.
    typedef typename BCM::little_block_type little_block_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType MT;

    out << "Testing Tpetra::Experimental::BlockCrsMatrix basic "
      "functionality" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const Comm<int> > comm = getDefaultComm ();
    RCP<Node> node = getNode<Node> ();
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
      rcp (new map_type (INVALID, numLocalMeshPoints, indexBase, comm, node));
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
        little_block_type tempBlock (tempBlockPtr, blockSize, blockSize, 1);
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


  // Test BlockCrsMatrix::setAllToScalar.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( ExpBlockCrsMatrix, SetAllToScalar, Scalar, LO, GO, Node )
  {
    typedef Tpetra::Experimental::BlockMultiVector<Scalar, LO, GO, Node> BMV;
    typedef Tpetra::Experimental::BlockCrsMatrix<Scalar, LO, GO, Node> BCM;
    typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;

    out << "Testing Tpetra::Experimental::BlockCrsMatrix::setAllToScalar" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const Comm<int> > comm = getDefaultComm ();
    RCP<Node> node = getNode<Node> ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    out << "Creating mesh row Map" << endl;

    const size_t numLocalMeshPoints = 12;
    const GO indexBase = 1;
    // mfh 16 May 2014: Tpetra::CrsGraph still needs the row Map as an
    // RCP.  Later interface changes will let us pass in the Map by
    // const reference and assume view semantics.
    RCP<const map_type> meshRowMapPtr =
      rcp (new map_type (INVALID, numLocalMeshPoints, indexBase, comm, node));
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
    map_type pointDomainMap = BMV::makePointMap (* (graph.getDomainMap ()), blockSize);
    TEST_ASSERT( ! blockMat.getDomainMap ().is_null () &&
                 pointDomainMap.isSameAs (* (blockMat.getDomainMap ())) );
    map_type pointRangeMap = BMV::makePointMap (* (graph.getRangeMap ()), blockSize);
    TEST_ASSERT( ! blockMat.getRangeMap ().is_null () &&
                 pointRangeMap.isSameAs (* (blockMat.getRangeMap ())) );

    // Test that the result of getGraph() has the same Maps.
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
    const Scalar three = STS::one () + STS::one () + STS::one ();
    blockMat.setAllToScalar (three);

    // Y := A*X, where X is a block multivector (with one column) full
    // of 1s.  Since there are two block entries per row, each of
    // which is all 3s, we know that each entry of the result Y will
    // be 6*blockSize.
    const Scalar requiredValue = static_cast<Scalar> (6 * blockSize);
    BMV X (* (graph.getDomainMap ()), pointDomainMap, blockSize, static_cast<LO> (1));
    X.putScalar (STS::one ());
    BMV Y (* (graph.getRangeMap ()), pointRangeMap, blockSize, static_cast<LO> (1));
    blockMat.applyBlock (X, Y, Teuchos::NO_TRANS, STS::one (), STS::zero ());

    const LO myMinLclMeshRow = Y.getMap ()->getMinLocalIndex ();
    const LO myMaxLclMeshRow = Y.getMap ()->getMaxLocalIndex ();
    for (LO lclMeshRow = myMinLclMeshRow; lclMeshRow <= myMaxLclMeshRow; ++lclMeshRow) {
      typename BMV::little_vec_type Y_lcl = Y.getLocalBlock (lclMeshRow, 0);
      for (LO i = 0; i < blockSize; ++i) {
        TEST_EQUALITY( Y_lcl(i), requiredValue );
      }
    }

    TEST_NOTHROW( blockMat.setAllToScalar (STS::zero ()) );
    blockMat.applyBlock (X, Y, Teuchos::NO_TRANS, STS::one (), STS::zero ());
    for (LO lclMeshRow = myMinLclMeshRow; lclMeshRow <= myMaxLclMeshRow; ++lclMeshRow) {
      typename BMV::little_vec_type Y_lcl = Y.getLocalBlock (lclMeshRow, 0);
      for (LO i = 0; i < blockSize; ++i) {
        TEST_EQUALITY( Y_lcl(i), STS::zero () );
      }
    }

    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
  }


  // Test BlockCrsMatrix Import for the same graphs.  This is really
  // just a test of the "copy" part of copyAndPermute.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( ExpBlockCrsMatrix, ImportCopy, Scalar, LO, GO, Node )
  {
    typedef Tpetra::Experimental::BlockMultiVector<Scalar, LO, GO, Node> BMV;
    typedef Tpetra::Experimental::BlockCrsMatrix<Scalar, LO, GO, Node> BCM;
    typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::Import<LO, GO, Node> import_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;

    out << "Testing Tpetra::Experimental::BlockCrsMatrix Import with same "
      "graphs" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const Comm<int> > comm = getDefaultComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();
    RCP<Node> node = getNode<Node> ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    out << "Creating mesh row Map" << endl;

    const size_t numLocalMeshPoints = 12;
    const GO indexBase = 1;
    // mfh 16 May 2014: Tpetra::CrsGraph still needs the row Map as an
    // RCP.  Later interface changes will let us pass in the Map by
    // const reference and assume view semantics.
    RCP<const map_type> meshRowMapPtr =
      rcp (new map_type (INVALID, numLocalMeshPoints, indexBase, comm, node));
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

    out << "The matrix A1, after construction:" << endl;
    A1.describe (out, Teuchos::VERB_EXTREME);

    // Fill all entries of the second matrix with -2.
    const Scalar minusTwo = -STS::one () - STS::one ();
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
          if (Y_lcl(i) != requiredValue) {
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
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( ExpBlockCrsMatrix, ExportDiffRowMaps, Scalar, LO, GO, Node )
  {
    // typedef Tpetra::Experimental::BlockMultiVector<Scalar, LO, GO, Node> BMV;
    typedef Tpetra::Experimental::BlockCrsMatrix<Scalar, LO, GO, Node> BCM;
    typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::Export<LO, GO, Node> export_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    int lclSuccess = 1;
    int gblSuccess = 1;

    out << "Testing Tpetra::Experimental::BlockCrsMatrix Import "
      "with different graphs with different row Maps" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const Comm<int> > comm = getDefaultComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();
    RCP<Node> node = getNode<Node> ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    out << "Create nonoverlapping mesh row Map" << endl;
    const size_t numLocalMeshPoints = 12;
    const GO indexBase = 1;
    // mfh 16 May 2014: Tpetra::CrsGraph still needs the row Map as an
    // RCP.  Later interface changes will let us pass in the Map by
    // const reference and assume view semantics.
    RCP<const map_type> meshRowMapPtr =
      rcp (new map_type (INVALID, numLocalMeshPoints, indexBase, comm, node));
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
      rcp (new map_type (INVALID, numOverlapMeshPoints, indexBase, comm, node));

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
    } catch (std::exception&) {
      lclSuccess = 0;
    }
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      out << "*** FAILED to get past fillComplete on nonoverlap graph! ***" << endl;
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
    } catch (std::exception&) {
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
    } catch (std::exception&) {
      lclSuccess = 0;
    }
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      out << "*** A_overlap->setAllToScalar(3) FAILED! ***" << endl;
      return;
    }

    out << "The matrix A_overlap, after construction:" << endl;
    A_overlap->describe (out, Teuchos::VERB_EXTREME);

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

    out << "The matrix A, after construction:" << endl;
    A->describe (out, Teuchos::VERB_EXTREME);

    out << "Export A_overlap into A" << endl;
    try {
      A->doExport (*A_overlap, *theExport, Tpetra::REPLACE);
    } catch (std::exception& e) {
      lclSuccess = 0;
      std::ostringstream os;
      os << "Proc " << myRank << ": A->doExport(...) threw an exception: "
         << e.what () << endl;
      std::cerr << os.str ();
    }

    if (A->localError ()) {
      lclSuccess = 0;
      for (int p = 0; p < numProcs; ++p) {
        if (myRank == p && A->localError ()) {
          std::ostringstream os;
          os << "Proc " << myRank << ": A reports local error: "
             << A->errorMessages () << endl;
          std::cerr << os.str ();
        }
        comm->barrier ();
        comm->barrier ();
        comm->barrier ();
      }
    }

    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      out << "*** Export FAILED!" << endl;
      for (int p = 0; p < numProcs; ++p) {
        if (p == myRank) {
          std::ostringstream os;
          os << "Process " << myRank << ": error messages from A_overlap: "
             << A_overlap->errorMessages () << endl
             << "Process " << myRank << ": error messages from A: "
             << A->errorMessages () << endl;
          std::cerr << os.str ();
        }
        comm->barrier (); // give time for output to complete
        comm->barrier ();
        comm->barrier ();
      }
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
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( ExpBlockCrsMatrix, localGSDiagonalMatrix, Scalar, LO, GO, Node )
  {
    typedef Tpetra::Experimental::BlockVector<Scalar, LO, GO, Node> BV;
    typedef Tpetra::Experimental::BlockCrsMatrix<Scalar, LO, GO, Node> BCM;
    typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;

    const Scalar two = STS::one () + STS::one ();
    const Scalar three = STS::one () + STS::one () + STS::one ();

    out << "Testing Tpetra::Experimental::BlockCrsMatrix::localGaussSeidel "
      "with a matrix whose graph is diagonal" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const Comm<int> > comm = getDefaultComm ();
    RCP<Node> node = getNode<Node> ();
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
      rcp (new map_type (INVALID, numLocalMeshPoints, indexBase, comm, node));
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


    Teuchos::Array<Scalar> basematrix(blockSize*blockSize, STS::zero());
    basematrix[0] = two;
    basematrix[2] = three;
    basematrix[3] = three;
    basematrix[4] = two;
    basematrix[7] = three;
    basematrix[8] = two;

    Teuchos::Array<Scalar> baseResidual(blockSize, STS::zero());
    baseResidual[0] = STS::one();
    baseResidual[1] = three;
    baseResidual[2] = -two;

    // FIXME (mfh 25 Aug 2014) This will likely only work with Scalar
    // = float or double.  On the other hand, the author of these
    // tests understood that and only instantiated them for
    // Scalar=double (see the instantiations list below).
    Teuchos::Array<Scalar> exactSolution(blockSize, STS::zero());
    exactSolution[0] = 43.0/35.0;
    exactSolution[1] = -12.0/35.0;
    exactSolution[2] = -17.0/35.0;

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
      blockMat.replaceLocalValues(lclRowInd, lclColInds.getRawPtr(), &basematrix[0], 1);
      residual.replaceLocalValues(lclRowInd, &baseResidual[0]);
      solution.replaceLocalValues(lclRowInd, &baseResidual[0]);
    }

    BCM diagonalMat(graph, blockSize);

    Teuchos::ArrayRCP<size_t> diagonalOffsets(numLocalMeshPoints);
    blockMat.getLocalDiagOffsets(diagonalOffsets);
    blockMat.getLocalDiagCopy(diagonalMat, diagonalOffsets());

    Scalar* blockVals;
    Scalar* diagVals;
    const LO* blkColInds;
    const LO* diagColInds;
    LO blkNumColInds;
    LO diagNumColInds;

    Teuchos::Array<int> pivots (blockSize*numLocalMeshPoints+1, 1);
    out << "pivots size = " << pivots.size() << endl;

    int* ipiv = pivots.getRawPtr ();
    for (LO lclRowInd = meshRowMap.getMinLocalIndex ();
         lclRowInd <= meshRowMap.getMaxLocalIndex (); ++lclRowInd) {
      blockMat.getLocalRowView (lclRowInd, blkColInds, blockVals, blkNumColInds);
      diagonalMat.getLocalRowView (lclRowInd, diagColInds, diagVals, diagNumColInds);
      for (LO k = 0; k < blockSize*blockSize; ++k) {
        TEST_EQUALITY( blockVals[k], diagVals[k] );
      }

      typename BCM::little_block_type diagBlock =
        diagonalMat.getLocalBlock(lclRowInd, lclRowInd);
      const LO pivotOffset = blockSize * (lclRowInd - meshRowMap.getMinLocalIndex ());
      int info = -5;
      diagBlock.factorize (&ipiv[pivotOffset], info);
    }

    blockMat.localGaussSeidel (residual, solution, diagonalMat, &pivots[0],
                               STS::one(), Tpetra::Forward);

    for (LO lclRowInd = meshRowMap.getMinLocalIndex ();
         lclRowInd <= meshRowMap.getMaxLocalIndex (); ++lclRowInd) {
      typename BV::little_vec_type xlcl = solution.getLocalBlock (lclRowInd);
      Scalar* x = xlcl.getRawPtr ();
      out << "row = " << lclRowInd << endl;
      for (LO k = 0; k < blockSize; ++k) {
        TEST_FLOATING_EQUALITY( x[k], exactSolution[k], 1e-12 );
        x[k] = -STS::one ();
      }
    }

    blockMat.localGaussSeidel (residual, solution, diagonalMat,
                               pivots.getRawPtr (), STS::one (),
                               Tpetra::Backward);
    for (LO lclRowInd = meshRowMap.getMinLocalIndex ();
         lclRowInd <= meshRowMap.getMaxLocalIndex (); ++lclRowInd) {
      typename BV::little_vec_type xlcl = solution.getLocalBlock (lclRowInd);
      Scalar* x = xlcl.getRawPtr ();
      for (LO k = 0; k < blockSize; ++k) {
        TEST_FLOATING_EQUALITY( x[k], exactSolution[k], 1e-12 );
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
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( ExpBlockCrsMatrix, localGSTriangularMatrices, Scalar, LO, GO, Node )
  {
    typedef Tpetra::Experimental::BlockVector<Scalar, LO, GO, Node> BV;
    typedef Tpetra::Experimental::BlockCrsMatrix<Scalar, LO, GO, Node> BCM;
    typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;

    const Scalar two = STS::one () + STS::one ();
    const Scalar three = STS::one () + STS::one () + STS::one ();

    out << "Testing Tpetra::Experimental::BlockCrsMatrix::localGaussSeidel "
      "with a triangular matrix ( ??? )" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const Comm<int> > comm = getDefaultComm ();
    RCP<Node> node = getNode<Node> ();
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
      rcp (new map_type (INVALID, numLocalMeshPoints, indexBase, comm, node));
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
    blockMat.computeDiagonalGraph();
    BCM diagonalMat;
    TEST_NOTHROW( diagonalMat = BCM (* (blockMat.getDiagonalGraph ()), blockSize) );


    Teuchos::Array<Scalar> basematrix (maxNumEntPerRow * maxNumEntPerRow,
                                       STS::zero ());
    basematrix[0] = two;
    basematrix[3] = three;
    basematrix[4] = two;
    basematrix[6] = STS::one ();
    basematrix[7] = three;
    basematrix[8] = two;

    Teuchos::Array<Scalar> baseResidual (maxNumEntPerRow, STS::zero ());
    baseResidual[0] = STS::one ();
    baseResidual[1] = three;
    baseResidual[2] = -two;

    // FIXME (mfh 25 Aug 2014) This will likely only work with Scalar
    // = float or double.  On the other hand, the author of these
    // tests understood that and only instantiated them for
    // Scalar=double (see the instantiations list below).
    Teuchos::Array<Scalar> exactSolution (maxNumEntPerRow, STS::zero ());
    exactSolution[0] = 0.5;
    exactSolution[1] = 0.75;
    exactSolution[2] = -19.0/8.0;

    Teuchos::Array<Scalar> assembleMatrix (blockSize * blockSize, STS::zero ());
    Teuchos::Array<Scalar> assembleResidual (blockSize, STS::zero ());

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

    Teuchos::ArrayRCP<size_t> diagonalOffsets(numLocalMeshPoints);
    blockMat.getLocalDiagOffsets(diagonalOffsets);
    blockMat.getLocalDiagCopy(diagonalMat, diagonalOffsets());

    Scalar* blockVals;
    Scalar* diagVals;

    Teuchos::Array<int> pivots(blockSize*numLocalMeshPoints+1, Teuchos::OrdinalTraits<int>::one());

    int* ipiv = pivots.getRawPtr ();
    for (LO lclRowInd = meshRowMap.getMinLocalIndex ();
         lclRowInd <= meshRowMap.getMaxLocalIndex (); ++lclRowInd) {
      typename BCM::little_block_type diagBlock =
        diagonalMat.getLocalBlock (lclRowInd, lclRowInd);
      typename BCM::little_block_type block =
        blockMat.getLocalBlock (lclRowInd, lclRowInd);

      diagVals = diagBlock.getRawPtr ();
      blockVals = block.getRawPtr ();
      for (LO k = 0; k < blockSize * blockSize; ++k) {
        TEST_EQUALITY( blockVals[k], diagVals[k] );
      }
      const LO pivotOffset = blockSize * (lclRowInd - meshRowMap.getMinLocalIndex ());
      int info = -5;
      diagBlock.factorize (&ipiv[pivotOffset], info);
    }

    blockMat.localGaussSeidel (residual, solution, diagonalMat, &pivots[0],
                               STS::one (), Tpetra::Forward);

    for (LO lclRowInd = meshRowMap.getMinLocalIndex ();
         lclRowInd <= meshRowMap.getMaxLocalIndex(); ++lclRowInd) {
      const LO rowOffset = lclRowInd - meshRowMap.getMinLocalIndex ();
      typename BV::little_vec_type xlcl = solution.getLocalBlock (lclRowInd);
      Scalar* x = xlcl.getRawPtr ();
      for (LO k = 0; k < blockSize; ++k) {
        TEST_FLOATING_EQUALITY( x[k], exactSolution[rowOffset], 1e-12 );
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

    blockMat.getLocalDiagCopy (diagonalMat, diagonalOffsets ());

    ipiv = pivots.getRawPtr();
    for (LO lclRowInd = meshRowMap.getMinLocalIndex ();
         lclRowInd <= meshRowMap.getMaxLocalIndex (); ++lclRowInd) {
      typename BCM::little_block_type diagBlock =
        diagonalMat.getLocalBlock(lclRowInd, lclRowInd);
      typename BCM::little_block_type block =
        blockMat.getLocalBlock(lclRowInd, lclRowInd);
      diagVals = diagBlock.getRawPtr ();
      blockVals = block.getRawPtr ();
      for (LO k = 0; k < blockSize*blockSize; ++k) {
        TEST_EQUALITY( blockVals[k], diagVals[k] );
      }
      const LO pivotOffset = blockSize * (lclRowInd - meshRowMap.getMinLocalIndex ());
      int info = -5;
      diagBlock.factorize (&ipiv[pivotOffset], info);
    }

    blockMat.localGaussSeidel (residual, solution, diagonalMat, &pivots[0],
                               STS::one (), Tpetra::Symmetric);

    for (LO lclRowInd = meshRowMap.getMinLocalIndex ();
         lclRowInd <= meshRowMap.getMaxLocalIndex(); ++lclRowInd) {
      const LO rowOffset = lclRowInd - meshRowMap.getMinLocalIndex ();
      typename BV::little_vec_type xlcl = solution.getLocalBlock (lclRowInd);
      Scalar* x = xlcl.getRawPtr ();
      for (LO k = 0; k < blockSize; ++k) {
        TEST_FLOATING_EQUALITY( x[k], exactSolution[rowOffset], 1e-12 );
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
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( ExpBlockCrsMatrix, ctor, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( ExpBlockCrsMatrix, basic, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( ExpBlockCrsMatrix, write, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( ExpBlockCrsMatrix, SetAllToScalar, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( ExpBlockCrsMatrix, ImportCopy, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( ExpBlockCrsMatrix, ExportDiffRowMaps, SCALAR, LO, GO, NODE )

# define UNIT_TEST_GROUP_LGN( LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( ExpBlockCrsMatrix, localGSDiagonalMatrix, double, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( ExpBlockCrsMatrix, localGSTriangularMatrices, double, LO, GO, NODE )


  TPETRA_ETI_MANGLING_TYPEDEFS()

  // FIXME (mfh 30 Jul 2014) We might like to use Teuchos::LAPACK in
  // the implementation of BlockCrsMatrix, so it wouldn't make sense
  // to do explicit instantiation for Scalar types that
  // Teuchos::LAPACK doesn't support.  It seems that the *_TESTMV
  // macro does whaat we want, but it would make sense to have a macro
  // like TPETRA_INSTANTIATE_FLOATINGPOINT or
  // TPETRA_INSTANTIATE_LAPACK_TYPES.
  //
  // NOTE (mfh 30 Jul 2014) BlockCrsMatrix's explicit instantiation
  // (in ../../src/Tpetra_Experimental_BlockCrsMatrix.cpp) must also
  // use this macro, so that the class itself and the tests are
  // explicitly instantiated over the same set of types.  I have also
  // put a note there that points here.

  TPETRA_INSTANTIATE_TESTMV( UNIT_TEST_GROUP )
  TPETRA_INSTANTIATE_LGN( UNIT_TEST_GROUP_LGN )

} // namespace (anonymous)


