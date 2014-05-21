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

// Some Macro Magic to ensure that if CUDA and KokkosCompat is enabled
// only the .cu version of this file is actually compiled
#include <Tpetra_config.h>
#ifdef HAVE_TPETRA_KOKKOSCOMPAT
#include <KokkosCore_config.h>
#ifdef KOKKOS_USE_CUDA_BUILD
  #define DO_COMPILATION
#else
  #ifndef KOKKOS_HAVE_CUDA
    #define DO_COMPILATION
  #endif
#endif
#else
  #define DO_COMPILATION
#endif

#ifdef DO_COMPILATION

#include <Tpetra_TestingUtilities.hpp>
#include <Tpetra_Experimental_BlockCrsMatrix.hpp>
#include <Tpetra_Experimental_BlockVector.hpp>

namespace {

  //
  // UNIT TESTS
  //

  // Test BlockCrsMatrix's constructors.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( ExpBlockCrsMatrix, ctor, Scalar, LO, GO, Node )
  {
    using Tpetra::TestingUtilities::getNode;
    using Tpetra::TestingUtilities::getDefaultComm;
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::endl;
    typedef Tpetra::Experimental::BlockMultiVector<Scalar, LO, GO, Node> BMV;
    typedef Tpetra::Experimental::BlockCrsMatrix<Scalar, LO, GO, Node> BCM;
    typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::global_size_t GST;

    out << "Testing Tpetra::Experimental::BlockCrsMatrix" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const Comm<int> > comm = getDefaultComm ();
    RCP<Node> node = getNode<Node> ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    out << "Creating mesh row Map" << endl;

    const size_t numLocalMeshPoints = 12;
    const GO indexBase = 0;
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

    // Test that the result of getGraph() has the same Maps.
    {
      graph_type graph2 = blockMat.getGraph ();
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
  }

  // Test creating a BlockCrsMatrix with a nontrivial graph.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( ExpBlockCrsMatrix, NontrivialGraph, Scalar, LO, GO, Node )
  {
    using Tpetra::TestingUtilities::getNode;
    using Tpetra::TestingUtilities::getDefaultComm;
    using Teuchos::Array;
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::endl;
    typedef Tpetra::Experimental::BlockMultiVector<Scalar, LO, GO, Node> BMV;
    typedef Tpetra::Experimental::BlockVector<Scalar, LO, GO, Node> BV;
    typedef Tpetra::Experimental::BlockCrsMatrix<Scalar, LO, GO, Node> BCM;
    typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::global_size_t GST;
    // The typedef below is also a test.  BlockCrsMatrix must have
    // this typedef, or this test won't compile.
    typedef typename BCM::little_block_type little_block_type;
    typedef typename BV::little_vec_type little_vec_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;

    const Scalar two = STS::one () + STS::one ();
    const Scalar three = STS::one () + STS::one () + STS::one ();

    out << "Testing Tpetra::Experimental::BlockCrsMatrix creation with a "
      "nontrivial graph" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const Comm<int> > comm = getDefaultComm ();
    RCP<Node> node = getNode<Node> ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    out << "Creating mesh row Map" << endl;

    const size_t numLocalMeshPoints = 12;
    const GO indexBase = 0;
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
        gblColInds[k] = (gblRowInd + static_cast<GO> (k + 1)) %
          static_cast<GO> (globalNumRows);
      }
      graph.insertGlobalIndices (gblRowInd, gblColInds ());
    }
    graph.fillComplete ();

    // Get the graph's column Map (the "mesh column Map").
    TEST_ASSERT( ! graph.getColMap ().is_null () );
    map_type meshColMap = * (graph.getColMap ());

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

    // Test that the result of getGraph() has the same Maps.
    {
      graph_type graph2 = blockMat.getGraph ();
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

    // // Fill the matrix with zeros.
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
        const GO gblColInd = (gblRowInd + static_cast<GO> (k + 1)) %
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
            tempBlock(i,j) = static_cast<Scalar> (j + i * blockSize);
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
            TEST_ASSERT( curBlk(i,j) == static_cast<Scalar> (j + i * blockSize) );
          }
        }
      } // for each entry in the row
    } // for each local row

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
          X_lcl(i) = static_cast<Scalar> (blockSize - i);
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
          expectedVal *= two;
          out << "Y_lcl(" << i << ") = " << Y_lcl(i)
              << "; expectedVal = " << expectedVal << std::endl;
          TEST_ASSERT( Y_lcl(i) == static_cast<Scalar> (expectedVal) );
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
          expectedVal *= (two * alpha);
          expectedVal += beta * STS::one ();
          out << "Y_lcl(" << i << ") = " << Y_lcl(i)
              << "; expectedVal = " << expectedVal << std::endl;
          TEST_ASSERT( Y_lcl(i) == static_cast<Scalar> (expectedVal) );
        }
      }
    } // done with single-vector applyBlock test


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
            X_lcl(i) = static_cast<Scalar> ((blockSize - i) * (j + 1));
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
            expectedVal *= two;
            expectedVal *= static_cast<Scalar> (col + 1);
            out << "Y_lcl(" << i << ") = " << Y_lcl(i)
                << "; expectedVal = " << expectedVal << std::endl;
            TEST_ASSERT( Y_lcl(i) == static_cast<Scalar> (expectedVal) );
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
            expectedVal *= (two * alpha);
            expectedVal *= static_cast<Scalar> (col + 1);
            expectedVal += beta * STS::one ();
            out << "Y_lcl(" << i << ") = " << Y_lcl(i)
                << "; expectedVal = " << expectedVal << std::endl;
            TEST_ASSERT( Y_lcl(i) == static_cast<Scalar> (expectedVal) );
          }
        } // for each column (vector) of the BlockMultiVector
      } // for each local (mesh) row of the output BlockMultiVector
    } // done with multiple-vector applyBlock test

    // Finishing with a barrier ensures that the test won't finish on
    // Process 0 (and therefore report a "pass") if there is deadlock
    // somewhere.
    comm->barrier ();
    out << "Hooray, got to the end of the BlockCrsMatrix test!" << std::endl;
  }



//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( ExpBlockCrsMatrix, ctor, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( ExpBlockCrsMatrix, NontrivialGraph, SCALAR, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV_NOGPU( UNIT_TEST_GROUP )

} // namespace (anonymous)

#endif  //DO_COMPILATION

