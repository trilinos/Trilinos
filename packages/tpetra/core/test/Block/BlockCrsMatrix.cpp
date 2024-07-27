// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_BlockCrsMatrix.hpp"
#include "Tpetra_BlockCrsMatrix_Helpers.hpp"
#include "Tpetra_BlockVector.hpp"
#include "Tpetra_BlockView.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"
#include "Tpetra_Details_gathervPrint.hpp"

namespace {
  using Tpetra::TestingUtilities::getDefaultComm;
  using Tpetra::TestingUtilities::arcp_from_view;
  using Tpetra::Details::gathervPrint;
  using Teuchos::Array;
  using Teuchos::Comm;
  using Teuchos::outArg;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ScalarTraits;
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
    graph_type graph (meshRowMapPtr, 0);
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

    using lids_type = typename graph_type::nonconst_local_inds_host_view_type;
    using gids_type = typename graph_type::nonconst_global_inds_host_view_type;
    using vals_type = typename BCM::nonconst_values_host_view_type;
    using local_inds_host_view_type = typename BCM::local_inds_host_view_type;
    using values_host_view_type = typename BCM::values_host_view_type;
    using impl_scalar_type = typename BCM::impl_scalar_type;

    typedef typename BCM::little_block_host_type little_block_host_type;
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
    graph_type graph (meshRowMapPtr, maxNumEntPerRow);

    // Fill the graph.
    gids_type gblColInds ("gblColIds",maxNumEntPerRow);
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
      graph.insertGlobalIndices (gblRowInd, gblColInds.extent(0),gblColInds.data());
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

    // KK: not meaningfule test
    // {
    //   auto val = blockMat.getValuesHost(); 
    //   static_assert (std::is_same<typename decltype (val)::execution_space::memory_space,
    //                  Kokkos::HostSpace>::value,
    //                  "Host View is not actually a host View.");
    // }
    // {
    //   auto val = blockMat.getValuesHostNonConst ();
    //   static_assert (std::is_same<typename decltype (val)::execution_space::memory_space,
    //                  Kokkos::HostSpace>::value,
    //                  "Host View is not actually a host View.");
    // }

    Array<Scalar> tempBlockSpace (maxNumEntPerRow * entriesPerBlock);

    // Test that getLocalRowView returns the right column indices.

    lids_type lclColInds ("lclColInds",maxNumEntPerRow);
    lids_type myLclColIndsCopy ("myLclColIndsCopy",maxNumEntPerRow);
    vals_type myValsCopy ("myValsCopy",maxNumEntPerRow*entriesPerBlock);
    lids_type myLclColIndsSorted ("myLclColIndsSorted",maxNumEntPerRow);
    for (LO lclRowInd = meshRowMap.getMinLocalIndex ();
         lclRowInd <= meshRowMap.getMaxLocalIndex (); ++lclRowInd) {
      local_inds_host_view_type myLclColInds;
      values_host_view_type myVals;
      LO numEnt = 0;
      blockMat.getLocalRowView (lclRowInd, myLclColInds, myVals); numEnt = myLclColInds.extent(0);

      TEST_ASSERT( numEnt == static_cast<LO> (maxNumEntPerRow) );
      TEST_ASSERT( myLclColInds.data() != NULL );
      TEST_ASSERT( myVals.data() != NULL );

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
      Tpetra::sort (lclColInds, lclColInds.extent(0));
      std::copy (myLclColInds.data(), myLclColInds.data() + 2, arcp_from_view(myLclColIndsSorted).begin());
      Tpetra::sort (myLclColIndsSorted, myLclColIndsSorted.extent(0));
      TEST_COMPARE_ARRAYS( lclColInds, myLclColIndsSorted );

      // Test that getLocalRowCopy works.
      size_t numEntries;
      blockMat.getLocalRowCopy (lclRowInd, myLclColIndsCopy, myValsCopy, numEntries);
      numEnt = static_cast<LO>(numEntries);
      TEST_ASSERT( numEnt == static_cast<LO> (maxNumEntPerRow) );

      // CrsGraph doesn't technically need to promise to sort by local
      // column indices, so we sort both arrays before comparing.
      Kokkos::deep_copy(myLclColIndsSorted,Kokkos::subview(myLclColIndsCopy,std::make_pair(0,2)));
      //      std::copy (myLclColIndsCopy.getRawPtr(), myLclColIndsCopy.getRawPtr() + 2, myLclColIndsSorted.begin ());
      Tpetra::sort (myLclColIndsSorted, myLclColIndsSorted.extent(0));
      TEST_COMPARE_ARRAYS( lclColInds, myLclColIndsSorted );

      // Fill the entries in the row with zeros.
      std::fill (tempBlockSpace.begin (), tempBlockSpace.end (), STS::zero ());
      int err = blockMat.replaceLocalValues (lclRowInd, lclColInds.data(),
                                         tempBlockSpace.getRawPtr (), numEnt);
      TEST_ASSERT( err == numEnt );
      // Make sure that the input Scalar values didn't change (are
      // still all zero).
      for (LO k = 0; k < numEnt; ++k) {
        Scalar* const tempBlockPtr = tempBlockSpace.getRawPtr () +
          k * blockSize * blockSize;
        little_block_host_type tempBlock ((typename little_block_host_type::value_type*) tempBlockPtr, blockSize, blockSize);
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
        little_block_host_type tempBlock ((typename little_block_host_type::value_type*) tempBlockPtr, blockSize, blockSize);
        for (LO j = 0; j < blockSize; ++j) {
          for (LO i = 0; i < blockSize; ++i) {
            tempBlock(i,j) = static_cast<Scalar> (static_cast<MT> (j + i * blockSize));
          }
        }
      } // for each entry in the row
      err = blockMat.replaceLocalValues (lclRowInd, lclColInds.data(),
                                         tempBlockSpace.getRawPtr (), numEnt);
      TEST_ASSERT( err == numEnt );

      // Get a view of the current row again, and test that the
      // entries were modified as expected.  This tests that the
      // method assumes that the input blocks are row major.
      blockMat.getLocalRowView (lclRowInd, myLclColInds, myVals); numEnt = static_cast<LO>(myLclColInds.extent(0));


      TEST_ASSERT( numEnt == static_cast<LO> (maxNumEntPerRow) );
      TEST_ASSERT( myLclColInds.data() != NULL );
      TEST_ASSERT( myVals.data() != NULL );
      for (LO k = 0; k < numEnt; ++k) {
        impl_scalar_type* curBlkPtr = const_cast<impl_scalar_type*>(reinterpret_cast<const impl_scalar_type*>(myVals.data())) + k * blockSize * blockSize;
        little_block_host_type curBlk ((typename little_block_host_type::value_type*) curBlkPtr, blockSize, blockSize);

        for (LO j = 0; j < blockSize; ++j) {
          for (LO i = 0; i < blockSize; ++i) {
            TEST_ASSERT( static_cast<Scalar> (curBlk(i,j)) ==
                         static_cast<Scalar> (static_cast<MT> (j + i * blockSize)) );
          }
        }
      } // for each entry in the row
    } // for each local row

    // KK: not meaningfule test; will be deprecated
//     {
//       auto val = blockMat.template getValues<typename Node::device_type::memory_space> ();
//       // "Device" View may live in CudaUVMSpace.
// #if defined(KOKKOS_ENABLE_CUDA)
//       constexpr bool testing_cuda =
//         std::is_same<typename Node::device_type::execution_space, Kokkos::Cuda>::value;
//       static_assert (! testing_cuda ||
//                      std::is_same<typename decltype (val)::execution_space, Kokkos::Cuda>::value,
//                      "Device View is not actually a Device View.");
//       auto val2 = blockMat.getValuesDevice ();
//       static_assert (! testing_cuda ||
//                      std::is_same<typename decltype (val2)::execution_space, Kokkos::Cuda>::value,
//                      "Device View is not actually a Device View.");
// #endif // defined(KOKKOS_ENABLE_CUDA)
//     }

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

      //Y.sync_host(); // X and Y are same map and write to X_lcl(i) needs fence

      const map_type& meshDomainMap = * (graph.getDomainMap ());
      for (LO lclDomIdx = meshDomainMap.getMinLocalIndex ();
           lclDomIdx <= meshDomainMap.getMaxLocalIndex (); ++lclDomIdx) {
        auto X_lcl = X.getLocalBlockHost (lclDomIdx, Tpetra::Access::OverwriteAll);
        TEST_ASSERT( X_lcl.data () != NULL );
        TEST_ASSERT( static_cast<size_t> (X_lcl.extent (0)) == static_cast<size_t> (blockSize) );
        for (LO i = 0; i < blockSize; ++i) {
          X_lcl(i) = static_cast<Scalar> (static_cast<MT> (blockSize - i));
        }
      }

      TEST_NOTHROW( blockMat.applyBlock (X, Y) );
      Kokkos::fence ();

      const map_type& meshRangeMap = * (graph.getRangeMap ());
      for (LO lclRanIdx = meshRangeMap.getMinLocalIndex ();
           lclRanIdx <= meshRangeMap.getMaxLocalIndex (); ++lclRanIdx) {
        auto Y_lcl = Y.getLocalBlockHost (lclRanIdx, Tpetra::Access::ReadOnly);
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
      Kokkos::fence ();

      for (LO lclRanIdx = meshRangeMap.getMinLocalIndex ();
           lclRanIdx <= meshRangeMap.getMaxLocalIndex (); ++lclRanIdx) {
        auto Y_lcl = Y.getLocalBlockHost (lclRanIdx, Tpetra::Access::ReadOnly);
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

      //Y.sync_host(); // X and Y are same map and write to X_lcl(i) needs fence

      const map_type& meshDomainMap = * (graph.getDomainMap ());
      for (LO lclDomIdx = meshDomainMap.getMinLocalIndex ();
           lclDomIdx <= meshDomainMap.getMaxLocalIndex (); ++lclDomIdx) {
        for (LO j = 0; j < numVecs; ++j) {
          auto X_lcl = X.getLocalBlockHost(lclDomIdx, j, Tpetra::Access::OverwriteAll);
          TEST_ASSERT( X_lcl.data () != NULL );
          TEST_ASSERT( static_cast<size_t> (X_lcl.extent (0)) == static_cast<size_t> (blockSize) );
          for (LO i = 0; i < blockSize; ++i) {
            X_lcl(i) = static_cast<Scalar> (static_cast<MT> ((blockSize - i) * (j + 1)));
          }
        }
      }

      TEST_NOTHROW( blockMat.applyBlock (X, Y) );
      Kokkos::fence ();

      const map_type& meshRangeMap = * (graph.getRangeMap ());
      for (LO lclRanIdx = meshRangeMap.getMinLocalIndex ();
           lclRanIdx <= meshRangeMap.getMaxLocalIndex (); ++lclRanIdx) {
        for (LO col = 0; col < numVecs; ++col) {
          auto Y_lcl = Y.getLocalBlockHost (lclRanIdx, col, Tpetra::Access::ReadOnly);
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
      Kokkos::fence ();

      for (LO lclRanIdx = meshRangeMap.getMinLocalIndex ();
           lclRanIdx <= meshRangeMap.getMaxLocalIndex (); ++lclRanIdx) {
        for (LO col = 0; col < numVecs; ++col) {
          auto Y_lcl = Y.getLocalBlockHost (lclRanIdx, col, Tpetra::Access::ReadOnly);
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

      //Y.sync_host(); // X and Y are same map and write to X_lcl(i) needs fence

      const map_type& meshDomainMap = * (graph.getDomainMap ());
      for (LO lclDomIdx = meshDomainMap.getMinLocalIndex ();
           lclDomIdx <= meshDomainMap.getMaxLocalIndex (); ++lclDomIdx) {
        auto X_lcl = X.getLocalBlockHost (lclDomIdx, Tpetra::Access::OverwriteAll);
        TEST_ASSERT( X_lcl.data () != NULL );
        TEST_ASSERT( static_cast<size_t> (X_lcl.extent (0)) == static_cast<size_t> (blockSize) );
        for (LO i = 0; i < blockSize; ++i) {
          X_lcl(i) = static_cast<Scalar> (static_cast<MT> (blockSize - i));
        }
      }

      vec_type X_vec = X.getVectorView ();
      vec_type Y_vec = Y.getVectorView ();

      TEST_NOTHROW( blockMat.apply (X_vec, Y_vec) );
      Kokkos::fence ();

      // This test also exercises whether getVectorView really does
      // return a view, since we access and test results using the
      // BlockVector.

      const map_type& meshRangeMap = * (graph.getRangeMap ());
      for (LO lclRanIdx = meshRangeMap.getMinLocalIndex ();
           lclRanIdx <= meshRangeMap.getMaxLocalIndex (); ++lclRanIdx) {
        auto Y_lcl = Y.getLocalBlockHost (lclRanIdx, Tpetra::Access::ReadOnly);
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
      Kokkos::fence ();

      for (LO lclRanIdx = meshRangeMap.getMinLocalIndex ();
           lclRanIdx <= meshRangeMap.getMaxLocalIndex (); ++lclRanIdx) {
        auto Y_lcl = Y.getLocalBlockHost (lclRanIdx, Tpetra::Access::ReadOnly);
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

      //Y.sync_host(); // X and Y are same map and write to X_lcl(i) needs fence

      const map_type& meshDomainMap = * (graph.getDomainMap ());
      for (LO lclDomIdx = meshDomainMap.getMinLocalIndex ();
           lclDomIdx <= meshDomainMap.getMaxLocalIndex (); ++lclDomIdx) {
        for (LO j = 0; j < numVecs; ++j) {
          auto X_lcl = X.getLocalBlockHost(lclDomIdx, j, Tpetra::Access::OverwriteAll);
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
      Kokkos::fence ();

      // This test also exercises whether getMultiVectorView really
      // does return a view, since we access and test results using
      // the BlockMultiVector.

      const map_type& meshRangeMap = * (graph.getRangeMap ());
      for (LO lclRanIdx = meshRangeMap.getMinLocalIndex ();
           lclRanIdx <= meshRangeMap.getMaxLocalIndex (); ++lclRanIdx) {
        for (LO col = 0; col < numVecs; ++col) {
          auto Y_lcl = Y.getLocalBlockHost (lclRanIdx, col, Tpetra::Access::ReadOnly);
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
      Kokkos::fence ();

      for (LO lclRanIdx = meshRangeMap.getMinLocalIndex ();
           lclRanIdx <= meshRangeMap.getMaxLocalIndex (); ++lclRanIdx) {
        for (LO col = 0; col < numVecs; ++col) {
          auto Y_lcl = Y.getLocalBlockHost (lclRanIdx, col, Tpetra::Access::ReadOnly);
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
    typedef typename BCM::little_block_host_type little_block_host_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType MT;

    using local_inds_host_view_type = typename BCM::local_inds_host_view_type;
    using values_host_view_type = typename BCM::values_host_view_type;

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
    graph_type graph (meshRowMapPtr, maxNumEntPerRow);

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

      local_inds_host_view_type myLclColInds;
      values_host_view_type myVals;
      LO numEnt = 0;
      blockMat.getLocalRowView (lclRowInd, myLclColInds, myVals); numEnt = myLclColInds.extent(0);

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
      std::copy (myLclColInds.data(), myLclColInds.data() + 2, myLclColIndsCopy.begin ());
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
        little_block_host_type tempBlock ((typename little_block_host_type::value_type*) tempBlockPtr, blockSize, blockSize);
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
    Tpetra::blockCrsMatrixWriter(blockMat, "savedBlockMatrix.m", pl);
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
    typedef typename BCM::little_block_host_type little_block_host_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;

    using local_inds_host_view_type = typename BCM::local_inds_host_view_type;
    using values_host_view_type = typename BCM::values_host_view_type;
    using nonconst_values_host_view_type = typename BCM::nonconst_values_host_view_type;

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
    graph_type graph (meshRowMapPtr, maxNumEntPerRow);

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
    Kokkos::fence ();
    try {
      graph.getLocalDiagOffsets (diagMeshOffsets);
    } catch (std::exception& e) {
      success = false;
      std::ostringstream os;
      os << "Proc " << comm->getRank () << ": getLocalDiagOffsets "
        "threw an exception: " << e.what () << endl;
      std::cerr << os.str ();
    }
    Kokkos::fence ();
    auto diagMeshOffsetsHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), diagMeshOffsets);

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
      auto localGraph = graph.getLocalGraphDevice ();
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

          local_inds_host_view_type lclColInds;
          values_host_view_type lclVals;
          LO numEnt = 0;
          blockMat.getLocalRowView (lclRowInd, lclColInds, lclVals); numEnt = lclColInds.extent(0);
          {
            const size_t offset = diagMeshOffsetsHost[lclRowInd];
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
      local_inds_host_view_type lclColInds;
      nonconst_values_host_view_type myVals;
      LO numEnt = 0;
      blockMat.getLocalRowViewNonConst (lclRowInd, lclColInds, myVals); numEnt = lclColInds.extent(0);
      
      // Fill the diagonal block D such that D(i,j) = (lclRowInd+1) *
      // (1 + i + j*blockSize).  Fill the off-diagonal block with -1.
      // This ensures that we can tell we got the right blocks, and
      // that we copied them in the correct order.
      for (LO k = 0; k < numEnt; ++k) {
        const LO offset = blockSize * blockSize * k;
        little_block_host_type curBlock (const_cast<IST*> (myVals.data()) + offset,
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
    Kokkos::fence ();
    blockMat.getLocalDiagCopy (diagBlocks, diagMeshOffsets);
    Kokkos::fence ();
    auto diagBlocksHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), diagBlocks);

    bool allBlocksGood = true;
    for (LO lclRowInd = 0; lclRowInd < static_cast<LO> (numLclMeshPoints); ++lclRowInd) {
      const GO gblRowInd = meshRowMap.getGlobalElement (lclRowInd);
      local_inds_host_view_type lclColInds;
      values_host_view_type myVals;
      LO numEnt = 0;
      blockMat.getLocalRowView (lclRowInd, lclColInds, myVals); numEnt = lclColInds.extent(0);

      // Make sure that the diagonal blocks from getLocalDiagCopy
      // match those in the matrix.
      for (LO k = 0; k < numEnt; ++k) {
        const LO offset = blockSize * blockSize * k;
        little_block_host_type curBlock (const_cast<IST*> (myVals.data()) + offset,
                                    blockSize, blockSize); // row major
        const GO gblColInd = meshColMap.getGlobalElement (lclColInds[k]);
        if (gblColInd == gblRowInd) { // the diagonal block
          auto diagBlock = subview (diagBlocksHost, lclRowInd, ALL (), ALL ());
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
    graph_type graph (meshRowMapPtr, 2);

    if (meshRowMapPtr->getLocalNumElements () > 0) {
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
    Kokkos::fence ();

    out << "Make sure applyBlock got the right answer" << endl;
    const LO myMinLclMeshRow = Y.getMap ()->getMinLocalIndex ();
    const LO myMaxLclMeshRow = Y.getMap ()->getMaxLocalIndex ();
    for (LO lclMeshRow = myMinLclMeshRow; lclMeshRow <= myMaxLclMeshRow; ++lclMeshRow) {
      auto Y_lcl = Y.getLocalBlockHost (lclMeshRow, 0, Tpetra::Access::ReadOnly);
      for (LO i = 0; i < blockSize; ++i) {
        TEST_EQUALITY( static_cast<Scalar> (Y_lcl(i)), requiredValue );
      }
    }

    TEST_NOTHROW( blockMat.setAllToScalar (STS::zero ()) );
    blockMat.applyBlock (X, Y, Teuchos::NO_TRANS, STS::one (), STS::zero ());
    Kokkos::fence ();

    for (LO lclMeshRow = myMinLclMeshRow; lclMeshRow <= myMaxLclMeshRow; ++lclMeshRow) {
      auto Y_lcl = Y.getLocalBlockHost (lclMeshRow, 0, Tpetra::Access::ReadOnly);
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
    using Tpetra::Details::gathervPrint;
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
    graph_type graph (meshRowMapPtr, 2);

    if (meshRowMapPtr->getLocalNumElements () > 0) {
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
    out << "Fill all entries of A2 with " << minusTwo << endl;
    A2.setAllToScalar (minusTwo);

    out << "The matrix A2, after construction:" << endl;
    A2.describe (out, Teuchos::VERB_EXTREME);

    out << "Create the Import" << endl;
    import_type imp (graph.getMap (), graph.getMap ());

    out << "Import A1 into A2" << endl;
    bool importFailedGlobally = false;
    bool importThrewLocally = false;
    std::string exceptionMessage;
    try {
      // The CombineMode doesn't matter for this example, since it
      // amounts to a matrix copy.  We use ADD arbitrarily.
      A2.doImport (A1, imp, Tpetra::ADD);
    }
    catch (std::exception& e) {
      importThrewLocally = true;
      exceptionMessage = e.what ();
    }

    int lclImportSuccess = importThrewLocally ? 0 : 1;
    int gblImportSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclImportSuccess, outArg (gblImportSuccess));
    TEST_ASSERT( gblImportSuccess == 1 );
    if (gblImportSuccess != 1) {
      importFailedGlobally = true;
      if (myRank == 0) {
        out << "Import FAILED: It threw an exception on at least one process." << endl;
      }
      gathervPrint (out, exceptionMessage, *comm);
    }

    lclImportSuccess = (A1.localError () || A2.localError ()) ? 0 : 1;
    gblImportSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclImportSuccess, outArg (gblImportSuccess));
    TEST_ASSERT( gblImportSuccess == 1 );
    if (gblImportSuccess != 1) {
      importFailedGlobally = true;
      if (myRank == 0) {
        out << "Import FAILED by reporting local error on at least one process." << endl
            << "Error messages from A1:" << endl;
      }
      gathervPrint (out, A1.errorMessages (), *comm);
      if (myRank == 0) {
        out << "Error messages from A2:" << endl;
      }
      gathervPrint (out, A2.errorMessages (), *comm);
    }

    if (importFailedGlobally) {
      return;
    }

    // doImport claims that it succeeded
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
    Kokkos::fence ();

    const LO myMinLclMeshRow = Y.getMap ()->getMinLocalIndex ();
    const LO myMaxLclMeshRow = Y.getMap ()->getMaxLocalIndex ();
    bool valsMatch = true;
    for (LO lclMeshRow = myMinLclMeshRow; lclMeshRow <= myMaxLclMeshRow; ++lclMeshRow) {
      auto Y_lcl = Y.getLocalBlockHost (lclMeshRow, 0, Tpetra::Access::ReadOnly);
      for (LO i = 0; i < blockSize; ++i) {
        if (static_cast<Scalar> (Y_lcl(i)) != requiredValue) {
          valsMatch = false;
        }
      }
    }
    TEST_ASSERT( valsMatch );

    //out << "The matrix A2, after Import (should be same as A1):" << endl;
    //A2.describe (out, Teuchos::VERB_EXTREME);

    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
  }

  // Test that two matrices' rows have the same entries.
  template<class BlockCrsMatrixType>
  bool matrices_are_same(const RCP<BlockCrsMatrixType>& A1,
                         const RCP<BlockCrsMatrixType>& A2)
  {
    // Loop through A1 and make sure each row has the same
    // entries as A2.  In the fully general case, the
    // redistribution may have added together values, resulting in
    // small rounding errors.  This is why we use an error tolerance
    // (with a little bit of wiggle room).

    int my_rank = A1->getRowMap()->getComm()->getRank();

    using Scalar = typename BlockCrsMatrixType::scalar_type;
    using LO = typename BlockCrsMatrixType::local_ordinal_type;
    using GO = typename BlockCrsMatrixType::global_ordinal_type;
    using lids_type = typename BlockCrsMatrixType::local_inds_host_view_type;
    using vals_type = typename BlockCrsMatrixType::values_host_view_type;

    using ST = ScalarTraits<Scalar>;
    using magnitude_type = typename ST::magnitudeType;
    const magnitude_type tol =
       Teuchos::as<magnitude_type> (10) * ScalarTraits<magnitude_type>::eps ();

    const LO blocksize = A1->getBlockSize();
    // Verify the blocksizes are identical
    if (blocksize != A2->getBlockSize()) {
      if (my_rank==0) std::cerr << "Error: Blocksizes are not the same!" << std::endl;
      return false;
    }

    lids_type A1LocalColInds;
    vals_type A1LocalRowVals;
    lids_type A2LocalColInds;
    vals_type A2LocalRowVals;
    for (LO localrow = A1->getRowMap()->getMinLocalIndex();
        localrow <= A1->getRowMap()->getMaxLocalIndex();
        ++localrow)
    {
      size_t A1NumEntries = A1->getNumEntriesInLocalRow (localrow);
      size_t A2NumEntries = A1->getNumEntriesInLocalRow (localrow);

      // Verify the same number of entries in each row
      if (A1NumEntries != A2NumEntries) {
	std::cerr << "Error: Matrices have different number of entries in at least one row!" << std::endl;
        return false;
      }

      A1->getLocalRowView (localrow, A1LocalColInds, A1LocalRowVals);
      A2->getLocalRowView (localrow, A2LocalColInds, A2LocalRowVals);

      // Verify the same number of values in each row
      if (A1LocalRowVals.extent(0) != A2LocalRowVals.extent(0)) {
	std::cerr << "Error: Matrices have different number of entries in at least one row!" << std::endl;
        return false;
      }

      // There's no guarantee the matrices have the same col map, so we compare global indices in
      // sets (because indices may be in a different order
      std::set<GO> a1_inds;
      std::set<GO> a2_inds;
      typedef typename Array<Scalar>::size_type size_type;
      for (size_type k = 0; k < static_cast<size_type> (A1NumEntries); ++k) {
        a1_inds.insert(A1->getColMap()->getGlobalElement(A1LocalColInds[k]));
        a2_inds.insert(A2->getColMap()->getGlobalElement(A2LocalColInds[k]));
      }
      if(a1_inds!=a2_inds) {
	std::cerr << "["<<localrow<<" ] Error: Matrices have different column indices!" << std::endl;
	std::cerr<< "A1_inds :";
	std::for_each(a1_inds.cbegin(), a1_inds.cend(), [&](int x) {
	    std::cerr << x << "("<<A1->getColMap()->getGlobalElement(x)<<") ";
	  });
	std::cerr<<std::endl;
	std::cerr<< "A2_inds :";
	std::for_each(a2_inds.cbegin(), a2_inds.cend(), [&](int x) {
	    std::cerr << x << "("<<A2->getColMap()->getGlobalElement(x)<<") ";
	  });
	std::cerr<<std::endl;


	return false;
      }

      // Loop over each local col entry of A1, find the corresponding col index of A2, and compare these value.
      LO INVALID = Teuchos::OrdinalTraits<LO>::invalid();
      for (size_type a1_k = 0; a1_k < static_cast<size_type> (A1NumEntries); ++a1_k) {
        LO a2_k = INVALID;
        for (size_type i = 0; i < static_cast<size_type> (A2NumEntries); ++i) {
	  if (A1->getColMap()->getGlobalElement(A1LocalColInds[a1_k]) == A2->getColMap()->getGlobalElement(A2LocalColInds[i])) {
            a2_k = i;
	    break;
	  }
        }
      
	if(a2_k == INVALID) {
	  std::cerr << "Cannot find corresponding column" << std::endl;
	  return false;
	}

        const int a1_start = a1_k*blocksize*blocksize;
        const int a2_start = a2_k*blocksize*blocksize;
        for (int b=0; b<blocksize*blocksize; ++b) {
          const magnitude_type rel_err = ST::magnitude(A1LocalRowVals[a1_start+b] - A2LocalRowVals[a2_start+b]);
          if(rel_err > tol) {
	    std::cerr << "Error: Matrices have different values!" << std::endl;
            return false;
          }
        }
      }
    }

    return true;
  }

  // Build lower diag matrix for test
  template<class BlockCrsMatrixType>
  void build_lower_diag_matrix (const RCP<BlockCrsMatrixType>& A) {

    using LO = typename BlockCrsMatrixType::local_ordinal_type;
    using GO = typename BlockCrsMatrixType::global_ordinal_type;
    using Scalar = typename BlockCrsMatrixType::scalar_type;

    const typename BlockCrsMatrixType::map_type row_map = *(A->getRowMap());
    const typename BlockCrsMatrixType::map_type col_map = *(A->getColMap());

    int my_rank = row_map.getComm()->getRank();

    if(A->getBlockSize() != 3) {
      if (my_rank==0) std::cerr << "Error: A->getBlockSize != 3!" << std::endl;
      return;
    }
    const int blocksize = 3;

    for (LO localrow = row_map.getMinLocalIndex();
         localrow <= row_map.getMaxLocalIndex();
         ++localrow) {

      const GO globalrow = row_map.getGlobalElement(localrow);

      if (globalrow == 0) {

        LO local_col_indices[1];
        local_col_indices[0] = col_map.getLocalElement(0);

        Scalar values[blocksize*blocksize];
        for (size_t b=0; b<blocksize*blocksize; ++b) {
          values[b] = 10*(globalrow+1);
        }
        A->replaceLocalValues(localrow,
                              local_col_indices,
                              values,
                              1);
      }
      else if (globalrow == 1) {

        LO local_col_indices[2];
        local_col_indices[0] = col_map.getLocalElement(0);
        local_col_indices[1] = col_map.getLocalElement(1);

        Scalar values[2*blocksize*blocksize];
        for (GO globalcol=0; globalcol<2; ++globalcol) {
          int start = globalcol*blocksize*blocksize;
          for (size_t b=0; b<blocksize*blocksize; ++b) {
            values[start+b] = 10*(globalrow+1)+globalcol;
          }
        }
        A->replaceLocalValues(localrow,
                              local_col_indices,
                              values,
                              2);
      } else {

        LO local_col_indices[3];
        local_col_indices[0] = col_map.getLocalElement(globalrow-2);
        local_col_indices[1] = col_map.getLocalElement(globalrow-1);
        local_col_indices[2] = col_map.getLocalElement(globalrow);

        Scalar values[3*blocksize*blocksize];
        int local_indx = 0;
        for (GO globalcol=globalrow-2; globalcol<=globalrow; ++globalcol) {
          int start = local_indx*blocksize*blocksize;
          for (size_t b=0; b<blocksize*blocksize; ++b) {
            values[start+b] = 10*(globalrow+1)+globalcol;
          }
          ++local_indx;
        }
        A->replaceLocalValues(localrow,
                              local_col_indices,
                              values,
                              3);
      }
    }

    return;
  }

  // Test BlockCrsMatrix importAndFillComplete
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockCrsMatrix, importAndFillComplete, Scalar, LO, GO, Node )
  {
    using Tpetra::Details::gathervPrint;
    typedef Tpetra::BlockCrsMatrix<Scalar, LO, GO, Node> block_crs_type;
    typedef Tpetra::CrsGraph<LO, GO, Node> crs_graph_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::Import<LO, GO, Node> import_type;
    using Teuchos::REDUCE_MAX;

    std::ostringstream err;
    int lclErr = 0;
    int gblErr = 0;

    out << "Testing Tpetra::BlockCrsMatrix importAndFillComplete" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const Comm<int> > comm = getDefaultComm ();
    const int myRank = comm->getRank ();
    const int numRanks = comm->getSize();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    out << "1st test: Import a diagonal BlockCrsMatrix from a source row Map "
           "that has all indices on Process 0, to a target row Map that is "
           "uniformly distributed over processes. Blocksize=3." << endl;
    try {
      Teuchos::OSTab tab1 (out);

      const GO indexBase = 0;
      const LO tgt_num_local_elements = 2;
      const LO src_num_local_elements = (myRank == 0) ?
        static_cast<LO> (numRanks*tgt_num_local_elements) :
        static_cast<LO> (0);

      const int blocksize = 3;

      // Create row Maps for the source and target
      RCP<const map_type> src_map =
        rcp (new map_type (INVALID,
                           src_num_local_elements,
                           indexBase, comm));
      RCP<const map_type> tgt_map =
        rcp (new map_type (INVALID,
                           tgt_num_local_elements,
                           indexBase, comm));

      // Build src graph.
      Teuchos::RCP<crs_graph_type> src_graph =
        Teuchos::rcp (new crs_graph_type (src_map, 1));
      for (LO localrow = src_map->getMinLocalIndex();
           localrow<=src_map->getMaxLocalIndex(); 
           ++localrow) {

        const GO globalrow = src_map->getGlobalElement(localrow);
        GO globalcol[1];
        globalcol[0] = globalrow;
        
        src_graph->insertGlobalIndices(globalrow, 1, globalcol);
      }
      src_graph->fillComplete();

      // Build src matrix. Simple block diagonal matrix with A(b,b) = [b*b*row,...,+b*b].
      RCP<block_crs_type> src_mat =
        rcp (new block_crs_type (*src_graph, blocksize)); 
      if (src_num_local_elements != 0) {
        for (LO localrow = src_map->getMinLocalIndex();
             localrow <= src_map->getMaxLocalIndex();
             ++localrow) {
          const GO globalrow = src_map->getGlobalElement(localrow);
          LO col_indices[1];  Scalar values[blocksize*blocksize];
          col_indices[0] = localrow; 
          for (size_t b=0; b<blocksize*blocksize; ++b) {
            values[b] = blocksize*blocksize*globalrow + b;
          }
          const LO actual_num_replaces = src_mat->replaceLocalValues(localrow,
                                                                     col_indices,
                                                                     values,
                                                                     1);
          TEST_EQUALITY_CONST(actual_num_replaces, 1);
        }
      }

      // Create the importer
      import_type importer (src_map, tgt_map);

      // Call importAndFillComplete to get the tgt matrix
      RCP<block_crs_type> tgt_mat =
        Tpetra::importAndFillCompleteBlockCrsMatrix<block_crs_type> (src_mat, importer);
     
      // Manually build the tgt matrix and test that it matches the returned matrix

      // Build tgt graph.
      Teuchos::RCP<crs_graph_type> tgt_graph_for_testing =
        Teuchos::rcp (new crs_graph_type (tgt_map, 1));
      for (LO localrow = tgt_map->getMinLocalIndex();
           localrow<=tgt_map->getMaxLocalIndex();
           ++localrow) {

        const GO globalrow = tgt_map->getGlobalElement(localrow);
        GO globalcol[1];
        globalcol[0] = globalrow;

        tgt_graph_for_testing->insertGlobalIndices(globalrow, 1, globalcol);
      }
      tgt_graph_for_testing->fillComplete();

      // Build tgt matrix
      RCP<block_crs_type> tgt_mat_for_testing =
        rcp (new block_crs_type (*tgt_graph_for_testing, blocksize));
      for (LO localrow = tgt_map->getMinLocalIndex();
           localrow <= tgt_map->getMaxLocalIndex();
           ++localrow) {
        const GO globalrow = tgt_map->getGlobalElement(localrow);
        LO col_indices[1];  Scalar values[blocksize*blocksize];
        col_indices[0] = localrow;
        for (size_t b=0; b<blocksize*blocksize; ++b) {
          values[b] = blocksize*blocksize*globalrow + b;
        }
        const LO actual_num_replaces = tgt_mat_for_testing->replaceLocalValues(localrow,
                                                                               col_indices,
                                                                               values,
                                                                               1);
        TEST_EQUALITY_CONST(actual_num_replaces, 1);
      }

      // Test that matrices are identical
      bool matrices_match = matrices_are_same<block_crs_type>(tgt_mat, tgt_mat_for_testing);
      TEST_ASSERT(matrices_match);
     }
     catch (std::exception& e) { // end of the first test
       err << "Proc " << myRank << ": " << e.what () << endl;
       lclErr = 1;
     }

     reduceAll<int, int> (*comm, REDUCE_MAX, lclErr, outArg (gblErr));
     TEST_EQUALITY_CONST( gblErr, 0 );
     if (gblErr != 0) {
       Tpetra::Details::gathervPrint (out, err.str (), *comm);
       out << "Above test failed; aborting further tests" << endl;
       return;
     }

     out << "2nd test: Import a lower triangular BlockCrsMatrix from a source row Map "
            "where even processors have 1 element and odd processors have 3 elements, "
            "to a target row Map where each processor have 2 elements. Blocksize=3." << endl;
     try {
       Teuchos::OSTab tab1 (out);

       // This test only makes sense for even number of ranks
       if (numRanks % 2 != 0) {
         return;
       }

       const GO indexBase = 0;
       LO src_num_local_elements;
       if (myRank % 2 == 0) src_num_local_elements = 1;
       else                 src_num_local_elements = 3;
       LO tgt_num_local_elements = 2;
       const int blocksize = 3;

       // Create row Maps for the source and target
       RCP<const map_type> src_map =
         rcp (new map_type (INVALID,
                            src_num_local_elements,
                            indexBase, comm));
       RCP<const map_type> tgt_map =
         rcp (new map_type (INVALID,
                            tgt_num_local_elements,
                            indexBase, comm));
       //src_map->describe(out, Teuchos::VERB_EXTREME);
       //tgt_map->describe(out, Teuchos::VERB_EXTREME);

       // Build src graph. Allow for up to 2 off-diagonal entries.
       Teuchos::RCP<crs_graph_type> src_graph =
         Teuchos::rcp (new crs_graph_type (src_map, 3));
       {
         Array<GO> cols(3);
         for (GO globalrow = src_map->getMinGlobalIndex ();
              globalrow <= src_map->getMaxGlobalIndex (); ++globalrow) {
           if      (globalrow==0) cols.resize(1);
           else if (globalrow==1) cols.resize(2);
           else                   cols.resize(3);
           for (GO col = 0; col < cols.size(); ++col) {
             cols[col] = globalrow - col;
           }
           src_graph->insertGlobalIndices (globalrow, cols());
         }
         src_graph->fillComplete();
         //src_graph->describe(out, Teuchos::VERB_EXTREME);
       }

       // Build src matrix. Simple block lower-diagonal matrix with
       // A(b1,b2) = [(b1)+10*(b2+1)].
       RCP<block_crs_type> src_mat =
         rcp (new block_crs_type (*src_graph, blocksize));
       build_lower_diag_matrix<block_crs_type>(src_mat);
       //src_mat->describe(out, Teuchos::VERB_EXTREME);

       // Create the importer
       import_type importer (src_map, tgt_map);

       // Call importAndFillComplete to get the tgt matrix
       RCP<block_crs_type> tgt_mat =
         Tpetra::importAndFillCompleteBlockCrsMatrix<block_crs_type> (src_mat, importer);
       //tgt_mat->describe(out, Teuchos::VERB_EXTREME);

       // Manually build the tgt matrix and test that it matches the returned matrix

       // Build tgt graph.
       Teuchos::RCP<crs_graph_type> tgt_graph_for_testing =
         Teuchos::rcp (new crs_graph_type (tgt_map, 3));
       {
         Array<GO> cols(3);
         for (GO globalrow = tgt_map->getMinGlobalIndex ();
              globalrow <= tgt_map->getMaxGlobalIndex (); ++globalrow) {
           if      (globalrow==0) cols.resize(1);
           else if (globalrow==1) cols.resize(2);
           else                   cols.resize(3);
           for (GO col = 0; col < cols.size(); ++col) {
             cols[col] = globalrow - col;
           }
           tgt_graph_for_testing->insertGlobalIndices (globalrow, cols());
         }
         tgt_graph_for_testing->fillComplete();
         //tgt_graph_for_testing->describe(out, Teuchos::VERB_EXTREME);
       }

       // Build tgt matrix
       RCP<block_crs_type> tgt_mat_for_testing =
         rcp (new block_crs_type (*tgt_graph_for_testing, blocksize));
       build_lower_diag_matrix<block_crs_type>(tgt_mat_for_testing);



       tgt_mat_for_testing->describe(out, Teuchos::VERB_EXTREME);
       tgt_mat->describe(out, Teuchos::VERB_EXTREME);

       


       // Test that matrices are identical
       bool matrices_match = matrices_are_same<block_crs_type>(tgt_mat, tgt_mat_for_testing);
       TEST_ASSERT(matrices_match);
      }
      catch (std::exception& e) { // end of the first test
        err << "Proc " << myRank << ": " << e.what () << endl;
        lclErr = 1;
      }

      reduceAll<int, int> (*comm, REDUCE_MAX, lclErr, outArg (gblErr));
      TEST_EQUALITY_CONST( gblErr, 0 );
      if (gblErr != 0) {
        Tpetra::Details::gathervPrint (out, err.str (), *comm);
        out << "Above test failed; aborting further tests" << endl;
        return;
      }
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
    graph_type overlapGraph (overlapMeshRowMapPtr, 2);
    if (overlapMeshRowMapPtr->getLocalNumElements () > 0) {
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
      const LO lclNumRowsOverlap = static_cast<LO> (overlapMeshRowMapPtr->getLocalNumElements ());
      for (LO lclRow = 0; lclRow < lclNumRowsOverlap; ++lclRow) {
        const LO numEnt = static_cast<LO> (overlapGraph.getNumEntriesInLocalRow (lclRow));
        TEST_EQUALITY( numEnt, static_cast<LO> (2) );
      }
      const LO maxNumRowEnt = static_cast<LO> (overlapGraph.getLocalMaxNumRowEntries ());
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
      graph = rcp (new graph_type (meshRowMapPtr, 0));

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
        static_cast<LO> (meshRowMapPtr->getLocalNumElements ());
      for (LO lclRow = 0; lclRow < lclNumRowsNonoverlap; ++lclRow) {
        const LO numEnt = static_cast<LO> (graph->getNumEntriesInLocalRow (lclRow));
        TEST_ASSERT( numEnt >= static_cast<LO> (2) );
      }
      const LO maxNumRowEnt = static_cast<LO> (graph->getLocalMaxNumRowEntries ());
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
      //   auto Y_lcl = Y.getLocalBlockHost(lclMeshRow, 0,Tpetra::Access::ReadOnly);
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
      blockMatrix = Tpetra::convertToBlockCrsMatrix(*pointMatrix,blockSize);
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

    const magnitude_type tol =
      magnitude_type(10.0) * Teuchos::ScalarTraits<magnitude_type>::eps();
    TEUCHOS_TEST_FOR_EXCEPTION(relativeError[0] > tol, std::runtime_error, "BlockCrsMatrix matvec does not produce same result as CrsMatrix matvec.");
  }




  //
  // Test conversion from BlockCrsMatrix to CrsMatrix.  This test is valid only on
  // 1, 2, or 4 processes.  (See note below.)
  //
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockCrsMatrix, block2point, ST, LO, GO, Node )
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
      blockMatrix = Tpetra::convertToBlockCrsMatrix(*pointMatrix,blockSize);
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


    // Now convert back

    RCP<crs_matrix_type> newPointMatrix;
    newPointMatrix = Tpetra::convertToCrsMatrix(*blockMatrix);


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
      std::cerr << "NewPointMatrix::apply" << endl;
    }
    else {
      out << "NewPointMatrix::apply" << endl;
    }
    try {
      newPointMatrix->apply(*randVec, *resultVec2, Teuchos::NO_TRANS, one, zero);
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



  }



  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockCrsMatrix, Transpose, ST, LO, GO, Node )
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

    Teuchos::OSTab tab0 (out);
    out << "Testing Transpose"<<std::endl;
    Teuchos::OSTab tab1 (out);

    std::ostringstream errStrm;

    RCP<const Comm<int> > comm = getDefaultComm();
    std::string matrixFile;
    if (STS::isComplex) {
      matrixFile = "blockA-complex.mm";
    }
    else {
      matrixFile = "blockA.mm";
    }

    out << "Read CrsMatrix from file \"" << matrixFile << "\"" << endl;
    
    RCP<crs_matrix_type> pointMatrix = reader_type::readSparseFile(matrixFile, comm);

    out << "Migrate input CrsMatrix to final parallel distribution" << endl;

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

    parPointMatrix->doImport(*pointMatrix, *importer, Tpetra::INSERT);
    parPointMatrix->fillComplete();
    pointMatrix.swap(parPointMatrix);

    
    out << "Convert CrsMatrix to BlockCrsMatrix" << endl;
    int blockSize = 3;
    RCP<block_matrix_type> blockMatrix = Tpetra::convertToBlockCrsMatrix(*pointMatrix,blockSize);

    // Now Transpose them both
    Tpetra::RowMatrixTransposer<ST, LO, GO, Node> Pt(pointMatrix);
    RCP<crs_matrix_type> newPointMatrix = Pt.createTranspose();


    Tpetra::BlockCrsMatrixTransposer<ST, LO, GO, Node> Bt(blockMatrix);
    RCP<block_matrix_type> newBlockMatrix = Bt.createTranspose();


    // Check the maps
    TEUCHOS_TEST_FOR_EXCEPTION(!blockMatrix->getDomainMap()->isSameAs(*newBlockMatrix->getRangeMap()),std::runtime_error,"Incorrect range map on transpose");
    TEUCHOS_TEST_FOR_EXCEPTION(!blockMatrix->getRangeMap()->isSameAs(*newBlockMatrix->getDomainMap()),std::runtime_error,"Incorrect domain map on transpose");

    //normalized pseudo-random vector
    RCP<mv_type> randVec = rcp(new mv_type(newPointMatrix->getDomainMap(),1));
    randVec->randomize();
    Teuchos::Array<magnitude_type> normVec1(1);
    randVec->norm2(normVec1);
    randVec->scale(1.0/normVec1[0]);

    RCP<mv_type> resultVec1 = rcp(new mv_type(newPointMatrix->getRangeMap(),1));
    out << "CrsMatrix::apply" << endl;
    newPointMatrix->apply(*randVec, *resultVec1, Teuchos::NO_TRANS, one, zero);
    
    out << "Compute norm of result" << endl;
    resultVec1->norm2(normVec1);

    RCP<mv_type> resultVec2 = rcp(new mv_type(newBlockMatrix->getRangeMap(),1));
    out << "NewBlockMatrix::apply" << endl;
    newBlockMatrix->apply(*randVec, *resultVec2, Teuchos::NO_TRANS, one, zero);

    Teuchos::Array<magnitude_type> normVec2(1);
    resultVec2->norm2(normVec2);

    resultVec2->update(-1.0,*resultVec1,1.0);
    Teuchos::Array<magnitude_type> normDelta(1);
    resultVec2->norm2(normDelta);
    Teuchos::Array<magnitude_type> relativeError(1);
    relativeError[0] = STS::magnitude(normDelta[0] / normVec1[0]);

    std::ostringstream normStr;
    normStr << "||CSR*xrand|| = " << normVec1[0] << ", ||CSR*xrand - BCSR*xrand|| / ||CSR*xrand|| = " << relativeError[0];
    out << normStr.str() << std::endl;

    // Correctness check
    auto tol  = Teuchos::ScalarTraits<magnitude_type>::squareroot(Teuchos::ScalarTraits<magnitude_type>::eps());
    TEUCHOS_TEST_FLOATING_EQUALITY(normVec1[0], normVec2[0], tol, out, success);

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
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockCrsMatrix, importAndFillComplete, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockCrsMatrix, ExportDiffRowMaps, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockCrsMatrix, point2block, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockCrsMatrix, block2point, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockCrsMatrix, Transpose, SCALAR, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  // NOTE (mfh 24 Sep 2015) It only makes sense to test over Scalar
  // types which have a Teuchos::LAPACK implementation.

  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP )

} // namespace (anonymous)


