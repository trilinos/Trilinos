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
#include <Tpetra_ConfigDefs.hpp>

#include <Tpetra_TestingUtilities.hpp>

#include <Tpetra_Experimental_BlockMultiVector.hpp>

namespace {

  //
  // UNIT TESTS
  //

  // Test BlockMultiVector's constructors.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( ExpBlockMultiVector, ctor, Scalar, LO, GO, Node )
  {
    using Tpetra::TestingUtilities::getNode;
    using Tpetra::TestingUtilities::getDefaultComm;
    using Teuchos::Comm;
    using Teuchos::RCP;
    typedef Tpetra::Experimental::BlockMultiVector<Scalar, LO, GO, Node> BMV;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::global_size_t GST;

    RCP<const Comm<int> > comm = getDefaultComm ();
    RCP<Node> node = getNode<Node> ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    const size_t numLocalMeshPoints = 12;
    const GO indexBase = 0;
    map_type meshMap (INVALID, numLocalMeshPoints, indexBase, comm, node);
    //RCP<const map_type> mapPtr = Teuchos::rcpFromRef (map); // nonowning RCP

    const LO blockSize = 4;
    const LO numVecs = 3;

    // Exercise the three-argument constructor (that makes a point Map).
    BMV X (meshMap, blockSize, numVecs);
    TEST_EQUALITY( X.getBlockSize (), blockSize );
    TEST_EQUALITY( X.getNumVectors (), numVecs );
    TEST_ASSERT( ! X.getMap ().is_null () );
    TEST_ASSERT( ! X.getMap ().is_null () && X.getMap ()->isSameAs (meshMap) );

    map_type pointMap = BMV::makePointMap (meshMap, blockSize);
    TEST_ASSERT( pointMap.isSameAs (X.getPointMap ()) );

    // Exercise the four-argument constructor (that uses an existing
    // point Map).
    BMV Y (meshMap, pointMap, blockSize, numVecs);
    TEST_EQUALITY( Y.getBlockSize (), blockSize );
    TEST_EQUALITY( Y.getNumVectors (), numVecs );
    TEST_ASSERT( ! Y.getMap ().is_null () );
    TEST_ASSERT( ! Y.getMap ().is_null () && Y.getMap ()->isSameAs (meshMap) );

    TEST_ASSERT( Y.getPointMap ().isSameAs (X.getPointMap ()) );

    // Exercise the default constructor.
    BMV Z;
    TEST_ASSERT( Z.getMap ().is_null () );
    TEST_EQUALITY_CONST( Z.getBlockSize (), static_cast<LO> (0) );
    TEST_EQUALITY_CONST( Z.getNumVectors (), static_cast<LO> (0) );

    // Exercise the other three-argument constructor (that views an
    // existing MultiVector).
    typename BMV::mv_type X_mv;
    TEST_NOTHROW( X_mv = X.getMultiVectorView () );
    TEST_ASSERT( ! X_mv.getMap ().is_null () );
    if (X_mv.getMap ().is_null ()) {
      TEST_ASSERT( ! X_mv.getMap ()->getNode ().is_null () );
    }
    // Make sure X_mv has view semantics, before we give it to W's ctor.
    TEST_ASSERT( X_mv.getCopyOrView () == Teuchos::View );
    BMV W (X_mv, meshMap, blockSize);
    // TEST_ASSERT( ! W.getMap ().is_null () );
    // TEST_ASSERT( ! W.getMap ().is_null () && W.getMap ()->isSameAs (meshMap) );
    // TEST_ASSERT( W.getPointMap ().isSameAs (X.getPointMap ()) );
  }

  // Test BlockMultiVector::getMultiVectorView.  It must return a
  // MultiVector with view semantics, and it must view the same data
  // that the BlockMultiVector view.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( ExpBlockMultiVector, MVView, Scalar, LO, GO, Node )
  {
    using Tpetra::TestingUtilities::getNode;
    using Tpetra::TestingUtilities::getDefaultComm;
    using Teuchos::Comm;
    using Teuchos::RCP;
    typedef Tpetra::Experimental::BlockMultiVector<Scalar, LO, GO, Node> BMV;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::global_size_t GST;
    typedef Teuchos::ScalarTraits<Scalar> STS;

    RCP<const Comm<int> > comm = getDefaultComm ();
    RCP<Node> node = getNode<Node> ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    const size_t numLocalMeshPoints = 12;
    const LO blockSize = 4;
    const LO numVecs = 3;

    const GO indexBase = 0;
    map_type meshMap (INVALID, numLocalMeshPoints, indexBase, comm, node);
    //RCP<const map_type> mapPtr = Teuchos::rcpFromRef (map); // nonowning RCP

    BMV X (meshMap, blockSize, numVecs);

    typedef typename BMV::mv_type mv_type;
    mv_type X_mv = X.getMultiVectorView ();
    TEST_EQUALITY_CONST( X_mv.getCopyOrView (), Teuchos::View );

    TEST_ASSERT( ! X_mv.getMap ().is_null () );
    // Keep the "not null" test in subsequent tests, so that we don't
    // get segfaults if getMap() really does return null.
    TEST_ASSERT( ! X_mv.getMap ().is_null () &&
                 X_mv.getMap ()->isSameAs (X.getPointMap ()) );
    map_type mv_map;
    if (! X_mv.getMap ().is_null ()) {
      mv_map = * (X_mv.getMap ()); // Map implements view semantics
    }
    TEST_EQUALITY_CONST( mv_map.getNodeNumElements (), static_cast<size_t> (48) );

    // X has three columns.  In the following tests, change only
    // blocks in column 1.  Entries in all other columns must remain
    // unmodified.
    const LO colToModify = 1;

    // Test that X_mv views the same data that the BMV views.  X_mv
    // has 48 rows on each process: 12 mesh points, and 4 degrees of
    // freedom per mesh point ("block size").  Rows 20-23 thus
    // correspond to local mesh point 5.
    typedef typename BMV::little_vec_type little_vec_type;
    little_vec_type X_5_1 = X.getLocalBlock (5, colToModify);

    // All entries of X_5_1 must be zero.  First make a block with all
    // zero entries, then test.  It's not worth testing the
    // corresponding entries of X_mv yet; we'll do that below.
    Teuchos::Array<Scalar> zeroArray (blockSize, STS::zero ());
    little_vec_type zeroLittleVector (zeroArray.getRawPtr (), blockSize, 1);
    TEST_ASSERT( X_5_1.equal (zeroLittleVector) && zeroLittleVector.equal (X_5_1) );

    // Put some data in the block.  This will help us test whether the
    // corresponding entries in the MultiVector, and only those
    // entries, got modified.
    for (LO i = 0; i < blockSize; ++i) {
      X_5_1(i) = static_cast<Scalar> (i + 1); // all are nonzero
    }
    TEST_ASSERT( ! X_5_1.equal (zeroLittleVector) && ! zeroLittleVector.equal (X_5_1) );

    // Make sure that getLocalBlock() returns a read-and-write view,
    // not a deep copy.  Do this by calling getLocalBlock(5,1) again,
    // and testing that changes to X_5_1 are reflected in the result.
    little_vec_type X_5_1_new = X.getLocalBlock (5, colToModify);
    TEST_ASSERT( X_5_1_new.equal (X_5_1) && X_5_1.equal (X_5_1_new) );
    TEST_ASSERT( ! X_5_1_new.equal (zeroLittleVector) && ! zeroLittleVector.equal (X_5_1_new) );

    // Make sure that all blocks other than block 5 still contain zeros.
    for (LO localMeshIndex = 0; localMeshIndex < static_cast<LO> (numLocalMeshPoints); ++localMeshIndex) {
      for (LO curCol = 0; curCol < numVecs; ++curCol) {
        little_vec_type X_cur = X.getLocalBlock (localMeshIndex, curCol);
        if (curCol != colToModify) {
          TEST_ASSERT( X_cur.equal (zeroLittleVector) && zeroLittleVector.equal (X_cur) );
          TEST_ASSERT( ! X_cur.equal (X_5_1) && ! X_5_1.equal (X_cur) );
        }
        if (localMeshIndex != 5) {
          TEST_ASSERT( X_cur.equal (zeroLittleVector) && zeroLittleVector.equal (X_cur) );
          TEST_ASSERT( ! X_cur.equal (X_5_1) && ! X_5_1.equal (X_cur) );
        }
      }
    }

    // Make sure that all entries of X_mv (the MultiVector) contain
    // zeros, except for rows 20-23 in column colToModify.  Those must
    // contain 1,2,3,4 (in order).
    const size_t numLocalDofs = mv_map.getNodeNumElements ();
    for (LO curCol = 0; curCol < numVecs; ++curCol) {
      Teuchos::ArrayRCP<const Scalar> X_mv_curCol =
        X_mv.getData (static_cast<size_t> (curCol));
      if (curCol == static_cast<size_t> (colToModify)) {
        for (size_t localDofIndex = 0; localDofIndex < numLocalDofs; ++localDofIndex) {
          if (localDofIndex == 20) {
            TEST_ASSERT( X_mv_curCol[20] == X_5_1(0) );
            TEST_ASSERT( X_mv_curCol[20] == static_cast<Scalar> (1) );
          }
          else if (localDofIndex == 21) {
            TEST_ASSERT( X_mv_curCol[21] == X_5_1(1) );
            TEST_ASSERT( X_mv_curCol[21] == static_cast<Scalar> (2) );
          }
          else if (localDofIndex == 22) {
            TEST_ASSERT( X_mv_curCol[22] == X_5_1(2) );
            TEST_ASSERT( X_mv_curCol[22] == static_cast<Scalar> (3) );
          }
          else if (localDofIndex == 23) {
            TEST_ASSERT( X_mv_curCol[23] == X_5_1(3) );
            TEST_ASSERT( X_mv_curCol[23] == static_cast<Scalar> (4) );
          }
          else {
            TEST_ASSERT( X_mv_curCol[localDofIndex] == STS::zero () );
          }
        }
      }
      else {
        for (size_t localDofIndex = 0; localDofIndex < numLocalDofs; ++localDofIndex) {
          TEST_ASSERT( X_mv_curCol[localDofIndex] == STS::zero () );
        }
      }
    }
  }

  // Make sure that Import works with BlockMultiVector.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( ExpBlockMultiVector, Import, Scalar, LO, GO, Node )
  {
    using Tpetra::TestingUtilities::getNode;
    using Tpetra::TestingUtilities::getDefaultComm;
    using Teuchos::Comm;
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using Teuchos::RCP;
    typedef Tpetra::Experimental::BlockMultiVector<Scalar, LO, GO, Node> BMV;
    typedef typename BMV::little_vec_type little_vec_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::Import<LO, GO, Node> import_type;
    typedef Tpetra::global_size_t GST;
    typedef Teuchos::ScalarTraits<Scalar> STS;

    RCP<const Comm<int> > comm = getDefaultComm ();
    // const int myRank = comm->getRank ();
    // const int numProcs = comm->getSize ();
    RCP<Node> node = getNode<Node> ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    const size_t numLocalMeshPoints = 2;
    const LO blockSize = 5;
    const LO numVecs = 3;
    const GO indexBase = 0;
    map_type meshMap (INVALID, numLocalMeshPoints, indexBase, comm, node);
    //RCP<const map_type> mapPtr = Teuchos::rcpFromRef (map); // nonowning RCP

    Teuchos::Array<GO> overlappingGIDs (numLocalMeshPoints + 1);
    for (LO lid = 0; lid < static_cast<LO> (numLocalMeshPoints); ++lid) {
      overlappingGIDs[lid] = meshMap.getGlobalElement (lid);
    }
    // Every process gets exactly one overlapping GID.  In the
    // (nonoverlapping) mesh Map on the calling process, my min GID
    // overlaps with exactly one process ((myRank + 1) % numProcs) in
    // the overlapping mesh Map.
    overlappingGIDs[numLocalMeshPoints] =
      overlappingGIDs[numLocalMeshPoints-1] %
      static_cast<GO> (meshMap.getGlobalNumElements ());
    map_type overlappingMeshMap (INVALID, overlappingGIDs (), indexBase, comm, node);

    BMV X (meshMap, blockSize, numVecs);
    BMV Y (overlappingMeshMap, blockSize, numVecs);

    //
    // Fill X with meaningful things to test Import with REPLACE combine mode.
    //
    // Fill just the "outgoing" blocks of X, in column 1 only, with
    // 1,2,3,4,5.  Due to how we constructed the overlapping mesh Map,
    // the first GID on the calling process in the mesh Map overlaps
    // with one GID on exactly one process.
    const LO colToModify = 1;
    little_vec_type X_overlap =
      X.getLocalBlock (meshMap.getLocalElement (meshMap.getMinGlobalIndex ()), colToModify);
    TEST_ASSERT( X_overlap.getRawPtr () != NULL );
    TEST_EQUALITY_CONST( X_overlap.getBlockSize (), blockSize );
    // FIXME (mfh 07 May 2014) May have to fix this if strides change
    // in the future.
    TEST_EQUALITY_CONST( X_overlap.getStride (), static_cast<LO> (1) );

    // {
    //   std::ostringstream os;
    //   os << "Proc " << myRank
    //      << ": X_overlap.getRawPtr() = " << X_overlap.getRawPtr ()
    //      << ", X_overlap.getBlockSize() = " << X_overlap.getBlockSize ()
    //      << ", meshMap.getMinGlobalIndex() = " << meshMap.getMinGlobalIndex ()
    //      << std::endl;
    //   std::cerr << os.str ();
    // }

    {
      const int lclOk = (X_overlap.getRawPtr () != NULL &&
                         X_overlap.getBlockSize () == blockSize) ? 1 : 0;
      int gblOk = 1;
      reduceAll<int, int> (*comm, REDUCE_MIN, lclOk, outArg (gblOk));
      TEUCHOS_TEST_FOR_EXCEPTION(
        gblOk == 0, std::logic_error, "Some process reported that X_overlap "
        "is an empty block.");
    }

    for (LO i = 0; i < blockSize; ++i) {
      X_overlap(i) = static_cast<Scalar> (i+1);
    }

    import_type import (rcpFromRef (meshMap), rcpFromRef (overlappingMeshMap));
    Y.doImport (X, import, Tpetra::REPLACE);

    little_vec_type Y_overlap = Y.getLocalBlock (overlappingMeshMap.getLocalElement (overlappingMeshMap.getMinGlobalIndex ()), colToModify);

    // Forestall compiler warning for unused variable.
    (void) Y_overlap;

    Teuchos::Array<Scalar> zeroArray (blockSize, STS::zero ());
    little_vec_type zeroLittleVector (zeroArray.getRawPtr (), blockSize, 1);

    for (LO col = 0; col < numVecs; ++col) {
      for (LO localMeshRow = meshMap.getMinLocalIndex ();
           localMeshRow < meshMap.getMaxLocalIndex (); ++localMeshRow) {
        little_vec_type Y_cur = Y.getLocalBlock (localMeshRow, col);
        if (col != colToModify) {
          TEST_ASSERT( Y_cur.equal (zeroLittleVector) && zeroLittleVector.equal (Y_cur) );
          TEST_ASSERT( ! Y_cur.equal (X_overlap) && ! X_overlap.equal (Y_cur) );
        }
        if (localMeshRow != meshMap.getMinLocalIndex ()) {
          TEST_ASSERT( Y_cur.equal (zeroLittleVector) && zeroLittleVector.equal (Y_cur) );
          TEST_ASSERT( ! Y_cur.equal (X_overlap) && ! X_overlap.equal (Y_cur) );
        }
        if (col == colToModify && localMeshRow == meshMap.getMinLocalIndex ()) {
          TEST_ASSERT( Y_cur.equal (X_overlap) && X_overlap.equal (Y_cur) );
        }
      }
    }
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( ExpBlockMultiVector, ctor, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( ExpBlockMultiVector, MVView, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( ExpBlockMultiVector, Import, SCALAR, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV_NOGPU( UNIT_TEST_GROUP )

} // namespace (anonymous)


