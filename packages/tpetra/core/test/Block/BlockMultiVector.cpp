// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_BlockMultiVector.hpp"
#include "Tpetra_BlockVector.hpp"
#include "Tpetra_Access.hpp"
namespace {

  template<class ViewType1,
           class ViewType2,
           const int rank1 = ViewType1::rank,
           const int rank2 = ViewType2::rank>
  struct LittleEqual {
    static bool equal (const ViewType1& X, const ViewType2& Y) {
      return false; // default
    }
  };

  template<class ViewType1,
           class ViewType2>
  struct LittleEqual<ViewType1, ViewType2, 2, 2> {
    static bool equal (const ViewType1& X, const ViewType2& Y)
    {
      if (X.extent (0) != Y.extent (0) ||
          X.extent (1) != Y.extent (1)) {
        return false;
      }
      for (int j = 0; j < static_cast<int> (X.extent (1)); ++j) {
        for (int i = 0; i < static_cast<int> (X.extent (0)); ++i) {
          if (X(i,j) != Y(i,j)) {
            return false;
          }
        }
      }
      return true;
    }
  };

  template<class ViewType1,
           class ViewType2>
  struct LittleEqual<ViewType1, ViewType2, 1, 1> {
    static bool equal (const ViewType1& X, const ViewType2& Y)
    {
      if (X.extent (0) != Y.extent (0)) {
        return false;
      }
      for (int i = 0; i < static_cast<int> (X.extent (0)); ++i) {
        if (X(i) != Y(i)) {
          return false;
        }
      }
      return true;
    }
  };

  template<class ViewType1,
           class ViewType2,
           const int rank1 = ViewType1::rank,
           const int rank2 = ViewType2::rank>
  bool equal (const ViewType1& X, const ViewType2& Y) {
    return LittleEqual<ViewType1, ViewType2, rank1, rank2>::equal (X, Y);
  }

  //
  // UNIT TESTS
  //

  // Test BlockMultiVector's constructors.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockMultiVector, ctor, Scalar, LO, GO, Node )
  {
    using Teuchos::Comm;
    using Teuchos::RCP;
    typedef Tpetra::BlockMultiVector<Scalar, LO, GO, Node> BMV;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::global_size_t GST;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    const size_t numLocalMeshPoints = 12;
    const GO indexBase = 0;
    map_type meshMap (INVALID, numLocalMeshPoints, indexBase, comm);

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

    // Test BlockMultiVector's two-argument "copy constructor" that
    // can make either a deep or a shallow copy.
    {
      BMV X2 (X, Teuchos::Copy);
      auto X2_map = X2.getMap ();
      const bool maps_ok = ! X2_map.is_null () &&
        ! X.getMap ().is_null () &
        X2_map->isSameAs (* (X.getMap ()));
      TEST_ASSERT( maps_ok );

      auto X_mv_h = X.getMultiVectorView ().getLocalViewHost(Tpetra::Access::ReadOnly);
      auto X2_mv_h = X2.getMultiVectorView ().getLocalViewHost(Tpetra::Access::ReadOnly);
      TEST_ASSERT( X_mv_h.extent (0) == X2_mv_h.extent (0) &&
                   X_mv_h.extent (1) == X2_mv_h.extent (1) &&
                   X_mv_h.data () != X2_mv_h.data () );
    }

    {
      BMV X2 (X, Teuchos::View);
      auto X2_map = X2.getMap ();
      const bool maps_ok = ! X2_map.is_null () &&
        ! X.getMap ().is_null () &
        X2_map->isSameAs (* (X.getMap ()));
      TEST_ASSERT( maps_ok );

      auto X_mv_h = X.getMultiVectorView ().getLocalViewHost(Tpetra::Access::ReadOnly);
      auto X2_mv_h = X2.getMultiVectorView ().getLocalViewHost(Tpetra::Access::ReadOnly);
      TEST_ASSERT( X_mv_h.extent (0) == X2_mv_h.extent (0) &&
                   X_mv_h.extent (1) == X2_mv_h.extent (1) &&
                   X_mv_h.data () == X2_mv_h.data () );
    }

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
    // Make sure X_mv has view semantics, before we give it to W's ctor.
    TEST_ASSERT( X_mv.getCopyOrView () == Teuchos::View );
    BMV W (X_mv, meshMap, blockSize);
    // TEST_ASSERT( ! W.getMap ().is_null () );
    // TEST_ASSERT( ! W.getMap ().is_null () && W.getMap ()->isSameAs (meshMap) );
    // TEST_ASSERT( W.getPointMap ().isSameAs (X.getPointMap ()) );

    // Test BlockVector constructor
    typedef Tpetra::BlockVector<Scalar, LO, GO, Node> BV;
    typedef Tpetra::Vector<Scalar, LO, GO, Node> vec_type;
    Teuchos::RCP<const vec_type> V1 = X_mv.getVector (0);
    BV V (*V1, meshMap, blockSize);
    vec_type V2 = V.getVectorView ();
    TEST_EQUALITY( V1->getLocalViewHost(Tpetra::Access::ReadOnly).data (),
                   V2.getLocalViewHost(Tpetra::Access::ReadOnly).data () );
    TEST_EQUALITY( V1->getLocalViewDevice(Tpetra::Access::ReadOnly).data (),
                   V2.getLocalViewDevice(Tpetra::Access::ReadOnly).data () );

    // Test BlockVector's two-argument "copy constructor" that can
    // make either a deep or a shallow copy.
    {
      BV V_a (V, Teuchos::Copy);
      auto V_a_map = V_a.getMap ();
      const bool maps_ok = ! V_a_map.is_null () &&
        ! V.getMap ().is_null () &
        V_a_map->isSameAs (* (V.getMap ()));
      TEST_ASSERT( maps_ok );

      auto V_v_h = V.getVectorView ().getLocalViewHost(Tpetra::Access::ReadOnly);
      auto V_a_v_h = V_a.getVectorView ().getLocalViewHost(Tpetra::Access::ReadOnly);
      TEST_ASSERT( V_v_h.extent (0) == V_a_v_h.extent (0) &&
                   V_v_h.extent (1) == V_a_v_h.extent (1) &&
                   V_v_h.data () != V_a_v_h.data () );
    }

    {
      BV V_a (V, Teuchos::View);
      auto V_a_map = V_a.getMap ();
      const bool maps_ok = ! V_a_map.is_null () &&
        ! V.getMap ().is_null () &
        V_a_map->isSameAs (* (V.getMap ()));
      TEST_ASSERT( maps_ok );

      auto V_v_h = V.getVectorView ().getLocalViewHost(Tpetra::Access::ReadOnly);
      auto V_a_v_h = V_a.getVectorView ().getLocalViewHost(Tpetra::Access::ReadOnly);
      TEST_ASSERT( V_v_h.extent (0) == V_a_v_h.extent (0) &&
                   V_v_h.extent (1) == V_a_v_h.extent (1) &&
                   V_v_h.data () == V_a_v_h.data () );
    }
  }

  // Test BlockMultiVector::getMultiVectorView.  It must return a
  // MultiVector with view semantics, and it must view the same data
  // that the BlockMultiVector view.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockMultiVector, MVView, Scalar, LO, GO, Node )
  {
    using Teuchos::Comm;
    using Teuchos::RCP;
    typedef Tpetra::BlockMultiVector<Scalar, LO, GO, Node> BMV;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::global_size_t GST;
    typedef Teuchos::ScalarTraits<Scalar> STS;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    const size_t numLocalMeshPoints = 12;
    const LO blockSize = 4;
    const LO numVecs = 3;

    const GO indexBase = 0;
    map_type meshMap (INVALID, numLocalMeshPoints, indexBase, comm);
    //RCP<const map_type> mapPtr = Teuchos::rcpFromRef (map); // nonowning RCP

    BMV X (meshMap, blockSize, numVecs);
    //X.sync_host();

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
    TEST_EQUALITY_CONST( mv_map.getLocalNumElements (), static_cast<size_t> (48) );

    // X has three columns.  In the following tests, change only
    // blocks in column 1.  Entries in all other columns must remain
    // unmodified.
    const LO colToModify = 1;

    // Test that X_mv views the same data that the BMV views.  X_mv
    // has 48 rows on each process: 12 mesh points, and 4 degrees of
    // freedom per mesh point ("block size").  Rows 20-23 thus
    // correspond to local mesh point 5.
    auto X_5_1 = X.getLocalBlockHost (5, colToModify, Tpetra::Access::ReadWrite);

    // All entries of X_5_1 must be zero.  First make a block with all
    // zero entries, then test.  It's not worth testing the
    // corresponding entries of X_mv yet; we'll do that below.
    typename BMV::little_host_vec_type zeroLittleVector ("zero", blockSize);
    Kokkos::deep_copy(zeroLittleVector, STS::zero());

    TEST_ASSERT( equal (X_5_1, zeroLittleVector) && equal (zeroLittleVector, X_5_1) );

    // Put some data in the block.  This will help us test whether the
    // corresponding entries in the MultiVector, and only those
    // entries, got modified.
    for (LO i = 0; i < blockSize; ++i) {
      X_5_1(i) = static_cast<Scalar> (i + 1); // all are nonzero
    }
    TEST_ASSERT( ! equal (X_5_1, zeroLittleVector) && ! equal (zeroLittleVector, X_5_1) );

    // Make sure that getLocalBlockHost() returns a read-and-write view,
    // not a deep copy.  Do this by calling getLocalBlockHost(5,1) again,
    // and testing that changes to X_5_1 are reflected in the result.
    auto X_5_1_new = X.getLocalBlockHost (5, colToModify, Tpetra::Access::ReadOnly);
    TEST_ASSERT( equal (X_5_1_new, X_5_1) && equal (X_5_1, X_5_1_new) );
    TEST_ASSERT( ! equal (X_5_1_new, zeroLittleVector) &&
                 ! equal (zeroLittleVector, X_5_1_new) );

    // Make sure that all blocks other than block 5 still contain zeros.
    for (LO localMeshIndex = 0;
         localMeshIndex < static_cast<LO> (numLocalMeshPoints);
         ++localMeshIndex) {
      for (LO curCol = 0; curCol < numVecs; ++curCol) {
        auto X_cur = X.getLocalBlockHost (localMeshIndex, curCol,
                                      Tpetra::Access::ReadOnly);
        if (curCol != colToModify) {
          TEST_ASSERT( equal (X_cur, zeroLittleVector) &&
                       equal (zeroLittleVector, X_cur) );
          TEST_ASSERT( ! equal (X_cur, X_5_1) && ! equal (X_5_1, X_cur) );
        }
        if (localMeshIndex != 5) {
          TEST_ASSERT( equal (X_cur, zeroLittleVector) &&
                       equal (zeroLittleVector, X_cur) );
          TEST_ASSERT( ! equal (X_cur, X_5_1) && ! equal (X_5_1, X_cur) );
        }
      }
    }

    // Make sure that all entries of X_mv (the MultiVector) contain
    // zeros, except for rows 20-23 in column colToModify.  Those must
    // contain 1,2,3,4 (in order).
    const size_t numLocalDofs = mv_map.getLocalNumElements ();
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
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockMultiVector, Import, Scalar, LO, GO, Node )
  {
    using Teuchos::Comm;
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using Teuchos::RCP;
    typedef Tpetra::BlockMultiVector<Scalar, LO, GO, Node> BMV;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::Import<LO, GO, Node> import_type;
    typedef Tpetra::global_size_t GST;
    typedef Teuchos::ScalarTraits<Scalar> STS;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    // const int myRank = comm->getRank ();
    // const int numProcs = comm->getSize ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    const size_t numLocalMeshPoints = 2;
    const LO blockSize = 5;
    const LO numVecs = 3;
    const GO indexBase = 0;
    map_type meshMap (INVALID, numLocalMeshPoints, indexBase, comm);
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
    map_type overlappingMeshMap (INVALID, overlappingGIDs (), indexBase, comm);

    BMV X (meshMap, blockSize, numVecs);
    BMV Y (overlappingMeshMap, blockSize, numVecs);
    //X.sync_host();
    //Y.sync_host();

    //
    // Fill X with meaningful things to test Import with REPLACE combine mode.
    //
    // Fill just the "outgoing" blocks of X, in column 1 only, with
    // 1,2,3,4,5.  Due to how we constructed the overlapping mesh Map,
    // the first GID on the calling process in the mesh Map overlaps
    // with one GID on exactly one process.
    const LO colToModify = 1;
    {
      auto X_overlap =
        X.getLocalBlockHost (meshMap.getLocalElement (meshMap.getMinGlobalIndex ()), 
                         colToModify, Tpetra::Access::ReadWrite);
      TEST_ASSERT( X_overlap.data () != NULL );
      TEST_EQUALITY_CONST( static_cast<size_t> (X_overlap.extent (0)),
                           static_cast<size_t> (blockSize) );

      const int lclOk = (X_overlap.data () != NULL &&
                         static_cast<size_t> (X_overlap.extent (0)) == static_cast<size_t> (blockSize)) ? 1 : 0;
      int gblOk = 1;
      reduceAll<int, int> (*comm, REDUCE_MIN, lclOk, outArg (gblOk));
      TEUCHOS_TEST_FOR_EXCEPTION(
        gblOk == 0, std::logic_error, "Some process reported that X_overlap "
        "is an empty block.");

      for (LO i = 0; i < blockSize; ++i) {
        X_overlap(i) = static_cast<Scalar> (i+1);
      }
    }

    { // BlockMultiVector relies on the point multivector infrastructure
      const auto pointMeshMap = BMV::makePointMap(meshMap, blockSize);
      const auto pointOverlappingMeshMap = BMV::makePointMap(overlappingMeshMap, blockSize);
      import_type pointImport (rcpFromRef (pointMeshMap), rcpFromRef (pointOverlappingMeshMap));
      Y.getMultiVectorView().doImport (X.getMultiVectorView(), pointImport, Tpetra::REPLACE);
    }

    {
      auto X_overlap =
           X.getLocalBlockHost (meshMap.getLocalElement(meshMap.getMinGlobalIndex()), 
                            colToModify, Tpetra::Access::ReadOnly);

      typename BMV::little_host_vec_type zeroLittleVector ("zero", blockSize);
      Kokkos::deep_copy(zeroLittleVector, STS::zero());

      for (LO col = 0; col < numVecs; ++col) {
        for (LO localMeshRow = meshMap.getMinLocalIndex ();
             localMeshRow < meshMap.getMaxLocalIndex (); ++localMeshRow) {
          auto Y_cur = Y.getLocalBlockHost (localMeshRow, col,
                                        Tpetra::Access::ReadOnly);

          if (col != colToModify) {
            TEST_ASSERT( equal (Y_cur, zeroLittleVector) &&
                         equal (zeroLittleVector, Y_cur) );
            TEST_ASSERT( ! equal (Y_cur, X_overlap) &&
                         ! equal (X_overlap, Y_cur) );
          }
          if (localMeshRow != meshMap.getMinLocalIndex ()) {
            TEST_ASSERT( equal (Y_cur, zeroLittleVector) &&
                         equal (zeroLittleVector, Y_cur) );
            TEST_ASSERT( ! equal (Y_cur, X_overlap) &&
                         ! equal (X_overlap, Y_cur) );
          }
          if (col == colToModify && localMeshRow == meshMap.getMinLocalIndex ()) {
            TEST_ASSERT( equal (Y_cur, X_overlap) && equal (X_overlap, Y_cur) );
          }
        }
      }
    }
  }
    

  //
  // Make sure that BlockMultiVector's "offset view" constructors work.
  //
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockMultiVector, OffsetView, Scalar, LO, GO, Node )
  {
    using Teuchos::Comm;
    using Teuchos::RCP;
    typedef Tpetra::BlockMultiVector<Scalar, LO, GO, Node> BMV;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef typename map_type::device_type device_type;
    typedef Tpetra::global_size_t GST;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef Tpetra::MultiVector<Scalar, LO, GO, Node> MV;
    typedef typename MV::mag_type MT;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    const int numProcs = comm->getSize ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const size_t numLocalMeshPoints = 20;
    const LO blockSize = 4;
    const LO numVecs = 3;
    const GO indexBase = 0;
    map_type meshMap (INVALID, numLocalMeshPoints, indexBase, comm);

    // Create a BlockMultiVector to view.  Fill it with 1s.
    BMV X (meshMap, blockSize, numVecs);
    const Scalar one = STS::one ();
    X.putScalar (one);

    // View X as two row blocks [X1; X2], where X1 has 11 rows and X2
    // has 9 rows.  Note that X2's Map has 9 rows, as X2 should have;
    // we account for X2's starting position with the offset argument
    // to the offset view constructor.
    map_type X1_meshMap (INVALID, (meshMap.getLocalElementList ()) (0, 11), indexBase, comm);
    map_type X2_meshMap (INVALID, (meshMap.getLocalElementList ()) (11, 9), indexBase+11, comm);

    TEST_EQUALITY_CONST( X1_meshMap.getLocalNumElements (), static_cast<size_t> (11) );
    TEST_EQUALITY_CONST( X2_meshMap.getLocalNumElements (), static_cast<size_t> (9) );

    BMV X1 (X, X1_meshMap);
    BMV X2 (X, X2_meshMap, static_cast<size_t> (11));

    TEST_EQUALITY( X1.getMultiVectorView ().getLocalLength (),
                   static_cast<size_t> (11 * blockSize) );
    TEST_EQUALITY( X2.getMultiVectorView ().getLocalLength (),
                   static_cast<size_t> (9 * blockSize) );
    TEST_EQUALITY( X1.getMultiVectorView ().getNumVectors (),
                   static_cast<size_t> (numVecs) );
    TEST_EQUALITY( X2.getMultiVectorView ().getNumVectors (),
                   static_cast<size_t> (numVecs) );

    Kokkos::View<MT*, device_type> norms ("norms", numVecs);
    auto norms_h = Kokkos::create_mirror_view (norms);

    // Fill X1 with 2s and X2 with 3s.
    // Since X1 and X2 are both views, this should affect X.
    const Scalar two = one + one;
    X1.putScalar (two);
    const Scalar three = two + one;
    X2.putScalar (three);

    // The one-norm of X1 should be numProcs*11*2*blockSize, ...
    const MT X1_expectedOneNorm =
      static_cast<MT> (numProcs) * static_cast<MT> (11) *
      static_cast<MT> (2) * static_cast<MT> (blockSize);
    X1.getMultiVectorView ().norm1 (norms);
    Kokkos::deep_copy (norms_h, norms);
    for (LO j = 0; j < numVecs; ++j) {
      TEST_EQUALITY( norms_h(j), X1_expectedOneNorm );
    }

    // ... and the one-norm of X2 should be numProcs*9*3*blockSize.
    const MT X2_expectedOneNorm =
      static_cast<MT> (numProcs) * static_cast<MT> (9) *
      static_cast<MT> (3) * static_cast<MT> (blockSize);
    X2.getMultiVectorView ().norm1 (norms);
    Kokkos::deep_copy (norms_h, norms);
    for (LO j = 0; j < numVecs; ++j) {
      TEST_EQUALITY( norms_h(j), X2_expectedOneNorm );
    }

    // The one-norm of X should now be numProcs*(11*2 + 9*3)*blockSize.
    const MT X_expectedOneNorm = static_cast<MT> (numProcs) *
      (static_cast<MT> (11) * static_cast<MT> (2) +
       static_cast<MT> (9) * static_cast<MT> (3)) *
      static_cast<MT> (blockSize);
    X.getMultiVectorView ().norm1 (norms);
    Kokkos::deep_copy (norms_h, norms);
    for (LO j = 0; j < numVecs; ++j) {
      TEST_EQUALITY( norms_h(j), X_expectedOneNorm );
    }

    // Repeat the test, but exercise the four-argument offset view
    // constructor instead.
    X.putScalar (one);
    X1 = BMV (X, X1_meshMap, BMV::makePointMap (X1_meshMap, blockSize));
    X2 = BMV (X, X2_meshMap, BMV::makePointMap (X2_meshMap, blockSize),
              static_cast<size_t> (11));

    TEST_EQUALITY( X1.getMultiVectorView ().getLocalLength (),
                   static_cast<size_t> (11 * blockSize) );
    TEST_EQUALITY( X2.getMultiVectorView ().getLocalLength (),
                   static_cast<size_t> (9 * blockSize) );
    TEST_EQUALITY( X1.getMultiVectorView ().getNumVectors (),
                   static_cast<size_t> (numVecs) );
    TEST_EQUALITY( X2.getMultiVectorView ().getNumVectors (),
                   static_cast<size_t> (numVecs) );

    // Fill X1 with 2s and X2 with 3s.
    // Since X1 and X2 are both views, this should affect X.
    X1.putScalar (two);
    X2.putScalar (three);

    // The one-norm of X1 should be numProcs*11*2*blockSize, ...
    X1.getMultiVectorView ().norm1 (norms);
    Kokkos::deep_copy (norms_h, norms);
    for (LO j = 0; j < numVecs; ++j) {
      TEST_EQUALITY( norms_h(j), X1_expectedOneNorm );
    }

    // ... and the one-norm of X2 should be numProcs*9*3*blockSize.
    X2.getMultiVectorView ().norm1 (norms);
    Kokkos::deep_copy (norms_h, norms);
    for (LO j = 0; j < numVecs; ++j) {
      TEST_EQUALITY( norms_h(j), X2_expectedOneNorm );
    }

    // The one-norm of X should now be numProcs*(11*2 + 9*3)*blockSize.
    X.getMultiVectorView ().norm1 (norms);
    Kokkos::deep_copy (norms_h, norms);
    for (LO j = 0; j < numVecs; ++j) {
      TEST_EQUALITY( norms_h(j), X_expectedOneNorm );
    }
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockMultiVector, ctor, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockMultiVector, MVView, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockMultiVector, Import, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockMultiVector, OffsetView, SCALAR, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP )

} // namespace (anonymous)
