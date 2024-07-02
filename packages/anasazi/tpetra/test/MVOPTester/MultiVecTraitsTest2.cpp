// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include <Teuchos_UnitTestHarness.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>

#include <AnasaziConfigDefs.hpp>
#include <AnasaziTpetraAdapter.hpp>

//
// Test the Tpetra specialization of MultiVecTraits::SetBlock.
//

namespace {
  using Tpetra::Map;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::Comm;
  using Teuchos::REDUCE_SUM;
  using Teuchos::reduceAll;
  using Teuchos::RCP;
  using std::endl;
  typedef Tpetra::MultiVector<> MV;
  typedef Tpetra::Vector<> V;
  typedef Tpetra::global_size_t GST;

  // We don't need "typename" here because MV is not templated.  All
  // of its template parameters take their default values, so it's a
  // concrete type.
  typedef MV::scalar_type scalar_type;
  typedef MV::local_ordinal_type LO;
  typedef MV::global_ordinal_type GO;
  typedef Anasazi::MultiVecTraits<scalar_type, MV> MVT;
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef STS::magnitudeType norm_type;
  typedef Teuchos::ScalarTraits<norm_type> STN;

  //
  // For 3-column MultiVectors X and Y, use SetBlock with the
  // std::vector<int> index vector argument to do
  //
  // Y(:, [1,2]) := X(:, [0,1]).
  //
  TEUCHOS_UNIT_TEST( MultiVecTraits, TpetraSetBlock1 )
  {
    const GST INV = Teuchos::OrdinalTraits<GST>::invalid ();
    const size_t numLclRows = 5;
    const size_t numCols = 3;
    const GO indexBase = 1; // just for a change

    out << "Test SetBlock with index vector: Y(:, [1,2]) := X(:, [0,1])"
        << endl;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    RCP<const Map<> > map (new Map<> (INV, numLclRows, indexBase, comm));
    const GST numGblRows = map->getGlobalNumElements ();

    out << "Creating and filling MultiVector X" << endl;

    // Create a MultiVector X, and fill it such that each entry has a
    // globally unique value.  This will help us test SetBlock both
    // interactively and automatically.
    MV X (map, numCols);
    // Keep track of what the one-norm of each column of X should be.
    Array<norm_type> X_norms (numCols);
    for (size_t j = 0; j < numCols; ++j) {
      ArrayRCP<scalar_type> X_j = X.getDataNonConst (j);

      X_norms[j] = STN::zero ();
      if (map->getLocalNumElements () != 0) {
        const LO myLclRowMin = map->getMinLocalIndex ();
        const LO myLclRowMax = map->getMaxLocalIndex ();

        for (LO i_lcl = myLclRowMin; i_lcl <= myLclRowMax; ++i_lcl) {
          const GO i_gbl = map->getGlobalElement (i_lcl);
          // Some value which is unique for each MultiVector entry.
          const scalar_type X_ij = static_cast<scalar_type> (i_gbl) +
            static_cast<scalar_type> (numGblRows * numCols);
          X_j[i_lcl] = X_ij;
          X_norms[j] += STS::magnitude (X_ij);
        }
      }
    }

    // Compute what the one-norm of each column of X should be.
    {
      Array<norm_type> X_norms_in (X_norms); // deep copy
      reduceAll<int, norm_type> (*comm, REDUCE_SUM, static_cast<int> (numCols),
                                 X_norms_in.getRawPtr (), X_norms.getRawPtr ());
    }

    // Test that the norms of the columns of X have their expected values.
    Array<norm_type> norms (numCols);
    X.norm1 (norms);
    for (size_t j = 0; j < numCols; ++j) {
      out << "Expected 1-norm of X_" << j << ": " << X_norms[j] << endl;
      out << "Actual 1-norm of X_" << j << ": " << norms[j] << endl;
      TEST_ASSERT( X_norms[j] != STN::zero () );
      TEST_FLOATING_EQUALITY( norms[j], X_norms[j], STS::eps() );
    }

    out << "Create MultiVector Y" << endl;

    // Create a MultiVector Y, and make it a deep copy of X.
    MV Y (map, numCols);
    Tpetra::deep_copy (Y, X);

    // Test that the norms of the columns of Y have their expected values.
    Y.norm1 (norms);
    for (size_t j = 0; j < numCols; ++j) {
      TEST_FLOATING_EQUALITY( norms[j], X_norms[j], STS::eps() );
    }

    out << "Call SetBlock(X, [1, 2], Y)" << endl;

    // Y(:, [1,2]) := X(:, [0,1])
    std::vector<int> index (2);
    index[0] = 1;
    index[1] = 2;
    MVT::SetBlock (X, index, Y);

    out << "Test that the norms of the columns of X have not changed" << endl;

    // Test that the norms of the columns of X have not changed.
    X.norm1 (norms);
    for (size_t j = 0; j < numCols; ++j) {
      TEST_FLOATING_EQUALITY( norms[j], X_norms[j], STS::eps() );
    }

    out << "Test the norms of the columns of Y" << endl;

    // Test that the norm of the first column of Y has not changed.
    Y.norm1 (norms);
    TEST_FLOATING_EQUALITY( norms[0], X_norms[0], STS::eps() );

    // Test that the norms of the remaining columns of Y have their
    // expected values.
    TEST_FLOATING_EQUALITY( norms[1], X_norms[0], STS::eps() );
    TEST_FLOATING_EQUALITY( norms[2], X_norms[1], STS::eps() );

    out << "Test the values in the columns of Y" << endl;

    // Test that after calling SetBlock, Y(:,1:2) == X(:,0:1) (where
    // 0:1 means [0, 1] (inclusive range) and 1:2 means [1, 2]).
    MV Z (map, static_cast<size_t> (2));
    RCP<const MV> X_view = X.subView (Teuchos::Range1D (0, 1));
    RCP<const MV> Y_view = Y.subView (Teuchos::Range1D (1, 2));
    // Z := X(:, 0:1) - Y(: 1:2)
    Z.update (STS::one (), *X_view, -STS::one (), *Y_view, STS::zero ());
    Z.norm1 (norms (0, 2));
    TEST_ASSERT( norms[0] == STN::zero () );
    TEST_ASSERT( norms[1] == STN::zero () );

    // Test that overwriting X doesn't affect the new values of Y.
    X.putScalar (STS::zero ());
    Y.norm1 (norms);
    TEST_FLOATING_EQUALITY( norms[0], X_norms[0], STS::eps() );
    TEST_FLOATING_EQUALITY( norms[1], X_norms[0], STS::eps() );
    TEST_FLOATING_EQUALITY( norms[2], X_norms[1], STS::eps() );
  }

  //
  // Repeat TpetraSetBlock1, but use Teuchos::Range1D instead of an
  // index vector: Y(:, [1,2]) := X(:, [0,1]).
  //
  TEUCHOS_UNIT_TEST( MultiVecTraits, TpetraSetBlock2 )
  {
    const GST INV = Teuchos::OrdinalTraits<GST>::invalid ();
    const size_t numLclRows = 5;
    const size_t numCols = 3;
    const GO indexBase = 1; // just for a change

    out << "Test SetBlock with Range1D: Y(:, [1,2]) := X(:, [0,1])"
        << endl;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    RCP<const Map<> > map (new Map<> (INV, numLclRows, indexBase, comm));
    const GST numGblRows = map->getGlobalNumElements ();

    out << "Creating and filling MultiVector X" << endl;

    // Create a MultiVector X, and fill it such that each entry has a
    // globally unique value.  This will help us test SetBlock both
    // interactively and automatically.
    MV X (map, numCols);
    // Keep track of what the one-norm of each column of X should be.
    Array<norm_type> X_norms (numCols);
    for (size_t j = 0; j < numCols; ++j) {
      ArrayRCP<scalar_type> X_j = X.getDataNonConst (j);

      X_norms[j] = STN::zero ();
      if (map->getLocalNumElements () != 0) {
        const LO myLclRowMin = map->getMinLocalIndex ();
        const LO myLclRowMax = map->getMaxLocalIndex ();

        for (LO i_lcl = myLclRowMin; i_lcl <= myLclRowMax; ++i_lcl) {
          const GO i_gbl = map->getGlobalElement (i_lcl);
          // Some value which is unique for each MultiVector entry.
          const scalar_type X_ij = static_cast<scalar_type> (i_gbl) +
            static_cast<scalar_type> (numGblRows * numCols);
          X_j[i_lcl] = X_ij;
          X_norms[j] += STS::magnitude (X_ij);
        }
      }
    }

    // Compute what the one-norm of each column of X should be.
    {
      Array<norm_type> X_norms_in (X_norms); // deep copy
      reduceAll<int, norm_type> (*comm, REDUCE_SUM, static_cast<int> (numCols),
                                 X_norms_in.getRawPtr (), X_norms.getRawPtr ());
    }

    // Test that the norms of the columns of X have their expected values.
    Array<norm_type> norms (numCols);
    X.norm1 (norms);
    for (size_t j = 0; j < numCols; ++j) {
      out << "Expected 1-norm of X_" << j << ": " << X_norms[j] << endl;
      out << "Actual 1-norm of X_" << j << ": " << norms[j] << endl;
      TEST_ASSERT( X_norms[j] != STN::zero () );
      TEST_FLOATING_EQUALITY( norms[j], X_norms[j], STS::eps() );
    }

    out << "Create MultiVector Y" << endl;

    // Create a MultiVector Y, and make it a deep copy of X.
    MV Y (map, numCols);
    Tpetra::deep_copy (Y, X);

    // Test that the norms of the columns of Y have their expected values.
    Y.norm1 (norms);
    for (size_t j = 0; j < numCols; ++j) {
      TEST_FLOATING_EQUALITY( norms[j], X_norms[j], STS::eps() );
    }

    out << "Call SetBlock(X, Range1D(1, 2), Y)" << endl;

    // Repeat TpetraSetBlock1, but use Teuchos::Range1D instead of an
    // index vector: Y(:, [1,2]) := X(:, [0,1])
    Teuchos::Range1D colRng (1, 2); // inclusive range
    MVT::SetBlock (X, colRng, Y);

    out << "Test that the norms of the columns of X have not changed" << endl;

    // Test that the norms of the columns of X have not changed.
    X.norm1 (norms);
    for (size_t j = 0; j < numCols; ++j) {
      TEST_FLOATING_EQUALITY( norms[j], X_norms[j], STS::eps() );
    }

    out << "Test the norms of the columns of Y" << endl;

    // Test that the norm of the first column of Y has not changed.
    Y.norm1 (norms);
    TEST_FLOATING_EQUALITY( norms[0], X_norms[0], STS::eps() );

    // Test that the norms of the remaining columns of Y have their
    // expected values.
    TEST_FLOATING_EQUALITY( norms[1], X_norms[0], STS::eps() );
    TEST_FLOATING_EQUALITY( norms[2], X_norms[1], STS::eps() );

    // Test that after calling SetBlock, Y(:,1:2) == X(:,0:1) (where
    // 0:1 means [0, 1] (inclusive range) and 1:2 means [1, 2]).
    MV Z (map, static_cast<size_t> (2));
    RCP<const MV> X_view, Y_view;
    try {
      X_view = X.subView (Teuchos::Range1D (0, 1));
    } catch (std::exception& e) {
      out << "*** Yikes!  X.subView(Range1D(0,1)) raised an exception!  "
          << e.what ();
    }
    try {
      Y_view = Y.subView (Teuchos::Range1D (1, 2));
    } catch (std::exception& e) {
      out << "*** Yikes!  Y.subView(Range1D(1,2)) raised an exception!  "
          << e.what ();
    }
    try {
      // Z := X(:, 0:1) - Y(: 1:2)
      Z.update (STS::one (), *X_view, -STS::one (), *Y_view, STS::zero ());
    } catch (std::exception& e) {
      out << "*** Yikes!  Z.update raised an exception!  " << e.what ();
    }
    try {
      Z.norm1 (norms (0, 2));
    } catch (std::exception& e) {
      out << "*** Yikes!  Z.norm1 raised an exception!  " << e.what ();
    }
    TEST_ASSERT( norms[0] == STN::zero () );
    TEST_ASSERT( norms[1] == STN::zero () );

    out << "Test that overwriting X doesn't affect the new values of Y" << endl;

    // Test that overwriting X doesn't affect the new values of Y.
    X.putScalar (STS::zero ());
    Y.norm1 (norms);
    TEST_FLOATING_EQUALITY( norms[0], X_norms[0], STS::eps() );
    TEST_FLOATING_EQUALITY( norms[1], X_norms[0], STS::eps() );
    TEST_FLOATING_EQUALITY( norms[2], X_norms[1], STS::eps() );
  }

  //
  // Repeat TpetraSetBlock1, but put the column indices in a different order:
  //
  // Y(:, [2,1]) := X(:, [0,1])
  //
  TEUCHOS_UNIT_TEST( MultiVecTraits, TpetraSetBlock3 )
  {
    const GST INV = Teuchos::OrdinalTraits<GST>::invalid ();
    const size_t numLclRows = 5;
    const size_t numCols = 3;
    const GO indexBase = 1; // just for a change

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    RCP<const Map<> > map (new Map<> (INV, numLclRows, indexBase, comm));
    const GST numGblRows = map->getGlobalNumElements ();

    // Create a MultiVector X, and fill it such that each entry has a
    // globally unique value.  This will help us test SetBlock both
    // interactively and automatically.
    MV X (map, numCols);
    // Keep track of what the one-norm of each column of X should be.
    Array<norm_type> X_norms (numCols);
    for (size_t j = 0; j < numCols; ++j) {
      ArrayRCP<scalar_type> X_j = X.getDataNonConst (j);

      X_norms[j] = STN::zero ();
      if (map->getLocalNumElements () != 0) {
        const LO myLclRowMin = map->getMinLocalIndex ();
        const LO myLclRowMax = map->getMaxLocalIndex ();

        for (LO i_lcl = myLclRowMin; i_lcl <= myLclRowMax; ++i_lcl) {
          const GO i_gbl = map->getGlobalElement (i_lcl);
          // Some value which is unique for each MultiVector entry.
          const scalar_type X_ij = static_cast<scalar_type> (i_gbl) +
            static_cast<scalar_type> (numGblRows * numCols);
          X_j[i_lcl] = X_ij;
          X_norms[j] += STS::magnitude (X_ij);
        }
      }
    }

    // Compute what the one-norm of each column of X should be.
    {
      Array<norm_type> X_norms_in (X_norms); // deep copy
      reduceAll<int, norm_type> (*comm, REDUCE_SUM, static_cast<int> (numCols),
                                 X_norms_in.getRawPtr (), X_norms.getRawPtr ());
    }

    // Test that the norms of the columns of X have their expected values.
    Array<norm_type> norms (numCols);
    X.norm1 (norms);
    for (size_t j = 0; j < numCols; ++j) {
      out << "Expected 1-norm of X_" << j << ": " << X_norms[j] << endl;
      out << "Actual 1-norm of X_" << j << ": " << norms[j] << endl;
      TEST_ASSERT( X_norms[j] != STN::zero () );
      TEST_FLOATING_EQUALITY( norms[j], X_norms[j], STS::eps() );
    }

    // Create a MultiVector Y, and make it a deep copy of X.
    MV Y (map, numCols);
    Tpetra::deep_copy (Y, X);

    // Test that the norms of the columns of Y have their expected values.
    Y.norm1 (norms);
    for (size_t j = 0; j < numCols; ++j) {
      TEST_FLOATING_EQUALITY( norms[j], X_norms[j], STS::eps() );
    }

    // Y(:, [2,1]) := X(:, [0,1])
    std::vector<int> index (2);
    index[0] = 2;
    index[1] = 1;
    MVT::SetBlock (X, index, Y);

    // Test that the norms of the columns of X have not changed.
    X.norm1 (norms);
    for (size_t j = 0; j < numCols; ++j) {
      TEST_FLOATING_EQUALITY( norms[j], X_norms[j], STS::eps() );
    }

    // Test that the norm of the first column of Y has not changed.
    Y.norm1 (norms);
    TEST_FLOATING_EQUALITY( norms[0], X_norms[0], STS::eps() );

    // Test that the norms of the remaining columns of Y have their
    // expected values.
    TEST_FLOATING_EQUALITY( norms[1], X_norms[1], STS::eps() );
    TEST_FLOATING_EQUALITY( norms[2], X_norms[0], STS::eps() );
  }

  //
  // Replicate what AnasaziMVOPTester.hpp does to test SetBlock.
  //
  TEUCHOS_UNIT_TEST( MultiVecTraits, TpetraSetBlock4 )
  {
    const GST INV = Teuchos::OrdinalTraits<GST>::invalid ();
    const size_t numLclRows = 5;
    const size_t numCols = 3;
    const GO indexBase = 1; // just for a change

    out << "Test SetBlock by imitating AnasaziMVOPTester.hpp" << endl;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    RCP<const Map<> > map (new Map<> (INV, numLclRows, indexBase, comm));
    const GST numGblRows = map->getGlobalNumElements ();

    out << "Creating a \"prototype\" MultiVector X" << endl;

    // Create a MultiVector X, and fill it such that each entry has a
    // globally unique value.  This will help us test SetBlock both
    // interactively and automatically.
    MV X (map, numCols);
    // Keep track of what the one-norm of each column of X should be.
    Array<norm_type> X_norms (numCols);
    for (size_t j = 0; j < numCols; ++j) {
      ArrayRCP<scalar_type> X_j = X.getDataNonConst (j);

      X_norms[j] = STN::zero ();
      if (map->getLocalNumElements () != 0) {
        const LO myLclRowMin = map->getMinLocalIndex ();
        const LO myLclRowMax = map->getMaxLocalIndex ();

        for (LO i_lcl = myLclRowMin; i_lcl <= myLclRowMax; ++i_lcl) {
          const GO i_gbl = map->getGlobalElement (i_lcl);
          // Some value which is unique for each MultiVector entry.
          const scalar_type X_ij = static_cast<scalar_type> (i_gbl) +
            static_cast<scalar_type> (numGblRows * numCols);
          X_j[i_lcl] = X_ij;
          X_norms[j] += STS::magnitude (X_ij);
        }
      }
    }

    // Compute what the one-norm of each column of X should be.
    {
      Array<norm_type> X_norms_in (X_norms); // deep copy
      reduceAll<int, norm_type> (*comm, REDUCE_SUM, static_cast<int> (numCols),
                                 X_norms_in.getRawPtr (), X_norms.getRawPtr ());
    }

    // Test that the norms of the columns of X have their expected values.
    Array<norm_type> norms (numCols);
    X.norm1 (norms);
    for (size_t j = 0; j < numCols; ++j) {
      out << "Expected 1-norm of X_" << j << ": " << X_norms[j] << endl;
      out << "Actual 1-norm of X_" << j << ": " << norms[j] << endl;
      TEST_ASSERT( X_norms[j] != STN::zero () );
      TEST_FLOATING_EQUALITY( norms[j], X_norms[j], STS::eps() );
    }

    out << "Create (non-copy) clones of X: B and C" << endl;

    const size_t numVecsB = 10;
    const size_t numVecsC = 5;

    RCP<MV> B = MVT::Clone (X, numVecsB);
    TEST_ASSERT( static_cast<size_t> (B->getNumVectors ()) == numVecsB );
    RCP<MV> C = MVT::Clone (X, numVecsC);
    TEST_ASSERT( static_cast<size_t> (C->getNumVectors ()) == numVecsC );

    out << "B has " << B->getNumVectors () << " vectors; "
      "C has " << C->getNumVectors () << " vectors" << endl;

    out << "Fill B and C with random numbers" << endl;
    MVT::MvRandom (*B);
    MVT::MvRandom (*C);

    std::vector<norm_type> normsB1 (numVecsB);
    std::vector<norm_type> normsB2 (numVecsB);
    std::vector<norm_type> normsC1 (numVecsC);
    std::vector<norm_type> normsC2 (numVecsC);

    out << "Compute norms of B and C before SetBlock:" << endl;
    {
      Teuchos::OSTab tab1 (out);
      MVT::MvNorm (*B, normsB1);
      MVT::MvNorm (*C, normsC1);
      out << "Norms of B: [";
      for (size_t j = 0; j < numVecsB; ++j) {
        out << normsB1[j];
        if (j + 1 != numVecsB) {
          out << ", ";
        }
      }
      out << "]" << endl << "Norms of C: [";
      for (size_t j = 0; j < numVecsC; ++j) {
        out << normsC1[j];
        if (j + 1 != numVecsC) {
          out << ", ";
        }
      }
      out << "]" << endl;
    }

    std::vector<int> ind (numVecsC);
    for (size_t j = 0; j < numVecsC; ++j) {
      ind[j] = static_cast<int> (2 * j);
    }

    // mfh 11 Sep 2014: SetBlock works by calling CloneViewNonConst on
    // the target MultiVector, and doing a deep copy from the source
    // MultiVector to the resulting view of the target MultiVector.
    // If there are problems with SetBlock, then CloneViewNonConst
    // might be a cause.
    {
      out << "Check B_view = CloneViewNonConst(B, ind):" << endl;
      Teuchos::OSTab tab1 (out);

      out << "ind: [";
      for (size_t j = 0; j < static_cast<size_t> (ind.size ()); ++j) {
        out << ind[j];
        if (j + 1 != static_cast<size_t> (ind.size ())) {
          out << ", ";
        }
      }
      out << "]" << endl;
      RCP<MV> B_view = MVT::CloneViewNonConst (*B, ind);
      TEST_EQUALITY( static_cast<size_t> (B_view->getNumVectors ()),
                     static_cast<size_t> (ind.size ()) );

      std::vector<norm_type> B_view_norms (ind.size ());
      MVT::MvNorm (*B_view, B_view_norms);
      out << "norms of CloneViewNonConst(B, ind): [";
      for (size_t j = 0; j < static_cast<size_t> (ind.size ()); ++j) {
        out << B_view_norms[j];
        if (j + 1 != static_cast<size_t> (ind.size ())) {
          out << ", ";
        }
      }
      out << "]" << endl;

      for (size_t j = 0; j < static_cast<size_t> (ind.size ()); ++j) {
        TEST_FLOATING_EQUALITY( B_view_norms[j], normsB1[ind[j]], STS::eps() );
      }

      out << "Check that modifying B_view modifies the corresponding columns "
        "of B (we'll put B back)" << endl;
      Teuchos::OSTab tab2 (out);

      // Make a temporary copy of B.  Use native Tpetra calls in case
      // MultiVecTraits is broken.
      out << "Make a temporary copy of B" << endl;
      MV B_copy (*B, Teuchos::Copy);
      B_copy.setCopyOrView (Teuchos::View);

      out << "Test that a view of column 1 of B_view views column 2 of B" << endl;
      RCP<V> B_view_1 = B_view->getVectorNonConst (1);
      RCP<const V> B_2 = B->getVector (2);
      const norm_type B_view_1_norm2_old = B_view_1->norm2 ();
      const norm_type B_2_norm2_old = B_2->norm2 ();
      TEST_EQUALITY( B_view_1_norm2_old, B_2_norm2_old );

      out << "Set all the entries of B_view_1 to one" << endl;
      B_view_1->putScalar (STS::one ());
      const norm_type B_view_1_norm2 = B_view_1->norm2 ();
      // Make sure that setting the entries worked.
      TEST_INEQUALITY( B_view_1_norm2, STN::zero () );
      // In theory, equality is possible, since B was filled with
      // random numbers.  However, it should be unlikely.
      TEST_INEQUALITY( B_view_1_norm2_old, B_view_1_norm2 );

      out << "The norm of B(:,2) should equal the norm of B_view(:,1)" << endl;
      const norm_type B_2_norm2 = B_2->norm2 ();
      TEST_EQUALITY( B_view_1_norm2, B_2_norm2 );
      // Make sure that no fishy stuff (e.g., copy instead of view) is
      // going on with getVector().
      B_2 = B->getVector (2);
      TEST_EQUALITY( B_2->norm2 (), B_2_norm2 );

      out << "Restore B from its copy" << endl;
      Tpetra::deep_copy (*B, B_copy);
      // Check that B was actually restored.
      Array<norm_type> B_norms (numVecsB);
      B->norm2 (B_norms);
      for (size_t j = 0; j < numVecsB; ++j) {
        TEST_FLOATING_EQUALITY( B_norms[j], normsB1[j], STS::eps() );
      }

      out << "Make sure that \\|B_view_1\\| still equals \\|B(:,2)\\|" << endl;
      const norm_type B_view_1_norm2_new = B_view_1->norm2 ();
      TEST_FLOATING_EQUALITY( normsB1[2], B_view_1_norm2_new, STS::eps() );
      const norm_type B_2_norm2_new = B_2->norm2 ();
      TEST_FLOATING_EQUALITY( normsB1[2], B_2_norm2_new, STS::eps() );
      // In theory, it's possible for these two numbers to be equal,
      // since B_2_norm2 comes from a random assignment of B.
      // However, it should be unlikely.
      TEST_INEQUALITY( B_2_norm2, B_2_norm2_new );
      // Make sure that no fishy stuff (e.g., copy instead of view) is
      // going on with getVector().
      B_2 = B->getVector (2);
      TEST_FLOATING_EQUALITY( normsB1[2], B_2->norm2 (), STS::eps() );

      // TODO: Make sure that changing B_copy doesn't change B.
    }

    //
    // Moment of truth: Call SetBlock(*C, ind, *B).
    //
    out << "Call SetBlock(C, ind, B) with ind = [";
    for (size_t j = 0; j < static_cast<size_t> (ind.size ()); ++j) {
      out << ind[j];
      if (j + 1 != static_cast<size_t> (ind.size ())) {
        out << ", ";
      }
    }
    out << "]" << endl;
    MVT::SetBlock (*C, ind, *B);

    out << "Compute norms of B and C after SetBlock:" << endl;
    MVT::MvNorm (*B, normsB2);
    MVT::MvNorm (*C, normsC2);
    out << "  Norms of B: [";
    for (size_t j = 0; j < numVecsB; ++j) {
      out << normsB2[j];
      if (j + 1 != numVecsB) {
        out << ", ";
      }
    }
    out << "]" << endl << "  Norms of C: [";
    for (size_t j = 0; j < numVecsC; ++j) {
      out << normsC2[j];
      if (j + 1 != numVecsC) {
        out << ", ";
      }
    }
    out << "]" << endl;

    out << "Check that C was not changed by SetBlock" << endl;
    for (size_t j = 0; j < numVecsC; ++j) {
      TEST_FLOATING_EQUALITY( normsC1[j], normsC2[j], STS::eps() );
    }
    out << "Check that only the vectors of B that _should_ have changed did"
        << endl;
    for (size_t j = 0; j < numVecsB; ++j) {
      if (j % static_cast<size_t> (2) == 0) { // should be a vector from C
        TEST_FLOATING_EQUALITY( normsB2[j], normsC1[j/2], STS::eps() );
      }
      else { // should be an original vector
        TEST_FLOATING_EQUALITY( normsB1[j], normsB2[j], STS::eps() );
      }
    }

    out << "Verify that we copied and didn't reference" << endl;
    MVT::MvInit (*C, STS::zero ());
    MVT::MvNorm (*B, normsB1);
    for (size_t j = 0; j < numVecsB; ++j) {
      TEST_FLOATING_EQUALITY( normsB1[j], normsB2[j], STS::eps() );
    }
  }

} // namespace (anonymous)
