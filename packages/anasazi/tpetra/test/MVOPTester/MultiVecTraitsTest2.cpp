#include <Teuchos_UnitTestHarness.hpp>
#include <Tpetra_DefaultPlatform.hpp>
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

    RCP<const Comm<int> > comm =
      Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
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
      if (map->getNodeNumElements () != 0) {
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
      TEST_ASSERT( norms[j] == X_norms[j] );
    }

    out << "Create MultiVector Y" << endl;

    // Create a MultiVector Y, and make it a deep copy of X.
    MV Y (map, numCols);
    Tpetra::deep_copy (Y, X);

    // Test that the norms of the columns of Y have their expected values.
    Y.norm1 (norms);
    for (size_t j = 0; j < numCols; ++j) {
      TEST_ASSERT( norms[j] == X_norms[j] );
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
      TEST_ASSERT( norms[j] == X_norms[j] );
    }

    out << "Test the norms of the columns of Y" << endl;

    // Test that the norm of the first column of Y has not changed.
    Y.norm1 (norms);
    TEST_ASSERT( norms[0] == X_norms[0] );

    // Test that the norms of the remaining columns of Y have their
    // expected values.
    TEST_ASSERT( norms[1] == X_norms[0] );
    TEST_ASSERT( norms[2] == X_norms[1] );

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
    TEST_ASSERT( norms[0] == X_norms[0] );
    TEST_ASSERT( norms[1] == X_norms[0] );
    TEST_ASSERT( norms[2] == X_norms[1] );
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

    RCP<const Comm<int> > comm =
      Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
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
      if (map->getNodeNumElements () != 0) {
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
      TEST_ASSERT( norms[j] == X_norms[j] );
    }

    out << "Create MultiVector Y" << endl;

    // Create a MultiVector Y, and make it a deep copy of X.
    MV Y (map, numCols);
    Tpetra::deep_copy (Y, X);

    // Test that the norms of the columns of Y have their expected values.
    Y.norm1 (norms);
    for (size_t j = 0; j < numCols; ++j) {
      TEST_ASSERT( norms[j] == X_norms[j] );
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
      TEST_ASSERT( norms[j] == X_norms[j] );
    }

    out << "Test the norms of the columns of Y" << endl;

    // Test that the norm of the first column of Y has not changed.
    Y.norm1 (norms);
    TEST_ASSERT( norms[0] == X_norms[0] );

    // Test that the norms of the remaining columns of Y have their
    // expected values.
    TEST_ASSERT( norms[1] == X_norms[0] );
    TEST_ASSERT( norms[2] == X_norms[1] );

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
    TEST_ASSERT( norms[0] == X_norms[0] );
    TEST_ASSERT( norms[1] == X_norms[0] );
    TEST_ASSERT( norms[2] == X_norms[1] );
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

    RCP<const Comm<int> > comm =
      Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
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
      if (map->getNodeNumElements () != 0) {
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
      TEST_ASSERT( norms[j] == X_norms[j] );
    }

    // Create a MultiVector Y, and make it a deep copy of X.
    MV Y (map, numCols);
    Tpetra::deep_copy (Y, X);

    // Test that the norms of the columns of Y have their expected values.
    Y.norm1 (norms);
    for (size_t j = 0; j < numCols; ++j) {
      TEST_ASSERT( norms[j] == X_norms[j] );
    }

    // Y(:, [2,1]) := X(:, [0,1])
    std::vector<int> index (2);
    index[0] = 2;
    index[1] = 1;
    MVT::SetBlock (X, index, Y);

    // Test that the norms of the columns of X have not changed.
    X.norm1 (norms);
    for (size_t j = 0; j < numCols; ++j) {
      TEST_ASSERT( norms[j] == X_norms[j] );
    }

    // Test that the norm of the first column of Y has not changed.
    Y.norm1 (norms);
    TEST_ASSERT( norms[0] == X_norms[0] );

    // Test that the norms of the remaining columns of Y have their
    // expected values.
    TEST_ASSERT( norms[1] == X_norms[1] );
    TEST_ASSERT( norms[2] == X_norms[0] );
  }

} // namespace (anonymous)
