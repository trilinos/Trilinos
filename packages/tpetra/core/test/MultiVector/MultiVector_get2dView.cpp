// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include <type_traits>

namespace { // (anonymous)

  using Tpetra::TestingUtilities::getDefaultComm;
  using Teuchos::outArg;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::endl;

  //
  // UNIT TESTS
  //

  // Test MultiVector's get2dView and get2dViewNonConst methods with a
  // noncontiguous MultiVector input.  Commit
  // 0552ebd216150e0d1d0821517e4db94b8fbdc7af on 17 May 2016 shows
  // that these methods lacked unit tests exercising that case.  This
  // relates to Github Issue #358.  Thus, this test is particularly
  // important for CUDA builds.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, get2dView_noncontig, Scalar, LocalOrdinal, GlobalOrdinal, Node )
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::MultiVector<Scalar, LO, GO, Node> MV;
    typedef Tpetra::Vector<Scalar, LO, GO, Node> V;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename MV::mag_type mag_type;
    typedef Teuchos::ScalarTraits<mag_type> STM;

    out << "Test Tpetra::MultiVector::get2dView and get2dViewNonConst" << endl;
    Teuchos::OSTab tab1 (out);

    Scalar curVal = STS::zero ();

    // Create a Map with enough rows that we can view a noncontiguous
    // subset of columns without the ratio of rows to columns being
    // too small for a good test.
    out << "Create Map" << endl;
    auto comm = getDefaultComm ();
    const LO lclNumRows = 31;
    const GO gblNumRows =
      static_cast<GO> (comm->getSize ()) * static_cast<GO> (lclNumRows);
    const size_t numCols = 11;
    const GO indexBase = 0;
    RCP<const map_type> map =
      rcp (new map_type (gblNumRows, lclNumRows, indexBase, comm));

    // Create the MultiVector to view.
    out << "Create \"original\" MultiVector X" << endl;
    MV X (map, numCols);
    {
      Teuchos::OSTab tab2 (out);
      TEST_EQUALITY( static_cast<size_t> (X.getNumVectors ()), numCols );
      TEST_EQUALITY( static_cast<LO> (X.getLocalLength ()), lclNumRows );
      TEST_EQUALITY( static_cast<GO> (X.getGlobalLength ()), gblNumRows );
    }

    // Fill each column with a different number.  Start with 1 instead
    // of 0, because 0 is the default fill value.
    out << "Fill X" << endl;
    curVal = STS::one ();
    for (size_t j = 0; j < numCols; ++j) {
      X.getVectorNonConst (j)->putScalar (curVal);
      curVal += STS::one ();
    }

    // Make sure that fill succeeded, by computing differences between
    // expected values.
    out << "Test fill of X" << endl;
    V diff (map);
    {
      Teuchos::OSTab tab2 (out);

      curVal = STS::one ();
      for (size_t j = 0; j < numCols; ++j) {
        auto X_j = X.getVector (j);
        diff.putScalar (curVal);
        diff.update (STS::one (), *X_j, -STS::one ());
        const auto diffNorm = diff.normInf ();
        TEST_EQUALITY( diffNorm, STM::zero () );
        curVal += STS::one ();
      }
    }

    // Make a MultiVector, X_noncontig, which views a subset of columns of X.
    out << "Create X_noncontig, which views a noncontigous subset "
      "of the columns of X" << endl;
    const size_t numColsSubset = 4;
    Teuchos::Array<size_t> colsSubset (4);
    // Some columns are adjacent; some aren't.
    colsSubset[0] = 2;
    colsSubset[1] = 5;
    colsSubset[2] = 6;
    colsSubset[3] = 9;
    RCP<MV> X_noncontig = X.subViewNonConst (colsSubset ());
    {
      Teuchos::OSTab tab2 (out);
      TEST_EQUALITY( static_cast<size_t> (X_noncontig->getNumVectors ()), numColsSubset );
      TEST_EQUALITY( static_cast<LO> (X_noncontig->getLocalLength ()), lclNumRows );
      TEST_EQUALITY( static_cast<GO> (X_noncontig->getGlobalLength ()), gblNumRows );
    }

    // Make sure that the columns viewed are correct.
    out << "Test whether X_noncontig views the correct columns of X" << endl;
    {
      Teuchos::OSTab tab2 (out);

      for (size_t j = 0; j < numColsSubset; ++j) {
        const size_t col = colsSubset[j];
        out << "j = " << j << ", colsSubset[j] = " << col << endl;
        Teuchos::OSTab tab3 (out);

        auto X_col_orig = X.getVector (col);
        auto X_col = X_noncontig->getVector (j);
        Tpetra::deep_copy (diff, *X_col);

        diff.update (STS::one (), *X_col_orig, -STS::one ());
        auto diffNorm = diff.normInf ();
        TEST_EQUALITY( diffNorm, STM::zero () );

        // If we only ever test differences without actually expecting
        // a particular value, then the test will pass even if all the
        // entries of all vectors are zero.  That's why we do the
        // tests below.

        // Does the max entry have the correct value?
        const auto X_col_norm = X_col->normInf ();
        TEST_EQUALITY( X_col_norm, static_cast<mag_type> (col+1) );

        // If we fill 'diff' with the correct value and take the
        // difference, is the max of the difference zero?
        //
        // See note above for why we can't always cast directly from
        // an integer type to Scalar.
        diff.putScalar (static_cast<Scalar> (static_cast<mag_type> (col+1)));
        diff.update (STS::one (), *X_col, -STS::one ());
        diffNorm = diff.normInf ();
        TEST_EQUALITY( diffNorm, STM::zero () );
      }
    }

    // Pick a Scalar value bigger than any value currently in X, to
    // use as a "flag" value to indicate that the column was modified.
    // Casts from integer types to Scalar don't always work (for
    // example, with Scalar = std::complex<double>), so we cast first
    // to mag_type, then to Scalar.
    Scalar flagVal = static_cast<Scalar> (static_cast<mag_type> (33));

    // Make sure that X_noncontig is actually a view, not a (deep) copy.
    out << "Make sure X_noncontig is a view of X, not a copy" << endl;
    {
      Teuchos::OSTab tab2 (out);

      MV X_copy (X, Teuchos::Copy); // make a "backup" of X
      curVal = STS::one ();
      for (size_t j = 0; j < numColsSubset; ++j) {
        const size_t col = colsSubset[j];
        auto X_col_orig = X.getVectorNonConst (col);
        auto X_col = X_noncontig->getVector (j);

        X_col_orig->putScalar (flagVal);
        diff.putScalar (flagVal);
        diff.update (STS::one (), *X_col, -STS::one ());
        const auto diffNorm = diff.normInf ();
        TEST_EQUALITY( diffNorm, STM::zero () );
      }

      // Restore the original X (and thus, X_noncontig too).
      Tpetra::deep_copy (X, X_copy);
    }

    // If we didn't have success this far, we might as well quit now.
    // We haven't even tested get2dView or get2dViewNonConst yet!
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      out << "Test FAILED on some process in the communicator!" << endl;
      success = false;
      return;
    }

    // At this point, X and X_noncontig are active on device, and have
    // been most recently modified on device.  get2dView and
    // get2dViewNonConst both must sync to host.

    // The sync state of a MultiVector is only updated if the memory spaces 
    // are not the same between host and device

    out << "Make sure X_noncontig->get2dView() syncs to host" << endl;
    {
      Teuchos::OSTab tab2 (out);

      Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > viewsConst =
        X_noncontig->get2dView ();
      if (! std::is_same<typename MV::dual_view_type::t_dev::memory_space,
	  typename MV::dual_view_type::t_host::memory_space>::value) {
	TEST_ASSERT( ! X_noncontig->need_sync_host () );
      }
    }

    out << "Test X_noncontig->get2dViewNonConst()" << endl;
    {
      Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > viewsNonConst =
        X_noncontig->get2dViewNonConst ();
      if (! std::is_same<typename MV::dual_view_type::t_dev::memory_space,
	  typename MV::dual_view_type::t_host::memory_space>::value) {
	TEST_ASSERT( ! X_noncontig->need_sync_host () );
      }

      // get2dViewNonConst is supposed to mark the host data as
      // modified.  Thus, the device data need a sync, if host and
      // device are actually different data.
      //
      // The sync state of a MultiVector is only updated if the memory spaces 
      // are not the same between host and device
      if (! std::is_same<typename MV::dual_view_type::t_dev::memory_space,
                         typename MV::dual_view_type::t_host::memory_space>::value) {
        TEST_ASSERT( X_noncontig->need_sync_device () );
      }

      TEST_EQUALITY( static_cast<size_t> (viewsNonConst.size ()), numColsSubset );
      if (static_cast<size_t> (viewsNonConst.size ()) == numColsSubset) {
        for (size_t j = 0; j < numColsSubset; ++j) {
          Teuchos::ArrayRCP<Scalar> view_j = viewsNonConst[j];

          for (LO i = 0; i < lclNumRows; ++i) {
            view_j[i] = flagVal;
          }
        }
      }

      // In the Chris Baker era (<= 2012, before the Kokkos 2.0
      // refactor), Teuchos::ArrayRCP views of CUDA device data were
      // in host memory, and were "generalized views."  That is, they
      // had the right to defer updating the device data until the
      // reference count of the ArrayRCP went to zero, at which point
      // (if the view was nonconst) they would copy back.  We don't
      // retain this behavior in the implementation, but still reserve
      // the right to do that.  Thus, we need to drop the reference
      // count to zero here.
      viewsNonConst = Teuchos::null;

      // Sync back to device, so we can test whether the modifications
      // took effect.
      for (size_t j = 0; j < numColsSubset; ++j) {
        auto X_col = X_noncontig->getVector (j);
        diff.putScalar (flagVal);
        diff.update (STS::one (), *X_col, -STS::one ());
        const auto diffNorm = diff.normInf ();
        TEST_EQUALITY( diffNorm, STM::zero () );
      }
    }

    // Make sure that none of the columns of X that the above code was
    // not supposed to modify actually were modified.
    curVal = STS::one ();
    for (size_t j = 0; j < numCols; ++j) {
      bool skip = false;
      for (size_t k = 0; k < numColsSubset; ++k) {
        if (colsSubset[k] == j) {
          skip = true;
          break;
        }
      }
      if (! skip) {
        auto X_j = X.getVector (j);
        diff.putScalar (curVal);
        diff.update (STS::one (), *X_j, -STS::one ());
        const auto diffNorm = diff.normInf ();
        TEST_EQUALITY( diffNorm, STM::zero () );
      }
      curVal += STS::one ();
    }

    out << "Make sure all processes had correct results" << endl;
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));

    if (gblSuccess != 1) {
      out << "Test FAILED on some process in the communicator!" << endl;
      success = false;
    }
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, get2dView_noncontig, SCALAR, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV( UNIT_TEST_GROUP )

} // namespace (anonymous)

