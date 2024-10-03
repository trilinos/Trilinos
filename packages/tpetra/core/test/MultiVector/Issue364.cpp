// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Kokkos_ArithTraits.hpp"

// CUDA doesn't like Kokkos functors that live in anonymous namespaces.
namespace TpetraTest {

template<class ViewType>
class VectorsEqual {
  static_assert (Kokkos::is_view<ViewType>::value,
                 "ViewType must be a Kokkos::View specialization.");
  static_assert (static_cast<int> (ViewType::rank) == 1,
                 "ViewType must be a 1-D Kokkos::View.");
public:
  typedef typename ViewType::size_type size_type;
  typedef int value_type; // as a bool; bool doesn't work in Kokkos

  VectorsEqual (const ViewType& x, const ViewType& y) :
    x_ (x), y_ (y)
  {}

  static bool run (const ViewType& x, const ViewType& y) {
    typedef typename ViewType::execution_space execution_space;
    typedef Kokkos::RangePolicy<execution_space, size_type> range_type;
    int result = 1;
    Kokkos::parallel_reduce (range_type (0, x.extent (0)),
                             VectorsEqual<ViewType> (x, y),
                             result);
    return result == 1;
  }

  KOKKOS_INLINE_FUNCTION void
  operator () (const size_type& i, value_type& dst) const
  {
    if (x_(i) != y_(i)) {
      dst = 0;
    }
  }

  //! Set the initial value of the reduction result.
  KOKKOS_INLINE_FUNCTION void init (value_type& dst) const
  {
    dst = 1;
  }

  KOKKOS_INLINE_FUNCTION void
  join (value_type& dst,
        const value_type& src) const
  {
    dst = (src == 1 && dst == 1) ? 1 : 0;
  }

private:
  ViewType x_;
  ViewType y_;
};

template<class ViewType>
class NegateAllEntries {
  static_assert (Kokkos::is_view<ViewType>::value,
                 "ViewType must be a Kokkos::View specialization.");
  static_assert (static_cast<int> (ViewType::rank) == 1,
                 "ViewType must be a 1-D Kokkos::View.");
public:
  typedef typename ViewType::size_type size_type;

  NegateAllEntries (const ViewType& x) : x_ (x) {}

  static void run (const ViewType& x) {
    typedef typename ViewType::execution_space execution_space;
    typedef Kokkos::RangePolicy<execution_space, size_type> range_type;
    Kokkos::parallel_for (range_type (0, x.extent (0)),
                          NegateAllEntries<ViewType> (x));
  }

  KOKKOS_INLINE_FUNCTION void operator () (const size_type& i) const
  {
    x_(i) = -x_(i);
  }

private:
  ViewType x_;
};

} // namespace TpetraTest

namespace { // (anonymous)

  // Are all entries of x and y equal?
  template<class ViewType>
  bool
  vectorsEqual (const ViewType& x, const ViewType& y)
  {
    static_assert (Kokkos::is_view<ViewType>::value,
                   "ViewType must be a Kokkos::View specialization.");
    static_assert (static_cast<int> (ViewType::rank) == 1,
                   "ViewType must be a 1-D Kokkos::View.");
    return TpetraTest::VectorsEqual<ViewType>::run (x, y);
  }

  template<class ViewType>
  void
  negateAllEntries (const ViewType& x)
  {
    static_assert (Kokkos::is_view<ViewType>::value,
                   "ViewType must be a Kokkos::View specialization.");
    static_assert (static_cast<int> (ViewType::rank) == 1,
                   "ViewType must be a 1-D Kokkos::View.");
    TpetraTest::NegateAllEntries<ViewType>::run (x);
  }

  //
  // UNIT TESTS
  //

  using std::endl;

  // These tests only take Node as a template parameter, because the
  // subview mechanism does not depend on types other than the
  // DeviceType (which is a function of Node).

  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiVector, DualViewNoncontig, Node )
  {
    using Scalar = typename Tpetra::MultiVector<>::scalar_type;
    using LO = typename Tpetra::MultiVector<>::local_ordinal_type;
    using GO = typename Tpetra::MultiVector<>::global_ordinal_type;
    using map_type = Tpetra::Map<LO, GO, Node>;
    using MV = Tpetra::MultiVector<Scalar, LO, GO, Node>;
    using IST = typename MV::impl_scalar_type;

    const IST ONE = Kokkos::ArithTraits<IST>::one ();

    out << "Test MultiVector dual view semantics with a "
      "view of a noncontiguous set of columns" << endl;
    Teuchos::OSTab tab1 (out);

    auto comm = Tpetra::TestingUtilities::getDefaultComm ();

    const LO lclNumRows = 17;
    const int numProcs = comm->getSize ();
    const GO gblNumRows =
      static_cast<GO> (lclNumRows) * static_cast<GO> (numProcs);
    const GO indexBase = 0;
    auto map = rcp (new map_type (gblNumRows, lclNumRows, indexBase, comm));
    const LO numVecs = 7;

    MV X (map, numVecs); // the MultiVector to view

    out << "Fill the master MultiVector X" << endl;
    {
      auto X_lcl = X.getLocalViewHost(Tpetra::Access::ReadWrite);

      TEST_EQUALITY( static_cast<LO> (X_lcl.extent (0)), lclNumRows );
      TEST_EQUALITY( static_cast<LO> (X_lcl.extent (1)), numVecs );
      if (! success) {
        return;
      }

      // Fill this MultiVector so that each entry has a unique value
      // (Scalar is generally big enough to hold this much data).  This
      // will help us identify whether subviews are correct.
      IST curVal = ONE;
      for (LO j = 0; j < numVecs; ++j) {
        for (LO i = 0; i < lclNumRows; ++i) {
          X_lcl(i,j) = curVal;
          curVal = curVal + Kokkos::ArithTraits<IST>::one ();
        }
      }
    }

    // Make a "backup" of X, for later comparison.
    MV X_copy (X, Teuchos::Copy);

    out << "Create a view of a noncontiguous subset of columns of X" << endl;
    Teuchos::Array<size_t> cols (4);
    cols[0] = 1;
    cols[1] = 3;
    cols[2] = 5;
    cols[3] = 6;
    Teuchos::RCP<MV> X_sub = X.subViewNonConst (cols);

    out << "Modify the entries of that view, on host" << endl;
    {
      // Use negative values to distinguish changes to the subview.
      for (LO k = 0; k < static_cast<LO> (cols.size ()); ++k) {
        const LO j = cols[k];

        TEST_ASSERT( j < static_cast<LO> (X.getNumVectors ()) );
        TEST_ASSERT( k < static_cast<LO> (X_sub->getNumVectors ()) );
        if (j >= static_cast<LO> (X.getNumVectors ()) ||
            k >= static_cast<LO> (X_sub->getNumVectors ())) {
          continue;
        }

        auto X_sub_j = X_sub->getVectorNonConst (k);
        auto X_sub_lcl_j_2d =
          X_sub_j->getLocalViewHost(Tpetra::Access::ReadWrite);
        auto X_sub_lcl_j_1d =
          Kokkos::subview (X_sub_lcl_j_2d, Kokkos::ALL (), 0);

        TEST_EQUALITY( static_cast<LO> (X_sub_lcl_j_1d.extent (0)),
                       lclNumRows );
        if (! success) {
          return;
        }
        for (LO i = 0; i < lclNumRows; ++i) {
          X_sub_lcl_j_1d(i) = -X_sub_lcl_j_1d(i);
        }
      }
    }

    // At this point, the host version of X_sub and the host version
    // of the corresponding columns of X should be the same.
    out << "Compare X_sub to X and X_copy, on host" << endl;
    {
      for (LO k = 0; k < static_cast<LO> (cols.size ()); ++k) {
        const LO j = cols[k];

        TEST_ASSERT( j < static_cast<LO> (X.getNumVectors ()) );
        TEST_ASSERT( k < static_cast<LO> (X_sub->getNumVectors ()) );
        if (j >= static_cast<LO> (X.getNumVectors ()) ||
            k >= static_cast<LO> (X_sub->getNumVectors ())) {
          continue;
        }

        auto X_sub_j = X_sub->getVector (k);
        TEST_ASSERT( ! X_sub_j.is_null () );
        if (X_sub_j.is_null ()) {
          continue;
        }
        auto X_sub_lcl_j_2d =
          X_sub_j->getLocalViewHost(Tpetra::Access::ReadOnly);
        auto X_sub_lcl_j_1d =
          Kokkos::subview (X_sub_lcl_j_2d, Kokkos::ALL (), 0);

        auto X_j = X.getVector (j);
        TEST_ASSERT( ! X_j.is_null () );
        if (X_j.is_null ()) {
          continue;
        }
        auto X_lcl_j_2d = X_j->getLocalViewHost(Tpetra::Access::ReadOnly);
        auto X_lcl_j_1d = Kokkos::subview (X_lcl_j_2d, Kokkos::ALL (), 0);

        auto X_copy_j = X_copy.getVector (j);
        TEST_ASSERT( ! X_copy_j.is_null () );
        if (X_copy_j.is_null ()) {
          continue;
        }
        auto X_copy_lcl_j_2d =
          X_copy_j->getLocalViewHost(Tpetra::Access::ReadOnly);
        auto X_copy_lcl_j_1d =
          Kokkos::subview (X_copy_lcl_j_2d, Kokkos::ALL (), 0);

        TEST_EQUALITY( static_cast<LO> (X_sub_lcl_j_1d.extent (0)),
                       lclNumRows );
        TEST_EQUALITY( static_cast<LO> (X_lcl_j_1d.extent (0)),
                       lclNumRows );
        TEST_EQUALITY( static_cast<LO> (X_copy_lcl_j_1d.extent (0)),
                       lclNumRows );
        if (! success) {
          return;
        }
        for (LO i = 0; i < lclNumRows; ++i) {
          TEST_EQUALITY( X_sub_lcl_j_1d(i), X_lcl_j_1d(i) );
          TEST_EQUALITY( X_sub_lcl_j_1d(i), -X_copy_lcl_j_1d(i) );
        }
      }
    }

    // Furthermore, columns of X that are NOT included in X_sub should
    // not have changed.
    out << "Check columns of X not included in X_sub, on host" << endl;
    {
      for (LO j = 0; j < numVecs; ++j) {
        auto iter = std::find (cols.begin (), cols.end (), j);
        if (iter == cols.end ()) { // not in the sequence
          auto X_j = X.getVector (j);
          auto X_lcl_j_2d = X_j->getLocalViewHost(Tpetra::Access::ReadOnly);
          auto X_lcl_j_1d = Kokkos::subview (X_lcl_j_2d, Kokkos::ALL (), 0);

          auto X_copy_j = X_copy.getVector (j);
          auto X_copy_lcl_j_2d =
            X_copy_j->getLocalViewHost(Tpetra::Access::ReadOnly);
          auto X_copy_lcl_j_1d =
            Kokkos::subview (X_copy_lcl_j_2d, Kokkos::ALL (), 0);

          for (LO i = 0; i < lclNumRows; ++i) {
            TEST_EQUALITY( X_lcl_j_1d(i), X_copy_lcl_j_1d(i) );
          }
        }
      }
    }

    out << "Sync X_sub to device" << endl;

    // At this point, the device version of X_sub and the device
    // version of the corresponding columns of X should be the same.
    out << "Compare X_sub to X and X_copy, on device" << endl;
    {
      for (LO k = 0; k < static_cast<LO> (cols.size ()); ++k) {
        const LO j = cols[k];
        auto X_sub_j = X_sub->getVector (k);
        auto X_sub_lcl_j_2d = X_sub_j->getLocalViewDevice(Tpetra::Access::ReadOnly);
        auto X_sub_lcl_j_1d = Kokkos::subview (X_sub_lcl_j_2d, Kokkos::ALL (), 0);

        auto X_j = X.getVector (j);
        auto X_lcl_j_2d = X_j->getLocalViewDevice(Tpetra::Access::ReadOnly);
        auto X_lcl_j_1d = Kokkos::subview (X_lcl_j_2d, Kokkos::ALL (), 0);

        // auto X_copy_j = X_copy.getVector (j);
        // auto X_copy_lcl_j_2d = X_copy_j->getLocalViewDevice ();
        // auto X_copy_lcl_j_1d = Kokkos::subview (X_copy_lcl_j_2d, Kokkos::ALL (), 0);

        const bool eq = vectorsEqual (X_sub_lcl_j_1d, X_lcl_j_1d);
        TEST_ASSERT( eq );
      }
    }

    // Furthermore, columns of X that are NOT included in X_sub should
    // not have changed.
    out << "Check columns of X not included in X_sub, on device" << endl;
    {
      for (LO j = 0; j < numVecs; ++j) {
        auto iter = std::find (cols.begin (), cols.end (), j);
        if (iter == cols.end ()) { // not in the sequence
          auto X_j = X.getVector (j);
          auto X_lcl_j_2d = X_j->getLocalViewDevice(Tpetra::Access::ReadOnly);
          auto X_lcl_j_1d = Kokkos::subview (X_lcl_j_2d, Kokkos::ALL (), 0);

          auto X_copy_j = X_copy.getVector (j);
          auto X_copy_lcl_j_2d = X_copy_j->getLocalViewDevice(Tpetra::Access::ReadOnly);
          auto X_copy_lcl_j_1d = Kokkos::subview (X_copy_lcl_j_2d, Kokkos::ALL (), 0);

          const bool eq = vectorsEqual (X_lcl_j_1d, X_copy_lcl_j_1d);
          TEST_ASSERT( eq );
        }
      }
    }
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiVector, DualViewTwoDisjoint, Node )
  {
    typedef typename Tpetra::MultiVector<>::scalar_type Scalar;
    typedef typename Tpetra::MultiVector<>::local_ordinal_type LO;
    typedef typename Tpetra::MultiVector<>::global_ordinal_type GO;
    typedef typename Tpetra::MultiVector<>::device_type DT;
    typedef typename DT::memory_space dev_memory_space;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::MultiVector<Scalar, LO, GO, Node> MV;
    typedef typename MV::impl_scalar_type IST;
    typedef typename Kokkos::View<IST**, DT>::HostMirror::memory_space
      host_memory_space;

    const IST ONE = Kokkos::ArithTraits<IST>::one ();
    const IST TWO = ONE + ONE;

    out << "Test MultiVector dual view semantics with two disjoint, "
      "noncontiguous column views of a single MultiVector" << endl;
    Teuchos::OSTab tab1 (out);
    // This test only mattters if the device memory space and the host
    // memory space are actually different.  We use an if-else rather
    // than an early return, to avoid warnings about "unreachable
    // code."
    if (std::is_same<dev_memory_space, host_memory_space>::value) {
      out << "This test only matters if the device memory space and the "
        "host memory space are actually different.  In this case, they "
        "are both the same.  (Note that the HostMirror::memory_space of "
        "a CudaUVMSpace View is also CudaUVMSpace.)  Thus, there is no "
        "point in continuing the test for the current device type." << endl;
    }
    else {
      // 1. Create and fill a MultiVector X, and make sure it is sync'd
      //    to device.
      // 2. Create views Y and Z of X, that view disjoint sets of
      //    columns of X.
      // 3. Sync Y to host, and modify it there.  Don't sync it to
      //    device.
      // 4. Modify Z on device, and sync it to host.
      // 5. Check that Z's sync to host did not cause all of X to get
      //    sync'd to host, thus clobbering Y's changes on host.

      auto comm = Tpetra::TestingUtilities::getDefaultComm ();

      const LO lclNumRows = 17;
      const int numProcs = comm->getSize ();
      const GO gblNumRows =
        static_cast<GO> (lclNumRows) * static_cast<GO> (numProcs);
      const GO indexBase = 0;
      auto map = rcp (new map_type (gblNumRows, lclNumRows, indexBase, comm));
      const LO numVecs = 7;

      MV X (map, numVecs); // the MultiVector to view

      out << "Fill the master MultiVector X" << endl;
      {
        auto X_lcl = X.getLocalViewHost(Tpetra::Access::ReadWrite);

        TEST_EQUALITY( static_cast<LO> (X_lcl.extent (0)), lclNumRows );
        TEST_EQUALITY( static_cast<LO> (X_lcl.extent (1)), numVecs );
        if (! success) {
          return;
        }

        // Fill this MultiVector so that each entry has a unique value
        // (Scalar is generally big enough to hold this much data).  This
        // will help us identify whether subviews are correct.
        IST curVal = ONE;
        for (LO j = 0; j < numVecs; ++j) {
          for (LO i = 0; i < lclNumRows; ++i) {
            X_lcl(i,j) = curVal;
            curVal = curVal + Kokkos::ArithTraits<IST>::one ();
          }
        }
      }

      // Make a "backup" of X, for later comparison.
      MV X_copy (X, Teuchos::Copy);

      out << "Create first view Y of a noncontiguous subset of columns of X" << endl;
      Teuchos::Array<size_t> Y_cols (4);
      Y_cols[0] = 1;
      Y_cols[1] = 3;
      Y_cols[2] = 5;
      Y_cols[3] = 6;
      Teuchos::RCP<MV> Y = X.subViewNonConst (Y_cols);
      TEST_EQUALITY( Y->getNumVectors (), static_cast<size_t> (4) );

      out << "Create second view Z of a noncontiguous subset of columns "
        "of X, disjoint from the columns that Y (see above) views" << endl;
      Teuchos::Array<size_t> Z_cols (3);
      Z_cols[0] = 0;
      Z_cols[1] = 2;
      Z_cols[2] = 4;
      Teuchos::RCP<MV> Z = X.subViewNonConst (Z_cols);
      TEST_EQUALITY( Z->getNumVectors (), static_cast<size_t> (3) );

      out << "Sync Y to host, and modify it there.  "
        "Don't sync it back to device." << endl;
      {
        // Multiply all entries by 2.
        for (LO k = 0; k < static_cast<LO> (Y_cols.size ()); ++k) {
          const LO j = Y_cols[k];
          TEST_ASSERT( j < static_cast<LO> (X.getNumVectors ()) );
          TEST_ASSERT( k < static_cast<LO> (Y->getNumVectors ()) );
          if (j >= static_cast<LO> (X.getNumVectors ()) ||
              k >= static_cast<LO> (Y->getNumVectors ())) {
            continue;
          }
          auto Y_j = Y->getVectorNonConst (k);
          auto Y_lcl_j_2d = Y_j->getLocalViewHost(Tpetra::Access::ReadWrite);
          auto Y_lcl_j_1d = Kokkos::subview (Y_lcl_j_2d, Kokkos::ALL (), 0);
          TEST_EQUALITY( static_cast<LO> (Y_lcl_j_1d.extent (0)), lclNumRows );
          if (! success) {
            return;
          }
          for (LO i = 0; i < lclNumRows; ++i) {
            Y_lcl_j_1d(i) = TWO * Y_lcl_j_1d(i);
          }
        }
      }

      out << "Modify Z on device, and sync it to host." << endl;
      {
        /// with multi non-contiguous memory space, syncing Z 
        /// would sync Y part (actually it sync the entire view of X).
        for (LO k = 0; k < static_cast<LO> (Z_cols.size ()); ++k) {
          const LO j = Z_cols[k];
          TEST_ASSERT( j < static_cast<LO> (X.getNumVectors ()) );
          TEST_ASSERT( k < static_cast<LO> (Z->getNumVectors ()) );
          if (j >= static_cast<LO> (X.getNumVectors ()) ||
              k >= static_cast<LO> (Z->getNumVectors ())) {
            continue;
          }
          auto Z_j = Z->getVectorNonConst (k);
          auto Z_lcl_j_2d = Z_j->getLocalViewDevice(Tpetra::Access::ReadWrite);
          auto Z_lcl_j_1d = Kokkos::subview (Z_lcl_j_2d, Kokkos::ALL (), 0);
          TEST_EQUALITY( static_cast<LO> (Z_lcl_j_1d.extent (0)), lclNumRows );
          if (! success) {
            return;
          }
          negateAllEntries (Z_lcl_j_1d);
        }

        // THIS is the thing that tests for Issue #364.
      }

      // Check that Z's sync to host did not cause all of X to get
      // sync'd to host, thus clobbering Y's changes on host.  If
      // MultiVector's dual view semantics are implemented correctly, we
      // should see that Y(:,k) == 2*X_copy(:,j), where j = Y_cols[k].

      {
        // We don't need to sync Y here; it's already sync'd to host.
        out<<"Checkpoint #1"<<std::endl;
        for (LO k = 0; k < static_cast<LO> (Y_cols.size ()); ++k) {
          const LO j = Y_cols[k];
          TEST_ASSERT( j < static_cast<LO> (X.getNumVectors ()) );
          TEST_ASSERT( j < static_cast<LO> (X_copy.getNumVectors ()) );
          TEST_ASSERT( k < static_cast<LO> (Y->getNumVectors ()) );
          if (j >= static_cast<LO> (X.getNumVectors ()) ||
              k >= static_cast<LO> (Y->getNumVectors ())) {
            continue;
          }
          out<<"Checkpoint #2"<<std::endl;
          auto X_j = X.getVector (j);
          TEST_ASSERT( ! X_j.is_null () );
          if (X_j.is_null ()) {
            continue;
          }
          out<<"Checkpoint #3"<<std::endl;
          auto X_lcl_j_2d = X_j->getLocalViewHost(Tpetra::Access::ReadOnly);
          auto X_lcl_j_1d = Kokkos::subview (X_lcl_j_2d, Kokkos::ALL (), 0);
          TEST_EQUALITY( static_cast<LO> (X_lcl_j_1d.extent (0)),
                         lclNumRows );

          out<<"Checkpoint #4"<<std::endl;
          auto X_copy_j = X_copy.getVector (j);
          if (X_copy_j.is_null ()) {
            continue;
          }
          auto X_copy_lcl_j_2d = X_copy_j->getLocalViewHost(Tpetra::Access::ReadOnly);
          auto X_copy_lcl_j_1d = Kokkos::subview (X_copy_lcl_j_2d, Kokkos::ALL (), 0);
          TEST_EQUALITY( static_cast<LO> (X_copy_lcl_j_1d.extent (0)),
                         lclNumRows );

          auto Y_j = Y->getVector (k);
          if (Y_j.is_null ()) {
            continue;
          }
          out<<"Checkpoint #5"<<std::endl;
          auto Y_lcl_j_2d = Y_j->getLocalViewHost(Tpetra::Access::ReadOnly);
          auto Y_lcl_j_1d = Kokkos::subview (Y_lcl_j_2d, Kokkos::ALL (), 0);
          TEST_EQUALITY( static_cast<LO> (Y_lcl_j_1d.extent (0)),
                         lclNumRows );
          if (! success) {
            return;
          }
          out<<"Checkpoint #6"<<std::endl;
          for (LO i = 0; i < lclNumRows; ++i) {
            TEST_EQUALITY( Y_lcl_j_1d(i), X_lcl_j_1d(i) );
            TEST_EQUALITY( Y_lcl_j_1d(i), TWO * X_copy_lcl_j_1d(i) );
          }
        }
      }
    }
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiVector, DualViewNoncontig, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MultiVector, DualViewTwoDisjoint, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_N( UNIT_TEST_GROUP )

} // namespace (anonymous)

