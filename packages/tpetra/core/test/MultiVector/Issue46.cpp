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

// Test for Issue #46 kindly contributed by Jonas Thies, and modified
// by Mark Hoemmen to work in the Teuchos unit test framework.

namespace { // (anonymous)

// FIXME (mfh 29 Sep 2016) Return type is only correct if LayoutLeft.
template<class VectorType>
auto getHostViewOfVector (const VectorType& X) ->
  decltype (Kokkos::subview (X.getLocalViewHost (Tpetra::Access::ReadOnly), Kokkos::ALL (), 0))
{
  // Don't change the sync status of the input Vector.  If it needs
  // sync to device, then the host version is the current one, so we
  // can use it without changes.  Else, if it needs sync to host, then
  // the device version is the current one, so we need to copy that to
  // host (don't sync the input!).  Otherwise, if it doesn't need sync
  // at all, we can use the host version.

  typedef typename VectorType::impl_scalar_type IST;
  typedef typename VectorType::dual_view_type::array_layout array_layout;
  typedef typename VectorType::device_type device_type;
  typedef typename device_type::memory_space dev_memory_space;
  typedef typename Kokkos::View<IST*, array_layout,
    device_type>::HostMirror host_view_type;
  typedef typename host_view_type::memory_space host_memory_space;

  if (std::is_same<dev_memory_space, host_memory_space>::value ||
      X.template need_sync<dev_memory_space> () ||
      ! X.template need_sync<host_memory_space> ()) {
    // Can use host version of the MultiVector directly.
    auto X_lcl_host = X.getLocalViewHost(Tpetra::Access::ReadOnly);
    return Kokkos::subview (X_lcl_host, Kokkos::ALL (), 0);
  }
  else { // need to copy from device to host (see above discussion)
    auto X_lcl_dev = X.getLocalViewDevice (Tpetra::Access::ReadOnly);
    auto X_lcl_dev_1d = Kokkos::subview (X_lcl_dev, Kokkos::ALL (), 0);
    host_view_type X_lcl_host_1d;
    X_lcl_host_1d = Kokkos::create_mirror_view (X_lcl_dev_1d);
    Kokkos::deep_copy (X_lcl_host_1d, X_lcl_dev_1d);
    return X_lcl_host_1d;
  }
}

template<class MultiVectorType>
auto getHostViewOfMultiVector (const MultiVectorType& X) ->
  decltype (X.getLocalViewHost (Tpetra::Access::ReadOnly))
{
  // Don't change the sync status of the input MultiVector.  If it
  // needs sync to device, then the host version is the current one,
  // so we can use it without changes.  Else, if it needs sync to
  // host, then the device version is the current one, so we need to
  // copy that to host (don't sync the input!).  Otherwise, if it
  // doesn't need sync at all, we can use the host version.

  typedef typename MultiVectorType::impl_scalar_type IST;
  typedef typename MultiVectorType::dual_view_type::array_layout array_layout;
  typedef typename MultiVectorType::device_type device_type;
  typedef typename device_type::memory_space dev_memory_space;
  typedef typename Kokkos::View<IST**, array_layout,
    device_type>::HostMirror::memory_space host_memory_space;
  typedef typename MultiVectorType::dual_view_type::t_host host_view_type;

  if (std::is_same<dev_memory_space, host_memory_space>::value ||
      X.template need_sync<dev_memory_space> () ||
      ! X.template need_sync<host_memory_space> ()) {
    // Can use host version of the MultiVector directly.
    return X.getLocalViewHost (Tpetra::Access::ReadOnly);
  }
  else { // need to copy from device to host (see above discussion)
    host_view_type X_lcl_host;
    auto X_lcl_dev = X.getLocalViewDevice (Tpetra::Access::ReadOnly);
    X_lcl_host = Kokkos::create_mirror_view (X_lcl_dev);
    Kokkos::deep_copy (X_lcl_host, X_lcl_dev);
    return X_lcl_host;
  }
}


template<class MultiVectorType>
bool
multiVectorsEqual (const MultiVectorType& X, const MultiVectorType& Y)
{
  typedef typename MultiVectorType::local_ordinal_type LO;

  const LO lclNumRows = static_cast<LO> (X.getLocalLength ());
  if (lclNumRows != static_cast<LO> (Y.getLocalLength ())) {
    return false;
  }
  const LO numVecs = static_cast<LO> (X.getNumVectors ());
  if (numVecs != static_cast<LO> (Y.getNumVectors ())) {
    return false;
  }

  // The one-column-at-a-time approach works even if either X or Y is
  // a view of noncontiguous columns of some other MultiVector.
  for (LO j = 0; j < numVecs; ++j) {
    auto X_j = X.getVector (j);
    auto Y_j = Y.getVector (j);
    auto X_j_host = getHostViewOfVector (*X_j);
    auto Y_j_host = getHostViewOfVector (*Y_j);

    for (LO i = 0; i < lclNumRows; ++i) {
      if (X_j_host(i) != Y_j_host(i)) {
        return false;
      }
    }
  }
  return true; // all entries equal
}

template<class MultiVectorType>
void
printMultiVector (std::ostream& out, const MultiVectorType& X)
{
  typedef typename MultiVectorType::local_ordinal_type LO;

  auto X_lcl_host = getHostViewOfMultiVector (X);
  const LO lclNumRows = static_cast<LO> (X.getLocalLength ());
  const LO numCols = static_cast<LO> (X.getNumVectors ());
  for (LO i = 0; i < lclNumRows; ++i) {
    for (LO j = 0; j < numCols; ++j) {
      out << std::setprecision (3)
          << std::setw (4)
          << X_lcl_host(i,j);
    }
    out << std::endl;
  }
}

// Work-around for #680: GCC 4.7.2 is broken; it does not let me say
// "X.template sync<Kokkos::HostSpace>()", even when X's type does not
// depend on template parameters.  However, other compilers require
// "template" there.  Thus, to work around this GCC 4.7.2 issue, I'll
// make the test templated on MultiVectorType, and fill in the
// template parameter when I use the Teuchos macro TEUCHOS_UNIT_TEST
// to define the unit test.
template<class MultiVectorType>
void issue46Test (bool& success, Teuchos::FancyOStream& out)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::endl;
  typedef MultiVectorType MV;
  typedef typename MV::scalar_type SC;
  typedef typename MV::local_ordinal_type LO;
  typedef typename MV::global_ordinal_type GO;
  typedef typename MV::node_type NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;

  out << "Test Github Issue #46 (offset view of an offset view of a "
    "Tpetra::MultiVector)" << endl;
  Teuchos::OSTab tab1 (out);
  auto comm = Tpetra::TestingUtilities::getDefaultComm ();

  const LO lclNumRows = 10;
  const GO indexBase = 0;
  const LO numVecs = 4;
  RCP<const map_type> map =
    rcp (new map_type (lclNumRows, indexBase, comm, Tpetra::LocallyReplicated));
  MV X0 (map, numVecs);

  // Fill the "parent" MultiVector X0 with entries, such that each
  // entry is unique, and a human can look at its value and tell where
  // it belongs.
  out << "Fill \"parent\" MultiVector X0 with entries" << endl;
  {
    auto X0_lcl = X0.getLocalViewHost (Tpetra::Access::OverwriteAll);
    for (LO j = 0; j < numVecs; ++j) {
      auto X0_lcl_j = Kokkos::subview (X0_lcl, Kokkos::ALL (), j);
      for (LO i = 0; i < lclNumRows; ++i) {
        X0_lcl_j(i) = static_cast<SC> (i + 0.1*j);
      }
    }
  }

  // Create X1, a view of X0 with row offset=offset1
  const size_t offset1 = 2;
  out << "Create X1, a view of X0 with row offset = " << offset1 << endl;
  const LO lclNumRows1 = 6;
  RCP<const map_type> map1 =
    rcp (new map_type (lclNumRows1, indexBase, comm, Tpetra::LocallyReplicated));
  RCP<MV> X1 = X0.offsetViewNonConst (map1, offset1);

  // Create X2, a view of X1 with row offset=offset2
  const size_t offset2 = 3;
  out << "Create X2, a view of X1 with row offset = " << offset2 << endl;
  const LO lclNumRows2 = 2;
  RCP<const map_type> map2 =
    rcp (new map_type (lclNumRows2, 0, comm, Tpetra::LocallyReplicated));
  RCP<MV> X2 = X1->offsetViewNonConst (map2, offset2);

  // Create X3, a view of X0 with row offset = offset1 + offset2.
  // It should thus have the same values as X2.
  const size_t offset3 = offset1 + offset2;
  out << "Create X3, a view of X0 with row offset = " << offset3 << endl;
  const LO lclNumRows3 = lclNumRows2;
  RCP<const map_type> map3 =
    rcp (new map_type (lclNumRows3, 0, comm, Tpetra::LocallyReplicated));
  RCP<MV> X3 = X0.offsetViewNonConst (map3, offset3);

  out << "Original MultiVector X0:" << endl;
  printMultiVector (out, X0);

  out << "X1, which is a view of X0 with row offset=" << offset1 << "" << endl;
  printMultiVector (out, *X1);

  out << "X2, which is a view of X1 with offset=" << offset2 << ", and thus "
    "should equal X3, a view of X0 with offset=" << (offset1 + offset2) << endl;
  printMultiVector (out, *X2);

  out << "X3, a view of x0 with offset=" << (offset1 + offset2) << endl;
  printMultiVector (out, *X3);

  out << "Check that X2 == X3" << endl;
  const bool X2_eq_X3 = multiVectorsEqual (*X2, *X3);
  TEST_ASSERT( X2_eq_X3 );
}

TEUCHOS_UNIT_TEST( MultiVector, Issue46 )
{
  typedef Tpetra::MultiVector<> MV;
  issue46Test<MV> (success, out);
}


} // namespace (anonymous)
