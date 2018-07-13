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
  decltype (Kokkos::subview (X.template getLocalView<Kokkos::HostSpace> (), Kokkos::ALL (), 0))
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

  host_view_type X_lcl_host_1d;
  if (std::is_same<dev_memory_space, host_memory_space>::value ||
      X.template need_sync<dev_memory_space> () ||
      ! X.template need_sync<host_memory_space> ()) {
    // Can use host version of the MultiVector directly.
    auto X_lcl_host = X.template getLocalView<host_memory_space> ();
    X_lcl_host_1d = Kokkos::subview (X_lcl_host, Kokkos::ALL (), 0);
  }
  else { // need to copy from device to host (see above discussion)
    auto X_lcl_dev = X.template getLocalView<dev_memory_space> ();
    auto X_lcl_dev_1d = Kokkos::subview (X_lcl_dev, Kokkos::ALL (), 0);
    X_lcl_host_1d = Kokkos::create_mirror_view (X_lcl_dev_1d);
    Kokkos::deep_copy (X_lcl_host_1d, X_lcl_dev_1d);
  }
  return X_lcl_host_1d;
}

template<class MultiVectorType>
typename Kokkos::View<typename MultiVectorType::impl_scalar_type**,
                      typename MultiVectorType::dual_view_type::t_host::array_layout,
                      typename MultiVectorType::device_type>::HostMirror
getHostViewOfMultiVector (const MultiVectorType& X)
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

  host_view_type X_lcl_host;
  if (std::is_same<dev_memory_space, host_memory_space>::value ||
      X.template need_sync<dev_memory_space> () ||
      ! X.template need_sync<host_memory_space> ()) {
    // Can use host version of the MultiVector directly.
    X_lcl_host = X.template getLocalView<host_memory_space> ();
  }
  else { // need to copy from device to host (see above discussion)
    auto X_lcl_dev = X.template getLocalView<dev_memory_space> ();
    X_lcl_host = Kokkos::create_mirror_view (X_lcl_dev);
    Kokkos::deep_copy (X_lcl_host, X_lcl_dev);
  }
  return X_lcl_host;
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
  typedef typename MV::impl_scalar_type IST;
  typedef typename MV::dual_view_type::array_layout array_layout;
  typedef typename MV::device_type device_type;
  typedef typename device_type::memory_space dev_memory_space;
  typedef typename Kokkos::View<IST**, array_layout,
    device_type>::HostMirror::memory_space host_memory_space;

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
    X0.template sync<host_memory_space> ();
    X0.template modify<host_memory_space> ();
    auto X0_lcl = X0.template getLocalView<host_memory_space> ();
    for (LO j = 0; j < numVecs; ++j) {
      auto X0_lcl_j = Kokkos::subview (X0_lcl, Kokkos::ALL (), j);
      for (LO i = 0; i < lclNumRows; ++i) {
        X0_lcl_j(i) = static_cast<SC> (i + 0.1*j);
      }
    }
    X0.template sync<dev_memory_space> ();
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
