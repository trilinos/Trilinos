// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Map.hpp"
#include "Tpetra_TestingUtilities.hpp"

// Bug 6437 test
//
// Exercise Tpetra::Map::getLocalMap().  Make sure that the returned
// "local Map" object behaves like Tpetra::Map, for the methods that
// the two have in common.

namespace { // (anonymous)

using Tpetra::TestingUtilities::getDefaultComm;
using Teuchos::Comm;
using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::REDUCE_MIN;
using Teuchos::reduceAll;
using std::endl;
typedef Tpetra::global_size_t GST;

//
// Test LocalMap with a uniform contiguous Map.
//
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( LocalMap, UniformContig, LO, GO, NT )
{
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef typename NT::execution_space execution_space;
  typedef Kokkos::RangePolicy<execution_space, int> range_type;

  int lclSuccess = 1;
  int gblSuccess = 1;

  Teuchos::OSTab tab0 (out);
  out << "LocalMap: Uniform Contiguous Map" << endl;
  Teuchos::OSTab tab1 (out);

  RCP<const Comm<int> > comm = getDefaultComm ();
  const LO numLclElts = 10;

  const GST numGblElts = static_cast<GST> (comm->getSize () * numLclElts);
  const GO indexBase = 0;

  map_type gblMap (numGblElts, indexBase, comm);
  auto lclMap = gblMap.getLocalMap ();

  try {
    // Sanity check on the global Map.
    TEST_EQUALITY( gblMap.getLocalNumElements (), static_cast<size_t> (numLclElts) );
    TEST_EQUALITY( gblMap.getMinLocalIndex (), static_cast<LO> (0) );
    TEST_EQUALITY( gblMap.getMaxLocalIndex (), static_cast<LO> (numLclElts - 1) );

    // Test constants.
    TEST_EQUALITY( gblMap.getLocalNumElements (), static_cast<size_t> (lclMap.getLocalNumElements ()) );
    TEST_EQUALITY( gblMap.getIndexBase (), lclMap.getIndexBase () );
    TEST_EQUALITY( gblMap.isContiguous (), lclMap.isContiguous () );
    TEST_EQUALITY( gblMap.getMinLocalIndex (), lclMap.getMinLocalIndex () );
    TEST_EQUALITY( gblMap.getMaxLocalIndex (), lclMap.getMaxLocalIndex () );
    TEST_EQUALITY( gblMap.getMinGlobalIndex (), lclMap.getMinGlobalIndex () );
    TEST_EQUALITY( gblMap.getMaxGlobalIndex (), lclMap.getMaxGlobalIndex () );

    // Test global-to-local and local-to-global index conversion.
    if (numLclElts > 0) {
      for (LO lid = 0; lid < numLclElts; ++lid) {
        const GO expectedGid = gblMap.getGlobalElement (lid);
        // Make sure that the (global) Map behaves as expected.
        TEST_INEQUALITY( expectedGid, Tpetra::Details::OrdinalTraits<GO>::invalid () );

        // gblMap.getLocalElement is host
        TEST_EQUALITY( gblMap.getLocalElement (expectedGid), lid );

        // lclMap.getGlobalElement is device only
        GO globalElement;
        Kokkos::parallel_reduce("read GO element", range_type (0, 1),
          KOKKOS_LAMBDA(int dummy, GO& ge) {
            ge = lclMap.getGlobalElement(lid);
          }, globalElement);
        TEST_EQUALITY( globalElement, expectedGid );

        // lclMap.getLocalElement is device only
        LO localElement;
        Kokkos::parallel_reduce("read LO element", range_type (0, 1),
          KOKKOS_LAMBDA(int dummy, LO& le) {
            le = lclMap.getLocalElement(expectedGid);
          }, localElement);
        TEST_EQUALITY( localElement, lid );
      }
    }
  } catch (...) {
    // None of the above methods should throw exceptions, but we catch
    // just in case.
    success = false;
  }

  // Make sure that all processes succeeded.
  lclSuccess = success ? 1 : 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY_CONST( gblSuccess, 1 );

  if (gblSuccess != 1) {
    out << "The comparison test between Map and LocalMap "
      "failed on at least one process." << endl;
    // Make sure that the test fails, even if Process 0 was fine.
    success = false;
    return;
  }
}

//
// Test LocalMap with a NONuniform contiguous Map.
//
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( LocalMap, NonuniformContig, LO, GO, NT )
{
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef typename NT::execution_space execution_space;
  typedef Kokkos::RangePolicy<execution_space, int> range_type;

  int lclSuccess = 1;
  int gblSuccess = 1;

  Teuchos::OSTab tab0 (out);
  out << "LocalMap: Nonuniform Contiguous Map" << endl;
  Teuchos::OSTab tab1 (out);

  RCP<const Comm<int> > comm = getDefaultComm ();
  const int myRank = comm->getRank ();

  // Process p gets 5+p indices.
  const LO numLclElts = static_cast<LO> (5 + myRank);
  const GST numGblElts = Tpetra::Details::OrdinalTraits<GST>::invalid ();
  const GO indexBase = 0;

  map_type gblMap (numGblElts, static_cast<size_t> (numLclElts), indexBase, comm);
  auto lclMap = gblMap.getLocalMap ();

  try {
    // Sanity check on the global Map.
    TEST_EQUALITY( gblMap.getLocalNumElements (), static_cast<size_t> (numLclElts) );
    TEST_EQUALITY( gblMap.getMinLocalIndex (), static_cast<LO> (0) );
    TEST_EQUALITY( gblMap.getMaxLocalIndex (), static_cast<LO> (numLclElts - 1) );

    // Test constants.
    TEST_EQUALITY( gblMap.getLocalNumElements (), static_cast<size_t> (lclMap.getLocalNumElements ()) );
    TEST_EQUALITY( gblMap.getIndexBase (), lclMap.getIndexBase () );
    TEST_EQUALITY( gblMap.isContiguous (), lclMap.isContiguous () );
    TEST_EQUALITY( gblMap.getMinLocalIndex (), lclMap.getMinLocalIndex () );
    TEST_EQUALITY( gblMap.getMaxLocalIndex (), lclMap.getMaxLocalIndex () );
    TEST_EQUALITY( gblMap.getMinGlobalIndex (), lclMap.getMinGlobalIndex () );
    TEST_EQUALITY( gblMap.getMaxGlobalIndex (), lclMap.getMaxGlobalIndex () );

    // Test global-to-local and local-to-global index conversion.
    if (numLclElts > 0) {
      for (LO lid = 0; lid < numLclElts; ++lid) {
        const GO expectedGid = gblMap.getGlobalElement (lid);

        // lclMap.getGlobalElement is device only
        GO globalElement;
        Kokkos::parallel_reduce("read GO element", range_type (0, 1),
          KOKKOS_LAMBDA(int dummy, GO& ge) {
            ge = lclMap.getGlobalElement(lid);
          }, globalElement);
        TEST_EQUALITY( globalElement, expectedGid );

        // lclMap.getLocalElement is device only
        LO localElement;
        Kokkos::parallel_reduce("read LO element", range_type (0, 1),
          KOKKOS_LAMBDA(int dummy, LO& le) {
            le = lclMap.getLocalElement(expectedGid);
          }, localElement);
        TEST_EQUALITY( localElement, lid );
      }
    }
  } catch (...) {
    // None of the above methods should throw exceptions, but we catch
    // just in case.
    success = false;
  }

  // Make sure that all processes succeeded.
  lclSuccess = success ? 1 : 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY_CONST( gblSuccess, 1 );

  if (gblSuccess != 1) {
    out << "The comparison test between Map and LocalMap "
      "failed on at least one process." << endl;
    // Make sure that the test fails, even if Process 0 was fine.
    success = false;
    return;
  }
}

//
// Test LocalMap with a NONcontiguous Map.
//
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( LocalMap, Noncontig, LO, GO, NT )
{
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef typename NT::execution_space execution_space;
  typedef Kokkos::RangePolicy<execution_space, int> range_type;

  int lclSuccess = 1;
  int gblSuccess = 1;

  Teuchos::OSTab tab0 (out);
  out << "LocalMap: Noncontiguous Map" << endl;
  Teuchos::OSTab tab1 (out);

  RCP<const Comm<int> > comm = getDefaultComm ();
  const int myRank = comm->getRank ();

  // Process p gets 5 indices.
  const LO numLclElts = static_cast<LO> (5);
  const GST numGblElts = Tpetra::Details::OrdinalTraits<GST>::invalid ();
  const GO indexBase = 0;

  Teuchos::Array<GO> eltList (numLclElts);
  const GO myGidStart = static_cast<GO> ((myRank + 1) * numLclElts  - 1);
  for (LO lid = 0; lid < numLclElts; ++lid) {
    eltList[lid] = static_cast<GO> (myGidStart - lid);
  }

  map_type gblMap (numGblElts, eltList (), indexBase, comm);
  auto lclMap = gblMap.getLocalMap ();

  try {
    // Sanity check on the global Map.
    TEST_EQUALITY( gblMap.getLocalNumElements (), static_cast<size_t> (numLclElts) );
    TEST_EQUALITY( gblMap.getMinLocalIndex (), static_cast<LO> (0) );
    TEST_EQUALITY( gblMap.getMaxLocalIndex (), static_cast<LO> (numLclElts - 1) );

    // Test constants.
    TEST_EQUALITY( gblMap.getLocalNumElements (), static_cast<size_t> (lclMap.getLocalNumElements ()) );
    TEST_EQUALITY( gblMap.getIndexBase (), lclMap.getIndexBase () );
    TEST_EQUALITY( gblMap.isContiguous (), lclMap.isContiguous () );
    TEST_EQUALITY( gblMap.getMinLocalIndex (), lclMap.getMinLocalIndex () );
    TEST_EQUALITY( gblMap.getMaxLocalIndex (), lclMap.getMaxLocalIndex () );
    TEST_EQUALITY( gblMap.getMinGlobalIndex (), lclMap.getMinGlobalIndex () );
    TEST_EQUALITY( gblMap.getMaxGlobalIndex (), lclMap.getMaxGlobalIndex () );

    // Test global-to-local and local-to-global index conversion.
    if (numLclElts > 0) {
      for (LO lid = 0; lid < numLclElts; ++lid) {
        const GO expectedGid = gblMap.getGlobalElement (lid);

        // lclMap.getGlobalElement is device only
        GO globalElement;
        Kokkos::parallel_reduce("read GO element", range_type (0, 1),
          KOKKOS_LAMBDA(int dummy, GO& ge) {
            ge = lclMap.getGlobalElement(lid);
          }, globalElement);
        TEST_EQUALITY( globalElement, expectedGid );

        // lclMap.getLocalElement is device only
        LO localElement;
        Kokkos::parallel_reduce("read LO element", range_type (0, 1),
          KOKKOS_LAMBDA(int dummy, LO& le) {
            le = lclMap.getLocalElement(expectedGid);
          }, localElement);
        TEST_EQUALITY( localElement, lid );
      }
    }
  } catch (...) {
    // None of the above methods should throw exceptions, but we catch
    // just in case.
    success = false;
  }

  // Make sure that all processes succeeded.
  lclSuccess = success ? 1 : 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY_CONST( gblSuccess, 1 );

  if (gblSuccess != 1) {
    out << "The comparison test between Map and LocalMap "
      "failed on at least one process." << endl;
    // Make sure that the test fails, even if Process 0 was fine.
    success = false;
    return;
  }
}

/**
 * @test This test ensures that the fix brought by trilinos/Trilinos#11218 is tested.
 *
 * Mainly, it creates a dual view of local maps, assigns on the host view and syncs.
 * Without the fix from trilinos/Trilinos#11218, the assignment would crash at runtime
 * when running with Cuda UVM.
 */
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( LocalMap, KokkosView, LO, GO, NT )
{
    using execution_space = typename NT::execution_space;
    using map_t           = Tpetra::Map<LO, GO, NT>;
    using local_map_t     = typename map_t::local_map_type;
    using dual_view_t     = Kokkos::DualView<local_map_t*, execution_space>;

    dual_view_t my_dual_view("test view with local maps",1);

    my_dual_view.h_view(0) = local_map_t();

    my_dual_view.sync_device();
}

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( LocalMap, UniformContig, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( LocalMap, NonuniformContig, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( LocalMap, Noncontig, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( LocalMap, KokkosView, LO, GO, NT )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LGN(UNIT_TEST_GROUP)

} // namespace (anonymous)


