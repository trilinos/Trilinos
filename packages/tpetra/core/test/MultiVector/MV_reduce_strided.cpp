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
#include "Teuchos_CommHelpers.hpp"
//#include <type_traits>

namespace { // (anonymous)

template<class MapType>
Teuchos::RCP<const MapType>
makeLocalMap (const MapType& inputMap,
              const typename MapType::local_ordinal_type lclNumRows)
{
  Teuchos::RCP<MapType> map (new MapType (lclNumRows,
                                          inputMap.getIndexBase (),
                                          inputMap.getComm (),
                                          Tpetra::LocallyReplicated));
  return Teuchos::rcp_const_cast<const MapType> (map);
}

template<class MV>
bool multiVectorsLocallyEqual (Teuchos::FancyOStream& out,
                               const MV& X, const MV& Y)
{
  using std::endl;
  out << "Test MultiVector local equality" << endl;
  Teuchos::OSTab tab1 (out);

  const size_t lclNumRows = X.getLocalLength ();
  if (Y.getLocalLength () != lclNumRows) {
    return false;
  }
  const size_t numCols = X.getNumVectors ();
  if (Y.getNumVectors () != numCols) {
    return false;
  }

  out << "Dimensions match" << endl;

  // We don't want to change the user's sync state, so if we find
  // ourselves needing to sync to host, we make a deep copy first.
  MV X_copy;
  MV Y_copy;
  if (X.need_sync_host ()) {
    X_copy = MV (X, Teuchos::Copy);
  }
  else {
    X_copy = X; // harmless shallow copy
  }
  if (Y.need_sync_host ()) {
    Y_copy = MV (Y, Teuchos::Copy);
  }
  else {
    Y_copy = Y; // harmless shallow copy
  }

  // Comparing a vector at a time avoids issues with noncontiguous MVs.
  for (size_t j = 0; j < numCols; ++j) {
    auto X_j = X_copy.getVector (j);
    auto Y_j = Y_copy.getVector (j);
    auto X_j_lcl_2d = X_j->getLocalViewHost(Tpetra::Access::ReadOnly);
    auto Y_j_lcl_2d = Y_j->getLocalViewHost(Tpetra::Access::ReadOnly);
    auto X_j_lcl = Kokkos::subview (X_j_lcl_2d, Kokkos::ALL (), 0);
    auto Y_j_lcl = Kokkos::subview (Y_j_lcl_2d, Kokkos::ALL (), 0);

    for (size_t i = 0; i < lclNumRows; ++i) {
      if (X_j_lcl(i) != Y_j_lcl(i)) {
        out << "Oh no! X(" << i << "," << j << ")=" << X_j_lcl(i) << " != Y(" << i << "," << j << ")=" << Y_j_lcl(i) << endl;
        return false;
      }
    }
  }

  out << "Values match" << endl;
  return true;
}

template<class MV>
bool multiVectorsEqual (Teuchos::FancyOStream& out,
                        const MV& X, const MV& Y)
{
  using std::endl;
  out << "Test MultiVector global equality" << endl;
  Teuchos::OSTab tab1 (out);

  const auto& X_map = * (X.getMap ());
  const auto& Y_map = * (Y.getMap ());
  if (! X_map.isSameAs (Y_map)) {
    return false;
  }

  const bool lclEqual = multiVectorsLocallyEqual (out, X, Y);
  if (lclEqual) {
    out << "X and Y are locally equal" << endl;
  }
  else {
    out << "X and Y are NOT locally equal" << endl;
  }
  const int lclEq = lclEqual ? 1 : 0;
  int gblEq = 0;

  using Teuchos::outArg;
  using Teuchos::reduceAll;
  using Teuchos::REDUCE_MIN;
  const auto& comm = * (X_map.getComm ());
  reduceAll<int, int> (comm, REDUCE_MIN, lclEq, outArg (gblEq));

  if (gblEq == 1) {
    out << "X and Y are globally equal" << endl;
  }
  else {
    out << "X and Y are NOT globally equal" << endl;
  }
  return gblEq == 1;
}

template<class MV>
void reduceMultiVector (MV& Z)
{
  // Make a non-strided MultiVector, reduce on it, and copy back.
  MV Z2 (Z, Teuchos::Copy);
  Z2.reduce ();
  Tpetra::deep_copy (Z, Z2);
}

template<class MV>
void reduceMultiVector2 (MV& Z)
{
  using Teuchos::outArg;
  using Teuchos::reduceAll;
  using Teuchos::REDUCE_SUM;

  const size_t numRows = Z.getLocalLength ();
  const size_t numCols = Z.getNumVectors ();
  MV Z2 (Z.getMap (), Z.getNumVectors ());
  const auto& comm = * (Z.getMap ()->getComm ());

  for (size_t j = 0; j < numCols; ++j) {
    auto Z_j = Z.getVectorNonConst (j);
    auto Z_j_lcl_2d = Z_j->getLocalViewHost(Tpetra::Access::OverwriteAll);
    auto Z_j_lcl = Kokkos::subview (Z_j_lcl_2d, Kokkos::ALL (), 0);

    auto Z2_j = Z2.getVectorNonConst (j);
    auto Z2_j_lcl_2d = Z2_j->getLocalViewHost(Tpetra::Access::ReadWrite);
    auto Z2_j_lcl = Kokkos::subview (Z2_j_lcl_2d, Kokkos::ALL (), 0);

    reduceAll (comm, REDUCE_SUM, static_cast<int> (numRows),
               Z_j_lcl.data (), Z2_j_lcl.data ());
    Kokkos::deep_copy (Z_j_lcl, Z2_j_lcl);
  }
}

//
// UNIT TESTS
//

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, reduce_strided, Scalar, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::endl;
  using SC = Scalar;
  using LO = LocalOrdinal;
  using GO = GlobalOrdinal;
  using NT = Node;
  using MV = Tpetra::MultiVector<SC, LO, GO, NT>;
  using impl_scalar_type = typename MV::impl_scalar_type;
  using map_type = typename MV::map_type;
  using STS = Teuchos::ScalarTraits<SC>;
  using mag_type = typename MV::mag_type;
  using dual_view_type = typename MV::dual_view_type;
  using pair_type = std::pair<size_t, size_t>;

  out << "Test Tpetra::MultiVector::reduce "
    "where stride is greater than local number of rows" << endl;
  Teuchos::OSTab tab1 (out);

  const auto comm = Tpetra::TestingUtilities::getDefaultComm ();
  constexpr LO lclNumRows = 5;
  constexpr LO numCols = 2;
  constexpr GO indexBase = 0;
  const GO gblNumRows = GO (comm->getSize ()) * GO (lclNumRows);

  RCP<const map_type> map =
    rcp (new map_type (gblNumRows, lclNumRows, indexBase, comm));
  const auto lclMap = makeLocalMap (*map, lclNumRows);

  MV Z0 (lclMap, numCols);

  // Fill each column of Z0 with a different number.  Start with 1
  // instead of 0, because 0 is the default fill value.
  {
    SC curVal = STS::one ();
    for (size_t j = 0; j < numCols; ++j) {
      Z0.getVectorNonConst (j)->putScalar (curVal);
      curVal += STS::one ();
    }
  }

  MV Z1 (Z0, Teuchos::Copy);
  Z1.reduce ();
  {

    auto Z1_lcl = Z1.getLocalViewHost(Tpetra::Access::ReadOnly);
    bool reduce_expected_result = true;
    for (size_t j = 0; j < numCols; ++j) {
      SC expectedVal = SC (mag_type (j+1)) *
        SC (mag_type (comm->getSize ()));
      for (size_t i = 0; i < lclNumRows; ++i) {
        if (STS::magnitude (expectedVal - SC (Z1_lcl(i,j))) >= mag_type (lclNumRows) * STS::eps ()) {
          reduce_expected_result = false;
        }
      }
    }
    TEST_ASSERT( reduce_expected_result );

  }
  Z1.describe (out, Teuchos::VERB_EXTREME);

  MV Z2 (Z0, Teuchos::Copy);
  reduceMultiVector (Z2);

  const bool Z1_Z2_equal = multiVectorsEqual (out, Z1, Z2);
  TEST_ASSERT( Z1_Z2_equal );

  MV Z3 (Z0, Teuchos::Copy);
  reduceMultiVector2 (Z3);

  const bool Z1_Z3_equal = multiVectorsEqual (out, Z1, Z3);
  TEST_ASSERT( Z1_Z3_equal );

  {
    // Make sure Z4 has a stride greater than its number of rows.
    const size_t Z4_stride = static_cast<size_t> (lclNumRows + 31);
    dual_view_type Z4_dv_extra ("Z4", Z4_stride, numCols);
    auto Z4_dv = Kokkos::subview (Z4_dv_extra,
                                  pair_type (0, lclNumRows),
                                  pair_type (0, numCols));
    TEST_ASSERT( LO (Z4_dv.extent (0) == lclNumRows ) );
    TEST_ASSERT( LO (Z4_dv.extent (1) == numCols ) );
    TEST_ASSERT( LO (Z4_dv.d_view.extent (0) == lclNumRows ) );
    TEST_ASSERT( LO (Z4_dv.d_view.extent (1) == numCols ) );
    TEST_ASSERT( LO (Z4_dv.h_view.extent (0) == lclNumRows ) );
    TEST_ASSERT( LO (Z4_dv.h_view.extent (1) == numCols ) );
    // Kokkos could in theory insert padding in the row dimension.
    TEST_ASSERT( size_t (Z4_dv.d_view.stride (1)) >= Z4_stride );
    TEST_ASSERT( size_t (Z4_dv.h_view.stride (1)) >= Z4_stride );

    MV Z4 (lclMap, Z4_dv, Z4_dv_extra);
    TEST_ASSERT( Z4_dv.d_view.data () == Z4.getLocalViewDevice(Tpetra::Access::ReadOnly).data () );
    TEST_ASSERT( Z4_dv.h_view.data () == Z4.getLocalViewHost(Tpetra::Access::ReadOnly).data () );
    TEST_ASSERT( Z4.isConstantStride () );
    if (Z4.isConstantStride ()) {
      TEST_ASSERT( size_t (Z4_dv.d_view.stride (1)) == Z4.getStride () );
      TEST_ASSERT( size_t (Z4_dv.h_view.stride (1)) == Z4.getStride () );
      // Kokkos could in theory insert padding in the row dimension.
      TEST_ASSERT( Z4.getStride () >= Z4_stride );
    }
    for (size_t j = 0; j < numCols; ++j) {
      out << "Column j=" << j << " of Z4" << endl;
      Teuchos::OSTab colTab (out);

      auto Z4_dv_j = Kokkos::subview (Z4_dv, Kokkos::ALL (), j);
      const impl_scalar_type* Z4_h_raw = nullptr;
      const impl_scalar_type* Z4_d_raw = nullptr;
      {
        auto Z4_j_h = Z4.getVectorNonConst (j);
        auto Z4_j_h_lcl_2d = Z4_j_h->getLocalViewHost(Tpetra::Access::OverwriteAll);
        auto Z4_j_h_lcl = Kokkos::subview (Z4_j_h_lcl_2d, Kokkos::ALL (), 0);
        Z4_h_raw = Z4_j_h_lcl.data();

        TEST_ASSERT( Z4_j_h_lcl.data () == Z4_dv_j.h_view.data () );
        if (Z4_j_h_lcl.data () != Z4_dv_j.h_view.data ()) {
          out << "Z4_j_h_lcl.data() = " << Z4_j_h_lcl.data ()
              << ", Z4_dv_j.h_view.data() = " << Z4_dv_j.h_view.data ()
              << endl;
        }
        TEST_ASSERT( Z4_j_h_lcl.extent (0) == Z4_dv_j.h_view.extent (0) );
      }

      {
        auto Z4_j_d = Z4.getVectorNonConst (j);
        auto Z4_j_d_lcl_2d = Z4_j_d->getLocalViewDevice(Tpetra::Access::ReadWrite);
        auto Z4_j_d_lcl = Kokkos::subview (Z4_j_d_lcl_2d, Kokkos::ALL (), 0);
        Z4_d_raw = Z4_j_d_lcl.data();

        TEST_ASSERT( Z4_j_d_lcl.data () == Z4_dv_j.d_view.data () );
        if (Z4_j_d_lcl.data () != Z4_dv_j.d_view.data ()) {
          out << "Z4_j_d_lcl.data() = " << Z4_j_d_lcl.data ()
              << ", Z4_dv_j.d_view.data() = " << Z4_dv_j.d_view.data ()
              << endl;
        }
        TEST_ASSERT( Z4_j_d_lcl.extent (0) == Z4_dv_j.d_view.extent (0) );
      }

      if (j == 0) {
        TEST_ASSERT( Z4_h_raw == Z4_dv.h_view.data () );
        TEST_ASSERT( Z4_d_raw == Z4_dv.d_view.data () );
      }
    }

    Tpetra::deep_copy (Z4, Z0);
    const bool Z0_Z4_equal_before = multiVectorsEqual (out, Z0, Z4);
    TEST_ASSERT( Z0_Z4_equal_before );

    Z4.reduce ();
    const bool Z1_Z4_equal_after_reduce = multiVectorsEqual (out, Z1, Z4);
    TEST_ASSERT( Z1_Z4_equal_after_reduce );
  }

  {
    // Make sure Z5 has a stride greater than its number of rows.
    const size_t Z5_stride = Z0.getLocalLength () + lclNumRows;
    dual_view_type Z5_dv_extra ("Z5", Z5_stride, numCols);
    auto Z5_dv = Kokkos::subview (Z5_dv_extra,
                                  pair_type (0, lclNumRows),
                                  pair_type (0, numCols));
    TEST_ASSERT( LO (Z5_dv.extent (0) == lclNumRows ) );
    TEST_ASSERT( LO (Z5_dv.extent (1) == numCols ) );
    // Kokkos could in theory insert padding in the row dimension.
    TEST_ASSERT( size_t (Z5_dv.d_view.stride (1)) >= Z5_stride );
    TEST_ASSERT( size_t (Z5_dv.h_view.stride (1)) >= Z5_stride );

    MV Z5 (lclMap, Z5_dv);
    Tpetra::deep_copy (Z5, Z0);

    const bool Z0_Z5_equal_before = multiVectorsEqual (out, Z0, Z5);
    TEST_ASSERT( Z0_Z5_equal_before );

    reduceMultiVector (Z5);
    const bool Z1_Z5_equal_after_reduceMultiVector =
      multiVectorsEqual (out, Z1, Z5);
    TEST_ASSERT( Z1_Z5_equal_after_reduceMultiVector );
  }

  {
    // Make sure Z6 has a stride greater than its number of rows.
    const size_t Z6_stride = Z0.getLocalLength () + lclNumRows;
    dual_view_type Z6_dv_extra ("Z6", Z6_stride, numCols);
    auto Z6_dv = Kokkos::subview (Z6_dv_extra,
                                  pair_type (0, lclNumRows),
                                  pair_type (0, numCols));
    TEST_ASSERT( LO (Z6_dv.extent (0) == lclNumRows ) );
    TEST_ASSERT( LO (Z6_dv.extent (1) == numCols ) );
    // Kokkos could in theory insert padding in the row dimension.
    TEST_ASSERT( size_t (Z6_dv.d_view.stride (1)) >= Z6_stride );
    TEST_ASSERT( size_t (Z6_dv.h_view.stride (1)) >= Z6_stride );

    MV Z6 (lclMap, Z6_dv, Z6_dv_extra);
    Tpetra::deep_copy (Z6, Z0);

    const bool Z0_Z6_equal_before = multiVectorsEqual (out, Z0, Z6);
    TEST_ASSERT( Z0_Z6_equal_before );

    reduceMultiVector2 (Z6);
    const bool Z1_Z6_equal_after_reduceMultiVector2 =
      multiVectorsEqual (out, Z1, Z6);
    TEST_ASSERT( Z1_Z6_equal_after_reduceMultiVector2 );

    out << endl << "Z1:" << endl;
    Z1.describe (out, Teuchos::VERB_EXTREME);
    out << endl << "Z6:" << endl;
    Z6.describe (out, Teuchos::VERB_EXTREME);
  }
}

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, reduce_strided, SCALAR, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV( UNIT_TEST_GROUP )

} // namespace (anonymous)
