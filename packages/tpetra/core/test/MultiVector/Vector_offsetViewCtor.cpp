// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Map.hpp"

namespace { // (anonymous)

  template<class ScalarType, class IntegerType>
  KOKKOS_INLINE_FUNCTION ScalarType toScalar (const IntegerType x)
  {
    using KAT = Kokkos::ArithTraits<ScalarType>;
    using mag_type = typename KAT::mag_type;
    // The double cast handles conversions like integer to
    // std::complex<double>, where no one-step conversion exists.
    return static_cast<ScalarType> (static_cast<mag_type> (x));
  }

  template<class VectorType>
  void restoreVectorEntries (VectorType& x)
  {
    using LO = typename VectorType::local_ordinal_type;

    auto x_lcl_d_2d = x.getLocalViewDevice(Tpetra::Access::ReadWrite);
    auto x_lcl_d = Kokkos::subview (x_lcl_d_2d, Kokkos::ALL (), 0);

    using execution_space = typename VectorType::execution_space;
    using range_type = Kokkos::RangePolicy<execution_space, LO>;
    using IST = typename VectorType::impl_scalar_type;
    auto lclNumRows = x.getLocalLength();
    Kokkos::parallel_for
      ("Initial Vector fill", range_type (0, lclNumRows),
       KOKKOS_LAMBDA (const LO lclRow) {
        x_lcl_d(lclRow) = toScalar<IST> (lclRow);
      });
    execution_space().fence ();
  }

  template<class VectorType>
  void
  checkVectorEntries (bool& success,
                      Teuchos::FancyOStream& out,
                      VectorType& x_offset,
                      const typename VectorType::local_ordinal_type rowOffset)
  {
    using LO = typename VectorType::local_ordinal_type;
    using IST = typename VectorType::impl_scalar_type;

    auto x_lcl_h_2d = x_offset.getLocalViewHost(Tpetra::Access::ReadOnly);
    auto x_lcl_h = Kokkos::subview (x_lcl_h_2d, Kokkos::ALL (), 0);

    const LO newLclNumRows = static_cast<LO> (x_offset.getLocalLength ());
    for (LO newLclRow = 0; newLclRow < newLclNumRows; ++newLclRow) {
      TEST_EQUALITY( x_lcl_h(newLclRow),
                     toScalar<IST> (newLclRow + rowOffset) );
    }
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( Vector, OffsetViewCtor, ST, LO, GO, NT )
  {
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::endl;
    using GST = Tpetra::global_size_t;
    using map_type = Tpetra::Map<LO, GO, NT>;
    using vector_type = Tpetra::Vector<ST, LO, GO, NT>;
    using IST = typename vector_type::impl_scalar_type;
    int lclSuccess = 1;
    int gblSuccess = 1;

    out << "Test Vector's \"offset view\" constructor" << endl;
    Teuchos::OSTab tab0 (out);

    const auto comm = Tpetra::getDefaultComm ();
    const LO lclNumRows = 13;
    const GO gblNumRows = static_cast<GO> (comm->getSize ()) *
      static_cast<GO> (lclNumRows);
    const GO indexBase = 0;

    RCP<const map_type> originalMap =
      rcp (new map_type (gblNumRows, lclNumRows, indexBase, comm));

    vector_type x (originalMap);
    restoreVectorEntries (x);

    std::vector<GO> myGblRowInds (lclNumRows);
    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      myGblRowInds[lclRow] = originalMap->getGlobalElement (lclRow);
    }

    const GST INV = Teuchos::OrdinalTraits<GST>::invalid ();
    for (LO rowOffset : { LO (0), LO (1), LO (3) }) {
      const LO newLclNumRows = lclNumRows - rowOffset;
      const GO newIndexBase = indexBase + rowOffset;
      RCP<const map_type> map_offset =
        rcp (new map_type (INV, myGblRowInds.data () + rowOffset,
                           newLclNumRows, newIndexBase, comm));
      vector_type x_offset (x, map_offset, rowOffset);

      const bool expectedMap = map_offset->isSameAs (* (x_offset.getMap ()));
      TEST_ASSERT( expectedMap );
      TEST_EQUALITY( static_cast<LO> (x_offset.getLocalLength ()),
                     newLclNumRows );
      if (success) {
        checkVectorEntries (success, out, x_offset, rowOffset);
      }

      lclSuccess = success ? 1 : 0;
      gblSuccess = 0; // output argument
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY_CONST( gblSuccess, 1 );
      if (! gblSuccess) {
        return; // no point in continuing
      }

      // If I modify the original Vector, does the new Vector see it?
      // That is, is the new Vector a view of the original Vector?
      {
        x.putScalar (Teuchos::ScalarTraits<ST>::one ());
        auto x_offset_lcl_h_2d = x_offset.getLocalViewHost(Tpetra::Access::ReadWrite);
        auto x_offset_lcl_h =
          Kokkos::subview (x_offset_lcl_h_2d, Kokkos::ALL (), 0);
        for (LO newLclRow = 0; newLclRow < newLclNumRows; ++newLclRow) {
          TEST_EQUALITY( x_offset_lcl_h(newLclRow),
                         Kokkos::ArithTraits<IST>::one () );
        }
      }
      restoreVectorEntries (x);
    }

    lclSuccess = success ? 1 : 0;
    gblSuccess = 0; // output argument
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (! gblSuccess) {
      return; // no point in continuing
    }
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Vector, OffsetViewCtor, SCALAR, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV( UNIT_TEST_GROUP )

} // namespace (anonymous)
