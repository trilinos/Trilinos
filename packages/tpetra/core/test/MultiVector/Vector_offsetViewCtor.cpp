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
    using vector_type = VectorType;
    using LO = typename vector_type::local_ordinal_type;

    x.modify_device ();
    auto x_lcl_d_2d = x.getLocalViewDevice ();
    auto x_lcl_d = Kokkos::subview (x_lcl_d_2d, Kokkos::ALL (), 0);

    using execution_space = typename vector_type::execution_space;
    using range_type = Kokkos::RangePolicy<execution_space, LO>;
    Kokkos::parallel_for
      ("Initial Vector fill", range_type (0, lclNumRows),
       KOKKOS_LAMBDA (const LO lclRow) {
        x_lcl_d(lclRow) = toScalar<IST> (lclRow+1);
      });
    execution_space().fence ();
    x.sync_host ();
  }

  template<class VectorType>
  void
  checkVectorEntries (bool& success,
                      Teuchos::FancyOStream& out,
                      VectorType& x_offset,
                      const typename VectorType::local_ordinal_type rowOffset)
  {
    using vector_type = VectorType;
    using LO = typename vector_type::local_ordinal_type;

    TEST_ASSERT( ! x_offset.need_sync_host () );

    auto x_lcl_h_2d = x_offset.getLocalViewHost ();
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
      myGblRowInds = originalMap->getGlobalElement (lclRow);
    }

    const GST INV = Teuchos::OrdinalTraits<GST>::invalid ();
    for (LO rowOffset : { LO (0), LO (1), LO (3) }) {
      const LO newLclNumRows = lclNumRows - rowOffset;
      const GO newIndexBase = indexBase + rowOffset;
      RCP<const map_type> map_offset =
        rcp (new map_type (INV, myGblRowInds.data () + rowOffset,
                           newLclNumRows, newIndexBase, comm));
      vector_type x_offset (x, map_offset);

      const bool expectedMap = map_offset->isSameAs (* (x_offset->getMap ()));
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
        x.sync_host ();
        auto x_offset_lcl_h_2d = x_offset.getLocalViewDevice ();
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
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( Vector, OffsetViewCtor, ST, LO, GO, NT )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV( UNIT_TEST_GROUP )

} // namespace (anonymous)
