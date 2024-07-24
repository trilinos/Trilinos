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
#include "Teuchos_ScalarTraits.hpp"
#include "Kokkos_ArithTraits.hpp"

namespace { // (anonymous)

  using Tpetra::TestingUtilities::getDefaultComm;
  using Teuchos::RCP;
  using std::endl;

  //
  // UNIT TESTS
  //

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, TwoArgRandomize, SC, LO, GO, NT )
  {
    typedef Tpetra::Map<LO, GO, NT> map_type;
    typedef Tpetra::MultiVector<SC, LO, GO, NT> MV;
    typedef Teuchos::ScalarTraits<SC> STS;
    typedef Kokkos::ArithTraits<SC> KAT;
    typedef typename MV::impl_scalar_type IST;

    out << "Test that two-argument Tpetra::MultiVector::randomize builds, "
      "and that it returns values not outside the given range" << endl;

    auto comm = getDefaultComm ();
    const int numProcs = comm->getSize ();
    // We need enough rows to make sure that we exercise the
    // pseudorandom number generator sufficiently.
    const LO lclNumRows = 1000;
    const GO gblNumRows =
      static_cast<GO> (lclNumRows) * static_cast<GO> (numProcs);
    const GO indexBase = 0;
    RCP<const map_type> map (new map_type (gblNumRows, lclNumRows,
                                           indexBase, comm));

    // Constrain X.randomize to produce numbers in the range (3,5).
    // Don't test the complex case, since it's less well defined what
    // "in a range" means in that case (see #567 for the motivation
    // for two-argument randomize; Chebyshev is more sensible for the
    // real case anyway).
    const SC ONE = STS::one ();
    const SC THREE = ONE + ONE + ONE;
    const SC FIVE = ONE + ONE + ONE + ONE + ONE;

    const LO numVecs = 3;
    MV X (map, numVecs);
    const SC minVal = THREE;
    const SC maxVal = FIVE;
    TEST_NOTHROW( X.randomize (minVal, maxVal) );

    bool allInRange = true;
    if (! STS::isComplex) {
      auto X_lcl = X.getLocalViewHost(Tpetra::Access::ReadOnly);

      for (LO j = 0; j < numVecs; ++j) {
        auto X_lcl_j = Kokkos::subview (X_lcl, Kokkos::ALL (), j);

        for (LO i = 0; i < lclNumRows; ++i) {
          const auto X_lcl_ij = X_lcl_j(i);

          if (KAT::real (X_lcl_ij) < KAT::real (static_cast<IST> (minVal)) ||
              KAT::real (X_lcl_ij) > KAT::real (static_cast<IST> (maxVal))) {
            allInRange = false;
            break;
          }
        } // for each local row

        if (! allInRange) {
          break;
        }
      } // for each column
    }

    TEST_ASSERT( allInRange );
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SC, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, TwoArgRandomize, SC, LO, GO, NT )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV( UNIT_TEST_GROUP )

} // namespace (anonymous)

