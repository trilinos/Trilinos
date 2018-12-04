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
    typedef Kokkos::Details::ArithTraits<SC> KAT;
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
      X.sync_host ();
      auto X_lcl = X.getLocalViewHost ();

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

