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

namespace { // (anonymous)
  using Tpetra::TestingUtilities::getDefaultComm;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Comm;
  using Teuchos::outArg;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::endl;
  typedef Tpetra::global_size_t GST;

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Map, LocallyFitted, LO, GO, NT )
  {
    typedef Tpetra::Map<LO, GO, NT> map_type;
    int lclSuccess = 1;
    int gblSuccess = 0; // output argument

    out << "Test Map::locallyFitted" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const Comm<int> > comm = getDefaultComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();
    const GST GSTI = Teuchos::OrdinalTraits<GST>::invalid ();
    const GO indexBase = 0;

    // Situations to test:
    //
    // - (map1, map2) is (contiguous, noncontiguous)
    // - map1 is (fitted, not fitted) to map2
    // - Different reasons for not being fitted:
    //     - Remote indices not all at the end
    //     - Non-remote indices not in the same order
    //     - Some indices in the domain Map, not the column Map (???)

    // Test the case where map1 is uniform contiguous
    {
      out << "map1 is uniform contiguous" << endl;
      Teuchos::OSTab tab1 (out);

      const LO lclNumInds1 = 3;
      const GO gblNumInds1 = numProcs * lclNumInds1;
      map_type map1 (static_cast<GST> (gblNumInds1), indexBase, comm);

      {
        out << "map2 is uniform contiguous, and the same as map1" << endl;
        Teuchos::OSTab tab2 (out);

        const LO lclNumInds2 = 3;
        const GO gblNumInds2 = numProcs * lclNumInds2;
        map_type map2 (static_cast<GST> (gblNumInds2), indexBase, comm);

        const bool compat = map1.isCompatible (map2);
        TEST_ASSERT( compat );
        const bool same = map1.isSameAs (map2);
        TEST_ASSERT( same );
        const bool fitted = map2.isLocallyFitted (map1);
        TEST_ASSERT( fitted );
      }

      {
        out << "map2 is uniform contiguous, and NOT the same as map1" << endl;
        Teuchos::OSTab tab2 (out);

        // Different total number of global indices, so can't be fitted.
        const LO lclNumInds2 = 4;
        const GO gblNumInds2 = numProcs * lclNumInds2;
        map_type map2 (static_cast<GST> (gblNumInds2), indexBase, comm);

        const bool compat = map1.isCompatible (map2);
        TEST_ASSERT( ! compat );
        const bool same = map1.isSameAs (map2);
        TEST_ASSERT( ! same );
        const bool fitted = map2.isLocallyFitted (map1);
        // Always fitted on Process 0, never on other processes.
        TEST_ASSERT( (fitted && myRank == 0) || (! fitted && myRank != 0) );
      }

      {
        out << "map2 is noncontiguous, yet the same as map1" << endl;
        Teuchos::OSTab tab2 (out);

        const LO lclNumInds2 = lclNumInds1;
        const GO gblNumInds2 = gblNumInds1;
        Teuchos::Array<GO> gblInds2 (lclNumInds2);
        for (LO lclInd = 0; lclInd < lclNumInds2; ++lclInd) {
          gblInds2[lclInd] = map1.getGlobalElement (lclInd);
        }
        map_type map2 (static_cast<GST> (gblNumInds2), gblInds2 (), indexBase, comm);

        const bool compat = map1.isCompatible (map2);
        TEST_ASSERT( compat );
        const bool same = map1.isSameAs (map2);
        TEST_ASSERT( same );
        const bool fitted = map2.isLocallyFitted (map1);
        TEST_ASSERT( fitted );
      }

      {
        out << "map2 is noncontiguous and compatible with map1, "
          "yet not the same" << endl;
        Teuchos::OSTab tab2 (out);

        const LO lclNumInds2 = lclNumInds1;
        const GO gblNumInds2 = gblNumInds1;
        Teuchos::Array<GO> gblInds2 (lclNumInds2);
        for (LO lclInd = 0; lclInd < lclNumInds2; ++lclInd) {
          // Different order.
          gblInds2[(lclNumInds2 - 1) - lclInd] = map1.getGlobalElement (lclInd);
        }
        map_type map2 (static_cast<GST> (gblNumInds2), gblInds2 (), indexBase, comm);

        const bool compat = map1.isCompatible (map2);
        TEST_ASSERT( compat );
        const bool same = map1.isSameAs (map2);
        TEST_ASSERT( ! same );
        const bool fitted = map2.isLocallyFitted (map1);
        TEST_ASSERT( ! fitted );
      }

      if (numProcs > 1) {
        out << "map2 is noncontiguous, overlapping, and fitted "
          "(numProcs > 1)" << endl;
        Teuchos::OSTab tab2 (out);

        const LO lclNumInds2 = lclNumInds1 + 1;
        Teuchos::Array<GO> gblInds2 (lclNumInds2);
        for (LO lclInd = 0; lclInd < lclNumInds2; ++lclInd) {
          if (lclInd < lclNumInds1) {
            gblInds2[lclInd] = map1.getGlobalElement (lclInd);
          }
          else {
            // This is only a remote index if numProcs > 1.
            // However, we exclude this test if numProcs == 1.
            const GO gblInd = (gblInds2[lclNumInds1 - 1] + 1) % gblNumInds1;
            gblInds2[lclInd] = gblInd;
          }
        }
        map_type map2 (GSTI, gblInds2 (), indexBase, comm);

        const bool compat = map1.isCompatible (map2);
        TEST_ASSERT( ! compat );
        const bool same = map1.isSameAs (map2);
        TEST_ASSERT( ! same );
        const bool fitted = map2.isLocallyFitted (map1);
        TEST_ASSERT( fitted );
      }

      if (numProcs > 1) {
        out << "map2 is noncontiguous, overlapping, and NOT fitted "
          "(numProcs > 1), because remotes not all at the end" << endl;
        Teuchos::OSTab tab2 (out);

        const LO lclNumInds2 = lclNumInds1 + 1;
        Teuchos::Array<GO> gblInds2 (lclNumInds2);

        // This is only a remote index if numProcs > 1.
        // However, we exclude this test if numProcs == 1.
        const GO gblRemoteInd = (map1.getGlobalElement (lclNumInds1 - 1) + 1) % gblNumInds1;
        gblInds2[0] = gblRemoteInd;
        for (LO lclInd = 0; lclInd < lclNumInds1; ++lclInd) {
          gblInds2[lclInd+1] = map1.getGlobalElement (lclInd);
        }
        map_type map2 (GSTI, gblInds2 (), indexBase, comm);

        const bool compat = map1.isCompatible (map2);
        TEST_ASSERT( ! compat );
        const bool same = map1.isSameAs (map2);
        TEST_ASSERT( ! same );
        const bool fitted = map2.isLocallyFitted (map1);
        TEST_ASSERT( ! fitted );
      }
    } // map1 is uniform contiguous

    // Test the case where map1 is noncontiguous
    {
      out << "map1 is noncontiguous" << endl;
      Teuchos::OSTab tab1 (out);

      const LO lclNumInds1 = 3;
      const GO gblNumInds1 = numProcs * lclNumInds1;
      Teuchos::Array<GO> gblInds1 (lclNumInds1);
      {
        const GO gblIndStart1 = static_cast<GO> (myRank) * static_cast<GO> (lclNumInds1);
        // Same global indices on each process as in a uniform
        // contiguous Map, but just in a different order (backwards,
        // in this case).
        for (LO lclInd = 0; lclInd < lclNumInds1; ++lclInd) {
          gblInds1[(lclNumInds1 - 1) - lclInd] = gblIndStart1 + static_cast<GO> (lclInd);
        }
      }
      map_type map1 (GSTI, gblInds1 (), indexBase, comm);

      {
        out << "map2 is noncontiguous and same as map1" << endl;
        Teuchos::OSTab tab2 (out);

        map_type map2 (GSTI, gblInds1 (), indexBase, comm);

        const bool compat = map1.isCompatible (map2);
        TEST_ASSERT( compat );
        const bool same = map1.isSameAs (map2);
        TEST_ASSERT( same );
        const bool fitted = map2.isLocallyFitted (map1);
        TEST_ASSERT( fitted );
      }

      {
        out << "map2 is noncontiguous and compatible with map1, "
          "yet not the same" << endl;
        Teuchos::OSTab tab2 (out);

        Teuchos::Array<GO> gblInds2 (lclNumInds1);
        for (LO lclInd = 0; lclInd < lclNumInds1; ++lclInd) {
          // Change the order
          gblInds2[lclInd] = map1.getGlobalElement ((lclNumInds1 - 1) - lclInd);
        }
        map_type map2 (GSTI, gblInds2 (), indexBase, comm);

        const bool compat = map1.isCompatible (map2);
        TEST_ASSERT( compat );
        const bool same = map1.isSameAs (map2);
        TEST_ASSERT( ! same );
        const bool fitted = map2.isLocallyFitted (map1);
        TEST_ASSERT( ! fitted );
      }

      if (numProcs > 1) {
        out << "map2 is noncontiguous, overlapping, and fitted "
          "(numProcs > 1)" << endl;
        Teuchos::OSTab tab2 (out);

        const LO lclNumInds2 = lclNumInds1 + 1;
        Teuchos::Array<GO> gblInds2 (lclNumInds2);
        for (LO lclInd = 0; lclInd < lclNumInds2; ++lclInd) {
          if (lclInd < lclNumInds1) {
            gblInds2[lclInd] = map1.getGlobalElement (lclInd);
          }
          else {
            // This is only a remote index if numProcs > 1.
            // However, we exclude this test if numProcs == 1.
            const GO gblInd = (gblInds2[lclNumInds1 - 1] + 1) % gblNumInds1;
            gblInds2[lclInd] = gblInd;
          }
        }
        map_type map2 (GSTI, gblInds2 (), indexBase, comm);

        const bool compat = map1.isCompatible (map2);
        TEST_ASSERT( ! compat );
        const bool same = map1.isSameAs (map2);
        TEST_ASSERT( ! same );
        const bool fitted = map2.isLocallyFitted (map1);
        TEST_ASSERT( fitted );
      }

      if (numProcs > 1) {
        out << "map2 is noncontiguous, overlapping, and NOT fitted "
          "(numProcs > 1), because remotes not all at the end" << endl;
        Teuchos::OSTab tab2 (out);

        const LO lclNumInds2 = lclNumInds1 + 1;
        Teuchos::Array<GO> gblInds2 (lclNumInds2);

        // This is only a remote index if numProcs > 1.
        // However, we exclude this test if numProcs == 1.
        const GO gblRemoteInd = (map1.getGlobalElement (lclNumInds1 - 1) + 1) % gblNumInds1;
        gblInds2[0] = gblRemoteInd;
        for (LO lclInd = 0; lclInd < lclNumInds1; ++lclInd) {
          gblInds2[lclInd+1] = map1.getGlobalElement (lclInd);
        }
        map_type map2 (GSTI, gblInds2 (), indexBase, comm);

        const bool compat = map1.isCompatible (map2);
        TEST_ASSERT( ! compat );
        const bool same = map1.isSameAs (map2);
        TEST_ASSERT( ! same );
        const bool fitted = map2.isLocallyFitted (map1);
        TEST_ASSERT( ! fitted );
      }
    } // map1 is noncontiguous

    // Make sure that the test passed on all processes.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0; // output argument
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
  }

  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP( LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Map, LocallyFitted, LO, GO, NT )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LGN(UNIT_TEST_GROUP)

} // namespace (anonymous)


