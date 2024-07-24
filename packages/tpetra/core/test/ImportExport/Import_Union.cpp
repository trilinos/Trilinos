// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"

#include <TpetraCore_ETIHelperMacros.h>
#include <Tpetra_Core.hpp>
#include <Tpetra_Distributor.hpp>
#include <Tpetra_Import.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Details_gathervPrint.hpp>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_OrdinalTraits.hpp>

#include <algorithm>


namespace {

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions (true);
  }

#define TPETRA_IMPORT_UNION_RUN_AND_CATCH_EXCEPTION( CODE, NAME ) \
  do { \
    std::ostringstream os; \
    try { \
      CODE; \
    } catch (std::exception& e) { \
      lclSuccess = 0; \
      os << "Proc " << comm->getRank () << ": " << NAME << " threw exception: " << e.what () << endl; \
    } \
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess)); \
    if (gblSuccess != 1) { \
      out << NAME << " FAILED!" << endl; \
      Tpetra::Details::gathervPrint (out, os.str (), *comm); \
      success = false; \
      return; \
    } \
  } while (false)

#define TPETRA_IMPORT_UNION_TEST_EQUALITY( THING1, THING2, FAIL_MSG, STOP_ON_FAIL )  \
  do { \
    const int R = comm->getRank (); \
    std::ostringstream os; \
    if (THING1 != THING2) { \
      lclSuccess = 0; \
      os << "Proc " << comm->getRank () << ": " << FAIL_MSG << endl; \
    } \
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess)); \
    if (gblSuccess != 1) { \
      out << "Equality test FAILED on one or more processes!" << endl; \
      Tpetra::Details::gathervPrint (out, os.str (), *comm); \
      success = false; \
      if (STOP_ON_FAIL) { \
        return; \
      } \
    } \
  } while (false)

  //
  // UNIT TESTS
  //

  // Test whether the two Import objects represent the same
  // communication pattern.
  //
  // We consider "same" loosely; for example, we allow the two Import
  // objects to have different numbers of "same" IDs, as long as they
  // correctly (that is, in the correct order) account for those IDs
  // in the "permute" IDs.
  template<class LO, class GO, class NT>
  void
  compareImports (bool& success, Teuchos::FancyOStream& out,
                  const Tpetra::Import<LO, GO, NT>& imp1,
                  const Tpetra::Import<LO, GO, NT>& imp2,
                  const Teuchos::Comm<int>& comm)
  {
    using Teuchos::ArrayView;
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;

    int lclSuccess = 1;
    int gblSuccess = 1;

    // Start by comparing the "same" and "permute" IDs.  We have to
    // consider them together, because they don't have to match
    // exactly in order to match functionally.  That is, the number of
    // "sames" may differ, but as long as the "missing" "sames" show
    // up in the "permutes" in the right order, then the two Import
    // objects behave in the same way.

    ArrayView<const LO> permFromLIDs1 = imp1.getPermuteFromLIDs ();
    ArrayView<const LO> permFromLIDs2 = imp2.getPermuteFromLIDs ();
    ArrayView<const LO> permToLIDs1 = imp1.getPermuteToLIDs ();
    ArrayView<const LO> permToLIDs2 = imp2.getPermuteToLIDs ();
    const size_t numPermIDs1 = imp1.getNumPermuteIDs ();
    const size_t numPermIDs2 = imp2.getNumPermuteIDs ();

    // Sanity check: for each Import, reported number of permutes
    // (method that returns size_t) should match the actual number of
    // permutes (method that returns ArrayView).
    TEST_EQUALITY( static_cast<size_t> (permFromLIDs1.size ()), numPermIDs1 );
    TEST_EQUALITY( static_cast<size_t> (permFromLIDs2.size ()), numPermIDs2 );
    TEST_EQUALITY( static_cast<size_t> (permToLIDs1.size ()), numPermIDs1 );
    TEST_EQUALITY( static_cast<size_t> (permToLIDs2.size ()), numPermIDs2 );

    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    success = success && gblSuccess == 1;
    if (! success) {
      // One of the two Imports is in an invalid state, so don't
      // bother beyond this point.
      return;
    }

    // Note that the counts of sames and permutes may differ across
    // processes.  This means that we have to be careful using
    // all-reduces to check consistency.

    const size_t numSame1 = imp1.getNumSameIDs ();
    const size_t numSame2 = imp2.getNumSameIDs ();

    if (numSame1 == numSame2) {
      // Number of "sames" match, so we can compare permutes directly.
      TEST_EQUALITY( numPermIDs1, numPermIDs2 );
      // Compare permute "from" LIDs.
      const bool permFromEq =
        std::equal (permFromLIDs1.begin (), permFromLIDs1.end (),
                    permFromLIDs2.begin ());
      TEST_EQUALITY( permFromEq, true );
      // Compare permute "to" LIDs.
      const bool permToEq =
        std::equal (permToLIDs1.begin (), permToLIDs1.end (),
                    permToLIDs2.begin ());
      TEST_EQUALITY( permToEq, true );
    }
    else { // numSame1 != numSame2
      // Some of the sameIDs could come from the permuteIDs.  Check
      // that the total of sames and permutes is correct, and that the
      // "missing" sames in the permutes are correct (order matters).
      const bool sameTotals = numSame1 + numPermIDs1 == numSame2 + numPermIDs2;
      TEST_EQUALITY_CONST( sameTotals, true );

      if (sameTotals) {
        const size_t numSameMin = std::min (numSame1, numSame2);
        const size_t numSameDiff = (numSame1 > numSame2) ?
          (numSame1 - numSame2) : (numSame2 - numSame1);
        ArrayView<const LO> permFromDiff, permToDiff;
        if (numSame1 > numSame2) {
          // imp2 must have imp1's "sames" in its "permutes"
          permFromDiff = permFromLIDs2 (numSame2, numSameDiff);
          permToDiff = permToLIDs2 (numSame2, numSameDiff);
        }
        else { // imp1 must have imp2's "sames" in its "permutes"
          permFromDiff = permFromLIDs1 (numSame1, numSameDiff);
          permToDiff = permToLIDs1 (numSame1, numSameDiff);
        }

        for (size_t k = 0; k < numSameDiff; ++k) {
          TEST_EQUALITY( static_cast<size_t> (permFromDiff[k]), numSameMin + k );
          // If the "permutes" in question belong with the "sames,"
          // then they shouldn't actually permute.
          TEST_EQUALITY( permFromDiff[k], permToDiff[k] );
        }

        // Test that the remaining permutes line up.
        ArrayView<const LO> permFromA, permFromB, permToA, permToB;
        if (numSame1 > numSame2) { // exclude "sames" from imp2's "permutes"
          permFromA = permFromLIDs2 (numSameDiff, numPermIDs2 - numSameDiff);
          permFromB = permFromLIDs1;
          permToA = permToLIDs2 (numSameDiff, numPermIDs2 - numSameDiff);
          permToB = permToLIDs1;
        }
        else { // exclude "sames" from imp1's "permutes"
          permFromA = permFromLIDs1 (numSameDiff, numPermIDs1 - numSameDiff);
          permFromB = permFromLIDs2;
          permToA = permToLIDs1 (numSameDiff, numPermIDs1 - numSameDiff);
          permToB = permToLIDs2;
        }
        TEST_EQUALITY( permFromA.size (), permFromB.size () ); // sanity check
        TEST_EQUALITY( permToA.size (), permToB.size () ); // sanity check

        for (size_t k = 0; k < static_cast<size_t> (permFromA.size ()); ++k) {
          TEST_EQUALITY( permFromA[k], permFromB[k] );
          TEST_EQUALITY( permToA[k], permToB[k] );
        }
      }
    }

    TEST_EQUALITY( imp1.getNumRemoteIDs (), imp2.getNumRemoteIDs () );
    TEST_EQUALITY( imp1.getNumExportIDs (), imp2.getNumExportIDs () );

    Tpetra::Distributor& dist1 = imp1.getDistributor ();
    Tpetra::Distributor& dist2 = imp2.getDistributor ();

    TEST_EQUALITY( dist1.getNumReceives (), dist2.getNumReceives () );
    TEST_EQUALITY( dist1.getNumSends (), dist2.getNumSends () );
    TEST_EQUALITY( dist1.hasSelfMessage (), dist2.hasSelfMessage () );
    TEST_EQUALITY( dist1.getMaxSendLength (), dist2.getMaxSendLength () );
    TEST_EQUALITY( dist1.getTotalReceiveLength (), dist2.getTotalReceiveLength () );

    ArrayView<const int> imagesFrom1 = dist1.getProcsFrom ();
    ArrayView<const int> imagesFrom2 = dist2.getProcsFrom ();
    const bool sameProcsFrom =
      imagesFrom1.size () == imagesFrom2.size () &&
      std::equal (imagesFrom1.begin (),
                  imagesFrom1.end (),
                  imagesFrom2.begin ());
    TEST_EQUALITY( sameProcsFrom, true );

    ArrayView<const int> imagesTo1 = dist1.getProcsTo ();
    ArrayView<const int> imagesTo2 = dist2.getProcsTo ();
    const bool sameProcsTo =
      imagesTo1.size () == imagesTo2.size () &&
      std::equal (imagesTo1.begin (), imagesTo1.end (),
                  imagesTo2.begin ());
    TEST_EQUALITY( sameProcsTo, true );

    ArrayView<const size_t> lengthsFrom1 = dist1.getLengthsFrom ();
    ArrayView<const size_t> lengthsFrom2 = dist2.getLengthsFrom ();
    const bool sameLengthsFrom =
      lengthsFrom1.size () == lengthsFrom2.size () &&
      std::equal (lengthsFrom1.begin (), lengthsFrom1.end (),
                  lengthsFrom2.begin ());
    TEST_EQUALITY( sameLengthsFrom, true );

    ArrayView<const size_t> lengthsTo1 = dist1.getLengthsTo ();
    ArrayView<const size_t> lengthsTo2 = dist2.getLengthsTo ();
    const bool sameLengthsTo =
      lengthsTo1.size () == lengthsTo2.size () &&
      std::equal (lengthsTo1.begin (), lengthsTo1.end (),
                  lengthsTo2.begin ());
    TEST_EQUALITY( sameLengthsTo, true );
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( ImportUnion, ContigPlusContig, LocalOrdinalType, GlobalOrdinalType, NodeType )
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Teuchos::Comm;
    using Teuchos::getFancyOStream;
    using Teuchos::FancyOStream;
    using Teuchos::OSTab;
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcpFromRef;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using std::cerr;
    using std::cout;
    using std::endl;
    typedef Tpetra::global_size_t GST;
    typedef LocalOrdinalType LO;
    typedef GlobalOrdinalType GO;
    typedef NodeType NT;
    typedef Tpetra::Map<LO, GO, NT> map_type;
    typedef Tpetra::Import<LO, GO, NT> import_type;
    typedef Tpetra::Vector<>::scalar_type ST;
    typedef Tpetra::Vector<ST, LO, GO, NT> vector_type;

    int lclSuccess = 1; // local error flag
    int gblSuccess = 1; // global error flag (result of all-reduce on lclSuccess)

    out << "Tpetra::Import::setUnion test" << endl;
    OSTab tab1 (out);
    out << "Both target Maps contiguous" << endl;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const GO indexBase = 0;
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    const GO n = 10;
    // Let r = comm->getRank() and P = comm->getSize().  Then:
    //
    // Source Map (for both): indexBase + {n*r, ..., n*(r + 1) - 1}.
    // Target Map 1: indexBase + {max(n*r - 1, 0), ..., min(n(r + 1), n*P - 1)}.
    // Target Map 2: indexBase + {max(n*r - 2, 0), ..., min(n(r + 1) + 1, n*P - 1)}.
    //
    // Expected union target Map happens to be target Map 2 in this
    // case, except that the "remote" GIDs go at the end of the GID
    // list on each process.

    Array<GO> srcMapGids;
    for (GO k = n * myRank; k < n*(myRank + 1); ++k) {
      srcMapGids.push_back (indexBase + k);
    }
    Array<GO> tgtMap1Gids;
    {
      // std::max(n*myRank - 1, 0) doesn't work if GO is unsigned.
      const GO lower = (n*myRank < 1) ?
        static_cast<GO> (0) :
        static_cast<GO> (n*myRank - 1);
      const GO upper = std::min (n*(myRank + 1) + 1, n * numProcs);
      for (GO k = lower; k < upper; ++k) {
        tgtMap1Gids.push_back (indexBase + k);
      }
    }
    Array<GO> tgtMap2Gids;
    {
      // std::max(n*myRank - 2, 0) doesn't work if GO is unsigned.
      const GO lower = (n*myRank < 2) ?
        static_cast<GO> (0) :
        static_cast<GO> (n*myRank - 2);
      const GO upper = std::min (n*(myRank + 1) + 2, n * numProcs);
      for (GO k = lower; k < upper; ++k) {
        tgtMap2Gids.push_back (indexBase + k);
      }
    }

    Array<GO> unionTgtMapGids;
    // Non-remote GIDs first.
    std::copy (srcMapGids.begin (), srcMapGids.end (),
               std::back_inserter (unionTgtMapGids));
    // Remote GIDs last.
    //
    // Don't test (n * myRank - 2 >= 0), because GO might be unsigned.
    if (n * myRank >= 2) {
      unionTgtMapGids.push_back (indexBase + n * myRank - 2);
    }
    // Don't test (n * myRank - 1 >= 0), because GO might be unsigned.
    if (n * myRank >= 1) {
      unionTgtMapGids.push_back (indexBase + n * myRank - 1);
    }
    if (n * myRank + n < n*numProcs) {
      unionTgtMapGids.push_back (indexBase + n * myRank + n);
    }
    if (n * myRank + n + 1 < n*numProcs) {
      unionTgtMapGids.push_back (indexBase + n * myRank + n + 1);
    }

    out << "Making the Maps" << endl;

    RCP<const map_type> srcMap (new map_type (INVALID, srcMapGids (), indexBase, comm));
    RCP<const map_type> tgtMap1 (new map_type (INVALID, tgtMap1Gids (), indexBase, comm));
    RCP<const map_type> tgtMap2 (new map_type (INVALID, tgtMap2Gids (), indexBase, comm));
    RCP<const map_type> expectedUnionMap =
      rcp (new map_type (INVALID, unionTgtMapGids (), indexBase, comm));

    out << "Making the Import objects" << endl;

    // The two Import objects for which to compute the union.
    RCP<const import_type> imp1 (new import_type (srcMap, tgtMap1));
    RCP<const import_type> imp2 (new import_type (srcMap, tgtMap2));

    // The "expected" union of the two Import objects is the Import
    // from the source Map to the union of their target Maps.
    RCP<const import_type> expectedUnionImp =
      rcp (new import_type (srcMap, expectedUnionMap));

    // Compute setUnion using two different orders, to make sure that
    // the results do not depend on order (set union should commute).

    out << "Computing setUnion using first, second" << endl;
    RCP<const import_type> unionImp1;
    TPETRA_IMPORT_UNION_RUN_AND_CATCH_EXCEPTION
      ( unionImp1 = imp1->setUnion (*imp2), "First setUnion call" );

    if (unionImp1.is_null ()) {
      lclSuccess = 0;
    }
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      out << "First setUnion call returned null on some process!" << endl;
      return; // no sense in continuing
    }

    out << "Computing setUnion using second, first" << endl;
    RCP<const import_type> unionImp2;
    TPETRA_IMPORT_UNION_RUN_AND_CATCH_EXCEPTION
      ( unionImp2 = imp2->setUnion (*imp1), "Second setUnion call" );

    if (unionImp2.is_null ()) {
      lclSuccess = 0;
    }
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    if (gblSuccess != 1) {
      out << "Second setUnion call returned null on some process!" << endl;
      return; // no sense in continuing
    }

    out << "Running tests" << endl;

    // RCP<FancyOStream> cerrWrapped = getFancyOStream (rcpFromRef (cerr));

    // cerr << "Target Map of setUnion (first, second) result:" << endl;
    // unionImp1->getTargetMap ()->describe (*cerrWrapped, Teuchos::VERB_EXTREME);

    // cerr << "Expected Target Map:" << endl;
    // expectedUnionImp->getTargetMap ()->describe (*cerrWrapped, Teuchos::VERB_EXTREME);

    out << "Testing whether target Map (1,2) is same as expected target Map"
        << endl;
    const bool tgtMapSame1 =
      expectedUnionMap->isSameAs (* (unionImp1->getTargetMap ()));
    TEST_EQUALITY( tgtMapSame1, true );

    out << "Testing whether target Map (2,1) is same as expected target Map"
        << endl;
    const bool tgtMapSame2 =
      expectedUnionMap->isSameAs (* (unionImp2->getTargetMap ()));
    TEST_EQUALITY( tgtMapSame2, true );

    // If either of the target Maps is wrong, it doesn't make sense to
    // continue.  Furthermore, in that case, it's possible that
    // actually using the union Imports in a doImport or doExport
    // operation will break.
    if (! tgtMapSame1 || ! tgtMapSame2) {
      return;
    }

    // Use the three Imports in a doImport operation on Vectors.
    // They should all produce the same results.

    vector_type x (srcMap);
    vector_type y_expected (expectedUnionMap);
    vector_type y_actual_12 (unionImp1->getTargetMap ());
    vector_type y_actual_21 (unionImp2->getTargetMap ());

    x.randomize ();
    y_expected.putScalar (0.0);
    y_actual_12.putScalar (0.0);
    y_actual_21.putScalar (0.0);

    y_expected.doImport (x, *expectedUnionImp, Tpetra::ADD);
    y_actual_12.doImport (x, *unionImp1, Tpetra::ADD);
    y_actual_21.doImport (x, *unionImp2, Tpetra::ADD);

    out << "Testing whether union Import (1,2) works like "
      "expected Union Import with a Vector" << endl;
    {
      vector_type z_12 (y_expected.getMap ());
      Tpetra::deep_copy (z_12, y_expected);
      z_12.update (1.0, y_actual_12, -1.0);
      const typename vector_type::mag_type z_12_norm = z_12.norm2 ();
      out << "||y_expected - y_actual_12||_2 = " << z_12_norm << endl;
      TEST_EQUALITY( z_12_norm, 0.0 );
    }
    out << "Testing whether union Import (2,1) works like "
      "expected Union Import with a Vector" << endl;
    {
      vector_type z_21 (y_expected.getMap ());
      Tpetra::deep_copy (z_21, y_expected);
      z_21.update (1.0, y_actual_21, -1.0);
      const typename vector_type::mag_type z_21_norm = z_21.norm2 ();
      out << "||y_expected - y_actual_21||_2 = " << z_21_norm << endl;
      TEST_EQUALITY( z_21_norm, 0.0 );
    }

    out << "Test whether the Imports actually represent the same "
      "communication pattern" << endl;
    // Print counts of the different kinds of IDs from each of the
    // three Import objects.
    {
      std::ostringstream os;
      os << "Proc " << myRank << ":" << endl
         << "  unionImp1: {numSameIDs: " << unionImp1->getNumSameIDs ()
         << ", numPermuteIDs: " << unionImp1->getNumPermuteIDs ()
         << ", numRemoteIDs: " << unionImp1->getNumRemoteIDs ()
         << ", numExportIDs: " << unionImp1->getNumExportIDs ()
         << ", permuteFromLIDs: " << unionImp1->getPermuteFromLIDs ()
         << ", permuteToLIDs: " << unionImp1->getPermuteToLIDs ()
         << "}" << endl
         << "  unionImp2: {numSameIDs: " << unionImp2->getNumSameIDs ()
         << ", numPermuteIDs: " << unionImp2->getNumPermuteIDs ()
         << ", numRemoteIDs: " << unionImp2->getNumRemoteIDs ()
         << ", numExportIDs: " << unionImp2->getNumExportIDs ()
         << ", permuteFromLIDs: " << unionImp2->getPermuteFromLIDs ()
         << ", permuteToLIDs: " << unionImp2->getPermuteToLIDs ()
         << "}" << endl
         << "  expected:  {numSameIDs: " << expectedUnionImp->getNumSameIDs ()
         << ", numPermuteIDs: " << expectedUnionImp->getNumPermuteIDs ()
         << ", numRemoteIDs: " << expectedUnionImp->getNumRemoteIDs ()
         << ", numExportIDs: " << expectedUnionImp->getNumExportIDs ()
         << ", permuteFromLIDs: " << expectedUnionImp->getPermuteFromLIDs ()
         << ", permuteToLIDs: " << expectedUnionImp->getPermuteToLIDs ()
         << "}" << endl;
      Tpetra::Details::gathervPrint (out, os.str (), *comm);
    }

    compareImports (success, out, *unionImp1, *expectedUnionImp, *comm);
    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "*** imp1->setUnion(imp2) != expected Import ***" << endl;
    }

    compareImports (success, out, *unionImp2, *expectedUnionImp, *comm);
    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "*** imp2->setUnion(imp1) != expected Import ***" << endl;
    }
  }

  //
  // INSTANTIATIONS (template tests must be instantiated in the same
  // anonymous namespace as where the tests were defined)
  //

#define UNIT_TEST_GROUP(LOCAL_ORDINAL, GLOBAL_ORDINAL, NODE_TYPE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( ImportUnion, ContigPlusContig, LOCAL_ORDINAL, GLOBAL_ORDINAL, NODE_TYPE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LGN( UNIT_TEST_GROUP )

} // namespace (anonymous)





