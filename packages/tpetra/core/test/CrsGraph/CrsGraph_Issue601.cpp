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

#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Map.hpp"
#include "Teuchos_CommHelpers.hpp" // REDUCE_MIN, reduceAll
#include "Tpetra_TestingUtilities.hpp"

namespace { // (anonymous)

  using Tpetra::TestingUtilities::getDefaultComm;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using Teuchos::tuple;
  using std::endl;
  typedef Tpetra::global_size_t GST;

  //
  // UNIT TESTS
  //

  // See Github Issue #601.  This test is really only useful if run with > 2 MPI processes.
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, GlobalAssembleOverlappingRowMap, LO, GO, NT )
  {
    typedef Tpetra::CrsGraph<LO, GO, NT> crs_graph_type;
    typedef Tpetra::Map<LO, GO, NT> map_type;
    int lclSuccess = 1;
    int gblSuccess = 0;

    out << "Test Tpetra::CrsGraph::globalAssemble with overlapping row Map "
      "(see Github Issue #601)" << endl;
    Teuchos::OSTab tab1 (out);

    auto comm = getDefaultComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    TEST_ASSERT( numProcs > 2 );
    if (! success) {
      out << "This test is not useful unless you run it with 3 or more MPI "
        "processes!" << endl;
      return;
    }
    out << "Create overlapping row Map" << endl;

    // Each process r gets 2 rows: r, and (r+1) % numProcs.
    // The resulting row Map is overlapping, because numProcs > 2.

    const GST INV = Teuchos::OrdinalTraits<GST>::invalid ();
    Teuchos::Array<GO> myGblRowInds (2);
    myGblRowInds[0] = static_cast<GO> (myRank);
    myGblRowInds[1] = static_cast<GO> ((myRank + 1) % numProcs);
    const GO indexBase = 0;
    RCP<const map_type> rowMap =
      rcp (new map_type (INV, myGblRowInds (), indexBase, comm));

#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    const Tpetra::ProfileType profileTypes[2] = {Tpetra::DynamicProfile, Tpetra::StaticProfile};
#else
    const Tpetra::ProfileType profileTypes[1] = {Tpetra::StaticProfile};
#endif // TPETRA_ENABLE_DEPRECATED_CODE
    bool insertLocalEntryValues[] = { true, false };

    // Test both dynamic and static profile.
    for (Tpetra::ProfileType profileType : profileTypes) {

#ifdef TPETRA_ENABLE_DEPRECATED_CODE
      out << "ProfileType: "
          << ((profileType == Tpetra::StaticProfile) ? "Static" : "Dynamic")
          << "Profile" << endl;
#endif // TPETRA_ENABLE_DEPRECATED_CODE
      Teuchos::OSTab tab2 (out);

      // Test both with no local entries before globalAssemble(), and
      // with some (one) local entries before globalAssemble().  This
      // exercises lazy allocation of entries, if that's what CrsGraph
      // does (it does this as of 09 Sep 2016, though I've been aiming
      // to get rid of lazy allocation for years).
      for (bool insertLocalEntry : insertLocalEntryValues) {
        out << "insertLocalEntry: " << (insertLocalEntry ? "true" : "false")
            << endl;
        Teuchos::OSTab tab3 (out);

        const size_t maxNumEntPerRow = static_cast<size_t> (insertLocalEntry ? 2 : 1);
        crs_graph_type G (rowMap, maxNumEntPerRow, profileType);

        const GO gblRow0 = static_cast<GO> (myRank);
        const GO gblRow1 = static_cast<GO> ((myRank + 1) % numProcs);

        if (insertLocalEntry) {
          // Insert the diagonal entry in each of the rows on this
          // process.  This results in duplication across processes.
          // globalAssemble() should do the right thing in that case.
          G.insertGlobalIndices (gblRow0, tuple<GO> (gblRow0));
          G.insertGlobalIndices (gblRow1, tuple<GO> (gblRow1));
        }
        // Insert an entry into row (myRank - 1) % numProcs.  Add
        // numProcs to (myRank - 1) before taking the mod, so it
        // doesn't turn out negative.  (C(++) doesn't promise that the
        // result of mod is nonnegative if the input is negative.)
        {
          const GO nonlocalGblRow =
            static_cast<GO> ((numProcs + myRank - 1) % numProcs);
          // Make sure the new column index is not a duplicate (on the
          // process that will receive it at globalAssemble).  This is
          // very much NOT a square graph (/ matrix).
          const GO gblCol = static_cast<GO> (numProcs) + nonlocalGblRow;
          TEST_NOTHROW( G.insertGlobalIndices (nonlocalGblRow, tuple<GO> (gblCol)) );
        }

        out << "Call G.globalAssemble()" << endl;
        // This is the moment of truth.
        TEST_NOTHROW( G.globalAssemble () );

        lclSuccess = success ? 1 : 0;
        gblSuccess = 0;
        reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
        TEST_EQUALITY_CONST( gblSuccess, 1 );
        if (! success) {
          out << "Test FAILED on some process; returning early." << endl;
          return;
        }

        out << "Now make sure each of our process' rows are correct" << endl;

        // Test gblRow0
        {
          Teuchos::ArrayView<const GO> gblInds;
          G.getGlobalRowView (gblRow0, gblInds);

          const LO expectedNumEnt = static_cast<LO> (maxNumEntPerRow);
          TEST_EQUALITY( static_cast<LO> (gblInds.size ()), expectedNumEnt );
          if (static_cast<LO> (gblInds.size ()) == expectedNumEnt) {
            if (insertLocalEntry) {
              auto lclEntIter = std::find (gblInds.begin (), gblInds.end (), gblRow0);
              TEST_ASSERT( lclEntIter != gblInds.end () );
            }
            const GO gblCol0 = gblRow0 + static_cast<GO> (numProcs);
            auto nonlclEntIter = std::find (gblInds.begin (), gblInds.end (), gblCol0);
            TEST_ASSERT( nonlclEntIter != gblInds.end () );
          }
        }

        // Test gblRow1
        {
          Teuchos::ArrayView<const GO> gblInds;
          G.getGlobalRowView (gblRow1, gblInds);

          const LO expectedNumEnt = static_cast<LO> (maxNumEntPerRow);
          TEST_EQUALITY( static_cast<LO> (gblInds.size ()), expectedNumEnt );
          if (static_cast<LO> (gblInds.size ()) == expectedNumEnt) {
            if (insertLocalEntry) {
              auto lclEntIter = std::find (gblInds.begin (), gblInds.end (), gblRow1);
              TEST_ASSERT( lclEntIter != gblInds.end () );
            }
            const GO gblCol1 = gblRow1 + static_cast<GO> (numProcs);
            auto nonlclEntIter = std::find (gblInds.begin (), gblInds.end (), gblCol1);
            TEST_ASSERT( nonlclEntIter != gblInds.end () );
          }
        }

        // Test across processes.
        lclSuccess = success ? 1 : 0;
        gblSuccess = 0;
        reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
        TEST_EQUALITY_CONST( gblSuccess, 1 );
        if (! success) {
          out << "Test FAILED on some process; returning early." << endl;
          return;
        }
      } // insertLocalEntry
    } // profileType
  }

  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP( LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, GlobalAssembleOverlappingRowMap, LO, GO, NT )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LGN( UNIT_TEST_GROUP )

} // namespace (anonymous)
