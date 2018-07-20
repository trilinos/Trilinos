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

#include <Tpetra_TestingUtilities.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_Tuple.hpp>
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Distributor.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Import_Util.hpp"
#include "Tpetra_Import_Util2.hpp"
#include "Tpetra_Details_packCrsMatrix.hpp"
#include "Tpetra_Details_unpackCrsMatrixAndCombine.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"
#include "TpetraExt_MatrixMatrix.hpp"
#include "Tpetra_Details_gathervPrint.hpp"

namespace {

  // Get a Teuchos::ArrayView which views the host Kokkos::View of the
  // input 1-D Kokkos::DualView.
  template<class DualViewType>
  Teuchos::ArrayView<typename DualViewType::t_dev::value_type>
  getArrayViewFromDualView (const DualViewType& x)
  {
    static_assert (static_cast<int> (DualViewType::t_dev::rank) == 1,
                   "The input DualView must have rank 1.");
    TEUCHOS_TEST_FOR_EXCEPTION
      (x.template need_sync<Kokkos::HostSpace> (), std::logic_error, "The "
       "input Kokkos::DualView was most recently modified on device, but this "
       "function needs the host view of the data to be the most recently "
       "modified.");

    auto x_host = x.template view<Kokkos::HostSpace> ();
    typedef typename DualViewType::t_dev::value_type value_type;
    return Teuchos::ArrayView<value_type> (x_host.data (),
                                           x_host.extent (0));
  }

  using Teuchos::as;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::arrayViewFromVector;
  using Teuchos::Comm;
  using Teuchos::FancyOStream;
  using Teuchos::getFancyOStream;
  using Teuchos::includesVerbLevel;
  using Teuchos::OrdinalTraits;
  using Teuchos::OSTab;
  using Teuchos::outArg;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::REDUCE_MAX;
  using Teuchos::reduceAll;
  using Teuchos::ScalarTraits;
  using Teuchos::tuple;

  using Tpetra::ADD;
  using Tpetra::createContigMap;
  using Tpetra::CrsGraph;
  using Tpetra::CrsMatrix;
  using Tpetra::DynamicProfile;
  using Tpetra::Export;
  using Tpetra::Import;
  using Tpetra::INSERT;
  using Tpetra::Map;
  using Tpetra::REPLACE;
  using Tpetra::StaticProfile;

  using std::cerr;
  using std::cout;
  using std::ostream_iterator;
  using std::endl;

  using Node = Tpetra::Map<>::node_type;

  // Command-line argument values (initially set to defaults).
  bool testMpi = true;
  double errorTolSlack = 1e+1;
  std::string distributorSendType ("Send");
  bool barrierBetween = true;
  bool verbose = false;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
    // FIXME (mfh 02 Apr 2012) It would be better to ask Distributor
    // for the list of valid send types, but setOption() needs a const
    // char[], not an std::string.
    clp.setOption ("distributor-send-type", &distributorSendType,
                   "In MPI tests, the type of send operation that the Tpetra::"
                   "Distributor will use.  Valid values include \"Isend\", "
                   "\"Rsend\", \"Send\", and \"Ssend\".");
    clp.setOption ("barrier-between", "no-barrier-between", &barrierBetween,
                   "In MPI tests, whether Tpetra::Distributor will execute a "
                   "barrier between posting receives and posting sends.");
    clp.setOption ("verbose", "quiet", &verbose, "Whether to print verbose "
                   "output.");
  }

  RCP<const Comm<int> >
  getDefaultComm()
  {
    if (testMpi) {
      return Tpetra::getDefaultComm ();
    }
    else {
      return rcp (new Teuchos::SerialComm<int> ());
    }
  }

  RCP<ParameterList>
  getDistributorParameterList ()
  {
    static RCP<ParameterList> plist; // NOT THREAD SAFE, but that's OK here.
    if (plist.is_null ()) {
      plist = parameterList ("Tpetra::Distributor");
      plist->set ("Send type", distributorSendType);
      plist->set ("Barrier between receives and sends", barrierBetween);

      if (verbose && getDefaultComm()->getRank() == 0) {
        cout << "ParameterList for Distributor: " << *plist << endl;
      }

      if (verbose) {
        // Tell Distributor to print verbose output.
        Teuchos::VerboseObject<Tpetra::Distributor>::setDefaultVerbLevel (Teuchos::VERB_EXTREME);
      }
    }
    return plist;
  }

  RCP<ParameterList> getImportParameterList () {
    return getDistributorParameterList (); // For now.
  }

  RCP<ParameterList> getExportParameterList () {
    return parameterList (* (getDistributorParameterList ())); // For now.
  }

  RCP<ParameterList> getCrsGraphParameterList () {
    static RCP<ParameterList> plist; // NOT THREAD SAFE, but that's OK here.
    if (plist.is_null ()) {
      plist = parameterList ("Tpetra::CrsGraph");
      plist->set ("Import", * (getImportParameterList ()));
      plist->set ("Export", * (getExportParameterList ()));

      if (verbose && getDefaultComm()->getRank() == 0) {
        cout << "ParameterList for CrsGraph: " << *plist << endl;
      }
    }
    return plist;
  }

  RCP<ParameterList> getCrsMatrixParameterList () {
    return parameterList (* (getCrsGraphParameterList ())); // For now.
  }



  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsGraphImportExport, doImport, LO, GO )
  {
    using Teuchos::VERB_EXTREME;
    using Teuchos::VERB_NONE;
    typedef Tpetra::global_size_t GST;

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const RCP<const Comm<int> > comm = getDefaultComm();
    const int numProcs = comm->getSize();
    const int myRank = comm->getRank();
    if (numProcs < 2) {
      out << "This test is only nontrivial if run with multiple MPI "
        "processes, but you have run it with only 1 process." << endl;
      return;
    }

    // Prepare for verbose output, if applicable.
    Teuchos::EVerbosityLevel verbLevel = verbose ? VERB_EXTREME : VERB_NONE;
    const bool doPrint = includesVerbLevel (verbLevel, VERB_EXTREME, true);
    if (doPrint) {
      out << "CrsGraphImportExport unit test" << endl;
    }
    OSTab tab1 (out); // Add one tab level

    std::ostringstream err;
    int lclErr = 0;

    out << "Create a Map that is evenly distributed, "
      "and another that has all elements on Proc 0" << endl;
    try {
      OSTab tab2 (out);
      const int tgt_num_local_elements = 3;
      const int src_num_local_elements =
        (myRank == 0 ? numProcs*tgt_num_local_elements : 0);

      // create Maps
      if (doPrint) {
        out << "Creating source and target Maps" << endl;
      }
      RCP<const Map<LO, GO> > src_map =
        createContigMap<LO, GO> (INVALID, src_num_local_elements, comm);
      RCP<const Map<LO, GO> > tgt_map =
        createContigMap<LO, GO> (INVALID, tgt_num_local_elements, comm);

      // create CrsGraph objects
      if (doPrint) {
        out << "Creating source and target CrsGraphs" << endl;
      }
      RCP<CrsGraph<LO, GO> > src_graph =
        rcp (new CrsGraph<LO, GO> (src_map, 1, DynamicProfile,
                                   getCrsGraphParameterList ()));
      RCP<CrsGraph<LO, GO> > tgt_graph =
        rcp (new CrsGraph<LO, GO> (tgt_map, 1, DynamicProfile,
                                   getCrsGraphParameterList ()));

      // Create a simple diagonal source graph.
      if (doPrint) {
        out << "Filling source CrsGraph" << endl;
      }
      Array<GO> diag(1);
      LO row = 0;
      for (size_t i = 0; i < src_map->getNodeNumElements (); ++i, ++row) {
        const GO globalrow = src_map->getGlobalElement (row);
        diag[0] = globalrow;
        src_graph->insertGlobalIndices (globalrow, diag ());
      }

      // Import from the source graph to the target graph.
      if (doPrint) {
        out << "Importing from source to target CrsGraph" << endl;
      }
      Import<LO, GO> importer (src_map, tgt_map, getImportParameterList ());
      tgt_graph->doImport (*src_graph, importer, INSERT);
      tgt_graph->fillComplete ();

      // Loop through the target graph and make sure it is diagonal.
      if (doPrint) {
        out << "Verifying target CrsGraph" << endl;
      }
      row = 0;
      for (size_t i = 0; i < tgt_map->getNodeNumElements (); ++i, ++row) {
        ArrayView<const LO> rowview;
        tgt_graph->getLocalRowView( row, rowview );
        TEST_EQUALITY(rowview.size(), 1);
        TEST_EQUALITY(rowview[0], row);
      }
    } catch (std::exception& e) {
      err << "Proc " << myRank << ": " << e.what () << endl;
      lclErr = 1;
    }

    int gblErr = 0;
    reduceAll<int, int> (*comm, REDUCE_MAX, lclErr, outArg (gblErr));
    TEST_EQUALITY_CONST( gblErr, 0 );
    if (gblErr != 0) {
      Tpetra::Details::gathervPrint (out, err.str (), *comm);
      out << "Above test failed; aborting further tests" << endl;
      return;
    }

    out << "Test with even number of processes (skip test otherwise)" << endl;
    try {
      OSTab tab2 (out);
      if (numProcs % 2 == 0) {
        // Create Maps that are distributed differently but have the
        // same global number of elements. The source map will have 3
        // elements on even-numbered processes and 5 on odd-numbered
        // processes. The target map will have 4 elements on each
        // process.
        LO src_num_local = 5;
        if (myRank % 2 == 0) {
          src_num_local = 3;
        }
        LO tgt_num_local = 4;

        RCP<const Map<LO, GO> > src_map =
          createContigMap<LO, GO> (INVALID, src_num_local, comm);
        RCP<const Map<LO, GO> > tgt_map =
          createContigMap<LO, GO> (INVALID, tgt_num_local, comm);

        RCP<CrsGraph<LO, GO> > src_graph =
          rcp (new CrsGraph<LO, GO> (src_map, 24, DynamicProfile,
                                     getCrsGraphParameterList ()));
        RCP<CrsGraph<LO, GO> > tgt_graph =
          rcp (new CrsGraph<LO, GO> (tgt_map, 24, DynamicProfile,
                                     getCrsGraphParameterList ()));

        // This time make src_graph be a full lower-triangular graph.
        // Each row of column indices will have length 'globalrow'+1,
        // and contain column indices 0 .. 'globalrow'.
        Array<GO> cols(1);
        for (GO globalrow = src_map->getMinGlobalIndex ();
             globalrow <= src_map->getMaxGlobalIndex (); ++globalrow) {
          cols.resize(globalrow+1);
          for (GO col = 0; col < globalrow+1; ++col) {
            cols[col] = col;
          }
          src_graph->insertGlobalIndices (globalrow, cols ());
        }

        Import<LO, GO> importer (src_map, tgt_map, getImportParameterList ());
        tgt_graph->doImport (*src_graph, importer, INSERT);

        src_graph->fillComplete ();
        tgt_graph->fillComplete ();

        // Loop through tgt_graph and make sure that each row has length
        // globalrow+1 and has the correct contents.
        RCP<const Map<LO, GO> > colmap = tgt_graph->getColMap ();

        for (GO globalrow = tgt_map->getMinGlobalIndex ();
             globalrow <= tgt_map->getMaxGlobalIndex ();
             ++globalrow) {
          LO localrow = tgt_map->getLocalElement (globalrow);
          ArrayView<const LO> rowview;
          tgt_graph->getLocalRowView (localrow, rowview);
          TEST_EQUALITY(rowview.size(), globalrow+1);

          // The target graph doesn't necessarily promise sorted
          // order.  Thus, we copy out the local row view, convert to
          // global, sort, and test the sorted result.
          Array<GO> curEntries (rowview.size ());
          for (GO j = static_cast<GO> (0); j < static_cast<GO> (rowview.size ()); ++j) {
            curEntries[j] = colmap->getGlobalElement(rowview[j]);
          }
          std::sort (curEntries.begin (), curEntries.end ());

          for (GO j = static_cast<GO> (0); j < static_cast<GO> (globalrow + 1); ++j) {
            TEST_EQUALITY(curEntries[j], j);
          }
        }
      }
    } catch (std::exception& e) {
      err << "Proc " << myRank << ": " << e.what () << endl;
      lclErr = 1;
    }

    gblErr = 0;
    reduceAll<int, int> (*comm, REDUCE_MAX, lclErr, outArg (gblErr));
    TEST_EQUALITY_CONST( gblErr, 0 );
    if (gblErr != 0) {
      Tpetra::Details::gathervPrint (out, err.str (), *comm);
      out << "Above test failed; aborting further tests" << endl;
      return;
    }
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrixImportExport, doImport, LO, GO, Scalar )
  {
    typedef Tpetra::global_size_t GST;
    typedef Map<LO, GO> map_type;

    out << "(CrsMatrixImportExport,doImport) test" << endl;
    OSTab tab1 (out); // Add one tab level

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();

    // Parameters for CrsMatrix.  We'll let all CrsMatrix instances
    // use the same parameter list.
    RCP<ParameterList> crsMatPlist = getCrsMatrixParameterList ();

    if (numImages < 2) {
      out << "This test is only meaningful if running with multiple MPI "
        "processes, but you ran it with only 1 process." << endl;
      return;
    }

    std::ostringstream err;
    int lclErr = 0;
    int gblErr = 0;

    out << "First test: Import a diagonal CrsMatrix from a source row Map "
      "that has all indices on Process 0, to a target row Map that is "
      "uniformly distributed over processes." << endl;
    try {
      OSTab tab2 (out);
      const GO indexBase = 0;
      const LO tgt_num_local_elements = 3;
      const LO src_num_local_elements = (myImageID == 0) ?
        static_cast<LO> (numImages*tgt_num_local_elements) :
        static_cast<LO> (0);

      // Create row Maps for the source and target
      // Create Maps.
      RCP<const map_type> src_map =
        rcp (new map_type (INVALID,
                           static_cast<size_t> (src_num_local_elements),
                           indexBase, comm));
      RCP<const map_type> tgt_map =
        rcp (new map_type (INVALID,
                           static_cast<size_t> (tgt_num_local_elements),
                           indexBase, comm));

      // Create CrsMatrix objects.
      RCP<CrsMatrix<Scalar, LO, GO> > src_mat =
        rcp (new CrsMatrix<Scalar, LO, GO> (src_map, 1, StaticProfile,
                                            crsMatPlist));
      RCP<CrsMatrix<Scalar, LO, GO> > tgt_mat =
        rcp (new CrsMatrix<Scalar, LO, GO> (tgt_map, 1, StaticProfile,
                                            crsMatPlist));

      // Create a simple diagonal source graph.
      if (src_num_local_elements != 0) {
        for (GO globalrow = src_map->getMinGlobalIndex();
             globalrow <= src_map->getMaxGlobalIndex();
             ++globalrow) {
          src_mat->insertGlobalValues (globalrow,
                                       tuple<GO> (globalrow),
                                       tuple<Scalar> (globalrow));
        }
      }
      src_mat->fillComplete ();

      // Create the importer
      Import<LO, GO> importer (src_map, tgt_map, getImportParameterList ());
      // Do the import, and fill-complete the target matrix.
      tgt_mat->doImport (*src_mat, importer, INSERT);
      tgt_mat->fillComplete ();

      // Loop through tgt_mat and make sure it is diagonal.
      if (tgt_num_local_elements != 0) {
        for (GO gblRow = tgt_map->getMinGlobalIndex ();
             gblRow <= tgt_map->getMaxGlobalIndex ();
             ++gblRow) {
          const LO lclRow = tgt_map->getLocalElement (gblRow);

          ArrayView<const LO> lclInds;
          ArrayView<const Scalar> lclVals;
          tgt_mat->getLocalRowView (lclRow, lclInds, lclVals);
          TEST_EQUALITY_CONST(lclInds.size(), 1);
          TEST_EQUALITY_CONST(lclVals.size(), 1);

          if (lclInds.size () != 0) { // don't segfault in error case
            TEST_EQUALITY(tgt_mat->getColMap ()->getGlobalElement (lclInds[0]), gblRow);
          }
          if (lclVals.size () != 0) { // don't segfault in error case
            TEST_EQUALITY(lclVals[0], as<Scalar> (gblRow));
          }
        }
      }

      // Test the all-in-one import and fill complete nonmember
      // constructor.  The returned matrix should also be diagonal and
      // should equal tgt_mat.
      Teuchos::ParameterList dummy;
      typedef CrsMatrix<Scalar, LO, GO> crs_type;
      RCP<crs_type> A_tgt2 =
        Tpetra::importAndFillCompleteCrsMatrix<crs_type> (src_mat, importer,
                                                          Teuchos::null,
                                                          Teuchos::null,
                                                          rcp(&dummy,false));

      // Make sure that A_tgt2's row Map is the same as tgt_map, and
      // is also the same as the Import's targetMap.  They should have
      // the same Map: not just in the sense of Map::isSameAs(), but
      // also in the sense of pointer equality.  (A Map isSameAs
      // itself, so we only need to test for pointer equality.)
      TEST_EQUALITY(A_tgt2->getRowMap ().getRawPtr (),
                    tgt_map.getRawPtr ());
      TEST_EQUALITY(A_tgt2->getRowMap ().getRawPtr (),
                    importer.getTargetMap ().getRawPtr ());

      // Loop through A_tgt2 and make sure each row has the same
      // entries as tgt_mat.  In the fully general case, the
      // redistribution may have added together values, resulting in
      // small rounding errors.  This is why we use an error tolerance
      // (with a little bit of wiggle room).
      typedef typename ScalarTraits<Scalar>::magnitudeType magnitude_type;
      // Include a little wiggle room in the error tolerance.  It
      // would be smarter to use an a posteriori error bound, but we
      // can't get inside the Import to see which values it's adding
      // together, so we make a rough guess.
      const magnitude_type tol =
          as<magnitude_type> (10) * ScalarTraits<magnitude_type>::eps ();

      Array<LO> tgtRowInds;
      Array<Scalar>  tgtRowVals;
      Array<LO> tgt2RowInds;
      Array<Scalar>  tgt2RowVals;
      for (LO localrow = tgt_map->getMinLocalIndex();
           localrow <= tgt_map->getMaxLocalIndex();
           ++localrow)
      {
        size_t tgtNumEntries = tgt_mat->getNumEntriesInLocalRow (localrow);
        size_t tgt2NumEntries = tgt_mat->getNumEntriesInLocalRow (localrow);

        // Same number of entries in each row?
        TEST_EQUALITY(tgtNumEntries, tgt2NumEntries);

        if (tgtNumEntries > as<size_t> (tgtRowInds.size ())) {
          tgtRowInds.resize (tgtNumEntries);
          tgtRowVals.resize (tgtNumEntries);
        }
        if (tgt2NumEntries > as<size_t> (tgt2RowInds.size ())) {
          tgt2RowInds.resize (tgt2NumEntries);
          tgt2RowVals.resize (tgt2NumEntries);
        }
        tgt_mat->getLocalRowCopy (localrow, tgtRowInds(), tgtRowVals(), tgtNumEntries);
        A_tgt2->getLocalRowCopy (localrow, tgt2RowInds(), tgt2RowVals(), tgt2NumEntries);

        // Entries should be sorted, but let's sort them by column
        // index just in case.  This is why we got a row copy instead
        // of a row view.
        Tpetra::sort2 (tgtRowInds.begin(), tgtRowInds.end(), tgtRowVals.begin());
        Tpetra::sort2 (tgt2RowInds.begin(), tgt2RowInds.end(), tgt2RowVals.begin());

        // Now that the entries are sorted, compare to make sure they
        // have the same column indices and values.  In the fully
        // general case, the redistribution may have added together
        // values, resulting in small rounding errors.
        typedef typename Array<Scalar>::size_type size_type;
        for (size_type k = 0; k < static_cast<size_type> (tgtNumEntries); ++k) {
          TEST_EQUALITY(tgtRowInds[k], tgt2RowInds[k]);
          // The "out" and "success" variables should have been
          // automatically defined by the unit test framework, in case
          // you're wondering where they came from.
          TEUCHOS_TEST_FLOATING_EQUALITY(tgtRowVals[k], tgt2RowVals[k], tol, out, success);
        } // for each entry in the current row
      } // for each row in the matrix
    }
    catch (std::exception& e) { // end of the first test
      err << "Proc " << myImageID << ": " << e.what () << endl;
      lclErr = 1;
    }

    reduceAll<int, int> (*comm, REDUCE_MAX, lclErr, outArg (gblErr));
    TEST_EQUALITY_CONST( gblErr, 0 );
    if (gblErr != 0) {
      Tpetra::Details::gathervPrint (out, err.str (), *comm);
      out << "Above test failed; aborting further tests" << endl;
      return;
    }

    out << "Test with even number of processes (skip test otherwise)" << endl;
    try {
      OSTab tab2 (out);
      if (numImages%2 == 0) {
        // Create Maps that are distributed differently but have the
        // same global number of elements. The source-map will have 3
        // elements on even-numbered processes and 5 on odd-numbered
        // processes. The target-map will have 4 elements on each
        // process.
        const LO src_num_local = (myImageID%2 == 0 ? 3 : 5);
        const LO tgt_num_local = 4;

        RCP<const map_type> src_map =
          createContigMap<LO, GO> (INVALID, src_num_local, comm);
        RCP<const map_type> tgt_map =
          createContigMap<LO, GO> (INVALID, tgt_num_local, comm);

        RCP<CrsMatrix<Scalar, LO, GO> > src_mat =
          rcp (new CrsMatrix<Scalar, LO, GO> (src_map, 24, DynamicProfile, crsMatPlist));
        RCP<CrsMatrix<Scalar, LO, GO> > tgt_mat =
          rcp (new CrsMatrix<Scalar, LO, GO> (tgt_map, 24, DynamicProfile, crsMatPlist));

        // This time make src_mat a full lower-triangular matrix.  Each
        // row of column-indices will have length 'globalrow', and
        // contain column-indices 0 .. 'globalrow'-1
        Array<GO> cols(1);
        Array<Scalar>  vals(1);
        for (GO globalrow = src_map->getMinGlobalIndex();
             globalrow <= src_map->getMaxGlobalIndex();
             ++globalrow) {
          if (globalrow > 0) {
            cols.resize(globalrow);
            vals.resize(globalrow);
            for (GO col=0; col<globalrow; ++col) {
              cols[col] = as<GO>(col);
              vals[col] = as<Scalar>(col);
            }
            src_mat->insertGlobalValues (globalrow, cols (), vals ());
          }
        }

        Import<LO, GO> importer (src_map, tgt_map, getImportParameterList ());
        tgt_mat->doImport(*src_mat, importer, Tpetra::INSERT);
        tgt_mat->fillComplete();

        // now we're going to loop through tgt_mat and make sure that
        // each row has length 'globalrow' and has the correct contents:
        const Teuchos::RCP<const map_type> colmap = tgt_mat->getColMap();

        for (GO globalrow=tgt_map->getMinGlobalIndex();
             globalrow<=tgt_map->getMaxGlobalIndex(); ++globalrow) {
          LO localrow = tgt_map->getLocalElement(globalrow);
          ArrayView<const LO> rowinds;
          ArrayView<const Scalar> rowvals;
          tgt_mat->getLocalRowView(localrow, rowinds, rowvals);
          TEST_EQUALITY(rowinds.size(), globalrow);
          TEST_EQUALITY(rowvals.size(), globalrow);

          // The target graph doesn't necessarily promise sorted
          // order.  Thus, we copy out the local row view, convert to
          // global, sort, and test the sorted result.  Since this is
          // a matrix, not just a graph, we have to sort the two views
          // jointly -- that is, sort curInds and apply the resulting
          // permutation to curVals.
          Array<GO> curInds (rowinds.size ());
          Array<Scalar> curVals (rowvals.size ());

          for (decltype (rowinds.size()) j=0; j<rowinds.size(); ++j) {
            curInds[j] = colmap->getGlobalElement (rowinds[j]);
            curVals[j] = rowvals[j];
          }
          Tpetra::sort2 (curInds.begin (), curInds.end (), curVals.begin ());

          for (decltype (rowinds.size()) j=0; j<rowinds.size(); ++j) {
            TEST_EQUALITY( curInds[j], as<GO> (j) );
            TEST_EQUALITY( curVals[j], as<Scalar> (j)  );
          }
        }
      }
    } catch (std::exception& e) {
      err << "Proc " << myImageID << ": " << e.what () << endl;
      lclErr = 1;
    }

    reduceAll<int, int> (*comm, REDUCE_MAX, lclErr, outArg (gblErr));
    TEST_EQUALITY_CONST( gblErr, 0 );
    if (gblErr != 0) {
      Tpetra::Details::gathervPrint (out, err.str (), *comm);
      out << "Above test failed; aborting further tests" << endl;
      return;
    }
  }


template<class Graph>
bool graphs_are_same(const RCP<Graph>& G1, const RCP<const Graph>& G2)
{
  typedef typename Graph::local_ordinal_type LO;

  int my_rank = G1->getRowMap()->getComm()->getRank();

  // Make sure each graph is fill complete before checking other properties
  if (! G1->isFillComplete()) {
    if (my_rank == 0)
      cerr << "***Error: Graph 1 is not fill complete!" << endl;
    return false;
  }
  if (! G2->isFillComplete()) {
    if (my_rank == 0)
      cerr << "***Error: Graph 2 is not fill complete!" << endl;
    return false;
  }

  int errors = 0;

  if (! G1->getRowMap()->isSameAs(*G2->getRowMap())) {
    if (my_rank == 0)
      cerr << "***Error: Graph 1's row map is different than Graph 2's" << endl;
    errors++;
  }
  if (! G1->getDomainMap()->isSameAs(*G2->getDomainMap())) {
    if (my_rank == 0)
      cerr << "***Error: Graph 1's domain map is different than Graph 2's" << endl;
    errors++;
  }
  if (! G1->getRangeMap()->isSameAs(*G2->getRangeMap())) {
    if (my_rank == 0)
      cerr << "***Error: Graph 1's range map is different than Graph 2's" << endl;
    errors++;
  }
  if (G1->getNodeNumEntries() != G2->getNodeNumEntries()) {
    cerr << "***Error: Graph 1 does not have the same number of entries as Graph 2 on Process "
         << my_rank << endl;
    errors++;
  }

  if (errors != 0) return false;

  for (LO i=0; i<static_cast<LO>(G1->getNodeNumRows()); i++) {
    ArrayView<const LO> V1, V2;
    G1->getLocalRowView(i, V1);
    G2->getLocalRowView(i, V2);
    if (V1.size() != V2.size()) {
      cerr << "***Error: Graph 1 and Graph 2 have different number of entries in local row "
           << i << " on Process " << my_rank << endl;
      errors++;
      continue;
    }
    int jerr = 0;
    for (LO j=0; static_cast<LO>(j<V1.size()); j++) {
      if (V1[j] != V2[j])
        jerr++;
    }
    if (jerr != 0) {
      cerr << "***Error: One or more entries in row " << i << " on Process " << my_rank
           << " Graphs 1 and 2 are not the same" << endl;
      errors++;
      continue;
    }
  }

  return (errors == 0);

}

// All the fused Import/export test stuff
// ===============================================================================
template <class MatrixType, class ImportType>
void build_matrix_unfused_import(const MatrixType & SourceMatrix, ImportType & RowImporter, RCP<MatrixType> & A){
  A->doImport(SourceMatrix, RowImporter, Tpetra::INSERT);
  A->fillComplete(SourceMatrix.getDomainMap(), SourceMatrix.getRangeMap());
}

// ===============================================================================
template <class CrsMatrixType>
double test_with_matvec(const CrsMatrixType &A, const CrsMatrixType &B){
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef typename CrsMatrixType::local_ordinal_type LO;
  typedef typename CrsMatrixType::global_ordinal_type GO;
  typedef typename CrsMatrixType::scalar_type Scalar;
  typedef typename CrsMatrixType::node_type NT;

  typedef Tpetra::Map<LO,GO, NT> map_type;
  typedef Tpetra::Vector<Scalar, LO, GO, NT> vector_type;
  typedef Tpetra::Import<LO, GO, NT> import_type;
  typedef Tpetra::Export<LO, GO, NT> export_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  RCP<const map_type> Xamap  = A.getDomainMap();
  RCP<const map_type> Yamap  = A.getRangeMap();
  RCP<const map_type> Xbmap  = B.getDomainMap();
  RCP<const map_type> Ybmap  = B.getRangeMap();

  vector_type Xa(Xamap), Xb(Xbmap),Ya(Yamap), Yb(Ybmap);

  // Start with a generic vector which has a 1-1 map
  const GO indexBase = Teuchos::OrdinalTraits<GO>::zero ();
  RCP<const map_type> Xgmap =
    rcp (new map_type (Xamap->getGlobalNumElements (), indexBase, Xamap->getComm ()));
  RCP<const map_type> Ygmap =
    rcp (new map_type (Yamap->getGlobalNumElements (), indexBase, Yamap->getComm ()));
  vector_type X0(Xgmap), Y0a(Ygmap), Diff(Ygmap);
  Teuchos::ScalarTraits< Scalar >::seedrandom(24601);
  X0.putScalar (STS::one ());

  // Handle domain map change
  if (! Xgmap->isSameAs (*Xamap)) {
    import_type Ximport (Xgmap, Xamap);
    Xa.doImport (X0, Ximport, Tpetra::INSERT);
  } else {
    Xa = X0;
  }

  if (! Xgmap->isSameAs (*Xbmap)) {
    import_type Ximport (Xgmap, Xbmap);
    Xb.doImport (X0, Ximport, Tpetra::INSERT);
  } else {
    Xb = X0;
  }

  Xa.putScalar (STS::one ());
  Xb.putScalar (STS::one ());

  // Do the multiplies
  A.apply(Xa,Ya);
  B.apply(Xb,Yb);

  // Handle range Map change
  if (! Ygmap->isSameAs (*Yamap)) {
    export_type Yexport (Yamap, Ygmap);
    Y0a.doExport (Ya, Yexport, Tpetra::ADD);
  } else {
    Y0a = Ya;
  }

  if (! Ygmap->isSameAs (*Ybmap)) {
    export_type Yexport (Ybmap, Ygmap);
    Diff.doExport (Yb, Yexport, Tpetra::ADD);
  } else {
    Diff = Yb;
  }

  // Check solution
  Diff.update (-STS::one (), Y0a, STS::one ());
  Teuchos::Array<typename STS::magnitudeType> norms (1);
  Diff.norm2 (norms);
  return Teuchos::as<double> (norms[0]);
}


// ===============================================================================
template <class CrsMatrixType, class map_type>
double test_with_matvec_reduced_maps(const CrsMatrixType &A, const CrsMatrixType &B, const map_type & Bfullmap){
  typedef typename CrsMatrixType::local_ordinal_type LO;
  typedef typename CrsMatrixType::global_ordinal_type GO;
  typedef typename CrsMatrixType::scalar_type Scalar;
  typedef typename CrsMatrixType::node_type NT;
  typedef Tpetra::MultiVector<Scalar, LO, GO, NT> vector_type;
  typedef Tpetra::Import<LO, GO, NT> import_type;

  RCP<const map_type>  Amap  = A.getDomainMap();
  vector_type Xa(Amap,1), Ya(Amap,1), Diff(Amap,1);
  RCP<const map_type> Bmap  = Bfullmap.getNodeNumElements() > 0 ? B.getDomainMap() : Teuchos::null;

  vector_type Xb_alias(rcp(&Bfullmap,false),1);
  vector_type Yb_alias(rcp(&Bfullmap,false),1);

  RCP<vector_type> Xb       = !Bmap.is_null() ? Xb_alias.offsetViewNonConst(Bmap,0) : Teuchos::null;
  RCP<vector_type> Yb       = !Bmap.is_null() ? Yb_alias.offsetViewNonConst(Bmap,0) : Teuchos::null;

  import_type Ximport(Amap,rcp(&Bfullmap,false));

  // Set the input vector
  Teuchos::ScalarTraits< Scalar >::seedrandom(24601);
  Xa.randomize();
  Xb_alias.doImport(Xa,Ximport,Tpetra::INSERT);

  // Do the multiplies
  A.apply(Xa,Ya);
  if(!Bmap.is_null()) B.apply(*Xb,*Yb);

  // Check solution
  import_type Yimport(rcp(&Bfullmap,false),Amap);
  Diff.doImport(Yb_alias,Yimport,Tpetra::INSERT);
  Diff.update(-1.0,Ya,1.0);
  Teuchos::Array< typename Teuchos::ScalarTraits< Scalar >::magnitudeType > norms(1);
  Diff.norm2(norms);

  return Teuchos::as<double>(norms[0]);
}


// ===============================================================================
template<class CrsMatrixType>
void build_test_matrix(RCP<const Teuchos::Comm<int> > & Comm, RCP<CrsMatrixType> & A){
  typedef typename CrsMatrixType::local_ordinal_type LO;
  typedef typename CrsMatrixType::global_ordinal_type GO;
  typedef typename CrsMatrixType::scalar_type Scalar;
  typedef typename CrsMatrixType::node_type NT;

  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  const Scalar ONE = STS::one ();
  const Scalar TWO = ONE + ONE;
  const int NumProc = Comm->getSize ();
  const int MyPID   = Comm->getRank ();

  // Case 1: Tridiagonal
  LO NumMyEquations = 100;
  GO NumGlobalEquations = (NumMyEquations * NumProc) + (NumProc < 3 ? NumProc : 3);
  if (MyPID < 3) {
    ++NumMyEquations;
  }

  // Construct a Map that puts approximately the same Number of equations on each processor
  RCP<const map_type> MyMap = rcp(new map_type(NumGlobalEquations, NumMyEquations, 0, Comm));

  // Create the matrix
  A = rcp(new CrsMatrixType(MyMap,0));

  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1
  Teuchos::Array<Scalar> Values(3);
  Values[0] = -ONE;
  Values[1] = -ONE;

  Teuchos::Array<GO> Indices(2);
  size_t NumEntries=1;

  for (LO i = 0; i < NumMyEquations; i++) {
    GO GID = MyMap->getGlobalElement(i);
    if(GID == 0){
      Indices[0] = 1;
      NumEntries = 1;
    }
    else if (GID == NumGlobalEquations-1) {
      Indices[0] = NumGlobalEquations-2;
      NumEntries = 1;
    }
    else {
      Indices[0] = GID-1;
      Indices[1] = GID+1;
      NumEntries = 2;
    }
    Values[0] = -ONE;
    Values[1] = -ONE;
    A->insertGlobalValues(GID,Indices.view(0,NumEntries),Values.view(0,NumEntries));

    Indices[0] = GID;
    Values[0] = TWO;
    A->insertGlobalValues(GID,Indices.view(0,1), Values.view(0,1));
  }

  A->fillComplete();
}

// ===============================================================================
template<class CrsMatrixType>
void build_test_matrix_wideband(RCP<const Teuchos::Comm<int> > & Comm, RCP<CrsMatrixType> & A){
  typedef typename CrsMatrixType::local_ordinal_type LO;
  typedef typename CrsMatrixType::global_ordinal_type GO;
  typedef typename CrsMatrixType::scalar_type Scalar;
  typedef typename CrsMatrixType::node_type NT;

  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  const Scalar ONE = STS::one ();
  const Scalar TWO = ONE + ONE;
  const int NumProc = Comm->getSize ();
  const int MyPID = Comm->getRank ();

  LO NumMyEquations = 1000;
  GO NumGlobalEquations = (NumMyEquations * NumProc) + (NumProc < 3 ? NumProc : 3);
  if (MyPID < 3) {
    ++NumMyEquations;
  }

  // Construct a Map that puts approximately the same Number of equations on each processor
  RCP<const map_type > MyMap = rcp(new map_type(NumGlobalEquations, NumMyEquations, 0, Comm));

  // Create the matrix
  A = rcp(new CrsMatrixType(MyMap,0));

  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1
  Teuchos::Array<Scalar> Values(10);
  Teuchos::Array<GO> Indices(10);
  size_t NumEntries=1;

  for (LO i = 0; i < NumMyEquations; i++) {
    GO GID = MyMap->getGlobalElement(i);
    // Super / subdiagonal
    if(GID == 0){
      Indices[0] = 1;
      NumEntries = 1;
    }
    else if (GID == NumGlobalEquations-1) {
      Indices[0] = NumGlobalEquations-2;
      NumEntries = 1;
    }
    else {
      Indices[0] = GID-1;
      Indices[1] = GID+1;
      NumEntries = 2;
    }

    // Wide stuff
    if(GID + NumMyEquations < NumGlobalEquations) {
      Indices[NumEntries]=GID + NumMyEquations;
      NumEntries++;
    }
    if(GID > Teuchos::as<GO>(NumMyEquations) ) { // Note: Unsigned integers are evil.
      Indices[NumEntries]=GID - NumMyEquations;
      NumEntries++;
    }

    // Double wide stuff.  Truck signage required
    if(GID + 2*NumMyEquations < NumGlobalEquations) {
      Indices[NumEntries]=GID + 2*NumMyEquations;
      NumEntries++;
    }
    if(GID > Teuchos::as<GO>(2*NumMyEquations) ) { // Note: Unsigned integers are evil.
      Indices[NumEntries]=GID - 2*NumMyEquations;
      NumEntries++;
    }

    Values[0] = -ONE;
    Values[1] = -ONE;
    Values[2] = -ONE;
    Values[3] = -ONE;
    Values[4] = -ONE;
    Values[5] = -ONE;
    Values[6] = -ONE;
    Values[7] = -ONE;
    A->insertGlobalValues(GID,Indices.view(0,NumEntries),Values.view(0,NumEntries));

    Indices[0] = GID;
    Values[0] = TWO;
    A->insertGlobalValues(GID,Indices.view(0,1), Values.view(0,1));
  }

  A->fillComplete();
}

// ===============================================================================
template<class CrsMatrixType>
void
build_test_matrix_with_row_overlap (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                    Teuchos::RCP<CrsMatrixType> & A)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef typename CrsMatrixType::scalar_type Scalar;
  typedef typename CrsMatrixType::local_ordinal_type LO;
  typedef typename CrsMatrixType::global_ordinal_type GO;
  typedef typename CrsMatrixType::node_type NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Tpetra::global_size_t GST;

  const Scalar ONE = Teuchos::ScalarTraits<Scalar>::one ();
  const Scalar TWO = ONE + ONE;
  const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
  const int NumProc = comm->getSize ();
  const int MyPID   = comm->getRank ();

  // Build the output matrix's overlapping row Map.
  const LO NumPrimaryEquations = 100;
  const LO NumLocalEquations = (NumProc > 1) ?
    (2*NumPrimaryEquations) : NumPrimaryEquations;
  GO start   = NumPrimaryEquations * MyPID;
  GO g_total = NumPrimaryEquations * NumProc;
  Teuchos::Array<GO> elementList (NumLocalEquations);
  for (LO i = 0; i < NumPrimaryEquations; ++i) {
    elementList[i] = start + i;
    if (NumProc > 1) {
      elementList[i+NumPrimaryEquations] =
        (start + NumPrimaryEquations + i) % g_total;
    }
  }
  RCP<const map_type> MyMap =
    rcp (new map_type (INVALID, elementList (), 0, comm));

  // Create the output matrix.
  A = rcp (new CrsMatrixType (MyMap, 0));

  // Fill the output matrix with entries.
  Teuchos::Array<Scalar> Values(1);
  Teuchos::Array<GO> Indices(1);
  for (LO i = 0; i < NumLocalEquations; ++i) {
    const GO GID = MyMap->getGlobalElement (i);
    Indices[0] = GID;
    if (i < NumPrimaryEquations) {
      Values[0] = TWO;
    } else {
      Values[0] = ONE;
    }
    A->insertGlobalValues (GID, Indices (), Values ());
  }

  // Create the domain/range Map.  Both of these must be one to one.
  RCP<const map_type> MyDRMap =
    rcp (new map_type (INVALID, NumPrimaryEquations, 0, comm));

  // Call fill complete on the output matrix.
  A->fillComplete (MyDRMap, MyDRMap);
}


// ===============================================================================
template<class CrsMatrixType>
void
build_test_prolongator (const Teuchos::RCP<const CrsMatrixType>& A,
                        Teuchos::RCP<CrsMatrixType>& P)
{
  typedef typename CrsMatrixType::local_ordinal_type LO;
  typedef typename CrsMatrixType::global_ordinal_type GO;
  typedef typename CrsMatrixType::scalar_type Scalar;
  typedef typename CrsMatrixType::node_type NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Tpetra::global_size_t GST;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
  RCP<const map_type> RowMap = A->getRowMap ();
  RCP<const map_type> DomainMap;

  // Create the matrix
  P = rcp(new CrsMatrixType(RowMap,0));

  // Make DomainMap
  Array<GO> gids;
  for (size_t i = 0; i < RowMap->getNodeNumElements (); ++i) {
    const GO gid = RowMap->getGlobalElement (i);
    if (gid % static_cast<GO> (3) == static_cast<GO> (0)) {
      gids.push_back (gid / 3);
    }
  }

  const GO indexBase = 0;
  DomainMap = rcp (new map_type (INVALID, gids (), indexBase, RowMap->getComm ()));

  Teuchos::Array<Scalar> Values(1);
  Teuchos::Array<GO> Indices(1);
  Values[0] = STS::one ();
  const GO minP = DomainMap->getMinGlobalIndex ();
  for (size_t i = 0; i < RowMap->getNodeNumElements (); ++i) {
    const GO GID = RowMap->getGlobalElement (i);
    Indices[0] = (static_cast<GO> (GID / 3.0) < minP) ? minP : static_cast<GO> (GID / 3.0);
    P->insertGlobalValues (GID, Indices (), Values ());
  }
  P->fillComplete (DomainMap, RowMap);

}

// ===============================================================================
template<class MapType>
void
build_test_map (const Teuchos::RCP<const MapType>& oldMap, Teuchos::RCP<MapType>& newMap)
{
  using Teuchos::rcp;
  typedef Tpetra::global_size_t GST;

  const int NumProc = oldMap->getComm()->getSize();
  const int MyPID   = oldMap->getComm()->getRank();

  if (NumProc < 3) {
    // Dump everything onto -proc 0
    GST num_global = oldMap->getGlobalNumElements();
    size_t num_local = MyPID==0 ? num_global : 0;
    newMap = rcp(new MapType(num_global,num_local,0,oldMap->getComm()));
  }
  else {
    // Split everything between procs 0 and 2 (leave proc 1 empty)
    GST num_global = oldMap->getGlobalNumElements();
    size_t num_local=0;
    if (MyPID == 0) {
      num_local = num_global/2;
    }
    else if (MyPID == 2) {
      num_local =  num_global - ((size_t)num_global/2);
    }
    newMap = rcp(new MapType(num_global,num_local,0,oldMap->getComm()));
  }
}

// ===============================================================================
template<class ImportType, class MapType>
void
build_remote_only_map (const Teuchos::RCP<const ImportType>& Import,
                       Teuchos::RCP<MapType>& newMap)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef typename MapType::local_ordinal_type LO;
  typedef typename MapType::global_ordinal_type GO;
  typedef Tpetra::global_size_t GST;

  if (Import.is_null ()) {
    return;
  }

  RCP<const MapType> targetMap = Import->getTargetMap ();
  const size_t NumRemotes = Import->getNumRemoteIDs ();
  ArrayView<const LO> oldRemoteLIDs = Import->getRemoteLIDs ();
  Array<GO> newRemoteGIDs (NumRemotes);

  for (size_t i=0; i < NumRemotes; ++i) {
    newRemoteGIDs[i] = targetMap->getGlobalElement (oldRemoteLIDs[i]);
  }

  const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
  newMap = rcp (new MapType (INVALID, newRemoteGIDs, targetMap->getIndexBase (),
                             targetMap->getComm ()));
}


TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( FusedImportExport, doImport, LO, GO, Scalar )
{
  typedef Tpetra::CrsGraph<LO, GO> Graph;
  typedef Tpetra::CrsMatrix<Scalar, LO, GO> CrsMatrixType;
  typedef Tpetra::Map<LO, GO> MapType;
  typedef Tpetra::Import<LO, GO> ImportType;
  typedef Tpetra::Export<LO, GO> ExportType;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MagType;
  typedef Tpetra::global_size_t GST;

  out << "Test importAndFillCompleteCrsMatrix and "
    "exportAndFillCompleteCrsMatrix (\"fused Import/Export + "
    "fillComplete\")" << endl;
  OSTab tab1 (out);
  RCP<const Comm<int> > Comm = getDefaultComm();

  RCP<CrsMatrixType> A, B, C;
  RCP<Graph> Bg;
  RCP<const MapType> Map1, Map2;
  RCP<MapType> Map3;

  RCP<ImportType> Import1;
  RCP<ExportType> Export1;
  const int MyPID = Comm->getRank();
  double diff;
  int total_err=0;
  MagType diff_tol = 1e4*Teuchos::ScalarTraits<Scalar>::eps();
  const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid();

  std::ostringstream err;
  int lclErr = 0;
  int gblErr = 0;

  try {
    build_test_matrix<CrsMatrixType> (Comm, A);
  } catch (std::exception& e) {
    err << "Proc " << MyPID << " threw an exception in build_test_matrix: "
        << e.what () << endl;
    lclErr = 1;
  }
  auto Ag = A->getCrsGraph();

  reduceAll<int, int> (*Comm, REDUCE_MAX, lclErr, outArg (gblErr));
  // The test fails if any (MPI) process had trouble.
  TEST_EQUALITY_CONST( gblErr, 0 );
  if (gblErr != 0) {
    Tpetra::Details::gathervPrint (out, err.str (), *Comm);
    out << "Above test failed; aborting further tests" << endl;
    return;
  }

  out << "Test #1: Tridiagonal Matrix; Migrate to Proc 0" << endl;
  try {
    OSTab tab2 (out);
    GST num_global = A->getRowMap()->getGlobalNumElements();

    // New map with all on Proc1
    if(MyPID==0) Map1 = rcp(new MapType(num_global,(size_t)num_global,0,Comm));
    else Map1 = rcp(new MapType(num_global,(size_t)0,0,Comm));

    // Execute fused import
    Import1 = rcp(new ImportType(A->getRowMap(),Map1));
    B = Tpetra::importAndFillCompleteCrsMatrix<CrsMatrixType>(A,*Import1);
    diff=test_with_matvec<CrsMatrixType>(*A,*B);
    if(diff > diff_tol) {
      if (MyPID == 0) {
        cerr << "FusedImport: Test #1 FAILED with norm diff = " << diff
             << "." << endl;
      }
      total_err--;
    }
    // Test the graph version
    Import1 = rcp(new ImportType(Ag->getRowMap(),Map1));
    Bg = Tpetra::importAndFillCompleteCrsGraph<Graph>(Ag,*Import1);
    if (!graphs_are_same(Bg, B->getCrsGraph())) {
      if (MyPID == 0) cerr << "FusedImport: CrsGraph test #1 FAILED." << endl;
      total_err--;
    }

    // Execute fused export
    Export1 = rcp(new ExportType(A->getRowMap(),Map1));
    B = Tpetra::exportAndFillCompleteCrsMatrix<CrsMatrixType>(A,*Export1);
    diff=test_with_matvec<CrsMatrixType>(*A,*B);
    if(diff > diff_tol) {
      if (MyPID == 0) {
        cerr << "FusedExport: Test #1 FAILED with norm diff = " << diff
             << "." << endl;
      }
      total_err--;
    }

    // Test the graph version
    Export1 = rcp(new ExportType(Ag->getRowMap(),Map1));
    Bg = Tpetra::exportAndFillCompleteCrsGraph<Graph>(Ag,*Export1);
    if (!graphs_are_same(Bg, B->getCrsGraph())) {
      if (MyPID == 0) cerr << "FusedExport: CrsGraph test #1 FAILED." << endl;
      total_err--;
    }

    Comm->barrier ();
  } catch (std::exception& e) {
    err << "Process " << MyPID << " threw an exception: " << e.what () << endl;
    lclErr = 1;
  }

  reduceAll<int, int> (*Comm, REDUCE_MAX, lclErr, outArg (gblErr));
  // The test fails if any (MPI) process had trouble.
  TEST_EQUALITY_CONST( gblErr, 0 );
  if (gblErr != 0) {
    Tpetra::Details::gathervPrint (out, err.str (), *Comm);
    out << "Above test failed; aborting further tests" << endl;
    return;
  }

  out << "Test #2: Tridiagonal Matrix; Locally Reversed Map" << endl;
  try {
    OSTab tab2 (out);
    size_t num_local = A->getRowMap()->getNodeNumElements();

    Teuchos::Array<GO> MyGIDs(num_local);
    for(size_t i=0; i<num_local; i++)
      MyGIDs[i] = A->getRowMap()->getGlobalElement(num_local-i-1);

    Map1 = rcp(new MapType(INVALID,MyGIDs(),0,Comm));

    // Execute fused import
    Import1 = rcp(new ImportType(A->getRowMap(),Map1));
    B = Tpetra::importAndFillCompleteCrsMatrix<CrsMatrixType>(A,*Import1);
    diff=test_with_matvec<CrsMatrixType>(*A,*B);
    if(diff > diff_tol) {
      if (MyPID == 0) {
        cerr << "FusedImport: Test #2 FAILED with norm diff = " << diff
             << "." << endl;
      }
      total_err--;
    }

    // Test the graph version
    Import1 = rcp(new ImportType(Ag->getRowMap(),Map1));
    Bg = Tpetra::importAndFillCompleteCrsGraph<Graph>(Ag,*Import1);
    if (!graphs_are_same(Bg, B->getCrsGraph())) {
      if (MyPID == 0) cerr << "FusedImport: CrsGraph test #2 FAILED." << endl;
      total_err--;
    }

    // Execute fused export
    Export1 = rcp(new ExportType(A->getRowMap(),Map1));
    B = Tpetra::exportAndFillCompleteCrsMatrix<CrsMatrixType>(A,*Export1);
    diff=test_with_matvec<CrsMatrixType>(*A,*B);
    if(diff > diff_tol) {
      if (MyPID == 0) {
        cerr << "FusedExport: Test #2 FAILED with norm diff = " << diff
             << "." << endl;
      }
      total_err--;
    }

    // Test the graph version
    Export1 = rcp(new ExportType(Ag->getRowMap(),Map1));
    Bg = Tpetra::exportAndFillCompleteCrsGraph<Graph>(Ag,*Export1);
    if (!graphs_are_same(Bg, B->getCrsGraph())) {
      if (MyPID == 0) cerr << "FusedExport: CrsGraph test #2 FAILED." << endl;
      total_err--;
    }

  } catch (std::exception& e) {
    err << "Process " << MyPID << " threw an exception: " << e.what ();
    lclErr = 1;
  }

  reduceAll<int, int> (*Comm, REDUCE_MAX, lclErr, outArg (gblErr));
  // The test fails if any (MPI) process had trouble.
  TEST_EQUALITY_CONST( gblErr, 0 );
  if (gblErr != 0) {
    Tpetra::Details::gathervPrint (out, err.str (), *Comm);
    out << "Above test failed; aborting further tests" << endl;
    return;
  }

  // mfh 28 Aug 2017: This test was skipped before; there was no code here.
  // I just made it print out that message.
  out << "Test #3: Tridiagonal Matrix; Globally Reversed Map (SKIPPED)" << endl;

  out << "Test #4: Tridiagonal Matrix; MMM style halo import" << endl;
  try {
    OSTab tab2 (out);
    // Assume we always own the diagonal
    size_t num_local = A->getNodeNumCols()-A->getNodeNumRows();
    Teuchos::Array<GO> MyGIDs(num_local);

    size_t idx=0;
    for (LO i_lcl = 0; static_cast<size_t> (i_lcl) < A->getNodeNumCols (); ++i_lcl) {
      const GO i_gbl = A->getColMap ()->getGlobalElement (i_lcl);
      if (A->getRowMap ()->getLocalElement (i_gbl) == Teuchos::OrdinalTraits<LO>::invalid ()) {
        MyGIDs[idx] = A->getColMap()->getGlobalElement (i_lcl);
        idx++;
      }
    }

    // New map & importer
    Map1=rcp(new MapType(INVALID,MyGIDs.view(0,idx),0,Comm));
    Import1 = rcp(new ImportType(A->getRowMap(),Map1));

    // Build unfused matrix to compare
    C = rcp(new CrsMatrixType(Map1,0));
    build_matrix_unfused_import<CrsMatrixType,ImportType>(*A,*Import1,C);

    // Execute fused import
    B = Tpetra::importAndFillCompleteCrsMatrix<CrsMatrixType>(A,*Import1);
    diff=test_with_matvec<CrsMatrixType>(*B,*C);
    if(diff > diff_tol) {
      if(MyPID==0) {
        cerr << "FusedImport: Test #4 FAILED with norm diff = " << diff
             << "." << endl;
      }
      total_err--;
    }

    // Test the graph version
    Import1 = rcp(new ImportType(Ag->getRowMap(),Map1));
    Bg = Tpetra::importAndFillCompleteCrsGraph<Graph>(Ag,*Import1);
    if (!graphs_are_same(Bg, B->getCrsGraph())) {
      if (MyPID == 0) cerr << "FusedImport: CrsGraph test #4 FAILED." << endl;
      total_err--;
    }

    // Execute fused export
    Export1 = rcp(new ExportType(A->getRowMap(),Map1));
    B = Tpetra::exportAndFillCompleteCrsMatrix<CrsMatrixType>(A,*Export1);
    diff=test_with_matvec<CrsMatrixType>(*B,*C);
    if(diff > diff_tol) {
      if(MyPID==0) {
        cerr << "FusedExport: Test #4 FAILED with norm diff = " << diff
             << "." << endl;
      }
      total_err--;
    }

    // Test the graph version
    Export1 = rcp(new ExportType(Ag->getRowMap(),Map1));
    Bg = Tpetra::exportAndFillCompleteCrsGraph<Graph>(Ag,*Export1);
    if (!graphs_are_same(Bg, B->getCrsGraph())) {
      if (MyPID == 0) cerr << "FusedExport: CrsGraph test #4 FAILED." << endl;
      total_err--;
    }

  } catch (std::exception& e) {
    err << "Process " << MyPID << " threw an exception: " << e.what () << endl;
    lclErr = 1;
  }

  reduceAll<int, int> (*Comm, REDUCE_MAX, lclErr, outArg (gblErr));
  // The test fails if any (MPI) process had trouble.
  TEST_EQUALITY_CONST( gblErr, 0 );
  if (gblErr != 0) {
    Tpetra::Details::gathervPrint (out, err.str (), *Comm);
    out << "Above test failed; aborting further tests" << endl;
    return;
  }

  out << "Test 5: Tridiagonal Matrix; Migrate to Proc 0, Replace Maps" << endl;
  try {
    OSTab tab2 (out);
    // New map with all on Proc 0 / 2
    build_test_map(A->getRowMap(),Map3);

    // Execute fused import
    Import1 = rcp(new ImportType(A->getRowMap(),Map3));
    B = Tpetra::importAndFillCompleteCrsMatrix<CrsMatrixType>(A,*Import1,Map3,Map3);
    diff=test_with_matvec<CrsMatrixType>(*A,*B);
    if(diff > diff_tol) {
      if(MyPID==0) {
        cerr << "FusedImport: Test #5 FAILED with norm diff = " << diff
             << "." << endl;
      }
      total_err--;
    }
    // Test the graph version
    Import1 = rcp(new ImportType(Ag->getRowMap(),Map3));
    Bg = Tpetra::importAndFillCompleteCrsGraph<Graph>(Ag,*Import1,Map3,Map3);
    if (!graphs_are_same(Bg, B->getCrsGraph())) {
      if (MyPID == 0) cerr << "FusedImport: CrsGraph test #5 FAILED." << endl;
      total_err--;
    }

    // Execute fused export
    Export1 = rcp(new ExportType(A->getRowMap(),Map3));
    B = Tpetra::exportAndFillCompleteCrsMatrix<CrsMatrixType>(A,*Export1,Map3,Map3);
    diff=test_with_matvec<CrsMatrixType>(*A,*B);
    if(diff > diff_tol) {
      if(MyPID==0) {
        cerr << "FusedExport: Test #5 FAILED with norm diff = " << diff
             << "." << endl;
      }
      total_err--;
    }

    // Test the graph version
    Export1 = rcp(new ExportType(Ag->getRowMap(),Map3));
    Bg = Tpetra::exportAndFillCompleteCrsGraph<Graph>(Ag,*Export1,Map3,Map3);
    if (!graphs_are_same(Bg, B->getCrsGraph())) {
      if (MyPID == 0) cerr << "FusedExport: CrsGraph test #5 FAILED." << endl;
      total_err--;
    }

  } catch (std::exception& e) {
    err << "Process " << MyPID << " threw an exception: " << e.what () << endl;
    lclErr = 1;
  }

  reduceAll<int, int> (*Comm, REDUCE_MAX, lclErr, outArg (gblErr));
  // The test fails if any (MPI) process had trouble.
  TEST_EQUALITY_CONST( gblErr, 0 );
  if (gblErr != 0) {
    Tpetra::Details::gathervPrint (out, err.str (), *Comm);
    out << "Above test failed; aborting further tests" << endl;
    return;
  }

  out << "Test 6: Tridiagonal Matrix; Migrate to Proc 0, Replace Comm" << endl;
  try {
    OSTab tab2 (out);
    // New map with all on Proc 0 / 2
    build_test_map(A->getRowMap(),Map3);

    // Parameters
    Teuchos::ParameterList params;
    params.set("Restrict Communicator",true);

    // Execute fused import constructor
    Import1 = rcp(new ImportType(A->getRowMap(),Map3));
    B = Tpetra::importAndFillCompleteCrsMatrix<CrsMatrixType>(A,*Import1,Map3,Map3,rcp(&params,false));
    diff=test_with_matvec_reduced_maps<CrsMatrixType,MapType>(*A,*B,*Map3);
    if(diff > diff_tol){
      if(MyPID==0) {
        cerr << "FusedImport: Test #6 FAILED with norm diff = " << diff
             << "." << endl;
      }
      total_err--;
    }

    // Test the graph version
    Import1 = rcp(new ImportType(Ag->getRowMap(),Map3));
    Bg = Tpetra::importAndFillCompleteCrsGraph<Graph>(Ag,*Import1,Map3,Map3,rcp(&params,false));
    if (Map3->getNodeNumElements() > 0) {
      if (!graphs_are_same(Bg, B->getCrsGraph())) {
        if (MyPID == 0) cerr << "FusedImport: CrsGraph test #6 FAILED." << endl;
        total_err--;
      }
    }

    // Execute fused export constructor
    Export1 = rcp(new ExportType(A->getRowMap(),Map3));
    B = Tpetra::exportAndFillCompleteCrsMatrix<CrsMatrixType>(A,*Export1,Map3,Map3,rcp(&params,false));
    diff=test_with_matvec_reduced_maps<CrsMatrixType,MapType>(*A,*B,*Map3);
    if(diff > diff_tol){
      if(MyPID==0) {
        cerr << "FusedExport: Test #6 FAILED with norm diff = " << diff
             << "." << endl;
      }
      total_err--;
    }

    // Test the graph version
    Export1 = rcp(new ExportType(Ag->getRowMap(),Map3));
    Bg = Tpetra::exportAndFillCompleteCrsGraph<Graph>(Ag,*Export1,Map3,Map3,rcp(&params,false));
    if (Map3->getNodeNumElements() > 0) {
      if (!graphs_are_same(Bg, B->getCrsGraph())) {
        if (MyPID == 0) cerr << "FusedExport: CrsGraph test #6 FAILED." << endl;
        total_err--;
      }
    }

  } catch (std::exception& e) {
    err << "Process " << MyPID << " threw an exception: " << e.what () << endl;
    lclErr = 1;
  }

  reduceAll<int, int> (*Comm, REDUCE_MAX, lclErr, outArg (gblErr));
  // The test fails if any (MPI) process had trouble.
  TEST_EQUALITY_CONST( gblErr, 0 );
  if (gblErr != 0) {
    Tpetra::Details::gathervPrint (out, err.str (), *Comm);
    out << "Above test failed; aborting further tests" << endl;
    return;
  }

  out << "Test 7: Tridiagonal Matrix; Migrate to Proc 0, Reverse Mode" << endl;
  try {
    OSTab tab2 (out);
    GST num_global = A->getRowMap()->getGlobalNumElements();

    // New map with all on Proc1
    if(MyPID==0) Map1 = rcp(new MapType(num_global,(size_t)num_global,0,Comm));
    else Map1 = rcp(new MapType(num_global,(size_t)0,0,Comm));

    // Parameters
    Teuchos::ParameterList params;
    params.set("Reverse Mode",true);

    // Execute fused import constructor
    Import1 = rcp(new ImportType(Map1,A->getRowMap()));
    B = Tpetra::importAndFillCompleteCrsMatrix<CrsMatrixType>(A,*Import1,Map1,Map1,rcp(&params,false));
    diff=test_with_matvec<CrsMatrixType>(*A,*B);
    if(diff > diff_tol){
      if(MyPID==0) {
        cerr << "FusedImport: Test #7 FAILED with norm diff = " << diff
             << "." << endl;
      }
      total_err--;
    }

    // Test the graph version
    Import1 = rcp(new ImportType(Map1,Ag->getRowMap()));
    Bg = Tpetra::importAndFillCompleteCrsGraph<Graph>(Ag,*Import1,Map1,Map1,rcp(&params,false));
    if (!graphs_are_same(Bg, B->getCrsGraph())) {
      if (MyPID == 0) cerr << "FusedImport: CrsGraph test #7 FAILED." << endl;
      total_err--;
    }

    // Execute fused export constructor
    Export1 = rcp(new ExportType(Map1,A->getRowMap()));
    B = Tpetra::exportAndFillCompleteCrsMatrix<CrsMatrixType>(A,*Export1,Map1,Map1,rcp(&params,false));
    diff=test_with_matvec<CrsMatrixType>(*A,*B);
    if(diff > diff_tol){
      if(MyPID==0) {
        cerr << "FusedExport: Test #7 FAILED with norm diff = " << diff
             << "." << endl;
      }
      total_err--;
    }

    // Test the graph version
    Export1 = rcp(new ExportType(Map1,Ag->getRowMap()));
    Bg = Tpetra::exportAndFillCompleteCrsGraph<Graph>(Ag,*Export1,Map1,Map1,rcp(&params,false));
    if (!graphs_are_same(Bg, B->getCrsGraph())) {
      if (MyPID == 0) cerr << "FusedExport: CrsGraph test #7 FAILED." << endl;
      total_err--;
    }

  } catch (std::exception& e) {
    err << "Process " << MyPID << " threw an exception: " << e.what () << endl;
    lclErr = 1;
  }

  reduceAll<int, int> (*Comm, REDUCE_MAX, lclErr, outArg (gblErr));
  // The test fails if any (MPI) process had trouble.
  TEST_EQUALITY_CONST( gblErr, 0 );
  if (gblErr != 0) {
    Tpetra::Details::gathervPrint (out, err.str (), *Comm);
    out << "Above test failed; aborting further tests" << endl;
    return;
  }

  out << "Test #8: Diagonal matrix w/ overlapping entries" << endl;
  try {
    OSTab tab2 (out);
    build_test_matrix_with_row_overlap<CrsMatrixType>(Comm,A);
    Ag = A->getCrsGraph();

    Teuchos::ArrayRCP<const size_t> rowptr;
    Teuchos::ArrayRCP<const LO> colind;
    Teuchos::ArrayRCP<const Scalar> vals;

    Map1 = A->getRangeMap();

    // Execute fused import constructor (reverse)
    Teuchos::ParameterList params;
    params.set("Reverse Mode",true);
    Import1 = rcp(new ImportType(Map1,A->getRowMap()));
    B = Tpetra::importAndFillCompleteCrsMatrix<CrsMatrixType>(A,*Import1,Map1,Map1,rcp(&params,false));
    diff=test_with_matvec<CrsMatrixType>(*B,*A);
    if(diff > diff_tol){
      if(MyPID==0) {
        cerr << "FusedImport: Test #8 FAILED with norm diff = " << diff
             << "." << endl;
      }
      total_err--;
    }

    // Test the graph version
    Import1 = rcp(new ImportType(Map1,Ag->getRowMap()));
    Bg = Tpetra::importAndFillCompleteCrsGraph<Graph>(Ag,*Import1,Map1,Map1,rcp(&params,false));
    if (!graphs_are_same(Bg, B->getCrsGraph())) {
      if (MyPID == 0) cerr << "FusedImport: CrsGraph test #8 FAILED." << endl;
      total_err--;
    }

    // Execute fused export constructor
    Export1 = rcp(new ExportType(A->getRowMap(),Map1));
    B = Tpetra::exportAndFillCompleteCrsMatrix<CrsMatrixType>(A,*Export1,Map1,Map1);
    diff=test_with_matvec<CrsMatrixType>(*B,*A);
    if(diff > diff_tol){
      if(MyPID==0) {
        cout << "FusedExport: Test #8 FAILED with norm diff = " << diff
             << "." << endl;
      }
      total_err--;
    }

    // Test the graph version
    Export1 = rcp(new ExportType(Ag->getRowMap(),Map1));
    Bg = Tpetra::exportAndFillCompleteCrsGraph<Graph>(Ag,*Export1,Map1,Map1);
    if (!graphs_are_same(Bg, B->getCrsGraph())) {
      if (MyPID == 0) cerr << "FusedExport: CrsGraph test #8 FAILED." << endl;
      total_err--;
    }

  } catch (std::exception& e) {
    err << "Process " << MyPID << " threw an exception: " << e.what () << endl;
    lclErr = 1;
  }

  reduceAll<int, int> (*Comm, REDUCE_MAX, lclErr, outArg (gblErr));
  // The test fails if any (MPI) process had trouble.
  TEST_EQUALITY_CONST( gblErr, 0 );
  if (gblErr != 0) {
    Tpetra::Details::gathervPrint (out, err.str (), *Comm);
    out << "Above test failed; aborting further tests" << endl;
    return;
  }

  TEST_EQUALITY(total_err,0);
}



TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( ReverseImportExport, doImport, LO, GO, Scalar )  {
  // NOTE: This test will fail.

  RCP<const Comm<int> > Comm = getDefaultComm();
  typedef Tpetra::Map<LO, GO> MapType;
  typedef Tpetra::Import<LO, GO> ImportType;
  typedef Tpetra::Export<LO, GO> ExportType;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MagType;
  typedef Tpetra::Vector<Scalar, LO, GO> VectorType;
  typedef Tpetra::CrsMatrix<Scalar, LO, GO> CrsMatrixType;
  typedef Tpetra::global_size_t GST;

  RCP<CrsMatrixType> A;

  RCP<VectorType> SourceVector, TargetVector, TestVector;
  RCP<const MapType> MapSource, MapTarget;
  RCP<ImportType> Import1;
  RCP<ExportType> Export1;
  int MyPID = Comm->getRank();
  int total_err=0;
  MagType diff_tol = 1e4*Teuchos::ScalarTraits<Scalar>::eps();
  int test_err=0;

  // Get the source map from a sample matrix
  build_test_matrix<CrsMatrixType>(Comm,A);
  MapSource=A->getRowMap();
  GST num_global = A->getRowMap()->getGlobalNumElements();

  // Target Map - all on one proc
  size_t num_local = MyPID==0 ? num_global : 0;
  MapTarget = rcp(new MapType(num_global,num_local,0,Comm));

  // Vectors
  SourceVector = rcp(new VectorType(MapSource));
  TargetVector = rcp(new VectorType(MapTarget));
  TestVector   = rcp(new VectorType(MapTarget));

  // Importer / Exporter
  Import1 = rcp(new ImportType(MapSource,MapTarget));
  Export1 = rcp(new ExportType(MapSource,MapTarget));
  Teuchos::Array< typename Teuchos::ScalarTraits< Scalar >::magnitudeType > norms(1);


  Teuchos::ScalarTraits< Scalar >::seedrandom(24601);
  SourceVector->randomize();
  TestVector->doImport(*SourceVector,*Import1,Tpetra::INSERT);

  /////////////////////////////////////////////////////////
  // Test #1: Use Exporter to create a reverse import
  /////////////////////////////////////////////////////////
  {
    TargetVector->putScalar(0.0);
    RCP<ImportType> Import2 = rcp(new ImportType(*Export1));

    TargetVector->doExport(*SourceVector,*Import2,Tpetra::ADD);

    TargetVector->update(-1.0,*TestVector,1.0);
    TargetVector->norm2(norms);
    if(norms[0] > diff_tol) {
      if(MyPID==0) cout<<"ReverseImport: Test #1 FAILED with norm diff = "<<norms[0]<<"."<<endl;
      total_err--;
    }
  }

  /////////////////////////////////////////////////////////
  // Test #2: Use Importer to create a reverse exporter
  /////////////////////////////////////////////////////////
  {
    TargetVector->putScalar(0.0);
    RCP<ExportType> Export2 = rcp(new ExportType(*Import1));

    TargetVector->doExport(*SourceVector,*Export2,Tpetra::ADD);

    TargetVector->update(-1.0,*TestVector,1.0);
    TargetVector->norm2(norms);

    if(norms[0] > diff_tol) {
      if(MyPID==0) cout<<"ReverseExport: Test #2 FAILED with norm diff = "<<norms[0]<<"."<<endl;
      total_err--;
    }
  }


  TEST_EQUALITY(test_err,0);
}



TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Import_Util, GetPids, LO, GO )  {
  // Unit Test the functionality in Tpetra_Import_Util
  RCP<const Comm<int> > Comm = getDefaultComm();
  typedef Tpetra::Map<LO, GO> MapType;
  typedef Tpetra::Import<LO, GO> ImportType;
  typedef Tpetra::Export<LO, GO> ExportType;
  typedef Tpetra::Vector<int, LO, GO> IntVectorType;
  typedef Tpetra::CrsMatrix<double, LO, GO> CrsMatrixType;
  typedef Tpetra::global_size_t GST;

  RCP<CrsMatrixType> A;

  RCP<IntVectorType> SourceVector, TargetVector, TestVector;
  RCP<const MapType> MapSource, MapTarget;
  RCP<ImportType> Import1;
  RCP<ExportType> Export1;
  int MyPID = Comm->getRank();
  int total_err=0;
  int test_err=0;

  // Get the source map from a sample matrix
  build_test_matrix<CrsMatrixType>(Comm,A);
  MapSource=A->getRowMap();
  const GST num_global = A->getRowMap()->getGlobalNumElements();

  // Target Map - all on one proc
  const size_t num_local = MyPID==0 ? num_global : 0;
  MapTarget = rcp (new MapType (num_global,num_local,0,Comm));

  // Vectors
  SourceVector = rcp(new IntVectorType(MapSource));
  TargetVector = rcp(new IntVectorType(MapTarget));
  TestVector   = rcp(new IntVectorType(MapTarget));

  // Importer / Exporter
  Import1 = rcp(new ImportType(MapSource,MapTarget));

  // Generate PID vector via explicit import
  SourceVector->putScalar(MyPID);
  TargetVector->putScalar(0);
  TargetVector->doImport(*SourceVector,*Import1,Tpetra::ADD);


  /////////////////////////////////////////////////////////
  // Test #1: getPids
  /////////////////////////////////////////////////////////
  {
    test_err=0;
    // Generate PID vector via getPids
    Teuchos::Array<int> pids;
    Tpetra::Import_Util::getPids<LO, GO, Node> (*Import1, pids, false);

    // Compare
    Teuchos::ArrayRCP<const int>  TargetView = TargetVector->get1dView();
    for(size_t i=0; i < TargetVector->getLocalLength(); ++i) {
      test_err += (pids[i] - TargetView[i] > 0) ?  1 : ((pids[i] - TargetView[i] < 0) ?  1 : 0);
    }
    total_err += test_err;
  }

  /////////////////////////////////////////////////////////
  // Test #2: getPidGidPairs
  /////////////////////////////////////////////////////////
  {
    test_err=0;
    // Generate PID vector via getPidGidPairs
    Teuchos::Array<std::pair<int, GO> > pgids;
    Tpetra::Import_Util::getPidGidPairs<LO, GO, Node> (*Import1, pgids, false);

    // Compare
    Teuchos::ArrayRCP<const int> TargetView = TargetVector->get1dView();
    for (size_t i=0; i < TargetVector->getLocalLength(); ++i) {
      test_err += (pgids[i].first - TargetView[i] > 0) ?  1 : ((pgids[i].first - TargetView[i] < 0) ?  1 : 0);
    }
    total_err += test_err;
  }

  /////////////////////////////////////////////////////////
  // Test #3: getRemotePids
  /////////////////////////////////////////////////////////
  {
    test_err=0;
    Teuchos::Array<int> RemotePIDs;
    Tpetra::Import_Util::getRemotePIDs<LO, GO, Node> (*Import1,RemotePIDs);

    // We can't easily test this, so let's at least make sure it doesn't crash.
  }

  TEST_EQUALITY(total_err,0);
}


TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Import_Util, PackAndPrepareWithOwningPIDs, LO, GO )  {
  // FIXME (mfh 28 Jan 2015) We've since changed the packing format,
  // so this test is broken for now.
  return;

#if 0
  // Unit Test the functionality in Tpetra_Import_Util
  RCP<const Comm<int> > Comm = getDefaultComm();
  typedef Tpetra::Import<LO, GO> ImportType;
  typedef Tpetra::CrsMatrix<double, LO, GO> CrsMatrixType;
  typedef typename CrsMatrixType::device_type device_type;
  using Teuchos::av_reinterpret_cast;

  RCP<CrsMatrixType> A;
  int total_err=0;
  int test_err=0;
  Teuchos::Array<char> exports1;
  Kokkos::DualView<char*, device_type> exports2;
  Teuchos::Array<size_t> numPackets1, numPackets2;
  size_t constantNumPackets1=0, constantNumPackets2=0;

  // Build sample matrix
  build_test_matrix<CrsMatrixType>(Comm,A);

  // Get Importer
  RCP<const ImportType> Importer = A->getCrsGraph()->getImporter();
  if(Importer == Teuchos::null)  {
    TEST_EQUALITY(0,0); // short circuit
  }
  else {
    /////////////////////////////////////////////////////////
    // Test #1: P&PWOPIDs
    /////////////////////////////////////////////////////////
    // Call Standard Pack & Prepare
    test_err=0;
    constantNumPackets1=0;
    numPackets1.resize(Importer->getExportLIDs().size());
    A->packAndPrepare(*A,Importer->getExportLIDs(),exports1,numPackets1(),constantNumPackets1,Importer->getDistributor());

    // Call the P&PWOPIDs
    Teuchos::Array<int> pids;
    Tpetra::Import_Util::getPids<LO, GO, Node>(*Importer,pids,false);
    constantNumPackets2=0;
    numPackets2.resize(Importer->getExportLIDs().size());
    Tpetra::Details::packCrsMatrixWithOwningPIDs<double, LO, GO, Node>(
        *A, exports2, InumPackets2(), Importer->getExportLIDs(), pids(),
        constantNumPackets2, Importer->getExportLIDs());

    // This test reads exports2 on the host, so sync there.
    exports2.template sync<Kokkos::HostSpace> ();
    Teuchos::ArrayView<char> exports2_av = getArrayViewFromDualView (exports2);

    // Loop through the parts that should be the same
    const size_t numExportLIDs = Importer->getExportLIDs().size();
    Teuchos::ArrayView<const LO> exportLIDs = Importer->getExportLIDs();

    const int sizeOfPacket1 = sizeof(double) + sizeof(GO);
    const int sizeOfPacket2 = sizeof(double) + sizeof(GO) + sizeof(int);

    size_t offset1=0,offset2=0;
    for (size_t i = 0; i < numExportLIDs; i++) {
      ArrayView<const LO> lidsView;
      ArrayView<const double>  valsView;
      A->getLocalRowView(exportLIDs[i], lidsView, valsView);
      const size_t curNumEntries = lidsView.size();

      ArrayView<char> gidsChar1 = exports1(offset1, curNumEntries*sizeof(GO));
      ArrayView<char> valsChar1 = exports1(offset1+curNumEntries*sizeof(GO), curNumEntries*sizeof(double));
      ArrayView<GO> gids1  = av_reinterpret_cast<GO>(gidsChar1);
      ArrayView<double>  vals1  = av_reinterpret_cast<double>(valsChar1);

      ArrayView<char> gidsChar2 = exports2_av(offset2, curNumEntries*sizeof(GO));
      //      ArrayView<char> pidsChar2 = exports2_av(offset2+curNumEntries*sizeof(GO), curNumEntries*sizeof(int));
      ArrayView<char> valsChar2 = exports2_av(offset2+curNumEntries*(sizeof(GO)+sizeof(int)), curNumEntries*sizeof(double));
      ArrayView<GO> gids2  = av_reinterpret_cast<GO>(gidsChar2);
      ArrayView<double>  vals2  = av_reinterpret_cast<double>(valsChar2);

      for (size_t k = 0; k < curNumEntries; ++k) {
        if(gids1[k] != gids2[k] || vals1[k] != vals2[k]) test_err++;
      }
      offset1 += sizeOfPacket1 * curNumEntries;
      offset2 += sizeOfPacket2 * curNumEntries;
    }
    total_err+=test_err;
  }

  TEST_EQUALITY(total_err,0);
#endif // 0
}



// Unit Test the functionality in Tpetra_Import_Util
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Import_Util, UnpackAndCombineWithOwningPIDs, LO, GO, Scalar)  {
  using Teuchos::av_reinterpret_cast;
  typedef Tpetra::Map<LO, GO> MapType;
  typedef Tpetra::Import<LO, GO> ImportType;
  typedef Tpetra::CrsMatrix<Scalar, LO, GO> CrsMatrixType;
  typedef typename Tpetra::CrsMatrix<Scalar, LO, GO>::packet_type PacketType;
  typedef typename MapType::device_type device_type;
  typedef Tpetra::global_size_t GST;
  typedef typename CrsMatrixType::impl_scalar_type IST;

  RCP<const Comm<int> > Comm = getDefaultComm();
  RCP<CrsMatrixType> A,B;
  int total_err=0;
  int test_err=0;
  int MyPID = Comm->getRank();
  RCP<const MapType> MapSource, MapTarget;

  // mfh 01 Aug 2017: Deal with fix for #1088, by not using
  // Kokkos::CudaUVMSpace for communication buffers.
#ifdef KOKKOS_ENABLE_CUDA
  typedef typename std::conditional<
  std::is_same<typename device_type::execution_space, Kokkos::Cuda>::value,
    Kokkos::CudaSpace,
    typename device_type::memory_space>::type buffer_memory_space;
#else
  typedef typename device_type::memory_space buffer_memory_space;
#endif // KOKKOS_ENABLE_CUDA
  typedef typename device_type::execution_space buffer_execution_space;
  typedef Kokkos::Device<buffer_execution_space, buffer_memory_space> buffer_device_type;

  Kokkos::DualView<char*, buffer_device_type> exports;
  Teuchos::Array<char> imports;
  Teuchos::Array<size_t> numImportPackets, numExportPackets;
  size_t constantNumPackets=0;

  // Build sample matrix & sourceMap
  build_test_matrix<CrsMatrixType>(Comm,A);
  GST num_global = A->getRowMap()->getGlobalNumElements();
  MapSource = A->getRowMap();

  // Target Map - all on one proc
  size_t num_local = MyPID==0 ? num_global : 0;
  MapTarget = rcp(new MapType(num_global,num_local,0,Comm));

  // Build Importer
  RCP<ImportType> Importer = rcp (new ImportType (MapSource, MapTarget));
  if (Importer != Teuchos::null)  {
    /////////////////////////////////////////////////////////
    // Test #1: Count test
    /////////////////////////////////////////////////////////
    // Do the traditional import
    test_err=0;
    B = rcp(new CrsMatrixType(MapTarget,0));
    B->doImport(*A, *Importer, Tpetra::INSERT);
    B->fillComplete(A->getDomainMap(),A->getRangeMap());
    size_t nnz1=B->getNodeNumEntries();

    // Call the P&PWOPIDs
    Teuchos::Array<int> SourcePids(A->getColMap()->getNodeNumElements());
    Tpetra::Distributor &distor = Importer->getDistributor();
    if(A->getGraph()->getImporter()==Teuchos::null) SourcePids.assign(SourcePids.size(),MyPID);
    else Tpetra::Import_Util::getPids<LO, GO, Node>(*A->getGraph()->getImporter(),SourcePids,false);
    numExportPackets.resize(Importer->getExportLIDs().size());
    numImportPackets.resize(Importer->getRemoteLIDs().size());

    Tpetra::Details::packCrsMatrixWithOwningPIDs<Scalar, LO, GO, Node>(
        *A, exports, numExportPackets(), Importer->getExportLIDs(),
        SourcePids(), constantNumPackets, distor);

    // This test reads exports on the host, so sync there.
    exports.template sync<Kokkos::HostSpace> ();
    Teuchos::ArrayView<char> exports_av = getArrayViewFromDualView (exports);

    // Do the moral equivalent of doTransfer
    distor.doPostsAndWaits<size_t>(numExportPackets().getConst(), 1,numImportPackets());
    size_t totalImportPackets = 0;
    for(size_t i = 0; i < (size_t)numImportPackets.size(); i++) {
      totalImportPackets += numImportPackets[i];
    }
    imports.resize(totalImportPackets);
    distor.doPostsAndWaits<PacketType>(exports_av,numExportPackets(),imports(),numImportPackets());

    // Run the count... which should get the same NNZ as the traditional import
    using Tpetra::Details::unpackAndCombineWithOwningPIDsCount;
    size_t nnz2 =
      unpackAndCombineWithOwningPIDsCount<Scalar, LO, GO, Node> (*A, Importer->getRemoteLIDs (),
                                                                 imports (), numImportPackets (),
                                                                 constantNumPackets, distor,
                                                                 Tpetra::INSERT,
                                                                 Importer->getNumSameIDs (),
                                                                 Importer->getPermuteToLIDs (),
                                                                 Importer->getPermuteFromLIDs ());
    if(nnz1!=nnz2) test_err++;
    total_err+=test_err;

    /////////////////////////////////////////////////////////
    // Test #2: Actual combine test
    /////////////////////////////////////////////////////////
    Teuchos::Array<size_t>  rowptr (MapTarget->getNodeNumElements () + 1);
    Teuchos::Array<GO>      colind (nnz2);
    Teuchos::Array<Scalar>  vals (nnz2);
    Teuchos::Array<int>     TargetPids;

    using Tpetra::Details::unpackAndCombineIntoCrsArrays;
    unpackAndCombineIntoCrsArrays<Scalar, LO, GO, Node> (
      *A,
      Importer->getRemoteLIDs (),
      imports (),
      numImportPackets (),
      constantNumPackets,
      distor,
      Tpetra::INSERT,
      Importer->getNumSameIDs (),
      Importer->getPermuteToLIDs (),
      Importer->getPermuteFromLIDs (),
      MapTarget->getNodeNumElements (),
      nnz2,
      MyPID,
      rowptr (),
      colind (),
      Teuchos::av_reinterpret_cast<IST> (vals ()),
      SourcePids (),
      TargetPids);
    // Do the comparison
    Teuchos::ArrayRCP<const size_t>  Browptr;
    Teuchos::ArrayRCP<const LO>      Bcolind;
    Teuchos::ArrayRCP<const Scalar>  Bvals;
    B->getAllValues (Browptr, Bcolind, Bvals);

    // Check the rowptrs
    if(Browptr.size()!= rowptr.size()) test_err++;
    if(!test_err) {
      for(size_t i=0; i < as<size_t>(rowptr.size()); i++) {
        if(Browptr[i]!=rowptr[i]) {
          test_err++;
          break;
        }
      }
    }
    // Check the indices / values... but sort first
    if(Bcolind.size()!=colind.size()) {test_err++; std::cout<<"--colind mismatch"<<std::endl;}
    if(Bvals.size()  !=vals.size())   {test_err++; std::cout<<"--vals mismatch"<<std::endl;}
    if(!test_err) {

      // Reindex colind to local indices
      const MapType & Bmap = *B->getColMap();
      for(size_t i=0; i<as<size_t>(colind.size()); i++){
        colind[i]=Bmap.getLocalElement(colind[i]);
      }

      // Sort the GIDs
      Tpetra::Import_Util::sortCrsEntries<Scalar, GO> (rowptr (), colind (), vals ());

      // Compare the gids
      for (size_t i=0; i < static_cast<size_t> (rowptr.size()-1); ++i) {
        for (size_t j=rowptr[i]; j < rowptr[i+1]; ++j) {
          if (colind[j] != Bcolind[j] || vals[j] != Bvals[j]) {
            test_err++;
          }
        }
      }
    }

    total_err+=test_err;
  }

  TEST_EQUALITY(total_err,0);
}


TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Import_Util,LowCommunicationMakeColMapAndReindex, LO, GO)  {
  // Test the colmap...
  RCP<const Comm<int> > Comm = getDefaultComm();
  typedef double Scalar;
  typedef Tpetra::Map<LO,GO> MapType;
  typedef Tpetra::Import<LO,GO> ImportType;
  typedef Tpetra::CrsMatrix<Scalar,LO,GO> CrsMatrixType;
  typedef Teuchos_Ordinal rsize_t;
  using Teuchos::av_reinterpret_cast;

  RCP<CrsMatrixType> A;
  RCP<const MapType> Acolmap, Adomainmap;
  RCP<const MapType> Bcolmap;
  int total_err=0;
  int test_err=0;

  // Build sample matrix
  build_test_matrix_wideband<CrsMatrixType>(Comm,A);

  // Get the matrix pointers / map
  ArrayRCP<const size_t> rowptr;
  ArrayRCP<const LO> colind;
  ArrayRCP<const Scalar> values;
  A->getAllValues(rowptr,colind,values);
  Acolmap = A->getColMap();
  Adomainmap = A->getDomainMap();

  // Get owning PID information
  size_t numMyCols = A->getColMap()->getNodeNumElements();
  RCP<const ImportType> Importer = A->getGraph()->getImporter();
  Teuchos::Array<int> AcolmapPIDs(numMyCols,-1);
  if(Importer!=Teuchos::null)
    Tpetra::Import_Util::getPids<LO, GO, Node>(*Importer,AcolmapPIDs,true);

  // Build a "gid" version of colind & colind-sized pid list
  Array<GO> colind_GID(colind.size());
  Array<int> colind_PID(colind.size());
  for(rsize_t i=0; i<colind.size(); i++) {
    colind_GID[i] = Acolmap->getGlobalElement(colind[i]);
    colind_PID[i] = AcolmapPIDs[colind[i]];
  }


  {
    /////////////////////////////////////////////////////////
    // Test #1: Pre-sorted colinds
    /////////////////////////////////////////////////////////
    Teuchos::Array<int> BcolmapPIDs;
    Teuchos::Array<LO> Bcolind_LID(colind.size());
    Tpetra::Import_Util::lowCommunicationMakeColMapAndReindex<LO, GO, Node>(rowptr(),Bcolind_LID(),colind_GID(),Adomainmap,colind_PID(),BcolmapPIDs,Bcolmap);

    // Since this was sorted to begin with, the outgoing maps should be
    // in an identical order.  So let's check.
    if(Acolmap->getNodeNumElements()!=Bcolmap->getNodeNumElements()) test_err++;
    else {
      for(size_t i=0; i<Acolmap->getNodeNumElements(); i++)
        if(Acolmap->getGlobalElement(i)!=Bcolmap->getGlobalElement(i)) test_err++;
    }

    // Now test the column indices
    if(colind.size()!=Bcolind_LID.size()) test_err++;
    else {
      for(rsize_t i=0; i<colind.size(); i++)
        if(colind[i] != Bcolind_LID[i]) test_err++;
    }

    total_err+=test_err;
  }

  TEST_EQUALITY(total_err,0);
}


TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Import, AdvancedConstructors, LO, GO )  {
  // Test the remotePIDs Tpetra::Import constructor
  RCP<const Comm<int> > Comm = getDefaultComm();
  typedef Tpetra::Map<LO, GO> MapType;
  typedef Tpetra::Import<LO, GO> ImportType;
  typedef Tpetra::CrsMatrix<double, LO, GO> CrsMatrixType;
  RCP<CrsMatrixType> A,B;

  RCP<const ImportType> Import1, Import2;
  RCP<MapType> Map0;
  int total_err=0;
  int test_err=0;

  // Build the sample matrix
  build_test_matrix<CrsMatrixType>(Comm,A);

  // Grab its importer
  Import1 = A->getGraph()->getImporter();

  // Only test in parallel
  if(Comm->getSize()==1) { TEST_EQUALITY(0,0); return;}


  // Rebalanced map (based on the *DomainMap* of A)
  build_test_map(A->getDomainMap(),Map0);

  /////////////////////////////////////////////////////////
  // Test #1: Constructor w/ remotePIDs test
  /////////////////////////////////////////////////////////
  {
    test_err=0;

    // Generate PID vector via getRemotePIDs
    Teuchos::Array<int> pids;
    Tpetra::Import_Util::getRemotePIDs<LO, GO, Node>(*Import1,pids);

    // Build a new (identical) importer via the other constructor
    Import2 = rcp(new ImportType(Import1->getSourceMap(),Import1->getTargetMap(),pids));

    // Compare
    if(Import1->getNumSameIDs() != Import2->getNumSameIDs()) test_err++;
    if(Import1->getNumPermuteIDs() != Import2->getNumPermuteIDs())
      test_err++;
    else {
      for(size_t i=0; i<Import1->getNumPermuteIDs(); i++) {
        test_err += (Import1->getPermuteFromLIDs()[i]!=Import2->getPermuteFromLIDs()[i]);
        test_err += (Import1->getPermuteToLIDs()[i]!=Import2->getPermuteToLIDs()[i]);
      }
    }
    if(Import1->getNumRemoteIDs() != Import2->getNumRemoteIDs()) {
      test_err++;
    }
    else {
      auto remoteLids1 = Import1->getRemoteLIDs ();
      auto remoteLids2 = Import2->getRemoteLIDs ();

      const size_t numRemoteIds = Import1->getNumRemoteIDs ();
      for(size_t i=0; i<numRemoteIds; i++) {
        test_err += (remoteLids1[i] != remoteLids2[i]);
      }
    }
    if(Import1->getNumExportIDs() != Import2->getNumExportIDs()) {
      test_err++;
    }
    else {
      auto exportLids1 = Import1->getExportLIDs ();
      auto exportLids2 = Import2->getExportLIDs ();
      auto exportPids1 = Import1->getExportPIDs ();
      auto exportPids2 = Import2->getExportPIDs ();

      const size_t numExportIds = Import1->getNumExportIDs ();
      for(size_t i=0; i<numExportIds; i++) {
        test_err += (exportLids1[i] != exportLids2[i]);
        test_err += (exportPids1[i] != exportPids2[i]);
      }
    }
    total_err += test_err;
  }

  /////////////////////////////////////////////////////////
  // Test #2: MueLu-style Transpose & Rebalance "R"
  /////////////////////////////////////////////////////////
  {
    // Create Transpose
    Tpetra::RowMatrixTransposer<double, LO, GO, Node> transposer(A);
    B = transposer.createTranspose();

    // Build Importer
    Import2 = rcp(new ImportType(B->getRowMap(),Map0));

    // Build the imported matrix
    Tpetra::importAndFillCompleteCrsMatrix<CrsMatrixType>(B,*Import2,Teuchos::null,Map0,Teuchos::null);

  }
  //


  TEST_EQUALITY(total_err,0);
}



 TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( FusedImportExport, MueLuStyle, LO, GO, Scalar )  {
  // Test a muelu-style SA build and rebalance.  Kind of like Gangnam Style, but with more cows.
  RCP<const Comm<int> > Comm = getDefaultComm();
  typedef Tpetra::Map<LO, GO> MapType;
  typedef Tpetra::Import<LO, GO> ImportType;
  typedef Tpetra::CrsMatrix<double, LO, GO> CrsMatrixType;
  RCP<CrsMatrixType> A,Ptent,P,R,AP,RAP,rebalancedP;
  RCP<MapType> Map0;
  RCP<const ImportType> Import0,ImportTemp;

  const int myRank = Comm->getRank ();
  int total_err=0;
  int test_err=0;

  std::ostringstream err;
  int lclErr = 0;
  int gblErr = 0;

  // Build sample matrix
  try {
    build_test_matrix_wideband<CrsMatrixType>(Comm,A);
  } catch (std::exception& e) {
    err << "Proc " << myRank << ": " << e.what ();
    lclErr = 1;
  }

  reduceAll<int, int> (*Comm, REDUCE_MAX, lclErr, outArg (gblErr));
  if (gblErr != 0) {
    for (int r = 0; r < Comm->getSize (); ++r) {
      if (r == myRank) {
        cerr << err.str () << endl;
      }
      Comm->barrier ();
      Comm->barrier ();
      Comm->barrier ();
    }
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Test failed!");
  }

  /////////////////////////////////////////////////////////
  // Test #1: Tentative P, no rebalance
  /////////////////////////////////////////////////////////
  try {
    // Build tentative prolongator
    build_test_prolongator<CrsMatrixType>(A,P);

    // Build R
    Tpetra::RowMatrixTransposer<double, LO, GO, Node> transposer(P);
    R = transposer.createTranspose();

    ArrayRCP<const size_t> rowptr;
    ArrayRCP<const LO> colind;
    ArrayRCP<const double> vals;
    R->getAllValues(rowptr,colind,vals);

    // Form AP
    AP = rcp (new CrsMatrixType(A->getRowMap(),0));
    Tpetra::MatrixMatrix::Multiply(*A,false,*P,false,*AP);

    // Form RAP
    RAP = rcp (new CrsMatrixType(R->getRangeMap(),0));
    Tpetra::MatrixMatrix::Multiply(*R,false,*AP,false,*RAP);

    total_err+=test_err;
  } catch (std::exception& e) {
    err << "Proc " << myRank << ": " << e.what ();
    lclErr = 1;
  }

  reduceAll<int, int> (*Comm, REDUCE_MAX, lclErr, outArg (gblErr));
  if (gblErr != 0) {
    for (int r = 0; r < Comm->getSize (); ++r) {
      if (r == myRank) {
        cerr << err.str () << endl;
      }
      Comm->barrier ();
      Comm->barrier ();
      Comm->barrier ();
    }
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Test failed!");
  }

  /////////////////////////////////////////////////////////
  // Test #2: SaP, plus rebalancing
  /////////////////////////////////////////////////////////
  try {
    // Build tentative prolongator
    build_test_prolongator<CrsMatrixType>(A,Ptent);

    // Build SA-stype P
    P = rcp (new CrsMatrixType(A->getRowMap(),0));
    Tpetra::MatrixMatrix::Multiply(*A,false,*Ptent,false,*P);

    // Build R
    Tpetra::RowMatrixTransposer<double, LO, GO, Node> transposer(P);
    R = transposer.createTranspose();

    // Form AP
    AP = rcp (new CrsMatrixType(A->getRowMap(),0));
    Tpetra::MatrixMatrix::Multiply(*A,false,*P,false,*AP);

    // Form RAP
    RAP = rcp (new CrsMatrixType(R->getRangeMap(),0));
    Tpetra::MatrixMatrix::Multiply(*R,false,*AP,false,*RAP);

    // "Rebalanced" Map (based on the *RowMap* of RAP)
    build_test_map(RAP->getRowMap(),Map0);

    // Build Importer
    Import0 = rcp(new ImportType(RAP->getRowMap(),Map0));

    // Rebalance P
    ImportTemp = rcp(new ImportType(Import0->getTargetMap(),P->getColMap()));
    P->replaceDomainMapAndImporter(Import0->getTargetMap(),ImportTemp);

    // Rebalance R
    R = Tpetra::importAndFillCompleteCrsMatrix<CrsMatrixType>(R,*Import0,Teuchos::null,Import0->getTargetMap(),Teuchos::null);

    total_err+=test_err;
  }
  catch (std::exception& e) {
    err << "Proc " << myRank << ": " << e.what ();
    lclErr = 1;
  }

  reduceAll<int, int> (*Comm, REDUCE_MAX, lclErr, outArg (gblErr));
  if (gblErr != 0) {
    for (int r = 0; r < Comm->getSize (); ++r) {
      if (r == myRank) {
        cerr << err.str () << endl;
      }
      Comm->barrier ();
      Comm->barrier ();
      Comm->barrier ();
    }
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Test failed!");
  }

  TEST_EQUALITY(total_err,0);
}



TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( RemoteOnlyImport, Basic, LO, GO )  {
// Test the remotePIDs Tpetra::Import constructor
  RCP<const Comm<int> > Comm = getDefaultComm();
  typedef Tpetra::Map<LO, GO> MapType;
  typedef Tpetra::Import<LO, GO> ImportType;
  typedef Tpetra::Vector<double, LO, GO> VectorType;
  typedef Tpetra::CrsMatrix<double, LO, GO> CrsMatrixType;
  RCP<CrsMatrixType> A,B;

  RCP<const ImportType> Import1, Import2;
  RCP<VectorType> SourceVector, TargetVector, TestVector;
  RCP<MapType> Map0;
  const int myRank = Comm->getRank ();
  int total_err=0;
  int test_err=0;

  std::ostringstream err;
  int lclErr = 0;
  int gblErr = 0;

  // Build the sample matrix
  try {
    build_test_matrix<CrsMatrixType>(Comm,A);
  } catch (std::exception& e) {
    err << "Proc " << myRank << ": " << e.what ();
    lclErr = 1;
  }

  reduceAll<int, int> (*Comm, REDUCE_MAX, lclErr, outArg (gblErr));
  if (gblErr != 0) {
    for (int r = 0; r < Comm->getSize (); ++r) {
      if (r == myRank) {
        cerr << err.str () << endl;
      }
      Comm->barrier ();
      Comm->barrier ();
      Comm->barrier ();
    }
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Test failed!");
  }

  // Grab its importer
  Import1 = A->getGraph()->getImporter();

  // Only test in parallel
  if(Comm->getSize()==1) { TEST_EQUALITY(0,0); return;}

  // Build the remote-only map
  try {
    build_remote_only_map<ImportType,MapType>(Import1,Map0);
  } catch (std::exception& e) {
    err << "Proc " << myRank << ": " << e.what ();
    lclErr = 1;
  }

  reduceAll<int, int> (*Comm, REDUCE_MAX, lclErr, outArg (gblErr));
  if (gblErr != 0) {
    for (int r = 0; r < Comm->getSize (); ++r) {
      if (r == myRank) {
        cerr << err.str () << endl;
      }
      Comm->barrier ();
      Comm->barrier ();
      Comm->barrier ();
    }
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Test failed!");
  }

  // Vectors
  SourceVector = rcp(new VectorType(Import1->getSourceMap()));
  TargetVector = rcp(new VectorType(Import1->getTargetMap()));
  TestVector   = rcp(new VectorType(Map0));

  /////////////////////////////////////////////////////////
  // Test #1: Import & Compare
  /////////////////////////////////////////////////////////
  try {
    // Import reference vector
    Teuchos::ScalarTraits< double >::seedrandom(24601);
    SourceVector->randomize();
    TargetVector->doImport(*SourceVector,*Import1,Tpetra::INSERT);

    // Build remote-only import
    Import2 = Import1->createRemoteOnlyImport(Map0);

    // Check validity
    bool vv = Tpetra::Import_Util::checkImportValidity(*Import2);
    if(!vv) lclErr=1;

    // Do remote-only import
    TestVector->doImport(*SourceVector,*Import2,Tpetra::INSERT);

    // Compare vector output
    Teuchos::ArrayRCP<const double> view1 = TargetVector->get1dView();
    Teuchos::ArrayRCP<const double> view2 = TestVector->get1dView();
    double diff=0;
    size_t NumComps = Map0->getNodeNumElements();
    for(size_t i=0; i < NumComps; i++) {
      const size_t j = (size_t) Import1->getTargetMap ()->getLocalElement (Map0->getGlobalElement (static_cast<LO> (i)));
      diff += std::abs (view1[j] - view2[i]);
    }
    test_err = (diff > 1e-10) ? 1 : 0;
    total_err+=test_err;
  } catch (std::exception& e) {
    err << "Proc " << myRank << ": " << e.what ();
    lclErr = 1;
  }

  reduceAll<int, int> (*Comm, REDUCE_MAX, lclErr, outArg (gblErr));
  if (gblErr != 0) {
    for (int r = 0; r < Comm->getSize (); ++r) {
      if (r == myRank) {
        cerr << err.str () << endl;
      }
      Comm->barrier ();
      Comm->barrier ();
      Comm->barrier ();
    }
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Test failed!");
  }

  TEST_EQUALITY(total_err,0);
}



TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Import_Util,GetTwoTransferOwnershipVector, LO, GO )  {
  RCP<const Comm<int> > Comm = getDefaultComm();
  typedef Tpetra::Map<LO, GO> MapType;
  typedef Tpetra::Import<LO, GO> ImportType;
  typedef Tpetra::Vector<int, LO, GO> IntVectorType;
  RCP<const ImportType> ImportOwn, ImportXfer;
  RCP<MapType> Map0, Map1, Map2;

  // Get Rank
  const int NumProc = Comm->getSize ();
  const int MyPID   = Comm->getRank ();
  if(NumProc==1) {TEST_EQUALITY(0,0); return;}

  typedef Tpetra::global_size_t GST;
  const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

  std::ostringstream err;
  int lclErr = 0;
  int gblErr = 0;
  const int num_per_proc = 10;

  // Map0  - Even
  GO num_global = num_per_proc * NumProc;
  Map0 = rcp(new MapType(num_global,0,Comm));


  // Map1 - Cycled
  Teuchos::Array<GO> map1_gids(num_per_proc);
  for(int i=0; i<num_per_proc; i++) {
    map1_gids[i] = MyPID + i*NumProc;
  }
  Map1 = rcp(new MapType(INVALID,map1_gids,0,Comm));

  // Map2 - Reversed
  Teuchos::Array<GO> map2_gids(num_per_proc);
  GO map2_base = (NumProc-MyPID-1)*num_per_proc;
  for(int i=0; i<num_per_proc; i++) {
    map2_gids[i] = map2_base + i;
  }
  Map2 = rcp(new MapType(INVALID,map2_gids,0,Comm));

  // Importer / Exporter
  ImportOwn = rcp(new ImportType(Map0,Map1));
  ImportXfer = rcp(new ImportType(Map1,Map2));

  // Get the owned PID's list
  IntVectorType ownership(Map2);
  Tpetra::Import_Util::getTwoTransferOwnershipVector(*ImportOwn,false,*ImportXfer,false,ownership);
  Teuchos::ArrayRCP<const int> odata = ownership.getData();

#if 0
  {
    std::ostringstream oss;
    oss<<"["<<MyPID<<"] Map0 = ";
    for(int i=0; i<num_per_proc; i++)
      oss<<Map0->getGlobalElement(i)<<" " ;
    oss<<"\n["<<MyPID<<"] Map1 = ";
    for(int i=0; i<num_per_proc; i++)
      oss<<Map1->getGlobalElement(i)<<" " ;
    oss<<"\n["<<MyPID<<"] Map2 = ";
    for(int i=0; i<num_per_proc; i++)
      oss<<Map2->getGlobalElement(i)<<" ";
    oss<<"\n["<<MyPID<<"] Ownership = ";
    for(int i=0; i<num_per_proc; i++)
      oss<<odata[i]<< " ";
    std::cout<<oss.str()<<std::endl;
  }
#endif


  // Check answer [ownership(GID i) should contain the owning PID in Map0]
  for(size_t i=0; i<Map2->getNodeNumElements(); i++) {
    GO GID  = Map2->getGlobalElement(i);
    int PID = (int)(GID / num_per_proc);
    if (odata[i]!=PID) {
      printf("For GID %d on PID %d expected owner %d; got %d\n",(int)GID,MyPID,PID,odata[i]);
      lclErr++;
    }
  }

  // Sum Errors
  reduceAll<int, int> (*Comm, REDUCE_MAX, lclErr, outArg (gblErr));
  if (gblErr != 0) {
    for (int r = 0; r <NumProc; ++r) {
      if (r == MyPID) {
        cerr << err.str () << endl;
      }
      Comm->barrier ();
      Comm->barrier ();
      Comm->barrier ();
    }
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Test failed!");
  }

  TEST_EQUALITY(gblErr,0);
}

  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP_LO_GO( LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Import_Util, PackAndPrepareWithOwningPIDs, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Import, AdvancedConstructors, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( RemoteOnlyImport, Basic, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Import_Util, LowCommunicationMakeColMapAndReindex, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsGraphImportExport, doImport, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Import_Util, GetPids, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Import_Util, GetTwoTransferOwnershipVector, LO, GO )

#define UNIT_TEST_GROUP_SC_LO_GO( SC, LO, GO )                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrixImportExport, doImport, LO, GO, SC ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( FusedImportExport, doImport, LO, GO, SC ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Import_Util, UnpackAndCombineWithOwningPIDs, LO, GO, SC ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( FusedImportExport, MueLuStyle, LO, GO, SC )

  // Note: This test fails.  Should fix later.
  //      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( ReverseImportExport, doImport, ORDINAL, SCALAR )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  // Test CrsGraph for all LO, GO template parameter combinations, and
  // the default Node type.
  TPETRA_INSTANTIATE_LG( UNIT_TEST_GROUP_LO_GO )

  // Test CrsMatrix for all Scalar, LO, GO template parameter
  // combinations, and the default Node type.
  TPETRA_INSTANTIATE_SLG_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP_SC_LO_GO )

}
