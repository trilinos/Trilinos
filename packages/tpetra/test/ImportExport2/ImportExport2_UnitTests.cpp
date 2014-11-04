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
// Some Macro Magic to ensure that if CUDA and KokkosCompat is enabled
// only the .cu version of this file is actually compiled
#include <Tpetra_config.h>

#include <Tpetra_TestingUtilities.hpp>
#include "Teuchos_UnitTestHarness.hpp"

#include <map>
#include <vector>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_Tuple.hpp>
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Distributor.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Import_Util.hpp"
#include "Tpetra_Import_Util2.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"
#include "TpetraExt_MatrixMatrix.hpp"

namespace {
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
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::ScalarTraits;
  using Teuchos::tuple;

  using Tpetra::ADD;
  using Tpetra::createContigMap;
  using Tpetra::CrsGraph;
  using Tpetra::CrsMatrix;
  using Tpetra::DefaultPlatform;
  using Tpetra::DynamicProfile;
  using Tpetra::Export;
  using Tpetra::global_size_t;
  using Tpetra::Import;
  using Tpetra::INSERT;
  using Tpetra::Map;
  using Tpetra::REPLACE;
  using Tpetra::StaticProfile;

  using std::cout;
  using std::ostream_iterator;
  using std::endl;

  typedef DefaultPlatform::DefaultPlatformType::NodeType Node;

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
      return DefaultPlatform::getDefaultPlatform ().getComm ();
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

  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( CrsGraphImportExport, doImport, Ordinal )
  {
    using Teuchos::VERB_EXTREME;
    using Teuchos::VERB_NONE;

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm and node
    const RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    if (numImages < 2) {
      // This test is more meaningful when run with multiple
      // processes.
      return;
    }


    // Prepare for verbose output, if applicable.
    Teuchos::EVerbosityLevel verbLevel = verbose ? VERB_EXTREME : VERB_NONE;
    const bool doPrint = includesVerbLevel (verbLevel, VERB_EXTREME, true);
    if (doPrint) {
      out << "CrsGraphImportExport unit test" << endl;
    }
    OSTab tab (rcpFromRef (out)); // Add one tab level

    //Create a Map that is evenly-distributed, and another that has all
    //elements on proc 0.
    {
      const int tgt_num_local_elements = 3;
      const int src_num_local_elements =
        (myImageID == 0 ? numImages*tgt_num_local_elements : 0);

      // create Maps
      if (doPrint) {
        out << "Creating source and target Maps" << endl;
      }
      RCP<const Map<Ordinal> > src_map =
        createContigMap<Ordinal,Ordinal> (INVALID, src_num_local_elements, comm);
      RCP<const Map<Ordinal> > tgt_map =
        createContigMap<Ordinal,Ordinal> (INVALID, tgt_num_local_elements, comm);

      // create CrsGraph objects
      if (doPrint) {
        out << "Creating source and target CrsGraphs" << endl;
      }
      RCP<CrsGraph<Ordinal> > src_graph =
        rcp (new CrsGraph<Ordinal>(src_map, 1, DynamicProfile,
                                   getCrsGraphParameterList ()));
      RCP<CrsGraph<Ordinal> > tgt_graph =
        rcp (new CrsGraph<Ordinal>(tgt_map, 1, DynamicProfile,
                                   getCrsGraphParameterList ()));

      // Create a simple diagonal source graph.
      if (doPrint) {
        out << "Filling source CrsGraph" << endl;
      }
      Array<Ordinal> diag(1);
      Ordinal row = 0;
      for (size_t i = 0; i < src_map->getNodeNumElements (); ++i, ++row) {
        Ordinal globalrow = src_map->getGlobalElement (row);
        diag[0] = globalrow;
        src_graph->insertGlobalIndices (globalrow, diag ());
      }

      // Import from the source graph to the target graph.
      if (doPrint) {
        out << "Importing from source to target CrsGraph" << endl;
      }
      Import<Ordinal> importer (src_map, tgt_map, getImportParameterList ());
      tgt_graph->doImport (*src_graph, importer, INSERT);
      tgt_graph->fillComplete ();

      // Loop through the target graph and make sure it is diagonal.
      if (doPrint) {
        out << "Verifying target CrsGraph" << endl;
      }
      row = 0;
      for (size_t i = 0; i < tgt_map->getNodeNumElements (); ++i, ++row) {
        ArrayView<const Ordinal> rowview;
        tgt_graph->getLocalRowView( row, rowview );
        TEST_EQUALITY(rowview.size(), 1);
        TEST_EQUALITY(rowview[0], row);
      }
    }

    // For the next test, we need an even number of processes.
    // Skip this test otherwise.
    if (numImages%2 == 0) {
      // Create Maps that are distributed differently but have the
      // same global number of elements. The source map will have 3
      // elements on even-numbered processes and 5 on odd-numbered
      // processes. The target map will have 4 elements on each
      // process.
      Ordinal src_num_local = 5;
      if (myImageID % 2 == 0) {
        src_num_local = 3;
      }
      Ordinal tgt_num_local = 4;

      RCP<const Map<Ordinal> > src_map =
        createContigMap<Ordinal,Ordinal> (INVALID, src_num_local, comm);
      RCP<const Map<Ordinal> > tgt_map =
        createContigMap<Ordinal,Ordinal> (INVALID, tgt_num_local, comm);

      RCP<CrsGraph<Ordinal> > src_graph =
        rcp (new CrsGraph<Ordinal> (src_map, 24, DynamicProfile,
                                    getCrsGraphParameterList ()));
      RCP<CrsGraph<Ordinal> > tgt_graph =
        rcp (new CrsGraph<Ordinal> (tgt_map, 24, DynamicProfile,
                                    getCrsGraphParameterList ()));

      // This time make src_graph be a full lower-triangular graph.
      // Each row of column indices will have length 'globalrow'+1,
      // and contain column indices 0 .. 'globalrow'.
      Array<Ordinal> cols(1);
      for (Ordinal globalrow = src_map->getMinGlobalIndex ();
           globalrow <= src_map->getMaxGlobalIndex (); ++globalrow) {
        cols.resize(globalrow+1);
        for (Ordinal col = 0; col < globalrow+1; ++col) {
          cols[col] = col;
        }
        src_graph->insertGlobalIndices (globalrow, cols ());
      }

      Import<Ordinal> importer (src_map, tgt_map, getImportParameterList ());
      tgt_graph->doImport (*src_graph, importer, INSERT);

      src_graph->fillComplete ();
      tgt_graph->fillComplete ();

      // Loop through tgt_graph and make sure that each row has length
      // globalrow+1 and has the correct contents.
      RCP<const Map<Ordinal> > colmap = tgt_graph->getColMap ();

      for (Ordinal globalrow = tgt_map->getMinGlobalIndex ();
           globalrow <= tgt_map->getMaxGlobalIndex ();
           ++globalrow) {
        Ordinal localrow = tgt_map->getLocalElement (globalrow);
        ArrayView<const Ordinal> rowview;
        tgt_graph->getLocalRowView (localrow, rowview);
        TEST_EQUALITY(rowview.size(), globalrow+1);
        for (Ordinal j = 0; j < globalrow+1; ++j) {
          TEST_EQUALITY(colmap->getGlobalElement(rowview[j]), j);
        }
      }
    }
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrixImportExport, doImport, Ordinal, Scalar )
  {
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();

    // Parameters for CrsMatrix.  We'll let all CrsMatrix instances
    // use the same parameter list.
    RCP<ParameterList> crsMatPlist = getCrsMatrixParameterList ();

    if (numImages < 2) {
      // Testing Import/Export is more meaningful with at least two processes.
      return;
    }

    // First test: Import from a source Map that has all elements on
    // process 0, to a target Map that is evenly distributed.
    {
      const Ordinal tgt_num_local_elements = 3;
      const Ordinal src_num_local_elements =
        (myImageID == 0) ? numImages*tgt_num_local_elements : 0;

      // Create Maps.
      RCP<const Map<Ordinal> > src_map =
        createContigMap<Ordinal,Ordinal> (INVALID,src_num_local_elements,comm);
      RCP<const Map<Ordinal> > tgt_map =
        createContigMap<Ordinal,Ordinal> (INVALID,tgt_num_local_elements,comm);

      // Create CrsMatrix objects.
      RCP<CrsMatrix<Scalar,Ordinal> > src_mat =
        rcp (new CrsMatrix<Scalar,Ordinal> (src_map, 1, StaticProfile,
                                            crsMatPlist));
      RCP<CrsMatrix<Scalar,Ordinal> > tgt_mat =
        rcp (new CrsMatrix<Scalar,Ordinal> (tgt_map, 1, StaticProfile,
                                            crsMatPlist));

      // Create a simple diagonal source graph.
      for (Ordinal globalrow = src_map->getMinGlobalIndex();
           globalrow <= src_map->getMaxGlobalIndex();
           ++globalrow) {
        src_mat->insertGlobalValues (globalrow,
                                     tuple<Ordinal> (globalrow),
                                     tuple<Scalar> (globalrow));
      }
      src_mat->fillComplete ();

      // Create the importer
      Import<Ordinal> importer (src_map, tgt_map, getImportParameterList ());
      // Do the import, and fill-complete the target matrix.
      tgt_mat->doImport (*src_mat, importer, INSERT);
      tgt_mat->fillComplete ();

      // Loop through tgt_mat and make sure it is diagonal.
      for (Ordinal localrow = tgt_map->getMinLocalIndex();
           localrow <= tgt_map->getMaxLocalIndex();
           ++localrow)
      {
        ArrayView<const Ordinal> rowinds;
        ArrayView<const Scalar>  rowvals;
        tgt_mat->getLocalRowView(localrow, rowinds, rowvals);
        TEST_EQUALITY_CONST(rowinds.size(), 1);
        TEST_EQUALITY(rowinds[0], as<Ordinal>(localrow));
        TEST_EQUALITY_CONST(rowvals.size(), 1);
        TEST_EQUALITY(rowvals[0], as<Scalar>(localrow));
      }

      // Test the all-in-one import and fill complete nonmember
      // constructor.  The returned matrix should also be diagonal and
      // should equal tgt_mat.
      Teuchos::ParameterList dummy;
      typedef CrsMatrix<Scalar, Ordinal> crs_type;
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

      Array<Ordinal> tgtRowInds;
      Array<Scalar>  tgtRowVals;
      Array<Ordinal> tgt2RowInds;
      Array<Scalar>  tgt2RowVals;
      for (Ordinal localrow = tgt_map->getMinLocalIndex();
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
        for (size_type k = 0; k < as<size_type> (tgtNumEntries); ++k) {
          TEST_EQUALITY(tgtRowInds[k], tgt2RowInds[k]);
          // The "out" and "success" variables should have been
          // automatically defined by the unit test framework, in case
          // you're wondering where they came from.
          TEUCHOS_TEST_FLOATING_EQUALITY(tgtRowVals[k], tgt2RowVals[k], tol, out, success);
        } // for each entry in the current row
      } // for each row in the matrix
    } // end of the first test


    // For the next test, we need an even number of processes:
    if (numImages%2 == 0) {
      // Create Maps that are distributed differently but have the
      // same global number of elements. The source-map will have 3
      // elements on even-numbered processes and 5 on odd-numbered
      // processes. The target-map will have 4 elements on each
      // process.
      const Ordinal src_num_local = (myImageID%2 == 0 ? 3 : 5);
      const Ordinal tgt_num_local = 4;

      RCP<const Map<Ordinal> > src_map =
        createContigMap<Ordinal,Ordinal> (INVALID, src_num_local, comm);
      RCP<const Map<Ordinal> > tgt_map =
        createContigMap<Ordinal,Ordinal> (INVALID, tgt_num_local, comm);

      RCP<CrsMatrix<Scalar,Ordinal> > src_mat =
        rcp (new CrsMatrix<Scalar,Ordinal> (src_map, 24, DynamicProfile, crsMatPlist));
      RCP<CrsMatrix<Scalar,Ordinal> > tgt_mat =
        rcp (new CrsMatrix<Scalar,Ordinal> (tgt_map, 24, DynamicProfile, crsMatPlist));

      // This time make src_mat a full lower-triangular matrix.  Each
      // row of column-indices will have length 'globalrow', and
      // contain column-indices 0 .. 'globalrow'-1
      Array<Ordinal> cols(1);
      Array<Scalar>  vals(1);
      for (Ordinal globalrow = src_map->getMinGlobalIndex();
           globalrow <= src_map->getMaxGlobalIndex();
           ++globalrow) {
        if (globalrow > 0) {
          cols.resize(globalrow);
          vals.resize(globalrow);
          for (Ordinal col=0; col<globalrow; ++col) {
            cols[col] = as<Ordinal>(col);
            vals[col] = as<Scalar>(col);
          }
          src_mat->insertGlobalValues (globalrow, cols (), vals ());
        }
      }

      Import<Ordinal> importer (src_map, tgt_map, getImportParameterList ());
      tgt_mat->doImport(*src_mat, importer, Tpetra::INSERT);
      tgt_mat->fillComplete();

      // now we're going to loop through tgt_mat and make sure that
      // each row has length 'globalrow' and has the correct contents:
      const Teuchos::RCP<const Map<Ordinal> > colmap = tgt_mat->getColMap();

      for (Ordinal globalrow=tgt_map->getMinGlobalIndex(); globalrow<=tgt_map->getMaxGlobalIndex(); ++globalrow)
      {
        Ordinal localrow = tgt_map->getLocalElement(globalrow);
        ArrayView<const Ordinal> rowinds;
        ArrayView<const Scalar> rowvals;
        tgt_mat->getLocalRowView(localrow, rowinds, rowvals);
        TEST_EQUALITY(rowinds.size(), globalrow);
        TEST_EQUALITY(rowvals.size(), globalrow);
        for (Teuchos_Ordinal j=0; j<rowinds.size(); ++j) {
          TEST_EQUALITY( colmap->getGlobalElement(rowinds[j]), as<Ordinal>(j) );
          TEST_EQUALITY( rowvals[j], as<Scalar>(j)  );
        }
      }
    }
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
void build_remote_only_map(const RCP<const ImportType> & Import, RCP<MapType> & newMap) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef typename MapType::local_ordinal_type LO;
  typedef typename MapType::local_ordinal_type GO;
  typedef Tpetra::global_size_t GST;

  if (Import.is_null ()) {
    return;
  }

  RCP<const MapType> targetMap = Import->getTargetMap();
  size_t NumRemotes = Import->getNumRemoteIDs();
  ArrayView<const LO> oldRemoteLIDs = Import->getRemoteLIDs();
  Array<GO> newRemoteGIDs(NumRemotes);

  for(size_t i=0; i < NumRemotes; i++)
    newRemoteGIDs[i] = targetMap->getGlobalElement(oldRemoteLIDs[i]);

  const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
  newMap = rcp (new MapType (INVALID, newRemoteGIDs, targetMap->getIndexBase (),
                             targetMap->getComm (), targetMap->getNode ()));
}


// ===============================================================================
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( FusedImportExport, doImport, Ordinal, Scalar )  {
   RCP<const Comm<int> > Comm = getDefaultComm();
   typedef Tpetra::CrsMatrix<Scalar,Ordinal,Ordinal> CrsMatrixType;
   typedef Tpetra::Map<Ordinal,Ordinal> MapType;
   typedef Tpetra::Import<Ordinal,Ordinal> ImportType;
   typedef Tpetra::Export<Ordinal,Ordinal> ExportType;
   typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MagType;

   RCP<CrsMatrixType> A, B, C;
   RCP<const MapType> Map1, Map2;
   RCP<MapType> Map3;

   RCP<ImportType> Import1;
   RCP<ExportType> Export1;
   int MyPID = Comm->getRank();
   double diff;
   int total_err=0;
   MagType diff_tol = 1e4*Teuchos::ScalarTraits<Scalar>::eps();

   Ordinal INVALID = Teuchos::OrdinalTraits<Ordinal>::invalid();
   build_test_matrix<CrsMatrixType>(Comm,A);

   /////////////////////////////////////////////////////////
   // Test #1: Tridiagonal Matrix; Migrate to Proc 0
   /////////////////////////////////////////////////////////
   {
     global_size_t num_global = A->getRowMap()->getGlobalNumElements();

     // New map with all on Proc1
     if(MyPID==0) Map1 = rcp(new MapType(num_global,(size_t)num_global,0,Comm));
     else Map1 = rcp(new MapType(num_global,(size_t)0,0,Comm));

     // Execute fused import
     Import1 = rcp(new ImportType(A->getRowMap(),Map1));
     B = Tpetra::importAndFillCompleteCrsMatrix<CrsMatrixType>(A,*Import1);
     diff=test_with_matvec<CrsMatrixType>(*A,*B);
     if(diff > diff_tol) {
       if(MyPID==0) cout<<"FusedImport: Test #1 FAILED with norm diff = "<<diff<<"."<<endl;
       total_err--;
     }

     // Execute fused export
     Export1 = rcp(new ExportType(A->getRowMap(),Map1));
     B = Tpetra::exportAndFillCompleteCrsMatrix<CrsMatrixType>(A,*Export1);
     diff=test_with_matvec<CrsMatrixType>(*A,*B);
     if(diff > diff_tol) {
       if(MyPID==0) cout<<"FusedExport: Test #1 FAILED with norm diff = "<<diff<<"."<<endl;
       total_err--;
     }

     Comm->barrier ();
   }

   /////////////////////////////////////////////////////////
   // Test #2: Tridiagonal Matrix; Locally Reversed Map
   /////////////////////////////////////////////////////////
   {
     size_t num_local = A->getRowMap()->getNodeNumElements();

     Teuchos::Array<Ordinal> MyGIDs(num_local);
     for(size_t i=0; i<num_local; i++)
       MyGIDs[i] = A->getRowMap()->getGlobalElement(num_local-i-1);

     Map1 = rcp(new MapType(INVALID,MyGIDs(),0,Comm));

     // Execute fused import
     Import1 = rcp(new ImportType(A->getRowMap(),Map1));
     B = Tpetra::importAndFillCompleteCrsMatrix<CrsMatrixType>(A,*Import1);
     diff=test_with_matvec<CrsMatrixType>(*A,*B);
     if(diff > diff_tol) {
       if(MyPID==0) cout<<"FusedImport: Test #2 FAILED with norm diff = "<<diff<<"."<<endl;
       total_err--;
     }

     // Execute fused export
     Export1 = rcp(new ExportType(A->getRowMap(),Map1));
     B = Tpetra::exportAndFillCompleteCrsMatrix<CrsMatrixType>(A,*Export1);
     diff=test_with_matvec<CrsMatrixType>(*A,*B);
     if(diff > diff_tol) {
       if(MyPID==0) cout<<"FusedExport: Test #2 FAILED with norm diff = "<<diff<<"."<<endl;
       total_err--;
     }
   }

   /////////////////////////////////////////////////////////
   // Test #3: Tridiagonal Matrix; Globally Reversed Map
   /////////////////////////////////////////////////////////

   // Skipped.


   /////////////////////////////////////////////////////////
   // Test #4: Tridiagonal Matrix; MMM style halo import
   /////////////////////////////////////////////////////////
   {
     // Assume we always own the diagonal
     size_t num_local = A->getNodeNumCols()-A->getNodeNumRows();
     Teuchos::Array<Ordinal> MyGIDs(num_local);

     size_t idx=0;
     for(Ordinal i=0; (size_t)i<A->getNodeNumCols(); i++)
       if(A->getRowMap()->getLocalElement(A->getColMap()->getGlobalElement(i)) == INVALID){
         MyGIDs[idx] = A->getColMap()->getGlobalElement(i);
         idx++;
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
       if(MyPID==0) cout<<"FusedImport: Test #4 FAILED with norm diff = "<<diff<<"."<<endl;
       total_err--;
     }

     // Execute fused export
     Export1 = rcp(new ExportType(A->getRowMap(),Map1));
     B = Tpetra::exportAndFillCompleteCrsMatrix<CrsMatrixType>(A,*Export1);
     diff=test_with_matvec<CrsMatrixType>(*B,*C);
     if(diff > diff_tol) {
       if(MyPID==0) cout<<"FusedExport: Test #4 FAILED with norm diff = "<<diff<<"."<<endl;
       total_err--;
     }
   }

  /////////////////////////////////////////////////////////
  // Test 5: Tridiagonal Matrix; Migrate to Proc 0, Replace Maps
  /////////////////////////////////////////////////////////
  {
    // New map with all on Proc 0 / 2
    build_test_map(A->getRowMap(),Map3);

    // Execute fused import
    Import1 = rcp(new ImportType(A->getRowMap(),Map3));
    B = Tpetra::importAndFillCompleteCrsMatrix<CrsMatrixType>(A,*Import1,Map3,Map3);
    diff=test_with_matvec<CrsMatrixType>(*A,*B);
    if(diff > diff_tol) {
      if(MyPID==0) cout<<"FusedImport: Test #5 FAILED with norm diff = "<<diff<<"."<<endl;
      total_err--;
    }

    // Execute fused export
    Export1 = rcp(new ExportType(A->getRowMap(),Map3));
    B = Tpetra::exportAndFillCompleteCrsMatrix<CrsMatrixType>(A,*Export1,Map3,Map3);
    diff=test_with_matvec<CrsMatrixType>(*A,*B);
    if(diff > diff_tol) {
      if(MyPID==0) cout<<"FusedExport: Test #5 FAILED with norm diff = "<<diff<<"."<<endl;
      total_err--;
    }
  }

  /////////////////////////////////////////////////////////
  // Test 6: Tridiagonal Matrix; Migrate to Proc 0, Replace Comm
  /////////////////////////////////////////////////////////
  {
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
      if(MyPID==0) cout<<"FusedImport: Test #6 FAILED with norm diff = "<<diff<<"."<<endl;
      total_err--;
    }

    // Execute fused export constructor
    Export1 = rcp(new ExportType(A->getRowMap(),Map3));
    B = Tpetra::exportAndFillCompleteCrsMatrix<CrsMatrixType>(A,*Export1,Map3,Map3,rcp(&params,false));

    diff=test_with_matvec_reduced_maps<CrsMatrixType,MapType>(*A,*B,*Map3);
    if(diff > diff_tol){
      if(MyPID==0) cout<<"FusedExport: Test #6 FAILED with norm diff = "<<diff<<"."<<endl;
      total_err--;
    }

  }


  /////////////////////////////////////////////////////////
  // Test 7: Tridiagonal Matrix; Migrate to Proc 0, Reverse Mode
  /////////////////////////////////////////////////////////
  {
    global_size_t num_global = A->getRowMap()->getGlobalNumElements();

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
      if(MyPID==0) cout<<"FusedImport: Test #7 FAILED with norm diff = "<<diff<<"."<<endl;
      total_err--;
    }

    // Execute fused export constructor
    Export1 = rcp(new ExportType(Map1,A->getRowMap()));
    B = Tpetra::exportAndFillCompleteCrsMatrix<CrsMatrixType>(A,*Export1,Map1,Map1,rcp(&params,false));

    diff=test_with_matvec<CrsMatrixType>(*A,*B);
    if(diff > diff_tol){
      if(MyPID==0) cout<<"FusedExport: Test #7 FAILED with norm diff = "<<diff<<"."<<endl;
      total_err--;
    }

  }

  /////////////////////////////////////////////////////////
  // Test #8: Diagonal matrix w/ overlapping entries
  /////////////////////////////////////////////////////////
  {
    build_test_matrix_with_row_overlap<CrsMatrixType>(Comm,A);

    Teuchos::ArrayRCP< const size_t > rowptr;
    Teuchos::ArrayRCP< const Ordinal > colind;
    Teuchos::ArrayRCP< const Scalar > vals;

    Map1 = A->getRangeMap();

    // Execute fused import constructor (reverse)
    Teuchos::ParameterList params;
    params.set("Reverse Mode",true);
    Import1 = rcp(new ImportType(Map1,A->getRowMap()));
    B = Tpetra::importAndFillCompleteCrsMatrix<CrsMatrixType>(A,*Import1,Map1,Map1,rcp(&params,false));

    diff=test_with_matvec<CrsMatrixType>(*B,*A);
    if(diff > diff_tol){
      if(MyPID==0) cout<<"FusedImport: Test #8 FAILED with norm diff = "<<diff<<"."<<endl;
      total_err--;
    }

    // Execute fused export constructor
    Export1 = rcp(new ExportType(A->getRowMap(),Map1));
    B = Tpetra::exportAndFillCompleteCrsMatrix<CrsMatrixType>(A,*Export1,Map1,Map1);
    diff=test_with_matvec<CrsMatrixType>(*B,*A);
    if(diff > diff_tol){
      if(MyPID==0) cout<<"FusedExport: Test #8 FAILED with norm diff = "<<diff<<"."<<endl;
      total_err--;
    }

  }

  TEST_EQUALITY(total_err,0);
}



TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( ReverseImportExport, doImport, Ordinal, Scalar )  {
  // NOTE: This test will fail.

  RCP<const Comm<int> > Comm = getDefaultComm();
  typedef Tpetra::Map<Ordinal,Ordinal> MapType;
  typedef Tpetra::Import<Ordinal,Ordinal> ImportType;
  typedef Tpetra::Export<Ordinal,Ordinal> ExportType;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MagType;
  typedef Tpetra::Vector<Scalar,Ordinal,Ordinal, Node> VectorType;
  typedef Tpetra::CrsMatrix<Scalar,Ordinal,Ordinal> CrsMatrixType;

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
  global_size_t num_global = A->getRowMap()->getGlobalNumElements();

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



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Import_Util, GetPids, Ordinal )  {
  // Unit Test the functionality in Tpetra_Import_Util
  RCP<const Comm<int> > Comm = getDefaultComm();
  typedef Tpetra::Map<Ordinal,Ordinal> MapType;
  typedef Tpetra::Import<Ordinal,Ordinal> ImportType;
  typedef Tpetra::Export<Ordinal,Ordinal> ExportType;
  typedef Tpetra::Vector<int,Ordinal,Ordinal, Node> IntVectorType;
  typedef Tpetra::CrsMatrix<double,Ordinal,Ordinal> CrsMatrixType;

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
  global_size_t num_global = A->getRowMap()->getGlobalNumElements();

  // Target Map - all on one proc
  size_t num_local = MyPID==0 ? num_global : 0;
  MapTarget = rcp(new MapType(num_global,num_local,0,Comm));

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
    Tpetra::Import_Util::getPids<Ordinal,Ordinal,Node>(*Import1,pids,false);

    // Compare
    Teuchos::ArrayRCP<const int>  TargetView = TargetVector->get1dView();
    for(size_t i=0; i < TargetVector->getLocalLength(); i++)
      test_err += (pids[i] - TargetView[i] > 0) ?  1 : ((pids[i] - TargetView[i] < 0) ?  1 : 0);
    total_err += test_err;
  }

  /////////////////////////////////////////////////////////
  // Test #2: getPidGidPairs
  /////////////////////////////////////////////////////////
  {
    test_err=0;
    // Generate PID vector via getPidGidPairs
    Teuchos::Array<std::pair<int,Ordinal> > pgids;
    Tpetra::Import_Util::getPidGidPairs<Ordinal,Ordinal,Node>(*Import1,pgids,false);

    // Compare
    Teuchos::ArrayRCP<const int>  TargetView = TargetVector->get1dView();
    for(size_t i=0; i < TargetVector->getLocalLength(); i++)
      test_err += (pgids[i].first - TargetView[i] > 0) ?  1 : ((pgids[i].first - TargetView[i] < 0) ?  1 : 0);
    total_err += test_err;
  }

  /////////////////////////////////////////////////////////
  // Test #3: getRemotePids
  /////////////////////////////////////////////////////////
  {
    test_err=0;
    Teuchos::Array<int> RemotePIDs;
    Tpetra::Import_Util::getRemotePIDs<Ordinal,Ordinal,Node>(*Import1,RemotePIDs);

    // We can't easily test this, so let's at least make sure it doesn't crash.
  }

  TEST_EQUALITY(total_err,0);
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Import_Util, PackAndPrepareWithOwningPIDs, Ordinal )  {
  // Unit Test the functionality in Tpetra_Import_Util
  RCP<const Comm<int> > Comm = getDefaultComm();
  typedef Tpetra::Import<Ordinal,Ordinal> ImportType;
  typedef Tpetra::CrsMatrix<double,Ordinal,Ordinal> CrsMatrixType;
  using Teuchos::av_reinterpret_cast;

  RCP<CrsMatrixType> A;
  int total_err=0;
  int test_err=0;
  Teuchos::Array<char> exports1, exports2;
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
    Tpetra::Import_Util::getPids<Ordinal,Ordinal,Node>(*Importer,pids,false);
    constantNumPackets2=0;
    numPackets2.resize(Importer->getExportLIDs().size());
    Tpetra::Import_Util::packAndPrepareWithOwningPIDs<double,Ordinal,Ordinal,Node>(*A,Importer->getExportLIDs(),exports2,numPackets2(),constantNumPackets2,Importer->getDistributor(),pids());

    // Loop through the parts that should be the same
    const size_t numExportLIDs = Importer->getExportLIDs().size();
    Teuchos::ArrayView<const Ordinal> exportLIDs = Importer->getExportLIDs();

    const int sizeOfPacket1 = sizeof(double) + sizeof(Ordinal);
    const int sizeOfPacket2 = sizeof(double) + sizeof(Ordinal) + sizeof(int);

    size_t offset1=0,offset2=0;
    for (size_t i = 0; i < numExportLIDs; i++) {
      ArrayView<const Ordinal> lidsView;
      ArrayView<const double>  valsView;
      A->getLocalRowView(exportLIDs[i], lidsView, valsView);
      const size_t curNumEntries = lidsView.size();

      ArrayView<char> gidsChar1 = exports1(offset1, curNumEntries*sizeof(Ordinal));
      ArrayView<char> valsChar1 = exports1(offset1+curNumEntries*sizeof(Ordinal), curNumEntries*sizeof(double));
      ArrayView<Ordinal> gids1  = av_reinterpret_cast<Ordinal>(gidsChar1);
      ArrayView<double>  vals1  = av_reinterpret_cast<double>(valsChar1);

      ArrayView<char> gidsChar2 = exports2(offset2, curNumEntries*sizeof(Ordinal));
      //      ArrayView<char> pidsChar2 = exports2(offset2+curNumEntries*sizeof(Ordinal), curNumEntries*sizeof(int));
      ArrayView<char> valsChar2 = exports2(offset2+curNumEntries*(sizeof(Ordinal)+sizeof(int)), curNumEntries*sizeof(double));
      ArrayView<Ordinal> gids2  = av_reinterpret_cast<Ordinal>(gidsChar2);
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
}




TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Import_Util, UnpackAndCombineWithOwningPIDs, Ordinal, Scalar)  {
  // Unit Test the functionality in Tpetra_Import_Util
  RCP<const Comm<int> > Comm = getDefaultComm();
  typedef Tpetra::Map<Ordinal,Ordinal> MapType;
  typedef Tpetra::Import<Ordinal,Ordinal> ImportType;
  typedef Tpetra::CrsMatrix<Scalar,Ordinal,Ordinal> CrsMatrixType;
  typedef typename Tpetra::CrsMatrix<Scalar,Ordinal,Ordinal>::packet_type PacketType;
  using Teuchos::av_reinterpret_cast;

  RCP<CrsMatrixType> A,B;
  int total_err=0;
  int test_err=0;
  int MyPID = Comm->getRank();
  RCP<const MapType> MapSource, MapTarget;
  Teuchos::Array<char> exports, imports;
  Teuchos::Array<size_t> numImportPackets, numExportPackets;
  size_t constantNumPackets=0;

  // Build sample matrix & sourceMap
  build_test_matrix<CrsMatrixType>(Comm,A);
  global_size_t num_global = A->getRowMap()->getGlobalNumElements();
  MapSource = A->getRowMap();

  // Target Map - all on one proc
  size_t num_local = MyPID==0 ? num_global : 0;
  MapTarget = rcp(new MapType(num_global,num_local,0,Comm));

  // Build Importer
  RCP<ImportType> Importer = rcp(new ImportType(MapSource,MapTarget));
  if(Importer != Teuchos::null)  {
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
    else Tpetra::Import_Util::getPids<Ordinal,Ordinal,Node>(*A->getGraph()->getImporter(),SourcePids,false);
    numExportPackets.resize(Importer->getExportLIDs().size());
    numImportPackets.resize(Importer->getRemoteLIDs().size());

    Tpetra::Import_Util::packAndPrepareWithOwningPIDs<Scalar,Ordinal,Ordinal,Node>(*A,Importer->getExportLIDs(),exports,numExportPackets(),constantNumPackets,distor,SourcePids());

    // Do the moral equivalent of doTransfer
    distor.doPostsAndWaits<size_t>(numExportPackets().getConst(), 1,numImportPackets());
    size_t totalImportPackets = 0;
    for(size_t i = 0; i < (size_t)numImportPackets.size(); i++) {
      totalImportPackets += numImportPackets[i];
    }
    imports.resize(totalImportPackets);
    distor.doPostsAndWaits<PacketType>(exports(),numExportPackets(),imports(),numImportPackets());

    // Run the count... which should get the same NNZ as the traditional import
    size_t nnz2=Tpetra::Import_Util::unpackAndCombineWithOwningPIDsCount<Scalar,Ordinal,Ordinal,Node>(*A,Importer->getRemoteLIDs(),imports(),numImportPackets(),
                                                                                                               constantNumPackets, distor,Tpetra::INSERT,Importer->getNumSameIDs(),
                                                                                                               Importer->getPermuteToLIDs(),Importer->getPermuteFromLIDs());
    if(nnz1!=nnz2) test_err++;
    total_err+=test_err;

    /////////////////////////////////////////////////////////
    // Test #2: Actual combine test
    /////////////////////////////////////////////////////////
    Teuchos::Array<size_t>  rowptr(MapTarget->getNodeNumElements()+1);
    Teuchos::Array<Ordinal> colind(nnz2);
    Teuchos::Array<Scalar>  vals(nnz2);
    Teuchos::Array<int>     TargetPids;

    Tpetra::Import_Util::unpackAndCombineIntoCrsArrays<Scalar,Ordinal,Ordinal,Node>(*A,Importer->getRemoteLIDs(),imports(),numImportPackets(),constantNumPackets,
                                                                                             distor,Tpetra::INSERT,Importer->getNumSameIDs(),Importer->getPermuteToLIDs(),
                                                                                             Importer->getPermuteFromLIDs(),MapTarget->getNodeNumElements(),nnz2,MyPID,rowptr(),
                                                                                             colind(),vals(),SourcePids(),TargetPids);
    // Do the comparison
    Teuchos::ArrayRCP<const size_t>  Browptr;
    Teuchos::ArrayRCP<const Ordinal> Bcolind;
    Teuchos::ArrayRCP<const Scalar>  Bvals;
    B->getAllValues(Browptr,Bcolind,Bvals);

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
      Tpetra::Import_Util::sortCrsEntries<Scalar,Ordinal>(rowptr(),colind(),vals());

      // Compare the gids
      for(size_t i=0; i < as<size_t>(rowptr.size()-1); i++) {
        for(size_t j=rowptr[i]; j < rowptr[i+1]; j++) {
          if(colind[j] != Bcolind[j] || vals[j] != Bvals[j]) {
            test_err++;
          }
        }
      }
    }

    total_err+=test_err;
  }

  TEST_EQUALITY(total_err,0);
}


TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Import_Util,LowCommunicationMakeColMapAndReindex, LocalOrdinal, GlobalOrdinal)  {
  // Test the colmap...
  RCP<const Comm<int> > Comm = getDefaultComm();
  typedef double Scalar;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal> MapType;
  typedef Tpetra::Import<LocalOrdinal,GlobalOrdinal> ImportType;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal> CrsMatrixType;
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
  ArrayRCP<const LocalOrdinal> colind;
  ArrayRCP<const Scalar> values;
  A->getAllValues(rowptr,colind,values);
  Acolmap = A->getColMap();
  Adomainmap = A->getDomainMap();

  // Get owning PID information
  size_t numMyCols = A->getColMap()->getNodeNumElements();
  RCP<const ImportType> Importer = A->getGraph()->getImporter();
  Teuchos::Array<int> AcolmapPIDs(numMyCols,-1);
  if(Importer!=Teuchos::null)
    Tpetra::Import_Util::getPids<LocalOrdinal,GlobalOrdinal,Node>(*Importer,AcolmapPIDs,true);

  // Build a "gid" version of colind & colind-sized pid list
  Array<GlobalOrdinal> colind_GID(colind.size());
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
    Teuchos::Array<LocalOrdinal> Bcolind_LID(colind.size());
    Tpetra::Import_Util::lowCommunicationMakeColMapAndReindex<LocalOrdinal,GlobalOrdinal,Node>(rowptr(),Bcolind_LID(),colind_GID(),Adomainmap,colind_PID(),BcolmapPIDs,Bcolmap);

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


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Import, AdvancedConstructors, Ordinal )  {
  // Test the remotePIDs Tpetra::Import constructor
  RCP<const Comm<int> > Comm = getDefaultComm();
  typedef Tpetra::Map<Ordinal,Ordinal> MapType;
  typedef Tpetra::Import<Ordinal,Ordinal> ImportType;
  typedef Tpetra::CrsMatrix<double,Ordinal,Ordinal> CrsMatrixType;
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
    Tpetra::Import_Util::getRemotePIDs<Ordinal,Ordinal,Node>(*Import1,pids);

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
    if(Import1->getNumRemoteIDs() != Import2->getNumRemoteIDs())
      test_err++;
    else {
      for(size_t i=0; i<Import1->getNumRemoteIDs(); i++)
        test_err += (Import1->getRemoteLIDs()[i]!=Import2->getRemoteLIDs()[i]);
    }
    if(Import1->getNumExportIDs() != Import2->getNumExportIDs())
      test_err++;
    else {
      for(size_t i=0; i<Import1->getNumExportIDs(); i++) {
        test_err += (Import1->getExportLIDs()[i]!=Import2->getExportLIDs()[i]);
        test_err += (Import1->getExportPIDs()[i]!=Import2->getExportPIDs()[i]);
      }
    }
    total_err += test_err;
  }

  /////////////////////////////////////////////////////////
  // Test #2: MueLu-style Transpose & Rebalance "R"
  /////////////////////////////////////////////////////////
  {
    // Create Transpose
    Tpetra::RowMatrixTransposer<double, Ordinal, Ordinal, Node> transposer(A);
    B = transposer.createTranspose();

    // Build Importer
    Import2 = rcp(new ImportType(B->getRowMap(),Map0));

    // Build the imported matrix
    Tpetra::importAndFillCompleteCrsMatrix<CrsMatrixType>(B,*Import2,Teuchos::null,Map0,Teuchos::null);

  }
  //


  TEST_EQUALITY(total_err,0);
}



TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( FusedImportExport, MueLuStyle, Ordinal, Scalar )  {
  // Test a muelu-style SA build and rebalance.  Kind of like Gangnam Style, but with more cows.
  RCP<const Comm<int> > Comm = getDefaultComm();
  typedef Tpetra::Map<Ordinal,Ordinal> MapType;
  typedef Tpetra::Import<Ordinal,Ordinal> ImportType;
  typedef Tpetra::CrsMatrix<double,Ordinal,Ordinal> CrsMatrixType;
  RCP<CrsMatrixType> A,Ptent,P,R,AP,RAP,rebalancedP;
  RCP<MapType> Map0;
  RCP<const ImportType> Import0,ImportTemp;

  int total_err=0;
  int test_err=0;

  // Build sample matrix
  build_test_matrix_wideband<CrsMatrixType>(Comm,A);

  /////////////////////////////////////////////////////////
  // Test #1: Tentative P, no rebalance
  /////////////////////////////////////////////////////////
  {
    // Build tentative prolongator
    build_test_prolongator<CrsMatrixType>(A,P);

    // Build R
    Tpetra::RowMatrixTransposer<double, Ordinal, Ordinal, Node> transposer(P);
    R = transposer.createTranspose();

    ArrayRCP< const size_t > rowptr;
    ArrayRCP< const Ordinal > colind;
    ArrayRCP< const double > vals;
    R->getAllValues(rowptr,colind,vals);

    // Form AP
    AP = rcp (new CrsMatrixType(A->getRowMap(),0));
    Tpetra::MatrixMatrix::Multiply(*A,false,*P,false,*AP);

    // Form RAP
    RAP = rcp (new CrsMatrixType(R->getRangeMap(),0));
    Tpetra::MatrixMatrix::Multiply(*R,false,*AP,false,*RAP);

    total_err+=test_err;
  }

  /////////////////////////////////////////////////////////
  // Test #2: SaP, plus rebalancing
  /////////////////////////////////////////////////////////
  {
    // Build tentative prolongator
    build_test_prolongator<CrsMatrixType>(A,Ptent);

    // Build SA-stype P
    P = rcp (new CrsMatrixType(A->getRowMap(),0));
    Tpetra::MatrixMatrix::Multiply(*A,false,*Ptent,false,*P);

    // Build R
    Tpetra::RowMatrixTransposer<double, Ordinal, Ordinal, Node> transposer(P);
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

  TEST_EQUALITY(total_err,0);
}



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RemoteOnlyImport, Basic, Ordinal )  {
// Test the remotePIDs Tpetra::Import constructor
  RCP<const Comm<int> > Comm = getDefaultComm();
  typedef Tpetra::Map<Ordinal,Ordinal> MapType;
  typedef Tpetra::Import<Ordinal,Ordinal> ImportType;
  typedef Tpetra::Vector<double,Ordinal,Ordinal, Node> VectorType;
  typedef Tpetra::CrsMatrix<double,Ordinal,Ordinal> CrsMatrixType;
  RCP<CrsMatrixType> A,B;

  RCP<const ImportType> Import1, Import2;
  RCP<VectorType> SourceVector, TargetVector, TestVector;
  RCP<MapType> Map0;
  int total_err=0;
  int test_err=0;

  // Build the sample matrix
  build_test_matrix<CrsMatrixType>(Comm,A);

  // Grab its importer
  Import1 = A->getGraph()->getImporter();

  // Only test in parallel
  if(Comm->getSize()==1) { TEST_EQUALITY(0,0); return;}

  // Build the remote-only map
  build_remote_only_map<ImportType,MapType>(Import1,Map0);

  // Vectors
  SourceVector = rcp(new VectorType(Import1->getSourceMap()));
  TargetVector = rcp(new VectorType(Import1->getTargetMap()));
  TestVector   = rcp(new VectorType(Map0));

  /////////////////////////////////////////////////////////
  // Test #1: Import & Compare
  /////////////////////////////////////////////////////////
  {
    // Import reference vector
    Teuchos::ScalarTraits< double >::seedrandom(24601);
    SourceVector->randomize();
    TargetVector->doImport(*SourceVector,*Import1,Tpetra::INSERT);

    // Build remote-only import
    Import2 = Import1->createRemoteOnlyImport(Map0);

    // Do remote-only import
    TestVector->doImport(*SourceVector,*Import2,Tpetra::INSERT);

    // Compare vector output
    Teuchos::ArrayRCP<const double> view1 = TargetVector->get1dView();
    Teuchos::ArrayRCP<const double> view2 = TestVector->get1dView();
    double diff=0;
    size_t NumComps = Map0->getNodeNumElements();
    for(size_t i=0; i < NumComps; i++) {
      size_t j = (size_t) Import1->getTargetMap()->getLocalElement(Map0->getGlobalElement((Ordinal)i));
      diff += std::abs( view1[j] - view2[i] );
    }
    test_err = (diff > 1e-10) ? 1 : 0;
    total_err+=test_err;
  }
  TEST_EQUALITY(total_err,0);
}



  //
  // INSTANTIATIONS
  //

#   define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( CrsGraphImportExport, doImport, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Import_Util, GetPids, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Import_Util, PackAndPrepareWithOwningPIDs, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Import, AdvancedConstructors, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( RemoteOnlyImport, Basic, ORDINAL )

#   define UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrixImportExport, doImport, ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( FusedImportExport, doImport, ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Import_Util, UnpackAndCombineWithOwningPIDs, ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( FusedImportExport,  MueLuStyle, ORDINAL, SCALAR )

#define UNIT_TEST_GROUP_LO_GO(LO, GO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Import_Util, LowCommunicationMakeColMapAndReindex, LO, GO)




  // Note: This test fails.  Should fix later.
  //      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( ReverseImportExport, doImport, ORDINAL, SCALAR )

  // test CrsGraph for <int,int,DefaultNode>
  // if explicit instantiation is enabled, this configuration is always built
  // if not, it is implicitly instantiated
  // therefore, this is always possible
  UNIT_TEST_GROUP_ORDINAL(int)

  // test CrsMatrix for some scalar
  // if explicit instantiation is enabled, test for all enabled types
  // if not, test for float
#if defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)
  #if defined(HAVE_TPETRA_INST_DOUBLE)
    UNIT_TEST_GROUP_ORDINAL_SCALAR(int,double)
  #elif defined(HAVE_TPETRA_INST_FLOAT)
    UNIT_TEST_GROUP_ORDINAL_SCALAR(int,float)
  #elif defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)
    typedef std::complex<float> ComplexFloat; \
    UNIT_TEST_GROUP_ORDINAL_SCALAR(int, ComplexFloat)
  #elif defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE)
    typedef std::complex<double> ComplexDouble; \
    UNIT_TEST_GROUP_ORDINAL_SCALAR(int, ComplexDouble)
  #endif
#else
  UNIT_TEST_GROUP_ORDINAL_SCALAR(int,float)
#endif

    TPETRA_ETI_MANGLING_TYPEDEFS()

    TPETRA_INSTANTIATE_LG( UNIT_TEST_GROUP_LO_GO )

}


