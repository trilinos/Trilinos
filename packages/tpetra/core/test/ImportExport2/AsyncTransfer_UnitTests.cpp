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
#include "Tpetra_MultiVector.hpp"
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
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include <memory>
#include <sstream>

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
  using Tpetra::DistObject;
  using Tpetra::createContigMap;
  using Tpetra::CrsGraph;
  using Tpetra::CrsMatrix;
  using Tpetra::MultiVector;
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
  // UNIT TEST FIXTURES
  //


  template <typename Scalar, typename LO, typename GO>
  class MultiVectorTransferFixture {
  private:
    using map_type = Map<LO, GO>;
    using mv_type = MultiVector<Scalar, LO, GO>;

  public:
    MultiVectorTransferFixture(Teuchos::FancyOStream& o, bool& s)
      : out(o),
        success(s),
        comm(getDefaultComm()),
        numProcs(comm->getSize()),
        myRank(comm->getRank())
    { }

    ~MultiVectorTransferFixture() { }

    bool shouldSkipTest() {
      return numProcs < 2;
    }

    void printSkippedTestMessage() {
      out << "This test is only meaningful if running with multiple MPI "
             "processes, but you ran it with only 1 process." << endl;
    }

    void setup(int collectRank) {
      setupMaps(collectRank);
      setupMultiVectors();
    }

    template <typename TransferMethod>
    void performTransfer(const TransferMethod& transfer) {
      transfer(sourceMV, targetMV);
    }

    template <typename ReferenceSolution>
    void checkResults(const ReferenceSolution& referenceSolution) {
      RCP<const mv_type> referenceMV = referenceSolution.generateWithClassicalCodePath(sourceMV, targetMap);
      compareMultiVectors(targetMV, referenceMV);
    }

  private:
    void setupMaps(int collectRank) {
      const GO indexBase = 0;
      const Tpetra::global_size_t INVALID = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();

      const size_t sourceNumLocalElements = 3;
      const size_t totalElements = numProcs*sourceNumLocalElements;
      const size_t targetNumLocalElements = (myRank == collectRank) ? totalElements : 0;

      sourceMap = rcp(new map_type(INVALID, sourceNumLocalElements, indexBase, comm));
      targetMap = rcp(new map_type(INVALID, targetNumLocalElements, indexBase, comm));
    }

    void setupMultiVectors() {
      sourceMV = rcp(new mv_type(sourceMap, 1));
      sourceMV->randomize();

      targetMV = rcp(new mv_type(targetMap, 1));
      targetMV->putScalar(Teuchos::ScalarTraits<Scalar>::zero());
    }

    void compareMultiVectors(RCP<const mv_type> resultMV, RCP<const mv_type> referenceMV) {
      auto data = resultMV->getLocalViewHost(Tpetra::Access::ReadOnly);
      auto referenceData = referenceMV->getLocalViewHost(Tpetra::Access::ReadOnly);

      for (GO globalRow = targetMap->getMinGlobalIndex(); globalRow <= targetMap->getMaxGlobalIndex(); ++globalRow) {
        const LO localRow = targetMap->getLocalElement(globalRow);
        TEST_EQUALITY(data(localRow, 0), referenceData(localRow, 0));
      }
    }

    Teuchos::FancyOStream& out;
    bool& success;

    RCP<const Comm<int>> comm;
    const int numProcs;
    const int myRank;

    RCP<const map_type> sourceMap;
    RCP<const map_type> targetMap;

    RCP<mv_type> sourceMV;
    RCP<mv_type> targetMV;
  };

  template <typename Scalar, typename LO, typename GO>
  class ReferenceImportMultiVector {
  private:
    using map_type = Map<LO, GO>;
    using mv_type = MultiVector<Scalar, LO, GO>;

  public:
    RCP<mv_type> generateWithClassicalCodePath(RCP<mv_type> sourceMV, RCP<const map_type> targetMap) const {
      RCP<const map_type> sourceMap = sourceMV->getMap();
      Import<LO, GO> importer(sourceMap, targetMap, getImportParameterList());
      expertSetRemoteLIDsContiguous(importer, false);
      TEUCHOS_ASSERT(!importer.areRemoteLIDsContiguous());

      RCP<mv_type> referenceMV = rcp(new mv_type(targetMap, 1));
      referenceMV->putScalar(Teuchos::ScalarTraits<Scalar>::zero());
      referenceMV->doImport(*sourceMV, importer, INSERT);
      TEUCHOS_ASSERT(!referenceMV->importsAreAliased());

      return referenceMV;
    }
  };


  template <typename Scalar, typename LO, typename GO>
  class DiagonalCrsMatrixTransferFixture {
  private:
    using map_type = Map<LO, GO>;
    using crs_type = CrsMatrix<Scalar, LO, GO>;

  public:
    DiagonalCrsMatrixTransferFixture(Teuchos::FancyOStream& o, bool& s)
      : out(o),
        success(s),
        comm(getDefaultComm()),
        numProcs(comm->getSize()),
        myRank(comm->getRank())
    { }

    ~DiagonalCrsMatrixTransferFixture() { }

    bool shouldSkipTest() {
      return numProcs < 2;
    }

    void printSkippedTestMessage() {
      out << "This test is only meaningful if running with multiple MPI "
             "processes, but you ran it with only 1 process." << endl;
    }

    void setup() {
      setupMaps();
      setupMatrices();
    }

    template <typename TransferMethod>
    void performTransfer(const TransferMethod& transfer) {
      transfer(sourceMat, targetMat);
      targetMat->fillComplete();
    }

    template <typename ReferenceSolution>
    void checkResults(const ReferenceSolution& referenceSolution) {
      RCP<const crs_type> referenceMat = referenceSolution.generateUsingAllInOneImport(sourceMat, targetMap);
      checkMatrixIsDiagonal(targetMat);
      checkMatrixIsDiagonal(referenceMat);
      compareMatrices(targetMat, referenceMat);
    }

  private:
    void setupMaps() {
      const GO indexBase = 0;
      const Tpetra::global_size_t INVALID = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();

      const size_t targetNumLocalElements = 3;
      const size_t totalElements = numProcs*targetNumLocalElements;
      const size_t sourceNumLocalElements = (myRank == 0) ? totalElements : 0;

      sourceMap = rcp(new map_type(INVALID, sourceNumLocalElements, indexBase, comm));
      targetMap = rcp(new map_type(INVALID, targetNumLocalElements, indexBase, comm));
    }

    void setupMatrices() {
      sourceMat = rcp(new crs_type(sourceMap, 1, StaticProfile, getCrsMatrixParameterList()));
      targetMat = rcp(new crs_type(targetMap, 1, StaticProfile, getCrsMatrixParameterList()));

      for (GO row = sourceMap->getMinGlobalIndex(); row <= sourceMap->getMaxGlobalIndex(); row++) {
        sourceMat->insertGlobalValues(row, tuple<GO>(row), tuple<Scalar>(row));
      }
      sourceMat->fillComplete();
    }

    void checkMatrixIsDiagonal(RCP<const crs_type> matrix) {
      for (GO globalRow = targetMap->getMinGlobalIndex(); globalRow <= targetMap->getMaxGlobalIndex(); ++globalRow) {
        const LO localRow = targetMap->getLocalElement(globalRow);

        typename crs_type::local_inds_host_view_type localInds;
        typename crs_type::values_host_view_type localVals;
        matrix->getLocalRowView(localRow, localInds, localVals);

        TEST_EQUALITY_CONST(localInds.size(), 1);
        if (localInds.size() == 1) {
          TEST_EQUALITY(matrix->getColMap()->getGlobalElement(localInds[0]), globalRow);
        }

        TEST_EQUALITY_CONST(localVals.size(), 1);
        if (localVals.size() == 1) {
          TEST_EQUALITY(localVals[0], as<Scalar>(globalRow));
        }
      }
    }

    void compareMatrices(RCP<const crs_type> resultMat, RCP<const crs_type> referenceMat) {
      using size_type = typename Array<Scalar>::size_type;
      using magnitude_type = typename ScalarTraits<Scalar>::magnitudeType;

      using lids_type = typename crs_type::nonconst_local_inds_host_view_type;
      using vals_type = typename crs_type::nonconst_values_host_view_type;

      const magnitude_type tol = as<magnitude_type>(10)*ScalarTraits<magnitude_type>::eps();

      lids_type resultRowIndices;
      vals_type resultRowValues;
      lids_type referenceRowIndices;
      vals_type referenceRowValues;

      for (LO localRow = targetMap->getMinLocalIndex(); localRow <= targetMap->getMaxLocalIndex(); ++localRow) {
        size_t resultNumEntries = resultMat->getNumEntriesInLocalRow(localRow);
        size_t referenceNumEntries = referenceMat->getNumEntriesInLocalRow(localRow);
        TEST_EQUALITY(resultNumEntries, referenceNumEntries);

        if (resultNumEntries > as<size_t>(resultRowIndices.size())) {
          Kokkos::resize(resultRowIndices, resultNumEntries);
          Kokkos::resize(resultRowValues, resultNumEntries);
        }
        if (referenceNumEntries > as<size_t>(referenceRowIndices.size())) {
          Kokkos::resize(referenceRowIndices, referenceNumEntries);
          Kokkos::resize(referenceRowValues, referenceNumEntries);
        }

        resultMat->getLocalRowCopy(localRow, resultRowIndices, resultRowValues, resultNumEntries);
        referenceMat->getLocalRowCopy(localRow, referenceRowIndices, referenceRowValues, referenceNumEntries);

        Tpetra::sort2(resultRowIndices, resultRowIndices.extent(0), resultRowValues);
        Tpetra::sort2(referenceRowIndices, referenceRowIndices.extent(0), referenceRowValues);

        for (size_type k = 0; k < static_cast<size_type>(resultNumEntries); ++k) {
          TEST_EQUALITY(resultRowIndices[k], referenceRowIndices[k]);
          TEST_FLOATING_EQUALITY(resultRowValues[k], referenceRowValues[k], tol);
        }
      }
    }

    Teuchos::FancyOStream& out;
    bool& success;

    RCP<const Comm<int>> comm;
    const int numProcs;
    const int myRank;

    RCP<const map_type> sourceMap;
    RCP<const map_type> targetMap;

    RCP<crs_type> sourceMat;
    RCP<crs_type> targetMat;
  };

  template <typename Scalar, typename LO, typename GO>
  class ReferenceImportMatrix {
  private:
    using map_type = Map<LO, GO>;
    using crs_type = CrsMatrix<Scalar, LO, GO>;

  public:
    RCP<crs_type> generateUsingAllInOneImport(RCP<crs_type> sourceMat, RCP<const map_type> targetMap) const {
      RCP<const map_type> sourceMap = sourceMat->getMap();
      Import<LO, GO> importer(sourceMap, targetMap, getImportParameterList());

      Teuchos::ParameterList dummy;
      RCP<crs_type> referenceMat = Tpetra::importAndFillCompleteCrsMatrix<crs_type>(
                                      sourceMat, importer, Teuchos::null, Teuchos::null, rcp(&dummy,false));
      return referenceMat;
    }
  };


  //
  // UNIT TESTS
  //


  template <typename Packet, typename LO, typename GO>
  class ForwardImport {
  private:
    using DistObjectRCP = RCP<DistObject<Packet, LO, GO>>;

  public:
    void operator()(DistObjectRCP source, DistObjectRCP target) const {
      Import<LO, GO> importer(source->getMap(), target->getMap(), getImportParameterList());
      target->doImport(*source, importer, INSERT);
    }
  };

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVectorTransfer, asyncImport, LO, GO, Scalar )
  {
    MultiVectorTransferFixture<Scalar, LO, GO> fixture(out, success);
    if (fixture.shouldSkipTest()) {
      fixture.printSkippedTestMessage();
      return;
    }

    fixture.setup(0);
    fixture.performTransfer(ForwardImport<Scalar, LO, GO>());
    fixture.checkResults(ReferenceImportMultiVector<Scalar, LO, GO>());
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( DiagonalCrsMatrixTransfer, asyncImport, LO, GO, Scalar )
  {
    DiagonalCrsMatrixTransferFixture<Scalar, LO, GO> fixture(out, success);
    if (fixture.shouldSkipTest()) {
      fixture.printSkippedTestMessage();
      return;
    }

    fixture.setup();
    fixture.performTransfer(ForwardImport<char, LO, GO>());
    fixture.checkResults(ReferenceImportMatrix<Scalar, LO, GO>());
  }


  template <typename Scalar, typename LO, typename GO>
  class CrsMatrixLowerTriangularTransferFixture {
  private:
    using map_type = Map<LO, GO>;
    using crs_type = CrsMatrix<Scalar, LO, GO>;

  public:
    CrsMatrixLowerTriangularTransferFixture(Teuchos::FancyOStream& o, bool& s)
      : out(o),
        success(s),
        comm(getDefaultComm()),
        numProcs(comm->getSize()),
        myRank(comm->getRank())
    { }

    ~CrsMatrixLowerTriangularTransferFixture() { }

    bool shouldSkipTest() {
      return numProcs%2 != 0;
    }

    void printSkippedTestMessage() {
      out << "This test is only meaningful if running with an even number of MPI processes." << endl;
    }

    void setupMatrices() {
      const GO indexBase = 0;
      const Tpetra::global_size_t INVALID = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();

      const size_t sourceNumLocalElements = (myRank%2 == 0) ? 3 : 5;
      const size_t targetNumLocalElements = 4;

      sourceMap = rcp(new map_type(INVALID, sourceNumLocalElements, indexBase, comm));
      sourceMat = rcp(new crs_type(sourceMap, 24, StaticProfile, getCrsMatrixParameterList()));

      Array<GO> cols(1);
      Array<Scalar>  vals(1);
      for (GO row = sourceMap->getMinGlobalIndex(); row <= sourceMap->getMaxGlobalIndex(); row++) {
        if (row > 0) {
          cols.resize(row);
          vals.resize(row);
          for (GO col=0; col<row; col++) {
            cols[col] = as<GO>(col);
            vals[col] = as<Scalar>(col);
          }
          sourceMat->insertGlobalValues(row, cols, vals);
        }
      }
      sourceMat->fillComplete();

      targetMap = rcp(new map_type(INVALID, targetNumLocalElements, indexBase, comm));
      targetMat = rcp(new crs_type(targetMap, 24, StaticProfile, getCrsMatrixParameterList()));
    }

    template <typename TransferMethod>
    void performTransfer(const TransferMethod& transfer) {
      transfer(sourceMat, targetMat);
      targetMat->fillComplete();
    }

    void checkResults() {
      using lids_type = typename crs_type::local_inds_host_view_type;
      using vals_type = typename crs_type::values_host_view_type;

      const RCP<const map_type> colMap = targetMat->getColMap();

      for (GO globalRow=targetMap->getMinGlobalIndex(); globalRow<=targetMap->getMaxGlobalIndex(); ++globalRow) {
        LO localRow = targetMap->getLocalElement(globalRow);
        lids_type rowIndices;
        vals_type rowValues;

        targetMat->getLocalRowView(localRow, rowIndices, rowValues);
        TEST_EQUALITY(rowIndices.extent(0), (size_t) globalRow);
        TEST_EQUALITY(rowValues.extent(0), (size_t) globalRow);

        Array<GO> indices(rowIndices.size());
        Array<Scalar> values(rowValues.size());

        for (decltype (rowIndices.size()) j=0; j<rowIndices.size(); ++j) {
          indices[j] = colMap->getGlobalElement(rowIndices[j]);
          values[j] = rowValues[j];
        }
        Tpetra::sort2(indices.begin(), indices.end(), values.begin());

        for (decltype (rowIndices.size()) j=0; j<rowIndices.size(); ++j) {
          TEST_EQUALITY(indices[j], as<GO>(j));
          TEST_EQUALITY(values[j], as<Scalar>(j));
        }
      }
    }

  private:
    Teuchos::FancyOStream& out;
    bool& success;

    RCP<const Comm<int>> comm;
    const int numProcs;
    const int myRank;

    RCP<const map_type> sourceMap;
    RCP<const map_type> targetMap;

    RCP<crs_type> sourceMat;
    RCP<crs_type> targetMat;
  };

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrixTransfer, asyncImport_lowerTriangular, LO, GO, Scalar )
  {
    CrsMatrixLowerTriangularTransferFixture<Scalar, LO, GO> fixture(out, success);
    if (fixture.shouldSkipTest()) {
      fixture.printSkippedTestMessage();
      return;
    }

    fixture.setupMatrices();
    fixture.performTransfer(ForwardImport<char, LO, GO>());
    fixture.checkResults();
  }


  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP_SC_LO_GO( SC, LO, GO )                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVectorTransfer, asyncImport, LO, GO, SC ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( DiagonalCrsMatrixTransfer, asyncImport, LO, GO, SC ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrixTransfer, asyncImport_lowerTriangular, LO, GO, SC ) \

  TPETRA_ETI_MANGLING_TYPEDEFS()

  // Test for all Scalar, LO, GO template parameter
  // combinations, and the default Node type.
  TPETRA_INSTANTIATE_SLG_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP_SC_LO_GO )
}
