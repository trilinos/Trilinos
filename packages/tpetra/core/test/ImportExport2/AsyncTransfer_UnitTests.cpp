// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_TestingUtilities.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <Teuchos_as.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_Tuple.hpp>

#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"

namespace {

  using std::endl;

  using Teuchos::as;
  using Teuchos::Array;
  using Teuchos::Comm;
  using Teuchos::FancyOStream;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::tuple;

  using Tpetra::DistObject;
  using Tpetra::CrsMatrix;
  using Tpetra::Export;
  using Tpetra::getDefaultComm;
  using Tpetra::global_size_t;
  using Tpetra::Import;
  using Tpetra::INSERT;
  using Tpetra::Map;
  using Tpetra::MultiVector;


  //
  // UNIT TEST FIXTURES
  //


  template <typename Scalar, typename LO, typename GO>
  class MultiVectorTransferFixture {
  private:
    using map_type = Map<LO, GO>;
    using mv_type = MultiVector<Scalar, LO, GO>;

  public:
    MultiVectorTransferFixture(FancyOStream& o, bool& s)
      : out(o),
        success(s),
        comm(getDefaultComm()),
        numProcs(comm->getSize()),
        myRank(comm->getRank())
    { }

    ~MultiVectorTransferFixture() { }

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
      const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();

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
      targetMV->putScalar(ScalarTraits<Scalar>::zero());
    }

    void compareMultiVectors(RCP<const mv_type> resultMV, RCP<const mv_type> referenceMV) {
      auto data = resultMV->getLocalViewHost(Tpetra::Access::ReadOnly);
      auto referenceData = referenceMV->getLocalViewHost(Tpetra::Access::ReadOnly);

      TEST_EQUALITY(data.size(), referenceData.size());
      for (LO localRow = 0; localRow < as<LO>(data.size()); localRow++) {
        TEST_EQUALITY(data(localRow, 0), referenceData(localRow, 0));
      }
    }

    FancyOStream& out;
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
      Import<LO, GO> importer(sourceMap, targetMap);
      expertSetRemoteLIDsContiguous(importer, false);
      TEUCHOS_ASSERT(!importer.areRemoteLIDsContiguous());

      RCP<mv_type> referenceMV = rcp(new mv_type(targetMap, 1));
      referenceMV->putScalar(ScalarTraits<Scalar>::zero());
      referenceMV->doImport(*sourceMV, importer, INSERT);
      TEUCHOS_ASSERT(!referenceMV->importsAreAliased());

      return referenceMV;
    }
  };

  template <typename Scalar, typename LO, typename GO>
  class ReferenceExportMultiVector {
  private:
    using map_type = Map<LO, GO>;
    using mv_type = MultiVector<Scalar, LO, GO>;

  public:
    RCP<mv_type> generateWithClassicalCodePath(RCP<mv_type> sourceMV, RCP<const map_type> targetMap) const {
      RCP<const map_type> sourceMap = sourceMV->getMap();
      Export<LO, GO> exporter(sourceMap, targetMap);
      expertSetRemoteLIDsContiguous(exporter, false);
      TEUCHOS_ASSERT(!exporter.areRemoteLIDsContiguous());

      RCP<mv_type> referenceMV = rcp(new mv_type(targetMap, 1));
      referenceMV->putScalar(ScalarTraits<Scalar>::zero());
      referenceMV->doExport(*sourceMV, exporter, INSERT);
      TEUCHOS_ASSERT(!referenceMV->importsAreAliased());

      return referenceMV;
    }
  };


  template <typename Scalar, typename LO, typename GO>
  class MultiVectorGroupTransferFixture {
  private:
    using map_type = Map<LO, GO>;
    using mv_type = MultiVector<Scalar, LO, GO>;

  public:
    MultiVectorGroupTransferFixture(FancyOStream& o, bool& s)
      : out(o),
        success(s),
        comm(getDefaultComm()),
        numProcs(comm->getSize()),
        myRank(comm->getRank()),
        numMVs(4)
    { }

    ~MultiVectorGroupTransferFixture() { }

    template <typename MapSetup>
    void setup(int collectRank) {
      setupMaps<MapSetup>(collectRank);
      setupMultiVectors();
    }

    template <typename TransferMethod>
    void performTransfer(const TransferMethod& transfer) {
      transfer(sourceMVs, targetMVs);
    }

    template <typename ReferenceSolution>
    void checkResults(const ReferenceSolution& referenceSolution) {
      for (int i=0; i<numMVs; i++) {
        RCP<const mv_type> referenceMV = referenceSolution.generateWithClassicalCodePath(sourceMVs[i], targetMap);
        compareMultiVectors(targetMVs[i], referenceMV);
      }
    }

  private:
    template <typename MapSetup>
    void setupMaps(int collectRank) {
      MapSetup maps(comm);
      maps.setup(collectRank);
      sourceMap = maps.getSourceMap();
      targetMap = maps.getTargetMap();
    }

    void setupMultiVectors() {
      for (int i=0; i<numMVs; i++) {
        sourceMVs.push_back(rcp(new mv_type(sourceMap, 1)));
        sourceMVs[i]->randomize();

        targetMVs.push_back(rcp(new mv_type(targetMap, 1)));
        targetMVs[i]->putScalar(ScalarTraits<Scalar>::zero());
      }
    }

    void compareMultiVectors(RCP<const mv_type> resultMV, RCP<const mv_type> referenceMV) {
      auto data = resultMV->getLocalViewHost(Tpetra::Access::ReadOnly);
      auto referenceData = referenceMV->getLocalViewHost(Tpetra::Access::ReadOnly);

      TEST_EQUALITY(data.size(), referenceData.size());
      for (LO localRow = 0; localRow < as<LO>(data.size()); localRow++) {
        TEST_EQUALITY(data(localRow, 0), referenceData(localRow, 0));
      }
    }

    FancyOStream& out;
    bool& success;

    RCP<const Comm<int>> comm;
    const int numProcs;
    const int myRank;

    const int numMVs;

    RCP<const map_type> sourceMap;
    RCP<const map_type> targetMap;

    std::vector<RCP<mv_type>> sourceMVs;
    std::vector<RCP<mv_type>> targetMVs;
  };

  template <typename Scalar, typename LO, typename GO>
  class DiagonalCrsMatrixTransferFixture {
  private:
    using map_type = Map<LO, GO>;
    using crs_type = CrsMatrix<Scalar, LO, GO>;

  public:
    DiagonalCrsMatrixTransferFixture(FancyOStream& o, bool& s)
      : out(o),
        success(s),
        comm(getDefaultComm()),
        numProcs(comm->getSize()),
        myRank(comm->getRank())
    { }

    ~DiagonalCrsMatrixTransferFixture() { }

    void setup() {
      setupMaps();
      setupMatrices();
    }

    template <typename TransferMethod>
    void performTransfer(const TransferMethod& transfer) {
      transfer(sourceMat, targetMat);
    }

    template <typename ReferenceSolution>
    void checkResults(const ReferenceSolution& referenceSolution) {
      targetMat->fillComplete();
      RCP<const crs_type> referenceMat = referenceSolution.generateUsingAllInOne(sourceMat, targetMap);
      checkMatrixIsDiagonal(targetMat);
      checkMatrixIsDiagonal(referenceMat);
      compareMatrices(targetMat, referenceMat);
    }

  private:
    void setupMaps() {
      const GO indexBase = 0;
      const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();

      const size_t targetNumLocalElements = 3;
      const size_t totalElements = numProcs*targetNumLocalElements;
      const size_t sourceNumLocalElements = (myRank == 0) ? totalElements : 0;

      sourceMap = rcp(new map_type(INVALID, sourceNumLocalElements, indexBase, comm));
      targetMap = rcp(new map_type(INVALID, targetNumLocalElements, indexBase, comm));
    }

    void setupMatrices() {
      sourceMat = rcp(new crs_type(sourceMap, 1));
      targetMat = rcp(new crs_type(targetMap, 1));

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

    FancyOStream& out;
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
    RCP<crs_type> generateUsingAllInOne(RCP<crs_type> sourceMat, RCP<const map_type> targetMap) const {
      RCP<const map_type> sourceMap = sourceMat->getMap();
      Import<LO, GO> importer(sourceMap, targetMap);

      Teuchos::ParameterList dummy;
      RCP<crs_type> referenceMat = Tpetra::importAndFillCompleteCrsMatrix<crs_type>(
                                      sourceMat, importer, Teuchos::null, Teuchos::null, rcp(&dummy,false));
      return referenceMat;
    }
  };

  template <typename Scalar, typename LO, typename GO>
  class ReferenceExportMatrix {
  private:
    using map_type = Map<LO, GO>;
    using crs_type = CrsMatrix<Scalar, LO, GO>;

  public:
    RCP<crs_type> generateUsingAllInOne(RCP<crs_type> sourceMat, RCP<const map_type> targetMap) const {
      RCP<const map_type> sourceMap = sourceMat->getMap();
      Export<LO, GO> exporter(sourceMap, targetMap);

      Teuchos::ParameterList dummy;
      RCP<crs_type> referenceMat = Tpetra::exportAndFillCompleteCrsMatrix<crs_type>(
                                      sourceMat, exporter, Teuchos::null, Teuchos::null, rcp(&dummy,false));
      return referenceMat;
    }
  };


  template <typename Scalar, typename LO, typename GO>
  class LowerTriangularCrsMatrixTransferFixture {
  private:
    using map_type = Map<LO, GO>;
    using crs_type = CrsMatrix<Scalar, LO, GO>;

  public:
    LowerTriangularCrsMatrixTransferFixture(FancyOStream& o, bool& s)
      : out(o),
        success(s),
        comm(getDefaultComm()),
        numProcs(comm->getSize()),
        myRank(comm->getRank())
    { }

    ~LowerTriangularCrsMatrixTransferFixture() { }

    void setup() {
      setupMaps();
      setupMatrices();
    }

    template <typename TransferMethod>
    void performTransfer(const TransferMethod& transfer) {
      transfer(sourceMat, targetMat);
    }

    void checkResults() {
      using lids_type = typename crs_type::local_inds_host_view_type;
      using vals_type = typename crs_type::values_host_view_type;

      targetMat->fillComplete();

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

        for (decltype(rowIndices.size()) j=0; j<rowIndices.size(); ++j) {
          indices[j] = colMap->getGlobalElement(rowIndices[j]);
          values[j] = rowValues[j];
        }
        Tpetra::sort2(indices.begin(), indices.end(), values.begin());

        for (decltype(rowIndices.size()) j=0; j<rowIndices.size(); ++j) {
          TEST_EQUALITY(indices[j], as<GO>(j));
          TEST_EQUALITY(values[j], as<Scalar>(j));
        }
      }
    }

  private:
    void setupMaps() {
      const GO indexBase = 0;
      const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();

      const size_t sourceNumLocalElements = isExtraRank() ? 4 : (myRank%2 == 0 ? 3 : 5);
      const size_t targetNumLocalElements = 4;

      sourceMap = rcp(new map_type(INVALID, sourceNumLocalElements, indexBase, comm));
      targetMap = rcp(new map_type(INVALID, targetNumLocalElements, indexBase, comm));
    }

    bool isExtraRank() {
      return numProcs%2 == 1 && myRank == numProcs-1;
    }

    void setupMatrices() {
      sourceMat = rcp(new crs_type(sourceMap, sourceMap->getGlobalNumElements()));
      targetMat = rcp(new crs_type(targetMap, targetMap->getGlobalNumElements()));

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
    }

    FancyOStream& out;
    bool& success;

    RCP<const Comm<int>> comm;
    const int numProcs;
    const int myRank;

    RCP<const map_type> sourceMap;
    RCP<const map_type> targetMap;

    RCP<crs_type> sourceMat;
    RCP<crs_type> targetMat;
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
      Import<LO, GO> importer(source->getMap(), target->getMap());
      target->beginImport(*source, importer, INSERT);
      target->endImport(*source, importer, INSERT);
    }
  };

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( AsyncForwardImport, MultiVector_rank0, Scalar, LO, GO )
  {
    MultiVectorTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.setup(0);
    fixture.performTransfer(ForwardImport<Scalar, LO, GO>());
    fixture.checkResults(ReferenceImportMultiVector<Scalar, LO, GO>());
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( AsyncForwardImport, MultiVector_rank1, Scalar, LO, GO )
  {
    MultiVectorTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.setup(1);
    fixture.performTransfer(ForwardImport<Scalar, LO, GO>());
    fixture.checkResults(ReferenceImportMultiVector<Scalar, LO, GO>());
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( AsyncForwardImport, DiagonalCrsMatrix, Scalar, LO, GO )
  {
    DiagonalCrsMatrixTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.setup();
    fixture.performTransfer(ForwardImport<char, LO, GO>());
    fixture.checkResults(ReferenceImportMatrix<Scalar, LO, GO>());
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( AsyncForwardImport, LowerTriangularCrsMatrix, Scalar, LO, GO )
  {
    LowerTriangularCrsMatrixTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.setup();
    fixture.performTransfer(ForwardImport<char, LO, GO>());
    fixture.checkResults();
  }


  template <typename Packet, typename LO, typename GO>
  class ReverseImport {
  private:
    using DistObjectRCP = RCP<DistObject<Packet, LO, GO>>;

  public:
    void operator()(DistObjectRCP source, DistObjectRCP target) const {
      Export<LO, GO> exporter(target->getMap(), source->getMap());
      target->beginImport(*source, exporter, INSERT);
      target->endImport(*source, exporter, INSERT);
    }
  };

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( AsyncReverseImport, MultiVector_rank0, Scalar, LO, GO )
  {
    MultiVectorTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.setup(0);
    fixture.performTransfer(ReverseImport<Scalar, LO, GO>());
    fixture.checkResults(ReferenceImportMultiVector<Scalar, LO, GO>());
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( AsyncReverseImport, MultiVector_rank1, Scalar, LO, GO )
  {
    MultiVectorTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.setup(1);
    fixture.performTransfer(ReverseImport<Scalar, LO, GO>());
    fixture.checkResults(ReferenceImportMultiVector<Scalar, LO, GO>());
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( AsyncReverseImport, DiagonalCrsMatrix, Scalar, LO, GO )
  {
    DiagonalCrsMatrixTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.setup();
    fixture.performTransfer(ReverseImport<char, LO, GO>());
    fixture.checkResults(ReferenceImportMatrix<Scalar, LO, GO>());
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( AsyncReverseImport, LowerTriangularCrsMatrix, Scalar, LO, GO )
  {
    LowerTriangularCrsMatrixTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.setup();
    fixture.performTransfer(ReverseImport<char, LO, GO>());
    fixture.checkResults();
  }


  template <typename Packet, typename LO, typename GO>
  class ForwardExport {
  private:
    using DistObjectRCP = RCP<DistObject<Packet, LO, GO>>;

  public:
    void operator()(DistObjectRCP source, DistObjectRCP target) const {
      Export<LO, GO> exporter(source->getMap(), target->getMap());
      target->beginExport(*source, exporter, INSERT);
      target->endExport(*source, exporter, INSERT);
    }
  };

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( AsyncForwardExport, MultiVector_rank0, Scalar, LO, GO )
  {
    MultiVectorTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.setup(0);
    fixture.performTransfer(ForwardExport<Scalar, LO, GO>());
    fixture.checkResults(ReferenceExportMultiVector<Scalar, LO, GO>());
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( AsyncForwardExport, MultiVector_rank1, Scalar, LO, GO )
  {
    MultiVectorTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.setup(1);
    fixture.performTransfer(ForwardExport<Scalar, LO, GO>());
    fixture.checkResults(ReferenceExportMultiVector<Scalar, LO, GO>());
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( AsyncForwardExport, DiagonalCrsMatrix, Scalar, LO, GO )
  {
    DiagonalCrsMatrixTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.setup();
    fixture.performTransfer(ForwardExport<char, LO, GO>());
    fixture.checkResults(ReferenceExportMatrix<Scalar, LO, GO>());
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( AsyncForwardExport, LowerTriangularCrsMatrix, Scalar, LO, GO )
  {
    LowerTriangularCrsMatrixTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.setup();
    fixture.performTransfer(ForwardExport<char, LO, GO>());
    fixture.checkResults();
  }


  template <typename Packet, typename LO, typename GO>
  class ReverseExport {
  private:
    using DistObjectRCP = RCP<DistObject<Packet, LO, GO>>;

  public:
    void operator()(DistObjectRCP source, DistObjectRCP target) const {
      Import<LO, GO> importer(target->getMap(), source->getMap());
      target->beginExport(*source, importer, INSERT);
      target->endExport(*source, importer, INSERT);
    }
  };

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( AsyncReverseExport, MultiVector_rank0, Scalar, LO, GO )
  {
    MultiVectorTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.setup(0);
    fixture.performTransfer(ReverseExport<Scalar, LO, GO>());
    fixture.checkResults(ReferenceExportMultiVector<Scalar, LO, GO>());
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( AsyncReverseExport, MultiVector_rank1, Scalar, LO, GO )
  {
    MultiVectorTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.setup(1);
    fixture.performTransfer(ReverseExport<Scalar, LO, GO>());
    fixture.checkResults(ReferenceExportMultiVector<Scalar, LO, GO>());
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( AsyncReverseExport, DiagonalCrsMatrix, Scalar, LO, GO )
  {
    DiagonalCrsMatrixTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.setup();
    fixture.performTransfer(ReverseExport<char, LO, GO>());
    fixture.checkResults(ReferenceExportMatrix<Scalar, LO, GO>());
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( AsyncReverseExport, LowerTriangularCrsMatrix, Scalar, LO, GO )
  {
    LowerTriangularCrsMatrixTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.setup();
    fixture.performTransfer(ReverseExport<char, LO, GO>());
    fixture.checkResults();
  }


  template <typename Packet, typename LO, typename GO>
  class CheckDefaultTransferStatus {
  private:
    using DistObjectRCP = RCP<DistObject<Packet, LO, GO>>;

  public:
    CheckDefaultTransferStatus(FancyOStream& o, bool& s)
      : out(o),
        success(s)
    { }

    void operator()(DistObjectRCP source, DistObjectRCP target) const {
      TEST_ASSERT(target->transferArrived());
    }

  private:
    FancyOStream& out;
    bool& success;
  };

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( TransferArrived, MultiVector_defaultTrue, Scalar, LO, GO )
  {
    MultiVectorTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.setup(0);
    fixture.performTransfer(CheckDefaultTransferStatus<Scalar, LO, GO>(out, success));
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( TransferArrived, CrsMatrix_defaultTrue, Scalar, LO, GO )
  {
    DiagonalCrsMatrixTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.setup();
    fixture.performTransfer(CheckDefaultTransferStatus<char, LO, GO>(out, success));
  }


  template <typename Packet, typename LO, typename GO>
  class CheckTransferArrivedForwardImport {
  private:
    using DistObjectRCP = RCP<DistObject<Packet, LO, GO>>;

  public:
    CheckTransferArrivedForwardImport(FancyOStream& o, bool& s)
      : out(o),
        success(s)
    { }

    void operator()(DistObjectRCP source, DistObjectRCP target) const {
      Import<LO, GO> importer(source->getMap(), target->getMap());

      target->beginImport(*source, importer, INSERT);
      while (!target->transferArrived());
      target->endImport(*source, importer, INSERT);

      TEST_ASSERT(target->transferArrived());
    }

  private:
    FancyOStream& out;
    bool& success;
  };

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( TransferArrived, MultiVector_forwardImportTrue, Scalar, LO, GO )
  {
    MultiVectorTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.setup(0);
    fixture.performTransfer(CheckTransferArrivedForwardImport<Scalar, LO, GO>(out, success));
    fixture.checkResults(ReferenceImportMultiVector<Scalar, LO, GO>());
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( TransferArrived, CrsMatrix_forwardImportTrue, Scalar, LO, GO )
  {
    DiagonalCrsMatrixTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.setup();
    fixture.performTransfer(CheckTransferArrivedForwardImport<char, LO, GO>(out, success));
    fixture.checkResults(ReferenceImportMatrix<Scalar, LO, GO>());
  }


  template <typename Packet, typename LO, typename GO>
  class CheckTransferArrivedForwardExport {
  private:
    using DistObjectRCP = RCP<DistObject<Packet, LO, GO>>;

  public:
    CheckTransferArrivedForwardExport(FancyOStream& o, bool& s)
      : out(o),
        success(s)
    { }

    void operator()(DistObjectRCP source, DistObjectRCP target) const {
      Export<LO, GO> exporter(source->getMap(), target->getMap());

      target->beginExport(*source, exporter, INSERT);
      while (!target->transferArrived());
      target->endExport(*source, exporter, INSERT);

      TEST_ASSERT(target->transferArrived());
    }

  private:
    FancyOStream& out;
    bool& success;
  };

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( TransferArrived, MultiVector_forwardExportTrue, Scalar, LO, GO )
  {
    MultiVectorTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.setup(0);
    fixture.performTransfer(CheckTransferArrivedForwardExport<Scalar, LO, GO>(out, success));
    fixture.checkResults(ReferenceExportMultiVector<Scalar, LO, GO>());
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( TransferArrived, CrsMatrix_forwardExportTrue, Scalar, LO, GO )
  {
    DiagonalCrsMatrixTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.setup();
    fixture.performTransfer(CheckTransferArrivedForwardExport<char, LO, GO>(out, success));
    fixture.checkResults(ReferenceExportMatrix<Scalar, LO, GO>());
  }


  template <typename Packet, typename LO, typename GO>
  class CheckTransferNotArrivedForwardImport {
  private:
    using DistObjectRCP = RCP<DistObject<Packet, LO, GO>>;

  public:
    CheckTransferNotArrivedForwardImport(FancyOStream& o, bool& s)
      : out(o),
        success(s)
    { }

    void operator()(DistObjectRCP source, DistObjectRCP target) const {
      auto comm = target->getMap()->getComm();
      int myRank = comm->getRank();
      int numProcs = comm->getSize();

      Import<LO, GO> importer(source->getMap(), target->getMap());

      if (numProcs > 1) {
        if (myRank == 0) {
          target->beginImport(*source, importer, INSERT);
          TEST_ASSERT(!target->transferArrived());
        }
        comm->barrier();
        if (myRank != 0) {
          target->beginImport(*source, importer, INSERT);
        }
        target->endImport(*source, importer, INSERT);
      }
    }

  private:
    FancyOStream& out;
    bool& success;
  };

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( TransferArrived, MultiVector_forwardImportFalse, Scalar, LO, GO )
  {
    MultiVectorTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.setup(0);
    fixture.performTransfer(CheckTransferNotArrivedForwardImport<Scalar, LO, GO>(out, success));
  }


  template <typename Packet, typename LO, typename GO>
  class CheckTransferNotArrivedForwardExport {
  private:
    using DistObjectRCP = RCP<DistObject<Packet, LO, GO>>;

  public:
    CheckTransferNotArrivedForwardExport(FancyOStream& o, bool& s)
      : out(o),
        success(s)
    { }

    void operator()(DistObjectRCP source, DistObjectRCP target) const {
      auto comm = target->getMap()->getComm();
      int myRank = comm->getRank();
      int numProcs = comm->getSize();

      Export<LO, GO> exporter(source->getMap(), target->getMap());

      if (numProcs > 1) {
        if (myRank == 0) {
          target->beginExport(*source, exporter, INSERT);
          TEST_ASSERT(!target->transferArrived());
        }
        comm->barrier();
        if (myRank != 0) {
          target->beginExport(*source, exporter, INSERT);
        }
        target->endExport(*source, exporter, INSERT);
      }
    }

  private:
    FancyOStream& out;
    bool& success;
  };

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( TransferArrived, MultiVector_forwardExportFalse, Scalar, LO, GO )
  {
    MultiVectorTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.setup(0);
    fixture.performTransfer(CheckTransferNotArrivedForwardExport<Scalar, LO, GO>(out, success));
  }


  template <typename Packet, typename LO, typename GO>
  class ForwardImportGroup {
  private:
    using DistObjectRCP = RCP<MultiVector<Packet, LO, GO>>;

  public:
    void operator()(std::vector<DistObjectRCP>& sources, std::vector<DistObjectRCP>& targets) const {
      Import<LO, GO> importer(sources[0]->getMap(), targets[0]->getMap());

      for (unsigned i=0; i<sources.size(); i++) {
        targets[i]->beginImport(*sources[i], importer, INSERT);
      }

      unsigned completedImports = 0;
      std::vector<bool> completedImport(sources.size(), false);
      while (completedImports < completedImport.size()) {
        for (unsigned i=0; i<sources.size(); i++) {
          if (completedImport[i]) continue;
          if (targets[i]->transferArrived()) {
            targets[i]->endImport(*sources[i], importer, INSERT);
            completedImport[i] = true;
            completedImports++;
          }
        }
      }
    }
  };

  template <typename LO, typename GO>
  class ContiguousMaps {
  private:
    using map_type = Map<LO, GO>;

  public:
    ContiguousMaps(RCP<const Comm<int>> c)
      : comm(c),
        numProcs(comm->getSize()),
        myRank(comm->getRank())
    { }

    void setup(int collectRank) {
      const GO indexBase = 0;
      const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();

      const size_t sourceNumLocalElements = 3;
      const size_t totalElements = numProcs*sourceNumLocalElements;
      const size_t targetNumLocalElements = (myRank == collectRank) ? totalElements : 0;

      sourceMap = rcp(new map_type(INVALID, sourceNumLocalElements, indexBase, comm));
      targetMap = rcp(new map_type(INVALID, targetNumLocalElements, indexBase, comm));
    }

    RCP<const map_type> getSourceMap() { return sourceMap; }
    RCP<const map_type> getTargetMap() { return targetMap; }

  private:
    RCP<const Comm<int>> comm;
    const int numProcs;
    const int myRank;

    RCP<const map_type> sourceMap;
    RCP<const map_type> targetMap;
  };

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( AsyncForwardImport, MultiVectorGroup_ContiguousMaps_rank0, Scalar, LO, GO )
  {
    MultiVectorGroupTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.template setup<ContiguousMaps<LO, GO>>(0);
    fixture.performTransfer(ForwardImportGroup<Scalar, LO, GO>());
    fixture.checkResults(ReferenceImportMultiVector<Scalar, LO, GO>());
  }

  template <typename LO, typename GO>
  class CyclicMaps {
  private:
    using map_type = Map<LO, GO>;

  public:
    CyclicMaps(RCP<const Comm<int>> c)
      : comm(c),
        numProcs(comm->getSize()),
        myRank(comm->getRank())
    { }

    void setup(int collectRank) {
      const GO indexBase = 0;
      const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();

      const size_t sourceNumLocalElements = 3;
      const size_t totalElements = numProcs*sourceNumLocalElements;
      const size_t targetNumLocalElements = (myRank == collectRank) ? totalElements : 0;

      Teuchos::Array<GO> sourceEntries(sourceNumLocalElements);
      for (size_t i=0; i<sourceNumLocalElements; i++) {
        sourceEntries[i] = i*numProcs + myRank;
      }

      sourceMap = rcp(new map_type(INVALID, sourceEntries, indexBase, comm));
      targetMap = rcp(new map_type(INVALID, targetNumLocalElements, indexBase, comm));
    }

    RCP<const map_type> getSourceMap() { return sourceMap; }
    RCP<const map_type> getTargetMap() { return targetMap; }

  private:
    RCP<const Comm<int>> comm;
    const int numProcs;
    const int myRank;

    RCP<const map_type> sourceMap;
    RCP<const map_type> targetMap;
  };

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( AsyncForwardImport, MultiVectorGroup_CyclicMaps_rank0, Scalar, LO, GO )
  {
    MultiVectorGroupTransferFixture<Scalar, LO, GO> fixture(out, success);

    fixture.template setup<CyclicMaps<LO, GO>>(0);
    fixture.performTransfer(ForwardImportGroup<Scalar, LO, GO>());
    fixture.checkResults(ReferenceImportMultiVector<Scalar, LO, GO>());
  }

  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP_SC_LO_GO( SC, LO, GO )                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( AsyncForwardImport, MultiVector_rank0, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( AsyncForwardImport, MultiVector_rank1, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( AsyncForwardImport, DiagonalCrsMatrix, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( AsyncForwardImport, LowerTriangularCrsMatrix, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( AsyncReverseImport, MultiVector_rank0, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( AsyncReverseImport, MultiVector_rank1, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( AsyncReverseImport, DiagonalCrsMatrix, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( AsyncReverseImport, LowerTriangularCrsMatrix, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( AsyncForwardExport, MultiVector_rank0, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( AsyncForwardExport, MultiVector_rank1, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( AsyncForwardExport, DiagonalCrsMatrix, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( AsyncForwardExport, LowerTriangularCrsMatrix, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( AsyncReverseExport, MultiVector_rank0, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( AsyncReverseExport, MultiVector_rank1, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( AsyncReverseExport, DiagonalCrsMatrix, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( AsyncReverseExport, LowerTriangularCrsMatrix, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( TransferArrived, MultiVector_defaultTrue, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( TransferArrived, CrsMatrix_defaultTrue, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( TransferArrived, MultiVector_forwardImportTrue, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( TransferArrived, CrsMatrix_forwardImportTrue, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( TransferArrived, MultiVector_forwardExportTrue, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( TransferArrived, CrsMatrix_forwardExportTrue, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( TransferArrived, MultiVector_forwardImportFalse, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( TransferArrived, MultiVector_forwardExportFalse, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( AsyncForwardImport, MultiVectorGroup_ContiguousMaps_rank0, SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( AsyncForwardImport, MultiVectorGroup_CyclicMaps_rank0, SC, LO, GO ) \

  TPETRA_ETI_MANGLING_TYPEDEFS()

  // Test for all Scalar, LO, GO template parameter
  // combinations, and the default Node type.
  TPETRA_INSTANTIATE_SLG_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP_SC_LO_GO )
}
