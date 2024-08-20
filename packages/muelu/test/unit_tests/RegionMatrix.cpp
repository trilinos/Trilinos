// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Galeri_XpetraParameters.hpp>

// Region MG headers
#include "SetupRegionUtilities.hpp"
#include "SetupRegionVector_def.hpp"
#include "SetupRegionMatrix_def.hpp"

namespace MueLuTests {

// createRegionMatrix is a helper function that allows us to easily
// generate a region matrix based on the corresponding composite
// matrix.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void createRegionMatrix(const Teuchos::ParameterList galeriList,
                        const int numDofsPerNode,
                        const RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > nodeMap,
                        const RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > dofMap,
                        const RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A,
                        RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >& rowMap,
                        RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >& colMap,
                        RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >& revisedRowMap,
                        RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >& revisedColMap,
                        RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> >& rowImport,
                        RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> >& colImport,
                        RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& regionMats,
                        Teuchos::ArrayRCP<LocalOrdinal>& regionMatVecLIDs,
                        Teuchos::RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> >& regionInterfaceImporter) {
#include <MueLu_UseShortNames.hpp>

  std::string matrixType = galeriList.get<std::string>("matrixType");
  int numDimensions      = 0;
  if (matrixType == "Laplace1D") {
    numDimensions = 1;
  } else if (matrixType == "Laplace2D" || matrixType == "Elasticity2D" ||
             matrixType == "BigStar2D" || matrixType == "Elasticity2D") {
    numDimensions = 2;
  } else if (matrixType == "Laplace3D" || matrixType == "Brick3D" || matrixType == "Elasticity3D") {
    numDimensions = 3;
  }

  // Set global geometric data
  Array<LO> lNodesPerDir(3);
  Array<GO> gNodesPerDir(3);
  Array<GO> procsPerDim(3);
  gNodesPerDir[0] = galeriList.get<GO>("nx");
  gNodesPerDir[1] = (1 < numDimensions) ? galeriList.get<GO>("ny") : 1;
  gNodesPerDir[2] = (2 < numDimensions) ? galeriList.get<GO>("nz") : 1;
  lNodesPerDir[0] = galeriList.get<LO>("lnx");
  lNodesPerDir[1] = (1 < numDimensions) ? galeriList.get<LO>("lny") : 1;
  lNodesPerDir[2] = (2 < numDimensions) ? galeriList.get<LO>("lnz") : 1;
  procsPerDim[0]  = galeriList.get<GO>("mx");
  procsPerDim[1]  = galeriList.get<GO>("my");
  procsPerDim[2]  = galeriList.get<GO>("mz");

  // std::cout << "p=" << nodeMap->getComm()->getRank() << " | numDimensions=" << numDimensions
  //           << ", useStructured=" << false << ", numDofsPerNode=" << numDofsPerNode
  //           << ", gNodesPerDir=" << gNodesPerDir << ", lNodesPerDir=" << lNodesPerDir
  //           << ", procsPerDim=" << procsPerDim << std::endl;

  Array<int> boundaryConditions;
  int maxRegPerGID       = 0;
  int numInterfaces      = 0;
  LO numLocalRegionNodes = 0;
  Array<GO> sendGIDs;
  Array<int> sendPIDs;
  Array<LO> rNodesPerDim(3);
  Array<LO> compositeToRegionLIDs(nodeMap->getLocalNumElements() * numDofsPerNode);
  Array<GO> quasiRegionGIDs;
  Array<GO> quasiRegionCoordGIDs;
  Array<GO> interfaceCompositeGIDs, interfaceRegionGIDs;
  Array<LO> interfaceRegionLIDs;
  createRegionData(numDimensions, false, numDofsPerNode,
                   gNodesPerDir(), lNodesPerDir(), procsPerDim(),
                   nodeMap, dofMap,
                   maxRegPerGID, numLocalRegionNodes, boundaryConditions,
                   sendGIDs, sendPIDs, numInterfaces, rNodesPerDim,
                   quasiRegionGIDs, quasiRegionCoordGIDs, compositeToRegionLIDs,
                   interfaceCompositeGIDs, interfaceRegionLIDs);

  // const int myRank = A->getRowMap()->getComm()->getRank();
  // const LO numSend = static_cast<LO>(sendGIDs.size());
  // std::cout << "p=" << myRank << " | numSend=" << numSend << std::endl;
  // std::cout << "p=" << myRank << " | sendGIDs: " << sendGIDs << std::endl;
  // std::cout << "p=" << myRank << " | sendPIDs: " << sendPIDs << std::endl;
  // std::cout << "p=" << myRank << " | compositeToRegionLIDs: " << compositeToRegionLIDs << std::endl;
  // std::cout << "p=" << myRank << " | quasiRegionGIDs: " << quasiRegionGIDs << std::endl;
  // std::cout << "p=" << myRank << " | interfaceCompositeGIDs" << interfaceCompositeGIDs << std::endl;
  // std::cout << "p=" << myRank << " | interfaceRegionLIDs" << interfaceRegionLIDs() << std::endl;

  rowMap = Xpetra::MapFactory<LO, GO, Node>::Build(A->getRowMap()->lib(),
                                                   Teuchos::OrdinalTraits<GO>::invalid(),
                                                   quasiRegionGIDs(),
                                                   A->getRowMap()->getIndexBase(),
                                                   A->getRowMap()->getComm());
  colMap = rowMap;

  revisedRowMap = Xpetra::MapFactory<LO, GO, Node>::Build(A->getRowMap()->lib(),
                                                          Teuchos::OrdinalTraits<GO>::invalid(),
                                                          quasiRegionGIDs.size(),
                                                          A->getRowMap()->getIndexBase(),
                                                          A->getRowMap()->getComm());
  revisedColMap = revisedRowMap;

  ExtractListOfInterfaceRegionGIDs(revisedRowMap, interfaceRegionLIDs, interfaceRegionGIDs);

  rowImport = ImportFactory::Build(dofMap, rowMap);
  colImport = ImportFactory::Build(dofMap, colMap);

  RCP<Xpetra::MultiVector<LO, LO, GO, NO> > regionsPerGIDWithGhosts;
  RCP<Xpetra::MultiVector<GO, LO, GO, NO> > interfaceGIDsMV;
  MakeRegionPerGIDWithGhosts(nodeMap, revisedRowMap, rowImport,
                             maxRegPerGID, numDofsPerNode,
                             lNodesPerDir, sendGIDs, sendPIDs, interfaceRegionLIDs,
                             regionsPerGIDWithGhosts, interfaceGIDsMV);

  SetupMatVec(interfaceGIDsMV, regionsPerGIDWithGhosts, revisedRowMap, rowImport,
              regionMatVecLIDs, regionInterfaceImporter);

  RCP<Matrix> quasiRegionGrpMats = Teuchos::null;
  MakeQuasiregionMatrices(Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A),
                          regionsPerGIDWithGhosts, rowMap, colMap, rowImport,
                          quasiRegionGrpMats, regionMatVecLIDs);

  MakeRegionMatrices(Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A), A->getRowMap(), rowMap,
                     revisedRowMap, revisedColMap,
                     rowImport, quasiRegionGrpMats, regionMats);

}  // createRegionMatrix

// Helper function that creates almost all the data needed to generate a unit-test
// numDofsPerNode [in]: number of degrees of freedom per grid point
// galeriParameters [in]: parameters passed to galeri to generate the composite problem
// comm [in]: the MPI communicator used with distributed objects
// A [out]: composite matrix
// regionMats [out]: the region matrix
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void createProblem(const LocalOrdinal numDofsPerNode,
                   Galeri::Xpetra::Parameters<GlobalOrdinal>& galeriParameters,
                   RCP<const Teuchos::Comm<int> > comm,
                   RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& A,
                   RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& regionMats,
                   RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >& revisedRowMap,
                   RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> >& rowImport,
                   Teuchos::ArrayRCP<LocalOrdinal>& regionMatVecLIDs,
                   RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> >& regionInterfaceImporter) {
#include <MueLu_UseShortNames.hpp>
  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;

  Teuchos::ParameterList galeriList = galeriParameters.GetParameterList();
  std::string matrixType            = galeriParameters.GetMatrixType();
  std::string mapType, coordinatesType;
  if ((matrixType == "Laplace2D") || (matrixType == "Elasticity2D")) {
    mapType         = "Cartesian2D";
    coordinatesType = "2D";
  } else if ((matrixType == "Laplace3D") || (matrixType == "Elasticity3D")) {
    mapType         = "Cartesian3D";
    coordinatesType = "3D";
  }

  // Build maps for the problem
  RCP<Map> nodeMap = Galeri::Xpetra::CreateMap<LO, GO, Node>(TestHelpers::Parameters::getLib(),
                                                             mapType, comm, galeriList);
  RCP<Map> dofMap  = Xpetra::MapFactory<LO, GO, Node>::Build(nodeMap, numDofsPerNode);

  // Build the Xpetra problem
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(galeriParameters.GetMatrixType(), dofMap, galeriList);

  // Generate the operator
  A = Pr->BuildMatrix();
  A->SetFixedBlockSize(numDofsPerNode);

  // Create auxiliary data for MG
  RCP<MultiVector> nullspace             = Pr->BuildNullspace();
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<double, LO, GO, Map, RealValuedMultiVector>(coordinatesType, nodeMap, galeriList);

  // create the region maps, importer and operator from composite counter parts
  RCP<const Map> rowMap        = Teuchos::null;
  RCP<const Map> colMap        = Teuchos::null;
  RCP<const Map> revisedColMap = Teuchos::null;
  RCP<Import> colImport        = Teuchos::null;
  createRegionMatrix(galeriList, numDofsPerNode, nodeMap, dofMap, A,
                     rowMap, colMap, revisedRowMap, revisedColMap,
                     rowImport, colImport, regionMats,
                     regionMatVecLIDs, regionInterfaceImporter);

  // Debug output
  // std::cout << "p=" << comm->getRank() << " | regionMatVecLIDs: " << regionMatVecLIDs << std::endl;
  // std::cout << "p=" << comm->getRank() << " | source map element list: "
  //           << regionInterfaceImporter->getSourceMap()->getLocalElementList() << std::endl;
  // std::cout << "p=" << comm->getRank() << " | target map element list: "
  //           << regionInterfaceImporter->getTargetMap()->getLocalElementList() << std::endl;

}  // createProblem

// test_matrix() is checking that performing a MatVec with composite A and region A
// yields the same vector. It also verifies that regionalToComposite(regA) returns
// the same matrix as composite A. It is a convenience function to perform common tests.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void test_matrix(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A,
                 RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regionMats,
                 RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > rowMap,
                 RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > colMap,
                 RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > revisedRowMap,
                 RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > rowImport,
                 Teuchos::FancyOStream& out,
                 bool& success) {
#include <MueLu_UseShortNames.hpp>
  using TST                           = Teuchos::ScalarTraits<SC>;
  using magnitude_type                = typename TST::magnitudeType;
  using TMT                           = Teuchos::ScalarTraits<magnitude_type>;
  RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();

  /************************************/
  /*                                  */
  /*  First test: MatVec comparison   */
  /*                                  */
  /************************************/
  RCP<Vector> X = VectorFactory::Build(A->getRowMap());
  RCP<Vector> B = VectorFactory::Build(A->getRowMap());

  // Generate a random composite X vector
  // we set seed for reproducibility
  Utilities::SetRandomSeed(*comm);
  X->randomize();

  // Perform composite MatVec
  A->apply(*X, *B, Teuchos::NO_TRANS, TST::one(), TST::zero());

  // Now build the region X vector
  RCP<Vector> quasiRegX = Teuchos::null;
  RCP<Vector> quasiRegB = Teuchos::null;
  RCP<Vector> regX      = Teuchos::null;
  RCP<Vector> regB      = Teuchos::null;
  compositeToRegional(X, quasiRegX, regX,
                      revisedRowMap, rowImport);
  regB = VectorFactory::Build(revisedRowMap, true);

  // Perform regional MatVec
  regionMats->apply(*regX, *regB, Teuchos::NO_TRANS, TST::one(), TST::zero());

  // Bring the result of the region MatVec
  // to composite format so it can be compared
  // with the original composite B vector.
  RCP<Vector> compB = VectorFactory::Build(A->getRowMap());
  regionalToComposite(regB, compB, rowImport);

  // Extract the data from B and compB to compare it
  ArrayRCP<const SC> dataB     = B->getData(0);
  ArrayRCP<const SC> dataCompB = compB->getData(0);
  for (size_t idx = 0; idx < B->getLocalLength(); ++idx) {
    TEST_FLOATING_EQUALITY(TST::magnitude(dataB[idx]),
                           TST::magnitude(dataCompB[idx]),
                           100 * TMT::eps());
  }

  /************************************/
  /*                                  */
  /*  Second test: composite matrix   */
  /*           calculation            */
  /*                                  */
  /************************************/
  RCP<Matrix> compositeMatrix = MatrixFactory::Build(A->getRowMap(), 10);
  // Transform region A into composite A.
  regionalToComposite(regionMats,
                      rowMap, colMap,
                      rowImport, Xpetra::INSERT,
                      compositeMatrix);

  // Extract the local data from the original and final matrices to compare them
  using local_matrix_type = typename Matrix::local_matrix_type;
  using local_graph_type  = typename local_matrix_type::staticcrsgraph_type;
  using entries_type      = typename local_graph_type::entries_type;
  using values_type       = typename local_matrix_type::values_type;

  local_matrix_type orignalA = A->getLocalMatrixDevice();  // Local matrix
  entries_type refEntries    = orignalA.graph.entries;     // view of local column indices
  values_type refValues      = orignalA.values;            // view of local values

  typename entries_type::HostMirror refEntries_h = Kokkos::create_mirror_view(refEntries);
  Kokkos::deep_copy(refEntries_h, refEntries);
  typename values_type::HostMirror refValues_h = Kokkos::create_mirror_view(refValues);
  Kokkos::deep_copy(refValues_h, refValues);

  local_matrix_type compositeA  = compositeMatrix->getLocalMatrixDevice();  // Local matrix
  entries_type compositeEntries = compositeA.graph.entries;                 // view of local column indices
  values_type compositeValues   = compositeA.values;                        // view of local values

  typename entries_type::HostMirror compositeEntries_h = Kokkos::create_mirror_view(compositeEntries);
  Kokkos::deep_copy(compositeEntries_h, compositeEntries);
  typename values_type::HostMirror compositeValues_h = Kokkos::create_mirror_view(compositeValues);
  Kokkos::deep_copy(compositeValues_h, compositeValues);

  TEST_EQUALITY(compositeEntries_h.extent(0), refEntries_h.extent(0));
  TEST_EQUALITY(compositeValues_h.extent(0), refValues_h.extent(0));
  for (LO idx = 0; idx < compositeEntries_h.extent_int(0); ++idx) {
    TEST_EQUALITY(compositeEntries_h(idx), refEntries_h(idx));
    TEST_FLOATING_EQUALITY(TST::magnitude(compositeValues_h(idx)),
                           TST::magnitude(refValues_h(idx)),
                           100 * TMT::eps());
  }
}  // test_matrix

// This test only checks the compositeToRegion() operation
// for matrices. More specifically we compute a region A
// based on a composite A and check the values in region A
// against know correct values.
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RegionMatrix, CompositeToRegionMatrix, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;

  // Get MPI parameter
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  LO numRanks                         = comm->getSize();
  LO myRank                           = comm->getRank();

  GO nx = 5, ny = 5, nz = 1;
  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP();
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D");
  Teuchos::ParameterList galeriList = galeriParameters.GetParameterList();
  std::string matrixType            = galeriParameters.GetMatrixType();

  // Build maps for the problem
  const LO numDofsPerNode = 1;
  RCP<Map> nodeMap        = Galeri::Xpetra::CreateMap<LO, GO, Node>(TestHelpers::Parameters::getLib(),
                                                             "Cartesian2D", comm, galeriList);
  RCP<Map> dofMap         = Xpetra::MapFactory<LO, GO, Node>::Build(nodeMap, numDofsPerNode);

  // Build the Xpetra problem
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(galeriParameters.GetMatrixType(), dofMap, galeriList);

  // Generate the operator
  RCP<Matrix> A = Pr->BuildMatrix();
  A->SetFixedBlockSize(numDofsPerNode);

  // Create the left and right handside vectors
  RCP<Vector> X = VectorFactory::Build(dofMap);
  RCP<Vector> B = VectorFactory::Build(dofMap);

  // Create auxiliary data for MG
  RCP<MultiVector> nullspace             = Pr->BuildNullspace();
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<double, LO, GO, Map, RealValuedMultiVector>("2D", nodeMap, galeriList);

  // Create the region version of A called regionMats
  RCP<const Map> rowMap        = Teuchos::null;
  RCP<const Map> colMap        = Teuchos::null;
  RCP<const Map> revisedRowMap = Teuchos::null;
  RCP<const Map> revisedColMap = Teuchos::null;
  RCP<Import> rowImport        = Teuchos::null;
  RCP<Import> colImport        = Teuchos::null;
  RCP<Matrix> regionMats       = Teuchos::null;
  Teuchos::ArrayRCP<LO> regionMatVecLIDs;
  RCP<Import> regionInterfaceImporter;
  createRegionMatrix(galeriList, numDofsPerNode, nodeMap, dofMap, A,
                     rowMap, colMap, revisedRowMap, revisedColMap,
                     rowImport, colImport, regionMats,
                     regionMatVecLIDs, regionInterfaceImporter);

  // Extract the local data from the region matrix
  using local_matrix_type = typename Matrix::local_matrix_type;
  using local_graph_type  = typename local_matrix_type::staticcrsgraph_type;
  using entries_type      = typename local_graph_type::entries_type;
  using values_type       = typename local_matrix_type::values_type;

  local_matrix_type myLocalA = regionMats->getLocalMatrixDevice();  // Local matrix
  entries_type myEntries     = myLocalA.graph.entries;              // view of local column indices
  values_type myValues       = myLocalA.values;                     // view of local values

  typename entries_type::HostMirror myEntries_h = Kokkos::create_mirror_view(myEntries);
  Kokkos::deep_copy(myEntries_h, myEntries);
  typename values_type::HostMirror myValues_h = Kokkos::create_mirror_view(myValues);
  Kokkos::deep_copy(myValues_h, myValues);

  // Now do a bunch of checks regarding the values stored in region A
  if (numRanks == 1) {
    TEST_EQUALITY(regionMats->getGlobalNumRows(), 25);
    TEST_EQUALITY(regionMats->getGlobalNumCols(), 25);
    TEST_EQUALITY(regionMats->getLocalNumRows(), 25);
    TEST_EQUALITY(regionMats->getGlobalNumEntries(), 105);
    TEST_EQUALITY(regionMats->getLocalNumEntries(), 105);

    // In the serial case we can just compare to the values in A
    entries_type refEntries                        = A->getLocalMatrixDevice().graph.entries;
    values_type refValues                          = A->getLocalMatrixDevice().values;
    typename entries_type::HostMirror refEntries_h = Kokkos::create_mirror_view(refEntries);
    Kokkos::deep_copy(refEntries_h, refEntries);
    typename values_type::HostMirror refValues_h = Kokkos::create_mirror_view(refValues);
    Kokkos::deep_copy(refValues_h, refValues);

    for (int idx = 0; idx < 105; ++idx) {
      TEST_EQUALITY(myEntries_h(idx), refEntries_h(idx));
      TEST_FLOATING_EQUALITY(TST::magnitude(myValues_h(idx)),
                             TST::magnitude(refValues_h(idx)),
                             100 * TMT::eps());
    }
  } else if (numRanks == 4) {
    // All ranks will have the same number of rows/cols/entries
    TEST_EQUALITY(regionMats->getGlobalNumRows(), 36);
    TEST_EQUALITY(regionMats->getGlobalNumCols(), 36);
    TEST_EQUALITY(regionMats->getLocalNumRows(), 9);
    TEST_EQUALITY(regionMats->getGlobalNumEntries(), 132);
    TEST_EQUALITY(regionMats->getLocalNumEntries(), 33);

    ArrayRCP<LO> refEntries;
    ArrayRCP<SC> refValues;
    refEntries.deepCopy(ArrayView<const LO>({0, 1, 3,
                                             0, 1, 2, 4,
                                             1, 2, 5,
                                             0, 3, 4, 6,
                                             1, 3, 4, 5, 7,
                                             2, 4, 5, 8,
                                             3, 6, 7,
                                             4, 6, 7, 8,
                                             5, 7, 8}));
    if (myRank == 0) {
      refValues.deepCopy(ArrayView<const SC>({4.0, -1.0, -1.0,
                                              -1.0, 4.0, -1.0, -1.0,
                                              -1.0, 2.0, -0.5,
                                              -1.0, 4.0, -1.0, -1.0,
                                              -1.0, -1.0, 4.0, -1.0, -1.0,
                                              -0.5, -1.0, 2.0, -0.5,
                                              -1.0, 2.0, -0.5,
                                              -1.0, -0.5, 2.0, -0.5,
                                              -0.5, -0.5, 1.0}));

    } else if (myRank == 1) {
      refValues.deepCopy(ArrayView<const SC>({2.0, -1.0, -0.5,
                                              -1.0, 4.0, -1.0, -1.0,
                                              -1.0, 4.0, -1.0,
                                              -0.5, 2.0, -1.0, -0.5,
                                              -1.0, -1.0, 4.0, -1.0, -1.0,
                                              -1.0, -1.0, 4.0, -1.0,
                                              -0.5, 1.0, -0.5,
                                              -1.0, -0.5, 2.0, -0.5,
                                              -1.0, -0.5, 2.0}));

    } else if (myRank == 2) {
      refValues.deepCopy(ArrayView<const SC>({2.0, -0.5, -1.0,
                                              -0.5, 2.0, -0.5, -1.0,
                                              -0.5, 1.0, -0.5,
                                              -1.0, 4.0, -1.0, -1.0,
                                              -1.0, -1.0, 4.0, -1.0, -1.0,
                                              -0.5, -1.0, 2.0, -0.5,
                                              -1.0, 4.0, -1.0,
                                              -1.0, -1.0, 4.0, -1.0,
                                              -0.5, -1.0, 2.0}));

    } else if (myRank == 3) {
      refValues.deepCopy(ArrayView<const SC>({1.0, -0.5, -0.5,
                                              -0.5, 2.0, -0.5, -1.0,
                                              -0.5, 2.0, -1.0,
                                              -0.5, 2.0, -1.0, -0.5,
                                              -1.0, -1.0, 4.0, -1.0, -1.0,
                                              -1.0, -1.0, 4.0, -1.0,
                                              -0.5, 2.0, -1.0,
                                              -1.0, -1.0, 4.0, -1.0,
                                              -1.0, -1.0, 4.0}));
    }

    // Loop over region matrix data and compare it to ref data
    for (int idx = 0; idx < 33; ++idx) {
      TEST_EQUALITY(myEntries_h(idx), refEntries[idx]);
      TEST_FLOATING_EQUALITY(TST::magnitude(myValues_h(idx)),
                             TST::magnitude(refValues[idx]),
                             100 * TMT::eps());
    }
  }

}  // CompositeToRegionMatrix

// Here the reverse operation: regionToComposite is tested.
// For that check the idea is to build a composite problem
// Then compute the equivalent region problem and finally
// use the region problem to go back to composite which should
// lead to a matrix identical to the original composite matrix
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RegionMatrix, RegionToCompositeMatrix, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;

  // Get MPI parameter
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  GO nx = 5, ny = 5, nz = 1;
  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP();
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D");
  Teuchos::ParameterList galeriList = galeriParameters.GetParameterList();
  std::string matrixType            = galeriParameters.GetMatrixType();

  // Build maps for the problem
  const LO numDofsPerNode = 1;
  RCP<Map> nodeMap        = Galeri::Xpetra::CreateMap<LO, GO, Node>(TestHelpers::Parameters::getLib(),
                                                             "Cartesian2D", comm, galeriList);
  RCP<Map> dofMap         = Xpetra::MapFactory<LO, GO, Node>::Build(nodeMap, numDofsPerNode);

  // Build the Xpetra problem
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(galeriParameters.GetMatrixType(), dofMap, galeriList);

  // Generate the operator
  RCP<Matrix> A = Pr->BuildMatrix();
  A->SetFixedBlockSize(numDofsPerNode);

  // Create the left and right handside vectors
  RCP<Vector> X = VectorFactory::Build(dofMap);
  RCP<Vector> B = VectorFactory::Build(dofMap);

  // Create auxiliary data for MG
  RCP<MultiVector> nullspace             = Pr->BuildNullspace();
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::
      CreateCartesianCoordinates<double, LO, GO, Map, RealValuedMultiVector>("2D", nodeMap, galeriList);

  // From the original composite matrix A, build the region equivalent: regionMats
  RCP<const Map> rowMap        = Teuchos::null;
  RCP<const Map> colMap        = Teuchos::null;
  RCP<const Map> revisedRowMap = Teuchos::null;
  RCP<const Map> revisedColMap = Teuchos::null;
  RCP<Import> rowImport        = Teuchos::null;
  RCP<Import> colImport        = Teuchos::null;
  RCP<Matrix> regionMats       = Teuchos::null;
  Teuchos::ArrayRCP<LO> regionMatVecLIDs;
  RCP<Import> regionInterfaceImporter;
  createRegionMatrix(galeriList, numDofsPerNode, nodeMap, dofMap, A,
                     rowMap, colMap, revisedRowMap, revisedColMap,
                     rowImport, colImport, regionMats,
                     regionMatVecLIDs, regionInterfaceImporter);

  // Finally do the revert operation: start with regionMats and bring it to composite format
  RCP<Matrix> compositeMatrix = MatrixFactory::Build(dofMap, 10);
  regionalToComposite(regionMats,
                      rowMap, colMap,
                      rowImport, Xpetra::INSERT,
                      compositeMatrix);

  // Now simply check that the original matrix A is identical to compositeMatrix.
  // This will ensure that createRegionMatrix() followed by regionalToComposite()
  // isequivalent to applying the identity operator.

  // Extract the local data from the original and final matrices to compare them
  using local_matrix_type = typename Matrix::local_matrix_type;
  using local_graph_type  = typename local_matrix_type::staticcrsgraph_type;
  using entries_type      = typename local_graph_type::entries_type;
  using values_type       = typename local_matrix_type::values_type;

  local_matrix_type orignalA = A->getLocalMatrixDevice();  // Local matrix
  entries_type refEntries    = orignalA.graph.entries;     // view of local column indices
  values_type refValues      = orignalA.values;            // view of local values

  typename entries_type::HostMirror refEntries_h = Kokkos::create_mirror_view(refEntries);
  Kokkos::deep_copy(refEntries_h, refEntries);
  typename values_type::HostMirror refValues_h = Kokkos::create_mirror_view(refValues);
  Kokkos::deep_copy(refValues_h, refValues);

  local_matrix_type compositeA  = compositeMatrix->getLocalMatrixDevice();  // Local matrix
  entries_type compositeEntries = compositeA.graph.entries;                 // view of local column indices
  values_type compositeValues   = compositeA.values;                        // view of local values

  typename entries_type::HostMirror compositeEntries_h = Kokkos::create_mirror_view(compositeEntries);
  Kokkos::deep_copy(compositeEntries_h, compositeEntries);
  typename values_type::HostMirror compositeValues_h = Kokkos::create_mirror_view(compositeValues);
  Kokkos::deep_copy(compositeValues_h, compositeValues);

  TEST_EQUALITY(compositeEntries_h.extent(0), refEntries_h.extent(0));
  TEST_EQUALITY(compositeValues_h.extent(0), refValues_h.extent(0));
  for (LO idx = 0; idx < compositeEntries_h.extent_int(0); ++idx) {
    TEST_EQUALITY(compositeEntries_h(idx), refEntries_h(idx));
    TEST_FLOATING_EQUALITY(TST::magnitude(compositeValues_h(idx)),
                           TST::magnitude(refValues_h(idx)),
                           100 * TMT::eps());
  }

}  // RegionToCompositeMatrix

// This test aims at checking that apply regionA to a region vector has the same effect as
// applying A to a composite vector. Of course the region vector needs to be brought back
// to composite formate before verifying the equivalence.
//
// Do this for 1 DOF per node
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RegionMatrix, FastMatVec, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;

  // Get MPI parameter
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  GO nx = 7, ny = 7, nz = 1;
  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP();
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D");
  Teuchos::ParameterList galeriList = galeriParameters.GetParameterList();
  std::string matrixType            = galeriParameters.GetMatrixType();

  // Build maps for the problem
  const LO numDofsPerNode = 1;
  RCP<Map> nodeMap        = Galeri::Xpetra::CreateMap<LO, GO, Node>(TestHelpers::Parameters::getLib(),
                                                             "Cartesian2D", comm, galeriList);
  RCP<Map> dofMap         = Xpetra::MapFactory<LO, GO, Node>::Build(nodeMap, numDofsPerNode);

  // Build the Xpetra problem
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(galeriParameters.GetMatrixType(), dofMap, galeriList);

  // Generate the operator
  RCP<Matrix> A = Pr->BuildMatrix();
  A->SetFixedBlockSize(numDofsPerNode);

  // Create auxiliary data for MG
  RCP<MultiVector> nullspace             = Pr->BuildNullspace();
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<double, LO, GO, Map, RealValuedMultiVector>("2D", nodeMap, galeriList);

  // create the region maps, importer and operator from composite counter parts
  RCP<const Map> rowMap        = Teuchos::null;
  RCP<const Map> colMap        = Teuchos::null;
  RCP<const Map> revisedRowMap = Teuchos::null;
  RCP<const Map> revisedColMap = Teuchos::null;
  RCP<Import> rowImport        = Teuchos::null;
  RCP<Import> colImport        = Teuchos::null;
  RCP<Matrix> regionMats       = Teuchos::null;
  Teuchos::ArrayRCP<LO> regionMatVecLIDs;
  RCP<Import> regionInterfaceImporter;
  createRegionMatrix(galeriList, numDofsPerNode, nodeMap, dofMap, A,
                     rowMap, colMap, revisedRowMap, revisedColMap,
                     rowImport, colImport, regionMats,
                     regionMatVecLIDs, regionInterfaceImporter);

  // Create initial vectors in composite format and apply composite A.
  // This will give a reference to compare with.
  RCP<Vector> X = VectorFactory::Build(dofMap);
  RCP<Vector> B = VectorFactory::Build(dofMap);

  // we set seed for reproducibility
  Utilities::SetRandomSeed(*comm);
  X->randomize();

  // Perform composite MatVec
  A->apply(*X, *B, Teuchos::NO_TRANS, TST::one(), TST::zero());

  // Create the region vectors and apply region A
  RCP<Vector> quasiRegX = Teuchos::null;
  RCP<Vector> quasiRegB = Teuchos::null;
  RCP<Vector> regX      = Teuchos::null;
  RCP<Vector> regB      = Teuchos::null;
  compositeToRegional(X, quasiRegX, regX, revisedRowMap, rowImport);
  regB = VectorFactory::Build(revisedRowMap, true);
  regionMats->apply(*regX, *regB, Teuchos::NO_TRANS, TST::one(), TST::zero());
  sumInterfaceValues(regB, revisedRowMap, rowImport);

  // Now create a refRegB vector using B as a starting point
  // and compare the result to regB
  RCP<Vector> refRegB = Teuchos::null;
  compositeToRegional(B, quasiRegB, refRegB, revisedRowMap, rowImport);

  // Extract the data from B and compB to compare it
  ArrayRCP<const SC> dataRegB    = regB->getData(0);
  ArrayRCP<const SC> dataRefRegB = refRegB->getData(0);
  for (size_t idx = 0; idx < refRegB->getLocalLength(); ++idx) {
    TEST_FLOATING_EQUALITY(TST::magnitude(dataRegB[idx]),
                           TST::magnitude(dataRefRegB[idx]),
                           100 * TMT::eps());
  }

  // Finally we perform the "fastMatVec" that does not require
  // to transform data from region to composite and back
  // it should perform faster and allow for easy customization
  // of the local MatVec
  RCP<const Map> regionMap = revisedRowMap;

  RCP<Vector> regC = Teuchos::null;
  regC             = VectorFactory::Build(revisedRowMap, true);
  regionMats->apply(*regX, *regC, Teuchos::NO_TRANS, TST::one(), TST::zero(), true, regionInterfaceImporter, regionMatVecLIDs);

  ArrayRCP<const SC> dataRegC = regC->getData(0);
  for (size_t idx = 0; idx < refRegB->getLocalLength(); ++idx) {
    TEST_FLOATING_EQUALITY(TST::magnitude(dataRegC[idx]),
                           TST::magnitude(dataRefRegB[idx]),
                           100 * TMT::eps());
  }

}  // FastMatVec

// This test aims at checking that apply regionA to a region vector has the same effect as
// applying A to a composite vector. Of course the region vector needs to be brought back
// to composite formate before verifying the equivalence.
//
// Do this for 1 DOF per node
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RegionMatrix, FastMatVec3D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;

  // Get MPI parameter
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  GO nx = 5, ny = 5, nz = 3;
  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP();
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace3D");
  Teuchos::ParameterList galeriList = galeriParameters.GetParameterList();
  std::string matrixType            = galeriParameters.GetMatrixType();
  const LO numDofsPerNode           = 1;

  // create the region maps, importer and operator from composite counter parts
  RCP<Matrix> regionMats = Teuchos::null;
  RCP<Matrix> A;
  RCP<const Map> revisedRowMap = Teuchos::null;
  RCP<Import> rowImport        = Teuchos::null;
  Teuchos::ArrayRCP<LocalOrdinal> regionMatVecLIDs;
  RCP<Import> regionInterfaceImporter;

  createProblem(numDofsPerNode, galeriParameters, comm,
                A, regionMats, revisedRowMap, rowImport,
                regionMatVecLIDs, regionInterfaceImporter);

  // Create initial vectors in composite format and apply composite A.
  // This will give a reference to compare with.
  RCP<Vector> X = VectorFactory::Build(A->getRowMap());
  RCP<Vector> B = VectorFactory::Build(A->getRowMap());

  // we set seed for reproducibility
  Utilities::SetRandomSeed(*comm);
  X->randomize();

  // Perform composite MatVec
  A->apply(*X, *B, Teuchos::NO_TRANS, TST::one(), TST::zero());

  // Create the region vectors and apply region A
  RCP<Vector> quasiRegX = Teuchos::null;
  RCP<Vector> quasiRegB = Teuchos::null;
  RCP<Vector> regX      = Teuchos::null;
  RCP<Vector> regB      = Teuchos::null;
  compositeToRegional(X, quasiRegX, regX,
                      revisedRowMap, rowImport);
  regB = VectorFactory::Build(revisedRowMap, true);
  regionMats->apply(*regX, *regB, Teuchos::NO_TRANS, TST::one(), TST::zero());
  sumInterfaceValues(regB, revisedRowMap, rowImport);

  // Now create a refRegB vector using B as a starting point
  // and compare the result to regB
  RCP<Vector> refRegB = Teuchos::null;
  compositeToRegional(B, quasiRegB, refRegB, revisedRowMap, rowImport);

  // Extract the data from B and compB to compare it
  ArrayRCP<const SC> dataRegB    = regB->getData(0);
  ArrayRCP<const SC> dataRefRegB = refRegB->getData(0);
  for (size_t idx = 0; idx < refRegB->getLocalLength(); ++idx) {
    TEST_FLOATING_EQUALITY(TST::magnitude(dataRegB[idx]),
                           TST::magnitude(dataRefRegB[idx]),
                           100 * TMT::eps());
  }

  // Finally we perform the "fastMatVec" that does not require
  // to transform data from region to composite and back
  // it should perform faster and allow for easy customization
  // of the local MatVec
  RCP<const Map> regionMap = revisedRowMap;

  RCP<Vector> regC = Teuchos::null;
  regC             = VectorFactory::Build(revisedRowMap, true);
  regionMats->apply(*regX, *regC, Teuchos::NO_TRANS, TST::one(), TST::zero(), true, regionInterfaceImporter, regionMatVecLIDs);

  ArrayRCP<const SC> dataRegC = regC->getData(0);
  for (size_t idx = 0; idx < refRegB->getLocalLength(); ++idx) {
    TEST_FLOATING_EQUALITY(TST::magnitude(dataRegC[idx]),
                           TST::magnitude(dataRefRegB[idx]),
                           100 * TMT::eps());
  }

}  // FastMatVec3D

// This test aims at checking that apply regionA to a region vector has the same effect as
// applying A to a composite vector. Of course the region vector needs to be brought back
// to composite formate before verifying the equivalence.
//
// Do this for 2 DOFs per node (two-dimensional elasticity)
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RegionMatrix, FastMatVec2D_Elasticity, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;

  // Get MPI parameter
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  const LO numDofsPerNode = 2;
  GO nx = 7, ny = 7, nz = 1;
  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP();
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Elasticity2D");
  RCP<Matrix> regionMats = Teuchos::null;
  RCP<Matrix> A;
  RCP<const Map> revisedRowMap = Teuchos::null;
  RCP<Import> rowImport        = Teuchos::null;
  Teuchos::ArrayRCP<LocalOrdinal> regionMatVecLIDs;
  RCP<Import> regionInterfaceImporter;

  createProblem(numDofsPerNode, galeriParameters, comm,
                A, regionMats, revisedRowMap, rowImport,
                regionMatVecLIDs, regionInterfaceImporter);

  // Create initial vectors in composite format and apply composite A.
  // This will give a reference to compare with.
  RCP<Vector> X = VectorFactory::Build(A->getRowMap());
  RCP<Vector> B = VectorFactory::Build(A->getRowMap());

  // we set seed for reproducibility
  Utilities::SetRandomSeed(*comm);
  X->randomize();

  // Perform composite MatVec
  A->apply(*X, *B, Teuchos::NO_TRANS, TST::one(), TST::zero());

  // Create the region vectors and apply region A
  RCP<Vector> quasiRegX = Teuchos::null;
  RCP<Vector> quasiRegB = Teuchos::null;
  RCP<Vector> regX      = Teuchos::null;
  RCP<Vector> regB      = Teuchos::null;
  compositeToRegional(X, quasiRegX, regX, revisedRowMap, rowImport);
  regB = VectorFactory::Build(revisedRowMap, true);
  regionMats->apply(*regX, *regB, Teuchos::NO_TRANS, TST::one(), TST::zero());
  sumInterfaceValues(regB, revisedRowMap, rowImport);

  // Now create a refRegB vector using B as a starting point
  // and compare the result to regB
  RCP<Vector> refRegB = Teuchos::null;
  compositeToRegional(B, quasiRegB, refRegB, revisedRowMap, rowImport);

  // Extract the data from B and compB to compare it
  ArrayRCP<const SC> dataRegB    = regB->getData(0);
  ArrayRCP<const SC> dataRefRegB = refRegB->getData(0);
  for (size_t idx = 0; idx < refRegB->getLocalLength(); ++idx) {
    TEST_FLOATING_EQUALITY(TST::magnitude(dataRegB[idx]),
                           TST::magnitude(dataRefRegB[idx]),
                           100 * TMT::eps());
  }

  // Finally we perform the "fastMatVec" that does not require
  // to transform data from region to composite and back
  // it should perform faster and allow for easy customization
  // of the local MatVec
  RCP<const Map> regionMap = revisedRowMap;

  RCP<Vector> regC = Teuchos::null;
  regC             = VectorFactory::Build(revisedRowMap, true);
  regionMats->apply(*regX, *regC, Teuchos::NO_TRANS, TST::one(), TST::zero(), true, regionInterfaceImporter, regionMatVecLIDs);

  ArrayRCP<const SC> dataRegC = regC->getData(0);
  for (size_t idx = 0; idx < refRegB->getLocalLength(); ++idx) {
    TEST_FLOATING_EQUALITY(TST::magnitude(dataRegC[idx]),
                           TST::magnitude(dataRefRegB[idx]),
                           100 * TMT::eps());
  }

}  // FastMatVec2D_Elasticity

// This test aims at checking that apply regionA to a region vector has the same effect as
// applying A to a composite vector. Of course the region vector needs to be brought back
// to composite formate before verifying the equivalence.
//
// Do this for 3 DOFs per node (three-dimensional elasticity)
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RegionMatrix, FastMatVec3D_Elasticity, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;

  // Get MPI parameter
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  const LO numDofsPerNode = 3;
  GO nx = 5, ny = 5, nz = 3;
  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP();
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Elasticity3D");
  RCP<Matrix> regionMats       = Teuchos::null;
  RCP<Matrix> A                = Teuchos::null;
  RCP<const Map> revisedRowMap = Teuchos::null;
  RCP<Import> rowImport        = Teuchos::null;
  Teuchos::ArrayRCP<LocalOrdinal> regionMatVecLIDs;
  RCP<Import> regionInterfaceImporter;

  createProblem(numDofsPerNode, galeriParameters, comm,
                A, regionMats, revisedRowMap, rowImport,
                regionMatVecLIDs, regionInterfaceImporter);

  // Create initial vectors in composite format and apply composite A.
  // This will give a reference to compare with.
  RCP<Vector> X = VectorFactory::Build(A->getRowMap());
  RCP<Vector> B = VectorFactory::Build(A->getRowMap());

  // we set seed for reproducibility
  Utilities::SetRandomSeed(*comm);
  X->randomize();

  // Perform composite MatVec
  A->apply(*X, *B, Teuchos::NO_TRANS, TST::one(), TST::zero());

  // Create the region vectors and apply region A
  RCP<Vector> quasiRegX = Teuchos::null;
  RCP<Vector> quasiRegB = Teuchos::null;
  RCP<Vector> regX      = Teuchos::null;
  RCP<Vector> regB      = Teuchos::null;
  compositeToRegional(X, quasiRegX, regX, revisedRowMap, rowImport);
  regB = VectorFactory::Build(revisedRowMap, true);
  regionMats->apply(*regX, *regB, Teuchos::NO_TRANS, TST::one(), TST::zero());
  sumInterfaceValues(regB, revisedRowMap, rowImport);

  // Now create a refRegB vector using B as a starting point
  // and compare the result to regB
  RCP<Vector> refRegB = Teuchos::null;
  compositeToRegional(B, quasiRegB, refRegB, revisedRowMap, rowImport);

  // Extract the data from B and compB to compare it
  ArrayRCP<const SC> dataRegB    = regB->getData(0);
  ArrayRCP<const SC> dataRefRegB = refRegB->getData(0);
  for (size_t idx = 0; idx < refRegB->getLocalLength(); ++idx) {
    TEST_FLOATING_EQUALITY(TST::magnitude(dataRegB[idx]),
                           TST::magnitude(dataRefRegB[idx]),
                           300 * TMT::eps());
  }

  // Finally we perform the "fastMatVec" that does not require
  // to transform data from region to composite and back
  // it should perform faster and allow for easy customization
  // of the local MatVec
  RCP<const Map> regionMap = revisedRowMap;

  RCP<Vector> regC = Teuchos::null;
  regC             = VectorFactory::Build(revisedRowMap, true);
  regionMats->apply(*regX, *regC, Teuchos::NO_TRANS, TST::one(), TST::zero(), true, regionInterfaceImporter, regionMatVecLIDs);

  ArrayRCP<const SC> dataRegC = regC->getData(0);
  for (size_t idx = 0; idx < refRegB->getLocalLength(); ++idx) {
    TEST_FLOATING_EQUALITY(TST::magnitude(dataRegC[idx]),
                           TST::magnitude(dataRefRegB[idx]),
                           300 * TMT::eps());
  }

}  // FastMatVec3D_Elasticity

// Here a Laplace 2D problem is tested for all the above checks mentioned:
//   1) the region operator is compared against know values
//   2) the action of the region MatVec is compared with the composite MatVec
//   3) compute the composite operator from the region operator leads to the original matrix
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RegionMatrix, Laplace2D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;

  // Get MPI parameter
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  GO nx = 6, ny = 5, nz = 1;
  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP();
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D");
  Teuchos::ParameterList galeriList = galeriParameters.GetParameterList();
  std::string matrixType            = galeriParameters.GetMatrixType();

  // Build maps for the problem
  const LO numDofsPerNode = 1;
  RCP<Map> nodeMap        = Galeri::Xpetra::CreateMap<LO, GO, Node>(TestHelpers::Parameters::getLib(),
                                                             "Cartesian2D", comm, galeriList);
  RCP<Map> dofMap         = Xpetra::MapFactory<LO, GO, Node>::Build(nodeMap, numDofsPerNode);

  // Build the Xpetra problem
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(galeriParameters.GetMatrixType(), dofMap, galeriList);

  // Generate the operator
  RCP<Matrix> A = Pr->BuildMatrix();
  A->SetFixedBlockSize(numDofsPerNode);

  // Create auxiliary data for MG
  RCP<MultiVector> nullspace             = Pr->BuildNullspace();
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<double, LO, GO, Map, RealValuedMultiVector>("2D", nodeMap, galeriList);

  RCP<const Map> rowMap        = Teuchos::null;
  RCP<const Map> colMap        = Teuchos::null;
  RCP<const Map> revisedRowMap = Teuchos::null;
  RCP<const Map> revisedColMap = Teuchos::null;
  RCP<Import> rowImport        = Teuchos::null;
  RCP<Import> colImport        = Teuchos::null;
  RCP<Matrix> regionMats       = Teuchos::null;
  Teuchos::ArrayRCP<LO> regionMatVecLIDs;
  RCP<Import> regionInterfaceImporter;
  createRegionMatrix(galeriList, numDofsPerNode, nodeMap, dofMap, A,
                     rowMap, colMap, revisedRowMap, revisedColMap,
                     rowImport, colImport, regionMats,
                     regionMatVecLIDs, regionInterfaceImporter);

  test_matrix(A, regionMats,
              rowMap, colMap, revisedRowMap, rowImport,
              out, success);

  // Extract the local data from the region matrix
  using local_matrix_type = typename Matrix::local_matrix_type;
  using local_graph_type  = typename local_matrix_type::staticcrsgraph_type;
  using entries_type      = typename local_graph_type::entries_type;
  using values_type       = typename local_matrix_type::values_type;

  local_matrix_type myLocalA = regionMats->getLocalMatrixDevice();  // Local matrix
  entries_type myEntries     = myLocalA.graph.entries;              // view of local column indices
  values_type myValues       = myLocalA.values;                     // view of local values

  typename entries_type::HostMirror myEntries_h = Kokkos::create_mirror_view(myEntries);
  Kokkos::deep_copy(myEntries_h, myEntries);
  typename values_type::HostMirror myValues_h = Kokkos::create_mirror_view(myValues);
  Kokkos::deep_copy(myValues_h, myValues);

  const int numRanks = comm->getSize();
  const int myRank   = comm->getRank();
  if (numRanks == 1) {
    TEST_EQUALITY(regionMats->getGlobalNumRows(), 30);
    TEST_EQUALITY(regionMats->getGlobalNumCols(), 30);
    TEST_EQUALITY(regionMats->getLocalNumRows(), 30);
    TEST_EQUALITY(regionMats->getGlobalNumEntries(), 128);
    TEST_EQUALITY(regionMats->getLocalNumEntries(), 128);

    // In the serial case we can just compare to the values in A
    entries_type refEntries                        = A->getLocalMatrixDevice().graph.entries;
    values_type refValues                          = A->getLocalMatrixDevice().values;
    typename entries_type::HostMirror refEntries_h = Kokkos::create_mirror_view(refEntries);
    Kokkos::deep_copy(refEntries_h, refEntries);
    typename values_type::HostMirror refValues_h = Kokkos::create_mirror_view(refValues);
    Kokkos::deep_copy(refValues_h, refValues);

    for (int idx = 0; idx < 128; ++idx) {
      TEST_EQUALITY(myEntries_h(idx), refEntries_h(idx));
      TEST_FLOATING_EQUALITY(TST::magnitude(myValues_h(idx)),
                             TST::magnitude(refValues_h(idx)),
                             100 * TMT::eps());
    }
  } else if (numRanks == 4) {
    // All ranks will have the same number of rows/cols/entries
    TEST_EQUALITY(regionMats->getGlobalNumRows(), 42);
    TEST_EQUALITY(regionMats->getGlobalNumCols(), 42);
    TEST_EQUALITY(regionMats->getGlobalNumEntries(), 158);

    ArrayRCP<SC> refValues;
    if (myRank == 0) {
      TEST_EQUALITY(regionMats->getLocalNumRows(), 9);
      TEST_EQUALITY(regionMats->getLocalNumEntries(), 33);
      refValues.deepCopy(ArrayView<const SC>({4.0, -1.0, -1.0,
                                              -1.0, 4.0, -1.0, -1.0,
                                              -1.0, 2.0, -0.5,
                                              -1.0, 4.0, -1.0, -1.0,
                                              -1.0, -1.0, 4.0, -1.0, -1.0,
                                              -0.5, -1.0, 2.0, -0.5,
                                              -1.0, 2.0, -0.5,
                                              -1.0, -0.5, 2.0, -0.5,
                                              -0.5, -0.5, 1.0}));

    } else if (myRank == 1) {
      TEST_EQUALITY(regionMats->getLocalNumRows(), 12);
      TEST_EQUALITY(regionMats->getLocalNumEntries(), 46);
      refValues.deepCopy(ArrayView<const SC>({2.0, -1.0, -0.5,
                                              -1.0, 4.0, -1.0, -1.0,
                                              -1.0, 4.0, -1.0, -1.0,
                                              -1.0, 4.0, -1.0,
                                              -0.5, 2.0, -1.0, -0.5,
                                              -1.0, -1.0, 4.0, -1.0, -1.0,
                                              -1.0, -1.0, 4.0, -1.0, -1.0,
                                              -1.0, -1.0, 4.0, -1.0,
                                              -0.5, 1.0, -0.5,
                                              -1.0, -0.5, 2.0, -0.5,
                                              -1.0, -0.5, 2.0, -0.5,
                                              -1.0, -0.5, 2.0}));

    } else if (myRank == 2) {
      TEST_EQUALITY(regionMats->getLocalNumRows(), 9);
      TEST_EQUALITY(regionMats->getLocalNumEntries(), 33);
      refValues.deepCopy(ArrayView<const SC>({2.0, -0.5, -1.0,
                                              -0.5, 2.0, -0.5, -1.0,
                                              -0.5, 1.0, -0.5,
                                              -1.0, 4.0, -1.0, -1.0,
                                              -1.0, -1.0, 4.0, -1.0, -1.0,
                                              -0.5, -1.0, 2.0, -0.5,
                                              -1.0, 4.0, -1.0,
                                              -1.0, -1.0, 4.0, -1.0,
                                              -0.5, -1.0, 2.0}));

    } else if (myRank == 3) {
      TEST_EQUALITY(regionMats->getLocalNumRows(), 12);
      TEST_EQUALITY(regionMats->getLocalNumEntries(), 46);
      refValues.deepCopy(ArrayView<const SC>({1.0, -0.5, -0.5,
                                              -0.5, 2.0, -0.5, -1.0,
                                              -0.5, 2.0, -0.5, -1.0,
                                              -0.5, 2.0, -1.0,
                                              -0.5, 2.0, -1.0, -0.5,
                                              -1.0, -1.0, 4.0, -1.0, -1.0,
                                              -1.0, -1.0, 4.0, -1.0, -1.0,
                                              -1.0, -1.0, 4.0, -1.0,
                                              -0.5, 2.0, -1.0,
                                              -1.0, -1.0, 4.0, -1.0,
                                              -1.0, -1.0, 4.0, -1.0,
                                              -1.0, -1.0, 4.0}));
    }

    // Loop over region matrix data and compare it to ref data
    for (int idx = 0; idx < 33; ++idx) {
      TEST_FLOATING_EQUALITY(TST::magnitude(myValues_h(idx)),
                             TST::magnitude(refValues[idx]),
                             100 * TMT::eps());
    }
  }
}  // Laplace2D

// Here a Laplace 3D problem is tested for all the above checks mentioned:
//   1) the region operator is compared against know values
//   2) the action of the region MatVec is compared with the composite MatVec
//   3) compute the composite operator from the region operator leads to the original matrix
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RegionMatrix, Laplace3D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;

  // Get MPI parameter
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  GO nx = 6, ny = 5, nz = 4;
  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP();
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace3D");
  Teuchos::ParameterList galeriList = galeriParameters.GetParameterList();
  std::string matrixType            = galeriParameters.GetMatrixType();

  // Build maps for the problem
  const LO numDofsPerNode = 1;
  RCP<Map> nodeMap        = Galeri::Xpetra::CreateMap<LO, GO, Node>(TestHelpers::Parameters::getLib(),
                                                             "Cartesian3D", comm, galeriList);
  RCP<Map> dofMap         = Xpetra::MapFactory<LO, GO, Node>::Build(nodeMap, numDofsPerNode);

  // Build the Xpetra problem
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(galeriParameters.GetMatrixType(), dofMap, galeriList);

  // Generate the operator
  RCP<Matrix> A = Pr->BuildMatrix();
  A->SetFixedBlockSize(numDofsPerNode);

  // Create auxiliary data for MG
  RCP<MultiVector> nullspace             = Pr->BuildNullspace();
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<double, LO, GO, Map, RealValuedMultiVector>("3D", nodeMap, galeriList);

  RCP<const Map> rowMap        = Teuchos::null;
  RCP<const Map> colMap        = Teuchos::null;
  RCP<const Map> revisedRowMap = Teuchos::null;
  RCP<const Map> revisedColMap = Teuchos::null;
  RCP<Import> rowImport        = Teuchos::null;
  RCP<Import> colImport        = Teuchos::null;
  RCP<Matrix> regionMats       = Teuchos::null;
  Teuchos::ArrayRCP<LO> regionMatVecLIDs;
  RCP<Import> regionInterfaceImporter;
  createRegionMatrix(galeriList, numDofsPerNode, nodeMap, dofMap, A,
                     rowMap, colMap, revisedRowMap, revisedColMap,
                     rowImport, colImport, regionMats,
                     regionMatVecLIDs, regionInterfaceImporter);

  test_matrix(A, regionMats,
              rowMap, colMap, revisedRowMap, rowImport,
              out, success);

  // Extract the local data from the region matrix
  using local_matrix_type = typename Matrix::local_matrix_type;
  using local_graph_type  = typename local_matrix_type::staticcrsgraph_type;
  using entries_type      = typename local_graph_type::entries_type;
  using values_type       = typename local_matrix_type::values_type;

  local_matrix_type myLocalA = regionMats->getLocalMatrixDevice();  // Local matrix
  entries_type myEntries     = myLocalA.graph.entries;              // view of local column indices
  values_type myValues       = myLocalA.values;                     // view of local values

  typename entries_type::HostMirror myEntries_h = Kokkos::create_mirror_view(myEntries);
  Kokkos::deep_copy(myEntries_h, myEntries);
  typename values_type::HostMirror myValues_h = Kokkos::create_mirror_view(myValues);
  Kokkos::deep_copy(myValues_h, myValues);

  const int numRanks = comm->getSize();
  const int myRank   = comm->getRank();
  if ((numRanks == 1) && (myRank == 0)) {
    TEST_EQUALITY(regionMats->getGlobalNumRows(), 120);
    TEST_EQUALITY(regionMats->getGlobalNumCols(), 120);
    TEST_EQUALITY(regionMats->getLocalNumRows(), 120);
    TEST_EQUALITY(regionMats->getGlobalNumEntries(), 692);
    TEST_EQUALITY(regionMats->getLocalNumEntries(), 692);

    // In the serial case we can just compare to the values in A
    entries_type refEntries                        = A->getLocalMatrixDevice().graph.entries;
    values_type refValues                          = A->getLocalMatrixDevice().values;
    typename entries_type::HostMirror refEntries_h = Kokkos::create_mirror_view(refEntries);
    Kokkos::deep_copy(refEntries_h, refEntries);
    typename values_type::HostMirror refValues_h = Kokkos::create_mirror_view(refValues);
    Kokkos::deep_copy(refValues_h, refValues);

    for (int idx = 0; idx < refEntries_h.extent_int(0); ++idx) {
      TEST_EQUALITY(myEntries_h(idx), refEntries_h(idx));
      TEST_FLOATING_EQUALITY(TST::magnitude(myValues_h(idx)),
                             TST::magnitude(refValues_h(idx)),
                             100 * TMT::eps());
    }
  } else if (numRanks == 4) {
    // All ranks will have the same number of rows/cols/entries
    TEST_EQUALITY(regionMats->getGlobalNumRows(), 168);
    TEST_EQUALITY(regionMats->getGlobalNumCols(), 168);
    TEST_EQUALITY(regionMats->getGlobalNumEntries(), 884);

    ArrayRCP<SC> refValues;
    if (myRank == 0) {
      TEST_EQUALITY(regionMats->getLocalNumRows(), 36);
      TEST_EQUALITY(regionMats->getLocalNumEntries(), 186);
      refValues.deepCopy(ArrayView<const SC>({6.0, -1.0, -1.0, -1.0,
                                              -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, 3.0, -0.5, -0.5,
                                              -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -0.5, -1.0, 3.0, -0.5, -0.5,
                                              -1.0, 3.0, -0.5, -0.5,
                                              -1.0, -0.5, 3.0, -0.5, -0.5,
                                              -0.5, -0.5, 1.5, -0.25,
                                              -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -0.5, -1.0, 3.0, -0.5, -0.5,
                                              -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -0.5, -0.5, -1.0, 3.0, -0.5, -0.5,
                                              -0.5, -1.0, 3.0, -0.5, -0.5,
                                              -0.5, -1.0, -0.5, 3.0, -0.5, -0.5,
                                              -0.25, -0.5, -0.5, 1.5, -0.25,
                                              -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -0.5, -1.0, 3.0, -0.5, -0.5,
                                              -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -0.5, -0.5, -1.0, 3.0, -0.5, -0.5,
                                              -0.5, -1.0, 3.0, -0.5, -0.5,
                                              -0.5, -1.0, -0.5, 3.0, -0.5, -0.5,
                                              -0.25, -0.5, -0.5, 1.5, -0.25,
                                              -1.0, 6.0, -1.0, -1.0,
                                              -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -0.5, -1.0, 3.0, -0.5,
                                              -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -0.5, -0.5, -1.0, 3.0, -0.5,
                                              -0.5, -1.0, 3.0, -0.5,
                                              -0.5, -1.0, -0.5, 3.0, -0.5,
                                              -0.25, -0.5, -0.5, 1.5}));

      // Loop over region matrix data and compare it to ref data
      for (int idx = 0; idx < 186; ++idx) {
        TEST_FLOATING_EQUALITY(TST::magnitude(myValues_h(idx)),
                               TST::magnitude(refValues[idx]),
                               100 * TMT::eps());
      }

    } else if (myRank == 1) {
      TEST_EQUALITY(regionMats->getLocalNumRows(), 48);
      TEST_EQUALITY(regionMats->getLocalNumEntries(), 256);
      refValues.deepCopy(ArrayView<const SC>({3.0, -1.0, -0.5, -0.5,
                                              -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, 6.0, -1.0, -1.0,
                                              -0.5, 3.0, -1.0, -0.5, -0.5,
                                              -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -0.5, 1.5, -0.5, -0.25,
                                              -1.0, -0.5, 3.0, -0.5, -0.5,
                                              -1.0, -0.5, 3.0, -0.5, -0.5,
                                              -1.0, -0.5, 3.0, -0.5,
                                              -0.5, 3.0, -1.0, -0.5, -0.5,
                                              -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -0.5, -0.5, 3.0, -1.0, -0.5, -0.5,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -0.25, -0.5, 1.5, -0.5, -0.25,
                                              -0.5, -1.0, -0.5, 3.0, -0.5, -0.5,
                                              -0.5, -1.0, -0.5, 3.0, -0.5, -0.5,
                                              -0.5, -1.0, -0.5, 3.0, -0.5,
                                              -0.5, 3.0, -1.0, -0.5, -0.5,
                                              -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -0.5, -0.5, 3.0, -1.0, -0.5, -0.5,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -0.25, -0.5, 1.5, -0.5, -0.25,
                                              -0.5, -1.0, -0.5, 3.0, -0.5, -0.5,
                                              -0.5, -1.0, -0.5, 3.0, -0.5, -0.5,
                                              -0.5, -1.0, -0.5, 3.0, -0.5,
                                              -0.5, 3.0, -1.0, -0.5,
                                              -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -1.0, -1.0, 6.0, -1.0,
                                              -0.5, -0.5, 3.0, -1.0, -0.5,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0,
                                              -0.25, -0.5, 1.5, -0.5,
                                              -0.5, -1.0, -0.5, 3.0, -0.5,
                                              -0.5, -1.0, -0.5, 3.0, -0.5,
                                              -0.5, -1.0, -0.5, 3.0}));

      // Loop over region matrix data and compare it to ref data
      for (int idx = 0; idx < 256; ++idx) {
        TEST_FLOATING_EQUALITY(TST::magnitude(myValues_h(idx)),
                               TST::magnitude(refValues[idx]),
                               100 * TMT::eps());
      }

    } else if (myRank == 2) {
      TEST_EQUALITY(regionMats->getLocalNumRows(), 36);
      TEST_EQUALITY(regionMats->getLocalNumEntries(), 186);
      refValues.deepCopy(ArrayView<const SC>({3.0, -0.5, -1.0, -0.5,
                                              -0.5, 3.0, -0.5, -1.0, -0.5,
                                              -0.5, 1.5, -0.5, -0.25,
                                              -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -0.5, -1.0, 3.0, -0.5, -0.5,
                                              -1.0, 6.0, -1.0, -1.0,
                                              -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -0.5, -1.0, 3.0, -0.5,
                                              -0.5, 3.0, -0.5, -1.0, -0.5,
                                              -0.5, -0.5, 3.0, -0.5, -1.0, -0.5,
                                              -0.25, -0.5, 1.5, -0.5, -0.25,
                                              -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -0.5, -0.5, -1.0, 3.0, -0.5, -0.5,
                                              -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -0.5, -0.5, -1.0, 3.0, -0.5,
                                              -0.5, 3.0, -0.5, -1.0, -0.5,
                                              -0.5, -0.5, 3.0, -0.5, -1.0, -0.5,
                                              -0.25, -0.5, 1.5, -0.5, -0.25,
                                              -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -0.5, -0.5, -1.0, 3.0, -0.5, -0.5,
                                              -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -0.5, -0.5, -1.0, 3.0, -0.5,
                                              -0.5, 3.0, -0.5, -1.0,
                                              -0.5, -0.5, 3.0, -0.5, -1.0,
                                              -0.25, -0.5, 1.5, -0.5,
                                              -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -0.5, -0.5, -1.0, 3.0, -0.5,
                                              -1.0, -1.0, 6.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0,
                                              -0.5, -0.5, -1.0, 3.0}));

      // Loop over region matrix data and compare it to ref data
      for (int idx = 0; idx < 186; ++idx) {
        TEST_FLOATING_EQUALITY(TST::magnitude(myValues_h(idx)),
                               TST::magnitude(refValues[idx]),
                               100 * TMT::eps());
      }

    } else if (myRank == 3) {
      TEST_EQUALITY(regionMats->getLocalNumRows(), 48);
      TEST_EQUALITY(regionMats->getLocalNumEntries(), 256);
      refValues.deepCopy(ArrayView<const SC>({1.5, -0.5, -0.5, -0.25,
                                              -0.5, 3.0, -0.5, -1.0, -0.5,
                                              -0.5, 3.0, -0.5, -1.0, -0.5,
                                              -0.5, 3.0, -1.0, -0.5,
                                              -0.5, 3.0, -1.0, -0.5, -0.5,
                                              -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -0.5, 3.0, -1.0, -0.5,
                                              -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -1.0, -1.0, 6.0, -1.0,
                                              -0.25, 1.5, -0.5, -0.5, -0.25,
                                              -0.5, -0.5, 3.0, -0.5, -1.0, -0.5,
                                              -0.5, -0.5, 3.0, -0.5, -1.0, -0.5,
                                              -0.5, -0.5, 3.0, -1.0, -0.5,
                                              -0.5, -0.5, 3.0, -1.0, -0.5, -0.5,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -0.5, -0.5, 3.0, -1.0, -0.5,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0,
                                              -0.25, 1.5, -0.5, -0.5, -0.25,
                                              -0.5, -0.5, 3.0, -0.5, -1.0, -0.5,
                                              -0.5, -0.5, 3.0, -0.5, -1.0, -0.5,
                                              -0.5, -0.5, 3.0, -1.0, -0.5,
                                              -0.5, -0.5, 3.0, -1.0, -0.5, -0.5,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -0.5, -0.5, 3.0, -1.0, -0.5,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0,
                                              -0.25, 1.5, -0.5, -0.5,
                                              -0.5, -0.5, 3.0, -0.5, -1.0,
                                              -0.5, -0.5, 3.0, -0.5, -1.0,
                                              -0.5, -0.5, 3.0, -1.0,
                                              -0.5, -0.5, 3.0, -1.0, -0.5,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0,
                                              -0.5, -0.5, 3.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0, -1.0,
                                              -1.0, -1.0, -1.0, 6.0}));

      // Loop over region matrix data and compare it to ref data
      for (int idx = 0; idx < 256; ++idx) {
        TEST_FLOATING_EQUALITY(TST::magnitude(myValues_h(idx)),
                               TST::magnitude(refValues[idx]),
                               100 * TMT::eps());
      }
    }
  }
}  // Laplace3D

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node)                                                       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RegionMatrix, CompositeToRegionMatrix, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RegionMatrix, RegionToCompositeMatrix, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RegionMatrix, FastMatVec, Scalar, LO, GO, Node)              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RegionMatrix, FastMatVec3D, Scalar, LO, GO, Node)            \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RegionMatrix, FastMatVec2D_Elasticity, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RegionMatrix, FastMatVec3D_Elasticity, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RegionMatrix, Laplace2D, Scalar, LO, GO, Node)               \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RegionMatrix, Laplace3D, Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
