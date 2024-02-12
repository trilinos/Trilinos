// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Galeri_XpetraParameters.hpp>

#include "MueLu_RegionRFactory.hpp"
#include "MueLu_RegionRFactory_kokkos.hpp"

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
                        Teuchos::Array<GlobalOrdinal>& quasiRegionCoordGIDs,
                        Teuchos::RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> >& regionInterfaceImporter,
                        Teuchos::Array<LocalOrdinal>& rNodesPerDim,
                        LocalOrdinal& numLocalRegionNodes) {
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
  int maxRegPerGID  = 0;
  int numInterfaces = 0;
  Array<GO> sendGIDs;
  Array<int> sendPIDs;
  Array<LO> compositeToRegionLIDs(nodeMap->getLocalNumElements() * numDofsPerNode);
  Array<GO> quasiRegionGIDs;
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
                   RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& regionNullspace,
                   RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node> >& regionCoordinates,
                   RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >& revisedRowMap,
                   RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> >& rowImport,
                   Teuchos::ArrayRCP<LocalOrdinal>& regionMatVecLIDs,
                   RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> >& regionInterfaceImporter,
                   Teuchos::Array<LocalOrdinal>& rNodesPerDim) {
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
  RCP<MultiVector> nullspace = Pr->BuildNullspace();
  RCP<RealValuedMultiVector> coordinates =
      Galeri::Xpetra::Utils::CreateCartesianCoordinates<double, LO, GO, Map, RealValuedMultiVector>(coordinatesType, nodeMap, galeriList);

  // create the region maps, importer and operator from composite counter parts
  RCP<const Map> rowMap        = Teuchos::null;
  RCP<const Map> colMap        = Teuchos::null;
  RCP<const Map> revisedColMap = Teuchos::null;
  RCP<Import> colImport        = Teuchos::null;
  Array<GO> quasiRegionCoordGIDs;
  LO numLocalRegionNodes = 0;
  createRegionMatrix(galeriList, numDofsPerNode, nodeMap, dofMap, A,
                     rowMap, colMap, revisedRowMap, revisedColMap,
                     rowImport, colImport, regionMats,
                     regionMatVecLIDs, quasiRegionCoordGIDs, regionInterfaceImporter,
                     rNodesPerDim, numLocalRegionNodes);

  // Build objects needed to construct the region coordinates
  RCP<Map> quasiRegCoordMap  = Teuchos::null;
  RCP<const Map> regCoordMap = Teuchos::null;
  RCP<Import> coordImporter  = Teuchos::null;

  quasiRegCoordMap = Xpetra::MapFactory<LO, GO, Node>::
      Build(nodeMap->lib(),
            Teuchos::OrdinalTraits<GO>::invalid(),
            quasiRegionCoordGIDs(),
            nodeMap->getIndexBase(),
            nodeMap->getComm());
  regCoordMap = Xpetra::MapFactory<LO, GO, Node>::
      Build(nodeMap->lib(),
            Teuchos::OrdinalTraits<GO>::invalid(),
            numLocalRegionNodes,
            nodeMap->getIndexBase(),
            nodeMap->getComm());

  coordImporter = ImportFactory::Build(nodeMap, quasiRegCoordMap);

  // create region coordinates vector
  regionCoordinates = Xpetra::MultiVectorFactory<real_type, LO, GO, NO>::Build(quasiRegCoordMap,
                                                                               coordinates->getNumVectors());
  regionCoordinates->doImport(*coordinates, *coordImporter, Xpetra::INSERT);
  regionCoordinates->replaceMap(regCoordMap);

  // Create regional nullspace and coordinates
  RCP<MultiVector> quasiRegionNullspace             = Teuchos::null;
  RCP<RealValuedMultiVector> quasiRegionCoordinates = Teuchos::null;

  compositeToRegional(nullspace, quasiRegionNullspace, regionNullspace,
                      revisedRowMap, rowImport);
  compositeToRegional(coordinates, quasiRegionCoordinates, regionCoordinates,
                      regCoordMap, coordImporter);

}  // createProblem

// This test aims at checking that RegionRFactory has a working constructor.
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RegionRFactory, RegionRFactCtor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  out << "version: " << MueLu::Version() << std::endl;

  RCP<RegionRFactory> myRFact = rcp(new RegionRFactory);
  TEST_EQUALITY(myRFact != Teuchos::null, true);

}  // RegionRFactCtor

// This test aims at checking that RegionRFactory has a working constructor.
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RegionRFactory, RegionRFactCtorKokkos, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  out << "version: " << MueLu::Version() << std::endl;

  RCP<RegionRFactory_kokkos> myRFact = rcp(new RegionRFactory_kokkos);
  TEST_EQUALITY(myRFact != Teuchos::null, true);

}  // RegionRFactCtor

// This test aims at checking that RegionRFactory produces a reasonable transfer
// operator on a simple Laplace 3D problem.
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RegionRFactory, RegionRFactLaplace3D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
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
  const int numRanks                  = comm->getSize();
  const int myRank                    = comm->getRank();

  const int numDimensions = 3;
  const LO numDofsPerNode = 1;
  GO nx = 7, ny = 7, nz = 4;
  Teuchos::Array<LO> lNodesPerDim({static_cast<LO>(nx), static_cast<LO>(ny), static_cast<LO>(nz)});
  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP();
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace3D");
  RCP<Matrix> regionMats                       = Teuchos::null;
  RCP<MultiVector> regionNullspace             = Teuchos::null;
  RCP<RealValuedMultiVector> regionCoordinates = Teuchos::null;
  RCP<Matrix> A                                = Teuchos::null;
  RCP<const Map> revisedRowMap                 = Teuchos::null;
  RCP<Import> rowImport                        = Teuchos::null;
  Teuchos::ArrayRCP<LocalOrdinal> regionMatVecLIDs;
  RCP<Import> regionInterfaceImporter = Teuchos::null;
  Teuchos::Array<LO> rNodesPerDim(3);

  createProblem(numDofsPerNode, galeriParameters, comm,
                A, regionMats, regionNullspace, regionCoordinates,
                revisedRowMap, rowImport,
                regionMatVecLIDs, regionInterfaceImporter, rNodesPerDim);

  // Generate levels for a two level hierarchy
  MueLu::Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);
  coarseLevel.SetFactoryManager(Teuchos::null);

  // Set requests and input data on fine level
  fineLevel.Request("A");
  fineLevel.Set("A", regionMats);
  fineLevel.Set("numDimensions", numDimensions);
  fineLevel.Set("lNodesPerDim", rNodesPerDim);
  fineLevel.Set("Nullspace", regionNullspace);
  fineLevel.Set("Coordinates", regionCoordinates);

  // Construct the region R factory
  RCP<RegionRFactory> myRFact                = rcp(new RegionRFactory);
  RCP<const Teuchos::ParameterList> myParams = myRFact->GetValidParameterList();

  // Set requests on coarse level to access
  // data generated by region R factory
  coarseLevel.Request("R", myRFact.get());  // request R
  coarseLevel.Request("P", myRFact.get());  // request P
  coarseLevel.Request("Nullspace", myRFact.get());
  coarseLevel.Request("Coordinates", myRFact.get());
  coarseLevel.Request(*myRFact);

  // Generate coarse level data with region R facotry
  myRFact->Build(fineLevel, coarseLevel);

  // Recover data from coarse level
  // and perform release mechanism
  // to free un-requested data.
  RCP<Matrix> R;
  coarseLevel.Get("R", R, myRFact.get());
  coarseLevel.Release("R", myRFact.get());
  TEST_EQUALITY(R != Teuchos::null, true);

  RCP<Matrix> P;
  coarseLevel.Get("P", P, myRFact.get());
  coarseLevel.Release("P", myRFact.get());
  TEST_EQUALITY(P != Teuchos::null, true);

  RCP<MultiVector> coarseNullspace;
  coarseLevel.Get("Nullspace", coarseNullspace, myRFact.get());
  coarseLevel.Release("Nullspace", myRFact.get());
  TEST_EQUALITY(coarseNullspace != Teuchos::null, true);

  RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<SC>::coordinateType, LO, GO, NO> > coarseCoordinates;
  coarseLevel.Get("Coordinates", coarseCoordinates, myRFact.get());
  coarseLevel.Release("Coordinates", myRFact.get());
  TEST_EQUALITY(coarseCoordinates != Teuchos::null, true);

  // R->describe(out, Teuchos::VERB_EXTREME);

  if (numRanks == 1) {
    TEST_EQUALITY(R->getGlobalNumRows(), 18);
    TEST_EQUALITY(R->getGlobalNumCols(), 196);
    TEST_EQUALITY(R->getLocalNumRows(), 18);
    TEST_EQUALITY(R->getCrsGraph()->getLocalNumCols(), 196);
    TEST_EQUALITY(R->getLocalNumEntries(), 726);

    Array<LO> rowLength = {{27, 45, 27, 45, 75, 45, 27, 45, 27,
                            27, 45, 27, 45, 75, 45, 27, 45, 27}};
    ArrayView<const LO> rowEntries;
    ArrayView<const SC> rowValues;
    for (int rowIdx = 0; rowIdx < static_cast<int>(R->getLocalNumRows()); ++rowIdx) {
      R->getLocalRowView(rowIdx, rowEntries, rowValues);
      TEST_EQUALITY(static_cast<LO>(rowEntries.size()), rowLength[rowIdx]);
    }

  } else {  // Running with 4 ranks
    TEST_EQUALITY(R->getGlobalNumRows(), 32);
    TEST_EQUALITY(R->getCrsGraph()->getLocalNumCols(), 64);
    TEST_EQUALITY(R->getLocalNumEntries(), 216);

    ArrayView<const LO> rowEntries;
    ArrayView<const SC> rowValues;
    if (myRank == 0) {
      TEST_EQUALITY(R->getLocalNumRows(), 8);
      TEST_EQUALITY(R->getCrsGraph()->getLocalNumCols(), 64);
      TEST_EQUALITY(R->getLocalNumEntries(), 216);

      Array<LO> rowLength = {{27, 27, 27, 27, 27, 27, 27, 27}};
      for (int rowIdx = 0; rowIdx < static_cast<int>(R->getLocalNumRows()); ++rowIdx) {
        R->getLocalRowView(rowIdx, rowEntries, rowValues);
        TEST_EQUALITY(static_cast<LO>(rowEntries.size()), rowLength[rowIdx]);
      }
    } else if (myRank == 1) {
      TEST_EQUALITY(R->getLocalNumRows(), 8);
      TEST_EQUALITY(R->getCrsGraph()->getLocalNumCols(), 64);
      TEST_EQUALITY(R->getLocalNumEntries(), 216);

      Array<LO> rowLength = {{27, 27, 27, 27, 27, 27, 27, 27}};
      for (int rowIdx = 0; rowIdx < static_cast<int>(R->getLocalNumRows()); ++rowIdx) {
        R->getLocalRowView(rowIdx, rowEntries, rowValues);
        TEST_EQUALITY(static_cast<LO>(rowEntries.size()), rowLength[rowIdx]);
      }
    } else if (myRank == 2) {
      TEST_EQUALITY(R->getLocalNumRows(), 8);
      TEST_EQUALITY(R->getCrsGraph()->getLocalNumCols(), 64);
      TEST_EQUALITY(R->getLocalNumEntries(), 216);

      Array<LO> rowLength = {{27, 27, 27, 27, 27, 27, 27, 27}};
      for (int rowIdx = 0; rowIdx < static_cast<int>(R->getLocalNumRows()); ++rowIdx) {
        R->getLocalRowView(rowIdx, rowEntries, rowValues);
        TEST_EQUALITY(static_cast<LO>(rowEntries.size()), rowLength[rowIdx]);
      }
    } else if (myRank == 3) {
      TEST_EQUALITY(R->getLocalNumRows(), 8);
      TEST_EQUALITY(R->getCrsGraph()->getLocalNumCols(), 64);
      TEST_EQUALITY(R->getLocalNumEntries(), 216);

      Array<LO> rowLength = {{27, 27, 27, 27, 27, 27, 27, 27}};
      for (int rowIdx = 0; rowIdx < static_cast<int>(R->getLocalNumRows()); ++rowIdx) {
        R->getLocalRowView(rowIdx, rowEntries, rowValues);
        TEST_EQUALITY(static_cast<LO>(rowEntries.size()), rowLength[rowIdx]);
      }
    }
  }

}  // RegionRFactLaplace3D

// This test aims at checking that RegionRFactory produces a reasonable transfer
// operator on a simple Elasticity 3D problem.
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RegionRFactory, RegionRFactElasticity3D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
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
  const int numRanks                  = comm->getSize();
  const int myRank                    = comm->getRank();

  const int numDimensions = 3;
  const LO numDofsPerNode = 3;
  GO nx = 7, ny = 7, nz = 4;
  Teuchos::Array<LO> lNodesPerDim({static_cast<LO>(nx), static_cast<LO>(ny), static_cast<LO>(nz)});
  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP();
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Elasticity3D");
  RCP<Matrix> regionMats                       = Teuchos::null;
  RCP<MultiVector> regionNullspace             = Teuchos::null;
  RCP<RealValuedMultiVector> regionCoordinates = Teuchos::null;
  RCP<Matrix> A                                = Teuchos::null;
  RCP<const Map> revisedRowMap                 = Teuchos::null;
  RCP<Import> rowImport                        = Teuchos::null;
  Teuchos::ArrayRCP<LocalOrdinal> regionMatVecLIDs;
  RCP<Import> regionInterfaceImporter;
  Teuchos::Array<LO> rNodesPerDim(4);

  createProblem(numDofsPerNode, galeriParameters, comm,
                A, regionMats, regionNullspace, regionCoordinates,
                revisedRowMap, rowImport,
                regionMatVecLIDs, regionInterfaceImporter, rNodesPerDim);

  // Generate levels for a two level hierarchy
  MueLu::Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);
  coarseLevel.SetFactoryManager(Teuchos::null);

  // Set requests and input data on fine level
  fineLevel.Request("A");
  fineLevel.Set("A", regionMats);
  fineLevel.Set("numDimensions", numDimensions);
  fineLevel.Set("lNodesPerDim", rNodesPerDim);
  fineLevel.Set("Nullspace", regionNullspace);
  fineLevel.Set("Coordinates", regionCoordinates);

  // Construct the region R factory
  RCP<RegionRFactory> myRFact                = rcp(new RegionRFactory);
  RCP<const Teuchos::ParameterList> myParams = myRFact->GetValidParameterList();

  // Set requests on coarse level to access
  // data generated by region R factory
  coarseLevel.Request("R", myRFact.get());  // request R
  coarseLevel.Request("P", myRFact.get());  // request P
  coarseLevel.Request("Nullspace", myRFact.get());
  coarseLevel.Request("Coordinates", myRFact.get());
  coarseLevel.Request(*myRFact);

  // Generate coarse level data with region R facotry
  myRFact->Build(fineLevel, coarseLevel);

  // Recover data from coarse level
  // and perform release mechanism
  // to free un-requested data.
  RCP<Matrix> R;
  coarseLevel.Get("R", R, myRFact.get());
  coarseLevel.Release("R", myRFact.get());
  TEST_EQUALITY(R != Teuchos::null, true);

  RCP<Matrix> P;
  coarseLevel.Get("P", P, myRFact.get());
  coarseLevel.Release("P", myRFact.get());
  TEST_EQUALITY(P != Teuchos::null, true);

  RCP<MultiVector> coarseNullspace;
  coarseLevel.Get("Nullspace", coarseNullspace, myRFact.get());
  coarseLevel.Release("Nullspace", myRFact.get());
  TEST_EQUALITY(coarseNullspace != Teuchos::null, true);

  RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<SC>::coordinateType, LO, GO, NO> > coarseCoordinates;
  coarseLevel.Get("Coordinates", coarseCoordinates, myRFact.get());
  coarseLevel.Release("Coordinates", myRFact.get());
  TEST_EQUALITY(coarseCoordinates != Teuchos::null, true);

  // R->describe(out, Teuchos::VERB_EXTREME);

  if (numRanks == 1) {
    TEST_EQUALITY(R->getGlobalNumRows(), 54);
    TEST_EQUALITY(R->getGlobalNumCols(), 588);
    TEST_EQUALITY(R->getLocalNumRows(), 54);
    TEST_EQUALITY(R->getCrsGraph()->getLocalNumCols(), 588);
    TEST_EQUALITY(R->getLocalNumEntries(), 2178);

    Array<LO> rowLength = {{27, 27, 27, 45, 45, 45, 27, 27, 27, 45, 45, 45,
                            75, 75, 75, 45, 45, 45, 27, 27, 27, 45, 45, 45,
                            27, 27, 27, 27, 27, 27, 45, 45, 45, 27, 27, 27,
                            45, 45, 45, 75, 75, 75, 45, 45, 45, 27, 27, 27,
                            45, 45, 45, 27, 27, 27}};
    ArrayView<const LO> rowEntries;
    ArrayView<const SC> rowValues;
    for (int rowIdx = 0; rowIdx < static_cast<int>(R->getLocalNumRows()); ++rowIdx) {
      R->getLocalRowView(rowIdx, rowEntries, rowValues);
      TEST_EQUALITY(static_cast<LO>(rowEntries.size()), rowLength[rowIdx]);
    }

  } else {  // Running with 4 ranks
    TEST_EQUALITY(R->getGlobalNumRows(), 96);
    TEST_EQUALITY(R->getCrsGraph()->getLocalNumCols(), 192);
    TEST_EQUALITY(R->getLocalNumEntries(), 648);

    ArrayView<const LO> rowEntries;
    ArrayView<const SC> rowValues;
    if (myRank == 0) {
      TEST_EQUALITY(R->getLocalNumRows(), 24);
      TEST_EQUALITY(R->getCrsGraph()->getLocalNumCols(), 192);
      TEST_EQUALITY(R->getLocalNumEntries(), 648);

      Array<LO> rowLength = {{27, 27, 27, 27, 27, 27, 27, 27,
                              27, 27, 27, 27, 27, 27, 27, 27,
                              27, 27, 27, 27, 27, 27, 27, 27}};
      for (int rowIdx = 0; rowIdx < static_cast<int>(R->getLocalNumRows()); ++rowIdx) {
        R->getLocalRowView(rowIdx, rowEntries, rowValues);
        TEST_EQUALITY(static_cast<LO>(rowEntries.size()), rowLength[rowIdx]);
      }
    } else if (myRank == 1) {
      TEST_EQUALITY(R->getLocalNumRows(), 24);
      TEST_EQUALITY(R->getCrsGraph()->getLocalNumCols(), 192);
      TEST_EQUALITY(R->getLocalNumEntries(), 648);

      Array<LO> rowLength = {{27, 27, 27, 27, 27, 27, 27, 27,
                              27, 27, 27, 27, 27, 27, 27, 27,
                              27, 27, 27, 27, 27, 27, 27, 27}};
      for (int rowIdx = 0; rowIdx < static_cast<int>(R->getLocalNumRows()); ++rowIdx) {
        R->getLocalRowView(rowIdx, rowEntries, rowValues);
        TEST_EQUALITY(static_cast<LO>(rowEntries.size()), rowLength[rowIdx]);
      }
    } else if (myRank == 2) {
      TEST_EQUALITY(R->getLocalNumRows(), 24);
      TEST_EQUALITY(R->getCrsGraph()->getLocalNumCols(), 192);
      TEST_EQUALITY(R->getLocalNumEntries(), 648);

      Array<LO> rowLength = {{27, 27, 27, 27, 27, 27, 27, 27,
                              27, 27, 27, 27, 27, 27, 27, 27,
                              27, 27, 27, 27, 27, 27, 27, 27}};
      for (int rowIdx = 0; rowIdx < static_cast<int>(R->getLocalNumRows()); ++rowIdx) {
        R->getLocalRowView(rowIdx, rowEntries, rowValues);
        TEST_EQUALITY(static_cast<LO>(rowEntries.size()), rowLength[rowIdx]);
      }
    } else if (myRank == 3) {
      TEST_EQUALITY(R->getLocalNumRows(), 24);
      TEST_EQUALITY(R->getCrsGraph()->getLocalNumCols(), 192);
      TEST_EQUALITY(R->getLocalNumEntries(), 648);

      Array<LO> rowLength = {{27, 27, 27, 27, 27, 27, 27, 27,
                              27, 27, 27, 27, 27, 27, 27, 27,
                              27, 27, 27, 27, 27, 27, 27, 27}};
      for (int rowIdx = 0; rowIdx < static_cast<int>(R->getLocalNumRows()); ++rowIdx) {
        R->getLocalRowView(rowIdx, rowEntries, rowValues);
        TEST_EQUALITY(static_cast<LO>(rowEntries.size()), rowLength[rowIdx]);
      }
    }
  }

}  // RegionRFactElasticity3D

// This test aims at checking that RegionRFactory produces a reasonable transfer
// operator on a simple Laplace 3D problem.
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RegionRFactory, RegionRFactLaplace3DKokkos, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
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
  const int numRanks                  = comm->getSize();
  const int myRank                    = comm->getRank();

  const int numDimensions = 3;
  const LO numDofsPerNode = 1;
  GO nx = 7, ny = 7, nz = 4;
  Teuchos::Array<LO> lNodesPerDim({static_cast<LO>(nx), static_cast<LO>(ny), static_cast<LO>(nz)});
  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP();
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace3D");
  RCP<Matrix> regionMats                       = Teuchos::null;
  RCP<MultiVector> regionNullspace             = Teuchos::null;
  RCP<RealValuedMultiVector> regionCoordinates = Teuchos::null;
  RCP<Matrix> A                                = Teuchos::null;
  RCP<const Map> revisedRowMap                 = Teuchos::null;
  RCP<Import> rowImport                        = Teuchos::null;
  Teuchos::ArrayRCP<LocalOrdinal> regionMatVecLIDs;
  RCP<Import> regionInterfaceImporter = Teuchos::null;
  Teuchos::Array<LO> rNodesPerDim(3);

  createProblem(numDofsPerNode, galeriParameters, comm,
                A, regionMats, regionNullspace, regionCoordinates,
                revisedRowMap, rowImport,
                regionMatVecLIDs, regionInterfaceImporter, rNodesPerDim);

  // Generate levels for a two level hierarchy
  MueLu::Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);
  coarseLevel.SetFactoryManager(Teuchos::null);

  // Set requests and input data on fine level
  fineLevel.Request("A");
  fineLevel.Set("A", regionMats);
  fineLevel.Set("numDimensions", numDimensions);
  fineLevel.Set("lNodesPerDim", rNodesPerDim);
  fineLevel.Set("Nullspace", regionNullspace);
  fineLevel.Set("Coordinates", regionCoordinates);

  // Construct the region R factory
  RCP<RegionRFactory_kokkos> myRFact         = rcp(new RegionRFactory_kokkos);
  RCP<const Teuchos::ParameterList> myParams = myRFact->GetValidParameterList();

  // Set requests on coarse level to access
  // data generated by region R factory
  coarseLevel.Request("R", myRFact.get());  // request R
  coarseLevel.Request("P", myRFact.get());  // request P
  coarseLevel.Request("Nullspace", myRFact.get());
  coarseLevel.Request("Coordinates", myRFact.get());
  coarseLevel.Request(*myRFact);

  // Generate coarse level data with region R facotry
  myRFact->Build(fineLevel, coarseLevel);

  // Recover data from coarse level
  // and perform release mechanism
  // to free un-requested data.
  RCP<Matrix> R;
  coarseLevel.Get("R", R, myRFact.get());
  coarseLevel.Release("R", myRFact.get());
  TEST_EQUALITY(R != Teuchos::null, true);

  RCP<Matrix> P;
  coarseLevel.Get("P", P, myRFact.get());
  coarseLevel.Release("P", myRFact.get());
  TEST_EQUALITY(P != Teuchos::null, true);

  RCP<MultiVector> coarseNullspace;
  coarseLevel.Get("Nullspace", coarseNullspace, myRFact.get());
  coarseLevel.Release("Nullspace", myRFact.get());
  TEST_EQUALITY(coarseNullspace != Teuchos::null, true);

  RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<SC>::coordinateType, LO, GO, NO> > coarseCoordinates;
  coarseLevel.Get("Coordinates", coarseCoordinates, myRFact.get());
  coarseLevel.Release("Coordinates", myRFact.get());
  TEST_EQUALITY(coarseCoordinates != Teuchos::null, true);

  // R->describe(out, Teuchos::VERB_EXTREME);

  if (numRanks == 1) {
    TEST_EQUALITY(R->getGlobalNumRows(), 18);
    TEST_EQUALITY(R->getGlobalNumCols(), 196);
    TEST_EQUALITY(R->getLocalNumRows(), 18);
    TEST_EQUALITY(R->getCrsGraph()->getLocalNumCols(), 196);
    TEST_EQUALITY(R->getLocalNumEntries(), 726);

    Array<LO> rowLength = {{27, 45, 27, 45, 75, 45, 27, 45, 27,
                            27, 45, 27, 45, 75, 45, 27, 45, 27}};
    ArrayView<const LO> rowEntries;
    ArrayView<const SC> rowValues;
    for (int rowIdx = 0; rowIdx < static_cast<int>(R->getLocalNumRows()); ++rowIdx) {
      R->getLocalRowView(rowIdx, rowEntries, rowValues);
      TEST_EQUALITY(static_cast<LO>(rowEntries.size()), rowLength[rowIdx]);
    }

  } else {  // Running with 4 ranks
    TEST_EQUALITY(R->getGlobalNumRows(), 32);
    TEST_EQUALITY(R->getCrsGraph()->getLocalNumCols(), 64);
    TEST_EQUALITY(R->getLocalNumEntries(), 216);

    ArrayView<const LO> rowEntries;
    ArrayView<const SC> rowValues;
    if (myRank == 0) {
      TEST_EQUALITY(R->getLocalNumRows(), 8);
      TEST_EQUALITY(R->getCrsGraph()->getLocalNumCols(), 64);
      TEST_EQUALITY(R->getLocalNumEntries(), 216);

      Array<LO> rowLength = {{27, 27, 27, 27, 27, 27, 27, 27}};
      for (int rowIdx = 0; rowIdx < static_cast<int>(R->getLocalNumRows()); ++rowIdx) {
        R->getLocalRowView(rowIdx, rowEntries, rowValues);
        TEST_EQUALITY(static_cast<LO>(rowEntries.size()), rowLength[rowIdx]);
      }
    } else if (myRank == 1) {
      TEST_EQUALITY(R->getLocalNumRows(), 8);
      TEST_EQUALITY(R->getCrsGraph()->getLocalNumCols(), 64);
      TEST_EQUALITY(R->getLocalNumEntries(), 216);

      Array<LO> rowLength = {{27, 27, 27, 27, 27, 27, 27, 27}};
      for (int rowIdx = 0; rowIdx < static_cast<int>(R->getLocalNumRows()); ++rowIdx) {
        R->getLocalRowView(rowIdx, rowEntries, rowValues);
        TEST_EQUALITY(static_cast<LO>(rowEntries.size()), rowLength[rowIdx]);
      }
    } else if (myRank == 2) {
      TEST_EQUALITY(R->getLocalNumRows(), 8);
      TEST_EQUALITY(R->getCrsGraph()->getLocalNumCols(), 64);
      TEST_EQUALITY(R->getLocalNumEntries(), 216);

      Array<LO> rowLength = {{27, 27, 27, 27, 27, 27, 27, 27}};
      for (int rowIdx = 0; rowIdx < static_cast<int>(R->getLocalNumRows()); ++rowIdx) {
        R->getLocalRowView(rowIdx, rowEntries, rowValues);
        TEST_EQUALITY(static_cast<LO>(rowEntries.size()), rowLength[rowIdx]);
      }
    } else if (myRank == 3) {
      TEST_EQUALITY(R->getLocalNumRows(), 8);
      TEST_EQUALITY(R->getCrsGraph()->getLocalNumCols(), 64);
      TEST_EQUALITY(R->getLocalNumEntries(), 216);

      Array<LO> rowLength = {{27, 27, 27, 27, 27, 27, 27, 27}};
      for (int rowIdx = 0; rowIdx < static_cast<int>(R->getLocalNumRows()); ++rowIdx) {
        R->getLocalRowView(rowIdx, rowEntries, rowValues);
        TEST_EQUALITY(static_cast<LO>(rowEntries.size()), rowLength[rowIdx]);
      }
    }
  }

}  // RegionRFactLaplace3DKokkos

// This test aims at checking that RegionRFactory produces a reasonable transfer
// operator on a simple Elasticity 3D problem.
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RegionRFactory, RegionRFactElasticity3DKokkos, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
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
  const int numRanks                  = comm->getSize();
  const int myRank                    = comm->getRank();

  const int numDimensions = 3;
  const LO numDofsPerNode = 3;
  GO nx = 7, ny = 7, nz = 4;
  Teuchos::Array<LO> lNodesPerDim({static_cast<LO>(nx), static_cast<LO>(ny), static_cast<LO>(nz)});
  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP();
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Elasticity3D");
  RCP<Matrix> regionMats                       = Teuchos::null;
  RCP<MultiVector> regionNullspace             = Teuchos::null;
  RCP<RealValuedMultiVector> regionCoordinates = Teuchos::null;
  RCP<Matrix> A                                = Teuchos::null;
  RCP<const Map> revisedRowMap                 = Teuchos::null;
  RCP<Import> rowImport                        = Teuchos::null;
  Teuchos::ArrayRCP<LocalOrdinal> regionMatVecLIDs;
  RCP<Import> regionInterfaceImporter;
  Teuchos::Array<LO> rNodesPerDim(4);

  createProblem(numDofsPerNode, galeriParameters, comm,
                A, regionMats, regionNullspace, regionCoordinates,
                revisedRowMap, rowImport,
                regionMatVecLIDs, regionInterfaceImporter, rNodesPerDim);

  // Generate levels for a two level hierarchy
  MueLu::Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);
  coarseLevel.SetFactoryManager(Teuchos::null);

  // Set requests and input data on fine level
  fineLevel.Request("A");
  fineLevel.Set("A", regionMats);
  fineLevel.Set("numDimensions", numDimensions);
  fineLevel.Set("lNodesPerDim", rNodesPerDim);
  fineLevel.Set("Nullspace", regionNullspace);
  fineLevel.Set("Coordinates", regionCoordinates);

  // Construct the region R factory
  RCP<RegionRFactory_kokkos> myRFact         = rcp(new RegionRFactory_kokkos);
  RCP<const Teuchos::ParameterList> myParams = myRFact->GetValidParameterList();

  // Set requests on coarse level to access
  // data generated by region R factory
  coarseLevel.Request("R", myRFact.get());  // request R
  coarseLevel.Request("P", myRFact.get());  // request P
  coarseLevel.Request("Nullspace", myRFact.get());
  coarseLevel.Request("Coordinates", myRFact.get());
  coarseLevel.Request(*myRFact);

  // Generate coarse level data with region R facotry
  myRFact->Build(fineLevel, coarseLevel);

  // Recover data from coarse level
  // and perform release mechanism
  // to free un-requested data.
  RCP<Matrix> R;
  coarseLevel.Get("R", R, myRFact.get());
  coarseLevel.Release("R", myRFact.get());
  TEST_EQUALITY(R != Teuchos::null, true);

  RCP<Matrix> P;
  coarseLevel.Get("P", P, myRFact.get());
  coarseLevel.Release("P", myRFact.get());
  TEST_EQUALITY(P != Teuchos::null, true);

  RCP<MultiVector> coarseNullspace;
  coarseLevel.Get("Nullspace", coarseNullspace, myRFact.get());
  coarseLevel.Release("Nullspace", myRFact.get());
  TEST_EQUALITY(coarseNullspace != Teuchos::null, true);

  RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<SC>::coordinateType, LO, GO, NO> > coarseCoordinates;
  coarseLevel.Get("Coordinates", coarseCoordinates, myRFact.get());
  coarseLevel.Release("Coordinates", myRFact.get());
  TEST_EQUALITY(coarseCoordinates != Teuchos::null, true);

  // R->describe(out, Teuchos::VERB_EXTREME);

  if (numRanks == 1) {
    TEST_EQUALITY(R->getGlobalNumRows(), 54);
    TEST_EQUALITY(R->getGlobalNumCols(), 588);
    TEST_EQUALITY(R->getLocalNumRows(), 54);
    TEST_EQUALITY(R->getCrsGraph()->getLocalNumCols(), 588);
    TEST_EQUALITY(R->getLocalNumEntries(), 2178);

    Array<LO> rowLength = {{27, 27, 27, 45, 45, 45, 27, 27, 27, 45, 45, 45,
                            75, 75, 75, 45, 45, 45, 27, 27, 27, 45, 45, 45,
                            27, 27, 27, 27, 27, 27, 45, 45, 45, 27, 27, 27,
                            45, 45, 45, 75, 75, 75, 45, 45, 45, 27, 27, 27,
                            45, 45, 45, 27, 27, 27}};
    ArrayView<const LO> rowEntries;
    ArrayView<const SC> rowValues;
    for (int rowIdx = 0; rowIdx < static_cast<int>(R->getLocalNumRows()); ++rowIdx) {
      R->getLocalRowView(rowIdx, rowEntries, rowValues);
      TEST_EQUALITY(static_cast<LO>(rowEntries.size()), rowLength[rowIdx]);
    }

  } else {  // Running with 4 ranks
    TEST_EQUALITY(R->getGlobalNumRows(), 96);
    TEST_EQUALITY(R->getCrsGraph()->getLocalNumCols(), 192);
    TEST_EQUALITY(R->getLocalNumEntries(), 648);

    ArrayView<const LO> rowEntries;
    ArrayView<const SC> rowValues;
    if (myRank == 0) {
      TEST_EQUALITY(R->getLocalNumRows(), 24);
      TEST_EQUALITY(R->getCrsGraph()->getLocalNumCols(), 192);
      TEST_EQUALITY(R->getLocalNumEntries(), 648);

      Array<LO> rowLength = {{27, 27, 27, 27, 27, 27, 27, 27,
                              27, 27, 27, 27, 27, 27, 27, 27,
                              27, 27, 27, 27, 27, 27, 27, 27}};
      for (int rowIdx = 0; rowIdx < static_cast<int>(R->getLocalNumRows()); ++rowIdx) {
        R->getLocalRowView(rowIdx, rowEntries, rowValues);
        TEST_EQUALITY(static_cast<LO>(rowEntries.size()), rowLength[rowIdx]);
      }
    } else if (myRank == 1) {
      TEST_EQUALITY(R->getLocalNumRows(), 24);
      TEST_EQUALITY(R->getCrsGraph()->getLocalNumCols(), 192);
      TEST_EQUALITY(R->getLocalNumEntries(), 648);

      Array<LO> rowLength = {{27, 27, 27, 27, 27, 27, 27, 27,
                              27, 27, 27, 27, 27, 27, 27, 27,
                              27, 27, 27, 27, 27, 27, 27, 27}};
      for (int rowIdx = 0; rowIdx < static_cast<int>(R->getLocalNumRows()); ++rowIdx) {
        R->getLocalRowView(rowIdx, rowEntries, rowValues);
        TEST_EQUALITY(static_cast<LO>(rowEntries.size()), rowLength[rowIdx]);
      }
    } else if (myRank == 2) {
      TEST_EQUALITY(R->getLocalNumRows(), 24);
      TEST_EQUALITY(R->getCrsGraph()->getLocalNumCols(), 192);
      TEST_EQUALITY(R->getLocalNumEntries(), 648);

      Array<LO> rowLength = {{27, 27, 27, 27, 27, 27, 27, 27,
                              27, 27, 27, 27, 27, 27, 27, 27,
                              27, 27, 27, 27, 27, 27, 27, 27}};
      for (int rowIdx = 0; rowIdx < static_cast<int>(R->getLocalNumRows()); ++rowIdx) {
        R->getLocalRowView(rowIdx, rowEntries, rowValues);
        TEST_EQUALITY(static_cast<LO>(rowEntries.size()), rowLength[rowIdx]);
      }
    } else if (myRank == 3) {
      TEST_EQUALITY(R->getLocalNumRows(), 24);
      TEST_EQUALITY(R->getCrsGraph()->getLocalNumCols(), 192);
      TEST_EQUALITY(R->getLocalNumEntries(), 648);

      Array<LO> rowLength = {{27, 27, 27, 27, 27, 27, 27, 27,
                              27, 27, 27, 27, 27, 27, 27, 27,
                              27, 27, 27, 27, 27, 27, 27, 27}};
      for (int rowIdx = 0; rowIdx < static_cast<int>(R->getLocalNumRows()); ++rowIdx) {
        R->getLocalRowView(rowIdx, rowEntries, rowValues);
        TEST_EQUALITY(static_cast<LO>(rowEntries.size()), rowLength[rowIdx]);
      }
    }
  }

}  // RegionRFactElasticity3DKokkos

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node)                                                            \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RegionRFactory, RegionRFactCtor, Scalar, LO, GO, Node)            \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RegionRFactory, RegionRFactCtorKokkos, Scalar, LO, GO, Node)      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RegionRFactory, RegionRFactLaplace3D, Scalar, LO, GO, Node)       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RegionRFactory, RegionRFactLaplace3DKokkos, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RegionRFactory, RegionRFactElasticity3D, Scalar, LO, GO, Node)    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RegionRFactory, RegionRFactElasticity3DKokkos, Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
