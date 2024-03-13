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
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <MueLu_FactoryManager.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_MapTransferFactory.hpp>
#include <MueLu_NoFactory.hpp>
#include <MueLu_TentativePFactory.hpp>

#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MapTransferFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  out << "version: " << MueLu::Version() << std::endl;

  RCP<MapTransferFactory> mapTransferFactory = rcp(new MapTransferFactory());
  TEST_ASSERT(!mapTransferFactory.is_null());
}  // Constructor

/* This tests coarsens the row map of the fine level operator, so the result from the MapTransferFactory
 * needs to match the domain map of the prolongator.
 *
 * Assume a 1D Poisson discretization.
 */
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MapTransferFactory, TransferFullMap1D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using test_factory = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;
  out << "Test transfer of a map with the MapTransferFactory" << std::endl;

  // Manual setup of a two-level hierarchy
  RCP<Level> fineLevel   = rcp(new Level());
  RCP<Level> coarseLevel = rcp(new Level());
  coarseLevel->SetPreviousLevel(fineLevel);
  fineLevel->SetLevelID(0);
  coarseLevel->SetLevelID(1);

  TEST_EQUALITY_CONST(fineLevel->GetLevelID(), 0);
  TEST_EQUALITY_CONST(coarseLevel->GetLevelID(), 1);

  // Create a dummy matrix needed to build a prolongator
  const GO nx   = 199;
  RCP<Matrix> A = test_factory::Build1DPoisson(nx);
  fineLevel->Set("A", A);

  const std::string mapName = "Dummy Map";
  fineLevel->Set(mapName, A->getRowMap());

  TEST_ASSERT(fineLevel->IsAvailable("A", MueLu::NoFactory::get()));
  TEST_ASSERT(fineLevel->IsAvailable(mapName, MueLu::NoFactory::get()));

  RCP<FactoryManager> factoryManager = rcp(new FactoryManager());
  factoryManager->SetKokkosRefactor(false);
  factoryManager->SetFactory(mapName, MueLu::NoFactory::getRCP());
  fineLevel->SetFactoryManager(factoryManager);
  coarseLevel->SetFactoryManager(factoryManager);

  RCP<TentativePFactory> tentativePFact = rcp(new TentativePFactory());

  RCP<MapTransferFactory> mapTransferFactory = rcp(new MapTransferFactory());
  mapTransferFactory->SetParameter("map: factory", Teuchos::ParameterEntry(mapName));
  mapTransferFactory->SetParameter("map: name", Teuchos::ParameterEntry(mapName));
  mapTransferFactory->SetFactory("P", tentativePFact);

  coarseLevel->Request(mapName, MueLu::NoFactory::get());
  coarseLevel->Request("P", tentativePFact.get(), mapTransferFactory.get());
  coarseLevel->Request(*mapTransferFactory);  // This calls DeclareInput() on mapTransferFactory

  TEST_ASSERT(coarseLevel->IsRequested(mapName, MueLu::NoFactory::get()));
  TEST_ASSERT(coarseLevel->IsRequested("P", tentativePFact.get()));

  RCP<Matrix> Ptent = coarseLevel->Get<RCP<Matrix>>("P", tentativePFact.get());
  TEST_ASSERT(!Ptent.is_null());

  TEST_ASSERT(!coarseLevel->IsAvailable(mapName, MueLu::NoFactory::get()));
  mapTransferFactory->Build(*fineLevel, *coarseLevel);
  TEST_ASSERT(coarseLevel->IsAvailable(mapName, MueLu::NoFactory::get()));

  RCP<const Map> coarsenedMap = coarseLevel->Get<RCP<const Map>>(mapName, MueLu::NoFactory::get());
  TEST_ASSERT(!coarsenedMap.is_null());

  TEST_ASSERT(coarsenedMap->isSameAs(*Ptent->getDomainMap()));
}

/* This tests coarsens the row map of the fine level operator, so the result from the MapTransferFactory
 * needs to match the domain map of the prolongator.
 *
 * Assume a 2D Poisson discretization.
 */
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MapTransferFactory, TransferFullMap2D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using test_factory = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;
  out << "Test transfer of a map with the MapTransferFactory" << std::endl;

  // Manual setup of a two-level hierarchy
  RCP<Level> fineLevel   = rcp(new Level());
  RCP<Level> coarseLevel = rcp(new Level());
  coarseLevel->SetPreviousLevel(fineLevel);
  fineLevel->SetLevelID(0);
  coarseLevel->SetLevelID(1);

  TEST_EQUALITY_CONST(fineLevel->GetLevelID(), 0);
  TEST_EQUALITY_CONST(coarseLevel->GetLevelID(), 1);

  // Create a dummy matrix needed to build a prolongator
  const GO nx   = 19;
  const GO ny   = 17;
  RCP<Matrix> A = test_factory::Build2DPoisson(nx, ny);
  fineLevel->Set("A", A);

  const std::string mapName = "Dummy Map";
  fineLevel->Set(mapName, A->getRowMap());

  TEST_ASSERT(fineLevel->IsAvailable("A", MueLu::NoFactory::get()));
  TEST_ASSERT(fineLevel->IsAvailable(mapName, MueLu::NoFactory::get()));

  RCP<FactoryManager> factoryManager = rcp(new FactoryManager());
  factoryManager->SetKokkosRefactor(false);
  factoryManager->SetFactory(mapName, MueLu::NoFactory::getRCP());
  fineLevel->SetFactoryManager(factoryManager);
  coarseLevel->SetFactoryManager(factoryManager);

  RCP<TentativePFactory> tentativePFact = rcp(new TentativePFactory());

  RCP<MapTransferFactory> mapTransferFactory = rcp(new MapTransferFactory());
  mapTransferFactory->SetParameter("map: factory", Teuchos::ParameterEntry(mapName));
  mapTransferFactory->SetParameter("map: name", Teuchos::ParameterEntry(mapName));
  mapTransferFactory->SetFactory("P", tentativePFact);

  coarseLevel->Request(mapName, MueLu::NoFactory::get());
  coarseLevel->Request("P", tentativePFact.get(), mapTransferFactory.get());
  coarseLevel->Request(*mapTransferFactory);  // This calls DeclareInput() on mapTransferFactory

  TEST_ASSERT(coarseLevel->IsRequested(mapName, MueLu::NoFactory::get()));
  TEST_ASSERT(coarseLevel->IsRequested("P", tentativePFact.get()));

  RCP<Matrix> Ptent = coarseLevel->Get<RCP<Matrix>>("P", tentativePFact.get());
  TEST_ASSERT(!Ptent.is_null());

  TEST_ASSERT(!coarseLevel->IsAvailable(mapName, MueLu::NoFactory::get()));
  mapTransferFactory->Build(*fineLevel, *coarseLevel);
  TEST_ASSERT(coarseLevel->IsAvailable(mapName, MueLu::NoFactory::get()));

  RCP<const Map> coarsenedMap = coarseLevel->Get<RCP<const Map>>(mapName, MueLu::NoFactory::get());
  TEST_ASSERT(!coarsenedMap.is_null());

  TEST_ASSERT(coarsenedMap->isSameAs(*Ptent->getDomainMap()));
}

/* This tests coarsens the row map of the fine level operator, so the result from the MapTransferFactory
 * needs to match the domain map of the prolongator.
 *
 * Assume a 3D elasticity discretization.
 */
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MapTransferFactory, TransferFullMap3DElasticity, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using test_factory = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;
  out << "Test transfer of a map with the MapTransferFactory" << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();

  // Manual setup of a two-level hierarchy
  RCP<Level> fineLevel   = rcp(new Level());
  RCP<Level> coarseLevel = rcp(new Level());
  coarseLevel->SetPreviousLevel(fineLevel);
  fineLevel->SetLevelID(0);
  coarseLevel->SetLevelID(1);

  TEST_EQUALITY_CONST(fineLevel->GetLevelID(), 0);
  TEST_EQUALITY_CONST(coarseLevel->GetLevelID(), 1);

  // Create a 3D elsasticity matrix needed to build a prolongator
  RCP<Matrix> A = Teuchos::null;
  {
    const GO nx                  = 5;
    const GO ny                  = 5;
    const GO nz                  = 5;
    const int numDofsPerNode     = 3;  // 3D elasticity
    const std::string matrixType = "Elasticity3D";

    Teuchos::ParameterList galeriList;
    galeriList.set("nx", nx);
    galeriList.set("ny", ny);
    galeriList.set("nz", nz);
    galeriList.set("matrixType", matrixType);
    RCP<const Map> nodeMap = Galeri::Xpetra::CreateMap<LO, GO, Node>(TestHelpers::Parameters::getLib(), "Cartesian3D", TestHelpers::Parameters::getDefaultComm(), galeriList);
    RCP<const Map> dofMap  = Xpetra::MapFactory<LO, GO, Node>::Build(nodeMap, numDofsPerNode);
    RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector>> Pr =
        Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(matrixType, dofMap, galeriList);
    A = Pr->BuildMatrix();

    TEST_ASSERT(!A.is_null());

    fineLevel->Set("A", A);
  }

  const std::string mapName = "Dummy Map";
  fineLevel->Set(mapName, A->getRowMap());

  TEST_ASSERT(fineLevel->IsAvailable("A", MueLu::NoFactory::get()));
  TEST_ASSERT(fineLevel->IsAvailable(mapName, MueLu::NoFactory::get()));

  RCP<FactoryManager> factoryManager = rcp(new FactoryManager());
  factoryManager->SetKokkosRefactor(false);
  factoryManager->SetFactory(mapName, MueLu::NoFactory::getRCP());
  fineLevel->SetFactoryManager(factoryManager);
  coarseLevel->SetFactoryManager(factoryManager);

  RCP<TentativePFactory> tentativePFact = rcp(new TentativePFactory());

  RCP<MapTransferFactory> mapTransferFactory = rcp(new MapTransferFactory());
  mapTransferFactory->SetParameter("map: factory", Teuchos::ParameterEntry(mapName));
  mapTransferFactory->SetParameter("map: name", Teuchos::ParameterEntry(mapName));
  mapTransferFactory->SetFactory("P", tentativePFact);

  coarseLevel->Request(mapName, MueLu::NoFactory::get());
  coarseLevel->Request("P", tentativePFact.get(), mapTransferFactory.get());
  coarseLevel->Request(*mapTransferFactory);  // This calls DeclareInput() on mapTransferFactory

  TEST_ASSERT(coarseLevel->IsRequested(mapName, MueLu::NoFactory::get()));
  TEST_ASSERT(coarseLevel->IsRequested("P", tentativePFact.get()));

  RCP<Matrix> Ptent = coarseLevel->Get<RCP<Matrix>>("P", tentativePFact.get());
  TEST_ASSERT(!Ptent.is_null());

  TEST_ASSERT(!coarseLevel->IsAvailable(mapName, MueLu::NoFactory::get()));
  mapTransferFactory->Build(*fineLevel, *coarseLevel);
  TEST_ASSERT(coarseLevel->IsAvailable(mapName, MueLu::NoFactory::get()));

  RCP<const Map> coarsenedMap = coarseLevel->Get<RCP<const Map>>(mapName, MueLu::NoFactory::get());
  TEST_ASSERT(!coarsenedMap.is_null());

  TEST_ASSERT(coarsenedMap->isSameAs(*Ptent->getDomainMap()));
}

/* This tests coarsens the row map of the fine level operator, so the result from the MapTransferFactory
 * needs to match the domain map of the prolongator. This particular test targets the special case of
 * limiting the number of columns from P to be used for the transfer (e.g. when excluding the rotations
 * from the nullspace when creating the MapTransfer.)
 *
 * Assume a 3D elasticity discretization.
 */
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MapTransferFactory, TransferFullMap3DElasticityReducedP, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using test_factory = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;
  out << "Test transfer of a map with the MapTransferFactory" << std::endl;

  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();

  // Manual setup of a two-level hierarchy
  RCP<Level> fineLevel   = rcp(new Level());
  RCP<Level> coarseLevel = rcp(new Level());
  coarseLevel->SetPreviousLevel(fineLevel);
  fineLevel->SetLevelID(0);
  coarseLevel->SetLevelID(1);

  TEST_EQUALITY_CONST(fineLevel->GetLevelID(), 0);
  TEST_EQUALITY_CONST(coarseLevel->GetLevelID(), 1);

  // Create a 3D elasticity matrix and its nullspace needed to build a prolongator
  RCP<Matrix> A = Teuchos::null;
  {
    const GO nx                  = 10;
    const GO ny                  = 10;
    const GO nz                  = 5;
    const int numDofsPerNode     = 3;  // 3D elasticity
    const std::string matrixType = "Elasticity3D";

    Teuchos::ParameterList galeriList;
    galeriList.set("nx", nx);
    galeriList.set("ny", ny);
    galeriList.set("nz", nz);
    galeriList.set("matrixType", matrixType);
    RCP<const Map> nodeMap = Galeri::Xpetra::CreateMap<LO, GO, Node>(TestHelpers::Parameters::getLib(), "Cartesian3D", comm, galeriList);
    RCP<const Map> dofMap  = Xpetra::MapFactory<LO, GO, Node>::Build(nodeMap, numDofsPerNode);
    RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector>> Pr =
        Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(matrixType, dofMap, galeriList);
    A                          = Pr->BuildMatrix();
    RCP<MultiVector> nullspace = Pr->BuildNullspace();

    TEST_ASSERT(!A.is_null());
    TEST_ASSERT(nullspace->getNumVectors() == 6);

    fineLevel->Set("A", A);
    fineLevel->Set("Nullspace", nullspace);
  }

  const std::string mapName = "Dummy Map";
  fineLevel->Set(mapName, A->getRowMap());

  TEST_ASSERT(fineLevel->IsAvailable("A", MueLu::NoFactory::get()));
  TEST_ASSERT(fineLevel->IsAvailable(mapName, MueLu::NoFactory::get()));

  RCP<FactoryManager> factoryManager = rcp(new FactoryManager());
  factoryManager->SetKokkosRefactor(false);
  factoryManager->SetFactory(mapName, MueLu::NoFactory::getRCP());
  fineLevel->SetFactoryManager(factoryManager);
  coarseLevel->SetFactoryManager(factoryManager);

  RCP<TentativePFactory> tentativePFact = rcp(new TentativePFactory());
  tentativePFact->SetParameter("Nullspace name", Teuchos::ParameterEntry(std::string("Nullspace")));

  RCP<MapTransferFactory> mapTransferFactory = rcp(new MapTransferFactory());
  mapTransferFactory->SetParameter("map: factory", Teuchos::ParameterEntry(mapName));
  mapTransferFactory->SetParameter("map: name", Teuchos::ParameterEntry(mapName));
  mapTransferFactory->SetFactory("P", tentativePFact);
  mapTransferFactory->SetParameter("nullspace vectors: limit to", Teuchos::ParameterEntry(std::string("translations")));

  coarseLevel->Request(mapName, MueLu::NoFactory::get());
  coarseLevel->Request("P", tentativePFact.get(), mapTransferFactory.get());
  coarseLevel->Request(*mapTransferFactory);  // This calls DeclareInput() on mapTransferFactory

  TEST_ASSERT(coarseLevel->IsRequested(mapName, MueLu::NoFactory::get()));
  TEST_ASSERT(coarseLevel->IsRequested("P", tentativePFact.get()));

  RCP<Matrix> Ptent = coarseLevel->Get<RCP<Matrix>>("P", tentativePFact.get());
  TEST_ASSERT(!Ptent.is_null());

  TEST_ASSERT(!coarseLevel->IsAvailable(mapName, MueLu::NoFactory::get()));
  mapTransferFactory->Build(*fineLevel, *coarseLevel);
  TEST_ASSERT(coarseLevel->IsAvailable(mapName, MueLu::NoFactory::get()));

  RCP<const Map> coarsenedMap = coarseLevel->Get<RCP<const Map>>(mapName, MueLu::NoFactory::get());
  TEST_ASSERT(!coarsenedMap.is_null());

  const global_size_t gNumTranslationalDOFs = Ptent->getDomainMap()->getGlobalNumElements() / 2;
  TEST_EQUALITY_CONST(coarsenedMap->getMinGlobalIndex(), Ptent->getDomainMap()->getMinGlobalIndex());
  TEST_EQUALITY_CONST(coarsenedMap->getGlobalNumElements(), gNumTranslationalDOFs);

  /* Manually extract translational part from domain map of prolongator and compare with coarsenedMap
   * Thereby, we can exploit the fact, that each nodes carries 6 unknowns on the coarse level: 3 translations
   * and 3 rotations. Here, we're only interested in the translations, i.e. the first 3 unknowns per node.
   */
  {
    ArrayView<const GO> domainMapGIDs = Ptent->getDomainMap()->getLocalElementList();
    Array<GO> gidsOfTranslations;
    const GO lNumCoarseNodes             = Teuchos::as<GO>(domainMapGIDs.size() / 6);
    const GO gNumCoarseTranslationalDOFs = Teuchos::as<GO>(Ptent->getDomainMap()->getGlobalNumElements() / 2);

    for (size_t nodeID = 0; nodeID < Teuchos::as<size_t>(lNumCoarseNodes); ++nodeID) {
      size_t idOfFirstDofOfThisNode = nodeID * 6;
      gidsOfTranslations.push_back(domainMapGIDs[idOfFirstDofOfThisNode]);
      gidsOfTranslations.push_back(domainMapGIDs[idOfFirstDofOfThisNode + 1]);
      gidsOfTranslations.push_back(domainMapGIDs[idOfFirstDofOfThisNode + 2]);
    }
    RCP<const Map> mapForComparison = MapFactory::Build(coarsenedMap->lib(),
                                                        gNumCoarseTranslationalDOFs, gidsOfTranslations, Teuchos::OrdinalTraits<GO>::zero(), comm);

    TEST_ASSERT(coarsenedMap->isSameAs(*mapForComparison));
  }
}

/* This tests coarsens the row map of the fine level operator, so the result from the MapTransferFactory
 * needs to match the domain map of the prolongator. This particular test targets the special case of
 * limiting the number of columns from P to be used for the transfer (e.g. when excluding the rotations
 * from the nullspace when creating the MapTransfer.)
 *
 * Assume a 2D elasticity discretization.
 */
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MapTransferFactory, TransferFullMap2DElasticityReducedP, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using test_factory = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;
  out << "Test transfer of a map with the MapTransferFactory" << std::endl;

  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();

  // Manual setup of a two-level hierarchy
  RCP<Level> fineLevel   = rcp(new Level());
  RCP<Level> coarseLevel = rcp(new Level());
  coarseLevel->SetPreviousLevel(fineLevel);
  fineLevel->SetLevelID(0);
  coarseLevel->SetLevelID(1);

  TEST_EQUALITY_CONST(fineLevel->GetLevelID(), 0);
  TEST_EQUALITY_CONST(coarseLevel->GetLevelID(), 1);

  // Create a 3D elasticity matrix and its nullspace needed to build a prolongator
  RCP<Matrix> A = Teuchos::null;
  {
    const GO nx                  = 10;
    const GO ny                  = 5;
    const int numDofsPerNode     = 2;  // 2D elasticity
    const std::string matrixType = "Elasticity2D";

    Teuchos::ParameterList galeriList;
    galeriList.set("nx", nx);
    galeriList.set("ny", ny);
    galeriList.set("matrixType", matrixType);
    RCP<const Map> nodeMap = Galeri::Xpetra::CreateMap<LO, GO, Node>(TestHelpers::Parameters::getLib(), "Cartesian2D", comm, galeriList);
    RCP<const Map> dofMap  = Xpetra::MapFactory<LO, GO, Node>::Build(nodeMap, numDofsPerNode);
    RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector>> Pr =
        Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(matrixType, dofMap, galeriList);
    A                          = Pr->BuildMatrix();
    RCP<MultiVector> nullspace = Pr->BuildNullspace();

    TEST_ASSERT(!A.is_null());
    TEST_ASSERT(nullspace->getNumVectors() == 3);

    fineLevel->Set("A", A);
    fineLevel->Set("Nullspace", nullspace);
  }

  const std::string mapName = "Dummy Map";
  fineLevel->Set(mapName, A->getRowMap());

  TEST_ASSERT(fineLevel->IsAvailable("A", MueLu::NoFactory::get()));
  TEST_ASSERT(fineLevel->IsAvailable(mapName, MueLu::NoFactory::get()));

  RCP<FactoryManager> factoryManager = rcp(new FactoryManager());
  factoryManager->SetKokkosRefactor(false);
  factoryManager->SetFactory(mapName, MueLu::NoFactory::getRCP());
  fineLevel->SetFactoryManager(factoryManager);
  coarseLevel->SetFactoryManager(factoryManager);

  RCP<TentativePFactory> tentativePFact = rcp(new TentativePFactory());
  tentativePFact->SetParameter("Nullspace name", Teuchos::ParameterEntry(std::string("Nullspace")));

  RCP<MapTransferFactory> mapTransferFactory = rcp(new MapTransferFactory());
  mapTransferFactory->SetParameter("map: factory", Teuchos::ParameterEntry(mapName));
  mapTransferFactory->SetParameter("map: name", Teuchos::ParameterEntry(mapName));
  mapTransferFactory->SetFactory("P", tentativePFact);
  mapTransferFactory->SetParameter("nullspace vectors: limit to", Teuchos::ParameterEntry(std::string("translations")));

  coarseLevel->Request(mapName, MueLu::NoFactory::get());
  coarseLevel->Request("P", tentativePFact.get(), mapTransferFactory.get());
  coarseLevel->Request(*mapTransferFactory);  // This calls DeclareInput() on mapTransferFactory

  TEST_ASSERT(coarseLevel->IsRequested(mapName, MueLu::NoFactory::get()));
  TEST_ASSERT(coarseLevel->IsRequested("P", tentativePFact.get()));

  RCP<Matrix> Ptent = coarseLevel->Get<RCP<Matrix>>("P", tentativePFact.get());
  TEST_ASSERT(!Ptent.is_null());

  TEST_ASSERT(!coarseLevel->IsAvailable(mapName, MueLu::NoFactory::get()));
  mapTransferFactory->Build(*fineLevel, *coarseLevel);
  TEST_ASSERT(coarseLevel->IsAvailable(mapName, MueLu::NoFactory::get()));

  RCP<const Map> coarsenedMap = coarseLevel->Get<RCP<const Map>>(mapName, MueLu::NoFactory::get());
  TEST_ASSERT(!coarsenedMap.is_null());

  const global_size_t gNumTranslationalDOFs = Ptent->getDomainMap()->getGlobalNumElements() * 2 / 3;
  TEST_EQUALITY_CONST(coarsenedMap->getMinGlobalIndex(), Ptent->getDomainMap()->getMinGlobalIndex());
  TEST_EQUALITY_CONST(coarsenedMap->getGlobalNumElements(), gNumTranslationalDOFs);

  /* Manually extract translational part from domain map of prolongator and compare with coarsenedMap
   * Thereby, we can exploit the fact, that each nodes carries 6 unknowns on the coarse level: 2 translations
   * and 1 rotation. Here, we're only interested in the translations, i.e. the first 2 unknowns per node.
   */
  {
    ArrayView<const GO> domainMapGIDs = Ptent->getDomainMap()->getLocalElementList();
    Array<GO> gidsOfTranslations;
    const GO lNumCoarseNodes             = Teuchos::as<GO>(domainMapGIDs.size() / 3);
    const GO gNumCoarseTranslationalDOFs = Teuchos::as<GO>(Ptent->getDomainMap()->getGlobalNumElements() * 2 / 3);

    for (size_t nodeID = 0; nodeID < Teuchos::as<size_t>(lNumCoarseNodes); ++nodeID) {
      size_t idOfFirstDofOfThisNode = nodeID * 3;
      gidsOfTranslations.push_back(domainMapGIDs[idOfFirstDofOfThisNode]);
      gidsOfTranslations.push_back(domainMapGIDs[idOfFirstDofOfThisNode + 1]);
    }
    RCP<const Map> mapForComparison = MapFactory::Build(coarsenedMap->lib(),
                                                        gNumCoarseTranslationalDOFs, gidsOfTranslations, Teuchos::OrdinalTraits<GO>::zero(), comm);

    TEST_ASSERT(coarsenedMap->isSameAs(*mapForComparison));
  }
}

/* This tests coarsens a subset of the row map of the fine level operator,
 * so the result from the MapTransferFactory needs to match a subset of the domain map of the prolongator.
 *
 * Assume a 1D Poisson discretization.
 */
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MapTransferFactory, TransferPartialMap1D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using MagnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  using test_factory  = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;
  out << "Test transfer of a map with the MapTransferFactory" << std::endl;

  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();

  // Manual setup of a two-level hierarchy
  RCP<Level> fineLevel   = rcp(new Level());
  RCP<Level> coarseLevel = rcp(new Level());
  coarseLevel->SetPreviousLevel(fineLevel);
  fineLevel->SetLevelID(0);
  coarseLevel->SetLevelID(1);

  TEST_EQUALITY_CONST(fineLevel->GetLevelID(), 0);
  TEST_EQUALITY_CONST(coarseLevel->GetLevelID(), 1);

  // Create a dummy matrix needed to build a prolongator
  const GO nx   = 49;
  RCP<Matrix> A = test_factory::Build1DPoisson(4 * nx);
  fineLevel->Set("A", A);

  // Extract a subset of A->getRowMap()'s GIDs to create the map to be transferred
  ArrayView<const GO> allRowGIDs = A->getRowMap()->getLocalElementList();
  Array<GO> myMapGIDs;
  for (LO lid = 0; lid < Teuchos::as<LO>(nx); ++lid) {
    if (lid % 3 == 0)
      myMapGIDs.push_back(allRowGIDs[lid]);
  }
  GO gNumFineEntries = 0;
  reduceAll(*comm, Teuchos::REDUCE_SUM, Teuchos::as<GO>(myMapGIDs.size()), Teuchos::outArg(gNumFineEntries));
  RCP<const Map> mapWithHoles = MapFactory::Build(TestHelpers::Parameters::getLib(), gNumFineEntries, myMapGIDs(), Teuchos::ScalarTraits<GO>::zero(), comm);

  const std::string mapName = "Dummy Map";
  fineLevel->Set(mapName, mapWithHoles);

  TEST_ASSERT(fineLevel->IsAvailable("A", MueLu::NoFactory::get()));
  TEST_ASSERT(fineLevel->IsAvailable(mapName, MueLu::NoFactory::get()));

  RCP<FactoryManager> factoryManager = rcp(new FactoryManager());
  factoryManager->SetKokkosRefactor(false);
  factoryManager->SetFactory(mapName, MueLu::NoFactory::getRCP());
  fineLevel->SetFactoryManager(factoryManager);
  coarseLevel->SetFactoryManager(factoryManager);

  RCP<TentativePFactory> tentativePFact = rcp(new TentativePFactory());

  RCP<MapTransferFactory> mapTransferFactory = rcp(new MapTransferFactory());
  mapTransferFactory->SetParameter("map: factory", Teuchos::ParameterEntry(mapName));
  mapTransferFactory->SetParameter("map: name", Teuchos::ParameterEntry(mapName));
  mapTransferFactory->SetFactory("P", tentativePFact);

  coarseLevel->Request(mapName, MueLu::NoFactory::get());
  coarseLevel->Request("P", tentativePFact.get(), mapTransferFactory.get());
  coarseLevel->Request(*mapTransferFactory);  // This calls DeclareInput() on mapTransferFactory

  TEST_ASSERT(coarseLevel->IsRequested(mapName, MueLu::NoFactory::get()));
  TEST_ASSERT(coarseLevel->IsRequested("P", tentativePFact.get()));

  RCP<Matrix> Ptent = coarseLevel->Get<RCP<Matrix>>("P", tentativePFact.get());
  TEST_ASSERT(!Ptent.is_null());

  TEST_ASSERT(!coarseLevel->IsAvailable(mapName, MueLu::NoFactory::get()));
  mapTransferFactory->Build(*fineLevel, *coarseLevel);
  TEST_ASSERT(coarseLevel->IsAvailable(mapName, MueLu::NoFactory::get()));

  RCP<const Map> coarsenedMap = coarseLevel->Get<RCP<const Map>>(mapName, MueLu::NoFactory::get());
  TEST_ASSERT(!coarsenedMap.is_null());

  /* Manually construct a coarse version of the mapWithHoles by
   * 1. populating a vector based on the mapWithHoles
   * 2. manually restricting that vector to the coarse level
   * 3. reconstructing the mapForComparison by only taking into accounts vector entries != 0
   */
  {
    // Populate a fine level vector w/ ones according to fineMap
    RCP<Vector> fullFineVec    = VectorFactory::Build(Ptent->getRangeMap(), true);
    RCP<Vector> partialFineVec = VectorFactory::Build(mapWithHoles, true);
    partialFineVec->putScalar(Teuchos::ScalarTraits<Scalar>::one());
    RCP<Import> fineImporter = ImportFactory::Build(mapWithHoles, fullFineVec->getMap());
    fullFineVec->doImport(*partialFineVec, *fineImporter, Xpetra::INSERT);

    // Restrict to coarse level manually
    RCP<Vector> fullCoarseVec = VectorFactory::Build(Ptent->getDomainMap(), true);
    Ptent->apply(*fullFineVec, *fullCoarseVec, Teuchos::TRANS);

    // Reconstruct coarse map for result checking
    ArrayRCP<const Scalar> coarseVecEntries = fullCoarseVec->getData(0);
    Array<GO> myCoarseGIDs;
    for (LO lid = 0; lid < static_cast<LO>(coarseVecEntries.size()); ++lid) {
      if (Teuchos::ScalarTraits<Scalar>::magnitude(coarseVecEntries[lid]) > Teuchos::ScalarTraits<MagnitudeType>::zero())
        myCoarseGIDs.push_back(Ptent->getDomainMap()->getGlobalElement(lid));
    }
    GO gNumCoarseEntries = 0;
    reduceAll(*comm, Teuchos::REDUCE_SUM, Teuchos::as<GO>(myCoarseGIDs.size()), Teuchos::outArg(gNumCoarseEntries));
    RCP<const Map> mapForComparison = MapFactory::Build(mapWithHoles->lib(), gNumCoarseEntries, myCoarseGIDs,
                                                        Teuchos::ScalarTraits<GO>::zero(), comm);

    TEST_ASSERT(coarsenedMap->isSameAs(*mapForComparison));
  }
}

/* This tests coarsens a subset of the row map of the fine level operator,
 * so the result from the MapTransferFactory needs to match a subset of the domain map of the prolongator.
 *
 * Assume a 2D Poisson discretization.
 */
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MapTransferFactory, TransferPartialMap2D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using MagnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  using test_factory  = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;
  out << "Test transfer of a map with the MapTransferFactory" << std::endl;

  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();

  // Manual setup of a two-level hierarchy
  RCP<Level> fineLevel   = rcp(new Level());
  RCP<Level> coarseLevel = rcp(new Level());
  coarseLevel->SetPreviousLevel(fineLevel);
  fineLevel->SetLevelID(0);
  coarseLevel->SetLevelID(1);

  TEST_EQUALITY_CONST(fineLevel->GetLevelID(), 0);
  TEST_EQUALITY_CONST(coarseLevel->GetLevelID(), 1);

  // Create a dummy matrix needed to build a prolongator
  const GO nx   = 49;
  RCP<Matrix> A = test_factory::Build2DPoisson(4 * nx);
  fineLevel->Set("A", A);

  // Extract a subset of A->getRowMap()'s GIDs to create the map to be transferred
  ArrayView<const GO> allRowGIDs = A->getRowMap()->getLocalElementList();
  Array<GO> myMapGIDs;
  for (LO lid = 0; lid < Teuchos::as<LO>(nx); ++lid) {
    if (lid % 3 == 0) {
      myMapGIDs.push_back(allRowGIDs[lid]);
      myMapGIDs.push_back(allRowGIDs[lid + 1]);
    }
  }
  GO gNumFineEntries = 0;
  reduceAll(*comm, Teuchos::REDUCE_SUM, Teuchos::as<GO>(myMapGIDs.size()), Teuchos::outArg(gNumFineEntries));
  RCP<const Map> mapWithHoles = MapFactory::Build(TestHelpers::Parameters::getLib(), gNumFineEntries, myMapGIDs(), Teuchos::ScalarTraits<GO>::zero(), comm);

  const std::string mapName = "Dummy Map";
  fineLevel->Set(mapName, mapWithHoles);

  TEST_ASSERT(fineLevel->IsAvailable("A", MueLu::NoFactory::get()));
  TEST_ASSERT(fineLevel->IsAvailable(mapName, MueLu::NoFactory::get()));

  RCP<FactoryManager> factoryManager = rcp(new FactoryManager());
  factoryManager->SetKokkosRefactor(false);
  factoryManager->SetFactory(mapName, MueLu::NoFactory::getRCP());
  fineLevel->SetFactoryManager(factoryManager);
  coarseLevel->SetFactoryManager(factoryManager);

  RCP<TentativePFactory> tentativePFact = rcp(new TentativePFactory());

  RCP<MapTransferFactory> mapTransferFactory = rcp(new MapTransferFactory());
  mapTransferFactory->SetParameter("map: factory", Teuchos::ParameterEntry(mapName));
  mapTransferFactory->SetParameter("map: name", Teuchos::ParameterEntry(mapName));
  mapTransferFactory->SetFactory("P", tentativePFact);

  coarseLevel->Request(mapName, MueLu::NoFactory::get());
  coarseLevel->Request("P", tentativePFact.get(), mapTransferFactory.get());
  coarseLevel->Request(*mapTransferFactory);  // This calls DeclareInput() on mapTransferFactory

  TEST_ASSERT(coarseLevel->IsRequested(mapName, MueLu::NoFactory::get()));
  TEST_ASSERT(coarseLevel->IsRequested("P", tentativePFact.get()));

  RCP<Matrix> Ptent = coarseLevel->Get<RCP<Matrix>>("P", tentativePFact.get());
  TEST_ASSERT(!Ptent.is_null());

  TEST_ASSERT(!coarseLevel->IsAvailable(mapName, MueLu::NoFactory::get()));
  mapTransferFactory->Build(*fineLevel, *coarseLevel);
  TEST_ASSERT(coarseLevel->IsAvailable(mapName, MueLu::NoFactory::get()));

  RCP<const Map> coarsenedMap = coarseLevel->Get<RCP<const Map>>(mapName, MueLu::NoFactory::get());
  TEST_ASSERT(!coarsenedMap.is_null());

  /* Manually construct a coarse version of the mapWithHoles by
   * 1. populating a vector based on the mapWithHoles
   * 2. manually restricting that vector to the coarse level
   * 3. reconstructing the mapForComparison by only taking into accounts vector entries != 0
   */
  {
    // Populate a fine level vector w/ ones according to fineMap
    RCP<Vector> fullFineVec    = VectorFactory::Build(Ptent->getRangeMap(), true);
    RCP<Vector> partialFineVec = VectorFactory::Build(mapWithHoles, true);
    partialFineVec->putScalar(Teuchos::ScalarTraits<Scalar>::one());
    RCP<Import> fineImporter = ImportFactory::Build(mapWithHoles, fullFineVec->getMap());
    fullFineVec->doImport(*partialFineVec, *fineImporter, Xpetra::INSERT);

    // Restrict to coarse level manually
    RCP<Vector> fullCoarseVec = VectorFactory::Build(Ptent->getDomainMap(), true);
    Ptent->apply(*fullFineVec, *fullCoarseVec, Teuchos::TRANS);

    // Reconstruct coarse map for result checking
    ArrayRCP<const Scalar> coarseVecEntries = fullCoarseVec->getData(0);
    Array<GO> myCoarseGIDs;
    for (LO lid = 0; lid < static_cast<LO>(coarseVecEntries.size()); ++lid) {
      if (Teuchos::ScalarTraits<Scalar>::magnitude(coarseVecEntries[lid]) > Teuchos::ScalarTraits<MagnitudeType>::zero())
        myCoarseGIDs.push_back(Ptent->getDomainMap()->getGlobalElement(lid));
    }
    GO gNumCoarseEntries = 0;
    reduceAll(*comm, Teuchos::REDUCE_SUM, Teuchos::as<GO>(myCoarseGIDs.size()), Teuchos::outArg(gNumCoarseEntries));
    RCP<const Map> mapForComparison = MapFactory::Build(mapWithHoles->lib(), gNumCoarseEntries, myCoarseGIDs,
                                                        Teuchos::ScalarTraits<GO>::zero(), comm);

    TEST_ASSERT(coarsenedMap->isSameAs(*mapForComparison));
  }
}

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node)                                                                         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MapTransferFactory, Constructor, Scalar, LO, GO, Node)                         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MapTransferFactory, TransferFullMap1D, Scalar, LO, GO, Node)                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MapTransferFactory, TransferFullMap2D, Scalar, LO, GO, Node)                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MapTransferFactory, TransferFullMap3DElasticity, Scalar, LO, GO, Node)         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MapTransferFactory, TransferFullMap3DElasticityReducedP, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MapTransferFactory, TransferFullMap2DElasticityReducedP, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MapTransferFactory, TransferPartialMap1D, Scalar, LO, GO, Node)                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MapTransferFactory, TransferPartialMap2D, Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests