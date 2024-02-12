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

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>

#include <Xpetra_IO.hpp>

#include "MueLu_Utilities.hpp"
#include "MueLu_SubBlockAFactory.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_RepartitionHeuristicFactory.hpp"
#include "MueLu_RepartitionFactory.hpp"
#include "MueLu_RepartitionBlockDiagonalFactory.hpp"
#include "MueLu_IsorropiaInterface.hpp"
#include "MueLu_RepartitionInterface.hpp"
#include "MueLu_CloneRepartitionInterface.hpp"
#include "MueLu_RebalanceBlockAcFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_BlockedRAPFactory.hpp"
#include "MueLu_BlockedPFactory.hpp"
#include "MueLu_RebalanceAcFactory.hpp"
#include "MueLu_RebalanceBlockAcFactory.hpp"
#include "MueLu_RebalanceBlockRestrictionFactory.hpp"
#include "MueLu_RebalanceBlockInterpolationFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_BlockedGaussSeidelSmoother.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_BlockedCoarseMapFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_BlockedCoordinatesTransferFactory.hpp"
#include "MueLu_CoordinatesTransferFactory.hpp"
#include "MueLu_RepartitionInterface.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"

#ifdef HAVE_MUELU_ZOLTAN
#include "MueLu_ZoltanInterface.hpp"
#endif
#ifdef HAVE_MUELU_ZOLTAN2
#include "MueLu_Zoltan2Interface.hpp"
#endif

#include "MueLu_Exceptions.hpp"

namespace MueLuTests {

/////////////////////////
// helper function
// note: we assume "domainmap" to be linear starting with GIDs from domainmap->getMinAllGlobalIndex() to
//       domainmap->getMaxAllGlobalIndex() and build a quadratic triangular matrix with the stencil (b,a,c)
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
GenerateProblemMatrix(const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > rangemap,
                      const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > domainmap,
                      Scalar a = 2.0, Scalar b = -1.0, Scalar c = -1.0) {
#include "MueLu_UseShortNames.hpp"
  Teuchos::RCP<CrsMatrixWrap> mtx = Galeri::Xpetra::MatrixTraits<Map, CrsMatrixWrap>::Build(rangemap, 3);

  LocalOrdinal NumMyRowElements = rangemap->getLocalNumElements();

  GlobalOrdinal minGColId       = domainmap->getMinAllGlobalIndex();  // minimum over all procs
  GlobalOrdinal maxGColId       = domainmap->getMaxAllGlobalIndex();  // maximum over all procs
  GlobalOrdinal numGColElements = domainmap->getGlobalNumElements();
  //    std::cout << maxGColId << " " << minGColId << " " << numGColElements <<std::endl;
  TEUCHOS_TEST_FOR_EXCEPTION(maxGColId - minGColId != numGColElements - 1, MueLu::Exceptions::RuntimeError, "GenerateProblemMatrix: incosistent number of map elements.");

  GlobalOrdinal minGRowId = rangemap->getMinAllGlobalIndex();  // minimum over all procs
  GlobalOrdinal maxGRowId = rangemap->getMaxAllGlobalIndex();  // maximum over all procs
  TEUCHOS_TEST_FOR_EXCEPTION(maxGRowId - minGRowId != maxGColId - minGColId, MueLu::Exceptions::RuntimeError, "GenerateProblemMatrix: incosistent number of map elements between range and domain maps.");

  GlobalOrdinal offset = minGColId - minGRowId;

  GlobalOrdinal NumEntries;
  LocalOrdinal nnz = 2;
  std::vector<Scalar> Values(nnz);
  std::vector<GlobalOrdinal> Indices(nnz);

  // loop over all local rows
  for (LocalOrdinal i = 0; i < NumMyRowElements; ++i) {
    GlobalOrdinal grid = rangemap->getGlobalElement(i);
    if (grid == minGRowId) {
      NumEntries = 1;
      Values[0]  = c;
      Indices[0] = minGColId + 1;
    } else if (grid == maxGRowId) {
      NumEntries = 1;
      Values[0]  = b;
      Indices[0] = maxGColId - 1;
    } else {
      NumEntries = 2;
      Indices[0] = offset + rangemap->getMinGlobalIndex() + i - 1;
      Indices[1] = offset + rangemap->getMinGlobalIndex() + i + 1;
      Values[0]  = b;
      Values[1]  = c;
    }
    // put the off-diagonal entries
    // Xpetra wants ArrayViews (sigh)
    Teuchos::ArrayView<Scalar> av(&Values[0], NumEntries);
    Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0], NumEntries);
    mtx->insertGlobalValues(rangemap->getGlobalElement(i), iv, av);

    // Put in the diagonal entry
    mtx->insertGlobalValues(grid,
                            Teuchos::tuple<GlobalOrdinal>(offset + rangemap->getMinGlobalIndex() + i),
                            Teuchos::tuple<Scalar>(a));
  }

  mtx->fillComplete(domainmap, rangemap);
  return mtx;

}  // GenerateProblemMatrix

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedRepartition, BlockedRAPFactory, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  if (comm->getSize() == 1) {
    out << "Skip BlockedRepartion test in serial case" << std::endl;
    return;
  }

  int nNumProcs = comm->getSize();  // number of procs used before rebalancing

  /**********************************************************************************/
  /* CREATE INITIAL MATRIX                                                          */
  /**********************************************************************************/
  RCP<const Map> bigMap;
  RCP<const Map> map1;
  RCP<const Map> map2;
  GO numElements  = 400;
  GO numElements1 = 200;
  GO numElements2 = 200;

  Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

  // the test matrix has to be a nxn block matrix with quadratic blocks
  // where the subblocks use consequent numbering of global DOF ids.
  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(1);

  map1 = StridedMapFactory::Build(lib, numElements1, 0, stridingInfo, comm);
  map2 = StridedMapFactory::Build(lib, numElements2, numElements1, stridingInfo, comm);

  std::vector<GlobalOrdinal> localGids;                                               // vector with all local GIDs on cur proc
  Teuchos::ArrayView<const GlobalOrdinal> map1eleList = map1->getLocalElementList();  // append all local gids from map1 and map2
  localGids.insert(localGids.end(), map1eleList.begin(), map1eleList.end());
  Teuchos::ArrayView<const GlobalOrdinal> map2eleList = map2->getLocalElementList();
  localGids.insert(localGids.end(), map2eleList.begin(), map2eleList.end());
  Teuchos::ArrayView<GlobalOrdinal> eleList(&localGids[0], localGids.size());
  bigMap = StridedMapFactory::Build(lib, numElements, eleList, 0, stridingInfo, comm);  // create full big map (concatenation of map1 and map2)
  std::vector<Teuchos::RCP<const Map> > maps;
  maps.push_back(map1);
  maps.push_back(map2);

  Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LO, GO, Node> > mapExtractor = Xpetra::MapExtractorFactory<Scalar, LO, GO, Node>::Build(bigMap, maps);

  RCP<CrsMatrixWrap> Op11 = GenerateProblemMatrix<Scalar, LO, GO, Node>(map1, map1, 2, -1, -1);
  RCP<CrsMatrixWrap> Op12 = GenerateProblemMatrix<Scalar, LO, GO, Node>(map1, map2, 1, 0, 0);
  RCP<CrsMatrixWrap> Op21 = GenerateProblemMatrix<Scalar, LO, GO, Node>(map2, map1, 1, 0, 0);
  RCP<CrsMatrixWrap> Op22 = GenerateProblemMatrix<Scalar, LO, GO, Node>(map2, map2, 3, -2, -1);

  // store output of simple MV products for OpIJ
  RCP<Vector> test11 = VectorFactory::Build(Op11->getDomainMap());
  test11->putScalar(1.0);
  RCP<Vector> test12 = VectorFactory::Build(Op12->getDomainMap());
  test12->putScalar(1.0);
  RCP<Vector> test21 = VectorFactory::Build(Op21->getDomainMap());
  test21->putScalar(1.0);
  RCP<Vector> test22 = VectorFactory::Build(Op22->getDomainMap());
  test22->putScalar(1.0);
  RCP<Vector> res11 = VectorFactory::Build(Op11->getRangeMap());
  RCP<Vector> res12 = VectorFactory::Build(Op12->getRangeMap());
  RCP<Vector> res21 = VectorFactory::Build(Op21->getRangeMap());
  RCP<Vector> res22 = VectorFactory::Build(Op22->getRangeMap());
  Op11->apply(*test11, *res11);
  Op12->apply(*test12, *res12);
  Op21->apply(*test21, *res21);
  Op22->apply(*test22, *res22);

  // build blocked operator
  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> > bOp = Teuchos::rcp(new Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>(mapExtractor, mapExtractor, 10));

  bOp->setMatrix(0, 0, Op11);
  bOp->setMatrix(0, 1, Op12);
  bOp->setMatrix(1, 0, Op21);
  bOp->setMatrix(1, 1, Op22);
  bOp->fillComplete();
  TEST_EQUALITY(bOp != Teuchos::null, true);
  TEST_EQUALITY(bOp->getGlobalNumEntries(), 2392);

  // build hierarchy
  RCP<Level> levelOne = rcp(new Level());
  levelOne->SetLevelID(-1);
  RCP<Level> levelTwo = rcp(new Level());
  levelTwo->SetPreviousLevel(levelOne);
  levelTwo->SetLevelID(0);
  levelTwo->Set("A", Teuchos::rcp_dynamic_cast<Matrix>(bOp));  // set blocked operator

  // define repartition heuristics
  RCP<RepartitionHeuristicFactory> RepHeuFact = Teuchos::rcp(new RepartitionHeuristicFactory);
  RepHeuFact->SetFactory("A", MueLu::NoFactory::getRCP());  // 2x2 blocked operator
  RepHeuFact->SetParameter("repartition: start level", Teuchos::ParameterEntry(0));
  RepHeuFact->SetParameter("repartition: min rows per proc", Teuchos::ParameterEntry(200));

  // define sub block factories for blocked operator "A"
  RCP<SubBlockAFactory> A11Fact = Teuchos::rcp(new SubBlockAFactory());
  A11Fact->SetFactory("A", MueLu::NoFactory::getRCP());
  A11Fact->SetParameter("block row", Teuchos::ParameterEntry(0));
  A11Fact->SetParameter("block col", Teuchos::ParameterEntry(0));
  RCP<SubBlockAFactory> A22Fact = Teuchos::rcp(new SubBlockAFactory());
  A22Fact->SetFactory("A", MueLu::NoFactory::getRCP());
  A22Fact->SetParameter("block row", Teuchos::ParameterEntry(1));
  A22Fact->SetParameter("block col", Teuchos::ParameterEntry(1));

  RCP<AmalgamationFactory> Amalg11Fact = Teuchos::rcp(new AmalgamationFactory());
  RCP<AmalgamationFactory> Amalg22Fact = Teuchos::rcp(new AmalgamationFactory());
  Amalg11Fact->SetFactory("A", A11Fact);
  Amalg22Fact->SetFactory("A", A22Fact);

#ifdef HAVE_MUELU_ISORROPIA
  RCP<Factory> Rep11Interface = Teuchos::null;
  if (TestHelpers::Parameters::getLib() == Xpetra::UseEpetra) {
    RCP<IsorropiaInterface> Iso11Interface = Teuchos::rcp(new IsorropiaInterface());
    Iso11Interface->SetFactory("A", A11Fact);
    Iso11Interface->SetFactory("number of partitions", RepHeuFact);
    Iso11Interface->SetFactory("UnAmalgamationInfo", Amalg11Fact);

    Rep11Interface = Teuchos::rcp(new RepartitionInterface());
    Rep11Interface->SetFactory("A", A11Fact);
    Rep11Interface->SetFactory("number of partitions", RepHeuFact);
    Rep11Interface->SetFactory("AmalgamatedPartition", Iso11Interface);
  } else {
    // we are in Tpetra mode (even though Isorropia would be available)
    // create dummy "Partition" array
    RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(Op11->getRowMap(), false);
    ArrayRCP<GO> decompEntries                         = decomposition->getDataNonConst(0);
    for (size_t r = 0; r < decomposition->getMap()->getLocalNumElements(); r++) {
      if (r % 2 == 0)
        decompEntries[r] = 0;
      else
        decompEntries[r] = 1;
    }
    Rep11Interface = Teuchos::rcp(new SubBlockAFactory());
    levelTwo->Request("Partition", Rep11Interface.get());
    levelTwo->Set("Partition", decomposition, Rep11Interface.get());
  }
#else
  // create dummy "Partition" array
  RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(Op11->getRowMap(), false);
  ArrayRCP<GO> decompEntries                         = decomposition->getDataNonConst(0);
  for (size_t r = 0; r < decomposition->getMap()->getLocalNumElements(); r++) {
    if (r % 2 == 0)
      decompEntries[r] = 0;
    else
      decompEntries[r] = 1;
  }
  RCP<SubBlockAFactory> Rep11Interface = Teuchos::rcp(new SubBlockAFactory());
  levelTwo->Request("Partition", Rep11Interface.get());
  levelTwo->Set("Partition", decomposition, Rep11Interface.get());
#endif

  RCP<RepartitionFactory> Rep11Factory = Teuchos::rcp(new RepartitionFactory);
  Rep11Factory->SetFactory("A", A11Fact);
  Rep11Factory->SetFactory("number of partitions", RepHeuFact);
  Rep11Factory->SetFactory("Partition", Rep11Interface);

  RCP<CloneRepartitionInterface> Rep22Interface = Teuchos::rcp(new CloneRepartitionInterface());
  Rep22Interface->SetFactory("A", A22Fact);
  Rep22Interface->SetFactory("Partition", Rep11Interface);

  RCP<RepartitionFactory> Rep22Factory = Teuchos::rcp(new RepartitionFactory);
  Rep22Factory->SetFactory("A", A22Fact);
  Rep22Factory->SetFactory("number of partitions", RepHeuFact);
  Rep22Factory->SetFactory("Partition", Rep22Interface);

  // set up factory manager
  RCP<FactoryManager> FC1 = rcp(new FactoryManager());
  FC1->SetFactory("A", A11Fact);
  FC1->SetFactory("Importer", Rep11Factory);
  FC1->SetIgnoreUserData(true);  // always use data from factories defined in factory manager

  RCP<FactoryManager> FC2 = rcp(new FactoryManager());
  FC2->SetFactory("A", A22Fact);
  FC2->SetFactory("Importer", Rep22Factory);
  FC2->SetIgnoreUserData(true);  // always use data from factories defined in factory manager

  /////////////////////////////////////////// define blocked transfer ops
  RCP<RebalanceBlockAcFactory> RebAFact = rcp(new RebalanceBlockAcFactory());
  RebAFact->SetParameter("repartition: use subcommunicators", Teuchos::ParameterEntry(false));
  RebAFact->SetFactory("A", MueLu::NoFactory::getRCP());
  RebAFact->AddFactoryManager(FC1);
  RebAFact->AddFactoryManager(FC2);

  /////////////////////////////////////////// request rebalanced coarse level matrix A
  levelTwo->Request("A", RebAFact.get(), MueLu::NoFactory::get());
  TEST_EQUALITY(levelTwo->IsRequested("A", RebAFact.get()), true);

  // request Partition data
  levelTwo->Request("Partition", Rep11Interface.get(), MueLu::NoFactory::get());
  levelTwo->Request("Partition", Rep22Interface.get(), MueLu::NoFactory::get());
  TEST_EQUALITY(levelTwo->IsRequested("Partition", Rep11Interface.get()), true);
  TEST_EQUALITY(levelTwo->IsRequested("Partition", Rep22Interface.get()), true);

  /////////////////////////////////////////// build rebalanced coarse level matrix A
  RCP<Matrix> rebA = levelTwo->Get<RCP<Matrix> >("A", RebAFact.get());
  TEST_EQUALITY(rebA != Teuchos::null, true);
  // get number of active processes used in rebalanced matrix
  std::vector<int> amActive  = std::vector<int>(comm->getSize(), 0);
  std::vector<int> areActive = std::vector<int>(comm->getSize(), 0);
  if (rebA->getLocalNumEntries() > 0) amActive[comm->getRank()] = 1;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, comm->getSize(), &amActive[0], &areActive[0]);
  int nNumProcsReb = 0;
  for (int p = 0; p < comm->getSize(); p++)
    nNumProcsReb += areActive[p];
  if (nNumProcsReb == nNumProcs) {
    out << "Skip detailed tests. Matrix was not rebalanced" << std::endl;
    return;
  }
  int expectedPartitions = 2;
  TEST_EQUALITY(nNumProcsReb, expectedPartitions);

  //////////////////////////////////////////////////
  // extract partitions
  RCP<Xpetra::Vector<GO, LO, GO, NO> > part1 = levelTwo->Get<RCP<Xpetra::Vector<GO, LO, GO, NO> > >("Partition", Rep11Interface.get());
  TEST_EQUALITY(part1 != Teuchos::null, true);
  RCP<Xpetra::Vector<GO, LO, GO, NO> > part2 = levelTwo->Get<RCP<Xpetra::Vector<GO, LO, GO, NO> > >("Partition", Rep22Interface.get());
  TEST_EQUALITY(part2 != Teuchos::null, true);
  TEST_EQUALITY(part1->getGlobalLength(), part2->getGlobalLength());
  TEST_EQUALITY(part1->getLocalLength(), part2->getLocalLength());

  Teuchos::ArrayRCP<const GO> part1_data = part1->getData(0);
  Teuchos::ArrayRCP<const GO> part2_data = part2->getData(0);
  for (size_t i = 0; i < part1->getLocalLength(); i++) {
    TEST_EQUALITY(part1->getMap()->getGlobalElement(i) + 200, part2->getMap()->getGlobalElement(i));
    TEST_EQUALITY(part1_data[i], part2_data[i]);
  }

  // check rebalanced operator
  RCP<Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> > bA = Teuchos::rcp_dynamic_cast<Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> >(rebA);
  TEST_EQUALITY(bA != Teuchos::null, true);
  TEST_EQUALITY(bA->Rows(), 2);
  TEST_EQUALITY(bA->Cols(), 2);

  RCP<Vector> rtest11 = VectorFactory::Build(bA->getMatrix(0, 0)->getDomainMap());
  rtest11->putScalar(1.0);
  RCP<Vector> rtest12 = VectorFactory::Build(bA->getMatrix(0, 1)->getDomainMap());
  rtest12->putScalar(1.0);
  RCP<Vector> rtest21 = VectorFactory::Build(bA->getMatrix(1, 0)->getDomainMap());
  rtest21->putScalar(1.0);
  RCP<Vector> rtest22 = VectorFactory::Build(bA->getMatrix(1, 1)->getDomainMap());
  rtest22->putScalar(1.0);
  RCP<Vector> rres11 = VectorFactory::Build(bA->getMatrix(0, 0)->getRangeMap());
  RCP<Vector> rres12 = VectorFactory::Build(bA->getMatrix(0, 1)->getRangeMap());
  RCP<Vector> rres21 = VectorFactory::Build(bA->getMatrix(1, 0)->getRangeMap());
  RCP<Vector> rres22 = VectorFactory::Build(bA->getMatrix(1, 1)->getRangeMap());
  bA->getMatrix(0, 0)->apply(*rtest11, *rres11);
  bA->getMatrix(0, 1)->apply(*rtest12, *rres12);
  bA->getMatrix(1, 0)->apply(*rtest21, *rres21);
  bA->getMatrix(1, 1)->apply(*rtest22, *rres22);
  TEST_EQUALITY(res11->norm1(), rres11->norm1());
  TEST_EQUALITY(res12->norm1(), rres12->norm1());
  TEST_EQUALITY(res21->norm1(), rres21->norm1());
  TEST_EQUALITY(res22->norm1(), rres22->norm1());
  TEST_EQUALITY(bA->getGlobalNumEntries(), 2392);
}  // BlockedRAPFactory

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedRepartition, BlockedRAPFactoryWithRestriction, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  if (comm->getSize() == 1) {
    out << "Skip BlockedRepartion test in serial case" << std::endl;
    return;
  }

  /**********************************************************************************/
  /* CREATE INITIAL MATRIX                                                          */
  /**********************************************************************************/
  RCP<const Map> bigMap;
  RCP<const Map> map1;
  RCP<const Map> map2;
  GO numElements  = 400;
  GO numElements1 = 200;
  GO numElements2 = 200;

  Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

  // the test matrix has to be a nxn block matrix with quadratic blocks
  // where the subblocks use consequent numbering of global DOF ids.
  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(1);

  map1 = StridedMapFactory::Build(lib, numElements1, 0, stridingInfo, comm);
  map2 = StridedMapFactory::Build(lib, numElements2, numElements1, stridingInfo, comm);

  std::vector<GlobalOrdinal> localGids;                                               // vector with all local GIDs on cur proc
  Teuchos::ArrayView<const GlobalOrdinal> map1eleList = map1->getLocalElementList();  // append all local gids from map1 and map2
  localGids.insert(localGids.end(), map1eleList.begin(), map1eleList.end());
  Teuchos::ArrayView<const GlobalOrdinal> map2eleList = map2->getLocalElementList();
  localGids.insert(localGids.end(), map2eleList.begin(), map2eleList.end());
  Teuchos::ArrayView<GlobalOrdinal> eleList(&localGids[0], localGids.size());
  bigMap = StridedMapFactory::Build(lib, numElements, eleList, 0, stridingInfo, comm);  // create full big map (concatenation of map1 and map2)
  std::vector<Teuchos::RCP<const Map> > maps;
  maps.push_back(map1);
  maps.push_back(map2);

  Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LO, GO, Node> > mapExtractor = Xpetra::MapExtractorFactory<Scalar, LO, GO, Node>::Build(bigMap, maps);

  RCP<CrsMatrixWrap> Op11 = GenerateProblemMatrix<Scalar, LO, GO, Node>(map1, map1, 2, -1, -1);
  RCP<CrsMatrixWrap> Op12 = GenerateProblemMatrix<Scalar, LO, GO, Node>(map1, map2, 1, 0, 0);
  RCP<CrsMatrixWrap> Op21 = GenerateProblemMatrix<Scalar, LO, GO, Node>(map2, map1, 1, 0, 0);
  RCP<CrsMatrixWrap> Op22 = GenerateProblemMatrix<Scalar, LO, GO, Node>(map2, map2, 3, -2, -1);

  // store output of simple MV products for OpIJ
  RCP<Vector> test11 = VectorFactory::Build(Op11->getDomainMap());
  test11->putScalar(1.0);
  RCP<Vector> test12 = VectorFactory::Build(Op12->getDomainMap());
  test12->putScalar(1.0);
  RCP<Vector> test21 = VectorFactory::Build(Op21->getDomainMap());
  test21->putScalar(1.0);
  RCP<Vector> test22 = VectorFactory::Build(Op22->getDomainMap());
  test22->putScalar(1.0);
  RCP<Vector> res11 = VectorFactory::Build(Op11->getRangeMap());
  RCP<Vector> res12 = VectorFactory::Build(Op12->getRangeMap());
  RCP<Vector> res21 = VectorFactory::Build(Op21->getRangeMap());
  RCP<Vector> res22 = VectorFactory::Build(Op22->getRangeMap());
  Op11->apply(*test11, *res11);
  Op12->apply(*test12, *res12);
  Op21->apply(*test21, *res21);
  Op22->apply(*test22, *res22);
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType res11norm = res11->norm1();
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType res12norm = res12->norm1();
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType res21norm = res21->norm1();
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType res22norm = res22->norm1();

  // build blocked operator
  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> > bOp = Teuchos::rcp(new Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>(mapExtractor, mapExtractor, 10));

  bOp->setMatrix(0, 0, Op11);
  bOp->setMatrix(0, 1, Op12);
  bOp->setMatrix(1, 0, Op21);
  bOp->setMatrix(1, 1, Op22);
  bOp->fillComplete();
  TEST_EQUALITY(bOp != Teuchos::null, true);
  TEST_EQUALITY(bOp->getGlobalNumEntries(), 2392);

  // build hierarchy
  RCP<Level> levelOne = rcp(new Level());
  levelOne->SetLevelID(-1);
  RCP<Level> levelTwo = rcp(new Level());
  levelTwo->SetPreviousLevel(levelOne);
  levelTwo->SetLevelID(0);
  levelTwo->Set("A", Teuchos::rcp_dynamic_cast<Matrix>(bOp));  // set blocked operator

  // define repartition heuristics
  RCP<RepartitionHeuristicFactory> RepHeuFact = Teuchos::rcp(new RepartitionHeuristicFactory);
  RepHeuFact->SetFactory("A", MueLu::NoFactory::getRCP());  // 2x2 blocked operator
  RepHeuFact->SetParameter("repartition: start level", Teuchos::ParameterEntry(0));
  RepHeuFact->SetParameter("repartition: min rows per proc", Teuchos::ParameterEntry(200));

  // define sub block factories for blocked operator "A"
  RCP<SubBlockAFactory> A11Fact = Teuchos::rcp(new SubBlockAFactory());
  A11Fact->SetFactory("A", MueLu::NoFactory::getRCP());
  A11Fact->SetParameter("block row", Teuchos::ParameterEntry(0));
  A11Fact->SetParameter("block col", Teuchos::ParameterEntry(0));
  RCP<SubBlockAFactory> A22Fact = Teuchos::rcp(new SubBlockAFactory());
  A22Fact->SetFactory("A", MueLu::NoFactory::getRCP());
  A22Fact->SetParameter("block row", Teuchos::ParameterEntry(1));
  A22Fact->SetParameter("block col", Teuchos::ParameterEntry(1));

  RCP<AmalgamationFactory> Amalg11Fact = Teuchos::rcp(new AmalgamationFactory());
  RCP<AmalgamationFactory> Amalg22Fact = Teuchos::rcp(new AmalgamationFactory());
  Amalg11Fact->SetFactory("A", A11Fact);
  Amalg22Fact->SetFactory("A", A22Fact);

#ifdef HAVE_MUELU_ISORROPIA
  RCP<Factory> Rep11Interface = Teuchos::null;
  if (TestHelpers::Parameters::getLib() == Xpetra::UseEpetra) {
    RCP<IsorropiaInterface> Iso11Interface = Teuchos::rcp(new IsorropiaInterface());
    Iso11Interface->SetFactory("A", A11Fact);
    Iso11Interface->SetFactory("number of partitions", RepHeuFact);
    Iso11Interface->SetFactory("UnAmalgamationInfo", Amalg11Fact);

    Rep11Interface = Teuchos::rcp(new RepartitionInterface());
    Rep11Interface->SetFactory("A", A11Fact);
    Rep11Interface->SetFactory("number of partitions", RepHeuFact);
    Rep11Interface->SetFactory("AmalgamatedPartition", Iso11Interface);
  } else {
    // we are in Tpetra mode (even though Isorropia would be available)
    // create dummy "Partition" array
    RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(Op11->getRowMap(), false);
    ArrayRCP<GO> decompEntries                         = decomposition->getDataNonConst(0);
    for (size_t r = 0; r < decomposition->getMap()->getLocalNumElements(); r++) {
      if (r % 2 == 0)
        decompEntries[r] = 0;
      else
        decompEntries[r] = 1;
    }
    Rep11Interface = Teuchos::rcp(new SubBlockAFactory());
    levelTwo->Request("Partition", Rep11Interface.get());
    levelTwo->Set("Partition", decomposition, Rep11Interface.get());
  }
#else
  // create dummy "Partition" array
  RCP<Xpetra::Vector<GO, LO, GO, NO> > decomposition = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(Op11->getRowMap(), false);
  ArrayRCP<GO> decompEntries                         = decomposition->getDataNonConst(0);
  for (size_t r = 0; r < decomposition->getMap()->getLocalNumElements(); r++) {
    if (r % 2 == 0)
      decompEntries[r] = 0;
    else
      decompEntries[r] = 1;
  }
  RCP<SubBlockAFactory> Rep11Interface = Teuchos::rcp(new SubBlockAFactory());
  levelTwo->Request("Partition", Rep11Interface.get());
  levelTwo->Set("Partition", decomposition, Rep11Interface.get());
#endif

  RCP<RepartitionFactory> Rep11Factory = Teuchos::rcp(new RepartitionFactory);
  Rep11Factory->SetFactory("A", A11Fact);
  Rep11Factory->SetFactory("number of partitions", RepHeuFact);
  Rep11Factory->SetFactory("Partition", Rep11Interface);

  RCP<CloneRepartitionInterface> Rep22Interface = Teuchos::rcp(new CloneRepartitionInterface());
  Rep22Interface->SetFactory("A", A22Fact);
  Rep22Interface->SetFactory("Partition", Rep11Interface);

  RCP<RepartitionFactory> Rep22Factory = Teuchos::rcp(new RepartitionFactory);
  Rep22Factory->SetFactory("A", A22Fact);
  Rep22Factory->SetFactory("number of partitions", RepHeuFact);
  Rep22Factory->SetFactory("Partition", Rep22Interface);

  // set up factory manager
  RCP<FactoryManager> FC1 = rcp(new FactoryManager());
  FC1->SetFactory("A", A11Fact);
  FC1->SetFactory("Importer", Rep11Factory);
  FC1->SetIgnoreUserData(true);  // always use data from factories defined in factory manager

  RCP<FactoryManager> FC2 = rcp(new FactoryManager());
  FC2->SetFactory("A", A22Fact);
  FC2->SetFactory("Importer", Rep22Factory);
  FC2->SetIgnoreUserData(true);  // always use data from factories defined in factory manager

  /////////////////////////////////////////// define blocked transfer ops
  RCP<RebalanceBlockAcFactory> RebAFact = rcp(new RebalanceBlockAcFactory());
  RebAFact->SetParameter("repartition: use subcommunicators", Teuchos::ParameterEntry(true));
  RebAFact->SetFactory("A", MueLu::NoFactory::getRCP());
  RebAFact->AddFactoryManager(FC1);
  RebAFact->AddFactoryManager(FC2);

  /////////////////////////////////////////// request rebalanced coarse level matrix A
  levelTwo->Request("A", RebAFact.get(), MueLu::NoFactory::get());
  TEST_EQUALITY(levelTwo->IsRequested("A", RebAFact.get()), true);

  // request Partition data
  levelTwo->Request("Partition", Rep11Interface.get(), MueLu::NoFactory::get());
  levelTwo->Request("Partition", Rep22Interface.get(), MueLu::NoFactory::get());
  TEST_EQUALITY(levelTwo->IsRequested("Partition", Rep11Interface.get()), true);
  TEST_EQUALITY(levelTwo->IsRequested("Partition", Rep22Interface.get()), true);

  /////////////////////////////////////////// build rebalanced coarse level matrix A
  RCP<Matrix> rebA = levelTwo->Get<RCP<Matrix> >("A", RebAFact.get());

  if (rebA != Teuchos::null) {
    TEST_EQUALITY(rebA->getRowMap()->getComm()->getSize(), 2);

    //////////////////////////////////////////////////
    // extract partitions
    RCP<Xpetra::Vector<GO, LO, GO, NO> > part1 = levelTwo->Get<RCP<Xpetra::Vector<GO, LO, GO, NO> > >("Partition", Rep11Interface.get());
    TEST_EQUALITY(part1 != Teuchos::null, true);
    RCP<Xpetra::Vector<GO, LO, GO, NO> > part2 = levelTwo->Get<RCP<Xpetra::Vector<GO, LO, GO, NO> > >("Partition", Rep22Interface.get());
    TEST_EQUALITY(part2 != Teuchos::null, true);
    TEST_EQUALITY(part1->getGlobalLength(), part2->getGlobalLength());
    TEST_EQUALITY(part1->getLocalLength(), part2->getLocalLength());

    Teuchos::ArrayRCP<const GO> part1_data = part1->getData(0);
    Teuchos::ArrayRCP<const GO> part2_data = part2->getData(0);
    for (size_t i = 0; i < part1->getLocalLength(); i++) {
      TEST_EQUALITY(part1->getMap()->getGlobalElement(i) + 200, part2->getMap()->getGlobalElement(i));
      TEST_EQUALITY(part1_data[i], part2_data[i]);
    }

    // check rebalanced operator
    RCP<Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> > bA = Teuchos::rcp_dynamic_cast<Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> >(rebA);
    TEST_EQUALITY(bA != Teuchos::null, true);
    TEST_EQUALITY(bA->Rows(), 2);
    TEST_EQUALITY(bA->Cols(), 2);

    RCP<Vector> rtest11 = VectorFactory::Build(bA->getMatrix(0, 0)->getDomainMap());
    rtest11->putScalar(1.0);
    RCP<Vector> rtest12 = VectorFactory::Build(bA->getMatrix(0, 1)->getDomainMap());
    rtest12->putScalar(1.0);
    RCP<Vector> rtest21 = VectorFactory::Build(bA->getMatrix(1, 0)->getDomainMap());
    rtest21->putScalar(1.0);
    RCP<Vector> rtest22 = VectorFactory::Build(bA->getMatrix(1, 1)->getDomainMap());
    rtest22->putScalar(1.0);
    RCP<Vector> rres11 = VectorFactory::Build(bA->getMatrix(0, 0)->getRangeMap());
    RCP<Vector> rres12 = VectorFactory::Build(bA->getMatrix(0, 1)->getRangeMap());
    RCP<Vector> rres21 = VectorFactory::Build(bA->getMatrix(1, 0)->getRangeMap());
    RCP<Vector> rres22 = VectorFactory::Build(bA->getMatrix(1, 1)->getRangeMap());
    bA->getMatrix(0, 0)->apply(*rtest11, *rres11);
    bA->getMatrix(0, 1)->apply(*rtest12, *rres12);
    bA->getMatrix(1, 0)->apply(*rtest21, *rres21);
    bA->getMatrix(1, 1)->apply(*rtest22, *rres22);
    TEST_EQUALITY(res11norm, rres11->norm1());
    TEST_EQUALITY(res12norm, rres12->norm1());
    TEST_EQUALITY(res21norm, rres21->norm1());
    TEST_EQUALITY(res22norm, rres22->norm1());
    TEST_EQUALITY(bA->getGlobalNumEntries(), 2392);
  }  // end if rebA != Teuchos::null
}  // BlockedRAPFactoryWithRestriction

#if defined(HAVE_MUELU_ISORROPIA) && defined(HAVE_MUELU_ZOLTAN)
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedRepartition, BlockedRAPFactoryWithDiagonal, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  if (comm->getSize() == 1) {
    out << "Skip BlockedRepartion test in serial case" << std::endl;
    return;
  }

  //    int NumProcs = comm->getSize(); // number of procs used before rebalancing

  /**********************************************************************************/
  /* CREATE INITIAL MATRIX                                                          */
  /**********************************************************************************/
  RCP<const Map> bigMap;
  RCP<const Map> map1;
  RCP<const Map> map2;
  GO numElements  = 400;
  GO numElements1 = 200;
  GO numElements2 = 200;

  Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

  // the test matrix has to be a nxn block matrix with quadratic blocks
  // where the subblocks use consequent numbering of global DOF ids.
  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(1);

  map1 = StridedMapFactory::Build(lib, numElements1, 0, stridingInfo, comm);
  map2 = StridedMapFactory::Build(lib, numElements2, 0, stridingInfo, comm, -1, numElements1);

  std::vector<GlobalOrdinal> localGids;                                               // vector with all local GIDs on cur proc
  Teuchos::ArrayView<const GlobalOrdinal> map1eleList = map1->getLocalElementList();  // append all local gids from map1 and map2
  localGids.insert(localGids.end(), map1eleList.begin(), map1eleList.end());
  Teuchos::ArrayView<const GlobalOrdinal> map2eleList = map2->getLocalElementList();
  localGids.insert(localGids.end(), map2eleList.begin(), map2eleList.end());
  Teuchos::ArrayView<GlobalOrdinal> eleList(&localGids[0], localGids.size());
  bigMap = StridedMapFactory::Build(lib, numElements, eleList, 0, stridingInfo, comm);  // create full big map (concatenation of map1 and map2)
  std::vector<Teuchos::RCP<const Map> > maps;
  maps.push_back(map1);
  maps.push_back(map2);

  Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LO, GO, Node> > mapExtractor = Xpetra::MapExtractorFactory<Scalar, LO, GO, Node>::Build(bigMap, maps);

  RCP<CrsMatrixWrap> Op11 = GenerateProblemMatrix<Scalar, LO, GO, Node>(map1, map1, 2, -1, -1);
  RCP<CrsMatrixWrap> Op12 = GenerateProblemMatrix<Scalar, LO, GO, Node>(map1, map2, 1, 0, 0);
  RCP<CrsMatrixWrap> Op21 = GenerateProblemMatrix<Scalar, LO, GO, Node>(map2, map1, 1, 0, 0);
  RCP<CrsMatrixWrap> Op22 = GenerateProblemMatrix<Scalar, LO, GO, Node>(map2, map2, 3, -2, -1);

  // store output of simple MV products for OpIJ
  RCP<Vector> test11 = VectorFactory::Build(Op11->getDomainMap());
  test11->putScalar(1.0);
  RCP<Vector> test12 = VectorFactory::Build(Op12->getDomainMap());
  test12->putScalar(1.0);
  RCP<Vector> test21 = VectorFactory::Build(Op21->getDomainMap());
  test21->putScalar(1.0);
  RCP<Vector> test22 = VectorFactory::Build(Op22->getDomainMap());
  test22->putScalar(1.0);
  RCP<Vector> res11 = VectorFactory::Build(Op11->getRangeMap());
  RCP<Vector> res12 = VectorFactory::Build(Op12->getRangeMap());
  RCP<Vector> res21 = VectorFactory::Build(Op21->getRangeMap());
  RCP<Vector> res22 = VectorFactory::Build(Op22->getRangeMap());
  Op11->apply(*test11, *res11);
  Op12->apply(*test12, *res12);
  Op21->apply(*test21, *res21);
  Op22->apply(*test22, *res22);

  // build blocked operator
  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> > bOp = Teuchos::rcp(new Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>(mapExtractor, mapExtractor, 10));

  bOp->setMatrix(0, 0, Op11);
  bOp->setMatrix(0, 1, Op12);
  bOp->setMatrix(1, 0, Op21);
  bOp->setMatrix(1, 1, Op22);
  bOp->fillComplete();
  TEST_EQUALITY(bOp != Teuchos::null, true);
  TEST_EQUALITY(bOp->getGlobalNumEntries(), 2392);

  // coordinates
  RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> > coord = Xpetra::MultiVectorFactory<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>::Build(bOp->getFullRangeMap(), 1);
  int PID                                                                                            = comm->getRank();
  // GH: scope data manipulation because we should not let a pointer to this data live once we call MueLu
  {
    Teuchos::ArrayRCP<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> coordData = coord->getDataNonConst(0);
    for (size_t i = 0; i < (size_t)coordData.size(); ++i)
      coordData[i] = PID + (typename Teuchos::ScalarTraits<Scalar>::magnitudeType)i / coordData.size();
  }
  // nullspace
  RCP<MultiVector> nullspace    = MultiVectorFactory::Build(bOp->getFullRangeMap(), 1);
  Teuchos::ArrayRCP<SC> nspData = nullspace->getDataNonConst(0);
  for (size_t i = 0; i < (size_t)nspData.size(); i++)
    nspData[i] = 1.0;

  // Grab sub-blocks for the Tobias-style goodies
  RCP<MultiVector> nullspace1 = mapExtractor->ExtractVector(nullspace, 0);
  RCP<MultiVector> nullspace2 = mapExtractor->ExtractVector(nullspace, 1);

  // build hierarchy
  Hierarchy H;
  H.SetMaxCoarseSize(10);  // FIXME: At present this does not work past two levels
  RCP<Level> levelOne = H.GetLevel();
  levelOne->Set("A", Teuchos::rcp_dynamic_cast<Matrix>(bOp));  // set blocked operator
  levelOne->Set("Coordinates", coord);
  levelOne->Set("Nullspace1", nullspace1);
  levelOne->Set("Nullspace2", nullspace2);

  // Define the operator's subblocks
  RCP<SubBlockAFactory> A11Fact = Teuchos::rcp(new SubBlockAFactory());
  A11Fact->SetFactory("A", MueLu::NoFactory::getRCP());
  A11Fact->SetParameter("block row", Teuchos::ParameterEntry(0));
  A11Fact->SetParameter("block col", Teuchos::ParameterEntry(0));
  RCP<SubBlockAFactory> A22Fact = Teuchos::rcp(new SubBlockAFactory());
  A22Fact->SetFactory("A", MueLu::NoFactory::getRCP());
  A22Fact->SetParameter("block row", Teuchos::ParameterEntry(1));
  A22Fact->SetParameter("block col", Teuchos::ParameterEntry(1));

  // Create smoothers for each block
  std::string ifpackType;
  Teuchos::ParameterList ifpackList;
  ifpackList.set("relaxation: sweeps", (LO)5);
  ifpackList.set("relaxation: damping factor", (SC)1.0);
  ifpackType = "RELAXATION";
  ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
  RCP<SmootherPrototype> smoProto11 = rcp(new TrilinosSmoother(ifpackType, ifpackList, 0));
  smoProto11->SetFactory("A", A11Fact);
  RCP<SmootherPrototype> smoProto22 = rcp(new TrilinosSmoother(ifpackType, ifpackList, 0));
  smoProto22->SetFactory("A", A22Fact);

  RCP<SmootherFactory> Smoo11Fact = rcp(new SmootherFactory(smoProto11));
  RCP<SmootherFactory> Smoo22Fact = rcp(new SmootherFactory(smoProto22));

  // Define how the transfer operators for block (0,0) is computed
  RCP<Factory> Amalg11  = rcp(new AmalgamationFactory());
  RCP<Factory> CDFact11 = rcp(new CoalesceDropFactory());
  CDFact11->SetFactory("UnAmalgamationInfo", Amalg11);
  RCP<Factory> AggFact11 = rcp(new UncoupledAggregationFactory());
  AggFact11->SetFactory("Graph", CDFact11);
  AggFact11->SetFactory("DofsPerNode", CDFact11);
  RCP<TentativePFactory> P11Fact = rcp(new TentativePFactory());
  P11Fact->SetFactory("Aggregates", AggFact11);
  P11Fact->SetFactory("UnAmalgamationInfo", Amalg11);
  P11Fact->SetParameter("tentative: build coarse coordinates", Teuchos::ParameterEntry(false));
  RCP<TransPFactory> R11Fact = rcp(new TransPFactory());
  RCP<Factory> Nullspace11   = rcp(new NullspaceFactory());
  Nullspace11->SetParameter("Fine level nullspace", Teuchos::ParameterEntry(std::string("Nullspace1")));
  Nullspace11->SetFactory("Nullspace1", P11Fact);
  RCP<Factory> Coord11 = rcp(new CoordinatesTransferFactory());
  RCP<Factory> Cmap11  = rcp(new CoarseMapFactory());
  Coord11->SetFactory("Aggregates", AggFact11);
  Coord11->SetFactory("CoarseMap", Cmap11);

  // Create a factory manager for block (0,0)
  RCP<FactoryManager> M11 = rcp(new FactoryManager());
  M11->SetFactory("A", A11Fact);
  M11->SetFactory("Graph", CDFact11);
  M11->SetFactory("DofsPerNode", CDFact11);
  M11->SetFactory("Aggregates", AggFact11);
  M11->SetFactory("UnAmalgamationInfo", Amalg11);
  M11->SetFactory("P", P11Fact);
  M11->SetFactory("R", R11Fact);
  M11->SetFactory("Ptent", P11Fact);  // for Nullspace
  M11->SetFactory("Coordinates", Coord11);
  M11->SetFactory("Smoother", Smoo11Fact);
  M11->SetFactory("CoarseMap", Cmap11);
  M11->SetFactory("Nullspace", Nullspace11);
  M11->SetIgnoreUserData(true);

  // Define how the transfer operators for block (1,1) is computed
  RCP<Factory> Amalg22  = rcp(new AmalgamationFactory());
  RCP<Factory> CDFact22 = rcp(new CoalesceDropFactory());
  CDFact22->SetFactory("UnAmalgamationInfo", Amalg22);
  RCP<Factory> AggFact22 = rcp(new UncoupledAggregationFactory());
  AggFact22->SetFactory("Graph", CDFact22);
  AggFact22->SetFactory("DofsPerNode", CDFact22);
  RCP<TentativePFactory> P22Fact = rcp(new TentativePFactory());
  P22Fact->SetFactory("Aggregates", AggFact22);
  P22Fact->SetFactory("UnAmalgamationInfo", Amalg22);
  P22Fact->SetParameter("tentative: build coarse coordinates", Teuchos::ParameterEntry(false));
  RCP<TransPFactory> R22Fact = rcp(new TransPFactory());
  RCP<Factory> Nullspace22   = rcp(new NullspaceFactory());
  Nullspace22->SetParameter("Fine level nullspace", Teuchos::ParameterEntry(std::string("Nullspace2")));
  Nullspace22->SetFactory("Nullspace2", P22Fact);
  RCP<Factory> Cmap22 = rcp(new BlockedCoarseMapFactory());
  Cmap22->SetFactory("Aggregates", AggFact22);
  Cmap22->SetFactory("CoarseMap", Cmap11);
  RCP<Factory> Coord22 = rcp(new CoordinatesTransferFactory());
  Coord22->SetFactory("Aggregates", AggFact22);
  Coord22->SetFactory("CoarseMap", Cmap22);

  // Create a factory manager for block (1,1)
  RCP<FactoryManager> M22 = rcp(new FactoryManager());
  M22->SetFactory("A", A22Fact);
  M22->SetFactory("Graph", CDFact22);
  M22->SetFactory("DofsPerNode", CDFact22);
  M22->SetFactory("Aggregates", AggFact22);
  M22->SetFactory("UnAmalgamationInfo", Amalg22);
  M22->SetFactory("P", P22Fact);
  M22->SetFactory("R", R22Fact);
  M22->SetFactory("Ptent", P22Fact);  // for Nullspace
  M22->SetFactory("Coordinates", P22Fact);
  M22->SetFactory("Smoother", Smoo22Fact);
  M22->SetFactory("CoarseMap", Cmap22);
  M22->SetFactory("Nullspace", Nullspace22);
  M22->SetIgnoreUserData(true);

  RCP<BlockedPFactory> PFact = rcp(new BlockedPFactory());
  // We might want to set a afactory for "A" here?
  PFact->AddFactoryManager(M11);
  PFact->AddFactoryManager(M22);

  RCP<GenericRFactory> RFact = rcp(new GenericRFactory());
  RFact->SetFactory("P", PFact);

  RCP<BlockedCoordinatesTransferFactory> CoordXfer = rcp(new BlockedCoordinatesTransferFactory());
  CoordXfer->AddFactory(Coord11);
  CoordXfer->AddFactory(Coord22);
  CoordXfer->SetFactory("CoarseMap", PFact);

  RCP<BlockedRAPFactory> AcFact = rcp(new BlockedRAPFactory());
  AcFact->AddTransferFactory(CoordXfer);
  AcFact->SetFactory("P", PFact);
  AcFact->SetFactory("R", RFact);

  // Smoothers
  RCP<BlockedGaussSeidelSmoother> smootherPrototype = rcp(new BlockedGaussSeidelSmoother());
  smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(2));
  smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<SC>(1.0)));
  smootherPrototype->AddFactoryManager(M11, 0);
  smootherPrototype->AddFactoryManager(M22, 1);
  RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

  // Coarse grid correction
  RCP<BlockedGaussSeidelSmoother> coarseSolverPrototype = rcp(new BlockedGaussSeidelSmoother());
  coarseSolverPrototype->AddFactoryManager(M11, 0);
  coarseSolverPrototype->AddFactoryManager(M22, 1);
  RCP<SmootherFactory> coarseSolverFact = rcp(new SmootherFactory(coarseSolverPrototype, Teuchos::null));

  // Repartitioning
  RCP<RepartitionHeuristicFactory> RepHeuFact = Teuchos::rcp(new RepartitionHeuristicFactory);
  RepHeuFact->SetFactory("A", AcFact);  // MueLu::NoFactory::getRCP());
  RepHeuFact->SetParameter("repartition: start level", Teuchos::ParameterEntry(0));
  RepHeuFact->SetParameter("repartition: min rows per proc", Teuchos::ParameterEntry(40));

  RCP<Factory> BDFact = rcp(new RepartitionBlockDiagonalFactory());
  BDFact->SetFactory("A", AcFact);

  RCP<ZoltanInterface> zoltan = rcp(new ZoltanInterface());
  zoltan->SetFactory("A", BDFact);
  zoltan->SetFactory("number of partitions", RepHeuFact);
  zoltan->SetFactory("Coordinates", CoordXfer);

  RCP<RepartitionFactory> repart = rcp(new RepartitionFactory());
  repart->SetFactory("A", BDFact);
  repart->SetFactory("number of partitions", RepHeuFact);
  repart->SetFactory("Partition", zoltan);

  RCP<RebalanceBlockAcFactory> RebAcFact = rcp(new RebalanceBlockAcFactory());
  RebAcFact->SetFactory("A", AcFact);
  RebAcFact->SetFactory("Importer", repart);
  RebAcFact->SetFactory("SubImporters", repart);

  RCP<RebalanceBlockRestrictionFactory> RebRFact = rcp(new RebalanceBlockRestrictionFactory());
  RebRFact->AddFactoryManager(M11);
  RebRFact->AddFactoryManager(M22);
  RebRFact->SetFactory("SubImporters", repart);
  RebRFact->SetFactory("R", RFact);
  RebRFact->SetSingleSourceImporters(true);

  RCP<RebalanceBlockInterpolationFactory> RebPFact = rcp(new RebalanceBlockInterpolationFactory());

  // main factory manager
  FactoryManager M;
  M.SetFactory("A", RebAcFact);
  M.SetFactory("P", PFact);
  M.SetFactory("R", RebRFact);
  M.SetFactory("Coordinates", CoordXfer);
  M.SetFactory("Smoother", smootherFact);  // TODO fix me
  M.SetFactory("CoarseSolver", coarseSolverFact);
  M.SetFactory("Nullspace", RebRFact);

  // setup hierarchy
  H.SetVerbLevel(MueLu::Test);
  H.EnableGraphDumping("dep_graph", 0);

  H.Setup(M);

}  // BlockedRAPFactoryWithDiagonal
#endif

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedRepartition, BlockedRAPFactoryCoords, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  if (comm->getSize() == 1) {
    out << "Skip BlockedRepartion test in serial case" << std::endl;
    return;
  }

  /**********************************************************************************/
  /* CREATE INITIAL MATRIX                                                          */
  /**********************************************************************************/
  RCP<const Map> bigMap;
  RCP<const Map> map1;
  RCP<const Map> map2;
  GO numElements  = 400;
  GO numElements1 = 200;
  GO numElements2 = 200;

  Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

  // the test matrix has to be a nxn block matrix with quadratic blocks
  // where the subblocks use consequent numbering of global DOF ids.
  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(2);

  map1 = StridedMapFactory::Build(lib, numElements1, 0, stridingInfo, comm);
  map2 = StridedMapFactory::Build(lib, numElements2, 0, stridingInfo, comm, -1, numElements1);

  std::vector<GlobalOrdinal> localGids;                                               // vector with all local GIDs on cur proc
  Teuchos::ArrayView<const GlobalOrdinal> map1eleList = map1->getLocalElementList();  // append all local gids from map1 and map2
  localGids.insert(localGids.end(), map1eleList.begin(), map1eleList.end());
  Teuchos::ArrayView<const GlobalOrdinal> map2eleList = map2->getLocalElementList();
  localGids.insert(localGids.end(), map2eleList.begin(), map2eleList.end());
  Teuchos::ArrayView<GlobalOrdinal> eleList(&localGids[0], localGids.size());
  bigMap = StridedMapFactory::Build(lib, numElements, eleList, 0, stridingInfo, comm);  // create full big map (concatenation of map1 and map2)
  std::vector<Teuchos::RCP<const Map> > maps;
  maps.push_back(map1);
  maps.push_back(map2);

  Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LO, GO, Node> > mapExtractor = Xpetra::MapExtractorFactory<Scalar, LO, GO, Node>::Build(bigMap, maps);

  RCP<CrsMatrixWrap> Op11 = GenerateProblemMatrix<Scalar, LO, GO, Node>(map1, map1, 2, -1, -1);
  RCP<CrsMatrixWrap> Op12 = GenerateProblemMatrix<Scalar, LO, GO, Node>(map1, map2, 1, 0, 0);
  RCP<CrsMatrixWrap> Op21 = GenerateProblemMatrix<Scalar, LO, GO, Node>(map2, map1, 1, 0, 0);
  RCP<CrsMatrixWrap> Op22 = GenerateProblemMatrix<Scalar, LO, GO, Node>(map2, map2, 3, -2, -1);

  // store output of simple MV products for OpIJ
  RCP<Vector> test11 = VectorFactory::Build(Op11->getDomainMap());
  test11->putScalar(1.0);
  RCP<Vector> test12 = VectorFactory::Build(Op12->getDomainMap());
  test12->putScalar(1.0);
  RCP<Vector> test21 = VectorFactory::Build(Op21->getDomainMap());
  test21->putScalar(1.0);
  RCP<Vector> test22 = VectorFactory::Build(Op22->getDomainMap());
  test22->putScalar(1.0);
  RCP<Vector> res11 = VectorFactory::Build(Op11->getRangeMap());
  RCP<Vector> res12 = VectorFactory::Build(Op12->getRangeMap());
  RCP<Vector> res21 = VectorFactory::Build(Op21->getRangeMap());
  RCP<Vector> res22 = VectorFactory::Build(Op22->getRangeMap());
  Op11->apply(*test11, *res11);
  Op12->apply(*test12, *res12);
  Op21->apply(*test21, *res21);
  Op22->apply(*test22, *res22);

  // build blocked operator
  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> > bOp = Teuchos::rcp(new Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>(mapExtractor, mapExtractor, 10));

  bOp->setMatrix(0, 0, Op11);
  bOp->setMatrix(0, 1, Op12);
  bOp->setMatrix(1, 0, Op21);
  bOp->setMatrix(1, 1, Op22);
  bOp->fillComplete();
  TEST_EQUALITY(bOp != Teuchos::null, true);
  TEST_EQUALITY(bOp->getGlobalNumEntries(), 2392);

  // coordinates
  RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> > coord = Xpetra::MultiVectorFactory<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>::Build(bOp->getFullRangeMap(), 1);
  int PID                                                                                            = comm->getRank();
  // GH: scope data manipulation because we should not let a pointer to this data live once we call MueLu
  {
    Teuchos::ArrayRCP<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> coordData = coord->getDataNonConst(0);
    for (size_t i = 0; i < (size_t)coordData.size(); ++i)
      coordData[i] = PID + (typename Teuchos::ScalarTraits<Scalar>::magnitudeType)i / coordData.size();
  }

  // nullspace
  RCP<MultiVector> nullspace    = MultiVectorFactory::Build(bOp->getFullRangeMap(), 1);
  Teuchos::ArrayRCP<SC> nspData = nullspace->getDataNonConst(0);
  for (size_t i = 0; i < (size_t)nspData.size(); i++)
    nspData[i] = 1.0;

  // Grab sub-blocks for the Tobias-style goodies
  RCP<MultiVector> nullspace1 = mapExtractor->ExtractVector(nullspace, 0);
  RCP<MultiVector> nullspace2 = mapExtractor->ExtractVector(nullspace, 1);

  // build hierarchy
  Hierarchy H;
  H.SetMaxCoarseSize(10);  // FIXME: At present this does not work past two levels
  RCP<Level> levelOne = H.GetLevel();
  levelOne->Set("A", Teuchos::rcp_dynamic_cast<Matrix>(bOp));  // set blocked operator
  levelOne->Set("Coordinates", coord);
  levelOne->Set("Nullspace1", nullspace1);
  levelOne->Set("Nullspace2", nullspace2);

  // Define the operator's subblocks
  RCP<SubBlockAFactory> A11Fact = Teuchos::rcp(new SubBlockAFactory());
  A11Fact->SetParameter("block row", Teuchos::ParameterEntry(0));
  A11Fact->SetParameter("block col", Teuchos::ParameterEntry(0));
  RCP<SubBlockAFactory> A22Fact = Teuchos::rcp(new SubBlockAFactory());
  A22Fact->SetParameter("block row", Teuchos::ParameterEntry(1));
  A22Fact->SetParameter("block col", Teuchos::ParameterEntry(1));

  // Create smoothers for each block
  std::string ifpackType;
  Teuchos::ParameterList ifpackList;
  ifpackList.set("relaxation: sweeps", (LO)5);
  ifpackList.set("relaxation: damping factor", (SC)1.0);
  ifpackType = "RELAXATION";
  ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
  RCP<SmootherPrototype> smoProto11 = rcp(new TrilinosSmoother(ifpackType, ifpackList, 0));
  smoProto11->SetFactory("A", A11Fact);
  RCP<SmootherPrototype> smoProto22 = rcp(new TrilinosSmoother(ifpackType, ifpackList, 0));
  smoProto22->SetFactory("A", A22Fact);

  RCP<SmootherFactory> Smoo11Fact = rcp(new SmootherFactory(smoProto11));
  RCP<SmootherFactory> Smoo22Fact = rcp(new SmootherFactory(smoProto22));

  // Define how the transfer operators for block (0,0) is computed
  RCP<Factory> Amalg11  = rcp(new AmalgamationFactory());
  RCP<Factory> CDFact11 = rcp(new CoalesceDropFactory());
  CDFact11->SetFactory("UnAmalgamationInfo", Amalg11);
  RCP<Factory> AggFact11 = rcp(new UncoupledAggregationFactory());
  AggFact11->SetFactory("Graph", CDFact11);
  AggFact11->SetFactory("DofsPerNode", CDFact11);
  RCP<Factory> Cmap11            = rcp(new CoarseMapFactory());
  RCP<TentativePFactory> P11Fact = rcp(new TentativePFactory());
  P11Fact->SetFactory("A", A11Fact);
  P11Fact->SetFactory("UnAmalgamationInfo", Amalg11);
  P11Fact->SetFactory("CoarseMap", Cmap11);
  P11Fact->SetFactory("Aggregates", AggFact11);
  P11Fact->SetParameter("tentative: build coarse coordinates", Teuchos::ParameterEntry(false));
  RCP<TransPFactory> R11Fact = rcp(new TransPFactory());
  R11Fact->SetFactory("P", P11Fact);
  RCP<Factory> Nullspace11 = rcp(new NullspaceFactory());
  Nullspace11->SetParameter("Fine level nullspace", Teuchos::ParameterEntry(std::string("Nullspace1")));
  Nullspace11->SetFactory("Nullspace1", P11Fact);
  RCP<Factory> Coord11 = rcp(new CoordinatesTransferFactory());
  Coord11->SetFactory("Aggregates", AggFact11);
  Coord11->SetFactory("CoarseMap", Cmap11);

  // Create a factory manager for block (0,0)
  RCP<FactoryManager> M11 = rcp(new FactoryManager());
  M11->SetFactory("A", A11Fact);
  M11->SetFactory("UnAmalgamationInfo", Amalg11);
  M11->SetFactory("Aggregates", AggFact11);
  M11->SetFactory("P", P11Fact);
  M11->SetFactory("R", R11Fact);
  M11->SetFactory("Coordinates", Coord11);
  M11->SetFactory("Smoother", Smoo11Fact);
  M11->SetFactory("CoarseMap", Cmap11);
  M11->SetFactory("Nullspace", Nullspace11);
  M11->SetIgnoreUserData(true);

  // Define how the transfer operators for block (1,1) is computed
  RCP<Factory> Amalg22  = rcp(new AmalgamationFactory());
  RCP<Factory> CDFact22 = rcp(new CoalesceDropFactory());
  CDFact22->SetFactory("UnAmalgamationInfo", Amalg22);
  RCP<Factory> AggFact22 = rcp(new UncoupledAggregationFactory());
  AggFact22->SetFactory("Graph", CDFact22);
  AggFact22->SetFactory("DofsPerNode", CDFact22);
  RCP<Factory> Cmap22            = rcp(new BlockedCoarseMapFactory());
  RCP<TentativePFactory> P22Fact = rcp(new TentativePFactory());
  P22Fact->SetFactory("A", A22Fact);
  P22Fact->SetFactory("UnAmalgamationInfo", Amalg22);
  P22Fact->SetFactory("CoarseMap", Cmap22);
  P22Fact->SetFactory("Aggregates", AggFact22);
  P22Fact->SetParameter("tentative: build coarse coordinates", Teuchos::ParameterEntry(false));
  RCP<TransPFactory> R22Fact = rcp(new TransPFactory());
  R22Fact->SetFactory("P", P22Fact);
  RCP<Factory> Nullspace22 = rcp(new NullspaceFactory());
  Nullspace22->SetParameter("Fine level nullspace", Teuchos::ParameterEntry(std::string("Nullspace2")));
  Nullspace22->SetFactory("Nullspace2", P22Fact);
  Cmap22->SetFactory("Aggregates", AggFact22);
  Cmap22->SetFactory("CoarseMap", Cmap11);
  RCP<Factory> Coord22 = rcp(new CoordinatesTransferFactory());
  Coord22->SetFactory("Aggregates", AggFact22);
  Coord22->SetFactory("CoarseMap", Cmap22);

  // Create a factory manager for block (1,1)
  RCP<FactoryManager> M22 = rcp(new FactoryManager());
  M22->SetFactory("A", A22Fact);
  M22->SetFactory("UnAmalgamationInfo", Amalg22);
  M22->SetFactory("Aggregates", AggFact22);
  M22->SetFactory("P", P22Fact);
  M22->SetFactory("R", R22Fact);
  M22->SetFactory("Coordinates", P22Fact);
  M22->SetFactory("Smoother", Smoo22Fact);
  M22->SetFactory("CoarseMap", Cmap22);
  M22->SetFactory("Nullspace", Nullspace22);
  M22->SetIgnoreUserData(true);

  RCP<BlockedPFactory> PFact = rcp(new BlockedPFactory());
  // We might want to set a afactory for "A" here?
  PFact->AddFactoryManager(M11);
  PFact->AddFactoryManager(M22);

  RCP<GenericRFactory> RFact = rcp(new GenericRFactory());
  RFact->SetFactory("P", PFact);

  RCP<BlockedCoordinatesTransferFactory> CoordXfer = rcp(new BlockedCoordinatesTransferFactory());
  CoordXfer->AddFactory(Coord11);
  CoordXfer->AddFactory(Coord22);
  CoordXfer->SetFactory("CoarseMap", PFact);

  RCP<BlockedRAPFactory> AcFact = rcp(new BlockedRAPFactory());
  AcFact->SetFactory("P", PFact);
  AcFact->SetFactory("R", RFact);
  AcFact->AddTransferFactory(CoordXfer);

  // Smoothers
  RCP<BlockedGaussSeidelSmoother> smootherPrototype = rcp(new BlockedGaussSeidelSmoother());
  smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(2));
  smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<SC>(1.0)));
  smootherPrototype->AddFactoryManager(M11, 0);
  smootherPrototype->AddFactoryManager(M22, 1);
  RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

  // Coarse grid correction
  RCP<BlockedGaussSeidelSmoother> coarseSolverPrototype = rcp(new BlockedGaussSeidelSmoother());
  coarseSolverPrototype->AddFactoryManager(M11, 0);
  coarseSolverPrototype->AddFactoryManager(M22, 1);
  RCP<SmootherFactory> coarseSolverFact = rcp(new SmootherFactory(coarseSolverPrototype, Teuchos::null));

  // Repartitioning
  RCP<SubBlockAFactory> AR11Fact = Teuchos::rcp(new SubBlockAFactory());
  AR11Fact->SetFactory("A", AcFact);
  AR11Fact->SetParameter("block row", Teuchos::ParameterEntry(0));
  AR11Fact->SetParameter("block col", Teuchos::ParameterEntry(0));
  RCP<SubBlockAFactory> AR22Fact = Teuchos::rcp(new SubBlockAFactory());
  AR22Fact->SetFactory("A", AcFact);
  AR22Fact->SetParameter("block row", Teuchos::ParameterEntry(1));
  AR22Fact->SetParameter("block col", Teuchos::ParameterEntry(1));

  RCP<RepartitionHeuristicFactory> RepHeuFact = Teuchos::rcp(new RepartitionHeuristicFactory);
  RepHeuFact->SetFactory("A", AcFact);
  RepHeuFact->SetParameter("repartition: start level", Teuchos::ParameterEntry(0));
  RepHeuFact->SetParameter("repartition: min rows per proc", Teuchos::ParameterEntry(40));

#ifdef HAVE_MUELU_ZOLTAN2
  RCP<Zoltan2Interface> zoltan11 = rcp(new Zoltan2Interface());
  zoltan11->SetFactory("A", AR11Fact);
  zoltan11->SetFactory("number of partitions", RepHeuFact);
  zoltan11->SetFactory("Coordinates", Coord11);
#endif

  RCP<RepartitionFactory> repart11 = rcp(new RepartitionFactory());
  repart11->SetFactory("A", AR11Fact);
  repart11->SetFactory("number of partitions", RepHeuFact);
#ifdef HAVE_MUELU_ZOLTAN2
  repart11->SetFactory("Partition", zoltan11);
#endif

#ifdef HAVE_MUELU_ZOLTAN2
  RCP<Zoltan2Interface> zoltan22 = rcp(new Zoltan2Interface());
  zoltan22->SetFactory("A", AR22Fact);
  zoltan22->SetFactory("number of partitions", RepHeuFact);
  zoltan22->SetFactory("Coordinates", Coord22);
#endif

  RCP<RepartitionFactory> repart22 = rcp(new RepartitionFactory());
  repart22->SetFactory("A", AR22Fact);
  repart22->SetFactory("number of partitions", RepHeuFact);
#ifdef HAVE_MUELU_ZOLTAN2
  repart22->SetFactory("Partition", zoltan22);
#endif

  RCP<FactoryManager> MReb11 = rcp(new FactoryManager());
  MReb11->SetFactory("A", AR11Fact);
  MReb11->SetFactory("number of partitions", RepHeuFact);
  MReb11->SetFactory("Importer", repart11);
  MReb11->SetFactory("Nullspace", Nullspace11);

  RCP<FactoryManager> MReb22 = rcp(new FactoryManager());
  MReb22->SetFactory("A", AR22Fact);
  MReb22->SetFactory("number of partitions", RepHeuFact);
  MReb22->SetFactory("Importer", repart22);
  MReb22->SetFactory("Nullspace", Nullspace22);

  RCP<RebalanceBlockInterpolationFactory> RebPFact = rcp(new RebalanceBlockInterpolationFactory());
  RebPFact->SetFactory("A", AcFact);
  RebPFact->SetFactory("P", PFact);
  RebPFact->AddFactoryManager(MReb11);
  RebPFact->AddFactoryManager(MReb22);

  RCP<RebalanceBlockRestrictionFactory> RebRFact = rcp(new RebalanceBlockRestrictionFactory());
  RebRFact->SetFactory("R", RFact);
  RebRFact->AddFactoryManager(MReb11);
  RebRFact->AddFactoryManager(MReb22);

  RCP<RebalanceBlockAcFactory> RebAcFact = rcp(new RebalanceBlockAcFactory());
  RebAcFact->SetFactory("A", AcFact);
  RebAcFact->AddFactoryManager(MReb11);
  RebAcFact->AddFactoryManager(MReb22);

  // main factory manager
  FactoryManager M;
  M.SetFactory("A", RebAcFact);
  M.SetFactory("P", RebPFact);
  M.SetFactory("R", RebRFact);
  M.SetFactory("Coordinates", CoordXfer);
  M.SetFactory("Smoother", smootherFact);  // TODO fix me
  M.SetFactory("CoarseSolver", coarseSolverFact);

  // setup hierarchy
  H.SetVerbLevel(MueLu::Test);
  H.EnableGraphDumping("dep_graph", 0);

  H.Setup(M);
}  // BlockedRAPFactoryCoords

#if defined(HAVE_MUELU_ZOLTAN2) && defined(HAVE_MUELU_ISORROPIA) && defined(HAVE_MUELU_ZOLTAN)
#define MUELU_ETI_GROUP(SC, LO, GO, Node)                                                                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedRepartition, BlockedRAPFactory, SC, LO, GO, Node)                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedRepartition, BlockedRAPFactoryWithRestriction, SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedRepartition, BlockedRAPFactoryWithDiagonal, SC, LO, GO, Node)    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedRepartition, BlockedRAPFactoryCoords, SC, LO, GO, Node)

#elif defined(HAVE_MUELU_ISORROPIA) && defined(HAVE_MUELU_ZOLTAN)
#define MUELU_ETI_GROUP(SC, LO, GO, Node)                                                                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedRepartition, BlockedRAPFactory, SC, LO, GO, Node)                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedRepartition, BlockedRAPFactoryWithRestriction, SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedRepartition, BlockedRAPFactoryWithDiagonal, SC, LO, GO, Node)

#else
#define MUELU_ETI_GROUP(SC, LO, GO, Node)                                                       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedRepartition, BlockedRAPFactory, SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedRepartition, BlockedRAPFactoryWithRestriction, SC, LO, GO, Node)
#endif

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
