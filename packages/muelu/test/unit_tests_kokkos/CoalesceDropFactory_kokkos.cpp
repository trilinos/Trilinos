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

#include "MueLu_TestHelpers_kokkos.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_CoalesceDropFactory_kokkos.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_LWGraph_kokkos.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<CoalesceDropFactory_kokkos> coalesceDropFact = rcp(new CoalesceDropFactory_kokkos());
  TEST_EQUALITY(coalesceDropFact != Teuchos::null, true);

  out << *coalesceDropFact << std::endl;
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, ClassicScalarWithoutFiltering, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory);
  CoalesceDropFactory_kokkos dropFact;
  dropFact.SetFactory("UnAmalgamationInfo", amalgFact);

  fineLevel.Request("Graph", &dropFact);
  fineLevel.Request("DofsPerNode", &dropFact);

  dropFact.Build(fineLevel);

  auto graph         = fineLevel.Get<RCP<LWGraph_kokkos> >("Graph", &dropFact);
  auto myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);

  TEST_EQUALITY(as<int>(myDofsPerNode) == 1, true);

  bool bCorrectGraph = false;
  int reduction_val  = 0;
  int comm_size = comm->getSize(), comm_rank = comm->getRank();
  auto lclLWGraph = graph->getLocalLWGraph();
  Kokkos::parallel_reduce(
      "MueLu:TentativePF:Build:compute_agg_sizes", Kokkos::RangePolicy<typename NO::execution_space, size_t>(0, 1),
      KOKKOS_LAMBDA(const LO i, int &correct) {
        if (comm_size == 1) {
          auto v0 = lclLWGraph.getNeighborVertices(0);
          auto v1 = lclLWGraph.getNeighborVertices(1);
          auto v2 = lclLWGraph.getNeighborVertices(2);
          if (v0.length == 2 && ((v0(0) == 0 && v0(1) == 1) || (v0(0) == 1 && v0(1) == 0)) &&
              v1.length == 3 && v2.length == 3)
            correct = true;
        } else {
          if (comm_rank == 0) {
            if (lclLWGraph.getNeighborVertices(0).length == 2)
              correct = true;

          } else {
            if (lclLWGraph.getNeighborVertices(0).length == 3)
              correct = true;
          }
        }
      },
      reduction_val);
  bCorrectGraph = reduction_val;
  TEST_EQUALITY(bCorrectGraph, true);

  auto myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  auto myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, ClassicScalarWithFiltering, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = TestHelpers_kokkos::Parameters::getLib();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  auto dofMap = MapFactory::Build(lib, 3 * comm->getSize(), 0, comm);
  auto mtx    = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildTridiag(dofMap, 1.0, -1.0, -0.0001);

  mtx->SetFixedBlockSize(1, 0);
  fineLevel.Set("A", mtx);

  RCP<AmalgamationFactory> amalgFact  = rcp(new AmalgamationFactory);
  CoalesceDropFactory_kokkos dropFact = CoalesceDropFactory_kokkos();
  dropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  dropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.5));

  fineLevel.Request("Graph", &dropFact);
  fineLevel.Request("DofsPerNode", &dropFact);

  dropFact.Build(fineLevel);

  auto graph         = fineLevel.Get<RCP<LWGraph_kokkos> >("Graph", &dropFact);
  auto myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);

  TEST_EQUALITY(as<int>(myDofsPerNode) == 1, true);
  TEST_EQUALITY(as<int>(graph->GetDomainMap()->getGlobalNumElements()) == 3 * comm->getSize(), true);

  bool bCorrectGraph = false;
  int reduction_val  = 0;
  int comm_size = comm->getSize(), comm_rank = comm->getRank();
  auto lclLWGraph = graph->getLocalLWGraph();
  Kokkos::parallel_reduce(
      "MueLu:TentativePF:Build:compute_agg_sizes", Kokkos::RangePolicy<typename NO::execution_space, size_t>(0, 1),
      KOKKOS_LAMBDA(const LO i, int &correct) {
        if (comm_size == 1) {
          auto v0 = lclLWGraph.getNeighborVertices(0);
          auto v1 = lclLWGraph.getNeighborVertices(1);
          auto v2 = lclLWGraph.getNeighborVertices(2);
          if (v0.length == 1 && v0(0) == 0 &&
              v1.length == 2 && ((v1(0) == 0 && v1(1) == 1) || (v1(0) == 1 && v1(1) == 0)) &&
              v2.length == 2 && ((v2(0) == 1 && v2(1) == 2) || (v2(0) == 2 && v2(1) == 1)))
            correct = true;
        } else {
          if (comm_rank == 0) {
            if (lclLWGraph.getNeighborVertices(0).length == 1)
              correct = true;

          } else {
            if (lclLWGraph.getNeighborVertices(0).length == 2)
              correct = true;
          }
        }
      },
      reduction_val);
  bCorrectGraph = reduction_val;
  TEST_EQUALITY(bCorrectGraph, true);

  auto myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  auto myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 3 * comm->getSize() - 1);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), as<size_t>(3 * comm->getSize() + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 3 * comm->getSize() - 1);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), as<size_t>(3 * comm->getSize()));
  TEST_EQUALITY(myDomainMap->getLocalNumElements(), 3);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, ClassicBlockWithoutFiltering, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = TestHelpers_kokkos::Parameters::getLib();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  int blockSize = 3;

  auto dofMap = MapFactory::Build(lib, blockSize * comm->getSize(), 0, comm);
  auto mtx    = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildTridiag(dofMap, 2.0, -1.0, -1.0);
  mtx->SetFixedBlockSize(blockSize, 0);
  fineLevel.Set("A", mtx);

  RCP<AmalgamationFactory> amalgFact  = rcp(new AmalgamationFactory);
  CoalesceDropFactory_kokkos dropFact = CoalesceDropFactory_kokkos();
  dropFact.SetFactory("UnAmalgamationInfo", amalgFact);

  fineLevel.Request("Graph", &dropFact);
  fineLevel.Request("DofsPerNode", &dropFact);

  dropFact.Build(fineLevel);

  auto graph         = fineLevel.Get<RCP<LWGraph_kokkos> >("Graph", &dropFact);
  auto myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);

  TEST_EQUALITY(as<int>(myDofsPerNode) == blockSize, true);
  TEST_EQUALITY(as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);

  bool bCorrectGraph = false;
  int reduction_val  = 0;
  int comm_size = comm->getSize(), comm_rank = comm->getRank();
  auto lclLWGraph = graph->getLocalLWGraph();
  Kokkos::parallel_reduce(
      "MueLu:TentativePF:Build:compute_agg_sizes", Kokkos::RangePolicy<typename NO::execution_space, size_t>(0, 1),
      KOKKOS_LAMBDA(const LO i, int &correct) {
        if (comm_size == 1 && lclLWGraph.getNeighborVertices(0).length == 1) {
          correct = true;
        } else {
          if (comm_rank == 0 || comm_rank == comm_size - 1) {
            if (lclLWGraph.getNeighborVertices(0).length == 2)
              correct = true;

          } else {
            if (static_cast<int>(lclLWGraph.getNeighborVertices(0).length) == blockSize)
              correct = true;
          }
        }
      },
      reduction_val);
  bCorrectGraph = reduction_val;
  TEST_EQUALITY(bCorrectGraph, true);

  auto myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  auto myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), as<size_t>(comm->getSize() + 2 * (comm->getSize() - 1)));

  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t numLocalImportElts = myImportMap->getLocalNumElements();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(numLocalImportElts, numLocalRowMapElts + 1);
    } else {
      TEST_EQUALITY(numLocalImportElts, numLocalRowMapElts + 2);
    }
  }
  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t maxLocalIndex      = myImportMap->getMaxLocalIndex();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(maxLocalIndex, numLocalRowMapElts * blockSize - 2);
    } else {
      TEST_EQUALITY(maxLocalIndex, numLocalRowMapElts * blockSize - 1);
    }
  }

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMaxLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), as<size_t>(comm->getSize()));
  TEST_EQUALITY(myDomainMap->getLocalNumElements(), 1);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, ClassicBlockWithFiltering, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = TestHelpers_kokkos::Parameters::getLib();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  auto dofMap = MapFactory::Build(lib, 3 * comm->getSize(), 0, comm);
  auto mtx    = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildTridiag(dofMap, 2.0, -1.0, 0.00001);

  mtx->SetFixedBlockSize(3, 0);
  fineLevel.Set("A", mtx);

  RCP<AmalgamationFactory> amalgFact  = rcp(new AmalgamationFactory);
  CoalesceDropFactory_kokkos dropFact = CoalesceDropFactory_kokkos();
  dropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  dropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(1.0));

  fineLevel.Request("Graph", &dropFact);
  fineLevel.Request("DofsPerNode", &dropFact);

  dropFact.Build(fineLevel);

  auto graph         = fineLevel.Get<RCP<LWGraph_kokkos> >("Graph", &dropFact);
  auto myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);

  TEST_EQUALITY(as<int>(myDofsPerNode) == 3, true);
  TEST_EQUALITY(as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);

  TEST_EQUALITY(graph->getNeighborVertices(0).size(), 1);

  auto myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  auto myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), as<size_t>(comm->getSize() + 2 * (comm->getSize() - 1)));
  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t numLocalImportElts = myImportMap->getLocalNumElements();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(numLocalImportElts, numLocalRowMapElts + 1);
    } else {
      TEST_EQUALITY(numLocalImportElts, numLocalRowMapElts + 2);
    }
  }

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMaxLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), as<size_t>(comm->getSize()));
  TEST_EQUALITY(myDomainMap->getLocalNumElements(), 1);
}

#if 0
  TEUCHOS_UNIT_TEST(CoalesceDropFactory_kokkos, LaplacianScalarWithoutFiltering)
  {
#include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    Level fineLevel;
    TestHelpers_kokkos::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

    RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC,LO,GO,NO>::Build1DPoisson(36);
    fineLevel.Set("A", A);

    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory);
    CoalesceDropFactory_kokkos dropFact;
    dropFact.SetFactory("UnAmalgamationInfo", amalgFact);
    dropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("distance laplacian")));

    fineLevel.Request("Graph",       &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    auto graph         = fineLevel.Get<RCP<LWGraph_kokkos> >("Graph",       &dropFact);
    auto myDofsPerNode = fineLevel.Get<LO>                  ("DofsPerNode", &dropFact);
    TEST_EQUALITY(as<int>(myDofsPerNode) == 1, true);

    bool bCorrectGraph = false;
    if (comm->getSize() == 1) {
      auto v0 = graph->getNeighborVertices(0);
      auto v1 = graph->getNeighborVertices(1);
      auto v2 = graph->getNeighborVertices(2);
      if (v0.size() == 2 && ((v0(0) == 0 && v0(1) == 1) || (v0(0) == 1 && v0(1) == 0)) &&
          v1.size() == 3 && v2.size() == 3)
        bCorrectGraph = true;
    } else {
      if (comm->getRank() == 0 ) {
        if (graph->getNeighborVertices(0).size() == 2)
          bCorrectGraph = true;

      } else {
        if (graph->getNeighborVertices(0).size() == 3)
          bCorrectGraph = true;
      }
    }
    TEST_EQUALITY(bCorrectGraph, true);

    auto myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
    auto myDomainMap = graph->GetDomainMap();

    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(),  35);
    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(),  0);
    TEST_EQUALITY(myImportMap->getMinLocalIndex(),      0);
    TEST_EQUALITY(myImportMap->getGlobalNumElements(),  as<size_t>(36 + (comm->getSize()-1)*2));

    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(),  35);
    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(),  0);
    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),      0);
    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),  36);
  }
#endif

#if 0
  TEUCHOS_UNIT_TEST(CoalesceDropFactory, AmalgamationStridedLW)
  {
#include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);
    out << "version: " << MueLu::Version() << std::endl;

    // unit test for block size 3 using a strided map

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    Level fineLevel;
    TestHelpers::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

    int blockSize=3;

    GO nx = blockSize*comm->getSize();
    RCP<Matrix> A = TestHelpers::TestFactory<SC,LO,GO,NO>::Build1DPoisson(nx);

    std::vector<size_t> stridingInfo;
    stridingInfo.push_back(as<size_t>(blockSize));
    LocalOrdinal stridedBlockId = -1;

    RCP<const Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node> > stridedRangeMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
                                                  A->getRangeMap(),
                                                  stridingInfo,
                                                  stridedBlockId,
                                                  0 /*offset*/
                                                  );
    RCP<const Map> stridedDomainMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
                                            A->getDomainMap(),
                                            stridingInfo,
                                            stridedBlockId,
                                            0 /*offset*/
                                            );

    if(A->IsView("stridedMaps") == true) A->RemoveView("stridedMaps");
    A->CreateView("stridedMaps", stridedRangeMap, stridedDomainMap);

    fineLevel.Set("A", A);
    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory);
    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetFactory("UnAmalgamationInfo", amalgFact);
    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    fineLevel.print(out);
    RCP<GraphBase> graph = fineLevel.Get<RCP<GraphBase> >("Graph", &dropFact);
    LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
    TEST_EQUALITY(as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
    TEST_EQUALITY(as<int>(myDofsPerNode) == blockSize, true);
    bool bCorrectGraph = false;
    if (comm->getSize() == 1 && graph->getNeighborVertices(0).size() == 1) {
      bCorrectGraph = true;
    } else {
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        if (graph->getNeighborVertices(0).size() == 2) bCorrectGraph = true;
      }
      else {
        if (graph->getNeighborVertices(0).size() == blockSize) bCorrectGraph = true;
      }
    }
    TEST_EQUALITY(bCorrectGraph, true);

    const RCP<const Map> myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
    const RCP<const Map> myDomainMap = graph->GetDomainMap();

    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myImportMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myImportMap->getGlobalNumElements(),as<size_t>(comm->getSize()+2*(comm->getSize()-1)));
    if (comm->getSize()>1) {
      size_t numLocalRowMapElts = graph->GetNodeNumVertices();
      size_t numLocalImportElts = myImportMap->getLocalNumElements();
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        TEST_EQUALITY(as<bool>(numLocalImportElts==numLocalRowMapElts+1), true);
      } else {
        TEST_EQUALITY(as<bool>(numLocalImportElts==numLocalRowMapElts+2), true);
      }
    }
    if (comm->getSize()>1) {
      size_t numLocalRowMapElts = graph->GetNodeNumVertices();
      size_t maxLocalIndex = myImportMap->getMaxLocalIndex();
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        TEST_EQUALITY(as<bool>(maxLocalIndex==numLocalRowMapElts*blockSize-2), true);
      } else {
        TEST_EQUALITY(as<bool>(maxLocalIndex==numLocalRowMapElts*blockSize-1), true);
      }
    }

    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getMaxLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),as<size_t>(comm->getSize()));
    TEST_EQUALITY(as<bool>(myDomainMap->getLocalNumElements()==1), true);
  } // AmalgamationStridedLW

  TEUCHOS_UNIT_TEST(CoalesceDropFactory, AmalgamationStrided2LW)
  {
#include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);
    out << "version: " << MueLu::Version() << std::endl;

    // unit test for block size 3 = (2,1). wrap block 0

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

    // create strided map information
    std::vector<size_t> stridingInfo;
    stridingInfo.push_back(as<size_t>(2));
    stridingInfo.push_back(as<size_t>(1));
    LocalOrdinal stridedBlockId = 0;

    int blockSize=3;

    RCP<const StridedMap> dofMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, blockSize*comm->getSize(), 0,
                                  stridingInfo, comm,
                                  stridedBlockId /*blockId*/, 0 /*offset*/);

    /////////////////////////////////////////////////////

    Teuchos::RCP<Matrix> mtx = TestHelpers::TestFactory<SC,LO,GO,NO>::BuildTridiag(dofMap, 2.0, -1.0, -1.0);

    Level fineLevel;
    TestHelpers::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

    RCP<const Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node> > stridedRangeMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
                                                  mtx->getRangeMap(),
                                                  stridingInfo,
                                                  stridedBlockId,
                                                  0 /*offset*/
                                                  );
    RCP<const Map> stridedDomainMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
                                            mtx->getDomainMap(),
                                            stridingInfo,
                                            stridedBlockId,
                                            0 /*offset*/
                                            );
    if(mtx->IsView("stridedMaps") == true) mtx->RemoveView("stridedMaps");
    mtx->CreateView("stridedMaps", stridedRangeMap, stridedDomainMap);

    fineLevel.Set("A", mtx);
    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory);
    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetFactory("UnAmalgamationInfo", amalgFact);
    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    fineLevel.print(out);
    RCP<GraphBase> graph = fineLevel.Get<RCP<GraphBase> >("Graph", &dropFact);

    LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
    TEST_EQUALITY(as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
    TEST_EQUALITY(as<int>(myDofsPerNode) == blockSize, true);
    bool bCorrectGraph = false;
    if (comm->getSize() == 1 && graph->getNeighborVertices(0).size() == 1) {
      bCorrectGraph = true;
    } else {
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        if (graph->getNeighborVertices(0).size() == 2) bCorrectGraph = true;
      }
      else {
        if (graph->getNeighborVertices(0).size() == blockSize) bCorrectGraph = true;
      }
    }
    TEST_EQUALITY(bCorrectGraph, true);

    const RCP<const Map> myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
    const RCP<const Map> myDomainMap = graph->GetDomainMap();

    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myImportMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myImportMap->getGlobalNumElements(),as<size_t>(comm->getSize()+2*(comm->getSize()-1)));
    if (comm->getSize()>1) {
      size_t numLocalRowMapElts = graph->GetNodeNumVertices();
      size_t numLocalImportElts = myImportMap->getLocalNumElements();
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        TEST_EQUALITY(as<bool>(numLocalImportElts==numLocalRowMapElts+1), true);
      } else {
        TEST_EQUALITY(as<bool>(numLocalImportElts==numLocalRowMapElts+2), true);
      }
    }
    if (comm->getSize()>1) {
      size_t numLocalRowMapElts = graph->GetNodeNumVertices();
      size_t maxLocalIndex = myImportMap->getMaxLocalIndex();
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        TEST_EQUALITY(as<bool>(maxLocalIndex==numLocalRowMapElts*blockSize-2), true);
      } else {
        TEST_EQUALITY(as<bool>(maxLocalIndex==numLocalRowMapElts*blockSize-1), true);
      }
    }

    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getMaxLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),as<size_t>(comm->getSize()));
    TEST_EQUALITY(as<bool>(myDomainMap->getLocalNumElements()==1), true);
  } // AmalgamationStrided2LW


  TEUCHOS_UNIT_TEST(CoalesceDropFactory, AmalgamationStridedOffsetDropping2LW)
  {
    // unit test for block size 9 = (2,3,4). wrap block 1.
    // drop small entries
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

    // create strided map information
    std::vector<size_t> stridingInfo;
    stridingInfo.push_back(as<size_t>(2));
    stridingInfo.push_back(as<size_t>(3));
    stridingInfo.push_back(as<size_t>(4));
    LocalOrdinal stridedBlockId = 1;
    GlobalOrdinal offset = 19;

    RCP<const StridedMap> dofMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, 9*comm->getSize(), 0,
                                  stridingInfo, comm,
                                  stridedBlockId, offset);

    /////////////////////////////////////////////////////

    Teuchos::RCP<Matrix> mtx = TestHelpers::TestFactory<SC,LO,GO,NO>::BuildTridiag(dofMap, 2.0, 1.0, 0.0001);

    Level fineLevel;
    TestHelpers::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

    RCP<const Map> stridedRangeMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
                                                  mtx->getRangeMap(),
                                                  stridingInfo,
                                                  stridedBlockId,
                                                  offset
                                                  );
    RCP<const Map> stridedDomainMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
                                            mtx->getDomainMap(),
                                            stridingInfo,
                                            stridedBlockId,
                                            offset
                                            );

    if(mtx->IsView("stridedMaps") == true) mtx->RemoveView("stridedMaps");
    mtx->CreateView("stridedMaps", stridedRangeMap, stridedDomainMap);

    fineLevel.Set("A", mtx);
    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory);
    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetFactory("UnAmalgamationInfo", amalgFact);
    dropFact.SetParameter("aggregation: drop tol",Teuchos::ParameterEntry(0.3));

    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    fineLevel.print(out);
    RCP<GraphBase> graph = fineLevel.Get<RCP<GraphBase> >("Graph", &dropFact);

    LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
    TEST_EQUALITY(as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
    TEST_EQUALITY(as<int>(myDofsPerNode) == 9, true);
    bool bCorrectGraph = false;
    if (comm->getSize() == 1 && graph->getNeighborVertices(0).size() == 1) {
      bCorrectGraph = true;
    } else {
      if (comm->getRank() == 0) {
        if (graph->getNeighborVertices(0).size() == 1) bCorrectGraph = true;
      }
      else {
        if (graph->getNeighborVertices(0).size() == 2) bCorrectGraph = true;
      }
    }
    TEST_EQUALITY(bCorrectGraph, true);

    const RCP<const Map> myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
    const RCP<const Map> myDomainMap = graph->GetDomainMap();

    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myImportMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myImportMap->getGlobalNumElements(),as<size_t>(comm->getSize()+2*(comm->getSize()-1)));
    if (comm->getSize()>1) {
      size_t numLocalRowMapElts = graph->GetNodeNumVertices();
      size_t numLocalImportElts = myImportMap->getLocalNumElements();
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        TEST_EQUALITY(as<bool>(numLocalImportElts==numLocalRowMapElts+1), true);
      } else {
        TEST_EQUALITY(as<bool>(numLocalImportElts==numLocalRowMapElts+2), true);
      }
    }
    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),as<size_t>(comm->getSize()));
    TEST_EQUALITY(as<bool>(myDomainMap->getLocalNumElements()==1), true);
  } // AmalgamationStridedOffsetDropping2LW
#endif

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, AggresiveDroppingIsMarkedAsBoundary, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  // Test that when everything but the diagonal is dropped, the node is marked as boundary
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = TestHelpers_kokkos::Parameters::getLib();

  RCP<const Map> dofMap    = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, 12 * comm->getSize(), 0, comm);
  Teuchos::RCP<Matrix> mtx = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildTridiag(dofMap, 2.0, -1.0, -1.0);

  {
    Level fineLevel;
    TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

    mtx->SetFixedBlockSize(1);
    fineLevel.Set("A", mtx);

    CoalesceDropFactory_kokkos dropFact = CoalesceDropFactory_kokkos();
    RCP<AmalgamationFactory> amalgFact  = rcp(new AmalgamationFactory());
    dropFact.SetFactory("UnAmalgamationInfo", amalgFact);
    dropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(4.1));
    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos> >("Graph", &dropFact);

    auto boundaryNodes     = graph->getLocalLWGraph().GetBoundaryNodeMap();
    auto boundaryNodesHost = Kokkos::create_mirror_view(boundaryNodes);
    Kokkos::deep_copy(boundaryNodesHost, boundaryNodes);
    bool allNodesAreOnBoundary = true;
    for (LO i = 0; i < Teuchos::as<LO>(boundaryNodesHost.size()); i++)
      allNodesAreOnBoundary &= boundaryNodesHost(i);
    TEST_EQUALITY(allNodesAreOnBoundary, true);
  }

  // {
  //   Level fineLevel;
  //   TestHelpers_kokkos::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

  //   mtx->SetFixedBlockSize(2);
  //   fineLevel.Set("A", mtx);

  //   CoalesceDropFactory_kokkos dropFact = CoalesceDropFactory_kokkos();
  //   RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  //   dropFact.SetFactory("UnAmalgamationInfo",amalgFact);
  //   dropFact.SetParameter("aggregation: drop tol",Teuchos::ParameterEntry(4.1));
  //   fineLevel.Request("Graph", &dropFact);
  //   fineLevel.Request("DofsPerNode", &dropFact);

  //   dropFact.Build(fineLevel);

  //   RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos> >("Graph", &dropFact);

  //   auto boundaryNodes = graph->GetBoundaryNodeMap();
  //   bool allNodesAreOnBoundary = true;
  //   for (LO i = 0; i < Teuchos::as<LO>(boundaryNodes.size()); i++)
  //     allNodesAreOnBoundary &= boundaryNodes[i];
  //   TEST_EQUALITY(allNodesAreOnBoundary, true);
  // }

  // {
  //   Level fineLevel;
  //   TestHelpers_kokkos::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

  //   mtx->SetFixedBlockSize(3);
  //   fineLevel.Set("A", mtx);

  //   CoalesceDropFactory_kokkos dropFact = CoalesceDropFactory_kokkos();
  //   RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  //   dropFact.SetFactory("UnAmalgamationInfo",amalgFact);
  //   dropFact.SetParameter("aggregation: drop tol",Teuchos::ParameterEntry(4.1));
  //   fineLevel.Request("Graph", &dropFact);
  //   fineLevel.Request("DofsPerNode", &dropFact);

  //   dropFact.Build(fineLevel);

  //   RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos> >("Graph", &dropFact);

  //   auto boundaryNodes = graph->GetBoundaryNodeMap();
  //   bool allNodesAreOnBoundary = true;
  //   for (LO i = 0; i < Teuchos::as<LO>(boundaryNodes.size()); i++)
  //     allNodesAreOnBoundary &= boundaryNodes[i];
  //   TEST_EQUALITY(allNodesAreOnBoundary, true);
  // }

}  // AggresiveDroppingIsMarkedAsBoundary

#define MUELU_ETI_GROUP(SC, LO, GO, NO)                                                                           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, Constructor, SC, LO, GO, NO)                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, ClassicScalarWithoutFiltering, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, ClassicScalarWithFiltering, SC, LO, GO, NO)    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, ClassicBlockWithoutFiltering, SC, LO, GO, NO)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, AggresiveDroppingIsMarkedAsBoundary, SC, LO, GO, NO)

// TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, ClassicBlockWithFiltering,     SC, LO, GO, NO) // not implemented yet

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
