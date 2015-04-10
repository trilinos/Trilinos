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

#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_PreDropFunctionConstVal.hpp"
#include "MueLu_Graph.hpp"

#include "MueLu_UseDefaultTypes.hpp"

namespace MueLuTests {

#include "MueLu_UseShortNames.hpp"

  TEUCHOS_UNIT_TEST(CoalesceDropFactory, Constructor)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<CoalesceDropFactory> coalesceDropFact = rcp(new CoalesceDropFactory());
    TEST_EQUALITY(coalesceDropFact != Teuchos::null, true);

  } //Constructor

  TEUCHOS_UNIT_TEST(CoalesceDropFactory, Build)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    Level fineLevel;
    TestHelpers::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

    RCP<Matrix> A = TestHelpers::TestFactory<SC,LO,GO,NO>::Build1DPoisson(36);
    fineLevel.Set("A", A);

    CoalesceDropFactory coalesceDropFact;
    fineLevel.Request("Graph",&coalesceDropFact);
    fineLevel.Request("DofsPerNode", &coalesceDropFact);

    coalesceDropFact.Build(fineLevel);

    RCP<GraphBase> graph = fineLevel.Get<RCP<GraphBase> >("Graph", &coalesceDropFact);
    LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
    TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

    const RCP<const Map> myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
    const RCP<const Map> myDomainMap = graph->GetDomainMap();

    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myImportMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myImportMap->getGlobalNumElements(),Teuchos::as<size_t>(36 + (comm->getSize()-1)*2));

    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),36);
  } //Build

  // TODO remove this
  TEUCHOS_UNIT_TEST(CoalesceDropFactory, PreDrop)
  {
    out << "version: " << MueLu::Version() << std::endl;

    Level fineLevel;
    TestHelpers::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

    RCP<Matrix> A = TestHelpers::TestFactory<SC,LO,GO,NO>::Build1DPoisson(3);
    fineLevel.Set("A", A);
    A->describe(out,Teuchos::VERB_EXTREME);

    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetVerbLevel(MueLu::Extreme);
    dropFact.SetPreDropFunction(rcp(new PreDropFunctionConstVal(0.00001)));

    fineLevel.Request("Graph", &dropFact);

    dropFact.Build(fineLevel);

    fineLevel.print(out);
    RCP<GraphBase> graph = fineLevel.Get<RCP<GraphBase> >("Graph", &dropFact);

    std::cout << graph->GetDomainMap()->getGlobalNumElements() << std::endl;
    graph->print(out, MueLu::Debug);

    //    TEST_EQUALITY(1 == 0, true);

  } //PreDrop

  TEUCHOS_UNIT_TEST(CoalesceDropFactory, AmalgamationBasic)
  {
    // unit test for block size 3.
    // lightweight wrap = false
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    Level fineLevel;
    TestHelpers::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

    int nx = 3*comm->getSize();
    RCP<Matrix> A = TestHelpers::TestFactory<SC,LO,GO,NO>::Build1DPoisson(nx);
    A->SetFixedBlockSize(3, 0);
    fineLevel.Set("A", A);
    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetParameter("lightweight wrap",Teuchos::ParameterEntry(false));
    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    fineLevel.print(out);
    RCP<GraphBase> graph = fineLevel.Get<RCP<GraphBase> >("Graph", &dropFact);
    LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
    TEST_EQUALITY(Teuchos::as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
    TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 3, true);
    bool bCorrectGraph = false;
    if (comm->getSize() == 1 && graph->getNeighborVertices(0).size() == 1) {
      bCorrectGraph = true;
    } else {
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        if (graph->getNeighborVertices(0).size() == 2) bCorrectGraph = true;
      }
      else {
        if (graph->getNeighborVertices(0).size() == 3) bCorrectGraph = true;
      }
    }
    TEST_EQUALITY(bCorrectGraph, true);

    const RCP<const Map> myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
    const RCP<const Map> myDomainMap = graph->GetDomainMap();

    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myImportMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myImportMap->getGlobalNumElements(),Teuchos::as<size_t>(comm->getSize()+2*(comm->getSize()-1)));
    TEST_EQUALITY(Teuchos::as<bool>(myImportMap->getNodeNumElements()==2) || comm->getSize() == 1, true);
    TEST_EQUALITY(Teuchos::as<bool>(myImportMap->getMaxLocalIndex()==1 || comm->getSize() == 1), true);

    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getMaxLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),Teuchos::as<size_t>(comm->getSize()));
    TEST_EQUALITY(Teuchos::as<bool>(myDomainMap->getNodeNumElements()==1), true);
  } // AmalgamationBasic

  TEUCHOS_UNIT_TEST(CoalesceDropFactory, AmalgamationStrided)
  {
    // unit test for block size 3 using a strided map
    // lightweight wrap = false
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    Level fineLevel;
    TestHelpers::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

    int nx = 3*comm->getSize();
    RCP<Matrix> A = TestHelpers::TestFactory<SC,LO,GO,NO>::Build1DPoisson(nx);

    std::vector<size_t> stridingInfo;
    stridingInfo.push_back(Teuchos::as<size_t>(3));
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
    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetParameter("lightweight wrap",Teuchos::ParameterEntry(false));
    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    fineLevel.print(out);
    RCP<GraphBase> graph = fineLevel.Get<RCP<GraphBase> >("Graph", &dropFact);
    LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
    TEST_EQUALITY(Teuchos::as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
    TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 3, true);
    bool bCorrectGraph = false;
    if (comm->getSize() == 1 && graph->getNeighborVertices(0).size() == 1) {
      bCorrectGraph = true;
    } else {
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        if (graph->getNeighborVertices(0).size() == 2) bCorrectGraph = true;
      }
      else {
        if (graph->getNeighborVertices(0).size() == 3) bCorrectGraph = true;
      }
    }
    TEST_EQUALITY(bCorrectGraph, true);

    const RCP<const Map> myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
    const RCP<const Map> myDomainMap = graph->GetDomainMap();

    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myImportMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myImportMap->getGlobalNumElements(),Teuchos::as<size_t>(comm->getSize()+2*(comm->getSize()-1)));
    TEST_EQUALITY(Teuchos::as<bool>(myImportMap->getNodeNumElements()==2) || comm->getSize() == 1, true);
    TEST_EQUALITY(Teuchos::as<bool>(myImportMap->getMaxLocalIndex()==1 || comm->getSize() == 1), true);

    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getMaxLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),Teuchos::as<size_t>(comm->getSize()));
    TEST_EQUALITY(Teuchos::as<bool>(myDomainMap->getNodeNumElements()==1), true);
  } // AmalgamationStrided

  TEUCHOS_UNIT_TEST(CoalesceDropFactory, AmalgamationStrided2)
  {
    // unit test for block size 3 = (2,1). wrap block 0
    // lightweight wrap = false
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

    // create strided map information
    std::vector<size_t> stridingInfo;
    stridingInfo.push_back(Teuchos::as<size_t>(2));
    stridingInfo.push_back(Teuchos::as<size_t>(1));
    LocalOrdinal stridedBlockId = 0;

    RCP<const StridedMap> dofMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, 3*comm->getSize(), 0,
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
    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetParameter("lightweight wrap",Teuchos::ParameterEntry(false));
    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    fineLevel.print(out);
    RCP<GraphBase> graph = fineLevel.Get<RCP<GraphBase> >("Graph", &dropFact);

    LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
    TEST_EQUALITY(Teuchos::as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
    TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 3, true);
    bool bCorrectGraph = false;
    if (comm->getSize() == 1 && graph->getNeighborVertices(0).size() == 1) {
      bCorrectGraph = true;
    } else {
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        if (graph->getNeighborVertices(0).size() == 2) bCorrectGraph = true;
      }
      else {
        if (graph->getNeighborVertices(0).size() == 3) bCorrectGraph = true;
      }
    }
    TEST_EQUALITY(bCorrectGraph, true);

    const RCP<const Map> myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
    const RCP<const Map> myDomainMap = graph->GetDomainMap();

    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myImportMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myImportMap->getGlobalNumElements(),Teuchos::as<size_t>(comm->getSize()+2*(comm->getSize()-1)));
    TEST_EQUALITY(Teuchos::as<bool>(myImportMap->getNodeNumElements()==2) || comm->getSize() == 1, true);
    TEST_EQUALITY(Teuchos::as<bool>(myImportMap->getMaxLocalIndex()==1 || comm->getSize() == 1), true);

    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getMaxLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),Teuchos::as<size_t>(comm->getSize()));
    TEST_EQUALITY(Teuchos::as<bool>(myDomainMap->getNodeNumElements()==1), true);
  } // AmalgamationStrided2

  TEUCHOS_UNIT_TEST(CoalesceDropFactory, AmalgamationStridedOffset)
  {
    // unit test for block size 9 = (2,3,4). wrap block 1.
    // lightweight wrap = false
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

    // create strided map information
    std::vector<size_t> stridingInfo;
    stridingInfo.push_back(Teuchos::as<size_t>(2));
    stridingInfo.push_back(Teuchos::as<size_t>(3));
    stridingInfo.push_back(Teuchos::as<size_t>(4));
    LocalOrdinal stridedBlockId = 1;
    GlobalOrdinal offset = 19;

    RCP<const StridedMap> dofMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, 9*comm->getSize(), 0,
                                  stridingInfo, comm,
                                  stridedBlockId, offset);

    /////////////////////////////////////////////////////

    Teuchos::RCP<Matrix> mtx = TestHelpers::TestFactory<SC,LO,GO,NO>::BuildTridiag(dofMap, 2.0, -1.0, -3.0);

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
    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetParameter("lightweight wrap",Teuchos::ParameterEntry(false));
    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    fineLevel.print(out);
    RCP<GraphBase> graph = fineLevel.Get<RCP<GraphBase> >("Graph", &dropFact);

    LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
    TEST_EQUALITY(Teuchos::as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
    TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 9, true);
    bool bCorrectGraph = false;
    if (comm->getSize() == 1 && graph->getNeighborVertices(0).size() == 1) {
      bCorrectGraph = true;
    } else {
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        if (graph->getNeighborVertices(0).size() == 2) bCorrectGraph = true;
      }
      else {
        if (graph->getNeighborVertices(0).size() == 3) bCorrectGraph = true;
      }
    }
    TEST_EQUALITY(bCorrectGraph, true);

    const RCP<const Map> myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
    const RCP<const Map> myDomainMap = graph->GetDomainMap();

    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myImportMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myImportMap->getGlobalNumElements(),Teuchos::as<size_t>(comm->getSize()+2*(comm->getSize()-1)));
    TEST_EQUALITY(Teuchos::as<bool>(myImportMap->getNodeNumElements()==2) || comm->getSize() == 1, true);
    TEST_EQUALITY(Teuchos::as<bool>(myImportMap->getMaxLocalIndex()==1 || comm->getSize() == 1), true);

    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getMaxLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),Teuchos::as<size_t>(comm->getSize()));
    TEST_EQUALITY(Teuchos::as<bool>(myDomainMap->getNodeNumElements()==1), true);
  } // AmalgamationStridedOffset

  TEUCHOS_UNIT_TEST(CoalesceDropFactory, AmalgamationLightweight)
  {
    // unit test for block size 3
    // lightweight wrap = true
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

    Level fineLevel;
    TestHelpers::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

    RCP<const Map> dofMap = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, 3*comm->getSize(), 0, comm);
    Teuchos::RCP<Matrix> mtx = TestHelpers::TestFactory<SC,LO,GO,NO>::BuildTridiag(dofMap, 2.0, -1.0, -1.0);
    mtx->SetFixedBlockSize(3, 0);
    fineLevel.Set("A", mtx);

    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetParameter("lightweight wrap",Teuchos::ParameterEntry(true));
    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    RCP<GraphBase> graph = fineLevel.Get<RCP<GraphBase> >("Graph", &dropFact);
    LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
    TEST_EQUALITY(Teuchos::as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
    TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 3, true);
    bool bCorrectGraph = false;
    if (comm->getSize() == 1 && graph->getNeighborVertices(0).size() == 1) {
      bCorrectGraph = true;
    } else {
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        if (graph->getNeighborVertices(0).size() == 2) bCorrectGraph = true;
      }
      else {
        if (graph->getNeighborVertices(0).size() == 3) bCorrectGraph = true;
      }
    }
    TEST_EQUALITY(bCorrectGraph, true);

    const RCP<const Map> myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
    const RCP<const Map> myDomainMap = graph->GetDomainMap();

    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myImportMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myImportMap->getGlobalNumElements(),Teuchos::as<size_t>(comm->getSize()+2*(comm->getSize()-1)));
    TEST_EQUALITY(Teuchos::as<bool>(myImportMap->getNodeNumElements()==2) || comm->getSize() == 1, true);
    TEST_EQUALITY(Teuchos::as<bool>(myImportMap->getMaxLocalIndex()==1 || comm->getSize() == 1), true);

    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getMaxLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),Teuchos::as<size_t>(comm->getSize()));
    TEST_EQUALITY(Teuchos::as<bool>(myDomainMap->getNodeNumElements()==1), true);
  } // AmalgamationLightweight

  TEUCHOS_UNIT_TEST(CoalesceDropFactory, AmalgamationLightweightDrop)
  {
    // unit test for block size 1
    // lightweight wrap = true
    // drop small values
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

    Level fineLevel;
    TestHelpers::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

    RCP<const Map> dofMap = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, 3*comm->getSize(), 0, comm);
    Teuchos::RCP<Matrix> mtx = TestHelpers::TestFactory<SC,LO,GO,NO>::BuildTridiag(dofMap, 1.0, -1.0, -0.0001);
    mtx->SetFixedBlockSize(1, 0);
    fineLevel.Set("A", mtx);

    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetParameter("lightweight wrap",Teuchos::ParameterEntry(true));
    dropFact.SetParameter("aggregation: drop tol",Teuchos::ParameterEntry(0.5));
    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    RCP<GraphBase> graph = fineLevel.Get<RCP<GraphBase> >("Graph", &dropFact);
    LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
    TEST_EQUALITY(Teuchos::as<int>(graph->GetDomainMap()->getGlobalNumElements()) == 3*comm->getSize(), true);
    TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);
    bool bCorrectGraph = false;
    if (comm->getSize() == 1 && graph->getNeighborVertices(0).size() == 1) {
      bCorrectGraph = true;
    } else {
      if (comm->getRank() == 0 ) {
        if (graph->getNeighborVertices(0).size() == 1) bCorrectGraph = true;
      }
      else {
        if (graph->getNeighborVertices(0).size() == 3) bCorrectGraph = true;
      }
    }
    TEST_EQUALITY(bCorrectGraph, true);

    const RCP<const Map> myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
    const RCP<const Map> myDomainMap = graph->GetDomainMap();

    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 3*comm->getSize()-1);
    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myImportMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myImportMap->getGlobalNumElements(),Teuchos::as<size_t>(3*comm->getSize()+(comm->getSize()-1)*2));

    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 3*comm->getSize()-1);
    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),Teuchos::as<size_t>(3*comm->getSize()));
    TEST_EQUALITY(Teuchos::as<bool>(myDomainMap->getNodeNumElements()==3), true);
  } // AmalgamationLightweightDrop

  TEUCHOS_UNIT_TEST(CoalesceDropFactory, AmalgamationStridedLW)
  {
    // unit test for block size 3 using a strided map
    // lightweight wrap = true
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    Level fineLevel;
    TestHelpers::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

    int nx = 3*comm->getSize();
    RCP<Matrix> A = TestHelpers::TestFactory<SC,LO,GO,NO>::Build1DPoisson(nx);

    std::vector<size_t> stridingInfo;
    stridingInfo.push_back(Teuchos::as<size_t>(3));
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
    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetParameter("lightweight wrap",Teuchos::ParameterEntry(true));
    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    fineLevel.print(out);
    RCP<GraphBase> graph = fineLevel.Get<RCP<GraphBase> >("Graph", &dropFact);
    LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
    TEST_EQUALITY(Teuchos::as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
    TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 3, true);
    bool bCorrectGraph = false;
    if (comm->getSize() == 1 && graph->getNeighborVertices(0).size() == 1) {
      bCorrectGraph = true;
    } else {
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        if (graph->getNeighborVertices(0).size() == 2) bCorrectGraph = true;
      }
      else {
        if (graph->getNeighborVertices(0).size() == 3) bCorrectGraph = true;
      }
    }
    TEST_EQUALITY(bCorrectGraph, true);

    const RCP<const Map> myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
    const RCP<const Map> myDomainMap = graph->GetDomainMap();

    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myImportMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myImportMap->getGlobalNumElements(),Teuchos::as<size_t>(comm->getSize()+2*(comm->getSize()-1)));
    TEST_EQUALITY(Teuchos::as<bool>(myImportMap->getNodeNumElements()==2) || comm->getSize() == 1, true);
    TEST_EQUALITY(Teuchos::as<bool>(myImportMap->getMaxLocalIndex()==1 || comm->getSize() == 1), true);

    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getMaxLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),Teuchos::as<size_t>(comm->getSize()));
    TEST_EQUALITY(Teuchos::as<bool>(myDomainMap->getNodeNumElements()==1), true);
  } // AmalgamationStridedLW

  TEUCHOS_UNIT_TEST(CoalesceDropFactory, AmalgamationStrided2LW)
  {
    // unit test for block size 3 = (2,1). wrap block 0
    // lightweight wrap = true
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

    // create strided map information
    std::vector<size_t> stridingInfo;
    stridingInfo.push_back(Teuchos::as<size_t>(2));
    stridingInfo.push_back(Teuchos::as<size_t>(1));
    LocalOrdinal stridedBlockId = 0;

    RCP<const StridedMap> dofMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, 3*comm->getSize(), 0,
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
    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetParameter("lightweight wrap",Teuchos::ParameterEntry(true));
    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    fineLevel.print(out);
    RCP<GraphBase> graph = fineLevel.Get<RCP<GraphBase> >("Graph", &dropFact);

    LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
    TEST_EQUALITY(Teuchos::as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
    TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 3, true);
    bool bCorrectGraph = false;
    if (comm->getSize() == 1 && graph->getNeighborVertices(0).size() == 1) {
      bCorrectGraph = true;
    } else {
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        if (graph->getNeighborVertices(0).size() == 2) bCorrectGraph = true;
      }
      else {
        if (graph->getNeighborVertices(0).size() == 3) bCorrectGraph = true;
      }
    }
    TEST_EQUALITY(bCorrectGraph, true);

    const RCP<const Map> myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
    const RCP<const Map> myDomainMap = graph->GetDomainMap();

    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myImportMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myImportMap->getGlobalNumElements(),Teuchos::as<size_t>(comm->getSize()+2*(comm->getSize()-1)));
    TEST_EQUALITY(Teuchos::as<bool>(myImportMap->getNodeNumElements()==2) || comm->getSize() == 1, true);
    TEST_EQUALITY(Teuchos::as<bool>(myImportMap->getMaxLocalIndex()==1 || comm->getSize() == 1), true);

    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getMaxLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),Teuchos::as<size_t>(comm->getSize()));
    TEST_EQUALITY(Teuchos::as<bool>(myDomainMap->getNodeNumElements()==1), true);
  } // AmalgamationStrided2LW

  TEUCHOS_UNIT_TEST(CoalesceDropFactory, AmalgamationStridedOffsetLW)
  {
    // unit test for block size 9 = (2,3,4). wrap block 1.
    // lightweight wrap = true
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

    // create strided map information
    std::vector<size_t> stridingInfo;
    stridingInfo.push_back(Teuchos::as<size_t>(2));
    stridingInfo.push_back(Teuchos::as<size_t>(3));
    stridingInfo.push_back(Teuchos::as<size_t>(4));
    LocalOrdinal stridedBlockId = 1;
    GlobalOrdinal offset = 19;

    RCP<const StridedMap> dofMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, 9*comm->getSize(), 0,
                                  stridingInfo, comm,
                                  stridedBlockId, offset);

    /////////////////////////////////////////////////////

    Teuchos::RCP<Matrix> mtx = TestHelpers::TestFactory<SC,LO,GO,NO>::BuildTridiag(dofMap, 2.0, -1.0, -3.0);

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
    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetParameter("lightweight wrap",Teuchos::ParameterEntry(true));
    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    fineLevel.print(out);
    RCP<GraphBase> graph = fineLevel.Get<RCP<GraphBase> >("Graph", &dropFact);

    LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
    TEST_EQUALITY(Teuchos::as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
    TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 9, true);
    bool bCorrectGraph = false;
    if (comm->getSize() == 1 && graph->getNeighborVertices(0).size() == 1) {
      bCorrectGraph = true;
    } else {
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        if (graph->getNeighborVertices(0).size() == 2) bCorrectGraph = true;
      }
      else {
        if (graph->getNeighborVertices(0).size() == 3) bCorrectGraph = true;
      }
    }
    TEST_EQUALITY(bCorrectGraph, true);

    const RCP<const Map> myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
    const RCP<const Map> myDomainMap = graph->GetDomainMap();

    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myImportMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myImportMap->getGlobalNumElements(),Teuchos::as<size_t>(comm->getSize()+2*(comm->getSize()-1)));
    TEST_EQUALITY(Teuchos::as<bool>(myImportMap->getNodeNumElements()==2) || comm->getSize() == 1, true);
    TEST_EQUALITY(Teuchos::as<bool>(myImportMap->getMaxLocalIndex()==1 || comm->getSize() == 1), true);

    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getMaxLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),Teuchos::as<size_t>(comm->getSize()));
    TEST_EQUALITY(Teuchos::as<bool>(myDomainMap->getNodeNumElements()==1), true);
  } // AmalgamationStridedOffsetLW

  TEUCHOS_UNIT_TEST(CoalesceDropFactory, AmalgamationDroppingLW)
  {
    // unit test for block size 3 = (3)
    // drop small entries
    // lightweight wrap = true
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

    Level fineLevel;
    TestHelpers::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

    RCP<const Map> dofMap = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, 3*comm->getSize(), 0, comm);
    Teuchos::RCP<Matrix> mtx = TestHelpers::TestFactory<SC,LO,GO,NO>::BuildTridiag(dofMap, 2.0, -1.0, 0.00001);
    mtx->SetFixedBlockSize(3, 0);
    fineLevel.Set("A", mtx);

    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetParameter("lightweight wrap",Teuchos::ParameterEntry(true));
    dropFact.SetParameter("aggregation: drop tol",Teuchos::ParameterEntry(1.0));
    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    fineLevel.print(out);
    RCP<GraphBase> graph = fineLevel.Get<RCP<GraphBase> >("Graph", &dropFact);

    LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
    TEST_EQUALITY(Teuchos::as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
    TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 3, true);
    TEST_EQUALITY(graph->getNeighborVertices(0).size(), 1);

    const RCP<const Map> myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
    const RCP<const Map> myDomainMap = graph->GetDomainMap();

    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myImportMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myImportMap->getGlobalNumElements(),Teuchos::as<size_t>(comm->getSize()+2*(comm->getSize()-1)));
    TEST_EQUALITY(Teuchos::as<bool>(myImportMap->getNodeNumElements()==2) || comm->getSize() == 1, true);
    TEST_EQUALITY(Teuchos::as<bool>(myImportMap->getMaxLocalIndex()==1 || comm->getSize() == 1), true);

    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myDomainMap->getMaxLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),Teuchos::as<size_t>(comm->getSize()));
    TEST_EQUALITY(Teuchos::as<bool>(myDomainMap->getNodeNumElements()==1), true);

  } // AmalgamationDroppingLW

  TEUCHOS_UNIT_TEST(CoalesceDropFactory, AmalgamationStridedOffsetDropping2LW)
  {
    // unit test for block size 9 = (2,3,4). wrap block 1.
    // drop small entries
    // lightweight wrap = true
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

    // create strided map information
    std::vector<size_t> stridingInfo;
    stridingInfo.push_back(Teuchos::as<size_t>(2));
    stridingInfo.push_back(Teuchos::as<size_t>(3));
    stridingInfo.push_back(Teuchos::as<size_t>(4));
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
    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetParameter("lightweight wrap",Teuchos::ParameterEntry(true));
    dropFact.SetParameter("aggregation: drop tol",Teuchos::ParameterEntry(0.3));

    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    fineLevel.print(out);
    RCP<GraphBase> graph = fineLevel.Get<RCP<GraphBase> >("Graph", &dropFact);

    LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
    TEST_EQUALITY(Teuchos::as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
    TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 9, true);
    bool bCorrectGraph = false;
    if (comm->getSize() == 1 && graph->getNeighborVertices(0).size() == 1) {
      bCorrectGraph = true;
    } else {
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        if (graph->getNeighborVertices(0).size() == 1) bCorrectGraph = true;
      }
      else {
        if (graph->getNeighborVertices(0).size() == 3) bCorrectGraph = true;
      }
    }
    TEST_EQUALITY(bCorrectGraph, true);

    const RCP<const Map> myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
    const RCP<const Map> myDomainMap = graph->GetDomainMap();

    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myImportMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myImportMap->getGlobalNumElements(),Teuchos::as<size_t>(comm->getSize()+2*(comm->getSize()-1)));
    TEST_EQUALITY(Teuchos::as<bool>(myImportMap->getNodeNumElements()==2) || comm->getSize() == 1, true);
    TEST_EQUALITY(Teuchos::as<bool>(myImportMap->getMaxLocalIndex()==1 || comm->getSize() == 1), true);

    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),Teuchos::as<size_t>(comm->getSize()));
    TEST_EQUALITY(Teuchos::as<bool>(myDomainMap->getNodeNumElements()==1), true);
  } // AmalgamationStridedOffsetDropping2LW

} // namespace MueLuTests

