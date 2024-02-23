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
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>

#include "MueLu_FactoryManager.hpp"
#include "MueLu_SubBlockAFactory.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLuTests {

/////////////////////////
// helper function

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
GenerateProblemMatrix(const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map,
                      Scalar a = 2.0, Scalar b = -1.0, Scalar c = -1.0) {
#include "MueLu_UseShortNames.hpp"

  Teuchos::RCP<CrsMatrixWrap> mtx = Galeri::Xpetra::MatrixTraits<Map, CrsMatrixWrap>::Build(map, 3);

  LocalOrdinal NumMyElements                               = map->getLocalNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getLocalElementList();
  GlobalOrdinal NumGlobalElements                          = map->getGlobalNumElements();
  GlobalOrdinal nIndexBase                                 = map->getIndexBase();

  GlobalOrdinal NumEntries;
  LocalOrdinal nnz = 2;
  std::vector<Scalar> Values(nnz);
  std::vector<GlobalOrdinal> Indices(nnz);

  for (LocalOrdinal i = 0; i < NumMyElements; ++i) {
    if (MyGlobalElements[i] == nIndexBase) {
      // off-diagonal for first row
      Indices[0] = nIndexBase;
      NumEntries = 1;
      Values[0]  = c;
    } else if (MyGlobalElements[i] == nIndexBase + NumGlobalElements - 1) {
      // off-diagonal for last row
      Indices[0] = nIndexBase + NumGlobalElements - 2;
      NumEntries = 1;
      Values[0]  = b;
    } else {
      // off-diagonal for internal row
      Indices[0] = MyGlobalElements[i] - 1;
      Values[1]  = b;
      Indices[1] = MyGlobalElements[i] + 1;
      Values[0]  = c;
      NumEntries = 2;
    }

    // put the off-diagonal entries
    // Xpetra wants ArrayViews (sigh)
    Teuchos::ArrayView<Scalar> av(&Values[0], NumEntries);
    Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0], NumEntries);
    mtx->insertGlobalValues(MyGlobalElements[i], iv, av);

    // Put in the diagonal entry
    mtx->insertGlobalValues(MyGlobalElements[i],
                            Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
                            Teuchos::tuple<Scalar>(a));
  }

  mtx->fillComplete(map, map);

  return mtx;
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SubBlockAFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);

  out << "version: " << MueLu::Version() << std::endl;

  RCP<SubBlockAFactory> subBlockAFactory = rcp(new SubBlockAFactory());
  TEST_EQUALITY(subBlockAFactory != Teuchos::null, true);
}  // Constructor

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SubBlockAFactory, ExtractMainDiagBlocks, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  /**********************************************************************************/
  /* CREATE INITIAL MATRIX                                                          */
  /**********************************************************************************/
  const GO numElements  = 500;
  const GO numElements1 = 400;
  const GO numElements2 = 100;

  Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(1);

  RCP<const Map> map1 = StridedMapFactory::Build(lib, numElements1, 0, stridingInfo, comm);
  RCP<const Map> map2 = StridedMapFactory::Build(lib, numElements2, numElements1, stridingInfo, comm);

  std::vector<GlobalOrdinal> localGids;                                               // vector with all local GIDs on cur proc
  Teuchos::ArrayView<const GlobalOrdinal> map1eleList = map1->getLocalElementList();  // append all local gids from map1 and map2
  localGids.insert(localGids.end(), map1eleList.begin(), map1eleList.end());
  Teuchos::ArrayView<const GlobalOrdinal> map2eleList = map2->getLocalElementList();
  localGids.insert(localGids.end(), map2eleList.begin(), map2eleList.end());
  Teuchos::ArrayView<GlobalOrdinal> eleList(&localGids[0], localGids.size());
  RCP<const Map> bigMap = MapFactory::Build(lib, numElements, eleList, 0, comm);  // create full big map (concatenation of map1 and map2)

  std::vector<Teuchos::RCP<const Map> > maps;
  maps.push_back(map1);
  maps.push_back(map2);

  RCP<const MapExtractor> mapExtractor = MapExtractorFactory::Build(bigMap, maps);

  RCP<CrsMatrixWrap> Op11 = GenerateProblemMatrix<Scalar, LO, GO, Node>(map1, 2, -1, -1);
  RCP<CrsMatrixWrap> Op22 = GenerateProblemMatrix<Scalar, LO, GO, Node>(map2, 3, -2, -1);

  // build blocked operator
  RCP<BlockedCrsMatrix> bOp = rcp(new BlockedCrsMatrix(mapExtractor, mapExtractor, 10));
  bOp->setMatrix(0, 0, Op11);
  bOp->setMatrix(1, 1, Op22);
  bOp->fillComplete();
  TEST_EQUALITY(bOp != Teuchos::null, true);

  // build hierarchy
  RCP<Level> levelOne = rcp(new Level());
  levelOne->Set("A", Teuchos::rcp_dynamic_cast<Matrix>(bOp));  // set blocked operator

  // define sub block factories for blocked operator "A"
  RCP<SubBlockAFactory> A11Fact = Teuchos::rcp(new SubBlockAFactory());
  A11Fact->SetFactory("A", MueLu::NoFactory::getRCP());
  A11Fact->SetParameter("block row", Teuchos::ParameterEntry(0));
  A11Fact->SetParameter("block col", Teuchos::ParameterEntry(0));
  RCP<SubBlockAFactory> A22Fact = Teuchos::rcp(new SubBlockAFactory());
  A22Fact->SetFactory("A", MueLu::NoFactory::getRCP());
  A22Fact->SetParameter("block row", Teuchos::ParameterEntry(1));
  A22Fact->SetParameter("block col", Teuchos::ParameterEntry(1));

  // Test the request mechanism
  levelOne->Request("A", A11Fact.get(), MueLu::NoFactory::get());
  levelOne->Request("A", A22Fact.get(), MueLu::NoFactory::get());
  TEST_EQUALITY(levelOne->IsRequested("A", A11Fact.get()), true);
  TEST_EQUALITY(levelOne->IsRequested("A", A22Fact.get()), true);

  // Test availability of subblocks of A
  RCP<Matrix> A11 = levelOne->Get<RCP<Matrix> >("A", A11Fact.get());
  RCP<Matrix> A22 = levelOne->Get<RCP<Matrix> >("A", A22Fact.get());
  TEST_EQUALITY(levelOne->IsAvailable("A", A11Fact.get()), true);
  TEST_EQUALITY(levelOne->IsAvailable("A", A22Fact.get()), true);

  // Test the release mechanism
  levelOne->Release("A", A11Fact.get());
  levelOne->Release("A", A22Fact.get());
  TEST_EQUALITY(levelOne->IsAvailable("A", A11Fact.get()), false);
  TEST_EQUALITY(levelOne->IsAvailable("A", A22Fact.get()), false);

  // A11 is supposed to match Op11
  TEST_EQUALITY(A11->getRowMap()->isSameAs(*(Op11->getRowMap())), true);
  TEST_EQUALITY(A11->getColMap()->isSameAs(*(Op11->getColMap())), true);
  TEST_EQUALITY(A11->getRangeMap()->isSameAs(*(Op11->getRangeMap())), true);
  TEST_EQUALITY(A11->getDomainMap()->isSameAs(*(Op11->getDomainMap())), true);
  TEST_EQUALITY(A11->getLocalNumEntries(), Op11->getLocalNumEntries());

  // A22 is supposed to match Op22
  TEST_EQUALITY(A22->getRowMap()->isSameAs(*(Op22->getRowMap())), true);
  TEST_EQUALITY(A22->getColMap()->isSameAs(*(Op22->getColMap())), true);
  TEST_EQUALITY(A22->getRangeMap()->isSameAs(*(Op22->getRangeMap())), true);
  TEST_EQUALITY(A22->getDomainMap()->isSameAs(*(Op22->getDomainMap())), true);
  TEST_EQUALITY(A22->getLocalNumEntries(), Op22->getLocalNumEntries());
}  // ExtractMainDiagBlocks

#define MUELU_ETI_GROUP(SC, LO, GO, Node)                                               \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SubBlockAFactory, Constructor, SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SubBlockAFactory, ExtractMainDiagBlocks, SC, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>
}  // namespace MueLuTests
