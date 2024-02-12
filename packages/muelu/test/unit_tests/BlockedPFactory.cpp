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

#include "MueLu_Utilities.hpp"
#include "MueLu_SubBlockAFactory.hpp"
#include "MueLu_BlockedPFactory.hpp"
#include "MueLu_FactoryManager.hpp"

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
  // std::cout << maxGColId << " " << minGColId << " " << numGColElements <<std::endl;
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
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedPFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);

  out << "version: " << MueLu::Version() << std::endl;

  RCP<BlockedPFactory> blockedPFactory = rcp(new BlockedPFactory());
  TEST_EQUALITY(blockedPFactory != Teuchos::null, true);
}  // Constructor

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedPFactory, CreateBlockedP, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  /**********************************************************************************/
  /* CREATE INITIAL MATRIX                                                          */
  /**********************************************************************************/
  const GO numElements  = 400;
  const GO numElements1 = 200;
  const GO numElements2 = 200;

  Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

  // the test matrix has to be a nxn block matrix with quadratic blocks
  // where the subblocks use consequent numbering of global DOF ids.
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
  RCP<const Map> bigMap = StridedMapFactory::Build(lib, numElements, eleList, 0, stridingInfo, comm);  // create full big map (concatenation of map1 and map2)
  std::vector<Teuchos::RCP<const Map> > maps;
  maps.push_back(map1);
  maps.push_back(map2);

  Teuchos::RCP<const MapExtractor> mapExtractor = MapExtractorFactory::Build(bigMap, maps);

  RCP<CrsMatrixWrap> Op11 = GenerateProblemMatrix<Scalar, LO, GO, Node>(map1, map1, 2, -1, -1);
  RCP<CrsMatrixWrap> Op12 = GenerateProblemMatrix<Scalar, LO, GO, Node>(map1, map2, 1, 0, 0);
  RCP<CrsMatrixWrap> Op21 = GenerateProblemMatrix<Scalar, LO, GO, Node>(map2, map1, 1, 0, 0);
  RCP<CrsMatrixWrap> Op22 = GenerateProblemMatrix<Scalar, LO, GO, Node>(map2, map2, 3, -2, -1);

  // build blocked operator
  Teuchos::RCP<BlockedCrsMatrix> bOp = Teuchos::rcp(new BlockedCrsMatrix(mapExtractor, mapExtractor, 10));
  bOp->setMatrix(0, 0, Op11);
  bOp->setMatrix(0, 1, Op12);
  bOp->setMatrix(1, 0, Op21);
  bOp->setMatrix(1, 1, Op22);
  bOp->fillComplete();
  TEST_EQUALITY(bOp != Teuchos::null, true);

  // build hierarchy
  RCP<Level> levelOne = rcp(new Level());
  RCP<Level> levelTwo = rcp(new Level());
  levelTwo->SetPreviousLevel(levelOne);
  levelOne->Set("A", Teuchos::rcp_dynamic_cast<Matrix>(bOp));  // set blocked operator

  // define sub block factories for blocked operator "A"
  RCP<SubBlockAFactory> A11Fact = Teuchos::rcp(new SubBlockAFactory());
  A11Fact->SetFactory("A", MueLu::NoFactory::getRCP());
  A11Fact->SetParameter("block row", Teuchos::ParameterEntry(0));
  A11Fact->SetParameter("block col", Teuchos::ParameterEntry(0));
  RCP<SubBlockAFactory> A12Fact = Teuchos::rcp(new SubBlockAFactory());
  A12Fact->SetFactory("A", MueLu::NoFactory::getRCP());
  A12Fact->SetParameter("block row", Teuchos::ParameterEntry(0));
  A12Fact->SetParameter("block col", Teuchos::ParameterEntry(1));
  RCP<SubBlockAFactory> A21Fact = Teuchos::rcp(new SubBlockAFactory());
  A21Fact->SetFactory("A", MueLu::NoFactory::getRCP());
  A21Fact->SetParameter("block row", Teuchos::ParameterEntry(1));
  A21Fact->SetParameter("block col", Teuchos::ParameterEntry(0));
  RCP<SubBlockAFactory> A22Fact = Teuchos::rcp(new SubBlockAFactory());
  A22Fact->SetFactory("A", MueLu::NoFactory::getRCP());
  A22Fact->SetParameter("block row", Teuchos::ParameterEntry(1));
  A22Fact->SetParameter("block col", Teuchos::ParameterEntry(1));

  // request subblocks of A
  levelOne->Request("A", A11Fact.get(), MueLu::NoFactory::get());
  levelOne->Request("A", A22Fact.get(), MueLu::NoFactory::get());
  TEST_EQUALITY(levelOne->IsRequested("A", A11Fact.get()), true);
  TEST_EQUALITY(levelOne->IsRequested("A", A22Fact.get()), true);

  RCP<Matrix> A11 = levelOne->Get<RCP<Matrix> >("A", A11Fact.get());
  RCP<Matrix> A22 = levelOne->Get<RCP<Matrix> >("A", A22Fact.get());
  TEST_EQUALITY(levelOne->IsAvailable("A", A11Fact.get()), true);
  TEST_EQUALITY(levelOne->IsAvailable("A", A22Fact.get()), true);

  // store data A11 and A22 as variable "P" on level 2

  levelTwo->Request("P", A11Fact.get(), MueLu::NoFactory::get());
  levelTwo->Request("P", A22Fact.get(), MueLu::NoFactory::get());
  TEST_EQUALITY(levelTwo->IsRequested("P", A11Fact.get()), true);
  TEST_EQUALITY(levelTwo->IsRequested("P", A22Fact.get()), true);

  levelTwo->Set("P", A11, A11Fact.get());
  levelTwo->Set("P", A22, A22Fact.get());
  TEST_EQUALITY(levelTwo->IsAvailable("P", A11Fact.get()), true);
  TEST_EQUALITY(levelTwo->IsAvailable("P", A22Fact.get()), true);

  levelOne->Release("A", A11Fact.get());
  levelOne->Release("A", A22Fact.get());

  // set up factory manager
  RCP<FactoryManager> FC1 = rcp(new FactoryManager());
  FC1->SetFactory("P", A11Fact);  // fool P to be generated by A11Fact
  FC1->SetIgnoreUserData(true);   // always use data from factories defined in factory manager

  RCP<FactoryManager> FC2 = rcp(new FactoryManager());
  FC2->SetFactory("P", A22Fact);  // fool P to be generated by A11Fact
  FC2->SetIgnoreUserData(true);   // always use data from factories defined in factory manager

  /////////////////////////////////////////// define blocked transfer ops
  RCP<BlockedPFactory> PFact = rcp(new BlockedPFactory());  // use row map index base from bOp
  PFact->AddFactoryManager(FC1);
  PFact->AddFactoryManager(FC2);

  levelTwo->Request("P", PFact.get(), MueLu::NoFactory::get());
  TEST_EQUALITY(levelTwo->IsRequested("P", PFact.get()), true);

  RCP<Matrix> P = levelTwo->Get<RCP<Matrix> >("P", PFact.get());
  TEST_EQUALITY(P != Teuchos::null, true);

  RCP<BlockedCrsMatrix> bP = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(P);
  TEST_EQUALITY(bP != Teuchos::null, true);

  TEST_EQUALITY(bP->Rows(), 2);
  TEST_EQUALITY(bP->Cols(), 2);

  // create test and rhs vector
  RCP<const Map> fullMap = mapExtractor->getFullMap();
  TEST_EQUALITY(fullMap == bigMap, true);
  RCP<Vector> iones = VectorFactory::Build(fullMap);
  RCP<Vector> rones = VectorFactory::Build(fullMap);
  iones->putScalar(1.0);
  bP->apply(*iones, *rones);  // the subblocks are chosen, such that bP*ones = zero (except for the first and row and the middle row)
  TEST_EQUALITY(rones->norm1(), 5.0);
  TEST_EQUALITY(rones->normInf(), 2.0);
}  // Constructor

#define MUELU_ETI_GROUP(SC, LO, GO, Node)                                              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedPFactory, Constructor, SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedPFactory, CreateBlockedP, SC, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
