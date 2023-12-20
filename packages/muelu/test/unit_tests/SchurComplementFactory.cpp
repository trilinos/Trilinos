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

#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_MapExtractor_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_BlockReorderManager.hpp>
#include <Xpetra_ReorderedBlockedCrsMatrix.hpp>

#include <MueLu_SchurComplementFactory.hpp>
#include <MueLu_InverseApproximationFactory.hpp>

namespace MueLuTests {

/////////////////////////
// helper function

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > CreateBlocked2x2MatrixThyra(const Teuchos::Comm<int>& comm, Xpetra::UnderlyingLib lib) {
  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix;
  typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractor;
  typedef Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> BlockedCrsMatrix;

  std::vector<RCP<const Map> > maps = std::vector<RCP<const Map> >(2, Teuchos::null);
  maps[0]                           = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildMap(100);
  maps[1]                           = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildMap(100);
  RCP<Matrix> A00                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[0], 4.0, -1.0, -1.0, lib);
  RCP<Matrix> A01                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[0], -1.0, 0.0, 0.0, lib);
  RCP<Matrix> A10                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[1], -1.0, 0.0, 0.0, lib);
  RCP<Matrix> A11                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[1], 4.0, -1.0, -1.0, lib);

  // create map extractor
  // To generate the Thyra style map extractor we do not need a full map but only the
  // information about the Map details (i.e. lib and indexBase). We can extract this
  // information from maps[0]
  Teuchos::RCP<const MapExtractor> rgMapExtractor =
      Teuchos::rcp(new MapExtractor(maps[0], maps, true));
  Teuchos::RCP<const MapExtractor> doMapExtractor =
      Teuchos::rcp(new MapExtractor(maps[0], maps, true));
  // build blocked operator
  Teuchos::RCP<BlockedCrsMatrix> bop = Teuchos::rcp(new BlockedCrsMatrix(rgMapExtractor, doMapExtractor, 5));
  bop->setMatrix(Teuchos::as<size_t>(0), Teuchos::as<size_t>(0), A00);
  bop->setMatrix(Teuchos::as<size_t>(0), Teuchos::as<size_t>(1), A01);
  bop->setMatrix(Teuchos::as<size_t>(1), Teuchos::as<size_t>(0), A10);
  bop->setMatrix(Teuchos::as<size_t>(1), Teuchos::as<size_t>(1), A11);
  bop->fillComplete();
  return bop;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > CreateBlocked3x3MatrixThyra(const Teuchos::Comm<int>& comm, Xpetra::UnderlyingLib lib) {
  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix;
  typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractor;
  typedef Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> BlockedCrsMatrix;

  std::vector<RCP<const Map> > maps = std::vector<RCP<const Map> >(3, Teuchos::null);
  maps[0]                           = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildMap(100);
  maps[1]                           = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildMap(100);
  maps[2]                           = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildMap(100);
  RCP<Matrix> A00                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[0], 4.0, -1.0, -1.0, lib);
  RCP<Matrix> A01                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[0], -1.0, 0.0, 0.0, lib);
  RCP<Matrix> A10                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[1], -1.0, 0.0, 0.0, lib);
  RCP<Matrix> A11                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[1], 4.0, -1.0, -1.0, lib);
  RCP<Matrix> A12                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[1], -1.0, 0.0, 0.0, lib);
  RCP<Matrix> A21                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[2], -1.0, 0.0, 0.0, lib);
  RCP<Matrix> A22                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[2], 4.0, -1.0, -1.0, lib);

  // create map extractor
  // To generate the Thyra style map extractor we do not need a full map but only the
  // information about the Map details (i.e. lib and indexBase). We can extract this
  // information from maps[0]
  Teuchos::RCP<const MapExtractor> rgMapExtractor =
      Teuchos::rcp(new MapExtractor(maps[0], maps, true));
  Teuchos::RCP<const MapExtractor> doMapExtractor =
      Teuchos::rcp(new MapExtractor(maps[0], maps, true));
  // build blocked operator
  Teuchos::RCP<BlockedCrsMatrix> bop = Teuchos::rcp(new BlockedCrsMatrix(rgMapExtractor, doMapExtractor, 5));
  bop->setMatrix(Teuchos::as<size_t>(0), Teuchos::as<size_t>(0), A00);
  bop->setMatrix(Teuchos::as<size_t>(0), Teuchos::as<size_t>(1), A01);
  bop->setMatrix(Teuchos::as<size_t>(1), Teuchos::as<size_t>(0), A10);
  bop->setMatrix(Teuchos::as<size_t>(1), Teuchos::as<size_t>(1), A11);
  bop->setMatrix(Teuchos::as<size_t>(1), Teuchos::as<size_t>(2), A12);
  bop->setMatrix(Teuchos::as<size_t>(2), Teuchos::as<size_t>(1), A21);
  bop->setMatrix(Teuchos::as<size_t>(2), Teuchos::as<size_t>(2), A22);
  bop->fillComplete();
  return bop;
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SchurComplementFactory, SchurConstructor2x2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

  Teuchos::RCP<const BlockedCrsMatrix> bop = MueLuTests::CreateBlocked2x2MatrixThyra<Scalar, LocalOrdinal, GlobalOrdinal, Node>(*comm, lib);
  Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
  Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

  Level level;
  TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
  level.Set("A", A);

  RCP<SchurComplementFactory> schurFact = rcp(new SchurComplementFactory());
  schurFact->SetFactory("A", MueLu::NoFactory::getRCP());
  schurFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));

  // request SchurComplement operator
  level.Request("A", schurFact.get());

  // generate Schur complement operator
  schurFact->Build(level);

  RCP<Matrix> sOp = level.Get<RCP<Matrix> >("A", schurFact.get());
  TEST_EQUALITY(sOp.is_null(), false);
  TEST_EQUALITY(sOp->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 100);
  TEST_EQUALITY(sOp->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * 100 + 99);
  TEST_EQUALITY(sOp->getDomainMap()->getMinGlobalIndex(), comm->getRank() * 100);
  TEST_EQUALITY(sOp->getDomainMap()->getMaxGlobalIndex(), comm->getRank() * 100 + 99);
  TEST_EQUALITY(Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(sOp), Teuchos::null);

  /*    size_t localEntries = 0;
  if(comm->getSize() > 2 && comm->getRank() == 0) localEntries = 299;
  else if(comm->getSize() > 2 && comm->getRank() == comm->getSize()-1) localEntries = 299;
  else localEntries = 300;
  if(comm->getSize() == 2) localEntries = 299;
  if(comm->getSize() == 1) localEntries = 298;*/
  //    TEST_EQUALITY(sOp->getLocalNumEntries(), localEntries);
  //    TEST_EQUALITY(sOp->getGlobalNumEntries(),  Teuchos::as<size_t>(comm->getSize() * 300 - 2));

  RCP<Vector> v = VectorFactory::Build(sOp->getRangeMap(), true);
  sOp->getLocalDiagCopy(*v);
  Teuchos::ArrayRCP<const Scalar> vdata = v->getData(0);
  bool bCheck                           = true;
  for (int i = 0; i < 100; i++)
    if (vdata[i] != Teuchos::as<Scalar>(3.75)) bCheck = false;
  TEST_EQUALITY(bCheck, true);

  // check effect of scaling

  schurFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0e12)));
  schurFact->Build(level);

  sOp = level.Get<RCP<Matrix> >("A", schurFact.get());
  TEST_EQUALITY(sOp.is_null(), false);
  TEST_EQUALITY(sOp->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 100);
  TEST_EQUALITY(sOp->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * 100 + 99);
  TEST_EQUALITY(sOp->getDomainMap()->getMinGlobalIndex(), comm->getRank() * 100);
  TEST_EQUALITY(sOp->getDomainMap()->getMaxGlobalIndex(), comm->getRank() * 100 + 99);
  TEST_EQUALITY(Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(sOp), Teuchos::null);

  /*localEntries = 0;
  if(comm->getSize() > 2 && comm->getRank() == 0) localEntries = 299;
  else if(comm->getSize() > 2 && comm->getRank() == comm->getSize()-1) localEntries = 299;
  else localEntries = 300;
  if(comm->getSize() == 2) localEntries = 299;
  if(comm->getSize() == 1) localEntries = 298;*/

  //    TEST_EQUALITY(sOp->getLocalNumEntries(), localEntries);
  //    TEST_EQUALITY(sOp->getGlobalNumEntries(), Teuchos::as<size_t>(comm->getSize() * 300 - 2));

  /*v = VectorFactory::Build(sOp->getRangeMap(),true);
  sOp->getLocalDiagCopy(*v);
  vdata = v->getData(0);
  bCheck = true;
  for(int i=0; i<100; i++) if(Teuchos::ScalarTraits<Scalar>::magnitude(vdata[i] - Teuchos::as<Scalar>(4.0)) < 1.0e-8) { bCheck = false; std::cout << vdata[i] << std::endl; }
  TEST_EQUALITY(bCheck, true);*/
}  // Constructor2x2

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SchurComplementFactory, SchurConstructorXpetra2x2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

  Teuchos::RCP<BlockedCrsMatrix> bop = MueLuTests::TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrix(lib, 2, comm);
  Teuchos::RCP<Matrix> A             = Teuchos::rcp_dynamic_cast<Matrix>(bop);

  Level level;
  TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
  level.Set("A", A);

  RCP<SchurComplementFactory> schurFact = rcp(new SchurComplementFactory());
  schurFact->SetFactory("A", MueLu::NoFactory::getRCP());
  schurFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));

  // request SchurComplement operator
  level.Request("A", schurFact.get());

  // generate Schur complement operator
  schurFact->Build(level);

  RCP<Matrix> sOp = level.Get<RCP<Matrix> >("A", schurFact.get());
  TEST_EQUALITY(sOp.is_null(), false);
  TEST_EQUALITY(sOp->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 10 + 5);
  TEST_EQUALITY(sOp->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * 10 + 9);
  TEST_EQUALITY(sOp->getDomainMap()->getMinGlobalIndex(), comm->getRank() * 10 + 5);
  TEST_EQUALITY(sOp->getDomainMap()->getMaxGlobalIndex(), comm->getRank() * 10 + 9);
  TEST_EQUALITY(Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(sOp), Teuchos::null);

  TEST_EQUALITY(sOp->getLocalNumEntries(), 5);
  TEST_EQUALITY(sOp->getGlobalNumEntries(), Teuchos::as<size_t>(comm->getSize() * 5));

  RCP<Vector> v = VectorFactory::Build(sOp->getRangeMap(), true);
  sOp->getLocalDiagCopy(*v);
  Teuchos::ArrayRCP<const Scalar> vdata = v->getData(0);
  bool bCheck                           = true;
  for (int i = 0; i < 5; i++)
    if (vdata[i] != Teuchos::as<Scalar>(2.0)) bCheck = false;
  TEST_EQUALITY(bCheck, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SchurComplementFactory, SchurConstructor3x3, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

  Teuchos::RCP<const BlockedCrsMatrix> bop = MueLuTests::CreateBlocked3x3MatrixThyra<Scalar, LocalOrdinal, GlobalOrdinal, Node>(*comm, lib);
  Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
  Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

  Level level;
  TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
  level.Set("A", A);

  RCP<InverseApproximationFactory> approxInvFactory = rcp(new InverseApproximationFactory());
  approxInvFactory->SetFactory("A", MueLu::NoFactory::getRCP());
  approxInvFactory->SetParameter("inverse: fixing", Teuchos::ParameterEntry(true));

  RCP<SchurComplementFactory> schurFact = rcp(new SchurComplementFactory());
  schurFact->SetFactory("A", MueLu::NoFactory::getRCP());
  schurFact->SetFactory("Ainv", approxInvFactory);
  schurFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(0.5)));

  // request SchurComplement operator
  level.Request("A", schurFact.get());

  // generate Schur complement operator
  TEST_THROW(schurFact->Build(level), MueLu::Exceptions::RuntimeError);

  // generate reordered blocked operator
  RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 0 [1 2] ]");
  RCP<const ReorderedBlockedCrsMatrix> brop =
      Teuchos::rcp_dynamic_cast<const ReorderedBlockedCrsMatrix>(buildReorderedBlockedCrsMatrix(brm, bop));
  RCP<const Matrix> rA = Teuchos::rcp_dynamic_cast<const Matrix>(brop);
  TEST_EQUALITY(brop->Rows(), 2);
  TEST_EQUALITY(brop->Cols(), 2);
  level.Set("A", Teuchos::rcp_const_cast<Matrix>(rA));

  // try to build Schur complement operator again
  schurFact->Build(level);

  // consecutive Thyra GIDs with offset comm->getSize() * 100!
  RCP<Matrix> sOp = level.Get<RCP<Matrix> >("A", schurFact.get());

  TEST_EQUALITY(sOp.is_null(), false);
  TEST_EQUALITY(sOp->getRangeMap()->getMinGlobalIndex(), comm->getSize() * 100 + comm->getRank() * 100);
  TEST_EQUALITY(sOp->getRangeMap()->getMaxGlobalIndex(), comm->getSize() * 200 + comm->getRank() * 100 + 99);
  TEST_EQUALITY(sOp->getDomainMap()->getMinGlobalIndex(), comm->getSize() * 100 + comm->getRank() * 100);
  TEST_EQUALITY(sOp->getDomainMap()->getMaxGlobalIndex(), comm->getSize() * 200 + comm->getRank() * 100 + 99);
  TEST_EQUALITY(sOp->getRangeMap()->getMinAllGlobalIndex(), comm->getSize() * 100);
  TEST_EQUALITY(sOp->getRangeMap()->getMaxAllGlobalIndex(), comm->getSize() * 300 - 1);
  TEST_EQUALITY(sOp->getDomainMap()->getMinAllGlobalIndex(), comm->getSize() * 100);
  TEST_EQUALITY(sOp->getDomainMap()->getMaxAllGlobalIndex(), comm->getSize() * 300 - 1);

  RCP<BlockedCrsMatrix> bsOp = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(sOp);
  TEST_EQUALITY(bsOp == Teuchos::null, false);
  TEST_EQUALITY(bsOp->Rows(), 2);
  TEST_EQUALITY(bsOp->Cols(), 2);
  TEST_EQUALITY(bsOp->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(bsOp->getDomainMapExtractor()->getThyraMode(), true);

  RCP<Vector> v = VectorFactory::Build(sOp->getRangeMap(), true);
  sOp->getLocalDiagCopy(*v);
  RCP<BlockedVector> bv = Teuchos::rcp_dynamic_cast<BlockedVector>(v);
  TEST_EQUALITY(bv->getBlockedMap()->getNumMaps(), 2);

  RCP<MultiVector> mergedv              = bv->Merge();
  Teuchos::ArrayRCP<const Scalar> vdata = mergedv->getData(0);
  bool bCheck                           = true;
  for (int i = 0; i < 100; i++)
    if (vdata[i] != Teuchos::as<Scalar>(3.50)) bCheck = false;
  for (int i = 100; i < 200; i++)
    if (vdata[i] != Teuchos::as<Scalar>(4.00)) bCheck = false;
  TEST_EQUALITY(bCheck, true);
}  // SchurConstructor3x3

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SchurComplementFactory, SchurConstructorXpetra4x4, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

  Teuchos::RCP<const BlockedCrsMatrix> bop = MueLuTests::TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrix(lib, 4, comm);
  Teuchos::RCP<const Matrix> cA            = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
  Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(cA);

  Level level;
  TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
  level.Set("A", A);

  RCP<SchurComplementFactory> schurFact = rcp(new SchurComplementFactory());
  schurFact->SetFactory("A", MueLu::NoFactory::getRCP());
  schurFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));

  // request SchurComplement operator
  level.Request("A", schurFact.get());

  // generate Schur complement operator
  TEST_THROW(schurFact->Build(level), MueLu::Exceptions::RuntimeError);

  // generate reordered blocked operator
  RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ [0 3] [1 2]]");
  RCP<const ReorderedBlockedCrsMatrix> brop =
      Teuchos::rcp_dynamic_cast<const ReorderedBlockedCrsMatrix>(buildReorderedBlockedCrsMatrix(brm, bop));
  RCP<const Matrix> rA = Teuchos::rcp_dynamic_cast<const Matrix>(brop);
  TEST_EQUALITY(brop->Rows(), 2);
  TEST_EQUALITY(brop->Cols(), 2);
  level.Set("A", Teuchos::rcp_const_cast<Matrix>(rA));

  // generate Schur complement operator
  schurFact->Build(level);
#if 1
  RCP<Matrix> sOp = level.Get<RCP<Matrix> >("A", schurFact.get());

  TEST_EQUALITY(sOp.is_null(), false);
  TEST_EQUALITY(sOp->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 40 + 5);
  TEST_EQUALITY(sOp->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 19);
  TEST_EQUALITY(sOp->getDomainMap()->getMinGlobalIndex(), comm->getRank() * 40 + 5);
  TEST_EQUALITY(sOp->getDomainMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 19);

  TEST_EQUALITY(sOp->getLocalNumEntries(), 15);
  TEST_EQUALITY(sOp->getGlobalNumEntries(), Teuchos::as<size_t>(comm->getSize() * 15));

  RCP<Vector> v = VectorFactory::Build(sOp->getRangeMap(), true);
  sOp->getLocalDiagCopy(*v);
  RCP<BlockedVector> bv = Teuchos::rcp_dynamic_cast<BlockedVector>(v);
  TEST_EQUALITY(bv.is_null(), false);
  TEST_EQUALITY(bv->getBlockedMap()->getNumMaps(), 2);

  RCP<MultiVector> mergedv              = bv->Merge();
  Teuchos::ArrayRCP<const Scalar> vdata = mergedv->getData(0);
  bool bCheck                           = true;
  for (int i = 0; i < 5; i++)
    if (vdata[i] != Teuchos::as<Scalar>(2.0)) bCheck = false;
  for (int i = 5; i < 15; i++)
    if (vdata[i] != Teuchos::as<Scalar>(3.0)) bCheck = false;
  TEST_EQUALITY(bCheck, true);
#endif
}

#define MUELU_ETI_GROUP(SC, LO, GO, Node)                                                                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SchurComplementFactory, SchurConstructor2x2, SC, LO, GO, Node)       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SchurComplementFactory, SchurConstructorXpetra2x2, SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SchurComplementFactory, SchurConstructor3x3, SC, LO, GO, Node)       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SchurComplementFactory, SchurConstructorXpetra4x4, SC, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>
}  // namespace MueLuTests
