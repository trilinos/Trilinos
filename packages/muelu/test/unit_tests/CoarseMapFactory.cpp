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
#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <Xpetra_MultiVectorFactory.hpp>

#include <MueLu_FactoryManagerBase.hpp>
#include <MueLu_CoalesceDropFactory.hpp>
#include <MueLu_AmalgamationFactory.hpp>
#include <MueLu_AmalgamationInfo.hpp>
#include <MueLu_Aggregates.hpp>
#include <MueLu_CoarseMapFactory.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoarseMapFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  out << "version: " << MueLu::Version() << std::endl;

  RCP<CoarseMapFactory> myCMF = rcp(new CoarseMapFactory());
  TEST_ASSERT(!myCMF.is_null());
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoarseMapFactory, StandardCase, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  Level myLevel;
  myLevel.SetLevelID(0);
  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(15);
  myLevel.Set("A", A);

  // build dummy aggregate structure
  Teuchos::RCP<Aggregates> aggs = Teuchos::rcp(new Aggregates(A->getRowMap()));
  aggs->SetNumAggregates(10);  // set (local!) number of aggregates
  myLevel.Set("Aggregates", aggs);

  // build dummy nullspace vector
  Teuchos::RCP<MultiVector> nsp = MultiVectorFactory::Build(A->getRowMap(), 1);
  nsp->putScalar(1.0);
  myLevel.Set("Nullspace", nsp);

  RCP<CoarseMapFactory> myCMF = Teuchos::rcp(new CoarseMapFactory());
  myLevel.Request("CoarseMap", myCMF.get());
  myCMF->SetFactory("Aggregates", MueLu::NoFactory::getRCP());
  myCMF->SetFactory("Nullspace", MueLu::NoFactory::getRCP());
  myCMF->Build(myLevel);
  Teuchos::RCP<const Map> myCoarseMap = myLevel.Get<Teuchos::RCP<const Map> >("CoarseMap", myCMF.get());

  TEST_EQUALITY(myCoarseMap->getMinAllGlobalIndex() == 0, true);
  TEST_EQUALITY(myCoarseMap->getMaxLocalIndex() == 9, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoarseMapFactory, NonStandardCaseA, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  Level myLevel;
  myLevel.SetLevelID(0);
  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(15);
  myLevel.Set("A", A);

  // build dummy aggregate structure
  Teuchos::RCP<Aggregates> aggs = Teuchos::rcp(new Aggregates(A->getRowMap()));
  aggs->SetNumAggregates(10);  // set (local!) number of aggregates
  myLevel.Set("Aggregates", aggs);

  // build dummy nullspace vector
  Teuchos::RCP<MultiVector> nsp = MultiVectorFactory::Build(A->getRowMap(), 1);
  nsp->putScalar(1.0);
  myLevel.Set("Nullspace", nsp);

  RCP<CoarseMapFactory> myCMF = Teuchos::rcp(new CoarseMapFactory());
  myLevel.Request("CoarseMap", myCMF.get());
  myCMF->SetParameter("Domain GID offsets", Teuchos::ParameterEntry(std::string("{100,50}")));
  myCMF->SetFactory("Aggregates", MueLu::NoFactory::getRCP());
  myCMF->SetFactory("Nullspace", MueLu::NoFactory::getRCP());
  myCMF->Build(myLevel);
  Teuchos::RCP<const Map> myCoarseMap = myLevel.Get<Teuchos::RCP<const Map> >("CoarseMap", myCMF.get());

  TEST_EQUALITY(myCoarseMap->getMinAllGlobalIndex() == 100, true);
  TEST_EQUALITY(myCoarseMap->getMaxLocalIndex() == 9, true);

  myLevel.Release("CoarseMap", myCMF.get());
  myLevel.SetLevelID(1);
  myLevel.Request("CoarseMap", myCMF.get());
  myCMF->SetParameter("Domain GID offsets", Teuchos::ParameterEntry(std::string("{100,50}")));
  myCMF->SetFactory("Aggregates", MueLu::NoFactory::getRCP());
  myCMF->SetFactory("Nullspace", MueLu::NoFactory::getRCP());
  myCMF->Build(myLevel);
  myCoarseMap = myLevel.Get<Teuchos::RCP<const Map> >("CoarseMap", myCMF.get());

  TEST_EQUALITY(myCoarseMap->getMinAllGlobalIndex() == 50, true);
  TEST_EQUALITY(myCoarseMap->getMaxLocalIndex() == 9, true);

  myLevel.Release("CoarseMap", myCMF.get());
  myLevel.SetLevelID(2);
  myLevel.Request("CoarseMap", myCMF.get());
  myCMF->SetParameter("Domain GID offsets", Teuchos::ParameterEntry(std::string("{100,50}")));
  myCMF->SetFactory("Aggregates", MueLu::NoFactory::getRCP());
  myCMF->SetFactory("Nullspace", MueLu::NoFactory::getRCP());
  myCMF->Build(myLevel);
  myCoarseMap = myLevel.Get<Teuchos::RCP<const Map> >("CoarseMap", myCMF.get());

  TEST_EQUALITY(myCoarseMap->getMinAllGlobalIndex() == 0, true);
  TEST_EQUALITY(myCoarseMap->getMaxLocalIndex() == 9, true);
}  // NonStandardCaseA

///////////////////////////////////////////////////////////////////////////

#define MUELU_ETI_GROUP(Scalar, LocalOrdinal, GlobalOrdinal, Node)                                                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoarseMapFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoarseMapFactory, StandardCase, Scalar, LocalOrdinal, GlobalOrdinal, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoarseMapFactory, NonStandardCaseA, Scalar, LocalOrdinal, GlobalOrdinal, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
