// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>

#include <MueLu_config.hpp>

#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <Xpetra_MultiVectorFactory.hpp>

#include <MueLu_Utilities.hpp>
#include <MueLu_RebalanceTransferFactory.hpp>
#include <MueLu_SaPFactory.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_TransPFactory.hpp>
#include <MueLu_NullspaceFactory.hpp>
#include <MueLu_CoalesceDropFactory.hpp>
#include <MueLu_UncoupledAggregationFactory.hpp>
#include <MueLu_FactoryManager.hpp>
#include <MueLu_ZoltanInterface.hpp>
#include <MueLu_RepartitionFactory.hpp>
#include <MueLu_MultiVectorTransferFactory.hpp>
#include <Galeri_XpetraUtils.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RebalanceTransfer, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<RebalanceTransferFactory> ptFactory = rcp(new RebalanceTransferFactory());
  TEST_EQUALITY(ptFactory != Teuchos::null, true);
}  // Constructor test

#ifdef NEVER_TESTED_TODO

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RebalanceTransfer, Build1, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  Level fineLevel, coarseLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);
  GO nx         = 199;
  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(nx);
  fineLevel.Set("A", A);

  // build coordinates
  Teuchos::ParameterList list;
  list.set("nx", nx);
  RCP<MultiVector> coordVector = Galeri::Xpetra::CreateCartesianCoordinates<SC, LO, GO, Map, MultiVector>("1D", A->getRowMap(), list);
  fineLevel.Set("Coordinates", coordVector);

  RCP<UncoupledAggregationFactory> UncoupledAggFact = rcp(new UncoupledAggregationFactory());
  RCP<TentativePFactory> Ptentfact                  = rcp(new TentativePFactory(UncoupledAggFact));
  RCP<SaPFactory> Pfact                             = rcp(new SaPFactory(Ptentfact));
  RCP<Factory> Rfact                                = rcp(new TransPFactory(Pfact));
  RCP<RAPFactory> Acfact                            = rcp(new RAPFactory(Pfact, Rfact));
  RCP<Factory> Rtentfact                            = rcp(new TransPFactory(Ptentfact));

  RCP<MultiVectorTransferFactory> mvTransFact = rcp(new MultiVectorTransferFactory("Coordinates", "R", Rtentfact));
  Acfact->AddTransferFactory(mvTransFact);
  RCP<ZoltanInterface> zoltan             = rcp(new ZoltanInterface(Acfact, mvTransFact));
  RCP<RepartitionFactory> RepartitionFact = rcp(new RepartitionFactory(zoltan, Acfact));

  coarseLevel.Request("A", Acfact.get());                     // kick off the DeclareInputs
  coarseLevel.Request("Permutation", RepartitionFact.get());  // request permutation matrix
  // coarseLevel.Request("P",Pfact.get());
  coarseLevel.Request("R", Rtentfact.get());
  coarseLevel.Request("Coordinates", mvTransFact.get());

  RCP<RebalanceTransferFactory> ptFactory = rcp(new RebalanceTransferFactory(RepartitionFact, Acfact, Pfact, MueLu::INTERPOLATION));
  ptFactory->SetParameter("type", ParameterEntry("Interpolation"));
  coarseLevel.Request("P", ptFactory.get());
  ptFactory->Build(fineLevel, coarseLevel);

  ptFactory = rcp(new RebalanceTransferFactory(RepartitionFact, Acfact, Rfact, MueLu::RESTRICTION));
  coarseLevel.Request("R", ptFactory.get());
  ptFactory->Build(fineLevel, coarseLevel);

}  // Constructor test

#endif

#define MUELU_ETI_GROUP(SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RebalanceTransfer, Constructor, SC, LO, GO, NO)
#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
