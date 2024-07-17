// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_IO.hpp>

#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_ClassicalMapFactory.hpp"
#include "MueLu_ClassicalPFactory.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(ClassicalPFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  out << "version: " << MueLu::Version() << std::endl;

  RCP<ClassicalPFactory> PFact = rcp(new ClassicalPFactory);
  TEST_EQUALITY(PFact != Teuchos::null, true);

}  // Constructor

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(ClassicalPFactory, BuildP_Direct, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real                  = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;

  Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
  coarseLevel.SetFactoryManager(Teuchos::null);

  GO nx         = 29;
  RCP<Matrix> A = test_factory::Build1DPoisson(nx);
  A->SetFixedBlockSize(1);
  fineLevel.Set("A", A);

  // This test only works in parallel if we have Zoltan2 & Tpetra
#ifndef HAVE_MUELU_ZOLTAN2
  if (A->getRowMap()->getComm()->getSize() > 1)
    return;
#else
  if ((A->getRowMap()->lib() == Xpetra::UseEpetra) &&
      (A->getRowMap()->getComm()->getSize() > 1))
    return;
#endif

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<real, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  LocalOrdinal NSdim         = 2;
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  nullSpace->randomize();
  fineLevel.Set("Nullspace", nullSpace);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  RCP<ClassicalMapFactory> cmFact = rcp(new ClassicalMapFactory());
  cmFact->SetFactory("Graph", dropFact);
  cmFact->SetFactory("UnAmalgamationInfo", amalgFact);

  Teuchos::ParameterList cp_params;
  cp_params.set("aggregation: classical scheme", "direct");
  RCP<ClassicalPFactory> PFact = rcp(new ClassicalPFactory());
  PFact->SetParameterList(cp_params);
  //    PFact->SetFactory("UnAmalgamationInfo", amalgFact);
  PFact->SetFactory("Graph", dropFact);
  PFact->SetFactory("DofsPerNode", dropFact);
  PFact->SetFactory("FC Splitting", cmFact);
  PFact->SetFactory("CoarseMap", cmFact);

  coarseLevel.Request("P", PFact.get());  // request Ptent
  coarseLevel.Request(*PFact);
  PFact->Build(fineLevel, coarseLevel);

  RCP<Matrix> P;
  coarseLevel.Get("P", P, PFact.get());

  // Check that the matrix is not zero.
  RCP<MultiVector> x = MultiVectorFactory::Build(P->getDomainMap(), 1);
  x->putScalar(TST::one());
  RCP<MultiVector> y = MultiVectorFactory::Build(P->getRangeMap(), 1);
  P->apply(*x, *y, Teuchos::NO_TRANS);
  const magnitude_type expected = TMT::one();
  TEST_FLOATING_EQUALITY(y->getVector(0)->normInf(), expected, 10 * TMT::eps());
}  // BuildP

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(ClassicalPFactory, BuildP_ClassicalModified, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real                  = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;

  Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
  coarseLevel.SetFactoryManager(Teuchos::null);

  GO nx         = 29;
  RCP<Matrix> A = test_factory::Build1DPoisson(nx);
  A->SetFixedBlockSize(1);
  fineLevel.Set("A", A);

  // This test only works in parallel if we have Zoltan2 & Tpetra
#ifndef HAVE_MUELU_ZOLTAN2
  if (A->getRowMap()->getComm()->getSize() > 1)
    return;
#else
  if ((A->getRowMap()->lib() == Xpetra::UseEpetra) &&
      (A->getRowMap()->getComm()->getSize() > 1))
    return;
#endif

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<real, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  LocalOrdinal NSdim         = 2;
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  nullSpace->randomize();
  fineLevel.Set("Nullspace", nullSpace);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  RCP<ClassicalMapFactory> cmFact = rcp(new ClassicalMapFactory());
  cmFact->SetFactory("Graph", dropFact);
  cmFact->SetFactory("UnAmalgamationInfo", amalgFact);

  Teuchos::ParameterList cp_params;
  cp_params.set("aggregation: classical scheme", "classical modified");
  RCP<ClassicalPFactory> PFact = rcp(new ClassicalPFactory());
  PFact->SetParameterList(cp_params);
  // PFact->SetFactory("UnAmalgamationInfo", amalgFact);
  PFact->SetFactory("Graph", dropFact);
  PFact->SetFactory("DofsPerNode", dropFact);
  PFact->SetFactory("FC Splitting", cmFact);
  PFact->SetFactory("CoarseMap", cmFact);

  coarseLevel.Request("P", PFact.get());  // request Ptent
  coarseLevel.Request(*PFact);
  PFact->Build(fineLevel, coarseLevel);

  RCP<Matrix> P;
  coarseLevel.Get("P", P, PFact.get());

  // Check that the matrix is not zero.
  RCP<MultiVector> x = MultiVectorFactory::Build(P->getDomainMap(), 1);
  x->putScalar(TST::one());
  RCP<MultiVector> y = MultiVectorFactory::Build(P->getRangeMap(), 1);
  P->apply(*x, *y, Teuchos::NO_TRANS);
  const magnitude_type expected = TMT::one();
  TEST_FLOATING_EQUALITY(y->getVector(0)->normInf(), expected, 10 * TMT::eps());

}  // BuildP

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(ClassicalPFactory, BuildP_Ext, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real                  = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;

  Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
  coarseLevel.SetFactoryManager(Teuchos::null);

  GO nx         = 29;
  RCP<Matrix> A = test_factory::Build1DPoisson(nx);
  A->SetFixedBlockSize(1);
  fineLevel.Set("A", A);

  // This test only works in parallel if we have Zoltan2 & Tpetra
#ifndef HAVE_MUELU_ZOLTAN2
  if (A->getRowMap()->getComm()->getSize() > 1)
    return;
#else
  if ((A->getRowMap()->lib() == Xpetra::UseEpetra) &&
      (A->getRowMap()->getComm()->getSize() > 1))
    return;
#endif

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<real, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  LocalOrdinal NSdim         = 2;
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  nullSpace->randomize();
  fineLevel.Set("Nullspace", nullSpace);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  RCP<ClassicalMapFactory> cmFact = rcp(new ClassicalMapFactory());
  cmFact->SetFactory("Graph", dropFact);
  cmFact->SetFactory("UnAmalgamationInfo", amalgFact);

  Teuchos::ParameterList cp_params;
  cp_params.set("aggregation: classical scheme", "ext+i");
  RCP<ClassicalPFactory> PFact = rcp(new ClassicalPFactory());
  PFact->SetParameterList(cp_params);
  //    PFact->SetFactory("UnAmalgamationInfo", amalgFact);
  PFact->SetFactory("Graph", dropFact);
  PFact->SetFactory("DofsPerNode", dropFact);
  PFact->SetFactory("FC Splitting", cmFact);
  PFact->SetFactory("CoarseMap", cmFact);

  coarseLevel.Request("P", PFact.get());  // request Ptent
  coarseLevel.Request(*PFact);
  PFact->Build(fineLevel, coarseLevel);

  RCP<Matrix> P;
  coarseLevel.Get("P", P, PFact.get());

  // Check that the matrix is not zero.
  RCP<MultiVector> x = MultiVectorFactory::Build(P->getDomainMap(), 1);
  x->putScalar(TST::one());
  RCP<MultiVector> y = MultiVectorFactory::Build(P->getRangeMap(), 1);
  P->apply(*x, *y, Teuchos::NO_TRANS);
  const magnitude_type expected = TMT::one();
  TEST_FLOATING_EQUALITY(y->getVector(0)->normInf(), expected, 10 * TMT::eps());

}  // BuildP

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node)                                                  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(ClassicalPFactory, Constructor, Scalar, LO, GO, Node)   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(ClassicalPFactory, BuildP_Direct, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(ClassicalPFactory, BuildP_ClassicalModified, Scalar, LO, GO, Node)
// Disabled until we actually have code to run these
#if 0

      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(ClassicalPFactory,BuildP_Ext,Scalar,LO,GO,Node)
#endif

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
