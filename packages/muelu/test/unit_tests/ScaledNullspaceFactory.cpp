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

#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_MatrixMatrix.hpp>

#include <MueLu_SaPFactory.hpp>
#include <MueLu_TrilinosSmoother.hpp>
#include <MueLu_ScaledNullspaceFactory.hpp>
#include <MueLu_NullspaceFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_TransPFactory.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_Utilities.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(ScaledNullspaceFactory, Test0, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_EPETRA_SCOPE_TPETRA_IS_DEFAULT(Scalar, GlobalOrdinal, Node);
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  Teuchos::Array<magnitude_type> results(2);
  Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();
  // Do not run this test for Epetra
  if (lib == Xpetra::UseEpetra) {
    TEST_EQUALITY(1, 1);
    return;
  }

  // Build a Matrix
  GlobalOrdinal nEle       = 20;
  const RCP<const Map> map = MapFactory::Build(lib, nEle, 0, comm);
  Teuchos::ParameterList matrixParameters;
  matrixParameters.set("nx", nEle);

  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, CrsMatrixWrap, MultiVector>("Laplace1D", map, matrixParameters);
  RCP<Matrix> Op = Pr->BuildMatrix();

  // build nullspace
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map, 1);
  nullSpace->putScalar((Scalar)1.0);
  Teuchos::Array<magnitude_type> norms(1);
  nullSpace->norm1(norms);
  if (comm->getRank() == 0)
    out << "||NS|| = " << norms[0] << std::endl;
  int maxLevels = 3;

  // fill hierarchy
  RCP<Hierarchy> H = rcp(new Hierarchy());
  H->setDefaultVerbLevel(Teuchos::VERB_HIGH);

  RCP<Level> Finest = H->GetLevel();
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  Finest->Set("A", Op);                 // set fine level matrix
  Finest->Set("Nullspace", nullSpace);  // set null space information for finest level

  // define transfer operators
  RCP<NullspaceFactory> NSFact        = rcp(new NullspaceFactory());
  RCP<ScaledNullspaceFactory> SNSFact = rcp(new ScaledNullspaceFactory());

  RCP<TentativePFactory> Ptent1 = rcp(new TentativePFactory());
  NSFact->SetFactory("Nullspace", Ptent1);
  RCP<TentativePFactory> Ptent2 = rcp(new TentativePFactory());
  Teuchos::ParameterList pt2_list;
  pt2_list.set("Nullspace name", "Scaled Nullspace");
  Ptent2->SetParameterList(pt2_list);

  RCP<Factory> Rfact = rcp(new TransPFactory());
  Rfact->SetFactory("P", Ptent2);
  RCP<RAPFactory> Acfact = rcp(new RAPFactory());
  H->SetMaxCoarseSize(1);

  // setup smoothers
  Teuchos::ParameterList smootherParamList;
  smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
  smootherParamList.set("relaxation: sweeps", (LocalOrdinal)1);
  smootherParamList.set("relaxation: damping factor", (Scalar)1.0);
  RCP<SmootherPrototype> smooProto = rcp(new TrilinosSmoother("RELAXATION", smootherParamList));
  RCP<SmootherFactory> SmooFact    = rcp(new SmootherFactory(smooProto));
  Acfact->setVerbLevel(Teuchos::VERB_HIGH);

  RCP<SmootherFactory> coarseSolveFact = rcp(new SmootherFactory(smooProto, Teuchos::null));

  FactoryManager M;
  M.SetKokkosRefactor(false);
  M.SetFactory("P", Ptent1);
  M.SetFactory("R", Rfact);
  M.SetFactory("Nullspace", NSFact);
  M.SetFactory("Scaled Nullspace", SNSFact);
  M.SetFactory("Ptent", Ptent1);
  M.SetFactory("Smoother", SmooFact);
  M.SetFactory("CoarseSolver", coarseSolveFact);

  H->Setup(M, 0, maxLevels);

  // test some basic multigrid data
  RCP<Level> coarseLevel = H->GetLevel(1);
  TEST_EQUALITY(coarseLevel->IsRequested("A", MueLu::NoFactory::get()), false);
  TEST_EQUALITY(coarseLevel->IsRequested("P", MueLu::NoFactory::get()), false);
  TEST_EQUALITY(coarseLevel->IsRequested("PreSmoother", MueLu::NoFactory::get()), false);
  TEST_EQUALITY(coarseLevel->IsRequested("PostSmoother", MueLu::NoFactory::get()), false);
  TEST_EQUALITY(coarseLevel->IsRequested("R", MueLu::NoFactory::get()), false);
  TEST_EQUALITY(coarseLevel->IsAvailable("A", MueLu::NoFactory::get()), true);
  TEST_EQUALITY(coarseLevel->IsAvailable("P", MueLu::NoFactory::get()), true);
  TEST_EQUALITY(coarseLevel->IsAvailable("PreSmoother", MueLu::NoFactory::get()), true);
  //    TEST_EQUALITY(coarseLevel->IsAvailable("PostSmoother",MueLu::NoFactory::get()), true);
  TEST_EQUALITY(coarseLevel->IsAvailable("R", MueLu::NoFactory::get()), true);
  TEST_EQUALITY(coarseLevel->GetKeepFlag("A", MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(coarseLevel->GetKeepFlag("P", MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(coarseLevel->GetKeepFlag("PreSmoother", MueLu::NoFactory::get()), MueLu::Final);
  //    TEST_EQUALITY(coarseLevel->GetKeepFlag("PostSmoother",MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(coarseLevel->GetKeepFlag("R", MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(coarseLevel->IsRequested("P", Ptent1.get()), false);
  TEST_EQUALITY(coarseLevel->IsRequested("PreSmoother", SmooFact.get()), false);
  TEST_EQUALITY(coarseLevel->IsRequested("PostSmoother", SmooFact.get()), false);
  TEST_EQUALITY(coarseLevel->IsRequested("R", Rfact.get()), false);
  TEST_EQUALITY(coarseLevel->IsRequested("A", Acfact.get()), false);
  TEST_EQUALITY(coarseLevel->IsAvailable("P", Ptent1.get()), false);
  TEST_EQUALITY(coarseLevel->IsAvailable("PreSmoother", SmooFact.get()), false);
  //    TEST_EQUALITY(coarseLevel->IsAvailable("PostSmoother",SmooFact.get()), false);
  TEST_EQUALITY(coarseLevel->IsAvailable("R", Rfact.get()), false);
  TEST_EQUALITY(coarseLevel->IsAvailable("A", Acfact.get()), false);
  TEST_EQUALITY(coarseLevel->GetKeepFlag("P", Ptent1.get()), 0);
  TEST_EQUALITY(coarseLevel->GetKeepFlag("PreSmoother", SmooFact.get()), 0);
  //    TEST_EQUALITY(coarseLevel->GetKeepFlag("PostSmoother",SmooFact.get()), 0);
  TEST_EQUALITY(coarseLevel->GetKeepFlag("R", Rfact.get()), 0);
  TEST_EQUALITY(coarseLevel->GetKeepFlag("A", Acfact.get()), 0);
}

#define MUELU_ETI_GROUP(SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(ScaledNullspaceFactory, Test0, SC, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
