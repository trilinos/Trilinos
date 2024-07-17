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

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_IteratorOps.hpp>
#include <Xpetra_IO.hpp>

#include "MueLu_TestHelpers_kokkos.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_SaPFactory.hpp"
#include "MueLu_Utilities.hpp"

#include "MueLu_UseDefaultTypes.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SaPFactory_kokkos, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<SaPFactory> sapFactory = rcp(new SaPFactory);
  TEST_EQUALITY(sapFactory != Teuchos::null, true);

  out << *sapFactory << std::endl;
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SaPFactory_kokkos, Build, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();

  // construct two levels
  Level fineLevel, coarseLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

  // construct matrices
  const typename Teuchos::ScalarTraits<Scalar>::magnitudeType lambdaMax = 5;
  RCP<Matrix> A                                                         = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build2DPoisson(27 * comm->getSize());
  A->SetMaxEigenvalueEstimate(lambdaMax);
  RCP<Matrix> Ptent = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build2DPoisson(27 * comm->getSize());

  // set level matrices
  fineLevel.Set("A", A);
  coarseLevel.Set("P", Ptent);

  // construct the factory to be tested
  const double dampingFactor = 0.5;
  RCP<SaPFactory> sapFactory = rcp(new SaPFactory);
  ParameterList Pparams;
  Pparams.set("sa: damping factor", dampingFactor);
  Pparams.set("use kokkos refactor", true);
  sapFactory->SetParameterList(Pparams);
  sapFactory->SetFactory("A", MueLu::NoFactory::getRCP());
  sapFactory->SetFactory("P", MueLu::NoFactory::getRCP());

  // build the data
  coarseLevel.Request("P", sapFactory.get());
  sapFactory->Build(fineLevel, coarseLevel);

  // fetch the data
  RCP<Matrix> Pfact = coarseLevel.Get<RCP<Matrix>>("P", sapFactory.get());

  // construct the data to compare
  SC omega                    = dampingFactor / lambdaMax;
  RCP<Vector> invDiag         = Utilities::GetMatrixDiagonalInverse(*A);
  RCP<ParameterList> APparams = rcp(new ParameterList);
  RCP<Matrix> Ptest           = Xpetra::IteratorOps<SC, LO, GO, NO>::Jacobi(omega, *invDiag, *A, *Ptent, Teuchos::null, out, "label", APparams);

  // compare matrices by multiplying them by a random vector
  RCP<MultiVector> X = MultiVectorFactory::Build(A->getDomainMap(), 1);
  X->setSeed(846930886);
  X->randomize();

  RCP<MultiVector> Bfact = MultiVectorFactory::Build(A->getRangeMap(), 1);
  RCP<MultiVector> Btest = MultiVectorFactory::Build(A->getRangeMap(), 1);

  typedef Teuchos::ScalarTraits<SC> STS;
  SC zero = STS::zero(), one = STS::one();

  Pfact->apply(*X, *Bfact, Teuchos::NO_TRANS, one, zero);
  Ptest->apply(*X, *Btest, Teuchos::NO_TRANS, one, zero);
  Btest->update(-one, *Bfact, one);

  Array<typename STS::magnitudeType> norms(1);
  Btest->norm2(norms);
  out << "|| B_factory - B_test || = " << norms[0] << std::endl;
  TEST_EQUALITY(norms[0] < 1e-12, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SaPFactory_kokkos, EnforceConstraints, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  typedef Teuchos::ScalarTraits<SC> STS;
  SC zero = STS::zero(), one = STS::one();

  RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();

  // construct two levels
  Level fineLevel, coarseLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

  // construct matrices
  SC lambdaMax = 2;

  Xpetra::UnderlyingLib lib = TestHelpers_kokkos::Parameters::getLib();
  Teuchos::ParameterList mp;
  mp.set("matrixType", "Star2D");
  const GO nx = 6;
  mp.set("nx", nx);
  mp.set("ny", nx);
  mp.set("a", 2.0);
  mp.set("b", -1.0);
  mp.set("c", -1.0);
  mp.set("d", 0.5);
  mp.set("e", 0.5);
  mp.set("z1", -0.25);
  mp.set("z2", -0.25);
  mp.set("z3", -0.25);
  mp.set("z4", -0.25);
  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildMatrix(mp, lib);
  A->SetMaxEigenvalueEstimate(lambdaMax);
  mp.set("a", 1.0);
  mp.set("b", 1.0);
  mp.set("c", 1.0);
  mp.set("d", 1.0);
  mp.set("e", 1.0);
  mp.set("z1", 1.0);
  mp.set("z2", 1.0);
  mp.set("z3", 1.0);
  mp.set("z4", 1.0);
  RCP<Matrix> Ptent = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildMatrix(mp, lib);

  // set level matrices
  fineLevel.Set("A", A);
  coarseLevel.Set("P", Ptent);

  // construct the factory to be tested
  const double dampingFactor = 0.5;
  RCP<SaPFactory> sapFactory = rcp(new SaPFactory);
  ParameterList Pparams;
  Pparams.set("sa: damping factor", dampingFactor);
  Pparams.set("sa: enforce constraints", true);
  Pparams.set("sa: max eigenvalue", (double)2.0);
  Pparams.set("tentative: calculate qr", false);
  Pparams.set("use kokkos refactor", true);

  sapFactory->SetParameterList(Pparams);
  sapFactory->SetFactory("A", MueLu::NoFactory::getRCP());
  sapFactory->SetFactory("P", MueLu::NoFactory::getRCP());

  // build the data
  coarseLevel.Request("P", sapFactory.get());
  sapFactory->Build(fineLevel, coarseLevel);

  // fetch the data
  RCP<Matrix> Pfact = coarseLevel.Get<RCP<Matrix>>("P", sapFactory.get());

  // check that row sums are all one by checking the norm of the vector
  RCP<MultiVector> X     = MultiVectorFactory::Build(A->getDomainMap(), 1);
  RCP<MultiVector> Bfact = MultiVectorFactory::Build(A->getRangeMap(), 1);
  X->putScalar(one);
  Pfact->apply(*X, *Bfact, Teuchos::NO_TRANS, one, zero);
  Array<typename STS::magnitudeType> norms(1);
  Bfact->norm2(norms);
  out << "|| B_factory ones || = " << norms[0] << std::endl;
  using magnitude_type = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  using TMT            = Teuchos::ScalarTraits<magnitude_type>;
  TEST_FLOATING_EQUALITY(STS::magnitude(norms[0]), STS::magnitude(as<Scalar>(6.0)), 100 * TMT::eps());
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SaPFactory_kokkos, ConstrainRowOptimalScalarPDE, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  typedef Teuchos::ScalarTraits<SC> STS;
  SC zero = STS::zero(), one = STS::one();

  RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();

  // Don't test for complex - matrix reader won't work
  if (STS::isComplex) {
    success = true;
    return;
  }

  RCP<Matrix> P;
  Xpetra::UnderlyingLib lib = TestHelpers_kokkos::Parameters::getLib();
  P                         = Xpetra::IO<SC, LO, GO, NO>::Read("../unit_tests/TestMatrices/SaP_constrainTest_P.mat", lib, comm);

  RCP<SaPFactory> sapFactory = rcp(new SaPFactory);
  ParameterList Pparams;
  Pparams.set("use kokkos refactor", true);
  sapFactory->SetParameterList(Pparams);
  sapFactory->optimalSatisfyPConstraintsForScalarPDEs(P);

  // check that row sums are all one by checking the norm of the vector
  // note: optimalSatisfyPConstraintsForScalarPDEs preserves row sum of original P (one in this case),
  //       but SatisfyPConstraints normalizes each row sum to one.
  RCP<MultiVector> X     = MultiVectorFactory::Build(P->getDomainMap(), 1);
  RCP<MultiVector> Bfact = MultiVectorFactory::Build(P->getRangeMap(), 1);
  X->putScalar(one);
  P->apply(*X, *Bfact, Teuchos::NO_TRANS, one, zero);
  Array<typename STS::magnitudeType> norms(1);
  Bfact->norm2(norms);
  out << "|| B_factory ones || = " << norms[0] << std::endl;
  using magnitude_type = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  using TMT            = Teuchos::ScalarTraits<magnitude_type>;
  TEST_FLOATING_EQUALITY(STS::magnitude(norms[0]), STS::magnitude(6.0), 100 * TMT::eps());

  // check that the min and max of each row are in [0,1]
  bool lowerViolation = false;
  bool upperViolation = false;
  for (size_t j = 0; j < as<size_t>(P->getRowMap()->getLocalNumElements()); j++) {
    Teuchos::ArrayView<const LocalOrdinal> indices;
    Teuchos::ArrayView<const Scalar> vals;
    P->getLocalRowView((LocalOrdinal)j, indices, vals);
    size_t nnz = indices.size();
    for (LO i = 0; i < (LO)nnz; i++) {
      if (Teuchos::ScalarTraits<SC>::real(vals[i]) < Teuchos::ScalarTraits<SC>::real(zero)) lowerViolation = true;
      if (Teuchos::ScalarTraits<SC>::real(vals[i]) > Teuchos::ScalarTraits<SC>::real(one)) upperViolation = true;
    }
  }
  TEST_EQUALITY(lowerViolation, false);
  TEST_EQUALITY(upperViolation, false);

}  // SaPFactory_ConstrainRowOptimalScalarPDE

// FIXME_KOKKOS: uncomment the test when we get all corresponding factories ported to kokkos
#if 0
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT) && defined(HAVE_MUELU_IFPACK) && defined(HAVE_MUELU_IFPACK2)
  TEUCHOS_UNIT_TEST(SaPFactory_kokkos, EpetraVsTpetra)
  {
#include "MueLu_UseShortNames.hpp"
    MueLu::VerboseObject::SetMueLuOStream(Teuchos::rcpFromRef(out));

    out << "version: " << MueLu::Version() << std::endl;
    out << "Compare results of Epetra and Tpetra" << std::endl;
    out << "for 3 level AMG solver using smoothed aggregation with" << std::endl;
    out << "one SGS sweep on each multigrid level as pre- and postsmoother" << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    typedef Teuchos::ScalarTraits<SC> STS;
    SC zero = STS::zero(), one = STS::one();

    Array<STS::magnitudeType> results(2);

    // run test only on 1 proc
    if(comm->getSize() == 1) {
      Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;

      // run Epetra and Tpetra test
      for (int run = 0; run < 2; run++) { //TODO: create a subfunction instead or Tuple of UnderlyingLib
        if (run == 0) lib = Xpetra::UseEpetra;
        else          lib = Xpetra::UseTpetra;

        // generate problem
        LO maxLevels = 3;
        LO its       = 10;
        GO nEle      = 63;
        const RCP<const Map> map = MapFactory::Build(lib, nEle, 0, comm);
        Teuchos::ParameterList matrixParameters;
        matrixParameters.set("nx", nEle);

        RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
          Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>("Laplace1D", map, matrixParameters);
        RCP<Matrix> Op = Pr->BuildMatrix();

        // build nullspace
        RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
        nullSpace->putScalar(one);
        Array<STS::magnitudeType> norms(1);
        nullSpace->norm1(norms);
        if (comm->getRank() == 0)
          out << "||NS|| = " << norms[0] << std::endl;

        // fill hierarchy
        RCP<Hierarchy> H = rcp( new Hierarchy() );
        H->setDefaultVerbLevel(Teuchos::VERB_HIGH);

        RCP<Level> Finest = H->GetLevel();
        Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
        Finest->Set("A",Op);                      // set fine level matrix
        Finest->Set("Nullspace",nullSpace);       // set null space information for finest level

        // define transfer operators
        RCP<UncoupledAggregationFactory> UncoupledAggFact = rcp(new UncoupledAggregationFactory());
        UncoupledAggFact->SetMinNodesPerAggregate(3);
        UncoupledAggFact->SetMaxNeighAlreadySelected(0);
        UncoupledAggFact->SetOrdering("natural");

        RCP<TentativePFactory> Ptentfact = rcp(new TentativePFactory());
        RCP<SaPFactory>        Pfact = rcp( new SaPFactory());
        RCP<Factory>      Rfact = rcp( new TransPFactory() );
        RCP<RAPFactory>        Acfact = rcp( new RAPFactory() );
        H->SetMaxCoarseSize(1);

        // setup smoothers
        Teuchos::ParameterList smootherParamList;
        smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
        smootherParamList.set("relaxation: sweeps", (LO) 1);
        smootherParamList.set("relaxation: damping factor", (SC) 1.0);
        RCP<SmootherPrototype> smooProto = rcp( new TrilinosSmoother("RELAXATION", smootherParamList) );
        RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smooProto) );
        Acfact->setVerbLevel(Teuchos::VERB_HIGH);

        RCP<SmootherFactory> coarseSolveFact = rcp(new SmootherFactory(smooProto, Teuchos::null));

        FactoryManager M;
        M.SetFactory("P", Pfact);
        M.SetFactory("R", Rfact);
        M.SetFactory("A", Acfact);
        M.SetFactory("Ptent", Ptentfact);
        M.SetFactory("Aggregates", UncoupledAggFact);
        M.SetFactory("Smoother", SmooFact);
        M.SetFactory("CoarseSolver", coarseSolveFact);

        H->Setup(M, 0, maxLevels);

        // test some basic multigrid data
        RCP<Level> coarseLevel = H->GetLevel(1);
        TEST_EQUALITY(coarseLevel->IsRequested("A",MueLu::NoFactory::get()), false);
        TEST_EQUALITY(coarseLevel->IsRequested("P",MueLu::NoFactory::get()), false);
        TEST_EQUALITY(coarseLevel->IsRequested("PreSmoother",MueLu::NoFactory::get()), false);
        TEST_EQUALITY(coarseLevel->IsRequested("PostSmoother",MueLu::NoFactory::get()), false);
        TEST_EQUALITY(coarseLevel->IsRequested("R",MueLu::NoFactory::get()), false);
        TEST_EQUALITY(coarseLevel->IsAvailable("A",MueLu::NoFactory::get()), true);
        TEST_EQUALITY(coarseLevel->IsAvailable("P",MueLu::NoFactory::get()), true);
        TEST_EQUALITY(coarseLevel->IsAvailable("PreSmoother",MueLu::NoFactory::get()), true);
        TEST_EQUALITY(coarseLevel->IsAvailable("PostSmoother",MueLu::NoFactory::get()), true);
        TEST_EQUALITY(coarseLevel->IsAvailable("R",MueLu::NoFactory::get()), true);
        TEST_EQUALITY(coarseLevel->GetKeepFlag("A",MueLu::NoFactory::get()), MueLu::Final);
        TEST_EQUALITY(coarseLevel->GetKeepFlag("P",MueLu::NoFactory::get()), MueLu::Final);
        TEST_EQUALITY(coarseLevel->GetKeepFlag("PreSmoother",MueLu::NoFactory::get()), MueLu::Final);
        TEST_EQUALITY(coarseLevel->GetKeepFlag("PostSmoother",MueLu::NoFactory::get()), MueLu::Final);
        TEST_EQUALITY(coarseLevel->GetKeepFlag("R",MueLu::NoFactory::get()), MueLu::Final);
        TEST_EQUALITY(coarseLevel->IsRequested("P",Pfact.get()), false);
        TEST_EQUALITY(coarseLevel->IsRequested("P",Ptentfact.get()), false);
        TEST_EQUALITY(coarseLevel->IsRequested("PreSmoother",SmooFact.get()), false);
        TEST_EQUALITY(coarseLevel->IsRequested("PostSmoother",SmooFact.get()), false);
        TEST_EQUALITY(coarseLevel->IsRequested("R",Rfact.get()), false);
        TEST_EQUALITY(coarseLevel->IsRequested("A",Acfact.get()), false);
        TEST_EQUALITY(coarseLevel->IsAvailable("P",Pfact.get()), false);
        TEST_EQUALITY(coarseLevel->IsAvailable("P",Ptentfact.get()), false);
        TEST_EQUALITY(coarseLevel->IsAvailable("PreSmoother",SmooFact.get()), false);
        TEST_EQUALITY(coarseLevel->IsAvailable("PostSmoother",SmooFact.get()), false);
        TEST_EQUALITY(coarseLevel->IsAvailable("R",Rfact.get()), false);
        TEST_EQUALITY(coarseLevel->IsAvailable("A",Acfact.get()), false);
        TEST_EQUALITY(coarseLevel->GetKeepFlag("P",Pfact.get()), 0);
        TEST_EQUALITY(coarseLevel->GetKeepFlag("P",Ptentfact.get()), 0);
        TEST_EQUALITY(coarseLevel->GetKeepFlag("PreSmoother",SmooFact.get()), 0);
        TEST_EQUALITY(coarseLevel->GetKeepFlag("PostSmoother",SmooFact.get()), 0);
        TEST_EQUALITY(coarseLevel->GetKeepFlag("R",Rfact.get()), 0);
        TEST_EQUALITY(coarseLevel->GetKeepFlag("A",Acfact.get()), 0);
        RCP<Matrix> P1 = coarseLevel->Get< RCP<Matrix> >("P");
        RCP<Matrix> R1 = coarseLevel->Get< RCP<Matrix> >("R");
        TEST_EQUALITY(P1->getGlobalNumRows(), 63);
        TEST_EQUALITY(P1->getGlobalNumCols(), 21);
        TEST_EQUALITY(R1->getGlobalNumRows(), 21);
        TEST_EQUALITY(R1->getGlobalNumCols(), 63);
        RCP<Level> coarseLevel2 = H->GetLevel(2);
        TEST_EQUALITY(coarseLevel2->IsRequested("A",MueLu::NoFactory::get()), false);
        TEST_EQUALITY(coarseLevel2->IsRequested("P",MueLu::NoFactory::get()), false);
        TEST_EQUALITY(coarseLevel2->IsRequested("R",MueLu::NoFactory::get()), false);
        TEST_EQUALITY(coarseLevel2->IsRequested("PreSmoother",MueLu::NoFactory::get()), false);
        TEST_EQUALITY(coarseLevel2->IsRequested("PostSmoother",MueLu::NoFactory::get()), false);
        TEST_EQUALITY(coarseLevel2->IsAvailable("A",MueLu::NoFactory::get()), true);
        TEST_EQUALITY(coarseLevel2->IsAvailable("P",MueLu::NoFactory::get()), true);
        TEST_EQUALITY(coarseLevel2->IsAvailable("PreSmoother",MueLu::NoFactory::get()), true);
        TEST_EQUALITY(coarseLevel2->IsAvailable("PostSmoother",MueLu::NoFactory::get()), false);
        TEST_EQUALITY(coarseLevel2->IsAvailable("R",MueLu::NoFactory::get()), true);
        TEST_EQUALITY(coarseLevel2->GetKeepFlag("A",MueLu::NoFactory::get()), MueLu::Final);
        TEST_EQUALITY(coarseLevel2->GetKeepFlag("P",MueLu::NoFactory::get()), MueLu::Final);
        TEST_EQUALITY(coarseLevel2->GetKeepFlag("PreSmoother",MueLu::NoFactory::get()), MueLu::Final);
        TEST_EQUALITY(coarseLevel2->GetKeepFlag("PostSmoother",MueLu::NoFactory::get()), 0);
        TEST_EQUALITY(coarseLevel2->GetKeepFlag("R",MueLu::NoFactory::get()), MueLu::Final);
        TEST_EQUALITY(coarseLevel2->IsRequested("P",Pfact.get()), false);
        TEST_EQUALITY(coarseLevel2->IsRequested("P",Ptentfact.get()), false);
        TEST_EQUALITY(coarseLevel2->IsRequested("R",Rfact.get()), false);
        TEST_EQUALITY(coarseLevel2->IsAvailable("P",Pfact.get()), false);
        TEST_EQUALITY(coarseLevel2->IsAvailable("P",Ptentfact.get()), false);
        TEST_EQUALITY(coarseLevel2->IsAvailable("PreSmoother",SmooFact.get()), false);
        TEST_EQUALITY(coarseLevel2->IsAvailable("PostSmoother",SmooFact.get()), false);
        TEST_EQUALITY(coarseLevel2->IsAvailable("R",Rfact.get()), false);
        TEST_EQUALITY(coarseLevel2->GetKeepFlag("P",Pfact.get()), 0);
        TEST_EQUALITY(coarseLevel2->GetKeepFlag("P",Ptentfact.get()), 0);
        TEST_EQUALITY(coarseLevel2->GetKeepFlag("PreSmoother",SmooFact.get()), 0);
        TEST_EQUALITY(coarseLevel2->GetKeepFlag("PostSmoother",SmooFact.get()), 0);
        TEST_EQUALITY(coarseLevel2->GetKeepFlag("R",Rfact.get()), 0);
        RCP<Matrix> P2 = coarseLevel2->Get< RCP<Matrix> >("P");
        RCP<Matrix> R2 = coarseLevel2->Get< RCP<Matrix> >("R");
        TEST_EQUALITY(P2->getGlobalNumRows(), 21);
        TEST_EQUALITY(P2->getGlobalNumCols(), 7);
        TEST_EQUALITY(R2->getGlobalNumRows(), 7);
        TEST_EQUALITY(R2->getGlobalNumCols(), 21);

        Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO> > PtentTPtent = Xpetra::MatrixMatrix<Scalar,LO,GO>::Multiply(*P1,true,*P1,false,out);
        TEST_EQUALITY(PtentTPtent->getGlobalMaxNumRowEntries()-3<1e-12, true);
        TEST_EQUALITY(P1->getGlobalMaxNumRowEntries()-2<1e-12, true);
        TEST_EQUALITY(P2->getGlobalMaxNumRowEntries()-2<1e-12, true);

        // Define RHS
        RCP<MultiVector> X = MultiVectorFactory::Build(map,1);
        RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);

        X->putScalar(1.0);
        X->norm2(norms);
        if (comm->getRank() == 0)
          out << "||X_true|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

        Op->apply(*X,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

        // Use AMG directly as an iterative method
        {
          X->putScalar( (SC) 0.0);

          H->Iterate(*RHS,*X,its);

          X->norm2(norms);
          if (comm->getRank() == 0)
            out << "||X_" << std::setprecision(2) << its << "|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;
          results[run] = norms[0];
        }
      }

      TEST_EQUALITY(results[0] - results[1] < 1e-10, true); // check results of EPETRA vs TPETRA
    } // comm->getSize == 1

  } //SaPFactory_EpetraVsTpetra
#endif
#endif

#define MUELU_ETI_GROUP(SC, LO, GO, NO)                                                       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SaPFactory_kokkos, Constructor, SC, LO, GO, NO)        \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SaPFactory_kokkos, Build, SC, LO, GO, NO)              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SaPFactory_kokkos, EnforceConstraints, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SaPFactory_kokkos, ConstrainRowOptimalScalarPDE, SC, LO, GO, NO)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
