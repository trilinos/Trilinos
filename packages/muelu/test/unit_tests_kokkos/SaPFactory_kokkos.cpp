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
#include <Xpetra_IO.hpp>

#include "MueLu_TestHelpers_kokkos.hpp"
#include "MueLu_IteratorOps.hpp"
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
  RCP<Matrix> Ptest           = MueLu::IteratorOps<SC, LO, GO, NO>::Jacobi(omega, *invDiag, *A, *Ptent, Teuchos::null, out, "label", APparams);

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
#endif

#define MUELU_ETI_GROUP(SC, LO, GO, NO)                                                       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SaPFactory_kokkos, Constructor, SC, LO, GO, NO)        \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SaPFactory_kokkos, Build, SC, LO, GO, NO)              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SaPFactory_kokkos, EnforceConstraints, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SaPFactory_kokkos, ConstrainRowOptimalScalarPDE, SC, LO, GO, NO)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
