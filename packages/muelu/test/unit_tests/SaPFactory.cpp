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
#include <Xpetra_IO.hpp>

#include <MueLu_SaPFactory.hpp>
#include <MueLu_TrilinosSmoother.hpp>
#include <MueLu_UncoupledAggregationFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_TransPFactory.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_Utilities.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SaPFactory, Test0, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<SaPFactory> sapFactory = rcp(new SaPFactory);
  TEST_EQUALITY(sapFactory != Teuchos::null, true);

  out << *sapFactory << std::endl;
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SaPFactory, EpetraVsTpetra, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  out << "Skipping test because some required packages are not enabled (Tpetra, Epetra, EpetraExt, Ifpack, Ifpack2)." << std::endl;

}  // SaPFactory_EpetraVsTpetra

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SaPFactory, ConstrainRowOptimalScalarPDE, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"

  typedef Teuchos::ScalarTraits<SC> STS;
  SC zero = STS::zero(), one = STS::one();

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // Don't test for complex - matrix reader won't work
  if (STS::isComplex) {
    success = true;
    return;
  }

  RCP<Matrix> P;
  Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();
  P                         = Xpetra::IO<SC, LO, GO, NO>::Read("TestMatrices/SaP_constrainTest_P.mat", lib, comm);

  RCP<SaPFactory> sapFactory = rcp(new SaPFactory);
  Teuchos::ParameterList Pparams;
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
  for (size_t i = 0; i < (size_t)(P->getRowMap()->getLocalNumElements()); i++) {
    Teuchos::ArrayView<const LO> indices;
    Teuchos::ArrayView<const SC> vals;
    P->getLocalRowView((LO)i, indices, vals);
    size_t nnz = indices.size();
    for (size_t j = 0; j < nnz; j++) {
      if (Teuchos::ScalarTraits<SC>::real(vals[j]) < Teuchos::ScalarTraits<SC>::real(zero)) lowerViolation = true;
      if (Teuchos::ScalarTraits<SC>::real(vals[j]) > Teuchos::ScalarTraits<SC>::real(one)) upperViolation = true;
      // if (STS::magnitude(vals[j]) < STS::magnitude(zero)) lowerViolation = true;
      // if (STS::magnitude(vals[j]) > STS::magnitude(one))  upperViolation = true;
    }
  }
  TEST_EQUALITY(lowerViolation, false);
  TEST_EQUALITY(upperViolation, false);

}  // SaPFactory_ConstrainRowOptimalScalarPDE

#define MUELU_ETI_GROUP(SC, LO, GO, Node)                                            \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SaPFactory, Test0, SC, LO, GO, Node)          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SaPFactory, EpetraVsTpetra, SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SaPFactory, ConstrainRowOptimalScalarPDE, SC, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
