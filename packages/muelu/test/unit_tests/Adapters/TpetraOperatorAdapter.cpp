// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_TpetraMultiVector.hpp>

#include <MueLu_FactoryManagerBase.hpp>
#include <MueLu_Hierarchy.hpp>
#include <MueLu_PFactory.hpp>
#include <MueLu_SaPFactory.hpp>
#include <MueLu_TransPFactory.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_TrilinosSmoother.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>

namespace MueLuTests {

namespace {  // helpers local to this translation unit, shared by the CreateTpetraPreconditioner
             // const-overload tests added together with MueLu::CreateTpetraPreconditioner's
             // RCP<const Tpetra::Operator> / RCP<CrsMatrix> overloads.

// Small 1D-Poisson based fixture that exposes the handles every const-overload test needs:
// the Xpetra matrix, its Tpetra::CrsMatrix view, and non-const / const Tpetra::Operator RCPs.
// Also provides a few one-liner utilities (pair creation, random RHS, residual reduction)
// so that each test is driven by a couple of lines plus its specific assertion tolerances.
template <class SC, class LO, class GO, class NO>
struct ConstOverloadFixture {
  using magnitude_type             = typename Teuchos::ScalarTraits<SC>::magnitudeType;
  using tpetra_crsmatrix_type      = Tpetra::CrsMatrix<SC, LO, GO, NO>;
  using tpetra_operator_type       = Tpetra::Operator<SC, LO, GO, NO>;
  using matrix_type                = Xpetra::Matrix<SC, LO, GO, NO>;
  using multivector_type           = Xpetra::MultiVector<SC, LO, GO, NO>;
  using multivector_factory_type   = Xpetra::MultiVectorFactory<SC, LO, GO, NO>;
  using muelu_tpetra_operator_type = MueLu::TpetraOperator<SC, LO, GO, NO>;
  using precond_pair               = std::pair<Teuchos::RCP<muelu_tpetra_operator_type>,
                                 Teuchos::RCP<muelu_tpetra_operator_type>>;

  Teuchos::RCP<matrix_type> Op;
  Teuchos::RCP<tpetra_crsmatrix_type> tpA;
  Teuchos::RCP<tpetra_operator_type> opNonConst;
  Teuchos::RCP<const tpetra_operator_type> opConst;

  explicit ConstOverloadFixture(GO rowsPerRank) {
    auto comm  = TestHelpers::Parameters::getDefaultComm();
    Op         = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(rowsPerRank * comm->getSize());
    tpA        = Xpetra::toTpetra(Op);
    opNonConst = Teuchos::rcp_implicit_cast<tpetra_operator_type>(tpA);
    opConst    = Teuchos::rcp_implicit_cast<const tpetra_operator_type>(tpA);
  }

  // CreateTpetraPreconditioner takes ParameterList by non-const reference and fills defaults
  // in place.  Both API paths therefore need their own fresh named lvalue copy of `base` so
  // they start from identical configuration (temporary ParameterList cannot bind to ParameterList&).
  precond_pair makeNonConstVsConstOperator(const Teuchos::ParameterList& base) const {
    Teuchos::ParameterList listA(base), listB(base);
    auto precA = MueLu::CreateTpetraPreconditioner<SC, LO, GO, NO>(opNonConst, listA);
    auto precB = MueLu::CreateTpetraPreconditioner<SC, LO, GO, NO>(opConst, listB);
    return std::make_pair(precA, precB);
  }

  // Compare the dedicated RCP<CrsMatrix> overload against the generic RCP<Operator> overload.
  precond_pair makeCrsMatrixVsOperator(const Teuchos::ParameterList& base) const {
    Teuchos::ParameterList listCrs(base), listOp(base);
    auto precCrs = MueLu::CreateTpetraPreconditioner<SC, LO, GO, NO>(tpA, listCrs);
    auto precOp  = MueLu::CreateTpetraPreconditioner<SC, LO, GO, NO>(opNonConst, listOp);
    return std::make_pair(precCrs, precOp);
  }

  // Same fine matrix, reuse entry point tested with non-const vs const CrsMatrix handle.
  void reuseOnConstNonConstCrsMatrix(muelu_tpetra_operator_type& precForNonConst,
                                     muelu_tpetra_operator_type& precForConst) const {
    MueLu::ReuseTpetraPreconditioner(tpA, precForNonConst);
    Teuchos::RCP<const tpetra_crsmatrix_type> tpAconst =
        Teuchos::rcp_implicit_cast<const tpetra_crsmatrix_type>(tpA);
    MueLu::ReuseTpetraPreconditioner(tpAconst, precForConst);
  }

  // Normalized random RHS (per-column unit 2-norm).
  Teuchos::RCP<multivector_type> makeRandomRHS(int numVec, unsigned int seed) const {
    auto RHS = multivector_factory_type::Build(Op->getRowMap(), numVec);
    RHS->setSeed(seed);
    RHS->randomize();
    Teuchos::Array<magnitude_type> colNorms(numVec);
    RHS->norm2(colNorms());
    for (int j = 0; j < numVec; ++j)
      RHS->getVectorNonConst(j)->scale(Teuchos::ScalarTraits<magnitude_type>::one() / colNorms[j]);
    return RHS;
  }

  Teuchos::RCP<multivector_type> makeZero(int numVec) const {
    auto X = multivector_factory_type::Build(Op->getRowMap(), numVec);
    X->putScalar(Teuchos::ScalarTraits<SC>::zero());
    return X;
  }

  // Per-column relative residual ||A*X - RHS||_2 / ||RHS||_2.
  Teuchos::Array<magnitude_type> residualReduction(const multivector_type& X,
                                                   const multivector_type& RHS) const {
    const int numVec = static_cast<int>(X.getNumVectors());
    auto r           = multivector_factory_type::Build(Op->getRowMap(), numVec);
    Op->apply(X, *r);
    r->update(Teuchos::ScalarTraits<SC>::one(), RHS, -Teuchos::ScalarTraits<SC>::one());
    Teuchos::Array<magnitude_type> rNorm(numVec), rhsNorm(numVec);
    r->norm2(rNorm());
    RHS.norm2(rhsNorm());
    Teuchos::Array<magnitude_type> red(numVec);
    for (int j = 0; j < numVec; ++j)
      red[j] = rNorm[j] / rhsNorm[j];
    return red;
  }
};

// Apply precA and precB to `RHS` and verify API equivalence:
//  (a) each preconditioner reduces the residual below `redMax`;
//  (b) the two API paths agree per-column within `gapMax`;
//  (c) both hierarchies have the same number of levels.
// Two CreateTpetraPreconditioner calls build independent hierarchies, so bit-identical
// behavior is NOT expected: UncoupledAggregation uses random graph coloring, so aggregate
// counts on coarse levels can differ by a handful (e.g. 721 vs 722 on level 2), and the
// SaPFactory eigenvalue estimate (PowerMethod) is nondeterministic across independent
// setups (slightly different Prolongator damping factor).  Tolerances reflect that:
// `redMax` is a "the preconditioner is not broken" bound, `gapMax` is a per-column
// disagreement bound that absorbs aggregation-level randomness, not numerical noise.
template <class SC, class LO, class GO, class NO>
void expectEquivalentPreconditioner(Teuchos::FancyOStream& out,
                                    bool& success,
                                    const ConstOverloadFixture<SC, LO, GO, NO>& fx,
                                    MueLu::TpetraOperator<SC, LO, GO, NO>& precA,
                                    MueLu::TpetraOperator<SC, LO, GO, NO>& precB,
                                    const Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>>& RHS,
                                    typename Teuchos::ScalarTraits<SC>::magnitudeType redMax,
                                    typename Teuchos::ScalarTraits<SC>::magnitudeType gapMax,
                                    const std::string& label) {
  using magnitude_type = typename Teuchos::ScalarTraits<SC>::magnitudeType;
  const int numVec     = static_cast<int>(RHS->getNumVectors());

  auto Xa = fx.makeZero(numVec);
  auto Xb = fx.makeZero(numVec);
  precA.apply(*Xpetra::toTpetra(RHS), *Xpetra::toTpetra(Xa));
  precB.apply(*Xpetra::toTpetra(RHS), *Xpetra::toTpetra(Xb));

  const auto redA = fx.residualReduction(*Xa, *RHS);
  const auto redB = fx.residualReduction(*Xb, *RHS);

  const magnitude_type zero = Teuchos::ScalarTraits<magnitude_type>::zero();
  magnitude_type maxRedA = zero, maxRedB = zero, maxGap = zero;
  for (int j = 0; j < numVec; ++j) {
    if (redA[j] > maxRedA) maxRedA = redA[j];
    if (redB[j] > maxRedB) maxRedB = redB[j];
    const magnitude_type gap = Teuchos::ScalarTraits<magnitude_type>::magnitude(redA[j] - redB[j]);
    if (gap > maxGap) maxGap = gap;
  }

  out << label << ": max residual reduction A=" << std::setiosflags(std::ios::fixed)
      << std::setprecision(10) << maxRedA << ", B=" << maxRedB
      << ", max per-column gap=" << maxGap << std::endl;

  TEUCHOS_TEST_EQUALITY(maxRedA < redMax, true, out, success);
  TEUCHOS_TEST_EQUALITY(maxRedB < redMax, true, out, success);
  TEUCHOS_TEST_EQUALITY(maxGap < gapMax, true, out, success);
  TEUCHOS_TEST_EQUALITY(precA.GetHierarchy()->GetNumLevels(),
                        precB.GetHierarchy()->GetNumLevels(), out, success);
}

template <class SC, class LO, class GO, class NO>
void expectEquivalentHierarchyLevels(Teuchos::FancyOStream& out,
                                     bool& success,
                                     MueLu::TpetraOperator<SC, LO, GO, NO>& precA,
                                     MueLu::TpetraOperator<SC, LO, GO, NO>& precB,
                                     const std::string& label) {
  const int levelsA = precA.GetHierarchy()->GetNumLevels();
  const int levelsB = precB.GetHierarchy()->GetNumLevels();
  out << label << ": numLevels A=" << levelsA << ", B=" << levelsB << std::endl;
  TEUCHOS_TEST_EQUALITY(levelsA, levelsB, out, success);
}

}  // namespace

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TpetraOperator, Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  typedef MueLu::TpetraOperator<SC, LO, GO, NO> muelu_tpetra_operator_type;
  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType magnitude_type;

  if (TestHelpers::Parameters::getLib() == Xpetra::UseTpetra) {
    // matrix
    RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();
    RCP<Matrix> Op                     = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(6561 * comm->getSize());  //=8*3^6
    RCP<const Map> map                 = Op->getRowMap();

    RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map, 1);
    nullSpace->putScalar((SC)1.0);
    // Teuchos::Array<Teuchos::ScalarTraits<SC>::magnitudeType> norms(1);
    Teuchos::Array<magnitude_type> norms(1);
    nullSpace->norm1(norms);

    RCP<Hierarchy> H = rcp(new Hierarchy());
    H->setDefaultVerbLevel(Teuchos::VERB_NONE);

    RCP<MueLu::Level> Finest = H->GetLevel();
    Finest->setDefaultVerbLevel(Teuchos::VERB_NONE);
    Finest->Set("A", Op);
    H->Setup();

    // ------------- test Tpetra Operator wrapping MueLu hierarchy ------------
    RCP<muelu_tpetra_operator_type> tH = rcp(new muelu_tpetra_operator_type(H));

    RCP<MultiVector> RHS1 = MultiVectorFactory::Build(Op->getRowMap(), 1);
    RCP<MultiVector> X1   = MultiVectorFactory::Build(Op->getRowMap(), 1);

    // normalized RHS, zero initial guess
    RHS1->setSeed(846930886);
    RHS1->randomize();
    RHS1->norm2(norms);
    RHS1->scale(1 / norms[0]);

    X1->putScalar((SC)0.0);

    tH->apply(*(Xpetra::toTpetra(RHS1)), *(Xpetra::toTpetra(X1)));

    X1->norm2(norms);
    out << "after apply, ||X1|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    // -------------- test MueLu Hierarchy directly -----------------------
    RCP<MultiVector> RHS2 = MultiVectorFactory::Build(Op->getRowMap(), 1);
    RCP<MultiVector> X2   = MultiVectorFactory::Build(Op->getRowMap(), 1);

    // normalized RHS, zero initial guess
    RHS2->setSeed(846930886);
    RHS2->randomize();
    RHS2->norm2(norms);
    RHS2->scale(1 / norms[0]);

    X2->putScalar((SC)0.0);

    int iterations = 1;
    H->Iterate(*RHS2, *X2, iterations);

    X2->norm2(norms);
    out << "after apply, ||X2|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    RCP<MultiVector> diff = MultiVectorFactory::Build(Op->getRowMap(), 1);
    diff->putScalar(0.0);

    diff->update(1.0, *X1, -1.0, *X2, 0.0);
    diff->norm2(norms);
    TEST_EQUALITY(norms[0] < 1e-10, true);

  } else {
    out << "This test is enabled only for linAlgebra=Tpetra." << std::endl;
  }
}  // Apply

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TpetraOperator, Getters, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  using TpetraOperatorType = MueLu::TpetraOperator<SC, LO, GO, NO>;

  if (TestHelpers::Parameters::getLib() == Xpetra::UseTpetra) {
    const auto Op = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(100);

    ////////////////////////////////////////
    //////////   WITH HIERARCHY   //////////
    ////////////////////////////////////////
    {
      const auto map = Op->getRowMap();
      auto H         = rcp(new Hierarchy());
      H->setDefaultVerbLevel(Teuchos::VERB_NONE);

      auto Finest = H->GetLevel();
      Finest->setDefaultVerbLevel(Teuchos::VERB_NONE);
      Finest->Set("A", Op);

      H->Setup();
      auto tH = rcp(new TpetraOperatorType(H));

      TEST_EQUALITY(tH->GetOperator(), Teuchos::null);
      TEST_INEQUALITY(tH->GetHierarchy(), Teuchos::null);

      TEST_INEQUALITY(tH->getRangeMap(), Teuchos::null);
      TEST_INEQUALITY(tH->getDomainMap(), Teuchos::null);

      // Hardcoded false
      TEST_EQUALITY(tH->hasTransposeApply(), false);
    }

    ///////////////////////////////////////
    //////////   WITH OPERATOR   //////////
    ///////////////////////////////////////
    {
      auto tO = rcp(new TpetraOperatorType((Teuchos::RCP<Xpetra::Operator<SC, LO, GO, NO>>)(Op)));

      TEST_INEQUALITY(tO->GetOperator(), Teuchos::null);
      TEST_EQUALITY(tO->GetHierarchy(), Teuchos::null);

      TEST_INEQUALITY(tO->getRangeMap(), Teuchos::null);
      TEST_INEQUALITY(tO->getDomainMap(), Teuchos::null);

      // Hardcoded false
      TEST_EQUALITY(tO->hasTransposeApply(), false);
    }
  }
}

// Reproducer for Trilinos issue #15062: callers that only have Teuchos::RCP<const Tpetra::Operator>
// could not use CreateTpetraPreconditioner without const_cast.  The overload taking
// RCP<const Operator> must produce a preconditioner equivalent to the non-const overload.
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TpetraOperator, CreatePreconditioner_ConstOperator, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

#if defined(HAVE_MUELU_IFPACK2) && defined(HAVE_MUELU_AMESOS2)
  using magnitude_type = typename Teuchos::ScalarTraits<SC>::magnitudeType;

  if (TestHelpers::Parameters::getLib() == Xpetra::UseTpetra) {
    ConstOverloadFixture<SC, LO, GO, NO> fx(/*rowsPerRank=*/6561);
    Teuchos::ParameterList baseList;
    baseList.set("aggregation: deterministic", true);
    auto precs = fx.makeNonConstVsConstOperator(baseList);
    auto RHS   = fx.makeRandomRHS(/*numVec=*/1, /*seed=*/846930886u);

    expectEquivalentPreconditioner(out, success, fx, *precs.first, *precs.second, RHS,
                                   static_cast<magnitude_type>(0.95),
                                   1000 * Teuchos::ScalarTraits<magnitude_type>::eps(),
                                   "CreatePreconditioner_ConstOperator");
  } else {
    out << "This test is enabled only for linAlgebra=Tpetra." << std::endl;
  }
#else
  out << "Skipping test because some required packages are not enabled (Tpetra, Ifpack2, Amesos2)." << std::endl;
#endif
}  // CreatePreconditioner_ConstOperator

// CreateTpetraPreconditioner(inA) with a default-constructed (empty) ParameterList must
// likewise be equivalent across the non-const / const Operator overloads.
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TpetraOperator, CreatePreconditioner_ConstOperator_NoParameterList, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

#if defined(HAVE_MUELU_IFPACK2) && defined(HAVE_MUELU_AMESOS2)
  using magnitude_type = typename Teuchos::ScalarTraits<SC>::magnitudeType;

  if (TestHelpers::Parameters::getLib() == Xpetra::UseTpetra) {
    ConstOverloadFixture<SC, LO, GO, NO> fx(/*rowsPerRank=*/6561);
    Teuchos::ParameterList baseList;
    baseList.set("aggregation: deterministic", true);
    auto precs = fx.makeNonConstVsConstOperator(baseList);
    auto RHS   = fx.makeRandomRHS(/*numVec=*/1, /*seed=*/846930886u);

    expectEquivalentPreconditioner(out, success, fx, *precs.first, *precs.second, RHS,
                                   static_cast<magnitude_type>(0.95),
                                   1000 * Teuchos::ScalarTraits<magnitude_type>::eps(),
                                   "CreatePreconditioner_ConstOperator_NoParameterList");
  } else {
    out << "This test is enabled only for linAlgebra=Tpetra." << std::endl;
  }
#else
  out << "Skipping test because some required packages are not enabled (Tpetra, Ifpack2, Amesos2)." << std::endl;
#endif
}  // CreatePreconditioner_ConstOperator_NoParameterList

// Smaller 1D problem (fewer rows per rank).  Residual reduction is looser on a shallow grid.
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TpetraOperator, CreatePreconditioner_ConstOperator_SmallGrid, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

#if defined(HAVE_MUELU_IFPACK2) && defined(HAVE_MUELU_AMESOS2)
  using magnitude_type = typename Teuchos::ScalarTraits<SC>::magnitudeType;

  if (TestHelpers::Parameters::getLib() == Xpetra::UseTpetra) {
    ConstOverloadFixture<SC, LO, GO, NO> fx(/*rowsPerRank=*/243);
    Teuchos::ParameterList baseList;
    baseList.set("aggregation: deterministic", true);
    auto precs = fx.makeNonConstVsConstOperator(baseList);
    auto RHS   = fx.makeRandomRHS(/*numVec=*/1, /*seed=*/123456789u);

    expectEquivalentPreconditioner(out, success, fx, *precs.first, *precs.second, RHS,
                                   static_cast<magnitude_type>(0.95),
                                   1000 * Teuchos::ScalarTraits<magnitude_type>::eps(),
                                   "CreatePreconditioner_ConstOperator_SmallGrid");
  } else {
    out << "This test is enabled only for linAlgebra=Tpetra." << std::endl;
  }
#else
  out << "Skipping test because some required packages are not enabled (Tpetra, Ifpack2, Amesos2)." << std::endl;
#endif
}  // CreatePreconditioner_ConstOperator_SmallGrid

// Hierarchy depth from CreateTpetraPreconditioner must not depend on const vs non-const Operator handle.
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TpetraOperator, CreatePreconditioner_ConstOperator_HierarchyLevels, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

#if defined(HAVE_MUELU_IFPACK2) && defined(HAVE_MUELU_AMESOS2)
  if (TestHelpers::Parameters::getLib() == Xpetra::UseTpetra) {
    ConstOverloadFixture<SC, LO, GO, NO> fx(/*rowsPerRank=*/2187);
    Teuchos::ParameterList baseList;
    baseList.set("aggregation: deterministic", true);
    auto precs = fx.makeNonConstVsConstOperator(baseList);

    expectEquivalentHierarchyLevels(out, success, *precs.first, *precs.second,
                                    "CreatePreconditioner_ConstOperator_HierarchyLevels");
  } else {
    out << "This test is enabled only for linAlgebra=Tpetra." << std::endl;
  }
#else
  out << "Skipping test because some required packages are not enabled (Tpetra, Ifpack2, Amesos2)." << std::endl;
#endif
}  // CreatePreconditioner_ConstOperator_HierarchyLevels

// Two-column multivector: const vs non-const CreateTpetraPreconditioner apply must agree column-wise.
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TpetraOperator, CreatePreconditioner_ConstOperator_TwoVectors, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

#if defined(HAVE_MUELU_IFPACK2) && defined(HAVE_MUELU_AMESOS2)
  using magnitude_type = typename Teuchos::ScalarTraits<SC>::magnitudeType;

  if (TestHelpers::Parameters::getLib() == Xpetra::UseTpetra) {
    ConstOverloadFixture<SC, LO, GO, NO> fx(/*rowsPerRank=*/729);
    Teuchos::ParameterList baseList;
    baseList.set("aggregation: deterministic", true);
    auto precs = fx.makeNonConstVsConstOperator(baseList);
    auto RHS   = fx.makeRandomRHS(/*numVec=*/2, /*seed=*/97531u);

    expectEquivalentPreconditioner(out, success, fx, *precs.first, *precs.second, RHS,
                                   static_cast<magnitude_type>(0.95),
                                   1000 * Teuchos::ScalarTraits<magnitude_type>::eps(),
                                   "CreatePreconditioner_ConstOperator_TwoVectors");
  } else {
    out << "This test is enabled only for linAlgebra=Tpetra." << std::endl;
  }
#else
  out << "Skipping test because some required packages are not enabled (Tpetra, Ifpack2, Amesos2)." << std::endl;
#endif
}  // CreatePreconditioner_ConstOperator_TwoVectors

// ReuseTpetraPreconditioner with RCP<const CrsMatrix> vs RCP<CrsMatrix> on the same fine matrix.
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TpetraOperator, ReuseTpetraPreconditioner_ConstCrsMatrix, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

#if defined(HAVE_MUELU_IFPACK2) && defined(HAVE_MUELU_AMESOS2)
  using magnitude_type = typename Teuchos::ScalarTraits<SC>::magnitudeType;

  if (TestHelpers::Parameters::getLib() == Xpetra::UseTpetra) {
    ConstOverloadFixture<SC, LO, GO, NO> fx(/*rowsPerRank=*/2187);
    Teuchos::ParameterList baseList;
    baseList.set("aggregation: deterministic", true);
    auto precs = fx.makeNonConstVsConstOperator(baseList);
    fx.reuseOnConstNonConstCrsMatrix(*precs.first, *precs.second);

    auto RHS = fx.makeRandomRHS(/*numVec=*/1, /*seed=*/314159265u);
    expectEquivalentPreconditioner(out, success, fx, *precs.first, *precs.second, RHS,
                                   static_cast<magnitude_type>(0.95),
                                   1000 * Teuchos::ScalarTraits<magnitude_type>::eps(),
                                   "ReuseTpetraPreconditioner_ConstCrsMatrix");
  } else {
    out << "This test is enabled only for linAlgebra=Tpetra." << std::endl;
  }
#else
  out << "Skipping test because some required packages are not enabled (Tpetra, Ifpack2, Amesos2)." << std::endl;
#endif
}  // ReuseTpetraPreconditioner_ConstCrsMatrix

// Regression: Teuchos::RCP<Tpetra::CrsMatrix> used to be ambiguous between the overloads taking
// RCP<Tpetra::Operator> and RCP<const Tpetra::Operator> (e.g. Zoltan2 Sphynx).  The dedicated
// RCP<CrsMatrix> overload must compile and agree with the RCP<Operator> entry point.
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TpetraOperator, CreatePreconditioner_RcpCrsMatrixOverload, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

#if defined(HAVE_MUELU_IFPACK2) && defined(HAVE_MUELU_AMESOS2)
  using magnitude_type = typename Teuchos::ScalarTraits<SC>::magnitudeType;

  if (TestHelpers::Parameters::getLib() == Xpetra::UseTpetra) {
    ConstOverloadFixture<SC, LO, GO, NO> fx(/*rowsPerRank=*/6561);
    Teuchos::ParameterList baseList;
    baseList.set("aggregation: deterministic", true);
    auto precs = fx.makeCrsMatrixVsOperator(baseList);
    auto RHS   = fx.makeRandomRHS(/*numVec=*/1, /*seed=*/846930886u);

    expectEquivalentPreconditioner(out, success, fx, *precs.first, *precs.second, RHS,
                                   static_cast<magnitude_type>(0.95),
                                   1000 * Teuchos::ScalarTraits<magnitude_type>::eps(),
                                   "CreatePreconditioner_RcpCrsMatrixOverload");
  } else {
    out << "This test is enabled only for linAlgebra=Tpetra." << std::endl;
  }
#else
  out << "Skipping test because some required packages are not enabled (Tpetra, Ifpack2, Amesos2)." << std::endl;
#endif
}  // CreatePreconditioner_RcpCrsMatrixOverload

#define MUELU_ETI_GROUP(Scalar, LocalOrdinal, GlobalOrdinal, Node)                                                                                    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TpetraOperator, Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node)                                              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TpetraOperator, Getters, Scalar, LocalOrdinal, GlobalOrdinal, Node)                                            \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TpetraOperator, CreatePreconditioner_ConstOperator, Scalar, LocalOrdinal, GlobalOrdinal, Node)                 \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TpetraOperator, CreatePreconditioner_ConstOperator_NoParameterList, Scalar, LocalOrdinal, GlobalOrdinal, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TpetraOperator, CreatePreconditioner_ConstOperator_SmallGrid, Scalar, LocalOrdinal, GlobalOrdinal, Node)       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TpetraOperator, CreatePreconditioner_ConstOperator_HierarchyLevels, Scalar, LocalOrdinal, GlobalOrdinal, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TpetraOperator, CreatePreconditioner_ConstOperator_TwoVectors, Scalar, LocalOrdinal, GlobalOrdinal, Node)      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TpetraOperator, ReuseTpetraPreconditioner_ConstCrsMatrix, Scalar, LocalOrdinal, GlobalOrdinal, Node)           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TpetraOperator, CreatePreconditioner_RcpCrsMatrixOverload, Scalar, LocalOrdinal, GlobalOrdinal, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
