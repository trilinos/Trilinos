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
#include "MueLu_TestHelpersSmoothers.hpp"

#include <MueLu_Ifpack2Smoother.hpp>
#include <MueLu_Utilities.hpp>

/*
   Comments about tests with hard coded results:
   1) Chebyshev smoothing must pass for any number of processors.
   2) Gauss-Seidel must pass for 1 and 4 processors.
   3) For any processor count except 1 and 4, the Gauss-Seidel test will
   report "passing", but this is only because the Teuchos test macro is skipped.
   */

namespace MueLuTests {

// this namespace already has:  #include "MueLu_UseShortNames.hpp"
using namespace TestHelpers::Smoothers;

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Ifpack2Smoother, NotSetup, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    Ifpack2Smoother smoother("RELAXATION", Teuchos::ParameterList());
    testApplyNoSetup(smoother, out, success);
  }
}

// Tests interface to Ifpack2's Jacobi preconditioner.
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Ifpack2Smoother, HardCodedResult_Jacobi, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    if (Teuchos::ScalarTraits<Scalar>::isComplex) {
      out << "Skipping Tpetra for SC type \"complex\"" << std::endl;
      return;
    }
    Teuchos::ParameterList paramList;
    paramList.set("relaxation: type", "Jacobi");
    paramList.set("relaxation: sweeps", (int)1);
    paramList.set("relaxation: damping factor", (double)0.8);
    paramList.set("relaxation: zero starting solution", false);

    Ifpack2Smoother smoother("RELAXATION", paramList);

    typename Teuchos::ScalarTraits<SC>::magnitudeType residualNorms = testApply_A125_X1_RHS0(smoother, out, success);

    RCP<const Teuchos::Comm<int> > comm                                  = TestHelpers::Parameters::getDefaultComm();
    const typename Teuchos::ScalarTraits<SC>::magnitudeType expectedNorm = 6.3245553203367577133e-01;
    switch (comm->getSize()) {
      case 1:
      case 4:
        TEST_FLOATING_EQUALITY(residualNorms, expectedNorm, 1e-12);
        break;
      default:
        out << "Pass/Fail is checked only for 1 and 4 processes." << std::endl;
        break;
    }  // switch

    auto ifpack2prec = smoother.getPreconditioner();

    // reuse setup & solve
    residualNorms = testApply_A125_X1_RHS0(smoother, out, success);
    switch (comm->getSize()) {
      case 1:
      case 4:
        TEST_FLOATING_EQUALITY(residualNorms, expectedNorm, 1e-12);
        break;
      default:
        out << "Pass/Fail is checked only for 1 and 4 processes." << std::endl;
        break;
    }  // switch

    // make sure the Ifpack2 preconditioner did not get replaced
    TEUCHOS_ASSERT_EQUALITY(ifpack2prec.ptr(), smoother.getPreconditioner().ptr());
  }
}  // Jacobi

// Tests interface to Ifpack2's Gauss-Seidel preconditioner.
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Ifpack2Smoother, HardCodedResult_GaussSeidel, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    if (Teuchos::ScalarTraits<Scalar>::isComplex) {
      out << "Skipping Tpetra for SC type \"complex\"" << std::endl;
      return;
    }
    Teuchos::ParameterList paramList;
    paramList.set("relaxation: type", "Gauss-Seidel");
    paramList.set("relaxation: sweeps", (int)1);
    paramList.set("relaxation: damping factor", (double)1.0);
    paramList.set("relaxation: zero starting solution", false);

    Ifpack2Smoother smoother("RELAXATION", paramList);

    typename Teuchos::ScalarTraits<SC>::magnitudeType residualNorms = testApply_A125_X1_RHS0(smoother, out, success);

    RCP<const Teuchos::Comm<int> > comm                                  = TestHelpers::Parameters::getDefaultComm();
    const typename Teuchos::ScalarTraits<SC>::magnitudeType expectedNorm = 5.773502691896257e-01;
    switch (comm->getSize()) {
      case 1:
      case 4:
        TEST_FLOATING_EQUALITY(residualNorms, expectedNorm, 1e-12);
        break;
      default:
        out << "Pass/Fail is checked only for 1 and 4 processes." << std::endl;
        break;
    }  // switch

    auto ifpack2prec = smoother.getPreconditioner();

    // reuse setup & solve
    residualNorms = testApply_A125_X1_RHS0(smoother, out, success);
    switch (comm->getSize()) {
      case 1:
      case 4:
        TEST_FLOATING_EQUALITY(residualNorms, expectedNorm, 1e-12);
        break;
      default:
        out << "Pass/Fail is checked only for 1 and 4 processes." << std::endl;
        break;
    }  // switch

    // make sure the Ifpack2 preconditioner did not get replaced
    TEUCHOS_ASSERT_EQUALITY(ifpack2prec.ptr(), smoother.getPreconditioner().ptr());
  }
}  // GaussSeidel

// Tests interface to Ifpack2's Gauss-Seidel preconditioner.
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Ifpack2Smoother, HardCodedResult_GaussSeidel2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    if (Teuchos::ScalarTraits<Scalar>::isComplex) {
      out << "Skipping Tpetra for SC type \"complex\"" << std::endl;
      return;
    }
    Teuchos::ParameterList paramList;
    paramList.set("relaxation: type", "Gauss-Seidel");
    paramList.set("relaxation: sweeps", (int)10);
    paramList.set("relaxation: damping factor", (double)1.0);
    paramList.set("relaxation: zero starting solution", false);

    Ifpack2Smoother smoother("RELAXATION", paramList);

    typename Teuchos::ScalarTraits<SC>::magnitudeType residualNorms = testApply_A125_X1_RHS0(smoother, out, success);

    RCP<const Teuchos::Comm<int> > comm                                   = TestHelpers::Parameters::getDefaultComm();
    const typename Teuchos::ScalarTraits<SC>::magnitudeType expectedNorm1 = 8.326553652741774e-02;
    const typename Teuchos::ScalarTraits<SC>::magnitudeType expectedNorm4 = 8.326553653078517e-02;
    switch (comm->getSize()) {
      case 1:
        TEST_FLOATING_EQUALITY(residualNorms, expectedNorm1, 1e-12);
        break;
      case 4:
        TEST_FLOATING_EQUALITY(residualNorms, expectedNorm4, 1e-12);
        break;
      default:
        out << "Pass/Fail is checked only for 1 and 4 processes." << std::endl;
        break;
    }  // switch

    auto ifpack2prec = smoother.getPreconditioner();

    // reuse setup & solve
    residualNorms = testApply_A125_X1_RHS0(smoother, out, success);
    switch (comm->getSize()) {
      case 1:
        TEST_FLOATING_EQUALITY(residualNorms, expectedNorm1, 1e-12);
        break;
      case 4:
        TEST_FLOATING_EQUALITY(residualNorms, expectedNorm4, 1e-12);
        break;
      default:
        out << "Pass/Fail is checked only for 1 and 4 processes." << std::endl;
        break;
    }  // switch

    // make sure the Ifpack2 preconditioner did not get replaced
    TEUCHOS_ASSERT_EQUALITY(ifpack2prec.ptr(), smoother.getPreconditioner().ptr());
  }
}  // GaussSeidel

// Tests interface to Ifpack2's Chebyshev preconditioner
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Ifpack2Smoother, HardCodedResult_Chebyshev, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    using magnitude_type = typename Teuchos::ScalarTraits<SC>::magnitudeType;
    using TMT            = Teuchos::ScalarTraits<magnitude_type>;

    if (Teuchos::ScalarTraits<Scalar>::isComplex) {
      out << "Skipping Tpetra for SC type \"complex\"" << std::endl;
      return;
    }

    Teuchos::ParameterList paramList;
    paramList.set("chebyshev: degree", (int)3);
    paramList.set("chebyshev: max eigenvalue", (double)1.98476);
    paramList.set("chebyshev: min eigenvalue", (double)1.0);
    paramList.set("chebyshev: ratio eigenvalue", (double)20);
    paramList.set("chebyshev: zero starting solution", false);
    Ifpack2Smoother smoother("CHEBYSHEV", paramList);

    typename Teuchos::ScalarTraits<SC>::magnitudeType residualNorms      = testApply_A125_X1_RHS0(smoother, out, success);
    const typename Teuchos::ScalarTraits<SC>::magnitudeType expectedNorm = 5.269156e-01;
    TEST_FLOATING_EQUALITY(residualNorms, expectedNorm, (1e-7 < TMT::eps() ? 10 * TMT::eps() : 1e-7));  // Compare to residual reported by ML

    auto ifpack2prec = smoother.getPreconditioner();

    // reuse setup & solve
    residualNorms = testApply_A125_X1_RHS0(smoother, out, success);
    TEST_FLOATING_EQUALITY(residualNorms, expectedNorm, (1e-7 < TMT::eps() ? 10 * TMT::eps() : 1e-7));  // Compare to residual reported by ML

    // make sure the Ifpack2 preconditioner did not get replaced
    TEUCHOS_ASSERT_EQUALITY(ifpack2prec.ptr(), smoother.getPreconditioner().ptr());
  }
}  // Chebyshev

// Tests interface to Ifpack2's ILU(0) preconditioner.
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Ifpack2Smoother, HardCodedResult_ILU, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    // FIXME this will probably fail in parallel b/c it becomes block Jacobi
    using magnitude_type = typename Teuchos::ScalarTraits<SC>::magnitudeType;

    Teuchos::ParameterList paramList;
    Ifpack2Smoother smoother("ILUT", paramList);

    magnitude_type residualNorms = testApply_A125_X0_RandomRHS(smoother, out, success);

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    if (comm->getSize() == 1) {
      TEST_EQUALITY(residualNorms < 100 * Teuchos::ScalarTraits<magnitude_type>::eps(), true);
    } else {
      out << "Pass/Fail is only checked in serial." << std::endl;
    }

    auto ifpack2prec = smoother.getPreconditioner();

    // reuse setup & solve
    residualNorms = testApply_A125_X1_RHS0(smoother, out, success);
    if (comm->getSize() == 1) {
      TEST_EQUALITY(residualNorms < 100 * Teuchos::ScalarTraits<magnitude_type>::eps(), true);
    } else {
      out << "Pass/Fail is only checked in serial." << std::endl;
    }

    // make sure the Ifpack2 preconditioner did not get replaced
    TEUCHOS_ASSERT_EQUALITY(ifpack2prec.ptr(), smoother.getPreconditioner().ptr());
  }
}  // ILU

// Tests two sweeps of ILUT in Ifpack2
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Ifpack2Smoother, ILU_TwoSweeps, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType magnitude_type;
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    // FIXME this will probably fail in parallel b/c it becomes block Jacobi

    Teuchos::ParameterList paramList;
    Ifpack2Smoother smoother("ILUT", paramList);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Teuchos::RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(125);
    Level level;
    TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
    level.Set("A", A);
    smoother.Setup(level);

    RCP<MultiVector> X   = MultiVectorFactory::Build(A->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(A->getRangeMap(), 1);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    A->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    smoother.Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    out << "solve with random initial guess" << std::endl;
    X->randomize();
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    smoother.Apply(*X, *RHS, false);  // nonzero initial guess

    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm2 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm2[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    if (comm->getSize() == 1) {
      // TEST_EQUALITY(residualNorms < 1e-10, true);
      TEST_EQUALITY(residualNorm1[0] != residualNorm2[0], true);
    } else {
      out << "Pass/Fail is only checked in serial." << std::endl;
    }
  }
}  // ILU

// Make sure Ifpack2's Banded relaxation actually gets called
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Ifpack2Smoother, BandedRelaxation, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType magnitude_type;
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    Teuchos::ParameterList paramList;
    paramList.set("partitioner: PDE equations", 5);  // Warning: This number has to be compatible with the problem size
    Ifpack2Smoother smoother("BANDED RELAXATION", paramList);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Teuchos::RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(125);
    Level level;
    TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
    level.Set("A", A);
    smoother.Setup(level);

    RCP<MultiVector> X   = MultiVectorFactory::Build(A->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(A->getRangeMap(), 1);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    A->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    smoother.Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    out << "solve with random initial guess" << std::endl;
    X->randomize();
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    smoother.Apply(*X, *RHS, false);  // nonzero initial guess

    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm2 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm2[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    if (comm->getSize() == 1) {
      // TEST_EQUALITY(residualNorms < 1e-10, true);
      TEST_EQUALITY(residualNorm1[0] != residualNorm2[0], true);
    } else {
      out << "Pass/Fail is only checked in serial." << std::endl;
    }
  }
}  // banded

// Make sure Ifpack2's TriDi relaxation actually gets called
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Ifpack2Smoother, TriDiRelaxation, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType magnitude_type;
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    Teuchos::ParameterList paramList;
    Ifpack2Smoother smoother("TRIDI RELAXATION", paramList);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Teuchos::RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(125);
    Level level;
    TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
    level.Set("A", A);
    smoother.Setup(level);

    RCP<MultiVector> X   = MultiVectorFactory::Build(A->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(A->getRangeMap(), 1);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    A->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    smoother.Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    out << "solve with random initial guess" << std::endl;
    X->randomize();
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    smoother.Apply(*X, *RHS, false);  // nonzero initial guess

    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm2 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm2[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    if (comm->getSize() == 1) {
      // TEST_EQUALITY(residualNorms < 1e-10, true);
      TEST_EQUALITY(residualNorm1[0] != residualNorm2[0], true);
    } else {
      out << "Pass/Fail is only checked in serial." << std::endl;
    }
  }
}  // tridi

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Ifpack2Smoother, BlockRelaxation_Autosize, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType magnitude_type;
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    Teuchos::ParameterList matrixParams, ifpack2Params;

    matrixParams.set("matrixType", "Laplace1D");
    matrixParams.set("nx", (GlobalOrdinal)20);  // needs to be even
    Teuchos::RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildMatrix(matrixParams, Xpetra::UseTpetra);
    A->SetFixedBlockSize(2);

    ifpack2Params.set("partitioner: type", "linear");
    Ifpack2Smoother smoother("BLOCK RELAXATION", ifpack2Params);

    Level level;
    TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
    level.Set("A", A);
    smoother.Setup(level);

    TEST_EQUALITY(1, 1);
  }
}  // banded

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Ifpack2Smoother, BlockCrsMatrix_Relaxation_ViaPoint, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    Teuchos::ParameterList matrixParams, ifpack2Params;

    matrixParams.set("matrixType", "Laplace1D");
    matrixParams.set("nx", (GlobalOrdinal)20);  // needs to be even

    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildBlockMatrixAsPoint(matrixParams, Xpetra::UseTpetra);
    ifpack2Params.set("smoother: use blockcrsmatrix storage", true);

    Ifpack2Smoother smoother("RELAXATION", ifpack2Params);

    Level level;
    TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
    level.Set("A", A);
    smoother.Setup(level);

    TEST_EQUALITY(1, 1);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Ifpack2Smoother, BlockCrsMatrix_Relaxation_AsBlock, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    Teuchos::ParameterList matrixParams, ifpack2Params;

    matrixParams.set("matrixType", "Laplace1D");
    matrixParams.set("nx", (GlobalOrdinal)20);  // needs to be even

    RCP<Matrix> A = TestHelpers::TpetraTestFactory<SC, LO, GO, NO>::BuildBlockMatrix(matrixParams, Xpetra::UseTpetra);
    ifpack2Params.set("smoother: use blockcrsmatrix storage", true);

    Ifpack2Smoother smoother("RELAXATION", ifpack2Params);

    Level level;
    TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
    level.Set("A", A);
    smoother.Setup(level);

    TEST_EQUALITY(1, 1);
  }
}

#define MUELU_ETI_GROUP(SC, LO, GO, NO)                                                                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Ifpack2Smoother, NotSetup, SC, LO, GO, NO)                           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Ifpack2Smoother, HardCodedResult_Jacobi, SC, LO, GO, NO)             \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Ifpack2Smoother, HardCodedResult_GaussSeidel, SC, LO, GO, NO)        \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Ifpack2Smoother, HardCodedResult_GaussSeidel2, SC, LO, GO, NO)       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Ifpack2Smoother, HardCodedResult_Chebyshev, SC, LO, GO, NO)          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Ifpack2Smoother, HardCodedResult_ILU, SC, LO, GO, NO)                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Ifpack2Smoother, ILU_TwoSweeps, SC, LO, GO, NO)                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Ifpack2Smoother, BandedRelaxation, SC, LO, GO, NO)                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Ifpack2Smoother, TriDiRelaxation, SC, LO, GO, NO)                    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Ifpack2Smoother, BlockRelaxation_Autosize, SC, LO, GO, NO)           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Ifpack2Smoother, BlockCrsMatrix_Relaxation_ViaPoint, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Ifpack2Smoother, BlockCrsMatrix_Relaxation_AsBlock, SC, LO, GO, NO)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
