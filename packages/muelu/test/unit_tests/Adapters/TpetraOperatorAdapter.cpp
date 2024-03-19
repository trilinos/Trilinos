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
#include <Teuchos_ScalarTraits.hpp>

#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <Xpetra_MultiVectorFactory.hpp>

#include <MueLu_FactoryManagerBase.hpp>
#include <MueLu_Hierarchy.hpp>
#include <MueLu_PFactory.hpp>
#include <MueLu_SaPFactory.hpp>
#include <MueLu_TransPFactory.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_AmesosSmoother.hpp>
#include <MueLu_TrilinosSmoother.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_AmesosSmoother.hpp>
#include <MueLu_Utilities.hpp>
#include <MueLu_TpetraOperator.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TpetraOperator, Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

#if defined(HAVE_MUELU_IFPACK2) && defined(HAVE_MUELU_AMESOS2)
  typedef MueLu::Utilities<SC, LO, GO, NO> Utils;
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

    tH->apply(*(Utils::MV2TpetraMV(RHS1)), *(Utils::MV2NonConstTpetraMV(X1)));

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
#else
  out << "Skipping test because some required packages are not enabled (Tpetra, Ifpack2, Amesos2)." << std::endl;
#endif
}  // Apply

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TpetraOperator, Getters, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

#if defined(HAVE_MUELU_IFPACK2) && defined(HAVE_MUELU_AMESOS2)
  using Utils              = MueLu::Utilities<SC, LO, GO, NO>;
  using TpetraOperatorType = MueLu::TpetraOperator<SC, LO, GO, NO>;
  using MagnitudeType      = typename Teuchos::ScalarTraits<SC>::magnitudeType;

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
#else
  out << "Skipping test because some required packages are not enabled (Tpetra, Ifpack2, Amesos2)." << std::endl;
#endif
}

#define MUELU_ETI_GROUP(Scalar, LocalOrdinal, GlobalOrdinal, Node)                                       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TpetraOperator, Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TpetraOperator, Getters, Scalar, LocalOrdinal, GlobalOrdinal, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
