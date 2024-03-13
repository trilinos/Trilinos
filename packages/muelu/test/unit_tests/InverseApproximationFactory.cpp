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
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_FancyOStream.hpp>

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_MapExtractor_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_IO.hpp>

#include <MueLu_InverseApproximationFactory.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(InverseApproximationFactory, InverseDiagonalConstructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

  {
    const int n            = 20;
    Teuchos::RCP<Matrix> A = MueLuTests::TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build1DPoisson(n, lib);

    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    RCP<InverseApproximationFactory> invapproxFact = rcp(new InverseApproximationFactory());
    invapproxFact->SetFactory("A", MueLu::NoFactory::getRCP());
    invapproxFact->SetParameter("inverse: approximation type", Teuchos::ParameterEntry(std::string("diagonal")));

    // request InverseApproximation operator
    level.Request("Ainv", invapproxFact.get());

    // generate Schur complement operator
    invapproxFact->Build(level);

    RCP<Matrix> Ainv = level.Get<RCP<Matrix> >("Ainv", invapproxFact.get());
    TEST_EQUALITY(Ainv.is_null(), false);
    TEST_EQUALITY(Ainv->getRangeMap()->getMinGlobalIndex(), comm->getRank() * int(n / comm->getSize()));
    TEST_EQUALITY(Ainv->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * int(n / comm->getSize()) + int(n / comm->getSize() - 1));
    TEST_EQUALITY(Ainv->getDomainMap()->getMinGlobalIndex(), comm->getRank() * int(n / comm->getSize()));
    TEST_EQUALITY(Ainv->getDomainMap()->getMaxGlobalIndex(), comm->getRank() * int(n / comm->getSize()) + int(n / comm->getSize() - 1));

    const RCP<Vector> AinvDiagonal = VectorFactory::Build(Ainv->getRangeMap(), true);
    Ainv->getLocalDiagCopy(*AinvDiagonal);
    Teuchos::ArrayRCP<const Scalar> AinvData = AinvDiagonal->getData(0);
    bool bCheck                              = true;
    for (int i = 0; i < int(n / comm->getSize()); i++)
      if (AinvData[i] != Teuchos::as<Scalar>(0.5)) bCheck = false;
    TEST_EQUALITY(bCheck, true);
  }

  {
    const int n            = 42;
    Teuchos::RCP<Matrix> A = MueLuTests::TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build2DPoisson(n, n, lib);

    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    RCP<InverseApproximationFactory> invapproxFact = rcp(new InverseApproximationFactory());
    invapproxFact->SetFactory("A", MueLu::NoFactory::getRCP());
    invapproxFact->SetParameter("inverse: approximation type", Teuchos::ParameterEntry(std::string("diagonal")));

    // request InverseApproximation operator
    level.Request("Ainv", invapproxFact.get());

    // generate Schur complement operator
    invapproxFact->Build(level);

    RCP<Matrix> Ainv = level.Get<RCP<Matrix> >("Ainv", invapproxFact.get());
    TEST_EQUALITY(Ainv.is_null(), false);
    TEST_EQUALITY(Ainv->getRangeMap()->getMinGlobalIndex(), comm->getRank() * int(n * n / comm->getSize()));
    TEST_EQUALITY(Ainv->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * int(n * n / comm->getSize()) + int(n * n / comm->getSize() - 1));
    TEST_EQUALITY(Ainv->getDomainMap()->getMinGlobalIndex(), comm->getRank() * int(n * n / comm->getSize()));
    TEST_EQUALITY(Ainv->getDomainMap()->getMaxGlobalIndex(), comm->getRank() * int(n * n / comm->getSize()) + int(n * n / comm->getSize() - 1));

    const RCP<Vector> AinvDiagonal = VectorFactory::Build(Ainv->getRangeMap(), true);
    Ainv->getLocalDiagCopy(*AinvDiagonal);
    Teuchos::ArrayRCP<const Scalar> AinvData = AinvDiagonal->getData(0);
    bool bCheck                              = true;
    for (int i = 0; i < int(n * n / comm->getSize()); i++)
      if (AinvData[i] != Teuchos::as<Scalar>(0.25)) bCheck = false;
    TEST_EQUALITY(bCheck, true);
  }

}  // InverseDiagonalConstructor

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(InverseApproximationFactory, InverseLumpingConstructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

  {
    const int n            = 20;
    Teuchos::RCP<Matrix> A = MueLuTests::TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build1DPoisson(n, lib);
    A->scale(0.5);

    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    RCP<InverseApproximationFactory> invapproxFact = rcp(new InverseApproximationFactory());
    invapproxFact->SetFactory("A", MueLu::NoFactory::getRCP());
    invapproxFact->SetParameter("inverse: approximation type", Teuchos::ParameterEntry(std::string("lumping")));

    // request InverseApproximation operator
    level.Request("Ainv", invapproxFact.get());

    // generate Schur complement operator
    invapproxFact->Build(level);

    RCP<Matrix> Ainv = level.Get<RCP<Matrix> >("Ainv", invapproxFact.get());
    TEST_EQUALITY(Ainv.is_null(), false);
    TEST_EQUALITY(Ainv->getRangeMap()->getMinGlobalIndex(), comm->getRank() * int(n / comm->getSize()));
    TEST_EQUALITY(Ainv->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * int(n / comm->getSize()) + int(n / comm->getSize() - 1));
    TEST_EQUALITY(Ainv->getDomainMap()->getMinGlobalIndex(), comm->getRank() * int(n / comm->getSize()));
    TEST_EQUALITY(Ainv->getDomainMap()->getMaxGlobalIndex(), comm->getRank() * int(n / comm->getSize()) + int(n / comm->getSize() - 1));

    const RCP<Vector> AinvDiagonal = VectorFactory::Build(Ainv->getRangeMap(), true);
    Ainv->getLocalDiagCopy(*AinvDiagonal);
    Teuchos::ArrayRCP<const Scalar> AinvData = AinvDiagonal->getData(0);
    bool bCheck                              = false;
    for (int i = 0; i < int(n / comm->getSize()); i++)
      if (std::abs(AinvData[i] - Teuchos::as<Scalar>(0.66666666667)) < std::abs(Teuchos::as<Scalar>(1e-8)) ||
          std::abs(AinvData[i] - Teuchos::as<Scalar>(0.5)) < std::abs(Teuchos::as<Scalar>(1e-8)))
        bCheck = true;
    TEST_EQUALITY(bCheck, true);
  }

  {
    const int n            = 42;
    Teuchos::RCP<Matrix> A = MueLuTests::TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build2DPoisson(n, n, lib);

    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    RCP<InverseApproximationFactory> invapproxFact = rcp(new InverseApproximationFactory());
    invapproxFact->SetFactory("A", MueLu::NoFactory::getRCP());
    invapproxFact->SetParameter("inverse: approximation type", Teuchos::ParameterEntry(std::string("lumping")));

    // request InverseApproximation operator
    level.Request("Ainv", invapproxFact.get());

    // generate Schur complement operator
    invapproxFact->Build(level);

    RCP<Matrix> Ainv = level.Get<RCP<Matrix> >("Ainv", invapproxFact.get());
    TEST_EQUALITY(Ainv.is_null(), false);
    TEST_EQUALITY(Ainv->getRangeMap()->getMinGlobalIndex(), comm->getRank() * int(n * n / comm->getSize()));
    TEST_EQUALITY(Ainv->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * int(n * n / comm->getSize()) + int(n * n / comm->getSize() - 1));
    TEST_EQUALITY(Ainv->getDomainMap()->getMinGlobalIndex(), comm->getRank() * int(n * n / comm->getSize()));
    TEST_EQUALITY(Ainv->getDomainMap()->getMaxGlobalIndex(), comm->getRank() * int(n * n / comm->getSize()) + int(n * n / comm->getSize() - 1));

    const RCP<Vector> AinvDiagonal = VectorFactory::Build(Ainv->getRangeMap(), true);
    Ainv->getLocalDiagCopy(*AinvDiagonal);
    Teuchos::ArrayRCP<const Scalar> AinvData = AinvDiagonal->getData(0);
    bool bCheck                              = false;
    for (int i = 0; i < int(n * n / comm->getSize()); i++)
      if (std::abs(AinvData[i] - Teuchos::as<Scalar>(0.166667)) < std::abs(Teuchos::as<Scalar>(1e-6)) ||
          std::abs(AinvData[i] - Teuchos::as<Scalar>(0.142857)) < std::abs(Teuchos::as<Scalar>(1e-6)))
        bCheck = true;
    TEST_EQUALITY(bCheck, true);
  }

}  // InverseLumpingConstructor

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(InverseApproximationFactory, InverseSpaiConstructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  using TST            = Teuchos::ScalarTraits<SC>;
  using magnitude_type = typename TST::magnitudeType;
  using TMT            = Teuchos::ScalarTraits<magnitude_type>;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

  {
    const int n            = 20;
    Teuchos::RCP<Matrix> A = MueLuTests::TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build1DPoisson(n, lib);

    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    RCP<InverseApproximationFactory> invapproxFact = rcp(new InverseApproximationFactory());
    invapproxFact->SetFactory("A", MueLu::NoFactory::getRCP());
    invapproxFact->SetParameter("inverse: approximation type", Teuchos::ParameterEntry(std::string("sparseapproxinverse")));

    // request InverseApproximation operator
    level.Request("Ainv", invapproxFact.get());

    // generate Schur complement operator
    invapproxFact->Build(level);

    RCP<Matrix> Ainv = level.Get<RCP<Matrix> >("Ainv", invapproxFact.get());
    TEST_EQUALITY(Ainv.is_null(), false);
    TEST_EQUALITY(Ainv->getRangeMap()->getMinGlobalIndex(), comm->getRank() * int(n / comm->getSize()));
    TEST_EQUALITY(Ainv->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * int(n / comm->getSize()) + int(n / comm->getSize() - 1));
    TEST_EQUALITY(Ainv->getDomainMap()->getMinGlobalIndex(), comm->getRank() * int(n / comm->getSize()));
    TEST_EQUALITY(Ainv->getDomainMap()->getMaxGlobalIndex(), comm->getRank() * int(n / comm->getSize()) + int(n / comm->getSize() - 1));
    TEST_FLOATING_EQUALITY(Ainv->getFrobeniusNorm(), 3.037251711528645, 1e2 * TMT::eps());
  }

  {
    const int n            = 42;
    Teuchos::RCP<Matrix> A = MueLuTests::TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build2DPoisson(n, n, lib);

    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    RCP<InverseApproximationFactory> invapproxFact = rcp(new InverseApproximationFactory());
    invapproxFact->SetFactory("A", MueLu::NoFactory::getRCP());
    invapproxFact->SetParameter("inverse: approximation type", Teuchos::ParameterEntry(std::string("sparseapproxinverse")));

    // request InverseApproximation operator
    level.Request("Ainv", invapproxFact.get());

    // generate Schur complement operator
    invapproxFact->Build(level);

    RCP<Matrix> Ainv = level.Get<RCP<Matrix> >("Ainv", invapproxFact.get());
    TEST_EQUALITY(Ainv.is_null(), false);
    TEST_EQUALITY(Ainv->getRangeMap()->getMinGlobalIndex(), comm->getRank() * int(n * n / comm->getSize()));
    TEST_EQUALITY(Ainv->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * int(n * n / comm->getSize()) + int(n * n / comm->getSize() - 1));
    TEST_EQUALITY(Ainv->getDomainMap()->getMinGlobalIndex(), comm->getRank() * int(n * n / comm->getSize()));
    TEST_EQUALITY(Ainv->getDomainMap()->getMaxGlobalIndex(), comm->getRank() * int(n * n / comm->getSize()) + int(n * n / comm->getSize() - 1));
    TEST_FLOATING_EQUALITY(Ainv->getFrobeniusNorm(), 12.41994595675205, 1e3 * TMT::eps());
  }

  // test with a highly nonsymmetric matrix
  {
    using STS = Teuchos::ScalarTraits<SC>;

    // Don't test for complex - matrix reader won't work
    if (STS::isComplex) {
      success = true;
      return;
    }
    RCP<Matrix> A = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read("TestMatrices/nonsym.mm", lib, comm);

    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    RCP<InverseApproximationFactory> invapproxFact = rcp(new InverseApproximationFactory());
    invapproxFact->SetFactory("A", MueLu::NoFactory::getRCP());
    invapproxFact->SetParameter("inverse: approximation type", Teuchos::ParameterEntry(std::string("sparseapproxinverse")));

    // request InverseApproximation operator
    level.Request("Ainv", invapproxFact.get());

    // generate Schur complement operator
    invapproxFact->Build(level);

    RCP<Matrix> Ainv = level.Get<RCP<Matrix> >("Ainv", invapproxFact.get());
    TEST_EQUALITY(Ainv.is_null(), false);
    TEST_FLOATING_EQUALITY(Ainv->getFrobeniusNorm(), 0.1235706050986417, 1e2 * TMT::eps());
    // check values of first row only on root
    if (comm->getRank() == 0) {
      ArrayView<const LocalOrdinal> indices;
      ArrayView<const Scalar> values;
      Ainv->getLocalRowView(0, indices, values);
      TEST_FLOATING_EQUALITY(values[0], Teuchos::as<Scalar>(1.0000000000000002e-01), 1e2 * TMT::eps());
      TEST_FLOATING_EQUALITY(values[1], Teuchos::as<Scalar>(-1.6666666666666673e-02), 1e2 * TMT::eps());
      TEST_FLOATING_EQUALITY(values[2], Teuchos::as<Scalar>(4.6666666666666688e-03), 1e2 * TMT::eps());
    }
  }

  // Test pre and post filtering of approximate inverse
  {
    using STS = Teuchos::ScalarTraits<SC>;

    // Don't test for complex - matrix reader won't work
    if (STS::isComplex) {
      success = true;
      return;
    }
    RCP<Matrix> A = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read("TestMatrices/beam.mm", lib, comm);

    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    RCP<InverseApproximationFactory> invapproxFact = rcp(new InverseApproximationFactory());
    invapproxFact->SetFactory("A", MueLu::NoFactory::getRCP());
    invapproxFact->SetParameter("inverse: drop tolerance", Teuchos::ParameterEntry(Scalar(1e-8)));
    invapproxFact->SetParameter("inverse: approximation type", Teuchos::ParameterEntry(std::string("sparseapproxinverse")));

    // request InverseApproximation operator
    level.Request("Ainv", invapproxFact.get());

    // generate Schur complement operator
    invapproxFact->Build(level);

    RCP<Matrix> Ainv = level.Get<RCP<Matrix> >("Ainv", invapproxFact.get());
    TEST_EQUALITY(Ainv.is_null(), false);
    TEST_EQUALITY(Ainv->getGlobalNumEntries(), 115760);
    // 8.31688788510637e+06
    TEST_FLOATING_EQUALITY(Ainv->getFrobeniusNorm(), 8.31688788510637e+06, 1e5 * TMT::eps());
  }

}  // InverseSpaiConstructor

#define MUELU_ETI_GROUP(SC, LO, GO, Node)                                                                         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(InverseApproximationFactory, InverseDiagonalConstructor, SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(InverseApproximationFactory, InverseLumpingConstructor, SC, LO, GO, Node)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(InverseApproximationFactory, InverseSpaiConstructor, SC, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>
}  // namespace MueLuTests
