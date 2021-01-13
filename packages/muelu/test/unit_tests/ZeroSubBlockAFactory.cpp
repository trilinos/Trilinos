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

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_BlockedCrsMatrix.hpp>

#include "MueLu_FactoryManager.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_ZeroSubBlockAFactory.hpp"

namespace MueLuTests {

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(ZeroSubBlockAFactory, ConstructZeroSubBlockAFactory, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);

    out << "version: " << MueLu::Version() << std::endl;

    RCP<ZeroSubBlockAFactory> subBlockAFactory = rcp(new ZeroSubBlockAFactory());
    TEST_ASSERT(!subBlockAFactory.is_null());
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(ZeroSubBlockAFactory, ExtractZeroMainDiagonalBlockFrom2x2SaddlePointSystem, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"

    using test_factory = TestHelpers::TestFactory<SC, LO, GO, NO>;

    MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra)

    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    // Let's create a 2x2 block matrix A, that represents meshtying of a 2D Poisson problem
    const GO numEleX = 10;
    RCP<BlockedCrsMatrix> saddlePointMatrix = test_factory::BuildPoissonSaddlePointMatrix(numEleX);

    TEST_ASSERT(!saddlePointMatrix.is_null());

    RCP<Level> myLevel = rcp(new Level());
    myLevel->Set("A", Teuchos::rcp_dynamic_cast<Matrix>(saddlePointMatrix));

    RCP<ZeroSubBlockAFactory> A11Fact = rcp(new ZeroSubBlockAFactory());
    A11Fact->SetFactory("A", MueLu::NoFactory::getRCP());
    A11Fact->SetParameter("block row", Teuchos::ParameterEntry(1));
    A11Fact->SetParameter("block col", Teuchos::ParameterEntry(1));

    myLevel->Request("A", A11Fact.get(), MueLu::NoFactory::get());
    TEST_ASSERT(myLevel->IsRequested("A", A11Fact.get()));

    A11Fact->Build(*myLevel);
    TEST_ASSERT(myLevel->IsAvailable("A", A11Fact.get()));

    RCP<Matrix> A11 = myLevel->Get<RCP<Matrix>>("A", A11Fact.get());

    myLevel->Release("A", A11Fact.get());
    TEST_ASSERT(!myLevel->IsAvailable("A", A11Fact.get()));

    TEST_ASSERT(A11->getRangeMap()->isSameAs(*saddlePointMatrix->getMatrix(1, 0)->getRangeMap()));
    TEST_ASSERT(A11->getDomainMap()->isSameAs(*saddlePointMatrix->getMatrix(0, 1)->getDomainMap()));
    TEST_ASSERT(A11->getRowMap()->isSameAs(*saddlePointMatrix->getMatrix(1, 0)->getRowMap()));
    TEST_ASSERT(A11->getColMap()->isSameAs(*saddlePointMatrix->getMatrix(0, 1)->getColMap()));
    TEST_EQUALITY_CONST(A11->GetFixedBlockSize(), saddlePointMatrix->getMatrix(1, 0)->GetFixedBlockSize());
    TEST_EQUALITY_CONST(A11->getGlobalNumEntries(), 0);
  }

#  define MUELU_ETI_GROUP(SC, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(ZeroSubBlockAFactory, ConstructZeroSubBlockAFactory, SC, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(ZeroSubBlockAFactory, ExtractZeroMainDiagonalBlockFrom2x2SaddlePointSystem, SC, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>
}


