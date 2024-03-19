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
#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_AmesosSmoother.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_AmesosSmoother.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"
#include "MueLu_AMGXOperator.hpp"
#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_TpetraOperator.hpp"

#include <ctime>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(AMGXOperator, Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  // NOTE: This test only works for double/int/int
  if (!TYPE_EQUAL(double, Scalar) || !TYPE_EQUAL(int, LocalOrdinal) || !TYPE_EQUAL(int, GlobalOrdinal)) {
    out << "This test is enabled only for double/int/int" << std::endl;
    return;
  }

  typedef MueLu::AMGXOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> AMGXOperator;
  out << "version: " << MueLu::Version() << std::endl;

  if (TestHelpers::Parameters::getLib() == Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    int nx;
    // disable amgx test in parallel
    if (comm->getSize() > 1)
      return;
    else
      nx = 91;

    // matrix
    RCP<Matrix> Op                                                         = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build2DPoisson(nx, -1, Xpetra::UseTpetra);
    RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tpA = MueLu::Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstTpetraCrs(Op);
    RCP<Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tOp  = tpA;
    Teuchos::ParameterList params, dummyList;
    params.set("use external multigrid package", "amgx");
    Teuchos::ParameterList subList = params.sublist("amgx:params", false);
    params.sublist("amgx:params").set("json file", "test.json");
    RCP<MueLu::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tH = MueLu::CreateTpetraPreconditioner(tOp, params);

    RCP<AMGXOperator> aH = Teuchos::rcp_dynamic_cast<AMGXOperator>(tH);
    TEST_EQUALITY(aH->sizeA() == nx * nx / comm->getSize(), true);

    RCP<MultiVector> RHS = MultiVectorFactory::Build(Op->getRowMap(), 1);
    RCP<MultiVector> X   = MultiVectorFactory::Build(Op->getRowMap(), 1);

    // RHS=1, zero initial guess
    RHS->putScalar((double)1.0);
    X->putScalar((double)0.0);

    aH->apply(*(MueLu::Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MV2TpetraMV(RHS)), *(MueLu::Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MV2NonConstTpetraMV(X)));
    // if(comm->getSize() == 1) TEST_EQUALITY(aH->iters()==16,true);
    TEST_EQUALITY(aH->getStatus() == 0, true);

  } else {
    out << "This test is enabled only for linAlgebra=Tpetra." << std::endl;
  }

}  // Apply

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(AMGXOperator, Apply, Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
