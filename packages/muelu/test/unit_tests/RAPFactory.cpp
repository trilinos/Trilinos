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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include "Teuchos_UnitTestHarness.hpp"

#include "MueLu_config.hpp"

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_Utilities.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_FactoryManager.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  TEUCHOS_UNIT_TEST(RAPFactory, Constructor)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<RAPFactory> rapFactory = rcp(new RAPFactory);
    TEST_EQUALITY(rapFactory != Teuchos::null, true);

    out << *rapFactory << std::endl;
  } // Constructor test

  TEUCHOS_UNIT_TEST(RAPFactory, Correctness)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    Level fineLevel, coarseLevel; TestHelpers::Factory<SC, LO, GO, NO, LMO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

    RCP<Matrix> Op = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(27*comm->getSize());
    fineLevel.Set("A",Op);

    TentativePFactory tentpFactory;
    SaPFactory sapFactory(rcpFromRef(tentpFactory));
    TransPFactory transPFactory(rcpFromRef(sapFactory)); //todo:rcpFromRef

    coarseLevel.Request("P",&sapFactory);
    coarseLevel.Request("R",&transPFactory);

    coarseLevel.Request(sapFactory);
    coarseLevel.Request(transPFactory);
    sapFactory.Build(fineLevel,coarseLevel);
    transPFactory.Build(fineLevel,coarseLevel);

    RAPFactory rap(rcpFromRef(sapFactory), rcpFromRef(transPFactory));
    coarseLevel.Request(rap);

    coarseLevel.Request("A",&rap);
    rap.Build(fineLevel,coarseLevel);

    RCP<Matrix> A = fineLevel.Get< RCP<Matrix> >("A");
    RCP<Matrix> P = coarseLevel.Get< RCP<Matrix> >("P", &sapFactory);
    RCP<Matrix> R = coarseLevel.Get< RCP<Matrix> >("R", &transPFactory);

    RCP<MultiVector> workVec1 = MultiVectorFactory::Build(P->getRangeMap(),1);
    RCP<MultiVector> workVec2 = MultiVectorFactory::Build(Op->getRangeMap(),1);
    RCP<MultiVector> result1 = MultiVectorFactory::Build(R->getRangeMap(),1);
    RCP<MultiVector> X = MultiVectorFactory::Build(P->getDomainMap(),1);
    X->randomize();

    //Calculate result1 = R*(A*(P*X))
    P->apply(*X,*workVec1,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
    Op->apply(*workVec1,*workVec2,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
    R->apply(*workVec2,*result1,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

    RCP<Matrix> coarseOp = coarseLevel.Get< RCP<Matrix> >("A", &rap);

    //Calculate result2 = (R*A*P)*X
    RCP<MultiVector> result2 = MultiVectorFactory::Build(R->getRangeMap(),1);
    coarseOp->apply(*X,*result2,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

    Teuchos::Array<ST::magnitudeType> normX(1), normResult1(1),normResult2(1);
    X->norm2(normX);
    out << "This test checks the correctness of the Galerkin triple "
        << "matrix product by comparing (RAP)*X to R(A(P*X))." << std::endl;
    out << "||X||_2 = " << normX << std::endl;
    result1->norm2(normResult1);
    result2->norm2(normResult2);
    TEST_FLOATING_EQUALITY(normResult1[0], normResult2[0], 1e-12);

  } // Correctness test

  TEUCHOS_UNIT_TEST(RAPFactory, ImplicitTranspose)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    if (comm->getSize() > 1 && TestHelpers::Parameters::getLib() == Xpetra::UseEpetra ) {
      out << "Skipping ImplicitTranspose test for Epetra and #proc>1" << std::endl;
      return;
    }

    // build test-specific default factory manager
    RCP<FactoryManager> defManager = rcp(new FactoryManager());
    defManager->SetFactory("A", rcp(MueLu::NoFactory::get(),false));         // dummy factory for A
    defManager->SetFactory("Nullspace", rcp(new NullspaceFactory()));        // real null space factory for Ptent
    defManager->SetFactory("Graph", rcp(new CoalesceDropFactory()));         // real graph factory for Ptent
    defManager->SetFactory("Aggregates", rcp(new UCAggregationFactory()));   // real aggregation factory for Ptent

    Level fineLevel, coarseLevel;
    TestHelpers::Factory<SC, LO, GO, NO, LMO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

    // overwrite default factory manager
    fineLevel.SetFactoryManager(defManager);
    coarseLevel.SetFactoryManager(defManager);

    RCP<Matrix> Op = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(19*comm->getSize());
    fineLevel.Set("A",Op);

    TentativePFactory tentpFactory;
    SaPFactory sapFactory(rcpFromRef(tentpFactory));
    TransPFactory transPFactory(rcpFromRef(sapFactory));

    coarseLevel.Request("P", &sapFactory);
    coarseLevel.Request("R", &transPFactory);

    coarseLevel.Request(sapFactory);
    coarseLevel.Request(transPFactory);
    sapFactory.Build(fineLevel, coarseLevel);
    transPFactory.Build(fineLevel,coarseLevel);
    RAPFactory rap(rcpFromRef(sapFactory), rcpFromRef(transPFactory));

    coarseLevel.Request("A", &rap);

    rap.SetImplicitTranspose(true);
    coarseLevel.Request(rap);
    rap.Build(fineLevel,coarseLevel);

    RCP<Matrix> A = fineLevel.Get< RCP<Matrix> >("A");
    RCP<Matrix> P = coarseLevel.Get< RCP<Matrix> >("P", &sapFactory);
    RCP<Matrix> R = coarseLevel.Get< RCP<Matrix> >("R", &transPFactory);

    //std::string filename = "A.dat";
    //Utils::Write(filename,Op);
    //filename = "P.dat";
    //Utils::Write(filename,P);

    RCP<MultiVector> workVec1 = MultiVectorFactory::Build(P->getRangeMap(),1);
    RCP<MultiVector> workVec2 = MultiVectorFactory::Build(Op->getRangeMap(),1);
    RCP<MultiVector> result1 = MultiVectorFactory::Build(P->getDomainMap(),1);
    RCP<MultiVector> X = MultiVectorFactory::Build(P->getDomainMap(),1);
    X->randomize();
    //out.precision(12);
    //out.setOutputToRootOnly(-1);
    //X->describe(out,Teuchos::VERB_EXTREME);

    //Calculate result1 = P^T*(A*(P*X))
    P->apply(*X,*workVec1,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
    Op->apply(*workVec1,*workVec2,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
    P->apply(*workVec2,*result1,Teuchos::TRANS,(SC)1.0,(SC)0.0);

    RCP<Matrix> coarseOp = coarseLevel.Get< RCP<Matrix> >("A", &rap);

    //Calculate result2 = (R*A*P)*X
    RCP<MultiVector> result2 = MultiVectorFactory::Build(P->getDomainMap(),1);
    coarseOp->apply(*X,*result2,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

    Teuchos::Array<ST::magnitudeType> normX(1), normResult1(1),normResult2(1);
    X->norm2(normX);
    out << "This test checks the correctness of the Galerkin triple "
        << "matrix product by comparing (RAP)*X to R(A(P*X)), where R is the implicit tranpose of P." << std::endl;
    out << "||X||_2 = " << normX << std::endl;
    result1->norm2(normResult1);
    result2->norm2(normResult2);
    TEST_FLOATING_EQUALITY(normResult1[0], normResult2[0], 1e-12);

  } // Correctness test

} // namespace MueLuTests

