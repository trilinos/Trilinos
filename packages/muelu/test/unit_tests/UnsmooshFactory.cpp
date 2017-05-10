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
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_TentativePFactory.hpp"
#include "MueLu_UnsmooshFactory.hpp"

namespace MueLuTests {


  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(UnsmooshFactory, UnsmooshTentativeP, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);
    if (!TYPE_EQUAL(GO, int)) { out << "Skipping test for GO != int"        << std::endl; return; }
    out << "version: " << MueLu::Version() << std::endl;

    //RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    //Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

    typedef Teuchos::ScalarTraits<Scalar> TST;
    typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;

    Level fineLevel, coarseLevel;
    test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);

    RCP<Matrix> A = test_factory::Build1DPoisson(199);
    fineLevel.Set("A", A);



    RCP<TentativePFactory> tentativePFact = rcp(new TentativePFactory());

    int maxDofPerNode = 2;
    Teuchos::Array<char> dofStatus = Teuchos::Array<char>(A->getRangeMap()->getNodeNumElements() * maxDofPerNode,'s');
    coarseLevel.Set("DofStatus", dofStatus);

    RCP<UnsmooshFactory> unsmooFact = Teuchos::rcp(new UnsmooshFactory());
    unsmooFact->SetFactory("P", tentativePFact);
    unsmooFact->SetFactory("DofStatus", MueLu::NoFactory::getRCP());
    unsmooFact->SetParameter("maxDofPerNode", Teuchos::ParameterEntry(maxDofPerNode));
    unsmooFact->SetParameter("fineIsPadded", Teuchos::ParameterEntry(true));

    coarseLevel.Request("P", unsmooFact.get());
    Teuchos::RCP<Matrix> test = coarseLevel.Get<RCP<Matrix> >("P",unsmooFact.get());


    //VariableDofLaplacianFactory lapFact;

    //l->Request("A",&lapFact);


    //lapA->describe(out, Teuchos::VERB_EXTREME);

    //TEST_EQUALITY(lapA->getRowMap()->isSameAs(*A->getRowMap()),true);
  } // VarLaplConstructor

#  define MUELU_ETI_GROUP(SC, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(UnsmooshFactory, UnsmooshTentativeP, SC, LO, GO, Node) \


#include <MueLu_ETI_4arg.hpp>
}
