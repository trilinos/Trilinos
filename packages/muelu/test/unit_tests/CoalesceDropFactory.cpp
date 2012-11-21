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
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_PreDropFunctionConstVal.hpp"
#include "MueLu_Graph.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  TEUCHOS_UNIT_TEST(CoalesceDropFactory, Constructor)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<CoalesceDropFactory> coalesceDropFact = rcp(new CoalesceDropFactory());
    TEST_EQUALITY(coalesceDropFact != Teuchos::null, true);

  } //Constructor

  TEUCHOS_UNIT_TEST(CoalesceDropFactory, Build)
  {
    out << "version: " << MueLu::Version() << std::endl;

    Level fineLevel;
    TestHelpers::TestFactory<SC,LO,GO,NO,LMO>::createSingleLevelHierarchy(fineLevel);

    RCP<Matrix> A = TestHelpers::TestFactory<SC,LO,GO,NO,LMO>::Build1DPoisson(36);
    fineLevel.Set("A", A);

    CoalesceDropFactory coalesceDropFact;
    fineLevel.Request("Graph",&coalesceDropFact);
    coalesceDropFact.Build(fineLevel);
    //FIXME how do we verify that this is correct?
  } //Build

  TEUCHOS_UNIT_TEST(CoalesceDropFactory, PreDrop)
  {
    out << "version: " << MueLu::Version() << std::endl;

    Level fineLevel;
    TestHelpers::TestFactory<SC,LO,GO,NO,LMO>::createSingleLevelHierarchy(fineLevel);

    RCP<Matrix> A = TestHelpers::TestFactory<SC,LO,GO,NO,LMO>::Build1DPoisson(3);
    fineLevel.Set("A", A);
    A->describe(out,Teuchos::VERB_EXTREME);

    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetVerbLevel(MueLu::Extreme);
    dropFact.SetPreDropFunction(rcp(new PreDropFunctionConstVal(0.00001)));

    fineLevel.Request("Graph", &dropFact);

    dropFact.Build(fineLevel);

    fineLevel.print(out);
    RCP<Graph> graph = fineLevel.Get<RCP<Graph> >("Graph", &dropFact);

    std::cout << graph->GetDomainMap()->getGlobalNumElements() << std::endl;
    graph->print(out, MueLu::Debug);

//    TEST_EQUALITY(1 == 0, true);

  } //PreDrop

} // namespace MueLuTests

