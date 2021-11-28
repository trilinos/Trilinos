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
// // Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_VerboseObject.hpp>

#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include <MueLu_config.hpp>
#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <MueLu_UncoupledAggregationFactory.hpp>

namespace MueLuTests {

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(UncoupledAggregationFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
    #   include <MueLu_UseShortNames.hpp>

      out << "version: " << MueLu::Version() << std::endl;
      RCP<UncoupledAggregationFactory> aggFact = rcp(new UncoupledAggregationFactory());
      TEST_EQUALITY(aggFact != Teuchos::null, true);
      TEST_THROW(aggFact->SetOrdering("unknown_ordering"), Teuchos::Exceptions::InvalidParameterValue);
      aggFact->SetOrdering("natural");
      TEST_EQUALITY(aggFact->GetOrdering() == "natural", true);
      aggFact->SetOrdering("graph");
      TEST_EQUALITY(aggFact->GetOrdering() == "graph", true);
      aggFact->SetOrdering("random");
      TEST_EQUALITY(aggFact->GetOrdering() == "random", true);

      aggFact->SetMaxNeighAlreadySelected(100);
      TEST_EQUALITY(aggFact->GetMaxNeighAlreadySelected() == 100, true);
      aggFact->SetMinNodesPerAggregate(10);
      TEST_EQUALITY(aggFact->GetMinNodesPerAggregate() == 10, true);

      auto validParamList = aggFact->GetValidParameterList();
      TEST_EQUALITY(validParamList != Teuchos::null, true);


  } // Constructor


  # define MUELU_ETI_GROUP(Scalar, LO, GO, Node) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(UncoupledAggregationFactory, Constructor, Scalar, LO, GO, Node)
# include <MueLu_ETI_4arg.hpp>

} // namespace MueLuTests

