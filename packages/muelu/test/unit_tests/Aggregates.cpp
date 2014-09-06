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

#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_CoupledAggregationFactory.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Aggregates.hpp"

#include "MueLu_UseDefaultTypes.hpp"

namespace MueLuTests {

#include "MueLu_UseShortNames.hpp"

  // Little utility to generate aggregates.
  RCP<Aggregates> gimmeAggregates(const RCP<Matrix> & A, RCP<AmalgamationInfo> & amalgInfo)
  {
    Level level;
    TestHelpers::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
    RCP<CoalesceDropFactory> dropFact = rcp(new CoalesceDropFactory());
    dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

    // Setup aggregation factory (use default factory for graph)
    RCP<CoupledAggregationFactory> aggFact = rcp(new CoupledAggregationFactory());
    aggFact->SetFactory("Graph", dropFact);
    aggFact->SetMinNodesPerAggregate(3);
    aggFact->SetMaxNeighAlreadySelected(0);
    aggFact->SetOrdering("natural");
    aggFact->SetPhase3AggCreation(0.5);

    level.Request("Aggregates", aggFact.get());
    level.Request("UnAmalgamationInfo", amalgFact.get());

    level.Request(*aggFact);
    aggFact->Build(level);
    RCP<Aggregates> aggregates = level.Get<RCP<Aggregates> >("Aggregates",aggFact.get()); // fix me
    amalgInfo = level.Get<RCP<AmalgamationInfo> >("UnAmalgamationInfo",amalgFact.get()); // fix me
    level.Release("UnAmalgamationInfo", amalgFact.get());
    level.Release("Aggregates", aggFact.get());
    return aggregates;
  }  // gimmeAggregates

  TEUCHOS_UNIT_TEST(Aggregates, JustAggregation)
  {
    out << "version: " << MueLu::Version() << std::endl;
    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(15);
    RCP<AmalgamationInfo> amalgInfo;
    RCP<Aggregates> aggregates = gimmeAggregates(A, amalgInfo);
    TEST_EQUALITY(aggregates != Teuchos::null, true);
  }

///////////////////////////////////////////////////////////////////////////

  TEUCHOS_UNIT_TEST(Aggregates, GetNumAggregates)
  {
      out << "version: " << MueLu::Version() << std::endl;

      RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
      RCP<const Map> rowmap = A->getRowMap();
      RCP<AmalgamationInfo> amalgInfo;
      RCP<Aggregates> aggregates = gimmeAggregates(A, amalgInfo);
      GO numAggs = aggregates->GetNumAggregates();
      RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

      ArrayRCP<LO> aggSizes = Teuchos::ArrayRCP<LO>(numAggs);
      ArrayRCP<LO> aggStart;
      ArrayRCP<GO> aggToRowMap;
      amalgInfo->UnamalgamateAggregates(*aggregates, aggStart, aggToRowMap);
      for (LO i = 0; i < numAggs; ++i)
        aggSizes[i] = aggStart[i+1] - aggStart[i];

      bool foundAggNotSize3=false;
      for (int i=0; i<aggSizes.size(); ++i)
        if (aggSizes[i] != 3) {
          foundAggNotSize3=true;
          break;
        }

      switch (comm->getSize()) {

        case 1 :
           TEST_EQUALITY(numAggs, 12);
           TEST_EQUALITY(foundAggNotSize3, false);
           break;

        case 2:
           TEST_EQUALITY(numAggs, 6);
           TEST_EQUALITY(foundAggNotSize3, false);
           break;

        case 3:
           TEST_EQUALITY(numAggs, 4);
           TEST_EQUALITY(foundAggNotSize3, false);
           break;

        case 4:
           TEST_EQUALITY(numAggs, 3);
           TEST_EQUALITY(foundAggNotSize3, false);
           break;

        default:
           std::string msg = "Only 1-4 MPI processes are supported.";
           //throw(MueLu::Exceptions::NotImplemented(msg));
           out << msg << std::endl;
           break;
      }

      //ArrayRCP< ArrayRCP<GO> > aggToRowMap(numAggs);
      int root = out.getOutputToRootOnly();
      out.setOutputToRootOnly(-1);
      for (int j=0; j<comm->getSize(); ++j) {
        if (comm->getRank() == j) {
            out << "++ pid " << j << " ++" << std::endl;
            out << "   num local DOFs = " << rowmap->getNodeNumElements() << std::endl;
          for (int i=0; i< numAggs; ++i) {
            out << "   aggregate " << i << ": ";
            for (int k=aggStart[i]; k< aggStart[i+1]; ++k)
              out << aggToRowMap[k] << " ";
            out << std::endl;
          }
        }
        comm->barrier();
      }
      out.setOutputToRootOnly(root);

  } //GetNumAggregates


} // namespace MueLuTests
