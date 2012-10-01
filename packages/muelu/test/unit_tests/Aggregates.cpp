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

#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Aggregates.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  // Little utility to generate aggregates.
  RCP<Aggregates> gimmeAggregates(RCP<Matrix> const &A, RCP<AmalgamationInfo> & amalgInfo)
  {
    Level level;
    TestHelpers::Factory<SC,LO,GO,NO,LMO>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    AmalgamationFactory amalgFact;
    CoalesceDropFactory dropFact(Teuchos::null, Teuchos::rcpFromRef(amalgFact));
    
    // Setup aggregation factory (use default factory for graph)
    UCAggregationFactory aggFact(Teuchos::rcpFromRef(dropFact));
    aggFact.SetMinNodesPerAggregate(3);
    aggFact.SetMaxNeighAlreadySelected(0);
    aggFact.SetOrdering(MueLu::AggOptions::NATURAL);
    aggFact.SetPhase3AggCreation(0.5);

    level.Request("Aggregates", &aggFact);
    level.Request("UnAmalgamationInfo", &amalgFact);

    level.Request(aggFact);
    aggFact.Build(level);
    RCP<Aggregates> aggregates = level.Get<RCP<Aggregates> >("Aggregates",&aggFact); // fix me
    amalgInfo = level.Get<RCP<AmalgamationInfo> >("UnAmalgamationInfo",&amalgFact); // fix me
    level.Release("UnAmalgamationInfo", &amalgFact);
    level.Release("Aggregates", &aggFact);
    return aggregates;
  }  // gimmeAggregates

  void ComputeAggregateSizes(const Aggregates & aggregates, const AmalgamationInfo & amalgInfo, Teuchos::ArrayRCP<LocalOrdinal> & aggSizes) {
    int myPid = aggregates.GetMap()->getComm()->getRank();
    Teuchos::ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);
    Teuchos::ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
    LO size = procWinner.size();

    for (LO i = 0; i< aggregates.GetNumAggregates(); ++i) aggSizes[i] = 0;
    for (LO lnode = 0; lnode < size; ++lnode) {
      LO myAgg = vertex2AggId[lnode];
      if (procWinner[lnode] == myPid) {
        GO gnodeid = aggregates.GetMap()->getGlobalElement(lnode);

        std::vector<GO> gDofIds = (*(amalgInfo.GetGlobalAmalgamationParams()))[gnodeid];
        aggSizes[myAgg] += Teuchos::as<LO>(gDofIds.size());
      }
    }
  }

  void ComputeAggregateToRowMap(const Aggregates& aggregates, const AmalgamationInfo& amalgInfo, const Teuchos::ArrayRCP<LocalOrdinal> & aggSizes, Teuchos::ArrayRCP<Teuchos::ArrayRCP<GlobalOrdinal> > & aggToRowMap) {
    int myPid = aggregates.GetMap()->getComm()->getRank();
    Teuchos::ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);
    Teuchos::ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
    LO size = procWinner.size();

    // initialize array aggToRowMap with empty arrays for each aggregate (with correct aggSize)
    LO t = 0;
    for (ArrayRCP<ArrayRCP<GO> >::iterator a2r = aggToRowMap.begin(); a2r!=aggToRowMap.end(); ++a2r) {
      *a2r = ArrayRCP<GO>(aggSizes[t++]);
    }

    // count, how many dofs have been recorded for each aggregate
    ArrayRCP<LO> numDofs(aggregates.GetNumAggregates(),0); // empty array with number of Dofs for each aggregate

    for (LO lnode = 0; lnode < size; ++lnode) {
      LO myAgg = vertex2AggId[lnode];
      if (procWinner[lnode] == myPid) {
        GO gnodeid = aggregates.GetMap()->getGlobalElement(lnode);
        std::vector<GO> gDofIds = (*(amalgInfo.GetGlobalAmalgamationParams()))[gnodeid];
        for (LO gDofId=0; gDofId < Teuchos::as<LO>(gDofIds.size()); gDofId++) {
          aggToRowMap[ myAgg ][ numDofs[myAgg] ] = gDofIds[gDofId]; // fill aggToRowMap structure
          ++(numDofs[myAgg]);
        }
      }
    }
  }
  
  TEUCHOS_UNIT_TEST(Aggregates, JustAggregation)
  {
    out << "version: " << MueLu::Version() << std::endl;
    RCP<Matrix> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(15);
    RCP<AmalgamationInfo> amalgInfo = Teuchos::rcp(new AmalgamationInfo);
    RCP<Aggregates> aggregates = gimmeAggregates(A, amalgInfo);
    TEST_EQUALITY(aggregates != Teuchos::null, true);
  }

///////////////////////////////////////////////////////////////////////////

  TEUCHOS_UNIT_TEST(Aggregates, GetNumAggregates)
  {
      out << "version: " << MueLu::Version() << std::endl;

      RCP<Matrix> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(36);
      RCP<const Map> rowmap = A->getRowMap();
      RCP<AmalgamationInfo> amalgInfo = Teuchos::rcp(new AmalgamationInfo);
      RCP<Aggregates> aggregates = gimmeAggregates(A, amalgInfo);
      GO numAggs = aggregates->GetNumAggregates();
      RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

      ArrayRCP<LO> aggSizes = Teuchos::ArrayRCP<LO>(numAggs);
      ComputeAggregateSizes(*aggregates, *amalgInfo, aggSizes);
      
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

      ArrayRCP< ArrayRCP<GO> > aggToRowMap(numAggs);
      ComputeAggregateToRowMap(*aggregates, *amalgInfo, aggSizes, aggToRowMap);
      int root = out.getOutputToRootOnly();
      out.setOutputToRootOnly(-1);
      for (int j=0; j<comm->getSize(); ++j) {
        if (comm->getRank() == j) {
            out << "++ pid " << j << " ++" << std::endl;
            out << "   num local DOFs = " << rowmap->getNodeNumElements() << std::endl;
          for (int i=0; i< aggToRowMap.size(); ++i) {
            out << "   aggregate " << i << ": ";
            for (int k=0; k< aggToRowMap[i].size(); ++k)
              out << aggToRowMap[i][k] << " ";
            out << std::endl;
          }
        }
        comm->barrier();
      }
      out.setOutputToRootOnly(root);

  } //GetNumAggregates


} // namespace MueLuTests
