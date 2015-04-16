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
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Aggregates.hpp"

#include "MueLu_UseDefaultTypes.hpp"

namespace MueLuTests {

#include "MueLu_UseShortNames.hpp"

  // Little utility to generate uncoupled aggregates.
  RCP<Aggregates> gimmeUncoupledAggregates(const RCP<Matrix> & A, RCP<AmalgamationInfo> & amalgInfo, bool bPhase1 = true, bool bPhase2a = true, bool bPhase2b = true, bool bPhase3 = true)
  {
    Level level;
    TestHelpers::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
    RCP<CoalesceDropFactory> dropFact = rcp(new CoalesceDropFactory());
    dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

    std::cout << "Phase2a=" << bPhase2a << std::endl;

    // Setup aggregation factory (use default factory for graph)
    RCP<UncoupledAggregationFactory> aggFact = rcp(new UncoupledAggregationFactory());
    aggFact->SetFactory("Graph", dropFact);
    aggFact->SetParameter("aggregation: mode",Teuchos::ParameterEntry(std::string("new")));
    aggFact->SetParameter("aggregation: max agg size",Teuchos::ParameterEntry(3));
    aggFact->SetParameter("aggregation: min agg size",Teuchos::ParameterEntry(3));
    aggFact->SetParameter("aggregation: max selected neighbors",Teuchos::ParameterEntry(0));
    aggFact->SetParameter("aggregation: ordering",Teuchos::ParameterEntry(std::string("natural")));
    aggFact->SetParameter("aggregation: enable phase 1",  Teuchos::ParameterEntry(bPhase1));
    aggFact->SetParameter("aggregation: enable phase 2a", Teuchos::ParameterEntry(bPhase2a));
    aggFact->SetParameter("aggregation: enable phase 2b", Teuchos::ParameterEntry(bPhase2b));
    aggFact->SetParameter("aggregation: enable phase 3",  Teuchos::ParameterEntry(bPhase3));

    level.Request("Aggregates", aggFact.get());
    level.Request("UnAmalgamationInfo", amalgFact.get());

    level.Request(*aggFact);
    aggFact->Build(level);
    RCP<Aggregates> aggregates = level.Get<RCP<Aggregates> >("Aggregates",aggFact.get()); // fix me
    amalgInfo = level.Get<RCP<AmalgamationInfo> >("UnAmalgamationInfo",amalgFact.get()); // fix me
    level.Release("UnAmalgamationInfo", amalgFact.get());
    level.Release("Aggregates", aggFact.get());
    return aggregates;
  }  // gimmeUncoupledAggregates

  // Little utility to generate coupled aggregates.
  RCP<Aggregates> gimmeCoupledAggregates(const RCP<Matrix> & A, RCP<AmalgamationInfo> & amalgInfo)
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
  }  // gimmeCoupledAggregates

  TEUCHOS_UNIT_TEST(Aggregates, JustAggregation)
  {
    out << "version: " << MueLu::Version() << std::endl;
    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(15);
    RCP<AmalgamationInfo> amalgInfo;
    RCP<Aggregates> aggregates = gimmeCoupledAggregates(A, amalgInfo);
    TEST_EQUALITY(aggregates != Teuchos::null, true);
    TEST_EQUALITY(aggregates->AggregatesCrossProcessors(),true);
  }

  ///////////////////////////////////////////////////////////////////////////

  TEUCHOS_UNIT_TEST(Aggregates, GetNumAggregates)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
    RCP<const Map> rowmap = A->getRowMap();
    RCP<AmalgamationInfo> amalgInfo;
    RCP<Aggregates> aggregates = gimmeCoupledAggregates(A, amalgInfo);
    GO numAggs = aggregates->GetNumAggregates();
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    TEST_EQUALITY(aggregates->AggregatesCrossProcessors(),true);

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

  TEUCHOS_UNIT_TEST(Aggregates, JustUncoupledAggregation)
  {
    out << "version: " << MueLu::Version() << std::endl;
    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(15);
    RCP<AmalgamationInfo> amalgInfo;
    RCP<Aggregates> aggregates = gimmeUncoupledAggregates(A, amalgInfo);
    TEST_EQUALITY(aggregates != Teuchos::null, true);
    TEST_EQUALITY(aggregates->AggregatesCrossProcessors(),false);
  }

  ///////////////////////////////////////////////////////////////////////////

  TEUCHOS_UNIT_TEST(Aggregates, GetNumUncoupledAggregates)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
    RCP<const Map> rowmap = A->getRowMap();
    RCP<AmalgamationInfo> amalgInfo;
    RCP<Aggregates> aggregates = gimmeUncoupledAggregates(A, amalgInfo);
    GO numAggs = aggregates->GetNumAggregates();
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    TEST_EQUALITY(aggregates->AggregatesCrossProcessors(),false);

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

  TEUCHOS_UNIT_TEST(Aggregates, UncoupledPhase1)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
    RCP<const Map> rowmap = A->getRowMap();
    RCP<AmalgamationInfo> amalgInfo;
    RCP<Aggregates> aggregates = gimmeUncoupledAggregates(A, amalgInfo,true,false,false,false);
    GO numAggs = aggregates->GetNumAggregates();
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    TEST_EQUALITY(aggregates->AggregatesCrossProcessors(),false);

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

  } //UncoupledPhase1

  TEUCHOS_UNIT_TEST(Aggregates, UncoupledPhase2)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
    RCP<const Map> rowmap = A->getRowMap();
    RCP<AmalgamationInfo> amalgInfo;
    bool bSuccess = true;
    TEUCHOS_TEST_THROW(gimmeUncoupledAggregates(A, amalgInfo,false,true,true,false),
      MueLu::Exceptions::RuntimeError, out, bSuccess);
    TEST_EQUALITY(bSuccess, true);
  } //UncoupledPhase2

  TEUCHOS_UNIT_TEST(Aggregates, UncoupledPhase3)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
    RCP<const Map> rowmap = A->getRowMap();
    RCP<AmalgamationInfo> amalgInfo;
    RCP<Aggregates> aggregates = gimmeUncoupledAggregates(A, amalgInfo,false,false,false,true);
    GO numAggs = aggregates->GetNumAggregates();
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    TEST_EQUALITY(aggregates->AggregatesCrossProcessors(),false);

    ArrayRCP<LO> aggSizes = Teuchos::ArrayRCP<LO>(numAggs);
    ArrayRCP<LO> aggStart;
    ArrayRCP<GO> aggToRowMap;
    amalgInfo->UnamalgamateAggregates(*aggregates, aggStart, aggToRowMap);
    for (LO i = 0; i < numAggs; ++i)
      aggSizes[i] = aggStart[i+1] - aggStart[i];

    bool foundAggNotSize2=false;
    for (int i=0; i<aggSizes.size(); ++i)
      if (aggSizes[i] != 2) {
        foundAggNotSize2=true;
        break;
      }

    switch (comm->getSize()) {

      case 1 :
        TEST_EQUALITY(numAggs, 18);
        TEST_EQUALITY(foundAggNotSize2, false);
        break;

      case 2:
        TEST_EQUALITY(numAggs, 9);
        TEST_EQUALITY(foundAggNotSize2, false);
        break;

      case 3:
        TEST_EQUALITY(numAggs, 6);
        TEST_EQUALITY(foundAggNotSize2, false);
        break;

      case 4:
        TEST_EQUALITY(numAggs, 4);
        TEST_EQUALITY(foundAggNotSize2, true);
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

  } //UncoupledPhase3
} // namespace MueLuTests
