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
#include <Teuchos_TestingHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <MueLu_FactoryManagerBase.hpp>
#include <MueLu_CoalesceDropFactory.hpp>
#include <MueLu_CoupledAggregationFactory.hpp>
#include <MueLu_StructuredAggregationFactory.hpp>
#include <MueLu_UncoupledAggregationFactory.hpp>
#include <MueLu_HybridAggregationFactory.hpp>
#include <MueLu_AmalgamationFactory.hpp>
#include <MueLu_AmalgamationInfo.hpp>
#include <MueLu_Aggregates.hpp>

namespace MueLuTests {

template
<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class AggregateGenerator {
  //typedef MueLu::AmalgamationInfo<LocalOrdinal,GlobalOrdinal,Node> amalgamation_info_type;
  //typedef MueLu::Aggregates<LocalOrdinal,GlobalOrdinal,Node> aggregates_type;
  //typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> xpetra_matrix_type;
#   include <MueLu_UseShortNames.hpp>

    public:

    // Little utility to generate uncoupled aggregates.
    static RCP<Aggregates>
    gimmeUncoupledAggregates(const RCP<Matrix> & A, RCP<AmalgamationInfo> & amalgInfo, bool bPhase1 = true, bool bPhase2a = true, bool bPhase2b = true, bool bPhase3 = true)
  //  RCP<MueLu::Aggregates> gimmeUncoupledAggregates(const RCP<xpetra_matrix_type> & A, RCP<AmalgamationInfo> & amalgInfo, bool bPhase1 = true, bool bPhase2a = true, bool bPhase2b = true, bool bPhase3 = true)
    {
      Level level;
      TestHelpers::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(level);
      level.Set("A", A);

      RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
      RCP<CoalesceDropFactory> dropFact = rcp(new CoalesceDropFactory());
      dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

      //std::cout << "Phase2a=" << bPhase2a << std::endl;

      // Setup aggregation factory (use default factory for graph)
      RCP<UncoupledAggregationFactory> aggFact = rcp(new UncoupledAggregationFactory());
      aggFact->SetFactory("Graph", dropFact);
      aggFact->SetParameter("aggregation: max agg size",Teuchos::ParameterEntry(3));
      aggFact->SetParameter("aggregation: min agg size",Teuchos::ParameterEntry(3));
      aggFact->SetParameter("aggregation: max selected neighbors",Teuchos::ParameterEntry(0));
      aggFact->SetParameter("aggregation: ordering",Teuchos::ParameterEntry(std::string("natural")));
      aggFact->SetParameter("aggregation: enable phase 1",  Teuchos::ParameterEntry(bPhase1));
      aggFact->SetParameter("aggregation: enable phase 2a", Teuchos::ParameterEntry(bPhase2a));
      aggFact->SetParameter("aggregation: enable phase 2b", Teuchos::ParameterEntry(bPhase2b));
      aggFact->SetParameter("aggregation: enable phase 3",  Teuchos::ParameterEntry(bPhase3));
      aggFact->SetParameter("aggregation: use interface aggregation",Teuchos::ParameterEntry(false));

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
    static RCP<Aggregates>
    gimmeCoupledAggregates(const RCP<Matrix> & A, RCP<AmalgamationInfo> & amalgInfo)
    //  RCP<Aggregates> gimmeCoupledAggregates(const RCP<xpetra_matrix_type> & A, RCP<AmalgamationInfo> & amalgInfo)
    {
      Level level;
      TestHelpers::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(level);
      level.Set("A", A);

      RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
      amalgFact->SetDefaultVerbLevel(MueLu::None);
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

    // Little utility to generate uncoupled aggregates.
    static RCP<Aggregates>
    gimmeStructuredAggregates(const RCP<Matrix>& A,
                              Array<GO> gNodesPerDir, Array<LO> lNodesPerDir, const bool coupled,
                              const LO numDimensions, const std::string meshLayout,
                              const Array<GO> meshData)
    {
      Level level;
      TestHelpers::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(level);

      const std::string coupling = (coupled ? "coupled" : "uncoupled");

      RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
      amalgFact->SetDefaultVerbLevel(MueLu::None);
      RCP<CoalesceDropFactory> dropFact = rcp(new CoalesceDropFactory());
      dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

      level.Set("A", A);
      level.Set("gNodesPerDim", gNodesPerDir);
      level.Set("lNodesPerDim", lNodesPerDir);
      level.Set("aggregation: mesh data", meshData);

      // Setup aggregation factory (use default factory for graph)
      RCP<StructuredAggregationFactory> aggFact = rcp(new StructuredAggregationFactory());
      aggFact->SetFactory("Graph", dropFact);
      aggFact->SetParameter("aggregation: coupling", Teuchos::ParameterEntry(coupling));
      aggFact->SetParameter("aggregation: mesh layout", Teuchos::ParameterEntry(meshLayout));
      aggFact->SetParameter("aggregation: number of spatial dimensions",
                            Teuchos::ParameterEntry(numDimensions));
      aggFact->SetParameter("aggregation: coarsening order", Teuchos::ParameterEntry(0));
      aggFact->SetParameter("aggregation: coarsening rate",
                            Teuchos::ParameterEntry(std::string("{3}")));

      level.Request("Aggregates", aggFact.get());

      level.Request(*aggFact);
      aggFact->Build(level);

      RCP<Aggregates> aggregates = level.Get<RCP<Aggregates> >("Aggregates",aggFact.get()); // fix me
      level.Release("Aggregates", aggFact.get());

      return aggregates;
    }  // gimmeStructuredAggregates

    // Little utility to generate uncoupled aggregates with some specified root nodes on interface.
    static RCP<Aggregates>
    gimmeInterfaceAggregates(const RCP<Matrix>&A, RCP<AmalgamationInfo> & amalgInfo, Teuchos::Array<LO>& nodeOnInterface)
    {
      Level level;
      TestHelpers::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(level);
      level.Set("A", A);

      RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
      RCP<CoalesceDropFactory> dropFact = rcp(new CoalesceDropFactory());
      dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

      // Setup aggregation factory (use default factory for graph)
      RCP<UncoupledAggregationFactory> aggFact = rcp(new UncoupledAggregationFactory());
      aggFact->SetFactory("Graph", dropFact);
      aggFact->SetParameter("aggregation: use interface aggregation",Teuchos::ParameterEntry(true));
      aggFact->SetParameter("aggregation: max agg size",Teuchos::ParameterEntry(3));
      aggFact->SetParameter("aggregation: min agg size",Teuchos::ParameterEntry(3));
      aggFact->SetParameter("aggregation: max selected neighbors",Teuchos::ParameterEntry(0));
      aggFact->SetParameter("aggregation: ordering",Teuchos::ParameterEntry(std::string("natural")));
      aggFact->SetParameter("aggregation: enable phase 1",  Teuchos::ParameterEntry(true));
      aggFact->SetParameter("aggregation: enable phase 2a", Teuchos::ParameterEntry(true));
      aggFact->SetParameter("aggregation: enable phase 2b", Teuchos::ParameterEntry(true));
      aggFact->SetParameter("aggregation: enable phase 3",  Teuchos::ParameterEntry(true));

      level.Set("nodeOnInterface",nodeOnInterface);



      level.Request("Aggregates", aggFact.get());
      level.Request("UnAmalgamationInfo", amalgFact.get());

      level.Request(*aggFact);
      aggFact->Build(level);

      RCP<Aggregates> aggregates = level.Get<RCP<Aggregates> >("Aggregates",aggFact.get());
      level.Release("Aggregates", aggFact.get());
      return aggregates;
    }  // gimmeInterfaceAggregates

    // Little utility to generate hybrid structured and uncoupled aggregates
    static RCP<Aggregates>
    gimmeHybridAggregates(const RCP<Matrix>&A, RCP<AmalgamationInfo> & amalgInfo,
                          const std::string regionType,
                          Array<GO> gNodesPerDir, Array<LO> lNodesPerDir,
                          const LO numDimensions, const std::string meshLayout,
                          const Array<GO> meshData)
    {
      Level level;
      TestHelpers::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(level);
      level.Set("A", A);
      level.Set("gNodesPerDim", gNodesPerDir);
      level.Set("lNodesPerDim", lNodesPerDir);
      level.Set("aggregation: mesh data", meshData);

      const std::string coupling = "uncoupled";

      RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
      amalgFact->SetDefaultVerbLevel(MueLu::None);
      RCP<CoalesceDropFactory> dropFact = rcp(new CoalesceDropFactory());
      dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

      // Setup aggregation factory (use default factory for graph)
      RCP<HybridAggregationFactory> aggFact = rcp(new HybridAggregationFactory());
      aggFact->SetFactory("Graph", dropFact);
      // Structured
      aggFact->SetParameter("aggregation: coupling",                     Teuchos::ParameterEntry(coupling));
      aggFact->SetParameter("aggregation: mesh layout",                  Teuchos::ParameterEntry(meshLayout));
      aggFact->SetParameter("aggregation: number of spatial dimensions", Teuchos::ParameterEntry(numDimensions));
      aggFact->SetParameter("aggregation: coarsening order",             Teuchos::ParameterEntry(0));
      aggFact->SetParameter("aggregation: coarsening rate",              Teuchos::ParameterEntry(std::string("{3}")));
      // Uncoupled
      aggFact->SetParameter("aggregation: use interface aggregation",    Teuchos::ParameterEntry(false));
      aggFact->SetParameter("aggregation: max agg size",                 Teuchos::ParameterEntry(3));
      aggFact->SetParameter("aggregation: min agg size",                 Teuchos::ParameterEntry(3));
      aggFact->SetParameter("aggregation: max selected neighbors",       Teuchos::ParameterEntry(0));
      aggFact->SetParameter("aggregation: ordering",                     Teuchos::ParameterEntry(std::string("natural")));
      aggFact->SetParameter("aggregation: enable phase 1",               Teuchos::ParameterEntry(true));
      aggFact->SetParameter("aggregation: enable phase 2a",              Teuchos::ParameterEntry(true));
      aggFact->SetParameter("aggregation: enable phase 2b",              Teuchos::ParameterEntry(true));
      aggFact->SetParameter("aggregation: enable phase 3",               Teuchos::ParameterEntry(true));

      // Hybrid
      level.Set("aggregationRegionType", regionType);


      level.Request("Aggregates", aggFact.get());
      level.Request("UnAmalgamationInfo", amalgFact.get());

      level.Request(*aggFact);
      aggFact->Build(level);

      RCP<Aggregates> aggregates = level.Get<RCP<Aggregates> >("Aggregates",aggFact.get());
      level.Release("Aggregates", aggFact.get());
      return aggregates;
    }  // gimmeHybridAggregates
};

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates, JustAggregation, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);

    out << "version: " << MueLu::Version() << std::endl;
    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(15);
    RCP<AmalgamationInfo> amalgInfo;
    RCP<Aggregates> aggregates = AggregateGenerator<SC,LO,GO,NO>::gimmeCoupledAggregates(A, amalgInfo);
    TEST_EQUALITY(aggregates != Teuchos::null, true);
    TEST_EQUALITY(aggregates->AggregatesCrossProcessors(),true);
  }

  ///////////////////////////////////////////////////////////////////////////

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates, GetNumAggregates, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;

    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
    RCP<const Map> rowmap = A->getRowMap();
    RCP<AmalgamationInfo> amalgInfo;
    RCP<Aggregates> aggregates = AggregateGenerator<SC,LO,GO,NO>::gimmeCoupledAggregates(A, amalgInfo);
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

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates, JustUncoupledAggregation, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;
    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(15);
    RCP<AmalgamationInfo> amalgInfo;
    RCP<Aggregates> aggregates = AggregateGenerator<SC,LO,GO,NO>::gimmeUncoupledAggregates(A, amalgInfo);
    TEST_EQUALITY(aggregates != Teuchos::null, true);
    TEST_EQUALITY(aggregates->AggregatesCrossProcessors(),false);
  }

  ///////////////////////////////////////////////////////////////////////////

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates, GetNumUncoupledAggregates, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;

    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
    RCP<const Map> rowmap = A->getRowMap();
    RCP<AmalgamationInfo> amalgInfo;
    RCP<Aggregates> aggregates = AggregateGenerator<SC,LO,GO,NO>::gimmeUncoupledAggregates(A, amalgInfo);
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

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates, UncoupledPhase1, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;

    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
    RCP<const Map> rowmap = A->getRowMap();
    RCP<AmalgamationInfo> amalgInfo;
    RCP<Aggregates> aggregates = AggregateGenerator<SC,LO,GO,NO>::gimmeUncoupledAggregates(A, amalgInfo,true,false,false,false);
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

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates, UncoupledPhase2, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;

    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
    RCP<const Map> rowmap = A->getRowMap();
    RCP<AmalgamationInfo> amalgInfo;
    bool bSuccess = true;
    TEUCHOS_TEST_THROW((AggregateGenerator<SC,LO,GO,NO>::gimmeUncoupledAggregates(A, amalgInfo,false,true,true,false)),
                        MueLu::Exceptions::RuntimeError, out, bSuccess);
    TEST_EQUALITY(bSuccess, true);
  } //UncoupledPhase2

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates, UncoupledPhase3, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;

    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
    RCP<const Map> rowmap = A->getRowMap();
    RCP<AmalgamationInfo> amalgInfo;
    RCP<Aggregates> aggregates = AggregateGenerator<SC,LO,GO,NO>::gimmeUncoupledAggregates(A, amalgInfo,false,false,false,true); GO numAggs = aggregates->GetNumAggregates();
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

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates, UncoupledInterface, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;

    // Get MPI parameters
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    LO numRanks = comm->getSize();
    LO myRank   = comm->getRank();

    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);


    RCP<const Map> rowmap = A->getRowMap();
    rowmap->describe(out,Teuchos::VERB_EXTREME);
    out<< "getNodeNumElements() gives: " << rowmap->getNodeNumElements() << std::endl;

    // Specify root nodes on interface
    Teuchos::Array<LO> nodeOnInterface(rowmap->getNodeNumElements(),0);
    if( rowmap->getMinAllGlobalIndex() == rowmap->getMinGlobalIndex() )
      nodeOnInterface[0] = 1;
    if( rowmap->getMaxAllGlobalIndex() == rowmap->getMaxGlobalIndex() )
      nodeOnInterface[rowmap->getNodeNumElements()-1] = 1;

    RCP<AmalgamationInfo> amalgInfo;
    RCP<Aggregates> aggregates = AggregateGenerator<SC,LO,GO,NO>::gimmeInterfaceAggregates(A, amalgInfo,nodeOnInterface);
    GO numAggs = aggregates->GetNumAggregates();


    // Check to see if specified nodes are root nodes
    for(LO i=0; i<nodeOnInterface.size(); i++){
      if (nodeOnInterface[i])
        TEST_EQUALITY(aggregates->IsRoot( i ), true);
    }

  } //UncoupledInterface

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates, JustStructuredAggregationGlobal, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;

    typedef typename Xpetra::MultiVector<double, LO, GO, NO> xdMV;

    // Get MPI parameters
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    LO numRanks = comm->getSize();
    LO myRank   = comm->getRank();

    // Set global geometric data
    const bool coupled = true;
    LO numDimensions = 3;
    Array<LO> lNodesPerDir(3);
    Array<GO> gNodesPerDir(3);
    Array<GO> meshData;
    const std::string meshLayout = "Global Lexicographic";
    for(int dim = 0; dim < 3; ++dim) {
      if(dim < numDimensions) {
        // Use more nodes in 1D to have a reasonable number of nodes per procs
        gNodesPerDir[dim] = ( (numDimensions == 1) ? 15 : 6 );
      } else {
        gNodesPerDir[dim] = 1;
      }
    }

    if(myRank == 0) std::cout << "About to build the coordinates" << std::endl;

    RCP<const Xpetra::MultiVector<double,LO,GO,NO> > Coordinates =
      TestHelpers::TestFactory<SC,LO,GO,NO>::BuildGeoCoordinates(numDimensions, gNodesPerDir,
                                                                 lNodesPerDir, meshData);

    Teuchos::ParameterList matrixList;
    matrixList.set("nx", gNodesPerDir[0]);
    matrixList.set("matrixType","Laplace1D");
    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr = Galeri::Xpetra::
      BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>("Laplace1D", Coordinates->getMap(),
                                                           matrixList);
    RCP<Matrix> A = Pr->BuildMatrix();

    if(myRank == 0) std::cout << "About to build the aggregates" << std::endl;

    RCP<Aggregates> aggregates = AggregateGenerator<SC,LO,GO,NO>::
      gimmeStructuredAggregates(A, gNodesPerDir, lNodesPerDir, coupled, numDimensions, meshLayout,
                                meshData);

    TEST_EQUALITY(aggregates != Teuchos::null, true);
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates, JustStructuredAggregationLocal, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;

    typedef typename Xpetra::MultiVector<double, LO, GO, NO> xdMV;

    // Get MPI parameter
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    LO numRanks = comm->getSize();
    LO myRank   = comm->getRank();

    // Set global geometric data
    const bool coupled = true;
    const std::string meshLayout = std::string("Local Lexicographic");
    LO numDimensions = 2;
    Array<LO> lNodesPerDir(3);
    Array<GO> gNodesPerDir(3);
    Array<GO> meshData(10*numRanks);
    for(int dim = 0; dim < 3; ++dim) {
      if(dim < numDimensions) {
        // Use more nodes in 1D to have a reasonable number of nodes per procs
        gNodesPerDir[dim] = ( (numDimensions == 1) ? 15 : 6 );
      } else {
        gNodesPerDir[dim] = 1;
      }
    }

    RCP<const Xpetra::MultiVector<double,LO,GO,NO> > Coordinates =
      TestHelpers::TestFactory<SC,LO,GO,NO>::BuildGeoCoordinates(numDimensions, gNodesPerDir,
                                                                 lNodesPerDir, meshData,
                                                                 meshLayout);

    Teuchos::ParameterList matrixList;
    matrixList.set("nx", gNodesPerDir[0]);
    matrixList.set("matrixType","Laplace1D");
    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr = Galeri::Xpetra::
      BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>("Laplace1D", Coordinates->getMap(),
                                                           matrixList);
    RCP<Matrix> A = Pr->BuildMatrix();

    RCP<Aggregates> aggregates = AggregateGenerator<SC,LO,GO,NO>::
      gimmeStructuredAggregates(A, gNodesPerDir, lNodesPerDir, coupled, numDimensions, meshLayout,
                                meshData);

    TEST_EQUALITY(aggregates != Teuchos::null, true);
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates, HybridAggregation, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;

    // Get MPI parameters
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    LO numRanks = comm->getSize();
    LO myRank   = comm->getRank();

    // Set hybrid region type
    std::string regionType;
    if( myRank % 2 == 0 ){
      regionType = "structured";
    } else {
      regionType = "uncoupled";
    }

    // Set global geometric data
    const std::string meshLayout = std::string("Global Lexicographic");
    LO numDimensions = 2;
    Array<LO> lNodesPerDir(3);
    Array<GO> gNodesPerDir(3);
    Array<GO> meshData(10*numRanks);
    for(int dim = 0; dim < 3; ++dim) {
      if(dim < numDimensions) {
        // Use more nodes in 1D to have a reasonable number of nodes per procs
        gNodesPerDir[dim] = ( (numDimensions == 1) ? 15 : 6 );
      } else {
        gNodesPerDir[dim] = 1;
      }
    }

    RCP<const Xpetra::MultiVector<double,LO,GO,NO> > Coordinates =
      TestHelpers::TestFactory<SC,LO,GO,NO>::BuildGeoCoordinates(numDimensions, gNodesPerDir,
                                                                 lNodesPerDir, meshData,
                                                                 meshLayout);
    // Build Problem
    Teuchos::ParameterList matrixList;
    matrixList.set("nx", gNodesPerDir[0]);
    matrixList.set("matrixType","Laplace1D");
    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr = Galeri::Xpetra::
      BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>("Laplace1D", Coordinates->getMap(),
                                                           matrixList);
    RCP<Matrix> A = Pr->BuildMatrix();

    //aggregates
    RCP<AmalgamationInfo> amalgInfo;
    RCP<Aggregates> aggregates = AggregateGenerator<SC,LO,GO,NO>::
      gimmeHybridAggregates(A, amalgInfo, regionType,
                            gNodesPerDir, lNodesPerDir,
                            numDimensions, meshLayout,
                            meshData);


    TEST_EQUALITY(aggregates != Teuchos::null, true);
    TEST_EQUALITY(aggregates->AggregatesCrossProcessors(),false);
  }
#define MUELU_ETI_GROUP(Scalar,LO,GO,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates,JustAggregation,Scalar,LO,GO,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates,GetNumAggregates,Scalar,LO,GO,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates,JustUncoupledAggregation,Scalar,LO,GO,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates,GetNumUncoupledAggregates,Scalar,LO,GO,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates,UncoupledPhase1,Scalar,LO,GO,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates,UncoupledPhase2,Scalar,LO,GO,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates,UncoupledPhase3,Scalar,LO,GO,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates,UncoupledInterface,Scalar,LO,GO,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates,JustStructuredAggregationGlobal,Scalar,LO,GO,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates,JustStructuredAggregationLocal,Scalar,LO,GO,Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates,HybridAggregation,Scalar,LO,GO,Node)

#include <MueLu_ETI_4arg.hpp>


} // namespace MueLuTests
