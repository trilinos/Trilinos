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
#include <MueLu_UncoupledAggregationFactory.hpp>
#include <MueLu_StructuredAggregationFactory.hpp>
#include <MueLu_UncoupledAggregationFactory.hpp>
#include <MueLu_HybridAggregationFactory.hpp>
#include <MueLu_AmalgamationFactory.hpp>
#include <MueLu_AmalgamationInfo.hpp>
#include <MueLu_Aggregates.hpp>
#include <MueLu_FilteredAFactory.hpp>

namespace MueLuTests {

template <class T>
void print_matrix(const char* label, const Teuchos::RCP<T>& Mat) {
  ArrayRCP<const size_t> rowptr;
  ArrayRCP<const typename T::local_ordinal_type> colind;
  ArrayRCP<const typename T::scalar_type> vals;
  Mat->getAllValues(rowptr, colind, vals);

  printf("*** %s ***\n", label);
  printf("rowptr = ");
  for (int i = 0; i < (int)rowptr.size(); i++)
    printf("%d ", (int)rowptr[i]);
  printf("\ncolind =");
  for (int i = 0; i < (int)colind.size(); i++)
    printf("%d ", (int)colind[i]);
  printf("\nvalues =");
  for (int i = 0; i < (int)vals.size(); i++)
    printf("%8.1e ", vals[i]);
  printf("\n");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class AggregateGenerator {
  // typedef MueLu::AmalgamationInfo<LocalOrdinal,GlobalOrdinal,Node> amalgamation_info_type;
  // typedef MueLu::Aggregates<LocalOrdinal,GlobalOrdinal,Node> aggregates_type;
  // typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> xpetra_matrix_type;
#include <MueLu_UseShortNames.hpp>

 public:
  // Little utility to generate uncoupled aggregates.
  static RCP<Aggregates>
  gimmeUncoupledAggregates(const RCP<Matrix>& A, RCP<AmalgamationInfo>& amalgInfo, bool bPhase1 = true, bool bPhase2a = true, bool bPhase2b = true, bool bPhase3 = true)
  //  RCP<MueLu::Aggregates> gimmeUncoupledAggregates(const RCP<xpetra_matrix_type> & A, RCP<AmalgamationInfo> & amalgInfo, bool bPhase1 = true, bool bPhase2a = true, bool bPhase2b = true, bool bPhase3 = true)
  {
    Level level;
    TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
    RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
    dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

    // Setup aggregation factory (use default factory for graph)
    RCP<UncoupledAggregationFactory> aggFact = rcp(new UncoupledAggregationFactory());
    aggFact->SetFactory("Graph", dropFact);
    aggFact->SetParameter("aggregation: max agg size", Teuchos::ParameterEntry(3));
    aggFact->SetParameter("aggregation: min agg size", Teuchos::ParameterEntry(3));
    aggFact->SetParameter("aggregation: max selected neighbors", Teuchos::ParameterEntry(0));
    aggFact->SetParameter("aggregation: ordering", Teuchos::ParameterEntry(std::string("natural")));
    aggFact->SetParameter("aggregation: allow user-specified singletons", Teuchos::ParameterEntry(true));

    aggFact->SetParameter("aggregation: enable phase 1", Teuchos::ParameterEntry(bPhase1));
    aggFact->SetParameter("aggregation: enable phase 2a", Teuchos::ParameterEntry(bPhase2a));
    aggFact->SetParameter("aggregation: enable phase 2b", Teuchos::ParameterEntry(bPhase2b));
    aggFact->SetParameter("aggregation: enable phase 3", Teuchos::ParameterEntry(bPhase3));
    aggFact->SetParameter("aggregation: phase3 avoid singletons", Teuchos::ParameterEntry(true));
    aggFact->SetParameter("aggregation: use interface aggregation", Teuchos::ParameterEntry(false));
    level.Request("Aggregates", aggFact.get());
    level.Request("UnAmalgamationInfo", amalgFact.get());

    level.Request(*aggFact);
    aggFact->Build(level);
    RCP<Aggregates> aggregates = level.Get<RCP<Aggregates>>("Aggregates", aggFact.get());                  // fix me
    amalgInfo                  = level.Get<RCP<AmalgamationInfo>>("UnAmalgamationInfo", amalgFact.get());  // fix me
    level.Release("UnAmalgamationInfo", amalgFact.get());
    level.Release("Aggregates", aggFact.get());
    return aggregates;
  }  // gimmeUncoupledAggregates

  // Little utility to generate uncoupled aggregates.
  static RCP<Aggregates>
  gimmeStructuredAggregates(const RCP<Matrix>& A,
                            Array<GO> gNodesPerDir, Array<LO> lNodesPerDir, const bool coupled,
                            const LO numDimensions, const std::string meshLayout,
                            const Array<GO> meshData) {
    Level level;
    TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);

    const std::string coupling = (coupled ? "coupled" : "uncoupled");

    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
    amalgFact->SetDefaultVerbLevel(MueLu::None);
    RCP<CoalesceDropFactory> dropFact = rcp(new CoalesceDropFactory());
    dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

    level.Set("A", A);
    level.Set("numDimensions", numDimensions);
    level.Set("gNodesPerDim", gNodesPerDir);
    level.Set("lNodesPerDim", lNodesPerDir);
    level.Set("aggregation: mesh data", meshData);

    // Setup aggregation factory (use default factory for graph)
    RCP<StructuredAggregationFactory> aggFact = rcp(new StructuredAggregationFactory());
    aggFact->SetFactory("Graph", dropFact);
    aggFact->SetParameter("aggregation: mode", Teuchos::ParameterEntry(coupling));
    aggFact->SetParameter("aggregation: mesh layout", Teuchos::ParameterEntry(meshLayout));
    aggFact->SetParameter("aggregation: coarsening order", Teuchos::ParameterEntry(0));
    aggFact->SetParameter("aggregation: coarsening rate",
                          Teuchos::ParameterEntry(std::string("{3}")));

    level.Request("Aggregates", aggFact.get());

    level.Request(*aggFact);
    aggFact->Build(level);

    RCP<Aggregates> aggregates = level.Get<RCP<Aggregates>>("Aggregates", aggFact.get());  // fix me
    level.Release("Aggregates", aggFact.get());

    return aggregates;
  }  // gimmeStructuredAggregates

  // Little utility to generate uncoupled aggregates with some specified root nodes on interface.
  static RCP<Aggregates>
  gimmeInterfaceAggregates(const RCP<Matrix>& A, RCP<AmalgamationInfo>& amalgInfo, Teuchos::Array<LO>& nodeOnInterface) {
    Level level;
    TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
    RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
    dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

    // Setup aggregation factory (use default factory for graph)
    RCP<UncoupledAggregationFactory> aggFact = rcp(new UncoupledAggregationFactory());
    aggFact->SetFactory("Graph", dropFact);
    aggFact->SetParameter("aggregation: use interface aggregation", Teuchos::ParameterEntry(true));
    aggFact->SetParameter("aggregation: max agg size", Teuchos::ParameterEntry(3));
    aggFact->SetParameter("aggregation: min agg size", Teuchos::ParameterEntry(3));
    aggFact->SetParameter("aggregation: max selected neighbors", Teuchos::ParameterEntry(0));
    aggFact->SetParameter("aggregation: ordering", Teuchos::ParameterEntry(std::string("natural")));
    aggFact->SetParameter("aggregation: enable phase 1", Teuchos::ParameterEntry(true));
    aggFact->SetParameter("aggregation: enable phase 2a", Teuchos::ParameterEntry(true));
    aggFact->SetParameter("aggregation: enable phase 2b", Teuchos::ParameterEntry(true));
    aggFact->SetParameter("aggregation: enable phase 3", Teuchos::ParameterEntry(true));

    level.Set("nodeOnInterface", nodeOnInterface);

    level.Request("Aggregates", aggFact.get());
    level.Request("UnAmalgamationInfo", amalgFact.get());

    level.Request(*aggFact);
    aggFact->Build(level);

    RCP<Aggregates> aggregates = level.Get<RCP<Aggregates>>("Aggregates", aggFact.get());
    level.Release("Aggregates", aggFact.get());
    return aggregates;
  }  // gimmeInterfaceAggregates

  // Little utility to generate hybrid structured and uncoupled aggregates
  static RCP<Aggregates>
  gimmeHybridAggregates(const RCP<Matrix>& A, RCP<AmalgamationInfo>& amalgInfo,
                        const std::string regionType,
                        Array<GO> gNodesPerDir, Array<LO> lNodesPerDir,
                        const LO numDimensions) {
    Level level;
    TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
    level.Set("A", A);
    level.Set("numDimensions", numDimensions);
    level.Set("lNodesPerDim", lNodesPerDir);

    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
    amalgFact->SetDefaultVerbLevel(MueLu::None);
    RCP<CoalesceDropFactory> dropFact = rcp(new CoalesceDropFactory());
    dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

    // Setup aggregation factory (use default factory for graph)
    RCP<HybridAggregationFactory> aggFact = rcp(new HybridAggregationFactory());
    aggFact->SetFactory("Graph", dropFact);
    // Structured
    aggFact->SetParameter("aggregation: coarsening order", Teuchos::ParameterEntry(0));
    aggFact->SetParameter("aggregation: coarsening rate", Teuchos::ParameterEntry(std::string("{3}")));
    // Uncoupled
    aggFact->SetParameter("aggregation: use interface aggregation", Teuchos::ParameterEntry(false));
    aggFact->SetParameter("aggregation: max agg size", Teuchos::ParameterEntry(3));
    aggFact->SetParameter("aggregation: min agg size", Teuchos::ParameterEntry(3));
    aggFact->SetParameter("aggregation: max selected neighbors", Teuchos::ParameterEntry(0));
    aggFact->SetParameter("aggregation: ordering", Teuchos::ParameterEntry(std::string("natural")));
    aggFact->SetParameter("aggregation: enable phase 1", Teuchos::ParameterEntry(true));
    aggFact->SetParameter("aggregation: enable phase 2a", Teuchos::ParameterEntry(true));
    aggFact->SetParameter("aggregation: enable phase 2b", Teuchos::ParameterEntry(true));
    aggFact->SetParameter("aggregation: enable phase 3", Teuchos::ParameterEntry(true));

    aggFact->SetParameter("aggregation: match ML phase2a", Teuchos::ParameterEntry(true));

    // Hybrid
    level.Set("aggregationRegionType", regionType);

    level.Request("Aggregates", aggFact.get());
    level.Request("UnAmalgamationInfo", amalgFact.get());

    level.Request(*aggFact);
    aggFact->Build(level);

    RCP<Aggregates> aggregates = level.Get<RCP<Aggregates>>("Aggregates", aggFact.get());
    level.Release("Aggregates", aggFact.get());
    return aggregates;
  }  // gimmeHybridAggregates
};

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates, JustAggregation, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  out << "version: " << MueLu::Version() << std::endl;
  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(15);
  RCP<AmalgamationInfo> amalgInfo;
  RCP<Aggregates> aggregates = AggregateGenerator<SC, LO, GO, NO>::gimmeUncoupledAggregates(A, amalgInfo);
  TEST_EQUALITY(aggregates != Teuchos::null, true);
  TEST_EQUALITY(aggregates->AggregatesCrossProcessors(), false);
}

///////////////////////////////////////////////////////////////////////////

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates, GetNumAggregates, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<Matrix> A         = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  RCP<const Map> rowmap = A->getRowMap();
  RCP<AmalgamationInfo> amalgInfo;
  RCP<Aggregates> aggregates         = AggregateGenerator<SC, LO, GO, NO>::gimmeUncoupledAggregates(A, amalgInfo);
  GO numAggs                         = aggregates->GetNumAggregates();
  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();

  TEST_EQUALITY(aggregates->AggregatesCrossProcessors(), false);

  ArrayRCP<LO> aggSizes = Teuchos::ArrayRCP<LO>(numAggs);
  ArrayRCP<LO> aggStart;
  ArrayRCP<GO> aggToRowMap;
  amalgInfo->UnamalgamateAggregates(*aggregates, aggStart, aggToRowMap);
  for (LO i = 0; i < numAggs; ++i)
    aggSizes[i] = aggStart[i + 1] - aggStart[i];

  bool foundAggNotSize3 = false;
  for (int i = 0; i < aggSizes.size(); ++i)
    if (aggSizes[i] != 3) {
      foundAggNotSize3 = true;
      break;
    }

  switch (comm->getSize()) {
    case 1:
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
      // throw(MueLu::Exceptions::NotImplemented(msg));
      out << msg << std::endl;
      break;
  }

  // ArrayRCP< ArrayRCP<GO> > aggToRowMap(numAggs);
  int root = out.getOutputToRootOnly();
  out.setOutputToRootOnly(-1);
  for (int j = 0; j < comm->getSize(); ++j) {
    if (comm->getRank() == j) {
      out << "++ pid " << j << " ++" << std::endl;
      out << "   num local DOFs = " << rowmap->getLocalNumElements() << std::endl;
      for (int i = 0; i < numAggs; ++i) {
        out << "   aggregate " << i << ": ";
        for (int k = aggStart[i]; k < aggStart[i + 1]; ++k)
          out << aggToRowMap[k] << " ";
        out << std::endl;
      }
    }
    comm->barrier();
  }
  out.setOutputToRootOnly(root);

}  // GetNumAggregates

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates, JustUncoupledAggregation, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(15);
  RCP<AmalgamationInfo> amalgInfo;
  RCP<Aggregates> aggregates = AggregateGenerator<SC, LO, GO, NO>::gimmeUncoupledAggregates(A, amalgInfo);
  TEST_EQUALITY(aggregates != Teuchos::null, true);
  TEST_EQUALITY(aggregates->AggregatesCrossProcessors(), false);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates, JustUncoupledAggregationFactory, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  // Setup aggregation factory (use default factory for graph)
  RCP<UncoupledAggregationFactory> aggFact = rcp(new UncoupledAggregationFactory);
  aggFact->SetOrdering("graph");
  TEST_EQUALITY(aggFact->GetOrdering() == "graph", true);
  aggFact->SetOrdering("natural");
  TEST_EQUALITY(aggFact->GetOrdering() == "natural", true);
  aggFact->SetOrdering("random");
  TEST_EQUALITY(aggFact->GetOrdering() == "random", true);

  aggFact->SetMaxNeighAlreadySelected(12);
  TEST_EQUALITY(aggFact->GetMaxNeighAlreadySelected() == 12, true);

  aggFact->SetMaxNeighAlreadySelected(0);
  TEST_EQUALITY(aggFact->GetMaxNeighAlreadySelected() == 0, true);

  aggFact->SetMinNodesPerAggregate(0);
  TEST_EQUALITY(aggFact->GetMinNodesPerAggregate() == 0, true);

  aggFact->SetMinNodesPerAggregate(3);
  TEST_EQUALITY(aggFact->GetMinNodesPerAggregate() == 3, true);
}
///////////////////////////////////////////////////////////////////////////

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates, GetNumUncoupledAggregates, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<Matrix> A         = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  RCP<const Map> rowmap = A->getRowMap();
  RCP<AmalgamationInfo> amalgInfo;
  RCP<Aggregates> aggregates         = AggregateGenerator<SC, LO, GO, NO>::gimmeUncoupledAggregates(A, amalgInfo);
  GO numAggs                         = aggregates->GetNumAggregates();
  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();

  TEST_EQUALITY(aggregates->AggregatesCrossProcessors(), false);

  ArrayRCP<LO> aggSizes = Teuchos::ArrayRCP<LO>(numAggs);
  ArrayRCP<LO> aggStart;
  ArrayRCP<GO> aggToRowMap;
  amalgInfo->UnamalgamateAggregates(*aggregates, aggStart, aggToRowMap);
  for (LO i = 0; i < numAggs; ++i)
    aggSizes[i] = aggStart[i + 1] - aggStart[i];

  bool foundAggNotSize3 = false;
  for (int i = 0; i < aggSizes.size(); ++i)
    if (aggSizes[i] != 3) {
      foundAggNotSize3 = true;
      break;
    }

  switch (comm->getSize()) {
    case 1:
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
      // throw(MueLu::Exceptions::NotImplemented(msg));
      out << msg << std::endl;
      break;
  }

  // ArrayRCP< ArrayRCP<GO> > aggToRowMap(numAggs);
  int root = out.getOutputToRootOnly();
  out.setOutputToRootOnly(-1);
  for (int j = 0; j < comm->getSize(); ++j) {
    if (comm->getRank() == j) {
      out << "++ pid " << j << " ++" << std::endl;
      out << "   num local DOFs = " << rowmap->getLocalNumElements() << std::endl;
      for (int i = 0; i < numAggs; ++i) {
        out << "   aggregate " << i << ": ";
        for (int k = aggStart[i]; k < aggStart[i + 1]; ++k)
          out << aggToRowMap[k] << " ";
        out << std::endl;
      }
    }
    comm->barrier();
  }
  out.setOutputToRootOnly(root);

}  // GetNumAggregates

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates, UncoupledPhase1, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<Matrix> A         = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  RCP<const Map> rowmap = A->getRowMap();
  RCP<AmalgamationInfo> amalgInfo;
  RCP<Aggregates> aggregates         = AggregateGenerator<SC, LO, GO, NO>::gimmeUncoupledAggregates(A, amalgInfo, true, false, false, false);
  GO numAggs                         = aggregates->GetNumAggregates();
  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();

  TEST_EQUALITY(aggregates->AggregatesCrossProcessors(), false);

  ArrayRCP<LO> aggSizes = Teuchos::ArrayRCP<LO>(numAggs);
  ArrayRCP<LO> aggStart;
  ArrayRCP<GO> aggToRowMap;
  amalgInfo->UnamalgamateAggregates(*aggregates, aggStart, aggToRowMap);
  for (LO i = 0; i < numAggs; ++i)
    aggSizes[i] = aggStart[i + 1] - aggStart[i];

  bool foundAggNotSize3 = false;
  for (int i = 0; i < aggSizes.size(); ++i)
    if (aggSizes[i] != 3) {
      foundAggNotSize3 = true;
      break;
    }

  switch (comm->getSize()) {
    case 1:
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
      // throw(MueLu::Exceptions::NotImplemented(msg));
      out << msg << std::endl;
      break;
  }

  // ArrayRCP< ArrayRCP<GO> > aggToRowMap(numAggs);
  int root = out.getOutputToRootOnly();
  out.setOutputToRootOnly(-1);
  for (int j = 0; j < comm->getSize(); ++j) {
    if (comm->getRank() == j) {
      out << "++ pid " << j << " ++" << std::endl;
      out << "   num local DOFs = " << rowmap->getLocalNumElements() << std::endl;
      for (int i = 0; i < numAggs; ++i) {
        out << "   aggregate " << i << ": ";
        for (int k = aggStart[i]; k < aggStart[i + 1]; ++k)
          out << aggToRowMap[k] << " ";
        out << std::endl;
      }
    }
    comm->barrier();
  }
  out.setOutputToRootOnly(root);

}  // UncoupledPhase1

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates, UncoupledPhase2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<Matrix> A         = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  RCP<const Map> rowmap = A->getRowMap();
  RCP<AmalgamationInfo> amalgInfo;
  bool bSuccess = true;
  TEUCHOS_TEST_THROW((AggregateGenerator<SC, LO, GO, NO>::gimmeUncoupledAggregates(A, amalgInfo, false, true, true, false)),
                     MueLu::Exceptions::RuntimeError, out, bSuccess);
  TEST_EQUALITY(bSuccess, true);
}  // UncoupledPhase2

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates, UncoupledPhase3, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<Matrix> A         = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  RCP<const Map> rowmap = A->getRowMap();
  RCP<AmalgamationInfo> amalgInfo;
  RCP<Aggregates> aggregates         = AggregateGenerator<SC, LO, GO, NO>::gimmeUncoupledAggregates(A, amalgInfo, false, false, false, true);
  GO numAggs                         = aggregates->GetNumAggregates();
  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();

  TEST_EQUALITY(aggregates->AggregatesCrossProcessors(), false);

  ArrayRCP<LO> aggSizes = Teuchos::ArrayRCP<LO>(numAggs);
  ArrayRCP<LO> aggStart;
  ArrayRCP<GO> aggToRowMap;
  amalgInfo->UnamalgamateAggregates(*aggregates, aggStart, aggToRowMap);
  for (LO i = 0; i < numAggs; ++i)
    aggSizes[i] = aggStart[i + 1] - aggStart[i];

  // LBV on 09/27/19: the aggregate size of 2 is only
  // expected because the problem is 1D and the nodes
  // are treated in a lexicographic order. Using a
  // different ordering will almost always result in
  // at least on aggregate of size 3 and one of size 1.
  bool foundAggNotSize2 = false;
  for (int i = 0; i < aggSizes.size(); ++i)
    if (aggSizes[i] != 2) {
      foundAggNotSize2 = true;
      break;
    }

  switch (comm->getSize()) {
    case 1:
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
      // throw(MueLu::Exceptions::NotImplemented(msg));
      out << msg << std::endl;
      break;
  }

  // ArrayRCP< ArrayRCP<GO> > aggToRowMap(numAggs);
  int root = out.getOutputToRootOnly();
  out.setOutputToRootOnly(-1);
  for (int j = 0; j < comm->getSize(); ++j) {
    if (comm->getRank() == j) {
      out << "++ pid " << j << " ++" << std::endl;
      out << "   num local DOFs = " << rowmap->getLocalNumElements() << std::endl;
      for (int i = 0; i < numAggs; ++i) {
        out << "   aggregate " << i << ": ";
        for (int k = aggStart[i]; k < aggStart[i + 1]; ++k)
          out << aggToRowMap[k] << " ";
        out << std::endl;
      }
    }
    comm->barrier();
  }
  out.setOutputToRootOnly(root);

}  // UncoupledPhase3

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates, UncoupledInterface, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  // Get MPI parameters
  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);

  RCP<const Map> rowmap = A->getRowMap();
  rowmap->describe(out, Teuchos::VERB_EXTREME);
  out << "getLocalNumElements() gives: " << rowmap->getLocalNumElements() << std::endl;

  // Specify root nodes on interface
  Teuchos::Array<LO> nodeOnInterface(rowmap->getLocalNumElements(), 0);
  if (rowmap->getMinAllGlobalIndex() == rowmap->getMinGlobalIndex())
    nodeOnInterface[0] = 1;
  if (rowmap->getMaxAllGlobalIndex() == rowmap->getMaxGlobalIndex())
    nodeOnInterface[rowmap->getLocalNumElements() - 1] = 1;

  RCP<AmalgamationInfo> amalgInfo;
  RCP<Aggregates> aggregates = AggregateGenerator<SC, LO, GO, NO>::gimmeInterfaceAggregates(A, amalgInfo, nodeOnInterface);

  // Check to see if specified nodes are root nodes
  for (LO i = 0; i < nodeOnInterface.size(); i++) {
    if (nodeOnInterface[i])
      TEST_EQUALITY(aggregates->IsRoot(i), true);
  }

}  // UncoupledInterface

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates, JustStructuredAggregationGlobal, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  typedef typename Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> xdMV;

  // Get MPI parameters
  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();
  LO myRank                          = comm->getRank();

  // Set global geometric data
  const bool coupled = true;
  LO numDimensions   = 3;
  Array<LO> lNodesPerDir(3);
  Array<GO> gNodesPerDir(3);
  Array<GO> meshData;
  const std::string meshLayout = "Global Lexicographic";
  for (int dim = 0; dim < 3; ++dim) {
    if (dim < numDimensions) {
      // Use more nodes in 1D to have a reasonable number of nodes per procs
      gNodesPerDir[dim] = ((numDimensions == 1) ? 15 : 6);
    } else {
      gNodesPerDir[dim] = 1;
    }
  }

  if (myRank == 0) std::cout << "About to build the coordinates" << std::endl;

  RCP<const Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>> Coordinates =
      TestHelpers::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions, gNodesPerDir,
                                                                    lNodesPerDir, meshData);

  Teuchos::ParameterList matrixList;
  matrixList.set("nx", gNodesPerDir[0]);
  matrixList.set("matrixType", "Laplace1D");
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector>> Pr = Galeri::Xpetra::
      BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace1D", Coordinates->getMap(),
                                                                matrixList);
  RCP<Matrix> A = Pr->BuildMatrix();

  if (myRank == 0) std::cout << "About to build the aggregates" << std::endl;

  RCP<Aggregates> aggregates = AggregateGenerator<SC, LO, GO, NO>::
      gimmeStructuredAggregates(A, gNodesPerDir, lNodesPerDir, coupled, numDimensions, meshLayout,
                                meshData);

  TEST_EQUALITY(aggregates != Teuchos::null, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates, JustStructuredAggregationLocal, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  typedef typename Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> xdMV;

  // Get MPI parameter
  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();
  LO numRanks                        = comm->getSize();

  // Set global geometric data
  const bool coupled           = true;
  const std::string meshLayout = std::string("Local Lexicographic");
  LO numDimensions             = 2;
  Array<LO> lNodesPerDir(3);
  Array<GO> gNodesPerDir(3);
  Array<GO> meshData(10 * numRanks);
  for (int dim = 0; dim < 3; ++dim) {
    if (dim < numDimensions) {
      // Use more nodes in 1D to have a reasonable number of nodes per procs
      gNodesPerDir[dim] = ((numDimensions == 1) ? 15 : 6);
    } else {
      gNodesPerDir[dim] = 1;
    }
  }

  RCP<const Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>> Coordinates =
      TestHelpers::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions, gNodesPerDir,
                                                                    lNodesPerDir, meshData,
                                                                    meshLayout);

  Teuchos::ParameterList matrixList;
  matrixList.set("nx", gNodesPerDir[0]);
  matrixList.set("matrixType", "Laplace1D");
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector>> Pr = Galeri::Xpetra::
      BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace1D", Coordinates->getMap(),
                                                                matrixList);
  RCP<Matrix> A = Pr->BuildMatrix();

  RCP<Aggregates> aggregates = AggregateGenerator<SC, LO, GO, NO>::
      gimmeStructuredAggregates(A, gNodesPerDir, lNodesPerDir, coupled, numDimensions, meshLayout,
                                meshData);

  TEST_EQUALITY(aggregates != Teuchos::null, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates, HybridAggregation, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  // Get MPI parameters
  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();
  LO numRanks                        = comm->getSize();
  LO myRank                          = comm->getRank();

  // Set hybrid region type
  std::string regionType;
  if (myRank % 2 == 0) {
    regionType = "structured";
  } else {
    regionType = "uncoupled";
  }

  // Set global geometric data
  const std::string meshLayout = std::string("Global Lexicographic");
  LO numDimensions             = 2;
  Array<LO> lNodesPerDir(3);
  Array<GO> gNodesPerDir(3);
  Array<GO> meshData(10 * numRanks);
  for (int dim = 0; dim < 3; ++dim) {
    if (dim < numDimensions) {
      // Use more nodes in 1D to have a reasonable number of nodes per procs
      gNodesPerDir[dim] = ((numDimensions == 1) ? 15 : 6);
    } else {
      gNodesPerDir[dim] = 1;
    }
  }

  RCP<const Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>> Coordinates =
      TestHelpers::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions, gNodesPerDir,
                                                                    lNodesPerDir, meshData,
                                                                    meshLayout);
  // Build Problem
  Teuchos::ParameterList matrixList;
  matrixList.set("nx", gNodesPerDir[0]);
  matrixList.set("matrixType", "Laplace1D");
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector>> Pr = Galeri::Xpetra::
      BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace1D", Coordinates->getMap(),
                                                                matrixList);
  RCP<Matrix> A = Pr->BuildMatrix();

  // aggregates
  RCP<AmalgamationInfo> amalgInfo;
  RCP<Aggregates> aggregates = AggregateGenerator<SC, LO, GO, NO>::
      gimmeHybridAggregates(A, amalgInfo, regionType,
                            gNodesPerDir, lNodesPerDir,
                            numDimensions);

  TEST_EQUALITY(aggregates != Teuchos::null, true);
  TEST_EQUALITY(aggregates->AggregatesCrossProcessors(), false);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates, GreedyDirichlet, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  out << "running greedy Dirichlet" << std::endl;

  typedef typename Teuchos::ScalarTraits<Scalar> TST;

  //    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(30);
  // Make a Matrix with multiple degrees of freedom per node
  GlobalOrdinal nx = 20, ny = 20;

  // Describes the initial layout of matrix rows across processors.
  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  galeriList.set("ny", ny);
  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();
  RCP<const Map> map                 = Galeri::Xpetra::CreateMap<LocalOrdinal, GlobalOrdinal, Node>(TestHelpers::Parameters::getLib(), "Cartesian2D", comm, galeriList);

  map = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(map, 2);  // expand map for 2 DOFs per node

  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector>> Pr =
      Galeri::Xpetra::BuildProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, CrsMatrixWrap, MultiVector>("Elasticity2D", map, galeriList);
  RCP<Matrix> A = Pr->BuildMatrix();
  A->SetFixedBlockSize(2);

  Teuchos::ArrayView<const LocalOrdinal> indices;
  Teuchos::ArrayView<const Scalar> values;

  // Create a dirichlet boundary row.
  LocalOrdinal localRowToZero = 5;  // Corresponds to a Dof on local graph node 2

  A->resumeFill();
  A->getLocalRowView(localRowToZero, indices, values);
  Array<Scalar> newvalues(values.size(), TST::zero());
  for (int j = 0; j < indices.size(); j++)
    // keep diagonal
    if (indices[j] == localRowToZero) newvalues[j] = values[j];
  A->replaceLocalValues(localRowToZero, indices, newvalues);

  A->fillComplete();

  ArrayRCP<const bool> drows = Utilities::DetectDirichletRows(*A);
  TEST_EQUALITY(drows[localRowToZero], true);
  TEST_EQUALITY(drows[localRowToZero - 1], false);

  RCP<AmalgamationInfo> amalgInfo;
  Level level;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
  level.Set("A", A);

  RCP<CoalesceDropFactory> dropFact;
  RCP<AmalgamationFactory> amalgFact;
  RCP<UncoupledAggregationFactory> aggFact;

  amalgFact = rcp(new AmalgamationFactory());
  dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetParameter("aggregation: greedy Dirichlet", Teuchos::ParameterEntry(false));
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  // Setup aggregation factory (use default factory for graph)
  aggFact = rcp(new UncoupledAggregationFactory());
  aggFact->SetFactory("Graph", dropFact);

  level.Request("Aggregates", aggFact.get());
  level.Request("UnAmalgamationInfo", amalgFact.get());

  level.Request(*aggFact);
  aggFact->Build(level);
  RCP<Aggregates> aggregates = level.Get<RCP<Aggregates>>("Aggregates", aggFact.get());
  typename Aggregates::LO_view aggPtr;
  typename Aggregates::LO_view aggNodes;
  typename Aggregates::LO_view unaggregated;

  aggregates->ComputeNodesInAggregate(aggPtr, aggNodes, unaggregated);

  typename Aggregates::LO_view::HostMirror unaggregated_h = Kokkos::create_mirror_view(unaggregated);
  Kokkos::deep_copy(unaggregated_h, unaggregated);

  //     Test to check that the dirichlet node is aggregated:
  for (LO i = 0; i < (LO)unaggregated_h.extent(0); i++) {
    TEST_EQUALITY(unaggregated_h(i) == 2, false);
  }
  // Repeat with greedy Dirichlet
  Level levelGreedyAndNoPreserve;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(levelGreedyAndNoPreserve);
  levelGreedyAndNoPreserve.Set("A", A);
  amalgFact = rcp(new AmalgamationFactory());
  dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetParameter("aggregation: greedy Dirichlet", Teuchos::ParameterEntry(true));
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  // Setup aggregation factory (use default factory for graph)
  aggFact = rcp(new UncoupledAggregationFactory());
  aggFact->SetFactory("Graph", dropFact);

  levelGreedyAndNoPreserve.Request("Aggregates", aggFact.get());
  levelGreedyAndNoPreserve.Request("UnAmalgamationInfo", amalgFact.get());

  levelGreedyAndNoPreserve.Request(*aggFact);
  aggFact->Build(levelGreedyAndNoPreserve);
  aggregates = levelGreedyAndNoPreserve.Get<RCP<Aggregates>>("Aggregates", aggFact.get());
  aggFact->SetParameter("aggregation: preserve Dirichlet points", Teuchos::ParameterEntry(false));

  aggregates->ComputeNodesInAggregate(aggPtr, aggNodes, unaggregated);
  unaggregated_h = Kokkos::create_mirror_view(unaggregated);
  Kokkos::deep_copy(unaggregated_h, unaggregated);

  // GH: loop over the unaggregated list and add to the counter each time node 2 appears
  int unaggregated_count_node2 = 0;
  for (LO i = 0; i < (LO)unaggregated_h.extent(0); i++)
    if (unaggregated_h(i) == (LO)2)
      unaggregated_count_node2++;

  TEST_EQUALITY(unaggregated_count_node2, 1);  // check that the node with the Dof flagged as dirichlet appears in the list once

  // GH: print this to out since the test can be sensitive to ordering issues
  out << "unaggregated_h = [";
  for (LO i = 0; i < (LO)unaggregated_h.extent(0) - 1; i++) {
    out << unaggregated_h(i) << ", ";
  }
  out << unaggregated_h(unaggregated_h.extent(0) - 1) << "]\n";

  // Repeat with greedy Dirichlet and preserve Dirichlet points
  Level levelGreedyAndPreserve;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(levelGreedyAndPreserve);
  levelGreedyAndPreserve.Set("A", A);

  amalgFact = rcp(new AmalgamationFactory());
  dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetParameter("aggregation: greedy Dirichlet", Teuchos::ParameterEntry(true));
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  // Setup aggregation factory (use default factory for graph)
  aggFact = rcp(new UncoupledAggregationFactory());
  aggFact->SetFactory("Graph", dropFact);
  aggFact->SetParameter("aggregation: preserve Dirichlet points", Teuchos::ParameterEntry(true));

  levelGreedyAndPreserve.Request("Aggregates", aggFact.get());
  levelGreedyAndPreserve.Request("UnAmalgamationInfo", amalgFact.get());

  levelGreedyAndPreserve.Request(*aggFact);
  aggFact->Build(levelGreedyAndPreserve);
  aggregates = levelGreedyAndPreserve.Get<RCP<Aggregates>>("Aggregates", aggFact.get());

  aggregates->ComputeNodesInAggregate(aggPtr, aggNodes, unaggregated);
  unaggregated_h = Kokkos::create_mirror_view(unaggregated);
  Kokkos::deep_copy(unaggregated_h, unaggregated);
  for (LO i = 0; i < (LO)unaggregated_h.extent(0); i++) {
    TEST_EQUALITY(unaggregated_h(i) == 2, false);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(FilteredA, RootStencil, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  // Make a "hot dog" matrix
  Teuchos::ParameterList matrixParams;
  matrixParams.set("matrixType", "Laplace2D");
  double factor = 10;
  matrixParams.set("stretchx", 1.0);
  matrixParams.set("stretchy", factor);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildMatrix(matrixParams, TestHelpers::Parameters::getLib());

  RCP<AmalgamationInfo> amalgInfo;
  Level level;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
  level.Set("A", A);

  RCP<CoalesceDropFactory> dropFact;
  RCP<AmalgamationFactory> amalgFact;
  RCP<UncoupledAggregationFactory> aggFact;
  double dropTol = 0.025;

  amalgFact = rcp(new AmalgamationFactory());
  dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(dropTol));
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  // Setup aggregation factory (use default factory for graph)
  aggFact = rcp(new UncoupledAggregationFactory());
  aggFact->SetFactory("Graph", dropFact);

  RCP<FilteredAFactory> filterFact = rcp(new FilteredAFactory);
  filterFact->SetFactory("UnAmalgamationInfo", amalgFact);
  filterFact->SetFactory("Aggregates", aggFact);
  filterFact->SetFactory("Graph", dropFact);
  filterFact->SetFactory("Filtering", dropFact);

  Teuchos::ParameterList params;
  params.set("filtered matrix: use lumping", true);
  params.set("filtered matrix: use root stencil", true);
  filterFact->SetParameterList(params);
  level.Request("A", filterFact.get());

  filterFact->Build(level);
  RCP<Matrix> Afiltered = level.Get<RCP<Matrix>>("A", filterFact.get());

  // Now check stuff
  //    print_matrix("A",MueLu::Utilities<SC,LO,GO,NO>::Op2NonConstTpetraCrs(A));
  //    print_matrix("Afiltered",MueLu::Utilities<SC,LO,GO,NO>::Op2NonConstTpetraCrs(Afiltered));

  // We use the full graph for the filtered matrix, so the notional nnz should be the same
  TEST_EQUALITY(A->getLocalNumRows() == Afiltered->getLocalNumRows(), true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(FilteredA, SpreadLumpingRootStencil, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  out << "running spread lumping" << std::endl;

  // Make a "hot dog" matrix
  Teuchos::ParameterList matrixParams;
  matrixParams.set("matrixType", "Laplace2D");
  double factor = 10;
  matrixParams.set("stretchx", 1.0);
  matrixParams.set("stretchy", factor);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildMatrix(matrixParams, TestHelpers::Parameters::getLib());

  RCP<AmalgamationInfo> amalgInfo;
  Level level;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
  level.Set("A", A);

  RCP<CoalesceDropFactory> dropFact;
  RCP<AmalgamationFactory> amalgFact;
  RCP<UncoupledAggregationFactory> aggFact;
  double dropTol = 0.025;

  amalgFact = rcp(new AmalgamationFactory());
  dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(dropTol));
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  // Setup aggregation factory (use default factory for graph)
  aggFact = rcp(new UncoupledAggregationFactory());
  aggFact->SetFactory("Graph", dropFact);

  RCP<FilteredAFactory> filterFact = rcp(new FilteredAFactory);
  filterFact->SetFactory("UnAmalgamationInfo", amalgFact);
  filterFact->SetFactory("Aggregates", aggFact);
  filterFact->SetFactory("Graph", dropFact);
  filterFact->SetFactory("Filtering", dropFact);

  Teuchos::ParameterList params;
  params.set("filtered matrix: use lumping", true);
  params.set("filtered matrix: use root stencil", true);
  params.set("filtered matrix: use spread lumping", true);
  params.set("filtered matrix: spread lumping diag dom growth factor", 1.1);
  params.set("filtered matrix: spread lumping diag dom cap", 2.0);
  filterFact->SetParameterList(params);
  level.Request("A", filterFact.get());

  filterFact->Build(level);
  RCP<Matrix> Afiltered = level.Get<RCP<Matrix>>("A", filterFact.get());

  // Now check stuff
  //    print_matrix("A",MueLu::Utilities<SC,LO,GO,NO>::Op2NonConstTpetraCrs(A));
  //    print_matrix("Afiltered",MueLu::Utilities<SC,LO,GO,NO>::Op2NonConstTpetraCrs(Afiltered));

  // We use the full graph for the filtered matrix, so the notional nnz should be the same
  TEST_EQUALITY(A->getLocalNumRows() == Afiltered->getLocalNumRows(), true);
}  // SpreadLumping

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(FilteredA, SpreadLumpingReuseGraph, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  out << "running spread lumping" << std::endl;

  // Make a "hot dog" matrix
  Teuchos::ParameterList matrixParams;
  matrixParams.set("matrixType", "Laplace2D");
  double factor = 10;
  matrixParams.set("stretchx", 1.0);
  matrixParams.set("stretchy", factor);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildMatrix(matrixParams, TestHelpers::Parameters::getLib());

  RCP<AmalgamationInfo> amalgInfo;
  Level level;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
  level.Set("A", A);

  RCP<CoalesceDropFactory> dropFact;
  RCP<AmalgamationFactory> amalgFact;
  RCP<UncoupledAggregationFactory> aggFact;
  double dropTol = 0.025;

  amalgFact = rcp(new AmalgamationFactory());
  dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(dropTol));
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  // Setup aggregation factory (use default factory for graph)
  aggFact = rcp(new UncoupledAggregationFactory());
  aggFact->SetFactory("Graph", dropFact);

  RCP<FilteredAFactory> filterFact = rcp(new FilteredAFactory);
  filterFact->SetFactory("UnAmalgamationInfo", amalgFact);
  filterFact->SetFactory("Aggregates", aggFact);
  filterFact->SetFactory("Graph", dropFact);
  filterFact->SetFactory("Filtering", dropFact);

  Teuchos::ParameterList params;
  params.set("filtered matrix: use lumping", true);
  params.set("filtered matrix: use spread lumping", false);
  params.set("filtered matrix: reuse graph", true);
  params.set("filtered matrix: spread lumping diag dom growth factor", 1.1);
  params.set("filtered matrix: spread lumping diag dom cap", 2.0);
  filterFact->SetParameterList(params);
  level.Request("A", filterFact.get());

  filterFact->Build(level);
  RCP<Matrix> Afiltered = level.Get<RCP<Matrix>>("A", filterFact.get());

  // We use the full graph for the filtered matrix, so the notional nnz should be the same
  TEST_EQUALITY(A->getLocalNumRows() == Afiltered->getLocalNumRows(), true);
}  // SpreadLumpingReuseGraph

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(FilteredA, SpreadLumpingNoStencilRootNoReuseGraph, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  out << "running spread lumping" << std::endl;

  // Make a "hot dog" matrix
  Teuchos::ParameterList matrixParams;
  matrixParams.set("matrixType", "Laplace2D");
  double factor = 10;
  matrixParams.set("stretchx", 1.0);
  matrixParams.set("stretchy", factor);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildMatrix(matrixParams, TestHelpers::Parameters::getLib());

  RCP<AmalgamationInfo> amalgInfo;
  Level level;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
  level.Set("A", A);

  RCP<CoalesceDropFactory> dropFact;
  RCP<AmalgamationFactory> amalgFact;
  RCP<UncoupledAggregationFactory> aggFact;
  double dropTol = 0.025;

  amalgFact = rcp(new AmalgamationFactory());
  dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(dropTol));
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  // Setup aggregation factory (use default factory for graph)
  aggFact = rcp(new UncoupledAggregationFactory());
  aggFact->SetFactory("Graph", dropFact);

  RCP<FilteredAFactory> filterFact = rcp(new FilteredAFactory);
  filterFact->SetFactory("UnAmalgamationInfo", amalgFact);
  filterFact->SetFactory("Aggregates", aggFact);
  filterFact->SetFactory("Graph", dropFact);
  filterFact->SetFactory("Filtering", dropFact);

  Teuchos::ParameterList params;
  params.set("filtered matrix: use lumping", true);
  params.set("filtered matrix: use spread lumping", false);
  params.set("filtered matrix: reuse graph", false);
  params.set("filtered matrix: spread lumping diag dom growth factor", 1.1);
  params.set("filtered matrix: spread lumping diag dom cap", 2.0);
  filterFact->SetParameterList(params);
  level.Request("A", filterFact.get());

  filterFact->Build(level);
  RCP<Matrix> Afiltered = level.Get<RCP<Matrix>>("A", filterFact.get());

  // We use the full graph for the filtered matrix, so the notional nnz should be the same
  TEST_EQUALITY(A->getLocalNumRows() == Afiltered->getLocalNumRows(), true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(FilteredA, SpreadNoLumpingNoStencilRootNoReuseGraph, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  out << "running spread lumping" << std::endl;

  // Make a "hot dog" matrix
  Teuchos::ParameterList matrixParams;
  matrixParams.set("matrixType", "Laplace2D");
  double factor = 10;
  matrixParams.set("stretchx", 1.0);
  matrixParams.set("stretchy", factor);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildMatrix(matrixParams, TestHelpers::Parameters::getLib());

  RCP<AmalgamationInfo> amalgInfo;
  Level level;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
  level.Set("A", A);

  RCP<CoalesceDropFactory> dropFact;
  RCP<AmalgamationFactory> amalgFact;
  RCP<UncoupledAggregationFactory> aggFact;
  double dropTol = 0.025;

  amalgFact = rcp(new AmalgamationFactory());
  dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(dropTol));
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  // Setup aggregation factory (use default factory for graph)
  aggFact = rcp(new UncoupledAggregationFactory());
  aggFact->SetFactory("Graph", dropFact);

  RCP<FilteredAFactory> filterFact = rcp(new FilteredAFactory);
  filterFact->SetFactory("UnAmalgamationInfo", amalgFact);
  filterFact->SetFactory("Aggregates", aggFact);
  filterFact->SetFactory("Graph", dropFact);
  filterFact->SetFactory("Filtering", dropFact);

  Teuchos::ParameterList params;
  params.set("filtered matrix: use lumping", false);
  params.set("filtered matrix: use spread lumping", false);
  params.set("filtered matrix: reuse graph", false);
  params.set("filtered matrix: spread lumping diag dom growth factor", 1.1);
  params.set("filtered matrix: spread lumping diag dom cap", 2.0);
  filterFact->SetParameterList(params);
  level.Request("A", filterFact.get());

  filterFact->Build(level);
  RCP<Matrix> Afiltered = level.Get<RCP<Matrix>>("A", filterFact.get());

  // We use the full graph for the filtered matrix, so the notional nnz should be the same
  TEST_EQUALITY(A->getLocalNumRows() == Afiltered->getLocalNumRows(), true);
}  // SpreadNoLumpingNoStencilRootNoReuseGraph

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates, AllowDroppingToCreateAdditionalDirichletRows, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  out << "Test option that allows dropping during aggregation to create new Dirichlet rows" << std::endl;

  typedef typename Teuchos::ScalarTraits<Scalar> TST;

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(15);
  A->resumeFill();

  // Create one row on every processor with small off-diagonal entries that will be dropped with
  // an appropriately chosen threshold.  Avoid domain boundaries.
  LocalOrdinal localRowToModify = 1;
  Teuchos::ArrayView<const LocalOrdinal> indices;
  Teuchos::ArrayView<const Scalar> values;
  A->getLocalRowView(localRowToModify, indices, values);
  Array<Scalar> newvalues(values.size(), TST::zero());
  for (int j = 0; j < indices.size(); j++) {
    if (indices[j] == localRowToModify)
      newvalues[j] = values[j];  // keep diagonal unmodified
    else
      newvalues[j] = -TST::eps();
  }
  A->replaceLocalValues(localRowToModify, indices, newvalues);
  A->fillComplete();

  // Dropping connections will not lead to the creation of new Dirichlet rows.
  RCP<AmalgamationInfo> amalgInfo;
  Level level;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
  level.Set("A", A);

  RCP<CoalesceDropFactory> dropFact = rcp(new CoalesceDropFactory());
  dropFact->SetParameter("aggregation: dropping may create Dirichlet", Teuchos::ParameterEntry(false));
  dropFact->SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(Teuchos::as<double>(-100 * TST::eps())));
  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
  level.Request("Graph", dropFact.get());

  // Setup aggregation factory (use default factory for graph)
  RCP<UncoupledAggregationFactory> aggFact = rcp(new UncoupledAggregationFactory());
  aggFact->SetFactory("Graph", dropFact);

  level.Request("Aggregates", aggFact.get());
  level.Request("UnAmalgamationInfo", amalgFact.get());

  level.Request(*aggFact);
  aggFact->Build(level);
  RCP<Aggregates> aggregates                = level.Get<RCP<Aggregates>>("Aggregates", aggFact.get());
  RCP<GraphBase> graph                      = level.Get<RCP<GraphBase>>("Graph", dropFact.get());
  ArrayRCP<const bool> dirichletBoundaryMap = graph->GetBoundaryNodeMap();
  int numDirichletRows                      = 0;
  LO numRows                                = graph->GetNodeNumVertices();
  for (LO i = 0; i < numRows; i++)
    if (dirichletBoundaryMap[i] == true)
      numDirichletRows++;
  TEST_EQUALITY(numDirichletRows, 0);

  Array<LO> aggPtr;
  Array<LO> aggNodes;
  Array<LO> unaggregated;

  // Repeat with "aggregation: dropping may create Dirichlet" = TRUE, i.e.,
  // dropping connections may lead to the creation of new Dirichlet rows.
  // The second row should be flagged as Dirichlet because all off-diagonal entries are dropped.
  amalgFact = rcp(new AmalgamationFactory());
  dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetParameter("aggregation: dropping may create Dirichlet", Teuchos::ParameterEntry(true));
  dropFact->SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(Teuchos::as<double>(-100 * TST::eps())));
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
  level.Request("Graph", dropFact.get());

  // Setup aggregation factory (use default factory for graph)
  aggFact = rcp(new UncoupledAggregationFactory());
  aggFact->SetFactory("Graph", dropFact);

  level.Request("Aggregates", aggFact.get());
  level.Request("UnAmalgamationInfo", amalgFact.get());

  level.Request(*aggFact);
  aggFact->Build(level);
  aggregates           = level.Get<RCP<Aggregates>>("Aggregates", aggFact.get());
  graph                = level.Get<RCP<GraphBase>>("Graph", dropFact.get());
  dirichletBoundaryMap = graph->GetBoundaryNodeMap();
  numDirichletRows     = 0;
  for (LO i = 0; i < numRows; i++)
    if (dirichletBoundaryMap[i] == true)
      numDirichletRows++;
  TEST_EQUALITY(numDirichletRows, 1);

}  // AllowDroppingToCreateAdditionalDirichletRows

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node)                                                                          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates, JustAggregation, Scalar, LO, GO, Node)                              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates, GetNumAggregates, Scalar, LO, GO, Node)                             \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates, JustUncoupledAggregation, Scalar, LO, GO, Node)                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates, JustUncoupledAggregationFactory, Scalar, LO, GO, Node)              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates, GetNumUncoupledAggregates, Scalar, LO, GO, Node)                    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates, UncoupledPhase1, Scalar, LO, GO, Node)                              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates, UncoupledPhase2, Scalar, LO, GO, Node)                              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates, UncoupledPhase3, Scalar, LO, GO, Node)                              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates, UncoupledInterface, Scalar, LO, GO, Node)                           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates, JustStructuredAggregationGlobal, Scalar, LO, GO, Node)              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates, JustStructuredAggregationLocal, Scalar, LO, GO, Node)               \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates, GreedyDirichlet, Scalar, LO, GO, Node)                              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates, AllowDroppingToCreateAdditionalDirichletRows, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(FilteredA, RootStencil, Scalar, LO, GO, Node)                                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(FilteredA, SpreadLumpingRootStencil, Scalar, LO, GO, Node)                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(FilteredA, SpreadLumpingReuseGraph, Scalar, LO, GO, Node)                       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(FilteredA, SpreadLumpingNoStencilRootNoReuseGraph, Scalar, LO, GO, Node)        \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(FilteredA, SpreadNoLumpingNoStencilRootNoReuseGraph, Scalar, LO, GO, Node)

//  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates,HybridAggregation,Scalar,LO,GO,Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
