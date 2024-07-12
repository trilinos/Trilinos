// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <Kokkos_Core.hpp>

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_LWGraph.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_PreDropFunctionConstVal.hpp"

namespace MueLuTests {

// Little utility to generate a LWGraph.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<MueLu::LWGraph<LocalOrdinal, GlobalOrdinal, Node> >
gimmeLWGraph(const Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& A) {
#include "MueLu_UseShortNames.hpp"

  Level level;

  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
  level.Set("A", A);

  RCP<CoalesceDropFactory> dropFact = rcp(new CoalesceDropFactory());
  dropFact->SetVerbLevel(MueLu::Extreme);
  dropFact->SetPreDropFunction(rcp(new PreDropFunctionConstVal(0.00001)));

  level.Request("Graph", dropFact.get());
  dropFact->Build(level);

  auto graph = level.Get<RCP<LWGraph> >("Graph", dropFact.get());
  level.Release("Graph", dropFact.get());
  return graph;
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(LWGraph, CreateLWGraph, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  if (TestHelpers::Parameters::getLib() == Xpetra::UseEpetra) {
    out << "skipping test for linAlgebra==UseEpetra" << std::endl;
    return;
  }

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(16);

  RCP<LWGraph> graph = gimmeLWGraph(A);

  auto comm          = graph->GetComm();
  const int numRanks = comm->getSize();

  graph->print(out, MueLu::MsgType::Extreme);
  //    TEST_EQUALITY( graph->description() == "LWGraph (graph of A)" ,  true);
  // TEST_EQUALITY(graph->description() == "MueLu.description()", true);  // To fix when LWGRaph description method will be fixed
  TEST_EQUALITY(graph != Teuchos::null, true);
  auto graphLWK = dynamic_cast<LWGraph*>(graph.get());
  auto rows     = graphLWK->getRowPtrs();

  if (numRanks == 1) {
    TEST_EQUALITY(graph->getLocalMaxNumRowEntries() == 3, true);
    TEST_EQUALITY(graph->GetNodeNumVertices() == 16, true);
    TEST_EQUALITY(graph->GetNodeNumEdges() == 46, true);
    TEST_EQUALITY(rows.size() == 17, true);

  } else if (numRanks == 4) {
    TEST_EQUALITY(graph->GetNodeNumVertices() == 4, true);
  }

  auto columns = graphLWK->getEntries();
  for (size_t i = 0; i < graph->GetNodeNumVertices(); ++i) {
    TEST_EQUALITY(graph->getNeighborVertices(i).length == Teuchos::as<LocalOrdinal>(Kokkos::subview(columns,
                                                                                                    Kokkos::make_pair(rows[i], rows[i + 1]))
                                                                                        .extent(0)),
                  true);
  }

  auto crsGraph = graphLWK->GetCrsGraph();
  TEST_EQUALITY(crsGraph->getLocalNumEntries() == graph->GetNodeNumEdges(), true);
  TEST_EQUALITY(crsGraph->getLocalNumRows() == graph->GetNodeNumVertices(), true);
  TEST_EQUALITY(crsGraph->getLocalMaxNumRowEntries() == graph->getLocalMaxNumRowEntries(), true);

}  // CreateLWGraph

#define MUELU_ETI_GROUP(SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(LWGraph, CreateLWGraph, SC, LO, GO, NO)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
