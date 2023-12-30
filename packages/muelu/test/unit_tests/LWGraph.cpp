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
#include <Teuchos_DefaultComm.hpp>

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
Teuchos::RCP<MueLu::GraphBase<LocalOrdinal, GlobalOrdinal, Node> >
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

  auto graph = level.Get<RCP<GraphBase> >("Graph", dropFact.get());
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

  RCP<GraphBase> graph = gimmeLWGraph(A);

  auto comm          = graph->GetComm();
  const int numRanks = comm->getSize();

  graph->print(out, MueLu::MsgType::Extreme);
  //    TEST_EQUALITY( graph->description() == "LWGraph (graph of A)" ,  true);
  TEST_EQUALITY(graph->description() == "MueLu.description()", true);  // To fix when LWGRaph description method will be fixed
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
    TEST_EQUALITY(graph->getNeighborVertices(i).size() == columns.view(rows[i], rows[i + 1] - rows[i]).size(), true);
  }

  auto crsGraph = graphLWK->GetCrsGraph();
  TEST_EQUALITY(crsGraph->getLocalNumEntries() == graph->GetNodeNumEdges(), true);
  TEST_EQUALITY(crsGraph->getLocalNumRows() == graph->GetNodeNumVertices(), true);
  TEST_EQUALITY(crsGraph->getLocalMaxNumRowEntries() == graph->getLocalMaxNumRowEntries(), true);

  TEST_THROW(graphLWK->SetBoundaryNodeMap(graphLWK->GetDomainMap()), MueLu::Exceptions::NotImplemented);

}  // CreateLWGraph

#define MUELU_ETI_GROUP(SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(LWGraph, CreateLWGraph, SC, LO, GO, NO)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
