// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_FE_GRAPH_TEST_UTILS_HPP
#define PANZER_FE_GRAPH_TEST_UTILS_HPP

// Utilities for testing the FE graphs of the Tpetra linear object factories: the
// reference graph built by insertion, and an exact comparison of two FECrsGraphs.

#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_ConnManager.hpp"
#include "Tpetra_FECrsGraph.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultMpiComm.hpp"

#include <algorithm>
#include <ostream>
#include <string>
#include <vector>

namespace panzer::fe_graph_test {

using LO      = int;
using GO      = panzer::GlobalOrdinal;
using NT      = panzer::TpetraNodeType;
using MapType = Tpetra::Map<LO, GO, NT>;
using FEGraph = Tpetra::FECrsGraph<LO, GO, NT>;
using Teuchos::RCP;
using Teuchos::rcp;

// The reference: an FECrsGraph built by insertion, exactly as the factories used to
// build their FE graphs.  Each element couples its row DOFs (rowIndexer) with its
// column DOFs (colIndexer) and its associated neighbors' column DOFs.
inline RCP<FEGraph> insertionFEGraph(const panzer::GlobalIndexer& rowIndexer, const panzer::GlobalIndexer& colIndexer,
                              const panzer::GlobalIndexer& elementIndexer,
                              const RCP<const panzer::ConnManager>& connMgr, const RCP<const MapType>& ownedRowMap,
                              const RCP<const MapType>& ownedPlusSharedRowMap,
                              const RCP<const MapType>& ownedPlusSharedDomainMap,
                              const RCP<const MapType>& ownedDomainMap) {
  const bool han = !connMgr.is_null() && connMgr->hasAssociatedNeighbors();
  std::vector<std::string> blockIds;
  elementIndexer.getElementBlockIds(blockIds);

  auto forEachElement = [&](auto&& f) {
    std::vector<GO> rowGIDs, colGIDs, neighborGIDs;
    for (const auto& blockId : blockIds)
      for (const LO element : elementIndexer.getElementBlock(blockId)) {
        rowIndexer.getElementGIDs(element, rowGIDs);
        colIndexer.getElementGIDs(element, colGIDs);
        if (han)
          for (const LO neighbor : connMgr->getAssociatedNeighbors(element)) {
            colIndexer.getElementGIDs(neighbor, neighborGIDs);
            colGIDs.insert(colGIDs.end(), neighborGIDs.begin(), neighborGIDs.end());
          }
        f(rowGIDs, colGIDs);
      }
  };

  std::vector<size_t> numEntries(ownedPlusSharedRowMap->getLocalNumElements(), 0);
  forEachElement([&](const std::vector<GO>& rows, const std::vector<GO>& cols) {
    for (const GO row : rows) numEntries[ownedPlusSharedRowMap->getLocalElement(row)] += cols.size();
  });
  size_t maxEntries = 1;
  for (const size_t n : numEntries) maxEntries = std::max(maxEntries, n);

  RCP<FEGraph> graph = rcp(new FEGraph(ownedRowMap, ownedPlusSharedRowMap, maxEntries, ownedPlusSharedDomainMap,
                                       Teuchos::null, ownedDomainMap));
  auto params = Teuchos::parameterList();
  params->set("Check Col GIDs In At Least One Owned Row", false);
  graph->setParameterList(params);
  graph->beginAssembly();
  forEachElement([&](const std::vector<GO>& rows, const std::vector<GO>& cols) {
    for (const GO row : rows) graph->insertGlobalIndices(row, Teuchos::ArrayView<const GO>(cols));
  });
  graph->endAssembly();
  return graph;
}

// Whether two maps have the same global indices in the same order on this process.
inline bool sameLocalGIDs(const MapType& a, const MapType& b) {
  if (a.getLocalNumElements() != b.getLocalNumElements()) return false;
  auto aGIDs = a.getMyGlobalIndices();
  auto bGIDs = b.getMyGlobalIndices();
  for (size_t i = 0; i < aGIDs.extent(0); ++i)
    if (aGIDs(i) != bGIDs(i)) return false;
  return true;
}

// Whether two fill-complete graphs are identical on this process.
inline bool identicalLocalGraphs(const Tpetra::CrsGraph<LO, GO, NT>& ref, const Tpetra::CrsGraph<LO, GO, NT>& test,
                          const std::string& label, std::ostream& out) {
  auto report = [&](const std::string& what) {
    out << label << ": " << what << std::endl;
    return false;
  };
  if (!ref.isFillComplete() || !test.isFillComplete()) return report("not fill complete");
  if (!sameLocalGIDs(*ref.getRowMap(), *test.getRowMap())) return report("row maps differ");
  if (!sameLocalGIDs(*ref.getColMap(), *test.getColMap())) return report("column maps differ");
  if (!sameLocalGIDs(*ref.getDomainMap(), *test.getDomainMap())) return report("domain maps differ");
  if (!sameLocalGIDs(*ref.getRangeMap(), *test.getRangeMap())) return report("range maps differ");
  auto refGraph  = ref.getLocalGraphHost();
  auto testGraph = test.getLocalGraphHost();
  if (refGraph.row_map.extent(0) != testGraph.row_map.extent(0)) return report("number of rows differs");
  for (size_t i = 0; i < refGraph.row_map.extent(0); ++i)
    if (refGraph.row_map(i) != testGraph.row_map(i)) return report("row offsets differ");
  const size_t nnz = refGraph.row_map.extent(0) ? refGraph.row_map(refGraph.row_map.extent(0) - 1) : 0;
  for (size_t k = 0; k < nnz; ++k)
    if (refGraph.entries(k) != testGraph.entries(k)) return report("column indices differ");
  return true;
}

// Whether two fully-assembled FECrsGraphs have identical owned and owned+shared graphs
// on every process.  Collective.
inline bool identicalFEGraphs(FEGraph& ref, FEGraph& test, const Teuchos::Comm<int>& comm, std::ostream& out) {
  bool ok = identicalLocalGraphs(ref, test, "owned graph", out);
  ref.switchActiveCrsGraph();
  test.switchActiveCrsGraph();
  ok = identicalLocalGraphs(ref, test, "owned+shared graph", out) && ok;
  ref.switchActiveCrsGraph();
  test.switchActiveCrsGraph();
  int lclOK = ok ? 1 : 0, gblOK = 0;
  Teuchos::reduceAll<int, int>(comm, Teuchos::REDUCE_MIN, lclOK, Teuchos::outArg(gblOK));
  return gblOK == 1;
}

} // namespace panzer::fe_graph_test

#endif
