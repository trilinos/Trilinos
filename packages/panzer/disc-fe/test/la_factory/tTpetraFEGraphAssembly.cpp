// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// The FE graphs of TpetraLinearObjFactory and BlockedTpetraLinearObjFactory are
// assembled on device (panzer::buildFEGraphGeometry + Tpetra::assembleFECrsGraphWithMaps).
// These tests check that they are identical to the graphs built by inserting each
// element's couplings into a Tpetra::FECrsGraph between beginAssembly() and
// endAssembly() (which is how the factories built them before): same row, column,
// domain and range maps (the column map compared index by index), and the same local
// row offsets and column indices, in the same order, for both the owned and the
// owned+shared graphs.
//
// They cover both ways the element connectivity is gathered: from the indexers'
// device local ID tables (DOFManager), and from global IDs on host (an indexer
// without a local ID table, or associated neighbors).

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <map>
#include <string>
#include <vector>

#include "Panzer_TpetraLinearObjFactory.hpp"
#include "Panzer_BlockedTpetraLinearObjFactory.hpp"
#include "Panzer_FEGraphGeometry.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"

#include "UnitTest_ConnManager.hpp"
#include "UnitTest_GlobalIndexer.hpp"

#include "FEGraphTestUtils.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

namespace {

using namespace panzer::fe_graph_test;
using LOFType   = panzer::TpetraLinearObjFactory<panzer::Traits, double, LO, GO>;
using BLOFType  = panzer::BlockedTpetraLinearObjFactory<panzer::Traits, double, LO, GO>;

RCP<const Teuchos::Comm<int>> getComm() {
  return rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));
}

template <typename Intrepid2Type>
RCP<const panzer::FieldPattern> buildFieldPattern() {
  RCP<Intrepid2::Basis<PHX::Device, double, double>> basis = rcp(new Intrepid2Type);
  return rcp(new panzer::Intrepid2FieldPattern(basis));
}

// A ConnManager whose elements each have an associated neighbor: the next element of
// the same element block (if any).
class NeighborConnManager : public panzer::unit_test::ConnManager {
 public:
  NeighborConnManager(int rank, int numProcs) : panzer::unit_test::ConnManager(rank, numProcs) {
    std::vector<std::string> blockIds;
    getElementBlockIds(blockIds);
    for (const auto& blockId : blockIds) {
      const auto& elements = getElementBlock(blockId);
      for (std::size_t k = 0; k + 1 < elements.size(); ++k) neighbors_[elements[k]].push_back(elements[k + 1]);
    }
  }
  const std::vector<int>& getAssociatedNeighbors(const panzer::LocalOrdinal& el) const override {
    auto it = neighbors_.find(el);
    return it == neighbors_.end() ? empty_ : it->second;
  }
  bool hasAssociatedNeighbors() const override { return true; }

 private:
  std::map<int, std::vector<int>> neighbors_;
  std::vector<int> empty_;
};

RCP<panzer::DOFManager> buildDOFManager(const RCP<panzer::ConnManager>& connManager,
                                        const std::vector<std::string>& fields) {
  RCP<const panzer::FieldPattern> patternC1 =
      buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::Device, double, double>>();
  RCP<panzer::DOFManager> dofManager = rcp(new panzer::DOFManager());
  dofManager->setConnManager(connManager, MPI_COMM_WORLD);
  for (const auto& field : fields) dofManager->addField(field, patternC1);
  dofManager->buildGlobalUnknowns();
  return dofManager;
}

// Check a TpetraLinearObjFactory's FE graph against the insertion reference, and that
// the connectivity was gathered with local indices iff expected.
void checkTpetraLOF(const RCP<const panzer::GlobalIndexer>& rowIndexer,
                    const RCP<const panzer::GlobalIndexer>& colIndexer, const bool expectLocalIndices,
                    Teuchos::FancyOStream& out, bool& success) {
  RCP<LOFType> factory = colIndexer.is_null() ? rcp(new LOFType(getComm(), rowIndexer, true))
                                              : rcp(new LOFType(getComm(), rowIndexer, colIndexer, true));
  const auto& cols   = colIndexer.is_null() ? *rowIndexer : *colIndexer;
  auto connMgr       = cols.getConnManager();

  auto geometry = panzer::buildFEGraphGeometry<LO, GO, NT>(*rowIndexer, cols, *rowIndexer, connMgr);
  TEST_EQUALITY(geometry.hasLocalIndices(), expectLocalIndices);

  RCP<FEGraph> ref = insertionFEGraph(*rowIndexer, cols, *rowIndexer, connMgr, factory->getMap(),
                                      factory->getGhostedMap(), factory->getGhostedColMap(), factory->getColMap());
  RCP<FEGraph> feGraph = factory->getFEGraph();
  TEST_ASSERT(identicalFEGraphs(*ref, *feGraph, *getComm(), out));
}

} // namespace

// panzer::unit_test::GlobalIndexer has no local ID table: global indices gathered on host.
TEUCHOS_UNIT_TEST(tTpetraFEGraphAssembly, test_indexer)
{
  auto comm = getComm();
  RCP<panzer::GlobalIndexer> indexer = rcp(new panzer::unit_test::GlobalIndexer(comm->getRank(), comm->getSize()));
  checkTpetraLOF(indexer, Teuchos::null, false, out, success);
}

// A DOFManager with two fields: local indices gathered on device.
TEUCHOS_UNIT_TEST(tTpetraFEGraphAssembly, dof_manager)
{
  auto comm = getComm();
  RCP<panzer::ConnManager> connManager = rcp(new panzer::unit_test::ConnManager(comm->getRank(), comm->getSize()));
  checkTpetraLOF(buildDOFManager(connManager, {"U", "V"}), Teuchos::null, true, out, success);
}

// Separate row and column DOFManagers (rectangular): local indices gathered on device.
TEUCHOS_UNIT_TEST(tTpetraFEGraphAssembly, row_and_column_dof_managers)
{
  auto comm = getComm();
  RCP<panzer::ConnManager> connManager = rcp(new panzer::unit_test::ConnManager(comm->getRank(), comm->getSize()));
  RCP<panzer::DOFManager> rows = buildDOFManager(connManager, {"U"});
  RCP<panzer::DOFManager> cols = buildDOFManager(connManager, {"U", "V"});
  checkTpetraLOF(rows, cols, true, out, success);
}

// Associated neighbors: global indices gathered on host, including the neighbors' columns.
TEUCHOS_UNIT_TEST(tTpetraFEGraphAssembly, associated_neighbors)
{
  auto comm = getComm();
  RCP<panzer::ConnManager> connManager = rcp(new NeighborConnManager(comm->getRank(), comm->getSize()));
  checkTpetraLOF(buildDOFManager(connManager, {"U", "V"}), Teuchos::null, false, out, success);
}

// BlockedTpetraLinearObjFactory: every block (i,j) of a three-field blocked DOFManager.
TEUCHOS_UNIT_TEST(tTpetraFEGraphAssembly, blocked)
{
  auto comm = getComm();
  RCP<const panzer::FieldPattern> patternC1 =
      buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::Device, double, double>>();
  RCP<panzer::ConnManager> connManager = rcp(new panzer::unit_test::ConnManager(comm->getRank(), comm->getSize()));
  RCP<panzer::BlockedDOFManager> blkIndexer = rcp(new panzer::BlockedDOFManager());
  blkIndexer->setConnManager(connManager, MPI_COMM_WORLD);
  // Block 0 has two fields, so its diagonal and off-diagonal blocks have different sizes.
  std::vector<std::vector<std::string>> fieldOrder = {{"U", "V"}, {"W"}, {"X"}};
  for (const auto& block : fieldOrder)
    for (const auto& field : block) blkIndexer->addField(field, patternC1);
  blkIndexer->setFieldOrder(fieldOrder);
  blkIndexer->buildGlobalUnknowns();

  const auto& indexers = blkIndexer->getFieldDOFManagers();
  const int numBlocks  = static_cast<int>(indexers.size());
  TEST_EQUALITY(numBlocks, 3);

  RCP<BLOFType> factory = rcp(new BLOFType(rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD)), blkIndexer, true));
  for (int i = 0; i < numBlocks; ++i)
    for (int j = 0; j < numBlocks; ++j) {
      out << "block (" << i << "," << j << ")" << std::endl;
      auto geometry = panzer::buildFEGraphGeometry<LO, GO, NT>(*indexers[i], *indexers[j], *indexers[0], Teuchos::null);
      TEST_ASSERT(geometry.hasLocalIndices());
      RCP<FEGraph> ref = insertionFEGraph(*indexers[i], *indexers[j], *indexers[0], Teuchos::null, factory->getMap(i),
                                          factory->getGhostedMap(i), factory->getGhostedMap(j), factory->getMap(j));
      RCP<FEGraph> feGraph = factory->getFEGraph(i, j);
      TEST_ASSERT(identicalFEGraphs(*ref, *feGraph, *getComm(), out));
    }
}
