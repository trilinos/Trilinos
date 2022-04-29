// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#include "Panzer_OrientationsInterface.hpp"

#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_ConnManager.hpp"
#include "PanzerDiscFE_config.hpp"
#include "Panzer_NodalFieldPattern.hpp"
#include "Panzer_LocalPartitioningUtilities.hpp"
#include "Panzer_NodeType.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_MultiVector.hpp"

namespace panzer {

namespace
{

void
buildIntrepidOrientation(const Teuchos::RCP<const Teuchos::Comm<int>> & comm,
                         panzer::ConnManager & conn,
                         std::vector<Intrepid2::Orientation> & orientations)
{

  using MVector = Tpetra::MultiVector<panzer::GlobalOrdinal, panzer::LocalOrdinal, panzer::GlobalOrdinal, panzer::TpetraNodeType>;
  using Map = Tpetra::Map<panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType>;
  using Importer = Tpetra::Import<panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType>;
  using NodeView = Kokkos::View<panzer::GlobalOrdinal*, Kokkos::DefaultHostExecutionSpace>;

  // First we need to build the indexing scheme
  PHX::View<panzer::GlobalOrdinal*> owned_cells, ghost_cells, virtual_cells;
  fillLocalCellIDs(comm, conn, owned_cells, ghost_cells, virtual_cells);

  // Build a map and importer for syncing the nodal connectivity
  auto owned_cell_map = Teuchos::rcp(new Map(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),owned_cells,0,comm));
  auto ghost_cell_map = Teuchos::rcp(new Map(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),ghost_cells,0,comm));

  // Build importer: imports from owned cells map to ghost cell map
  auto importer = Teuchos::rcp(new Importer(owned_cell_map,ghost_cell_map));

  // Grab the cell topology from the conn manager
  shards::CellTopology topology;
  {
    // Retrive element blocks and its meta data
    const int numElementBlocks = conn.numElementBlocks();

    std::vector<std::string> elementBlockIds;
    std::vector<shards::CellTopology> elementBlockTopologies;

    conn.getElementBlockIds(elementBlockIds);
    conn.getElementBlockTopologies(elementBlockTopologies);

    TEUCHOS_TEST_FOR_EXCEPTION(numElementBlocks <= 0 &&
                               numElementBlocks != static_cast<int>(elementBlockIds.size()) &&
                               numElementBlocks != static_cast<int>(elementBlockTopologies.size()),
                               std::logic_error,
                               "panzer::buildIntrepidOrientation: Number of element blocks does not match to element block meta data");

    topology = elementBlockTopologies.at(0);
  }
  const int num_nodes_per_cell = topology.getVertexCount();

  // Create Tpetra multivectors for storing global node ids
  auto owned_nodes_vector = Teuchos::rcp(new MVector(owned_cell_map,num_nodes_per_cell));
  auto ghost_nodes_vector = Teuchos::rcp(new MVector(ghost_cell_map,num_nodes_per_cell));

  // Make sure the conn is setup for a nodal connectivity
  panzer::NodalFieldPattern pattern(topology);
  conn.buildConnectivity(pattern);

  const int num_owned_cells = owned_cells.extent(0);
  const int num_ghost_cells = ghost_cells.extent(0);

  // Initialize the orientations vector
  orientations.clear();
  orientations.resize(num_owned_cells+num_ghost_cells);

  // Fill the owned vector with the nodal connectivity of the cells on this processor
  {
    auto vector_view = owned_nodes_vector->getLocalViewHost(Tpetra::Access::OverwriteAll);
    for(int cell=0; cell<owned_cells.extent_int(0); ++cell){
      const GlobalOrdinal * nodes = conn.getConnectivity(cell);
      for(int node=0; node<num_nodes_per_cell; ++node)
        vector_view(cell,node) = nodes[node];
    }
  }

  // Import into the ghost vector
  ghost_nodes_vector->doImport(*owned_nodes_vector,*importer,Tpetra::CombineMode::REPLACE);

  // Add owned orientations
  {
    auto vector_view = owned_nodes_vector->getLocalViewHost(Tpetra::Access::ReadOnly);
    for(int cell=0; cell<num_owned_cells; ++cell){
      NodeView nodes("nodes",num_nodes_per_cell);
      for(int node=0; node<num_nodes_per_cell; ++node)
        nodes(node) = vector_view(cell,node);
      orientations[cell] = Intrepid2::Orientation::getOrientation(topology, nodes);
    }
  }

  // Add ghost orientations
  {
    auto vector_view = ghost_nodes_vector->getLocalViewHost(Tpetra::Access::ReadOnly);
    for(int ghost_cell=0; ghost_cell<num_ghost_cells; ++ghost_cell){
      const int cell = num_owned_cells + ghost_cell;
      NodeView nodes("nodes",num_nodes_per_cell);
      for(int node=0; node<num_nodes_per_cell; ++node)
        nodes(node) = vector_view(ghost_cell,node);
      orientations[cell] = Intrepid2::Orientation::getOrientation(topology, nodes);
    }
  }

}

Teuchos::RCP<std::vector<Intrepid2::Orientation> >
buildIntrepidOrientation(const Teuchos::RCP<const GlobalIndexer> globalIndexer)
{
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::RCP;
  using Teuchos::rcp;

  auto orientation = rcp(new std::vector<Intrepid2::Orientation>);

  auto comm = globalIndexer->getComm();
  auto conn = globalIndexer->getConnManager()->noConnectivityClone();

  TEUCHOS_TEST_FOR_EXCEPTION(conn == Teuchos::null,std::logic_error,
                             "panzer::buildIntrepidOrientation: Could not cast ConnManagerBase");

  buildIntrepidOrientation(comm, *conn, *orientation);
  return orientation;

}

}

OrientationsInterface::
OrientationsInterface(const Teuchos::RCP<const panzer::GlobalIndexer> & indexer)
{
  orientations_ = buildIntrepidOrientation(indexer);
}


Teuchos::RCP<const std::vector<Intrepid2::Orientation> >
OrientationsInterface::
getOrientations() const
{
  TEUCHOS_ASSERT(not orientations_.is_null());
  return orientations_;
}

}
