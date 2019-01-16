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

#include "Panzer_LocalMeshInfo.hpp"

// Phalanx includes
#include "Phalanx_KokkosDeviceTypes.hpp"

// Teuchos includes
#include "Teuchos_Assert.hpp"

// Tpetra includes
#include "Tpetra_Import.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"

// Panzer includes
#include "Panzer_NodeType.hpp"
#include "Panzer_HashUtils.hpp"
#include "Panzer_LocalMeshInfo.hpp"
#include "Panzer_LocalPartitioningUtilities.hpp"
#include "Panzer_FaceToElement.hpp"
#include "Panzer_ConnManager.hpp"
#include "Panzer_FieldPattern.hpp"
#include "Panzer_NodalFieldPattern.hpp"
#include "Panzer_EdgeFieldPattern.hpp"
#include "Panzer_FaceFieldPattern.hpp"
#include "Panzer_ElemFieldPattern.hpp"

// STL includes
#include <string>
#include <map>
#include <set>
#include <vector>
#include <unordered_set>

namespace panzer
{

namespace
{

/** Build a Kokkos array of all the global cell IDs from a connection manager.
  * Note that this is mapping between local IDs to global IDs.
  */
PHX::View<const GlobalOrdinal*>
buildOwnedCellsGlobalIDs(const panzer::ConnManager & cell_conn)
{

  // calculate total number of local cells
  std::vector<std::string> element_blocks;
  cell_conn.getElementBlockIds(element_blocks);

  // Note: getElementBlock only returns OWNED cells
  std::size_t num_owned_cells = 0;
  for(const auto & element_block : element_blocks)
    num_owned_cells += cell_conn.getElementBlock(element_block).size();

  Kokkos::View<panzer::GlobalOrdinal*> global_cell_ids("global_cell_ids",num_owned_cells);

  for(const auto & element_block : element_blocks){
    const auto & owned_cell_local_ids = cell_conn.getElementBlock(element_block);

    for(size_t i=0; i<owned_cell_local_ids.size(); ++i){
      const auto local_cell_id = owned_cell_local_ids[i];

      // Sanity check
      TEUCHOS_ASSERT(cell_conn.getConnectivitySize(local_cell_id)==1);

      // Set the global id for this local id
      global_cell_ids(local_cell_id) = cell_conn.getConnectivity(local_cell_id)[0];
    }
  }

  return global_cell_ids;
}

/** Build a Kokkos array mapping local cells to global node IDs.
  */
PHX::View<const GlobalOrdinal**>
buildOwnedCellsToGlobalNodes(const panzer::ConnManager & node_conn)
{

  // calculate total number of local cells
  std::vector<std::string> element_blocks;
  node_conn.getElementBlockIds(element_blocks);

  // compute total owned cells and maximum nodes per cell
  size_t num_local_cells=0, max_num_nodes_per_cell=0;
  for (const auto & element_block : element_blocks) {

    // Get local cells for this block
    const auto & local_cell_ids = node_conn.getElementBlock(element_block);
    if ( local_cell_ids.size() == 0 )
      continue;

    // Safe bet that all cells have the same number of nodes in an element block
    const size_t num_nodes_per_cell = node_conn.getConnectivitySize(local_cell_ids[0]);

    num_local_cells += local_cell_ids.size();
    max_num_nodes_per_cell = std::max(max_num_nodes_per_cell, num_nodes_per_cell);
  }

  // Note that not all owned cells are fully made up of 'owned' nodes, that's how we build the connectivity in the first place
  PHX::View<GlobalOrdinal**> global_nodes("owned_cell_to_global_vertex",num_local_cells,max_num_nodes_per_cell);

  for (const auto & element_block : element_blocks) {
    const auto & local_cell_ids = node_conn.getElementBlock(element_block);

    for(const auto local_cell_id : local_cell_ids){
      const panzer::GlobalOrdinal * connectivity = node_conn.getConnectivity(local_cell_id);
      const int num_nodes_per_cell = node_conn.getConnectivitySize(local_cell_id);

      // This is a bit tricky, because GlobalOrdinal is >= 0, we can't fill empty space with a "void" value (e.g. -1)
      TEUCHOS_ASSERT(num_nodes_per_cell == max_num_nodes_per_cell);

      // Set the global node ids for
      for(int node=0;node<num_nodes_per_cell;++node)
        global_nodes(local_cell_id,node) = connectivity[node];
    }

  }

  return global_nodes;
}

Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,TpetraNodeType> >
buildNodeMap(const Teuchos::RCP<const Teuchos::Comm<int> > & comm,
             PHX::View<const GlobalOrdinal**> cells_to_nodes)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  /*

   This function converts

   cells_to_nodes(local cell, local node) = global node index

   to a map describing which global nodes are found on this process

   Note that we have to ensure that that the global nodes found on this process are unique.

   */

  typedef Tpetra::Map<panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType> map_type;

  // get locally unique global ids
  std::set<GlobalOrdinal> global_nodes;
  for(unsigned int i=0;i<cells_to_nodes.extent(0);i++)
    for(unsigned int j=0;j<cells_to_nodes.extent(1);j++)
      global_nodes.insert(cells_to_nodes(i,j));

  // build local vector contribution
  Kokkos::View<panzer::GlobalOrdinal*> node_ids("global_nodes",global_nodes.size());
  int i = 0;
  for(auto itr=global_nodes.begin();itr!=global_nodes.end();++itr,++i)
    node_ids(i) = *itr;

  return rcp(new map_type(-1,node_ids,0,comm));
}

/** Given a cell to node map in a Kokkos array, build the node
  * to cell map using a transpose operation.
  */
Teuchos::RCP<Tpetra::CrsMatrix<panzer::LocalOrdinal,panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType> >
buildNodeToCellMatrix(const Teuchos::RCP<const Teuchos::Comm<int> > & comm,
                      PHX::View<const panzer::GlobalOrdinal*> owned_cells,
                      PHX::View<const panzer::GlobalOrdinal**> owned_cells_to_nodes)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef Tpetra::Map<panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType> map_type;
  typedef Tpetra::CrsMatrix<panzer::LocalOrdinal,panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType> crs_type;
  typedef Tpetra::Import<panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType> import_type;

  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::buildNodeToCellMatrix",BNTCM);

  TEUCHOS_ASSERT(owned_cells.extent(0)==owned_cells_to_nodes.extent(0));

  // build a unque node map to use with fill complete

  // This map identifies all nodes linked to cells on this process
  auto node_map = buildNodeMap(comm,owned_cells_to_nodes);

  // This map identifies nodes owned by this process
  auto unique_node_map  = Tpetra::createOneToOne(node_map);

  // This map identifies the cells owned by this process
  RCP<map_type> cell_map = rcp(new map_type(-1,owned_cells,0,comm));

  // Create a CRS matrix that stores a pointless value for every global node that belongs to a global cell
  // This is essentially another way to store cells_to_nodes
  RCP<crs_type> cell_to_node;
  {
    PANZER_FUNC_TIME_MONITOR_DIFF("Build matrix",BuildMatrix);
    // The matrix is indexed by (global cell, global node) = local node
    cell_to_node = rcp(new crs_type(cell_map,0));

    // fill in the cell to node matrix
    const unsigned int num_local_cells = owned_cells_to_nodes.extent(0);
    const unsigned int num_nodes_per_cell = owned_cells_to_nodes.extent(1);
    std::vector<LocalOrdinal> local_node_indexes(num_nodes_per_cell);
    std::vector<GlobalOrdinal> global_node_indexes(num_nodes_per_cell);
    for(unsigned int i=0;i<num_local_cells;i++) {
      const panzer::GlobalOrdinal global_cell_index = owned_cells(i);
      for(unsigned int j=0;j<num_nodes_per_cell;j++) {
        local_node_indexes[j] = Teuchos::as<panzer::LocalOrdinal>(j);
        global_node_indexes[j] = owned_cells_to_nodes(i,j);
      }

  //    cell_to_node_mat->insertGlobalValues(cells(i),cols,vals);
      cell_to_node->insertGlobalValues(global_cell_index,global_node_indexes,local_node_indexes);
    }
    cell_to_node->fillComplete(unique_node_map,cell_map);

  }

  // Transpose the cell_to_node array to create the node_to_cell array
  RCP<crs_type> node_to_cell;
  {
    PANZER_FUNC_TIME_MONITOR_DIFF("Tranpose matrix",TransposeMatrix);
    // Create an object designed to transpose the (global cell, global node) matrix to give
    // a (global node, global cell) matrix
    Tpetra::RowMatrixTransposer<panzer::LocalOrdinal,panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType> transposer(cell_to_node);

    // Create the transpose crs matrix
    auto trans = transposer.createTranspose();

    // We want to import the portion of the transposed matrix relating to all nodes on this process
    // The importer must import nodes required by this process (node_map) from the unique nodes (nodes living on a process)
    RCP<import_type> import = rcp(new import_type(unique_node_map,node_map));

    // Create the crs matrix to store (ghost node, global cell) array
    // This CRS matrix will have overlapping rows for shared global nodes
    node_to_cell = rcp(new crs_type(node_map,0));

    // Transfer data from the transpose array (associated with unique_node_map) to node_to_cell (associated with node_map)
    node_to_cell->doImport(*trans,*import,Tpetra::INSERT);

    // Set the fill - basicially locks down the structured of the CRS matrix - required before doing some operations
    //node_to_cell->fillComplete();
    node_to_cell->fillComplete(cell_map,unique_node_map);
  }

  // Return the node to cell array
  return node_to_cell;
}

/** Build ghstd cell one ring based on shared nodes
  */
Kokkos::View<GlobalOrdinal*>
buildGhostedCellOneRing(const Teuchos::RCP<const Teuchos::Comm<int> > & comm,
                        PHX::View<const panzer::GlobalOrdinal*> owned_cells,
                        PHX::View<const panzer::GlobalOrdinal**> owned_cells_to_global_nodes)
{

  PANZER_FUNC_TIME_MONITOR_DIFF("panzer_stk::buildGhostedCellOneRing",BGCOR);
  typedef Tpetra::CrsMatrix<int,int,panzer::GlobalOrdinal,panzer::TpetraNodeType> crs_type;

  // cells : (local cell index) -> global cell index
  // cells_to_nodes : (local cell index, local node_index) -> global node index

  /*

   This function creates a list of global indexes relating to ghost cells

   It is a little misleading how it does this, but the idea is to use the indexing of a CRS matrix to identify
   what global cells are connected to which global nodes. The values in the CRS matrix are meaningless, however,
   the fact that they are filled allows us to ping what index combinations exist.

   To do this we are going to use cell 'nodes' which could also be cell 'vertices'. It is unclear.

   First we construct an array that stores that global node indexes make up the connectivity for a given global cell index (order doesn't matter)

   cell_to_node : cell_to_node[global cell index][global node index] = some value (not important, just has to exist)

   We are then going to transpose that array

   node_to_cell : node_to_cell[global node index][global cell index] = some value (not important, just has to exist)

   The coloring of the node_to_cell array tells us what global cells are connected to a given global node.


   */

  // the node to cell matrix: Row = Global Node Id, Cell = Global Cell Id, Value = Cell Local Node Id
  Teuchos::RCP<crs_type> node_to_cell = buildNodeToCellMatrix(comm,owned_cells,owned_cells_to_global_nodes);

  // the set of cells already known
  std::unordered_set<panzer::GlobalOrdinal> unique_cells;

  // mark all owned cells as already known, e.g. and not in the list of
  // ghstd cells to be constructed

  for(size_t i=0;i<owned_cells.extent(0);i++)
    unique_cells.insert(owned_cells(i));

  // The set of ghost cells that share a global node with an owned cell
  std::set<GlobalOrdinal> ghstd_cells_set;

  // Get a list of global node indexes associated with the cells owned by this process
  auto node_map = node_to_cell->getMap()->getMyGlobalIndices();

  // Iterate through the global node indexes associated with this process
  for(size_t i=0;i<node_map.extent(0);i++) {
    const panzer::GlobalOrdinal global_node_index = node_map(i);
    size_t numEntries = node_to_cell->getNumEntriesInGlobalRow(node_map(i));
    Teuchos::Array<panzer::GlobalOrdinal> indices(numEntries);
    Teuchos::Array<int> values(numEntries);

    // Copy the row for a global node index into a local vector
    node_to_cell->getGlobalRowCopy(global_node_index,indices,values,numEntries);

    for(auto index : indices) {
      // if this is a new index (not owned, not previously found ghstd index
      // add it to the list of ghstd cells
      if(unique_cells.find(index)==unique_cells.end()) {
        ghstd_cells_set.insert(index);
        unique_cells.insert(index); // make sure you don't find it again
      }
    }
  }

  // build an array containing only the ghstd cells
  int indx = 0;
  Kokkos::View<panzer::GlobalOrdinal*> ghstd_cells("ghstd_cells",ghstd_cells_set.size());
  for(auto global_cell_index : ghstd_cells_set) {
    ghstd_cells(indx) = global_cell_index;
    indx++;
  }

  return ghstd_cells;
}

/**
 * Shares vertices based on global ids of cells
 */
PHX::View<const double ***>
buildGhostCellVertices(Teuchos::RCP<const Teuchos::Comm<int> > comm,
                       PHX::View<const GlobalOrdinal*> owned_cell_global_ids,
                       PHX::View<const GlobalOrdinal*> ghost_cell_global_ids,
                       PHX::View<const double ***> owned_cell_vertices)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef Tpetra::Map<panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType> map_type;
  typedef Tpetra::Import<panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType> import_type;
  typedef Tpetra::MultiVector<double,int,panzer::GlobalOrdinal,panzer::TpetraNodeType> mvec_type;
  typedef typename mvec_type::dual_view_type dual_view_type;

  RCP<map_type> owned_cell_map = rcp(new map_type(-1,owned_cell_global_ids, 0,comm));
  RCP<map_type> ghstd_cell_map = rcp(new map_type(-1,ghost_cell_global_ids,0,comm));

  // build importer: cell importer, owned to ghstd
  RCP<import_type> importer = rcp(new import_type(owned_cell_map,ghstd_cell_map));

  const size_t owned_cell_cnt = importer->getSourceMap()->getNodeNumElements();
  const size_t ghstd_cell_cnt = importer->getTargetMap()->getNodeNumElements();
  const int vertices_per_cell = owned_cell_vertices.extent(1);
  const int space_dim         = owned_cell_vertices.extent(2);

  TEUCHOS_ASSERT(owned_cell_vertices.extent(0)==owned_cell_cnt);

  // build vertex multivector
  RCP<mvec_type> owned_vertices_mv   = rcp(new mvec_type(importer->getSourceMap(),vertices_per_cell*space_dim));
  RCP<mvec_type> ghstd_vertices_mv = rcp(new mvec_type(importer->getTargetMap(),vertices_per_cell*space_dim));

  {
    auto owned_vertices_view = owned_vertices_mv->template getLocalView<dual_view_type>();
    for(size_t i=0;i<owned_cell_cnt;i++) {
      int l = 0;
      for(int j=0;j<vertices_per_cell;j++)
        for(int k=0;k<space_dim;k++,l++)
          owned_vertices_view(i,l) = owned_cell_vertices(i,j,k);
    }
  }

  // communicate ghstd vertices
  ghstd_vertices_mv->doImport(*owned_vertices_mv,*importer,Tpetra::INSERT);

  // copy multivector into ghstd vertex structure
  PHX::View<double ***> ghost_cell_vertices("ghost_cell_vertices",ghstd_cell_cnt,vertices_per_cell,space_dim);
  {
    auto ghstd_vertices_view = ghstd_vertices_mv->template getLocalView<dual_view_type>();
    for(size_t i=0;i<ghstd_cell_cnt;i++) {
      int l = 0;
      for(int j=0;j<vertices_per_cell;j++)
        for(int k=0;k<space_dim;k++,l++)
          ghost_cell_vertices(i,j,k) = ghstd_vertices_view(i,l);
    }
  }

  return ghost_cell_vertices;
} // end build ghstd vertices


PHX::View<const int*>
buildCellOwnership(const Teuchos::RCP<const Teuchos::Comm<int> > & comm,
                   const panzer::ConnManager & cell_conn,
                   const std::map<std::string, ElementSets::SetIndex> & block_id_map,
                   PHX::View<const GlobalOrdinal*> owned_cell_global_ids,
                   PHX::View<const GlobalOrdinal*> ghost_cell_global_ids,
                   PHX::View<const GlobalOrdinal*> virtual_cell_global_ids)
{

  // TODO: This call assumes that all ranks see the same set of element blocks.

  const size_t num_cells = owned_cell_global_ids.extent(0)
                         + ghost_cell_global_ids.extent(0)
                         + virtual_cell_global_ids.extent(0);

  PHX::View<ElementSets::SetIndex*> cell_ownership("cell_ownership",num_cells);

  // Get element blocks
  std::vector<std::string> element_blocks;
  cell_conn.getElementBlockIds(element_blocks);

  for(const auto & element_block : element_blocks){
    // Make sure our inputs are correct
    TEUCHOS_ASSERT(block_id_map.find(element_block) != block_id_map.end());
    const int block_id = block_id_map.at(element_block);
    for(const auto & local_id : cell_conn.getElementBlock(element_block))
      cell_ownership(local_id) = block_id;
  }

  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,TpetraNodeType> map_type;
  typedef Tpetra::Import<LocalOrdinal,GlobalOrdinal,TpetraNodeType> import_type;
  typedef Tpetra::MultiVector<int,int,GlobalOrdinal,panzer::TpetraNodeType> mvec_type;
  typedef typename mvec_type::dual_view_type dual_view_type;

  Teuchos::RCP<map_type> owned_cell_map = Teuchos::rcp(new map_type(-1,owned_cell_global_ids,0,comm));
  Teuchos::RCP<map_type> ghstd_cell_map = Teuchos::rcp(new map_type(-1,ghost_cell_global_ids,0,comm));

  // build importer: cell importer, owned to ghstd
  Teuchos::RCP<import_type> importer = rcp(new import_type(owned_cell_map,ghstd_cell_map));

  const size_t owned_cell_cnt = importer->getSourceMap()->getNodeNumElements();
  const size_t ghost_cell_cnt = importer->getTargetMap()->getNodeNumElements();

  TEUCHOS_ASSERT(owned_cell_cnt == owned_cell_global_ids.extent(0));
  TEUCHOS_ASSERT(ghost_cell_cnt == ghost_cell_global_ids.extent(0));

  // build vertex multivector
  Teuchos::RCP<mvec_type> owned_mv = Teuchos::rcp(new mvec_type(importer->getSourceMap(),1));
  Teuchos::RCP<mvec_type> ghost_mv = Teuchos::rcp(new mvec_type(importer->getTargetMap(),1));

  {
    auto cell_owner = owned_mv->template getLocalView<dual_view_type>();

    for(size_t i=0; i<owned_cell_cnt; ++i)
      cell_owner(i,0) = cell_ownership(i);
  }

  // communicate ghstd vertices
  ghost_mv->doImport(*owned_mv,*importer,Tpetra::INSERT);

  // Read ghost ownership
  {
    auto ghost_owner_view = ghost_mv->template getLocalView<dual_view_type>();
    for(size_t i=0;i<ghost_cell_cnt;i++)
      cell_ownership(owned_cell_cnt + i) = ghost_owner_view(i,0);
  }

  // Lastly we add in the virtual cells
  {
    const size_t num_real_cells = owned_cell_cnt + ghost_cell_cnt;
    const auto virtual_block_index = block_id_map.at("<<virtual>>");
    for(size_t i=0; i<virtual_cell_global_ids.extent(0); ++i)
      cell_ownership(num_real_cells + i) = virtual_block_index;
  }

  return cell_ownership;
}


PHX::View<const GlobalOrdinal*>
buildVirtualCellsGlobalIDs(const Teuchos::RCP<const Teuchos::Comm<int> > & comm,
                          const size_t num_owned_cells,
                          const size_t num_virtual_cells)
{

  PANZER_FUNC_TIME_MONITOR_DIFF("Initial global index creation",InitialGlobalIndexCreation);

  PHX::View<GlobalOrdinal*> virtual_cells("virtual_cells_global_ids",num_virtual_cells);

  const int num_ranks = comm->getSize();
  const int rank = comm->getRank();

  std::vector<GlobalOrdinal> owned_cell_distribution(num_ranks,0);
  {
    std::vector<GlobalOrdinal> my_owned_cell_distribution(num_ranks,0);
    my_owned_cell_distribution[rank] = num_owned_cells;

    Teuchos::reduceAll(*comm,Teuchos::REDUCE_SUM, num_ranks, my_owned_cell_distribution.data(),owned_cell_distribution.data());
  }

  std::vector<GlobalOrdinal> virtual_cell_distribution(num_ranks,0);
  {
    std::vector<GlobalOrdinal> my_virtual_cell_distribution(num_ranks,0);
    my_virtual_cell_distribution[rank] = num_virtual_cells;

    Teuchos::reduceAll(*comm,Teuchos::REDUCE_SUM, num_ranks, my_virtual_cell_distribution.data(),virtual_cell_distribution.data());
  }

  GlobalOrdinal num_global_real_cells=0;
  for(int i=0;i<num_ranks;++i)
    num_global_real_cells+=owned_cell_distribution[i];

  GlobalOrdinal global_virtual_start_idx = num_global_real_cells;
  for(int i=0;i<rank;++i)
    global_virtual_start_idx += virtual_cell_distribution[i];

  for(size_t i=0;i<num_virtual_cells;++i)
    virtual_cells(i) = global_virtual_start_idx + GlobalOrdinal(i);

  return virtual_cells;
}


} // namespace


void
LocalMeshInfo::
initialize(const Teuchos::RCP<const Teuchos::Comm<int> > & comm,
           ConnManager & conn,
           const PHX::View<const double ***> & owned_cell_vertices)
{

  TEUCHOS_FUNC_TIME_MONITOR_DIFF("LocalMeshInfo::LocalMeshInfo",GhostedMeshInfo);

  // Names of element blocks
  std::vector<std::string> element_blocks;
  conn.getElementBlockIds(element_blocks);

  // Before we begin, we need to construct two connectivity objects
  Teuchos::RCP<ConnManager> cell_conn_rcp = conn.noConnectivityClone(),
                            node_conn_rcp = conn.noConnectivityClone();

  std::map<std::string, Teuchos::RCP<const shards::CellTopology>> cell_topology_map;
  {
    std::vector<shards::CellTopology> element_block_topologies;
    conn.getElementBlockTopologies(element_block_topologies);

    // Sanity check
    TEUCHOS_ASSERT(element_block_topologies.size() == element_blocks.size());

    // We will need this map later
    for(unsigned int i=0; i<element_blocks.size(); ++i)
      cell_topology_map[element_blocks[i]] = Teuchos::rcp(new shards::CellTopology(element_block_topologies[i]));

    const shards::CellTopology & topology = element_block_topologies[0];

    // FIXME: We assume that all element blocks have the same topology.
    for(const auto & other_topology : element_block_topologies){
      TEUCHOS_ASSERT(other_topology.getKey() == topology.getKey());
    }

    Teuchos::RCP<panzer::FieldPattern> cell_pattern;
    if(topology.getDimension() == 1){
      cell_pattern = Teuchos::rcp(new EdgeFieldPattern(topology));
    } else if(topology.getDimension() == 2){
      cell_pattern = Teuchos::rcp(new FaceFieldPattern(topology));
    } else if(topology.getDimension() == 3){
      cell_pattern = Teuchos::rcp(new ElemFieldPattern(topology));
    }

    cell_conn_rcp->buildConnectivity(*cell_pattern);
    node_conn_rcp->buildConnectivity(NodalFieldPattern(topology));
  }

  const ConnManager & cell_conn = * cell_conn_rcp;
  const ConnManager & node_conn = * node_conn_rcp;


  // Build the global ids for local cells
  const auto owned_cells_global_ids = buildOwnedCellsGlobalIDs(cell_conn);
  num_owned_cells = owned_cells_global_ids.extent(0);

  // Sanity check
  TEUCHOS_ASSERT(owned_cell_vertices.extent_int(0) == num_owned_cells);

  // Build the global node ids for local cells
  const auto owned_cells_to_global_nodes = buildOwnedCellsToGlobalNodes(node_conn);

  // Build the global ids for ghost cells
  const auto ghost_cells_global_ids = buildGhostedCellOneRing(comm,owned_cells_global_ids,owned_cells_to_global_nodes);
  num_ghost_cells = ghost_cells_global_ids.extent(0);

  // Build the ghost cell vertices
  const auto ghost_cell_vertices = buildGhostCellVertices(comm,owned_cells_global_ids,ghost_cells_global_ids,owned_cell_vertices);

  // Virtual cells require two additional components

  // Face connectivities allow us to find faces at the boundary of the domain
  FaceToElement<LocalOrdinal,GlobalOrdinal> face_connectivity;
  face_connectivity.initialize(conn);

  // The face connectivity is designed for GlobalOrdinal so we need a LocalOrdinal lookup
  std::unordered_map<GlobalOrdinal,LocalOrdinal> global_ids_to_local_ids;
  for(size_t i=0; i<owned_cells_global_ids.extent(0); ++i)
    global_ids_to_local_ids[owned_cells_global_ids(i)] = i;
  for(size_t i=0; i<ghost_cells_global_ids.extent(0); ++i)
    global_ids_to_local_ids[ghost_cells_global_ids(i)] = i;

  // Iterate through the faces and add up virtual cells
  std::vector<GlobalOrdinal> virtual_reflection_cells;
  {
    const auto elems_by_face = face_connectivity.getFaceToElementsMap();
    for(size_t face=0;face<elems_by_face.extent(0);++face){
      if(elems_by_face(face,0) < 0)
        virtual_reflection_cells.push_back(elems_by_face(face,1));
      else if(elems_by_face(face,1) < 0)
        virtual_reflection_cells.push_back(elems_by_face(face,0));
    }
  }
  num_virtual_cells = virtual_reflection_cells.size();

  // Build the virtual cell global ids
  PHX::View<GlobalOrdinal*> virtual_cells_global_ids("virtual_cells_global_ids",num_virtual_cells);
  {
    const int num_ranks = comm->getSize();
    const int rank = comm->getRank();

    std::vector<GlobalOrdinal> owned_cell_distribution(num_ranks,0);
    {
      std::vector<GlobalOrdinal> my_owned_cell_distribution(num_ranks,0);
      my_owned_cell_distribution[rank] = num_owned_cells;

      Teuchos::reduceAll(*comm,Teuchos::REDUCE_SUM, num_ranks, my_owned_cell_distribution.data(),owned_cell_distribution.data());
    }

    std::vector<GlobalOrdinal> virtual_cell_distribution(num_ranks,0);
    {
      std::vector<GlobalOrdinal> my_virtual_cell_distribution(num_ranks,0);
      my_virtual_cell_distribution[rank] = num_virtual_cells;

      Teuchos::reduceAll(*comm,Teuchos::REDUCE_SUM, num_ranks, my_virtual_cell_distribution.data(),virtual_cell_distribution.data());
    }

    GlobalOrdinal num_global_real_cells=0;
    for(int i=0;i<num_ranks;++i)
      num_global_real_cells+=owned_cell_distribution[i];

    GlobalOrdinal global_virtual_start_idx = num_global_real_cells;
    for(int i=0;i<rank;++i)
      global_virtual_start_idx += virtual_cell_distribution[i];

    for(size_t i=0;i<num_virtual_cells;++i)
      virtual_cells_global_ids(i) = global_virtual_start_idx + GlobalOrdinal(i);

  }

  // It is time to allocate our arrays
  const LocalOrdinal num_real_cells = num_owned_cells + num_ghost_cells;
  const LocalOrdinal num_cells = num_real_cells + num_virtual_cells;
  const int num_vertices_per_cell = owned_cell_vertices.extent_int(1);
  const int num_dims_per_vertex = owned_cell_vertices.extent_int(2);

  global_cells = PHX::View<GlobalOrdinal*>("global_cell_ids",num_cells);
  cell_vertices = PHX::View<double ***>("cell_vertices",num_cells,num_vertices_per_cell,num_dims_per_vertex);

  // This will be used later for generating virtual cell geometry
  std::map<GlobalOrdinal,LocalOrdinal> owned_global_to_local;

  // Start with owned cells
  for(LocalOrdinal local_cell_id=0; local_cell_id<num_owned_cells; ++local_cell_id){

    const GlobalOrdinal global_cell_id = owned_cells_global_ids(local_cell_id);

    // Local cell ids are aligned for owned cells
    global_cells(local_cell_id) = global_cell_id;
    owned_global_to_local[global_cell_id] = local_cell_id;

    for(int v=0; v<num_vertices_per_cell; ++v)
      for(int d=0; d<num_dims_per_vertex; ++d)
        cell_vertices(local_cell_id,v,d) = owned_cell_vertices(local_cell_id,v,d);

  }

  // Ghost cells are a little more complicated
  {
    // Because the ConnManager has an incomplete local representation of the mesh (e.g. no ghost cells for periodic boundary conditions)
    // we have to do a partial lookup on all known global ids for 'neighbor' cells to make sure the local ids line up

    std::unordered_map<GlobalOrdinal, LocalOrdinal> existing_ghost_global_to_local;
    for(const auto & element_block : element_blocks){
      for(const auto & local_cell_id : cell_conn.getNeighborElementBlock(element_block)){

        // Sanity check
        TEUCHOS_ASSERT(cell_conn.getConnectivitySize(local_cell_id)==1);

        const auto global_cell_id = cell_conn.getConnectivity(local_cell_id)[0];

        existing_ghost_global_to_local[global_cell_id] = local_cell_id;
      }
    }

    LocalOrdinal new_ghost_index = num_owned_cells + existing_ghost_global_to_local.size();

    // We use this lookup table to properly orient our ghost cells
    for(size_t i=0; i<ghost_cells_global_ids.extent(0); ++i){

      const GlobalOrdinal global_cell_id = ghost_cells_global_ids(i);

      const auto itr = existing_ghost_global_to_local.find(global_cell_id);

      // Figure out the local cell id for this global cell id
      LocalOrdinal local_cell_id = -1;
      if(itr == existing_ghost_global_to_local.end())
        local_cell_id = new_ghost_index++;
      else
        local_cell_id = itr->second;

      // Local cell ids are aligned for owned cells
      global_cells(local_cell_id) = global_cell_id;

      for(int v=0; v<num_vertices_per_cell; ++v)
        for(int d=0; d<num_dims_per_vertex; ++d)
          cell_vertices(local_cell_id,v,d) = ghost_cell_vertices(i,v,d);

    }

    // Since we just reordered our ghost cells, we need to make sure ghost_cells_global_ids is correct
    // This will be used when we set ghost ownership
    for(size_t i=0; i<ghost_cells_global_ids.extent(0); ++i)
      ghost_cells_global_ids(i) = global_cells(num_owned_cells + i);

  }

  // Finally we have virtual cells
  {

    for(size_t i=0; i<num_virtual_cells; ++i){

      const auto global_cell_id = virtual_cells_global_ids(i);
      const auto local_cell_id = num_real_cells + i;

      // Local cell ids are aligned for owned cells
      global_cells(local_cell_id) = global_cell_id;

      // Note we just copy the owned cell geometry corresponding to the virtual cell
      // TODO: We need to make this an extrusion, but it will require a bit more work
      const auto reflected_global_cell_id = virtual_reflection_cells[i];
      const auto reflected_local_cell_id = owned_global_to_local.at(reflected_global_cell_id);

      for(int v=0; v<num_vertices_per_cell; ++v)
        for(int d=0; d<num_dims_per_vertex; ++d)
          cell_vertices(local_cell_id,v,d) = cell_vertices(reflected_local_cell_id,v,d);
    }

  }

  // Local cells are initialized here as a simple 1-N list, which will be used in partitioning later
  local_cells = Kokkos::View<LocalOrdinal*>("local_cell_ids",num_cells);
  for(LocalOrdinal i=0; i<num_cells; ++i)
    local_cells(i) = i;

  // Now that the global ids are set, we need to define ownership of the cells
  // i.e. what blocks own what cells
  // As with ghost/virtual cells, the ConnManager has an incomplete definition of what is owned where

  std::map<std::string, ElementSets::SetIndex> block_id_map;
  block_id_map["<<virtual>>"] = 0;
  for(unsigned int i=0; i<element_blocks.size(); ++i)
    block_id_map[element_blocks[i]] = i+1;

  auto cell_ownership = buildCellOwnership(comm,
                                           cell_conn,
                                           block_id_map,
                                           owned_cells_global_ids,
                                           ghost_cells_global_ids,
                                           virtual_cells_global_ids);

  cell_sets = Teuchos::rcp(new ElementSets(block_id_map,cell_ownership));

  ////////////////////////////////////////////////////////////////////////////////////
  // Face connectivity

  // Transfer information from 'faceToElement' datasets to local arrays
  {
    PANZER_FUNC_TIME_MONITOR_DIFF("Transer faceToElement to local",TransferFaceToElementLocal);

    has_connectivity = true;

    // To work with the face_connectivity object, we have to be able to link global ids to local ids
    // Note we only need this lookup table to cover real cells since FaceToElement doesn't know about virtual cells
    std::unordered_map<GlobalOrdinal,LocalOrdinal> global_to_local;
    for(LocalOrdinal local_cell_id=0; local_cell_id < num_real_cells; ++local_cell_id)
      global_to_local[global_cells(local_cell_id)] = local_cell_id;

    // TODO: This will eventually be converted into a SubcellConnectivity object
    // which will be generated directly from the ConnManager once the ConnManager
    // handles periodic boundary conditions properly

    const auto elems_by_face = face_connectivity.getFaceToElementsMap();
    const auto face_to_localidx  = face_connectivity.getFaceToCellLocalIdxMap();

    const size_t num_total_faces = elems_by_face.extent(0);

    // TODO: I'm not sure how to make this universal
    const size_t num_faces_per_cell = cell_topology_map.begin()->second->getSubcellCount(num_dims_per_vertex-1);

    // Arrays for storing face connectivity information
    face_to_cells = Kokkos::View<panzer::LocalOrdinal*[2]>("face_to_cells",num_total_faces);
    face_to_lidx = Kokkos::View<panzer::LocalOrdinal*[2]>("face_to_lidx",num_total_faces);
    cell_to_faces = Kokkos::View<panzer::LocalOrdinal**>("cell_to_face",num_cells,num_faces_per_cell);

    // initialize with negative one cells that are not associated with a face
    Kokkos::deep_copy(cell_to_faces,-1);

    // Note that this iteration is designed to
    int virtual_cell_index = num_real_cells;
    for(size_t f=0;f<elems_by_face.extent(0);f++) {

      const GlobalOrdinal global_c0 = elems_by_face(f,0);
      const GlobalOrdinal global_c1 = elems_by_face(f,1);

      // Check to see if the global ids are in our map
      const auto g0itr = global_to_local.find(global_c0);
      const auto g1itr = global_to_local.find(global_c1);

      // If the global id was found, then we use it, otherwise we assume this is a virtual cell
      const LocalOrdinal cell0 = (g0itr == global_to_local.end()) ? virtual_cell_index++ : g0itr->second;
      const LocalOrdinal cell1 = (g1itr == global_to_local.end()) ? virtual_cell_index++ : g1itr->second;

      // We need to make sure both cells are not virtual
      TEUCHOS_ASSERT(not ((cell0 >= num_real_cells) and (cell1 >= num_real_cells)));

      // The subcell indexes of the faces with respect to the elements
      int lidx0 = face_to_localidx(f,0);
      int lidx1 = face_to_localidx(f,1);

      // Test for virtual cells

      // Left cell
      if(cell0 >= num_real_cells){
        // We make the subcell indexes line up - TODO: Why???
        lidx0 = lidx1;
      }
      cell_to_faces(cell0,lidx0) = f;

      // Right cell
      if(cell1 >= num_real_cells){
        // We make the subcell indexes line up - TODO: Why???
        lidx1 = lidx0;
      }
      cell_to_faces(cell1,lidx1) = f;

      // Faces point from low cell index to high cell index (i.e. always from owned TO ghost and owned TO virtual)
      if(cell0<cell1){
        face_to_cells(f,0) = cell0;
        face_to_lidx(f,0) = lidx0;
        face_to_cells(f,1) = cell1;
        face_to_lidx(f,1) = lidx1;
      } else {
        face_to_cells(f,0) = cell1;
        face_to_lidx(f,0) = lidx1;
        face_to_cells(f,1) = cell0;
        face_to_lidx(f,1) = lidx0;
      }

    }
  }


  ////////////////////////////////////////////////////////////////////////////////////
  // At this point all the data structures have been built.
  // Let's store it in the mesh_info object
  {
    PANZER_FUNC_TIME_MONITOR_DIFF("Setup LocalMeshBlockInfo",SetupLocalMeshBlockInfo);

    // Map from local cell id to element index (not sure if this is necessary)
    std::unordered_map<size_t, int> local_id_map;
    for(size_t i=0; i<num_real_cells; ++i)
      local_id_map[local_cells(i)] = i;

    for(const auto & element_block : element_blocks){

      // Build local block info for element block - note these are local indexes which align with the indexing of LocalMeshInfo
      const auto & element_block_cells = cell_conn.getElementBlock(element_block);

      if(element_block_cells.size() == 0)
        continue;

      // We have to convert to the datatype used by the partitioner
      std::vector<size_t> owned_block_cells(element_block_cells.size());
      for(size_t i=0; i<element_block_cells.size(); ++i)
        owned_block_cells[i] = element_block_cells[i];

      // Create the block info object
      auto & block_info = this->element_blocks[element_block];
      block_info.num_owned_cells = owned_block_cells.size();
      block_info.element_block_name = element_block;
      block_info.cell_topology = cell_topology_map.at(element_block);
      {
        PANZER_FUNC_TIME_MONITOR("panzer::partitioning_utilities::setupSubLocalMeshInfoWithConnectivity");
        panzer::partitioning_utilities::setupSubLocalMeshInfoWithConnectivity(*this, owned_block_cells, block_info);
      }

      // Setup sidesets - REQUIRES STK
//      for(const auto & sideset : sidesets){
//
//        PANZER_FUNC_TIME_MONITOR_DIFF("Setup LocalMeshSidesetInfo",SetupLocalMeshSidesetInfo);
//        setupLocalMeshSidesetInfo(mesh, conn, mesh_info, element_block, sideset, mesh_info.sidesets);
//      }
    }
  }

}
}
