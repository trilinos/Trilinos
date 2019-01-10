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

#include "Panzer_NodeType.hpp"
#include "Panzer_STK_LocalMeshUtilities.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STKConnManager.hpp"

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

#include "Phalanx_KokkosDeviceTypes.hpp"

#include "Teuchos_Assert.hpp"

#include "Tpetra_Import.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"

#include <string>
#include <map>
#include <vector>
#include <unordered_set>

namespace panzer_stk
{

// No external access
namespace
{

/** Build a Kokkos array of all the global cell IDs from a connection manager.
  * Note that this is mapping between local IDs to global IDs.
  */
template <typename LO,typename GO>
void
buildCellGlobalIDs(panzer::ConnManager<LO,GO> & conn,
                   Kokkos::View<GO*> & globals)
{
  // extract topologies, and build global connectivity...currently assuming only one topology
  std::vector<shards::CellTopology> elementBlockTopologies;
  conn.getElementBlockTopologies(elementBlockTopologies);
  const shards::CellTopology & topology = elementBlockTopologies[0];

  // FIXME: We assume that all element blocks have the same topology.
  for(const auto & other_topology : elementBlockTopologies){
    TEUCHOS_ASSERT(other_topology.getKey() == topology.getKey());
  }

  Teuchos::RCP<panzer::FieldPattern> cell_pattern;
  if(topology.getDimension() == 1){
    cell_pattern = Teuchos::rcp(new panzer::EdgeFieldPattern(topology));
  } else if(topology.getDimension() == 2){
    cell_pattern = Teuchos::rcp(new panzer::FaceFieldPattern(topology));
  } else if(topology.getDimension() == 3){
    cell_pattern = Teuchos::rcp(new panzer::ElemFieldPattern(topology));
  }

//  panzer::EdgeFieldPattern cell_pattern(elementBlockTopologies[0]);
  conn.buildConnectivity(*cell_pattern);

  // calculate total number of local cells
  std::vector<std::string> block_ids;
  conn.getElementBlockIds(block_ids);

  std::size_t totalSize = 0;
  for (std::size_t which_blk=0;which_blk<block_ids.size();which_blk++) {
    // get the elem to face mapping
    const std::vector<LO> & localIDs = conn.getElementBlock(block_ids[which_blk]);
    totalSize += localIDs.size();
  }
  globals = Kokkos::View<GO*>("global_cells",totalSize);

  for (std::size_t id=0;id<totalSize; ++id) {
    // sanity check
    int n_conn = conn.getConnectivitySize(id);
    TEUCHOS_ASSERT(n_conn==1);

    const GO * connectivity = conn.getConnectivity(id);
    globals(id) = connectivity[0];
  }

//  print_view_1D("buildCellGlobalIDs : globals",globals);
}

/** Build a Kokkos array mapping local cells to global node IDs.
  * Note that these are 'vertex nodes' and not 'basis nodes', 'quad nodes', or 'dof nodes'
  */
template <typename LO,typename GO>
void
buildCellToNodes(panzer::ConnManager<LO,GO> & conn, Kokkos::View<GO**> & globals)
{
  // extract topologies, and build global connectivity...currently assuming only one topology
  std::vector<shards::CellTopology> elementBlockTopologies;
  conn.getElementBlockTopologies(elementBlockTopologies);
  const shards::CellTopology & topology = elementBlockTopologies[0];

  // FIXME: We assume that all element blocks have the same topology.
  for(const auto & other_topology : elementBlockTopologies){
    TEUCHOS_ASSERT(other_topology.getKey() == topology.getKey());
  }

  panzer::NodalFieldPattern pattern(topology);
  conn.buildConnectivity(pattern);

  // calculate total number of local cells
  std::vector<std::string> block_ids;
  conn.getElementBlockIds(block_ids);

  // compute total cells and maximum nodes
  std::size_t totalCells=0, maxNodes=0;
  for (std::size_t which_blk=0;which_blk<block_ids.size();which_blk++) {
    // get the elem to face mapping
    const std::vector<LO> & localIDs = conn.getElementBlock(block_ids[which_blk]);
    if ( localIDs.size() == 0 )
      continue;
    LO thisSize = conn.getConnectivitySize(localIDs[0]);

    totalCells += localIDs.size();
    maxNodes = maxNodes<Teuchos::as<std::size_t>(thisSize) ? Teuchos::as<std::size_t>(thisSize) : maxNodes;
  }
  globals = Kokkos::View<GO**>("cell_to_node",totalCells,maxNodes);

  // build connectivity array
  for (std::size_t id=0;id<totalCells; ++id) {
    const GO * connectivity = conn.getConnectivity(id);
    LO nodeCnt = conn.getConnectivitySize(id);

    for(int n=0;n<nodeCnt;n++)
      globals(id,n) = connectivity[n];
  }

//  print_view("buildCellToNodes : globals",globals);
}

template <typename LO,typename GO>
Teuchos::RCP<const Tpetra::Map<LO,GO,panzer::TpetraNodeType> >
buildNodeMap(const Teuchos::RCP<const Teuchos::Comm<int> > & comm,
                    Kokkos::View<const GO**> cells_to_nodes)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  /*

   This function converts

   cells_to_nodes(local cell, local node) = global node index

   to a map describing which global nodes are found on this process

   Note that we have to ensure that that the global nodes found on this process are unique.

   */

  typedef Tpetra::Map<LO,GO,panzer::TpetraNodeType> map_type;

  // get locally unique global ids
  std::set<GO> global_nodes;
  for(unsigned int i=0;i<cells_to_nodes.extent(0);i++)
    for(unsigned int j=0;j<cells_to_nodes.extent(1);j++)
      global_nodes.insert(cells_to_nodes(i,j));

  // build local vector contribution
  Kokkos::View<GO*> node_ids("global_nodes",global_nodes.size());
  int i = 0;
  for(auto itr=global_nodes.begin();itr!=global_nodes.end();++itr,++i)
    node_ids(i) = *itr;

//  print_view("buildNodeMap : cells_to_nodes",cells_to_nodes);
//  print_view_1D("buildNodeMap : node_ids",node_ids);

  return rcp(new map_type(-1,node_ids,0,comm));
}

/** Given a cell to node map in a Kokkos array, build the node
  * to cell map using a transpose operation.
  */
template <typename LO,typename GO>
Teuchos::RCP<Tpetra::CrsMatrix<LO,LO,GO,panzer::TpetraNodeType> >
buildNodeToCellMatrix(const Teuchos::RCP<const Teuchos::Comm<int> > & comm,
                      Kokkos::View<const GO*> owned_cells,
                      Kokkos::View<const GO**> owned_cells_to_nodes)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef Tpetra::Map<LO,GO,panzer::TpetraNodeType> map_type;
  typedef Tpetra::CrsMatrix<LO,LO,GO,panzer::TpetraNodeType> crs_type;
  typedef Tpetra::Import<LO,GO,panzer::TpetraNodeType> import_type;


  PANZER_FUNC_TIME_MONITOR_DIFF("panzer_stk::buildNodeToCellMatrix",BNTCM);

  TEUCHOS_ASSERT(owned_cells.extent(0)==owned_cells_to_nodes.extent(0));

  // build a unque node map to use with fill complete

  // This map identifies all nodes linked to cells on this process
  auto node_map = buildNodeMap<LO,GO>(comm,owned_cells_to_nodes);

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
    std::vector<LO> local_node_indexes(num_nodes_per_cell);
    std::vector<GO> global_node_indexes(num_nodes_per_cell);
    for(unsigned int i=0;i<num_local_cells;i++) {
      const GO global_cell_index = owned_cells(i);
  //    std::vector<LO> vals(cells_to_nodes.extent(1));
  //    std::vector<GO> cols(cells_to_nodes.extent(1));
      for(unsigned int j=0;j<num_nodes_per_cell;j++) {
  //      vals[j] = Teuchos::as<LO>(j);
  //      cols[j] = cells_to_nodes(i,j);
        local_node_indexes[j] = Teuchos::as<LO>(j);
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
    Tpetra::RowMatrixTransposer<LO,LO,GO,panzer::TpetraNodeType> transposer(cell_to_node);

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
template <typename GO>
Kokkos::View<const GO*>
buildGhostedCellOneRing(const Teuchos::RCP<const Teuchos::Comm<int> > & comm,
                        Kokkos::View<const GO*> cells,
                        Kokkos::View<const GO**> cells_to_nodes)
{

  PANZER_FUNC_TIME_MONITOR_DIFF("panzer_stk::buildGhostedCellOneRing",BGCOR);
  typedef Tpetra::CrsMatrix<int,int,GO,panzer::TpetraNodeType> crs_type;

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
  Teuchos::RCP<crs_type> node_to_cell = buildNodeToCellMatrix<int,GO>(comm,cells,cells_to_nodes);

  // the set of cells already known
  std::unordered_set<GO> unique_cells;

  // mark all owned cells as already known, e.g. and not in the list of
  // ghstd cells to be constructed
  for(size_t i=0;i<cells.extent(0);i++) {
    unique_cells.insert(cells(i));
  }

  // The set of ghost cells that share a global node with an owned cell
  std::set<GO> ghstd_cells_set;

  // Get a list of global node indexes associated with the cells owned by this process
//  auto node_map = node_to_cell->getRangeMap()->getMyGlobalIndices();
  auto node_map = node_to_cell->getMap()->getMyGlobalIndices();

  // Iterate through the global node indexes associated with this process
  for(size_t i=0;i<node_map.extent(0);i++) {
    const GO global_node_index = node_map(i);
    size_t numEntries = node_to_cell->getNumEntriesInGlobalRow(node_map(i));
    Teuchos::Array<GO> indices(numEntries);
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
  Kokkos::View<GO*> ghstd_cells("ghstd_cells",ghstd_cells_set.size());
  for(auto global_cell_index : ghstd_cells_set) {
    ghstd_cells(indx) = global_cell_index;
    indx++;
  }

//  print_view_1D("ghstd_cells",ghstd_cells);

  return ghstd_cells;
}

/** This method takes a cell importer (owned to ghstd) and communicates vertices
  * of the ghstd elements.
  */
template <typename GO>
Kokkos::DynRankView<double,PHX::Device>
buildGhostedVertices(const Tpetra::Import<int,GO,panzer::TpetraNodeType> & importer,
                     Kokkos::DynRankView<const double,PHX::Device> owned_vertices)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef Tpetra::MultiVector<double,int,GO,panzer::TpetraNodeType> mvec_type;
  typedef typename mvec_type::dual_view_type dual_view_type;

  size_t owned_cell_cnt = importer.getSourceMap()->getNodeNumElements();
  size_t ghstd_cell_cnt = importer.getTargetMap()->getNodeNumElements();
  int vertices_per_cell = owned_vertices.extent(1);
  int space_dim         = owned_vertices.extent(2);

  TEUCHOS_ASSERT(owned_vertices.extent(0)==owned_cell_cnt);

  // build vertex multivector
  RCP<mvec_type> owned_vertices_mv   = rcp(new mvec_type(importer.getSourceMap(),vertices_per_cell*space_dim));
  RCP<mvec_type> ghstd_vertices_mv = rcp(new mvec_type(importer.getTargetMap(),vertices_per_cell*space_dim));

  {
    auto owned_vertices_view = owned_vertices_mv->template getLocalView<dual_view_type>();
    for(size_t i=0;i<owned_cell_cnt;i++) {
      int l = 0;
      for(int j=0;j<vertices_per_cell;j++)
        for(int k=0;k<space_dim;k++,l++)
          owned_vertices_view(i,l) = owned_vertices(i,j,k);
    }
  }

  // communicate ghstd vertices
  ghstd_vertices_mv->doImport(*owned_vertices_mv,importer,Tpetra::INSERT);

  // copy multivector into ghstd vertex structure
  Kokkos::DynRankView<double,PHX::Device> ghstd_vertices("ghstd_vertices",ghstd_cell_cnt,vertices_per_cell,space_dim);
  {
    auto ghstd_vertices_view = ghstd_vertices_mv->template getLocalView<dual_view_type>();
    for(size_t i=0;i<ghstd_cell_cnt;i++) {
      int l = 0;
      for(int j=0;j<vertices_per_cell;j++)
        for(int k=0;k<space_dim;k++,l++)
          ghstd_vertices(i,j,k) = ghstd_vertices_view(i,l);
    }
  }

  return ghstd_vertices;
} // end build ghstd vertices

template<typename LO, typename GO>
void
setupLocalMeshBlockInfo(const panzer_stk::STK_Interface & mesh,
                        panzer::ConnManager<LO,GO> & conn,
                        const panzer::LocalMeshInfo<LO,GO> & mesh_info,
                        const std::string & element_block_name,
                        panzer::LocalMeshBlockInfo<LO,GO> & block_info)
{

  // This function identifies all cells in mesh_info that belong to element_block_name
  // and creates a block_info from it.

  const int num_parent_owned_cells = mesh_info.num_owned_cells;

  // Make sure connectivity is setup for interfaces between cells
  {
    const shards::CellTopology & topology = *(mesh.getCellTopology(element_block_name));
    Teuchos::RCP<panzer::FieldPattern> cell_pattern;
    if(topology.getDimension() == 1){
      cell_pattern = Teuchos::rcp(new panzer::EdgeFieldPattern(topology));
    } else if(topology.getDimension() == 2){
      cell_pattern = Teuchos::rcp(new panzer::FaceFieldPattern(topology));
    } else if(topology.getDimension() == 3){
      cell_pattern = Teuchos::rcp(new panzer::ElemFieldPattern(topology));
    }

    {
      PANZER_FUNC_TIME_MONITOR("Build connectivity");
      conn.buildConnectivity(*cell_pattern);
    }
  }

  std::vector<LO> owned_block_cells;
  for(int parent_owned_cell=0;parent_owned_cell<num_parent_owned_cells;++parent_owned_cell){
    const LO local_cell = mesh_info.local_cells(parent_owned_cell);
    const bool is_in_block = conn.getBlockId(local_cell) == element_block_name;

    if(is_in_block){
      owned_block_cells.push_back(parent_owned_cell);
    }

  }

  if ( owned_block_cells.size() == 0 )
    return;
  block_info.num_owned_cells = owned_block_cells.size();
  block_info.element_block_name = element_block_name;
  block_info.cell_topology = mesh.getCellTopology(element_block_name);
  {
    PANZER_FUNC_TIME_MONITOR("panzer::partitioning_utilities::setupSubLocalMeshInfo");
    panzer::partitioning_utilities::setupSubLocalMeshInfo<LO,GO>(mesh_info, owned_block_cells, block_info);
  }
}


template<typename LO, typename GO>
void
setupLocalMeshSidesetInfo(const panzer_stk::STK_Interface & mesh,
                          panzer::ConnManager<LO,GO>& /* conn */,
                          const panzer::LocalMeshInfo<LO,GO> & mesh_info,
                          const std::string & element_block_name,
                          const std::string & sideset_name,
                          panzer::LocalMeshSidesetInfo<LO,GO> & sideset_info)
{

  // This function identifies all cells in mesh_info that belong to element_block_name
  // and creates a block_info from it.

  const size_t face_subcell_dimension = mesh.getCellTopology(element_block_name)->getDimension()-1;

  Kokkos::DynRankView<double,PHX::Device> data_allocation;

  // This is a list of all entities that make up the 'side'
  std::vector<stk::mesh::Entity> side_entities;

  // Grab the side entities associated with this sideset on the mesh
  // Note: Throws exception if element block or sideset doesn't exist
  try{

    mesh.getMySides(sideset_name, element_block_name, side_entities);

  } catch(STK_Interface::SidesetException & e) {
     std::stringstream ss;
     std::vector<std::string> sideset_names;
     mesh.getSidesetNames(sideset_names);

     // build an error message
     ss << e.what() << "\nChoose existing sideset:\n";
     for(const auto & optional_sideset_name : sideset_names){
        ss << "\t\"" << optional_sideset_name << "\"\n";
     }

     TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(true,std::logic_error,ss.str());

  } catch(STK_Interface::ElementBlockException & e) {
     std::stringstream ss;
     std::vector<std::string> element_block_names;
     mesh.getElementBlockNames(element_block_names);

     // build an error message
     ss << e.what() << "\nChoose existing element block:\n";
     for(const auto & optional_element_block_name : element_block_names){
        ss << "\t\"" << optional_element_block_name << "\"\n";
     }

     TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(true,std::logic_error,ss.str());

  } catch(std::logic_error & e) {
     std::stringstream ss;
     ss << e.what() << "\nUnrecognized logic error.\n";

     TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(true,std::logic_error,ss.str());

  }

  // We now have a list of sideset entities, lets unwrap them and create the sideset_info!

  // (subcell_dimension, subcell_index) -> list of cells
  std::map<std::pair<size_t,size_t>, std::vector<size_t> > local_cell_indexes_by_subcell;

  // List of local subcell indexes indexed by element:
  // For example: a Tet (element) would have
  //  - 4 triangular faces (subcell_index 0-3, subcell_dimension=2)
  //  - 6 edges (subcell_index 0-5, subcell_dimension=1)
  //  - 4 vertices (subcell_index 0-3, subcell_dimension=0)
  // Another example: a Line (element) would have
  //  - 2 vertices (subcell_index 0-1, subcell_dimension=0)
  std::vector<stk::mesh::Entity> elements;
  std::vector<size_t> subcell_indexes;
  std::vector<size_t> subcell_dimensions;
  panzer_stk::workset_utils::getSideElementCascade(mesh, element_block_name, side_entities, subcell_dimensions, subcell_indexes, elements);
  const LO num_elements = subcell_dimensions.size();

  // build local cell_ids, mapped by local side id
  for(LO element_index=0; element_index<num_elements; ++element_index) {
    const size_t subcell_dimension = subcell_dimensions[element_index];
    const size_t subcell_index = subcell_indexes[element_index];
    const size_t element_local_index = mesh.elementLocalId(elements[element_index]);

    // Add subcell to map
    local_cell_indexes_by_subcell[std::pair<size_t,size_t>(subcell_dimension, subcell_index)].push_back(element_local_index);
  }

  // We now have a list of all local cells 'touching' the side set, as well as the faces, edges, and vertices that define the 'touch'

  // For our purposes we only really want the cells and faces associated with the side

  // Each cell can have multiple faces on the sideset
  std::unordered_set<LO> owned_parent_cells_set;
  std::map<LO,std::vector<LO> > owned_parent_cell_index_map;
  {
    // We use a set to avoid duplicates
    for(const auto & cell_subcell_pair : local_cell_indexes_by_subcell){
      const std::pair<size_t,size_t> & subcell_definition = cell_subcell_pair.first;
      const LO subcell_index = LO(subcell_definition.second);
      if(subcell_definition.first == face_subcell_dimension){
        const std::vector<size_t> & local_cell_indexes_for_subcell = cell_subcell_pair.second;
        for(const size_t & local_cell_index : local_cell_indexes_for_subcell){

          // TODO: Super slow - fix this!!
          for(LO i=0;i<LO(mesh_info.local_cells.extent(0)); ++i){
            if(mesh_info.local_cells(i) == LO(local_cell_index)){
              owned_parent_cell_index_map[i].push_back(subcell_index);
            }
          }
        }
      }
    }

    for(const auto & pair : owned_parent_cell_index_map){
      owned_parent_cells_set.insert(pair.first);
    }

  }

  const LO num_owned_cells = owned_parent_cell_index_map.size();

  sideset_info.element_block_name = element_block_name;
  sideset_info.sideset_name = sideset_name;
  sideset_info.cell_topology = mesh.getCellTopology(element_block_name);

  sideset_info.num_owned_cells = num_owned_cells;

  struct face_t{
    face_t(LO c0, LO c1, LO sc0, LO sc1)
    {
      cell_0=c0;
      cell_1=c1;
      subcell_index_0=sc0;
      subcell_index_1=sc1;
    }
    LO cell_0;
    LO cell_1;
    LO subcell_index_0;
    LO subcell_index_1;
  };


  // Figure out how many cells on the other side of the sideset are ghost or virtual
  std::unordered_set<LO> ghstd_parent_cells_set, virtual_parent_cells_set;
  std::vector<face_t> faces;
  {
    LO parent_virtual_cell_offset = mesh_info.num_owned_cells + mesh_info.num_ghstd_cells;
    for(const auto & local_cell_index_pair : owned_parent_cell_index_map){
      const LO local_cell = local_cell_index_pair.first;
      const std::vector<LO> & subcell_indexes_vec = local_cell_index_pair.second;

      for(const LO & subcell_index : subcell_indexes_vec){

        const LO face = mesh_info.cell_to_faces(local_cell, subcell_index);
        const LO face_other_side = (mesh_info.face_to_cells(face,0) == local_cell) ? 1 : 0;

        TEUCHOS_ASSERT(subcell_index == mesh_info.face_to_lidx(face, 1-face_other_side));

        const LO other_side_cell = mesh_info.face_to_cells(face, face_other_side);
        const LO other_side_subcell_index = mesh_info.face_to_lidx(face, face_other_side);

        faces.push_back(face_t(local_cell, other_side_cell, subcell_index, other_side_subcell_index));

        if(other_side_cell >= parent_virtual_cell_offset){
          virtual_parent_cells_set.insert(other_side_cell);
        } else {
          ghstd_parent_cells_set.insert(other_side_cell);
        }
      }
    }
  }

  std::vector<LO> all_cells;
  all_cells.insert(all_cells.end(),owned_parent_cells_set.begin(),owned_parent_cells_set.end());
  all_cells.insert(all_cells.end(),ghstd_parent_cells_set.begin(),ghstd_parent_cells_set.end());
  all_cells.insert(all_cells.end(),virtual_parent_cells_set.begin(),virtual_parent_cells_set.end());

  sideset_info.num_ghstd_cells = ghstd_parent_cells_set.size();
  sideset_info.num_virtual_cells = virtual_parent_cells_set.size();

  const LO num_real_cells = sideset_info.num_owned_cells + sideset_info.num_ghstd_cells;
  const LO num_total_cells = num_real_cells + sideset_info.num_virtual_cells;
  const LO num_vertices_per_cell = mesh_info.cell_vertices.extent(1);
  const LO num_dims = mesh_info.cell_vertices.extent(2);

  sideset_info.global_cells = Kokkos::View<GO*>("global_cells", num_total_cells);
  sideset_info.local_cells = Kokkos::View<LO*>("local_cells", num_total_cells);
  sideset_info.cell_vertices = Kokkos::View<double***,PHX::Device>("cell_vertices", num_total_cells, num_vertices_per_cell, num_dims);
  Kokkos::deep_copy(sideset_info.cell_vertices,0.);

  for(LO i=0; i<num_total_cells; ++i){
    const LO parent_cell = all_cells[i];
    sideset_info.local_cells(i) = mesh_info.local_cells(parent_cell);
    sideset_info.global_cells(i) = mesh_info.global_cells(parent_cell);
    for(LO j=0; j<num_vertices_per_cell; ++j){
      for(LO k=0; k<num_dims; ++k){
        sideset_info.cell_vertices(i,j,k) = mesh_info.cell_vertices(parent_cell,j,k);
      }
    }
  }

  // Now we have to set the connectivity for the faces.

  const LO num_faces = faces.size();
  const LO num_faces_per_cell = mesh_info.cell_to_faces.extent(1);

  sideset_info.face_to_cells = Kokkos::View<LO*[2]>("face_to_cells", num_faces);
  sideset_info.face_to_lidx = Kokkos::View<LO*[2]>("face_to_lidx", num_faces);
  sideset_info.cell_to_faces = Kokkos::View<LO**>("cell_to_faces", num_total_cells, num_faces_per_cell);

  // Default the system with invalid cell index - this will be most of the entries
  Kokkos::deep_copy(sideset_info.cell_to_faces, -1);

  for(int face_index=0;face_index<num_faces;++face_index){
    const face_t & face = faces[face_index];
    const LO & cell_0 = face.cell_0;
    const LO & cell_1 = face.cell_1;

    // TODO: Super expensive... need to find a better way
    const LO sideset_cell_0 = std::distance(all_cells.begin(), std::find(all_cells.begin(), all_cells.end(), cell_0));
    const LO sideset_cell_1 = std::distance(all_cells.begin(), std::find(all_cells.begin()+num_owned_cells, all_cells.end(), cell_1));

    sideset_info.face_to_cells(face_index,0) = sideset_cell_0;
    sideset_info.face_to_cells(face_index,1) = sideset_cell_1;

    sideset_info.face_to_lidx(face_index,0) = face.subcell_index_0;
    sideset_info.face_to_lidx(face_index,1) = face.subcell_index_1;

    sideset_info.cell_to_faces(sideset_cell_0,face.subcell_index_0) = face_index;
    sideset_info.cell_to_faces(sideset_cell_1,face.subcell_index_1) = face_index;

  }

}

}

template <typename LO, typename GO>
void
generateLocalMeshInfo(const panzer_stk::STK_Interface & mesh,
                      panzer::LocalMeshInfo<LO,GO> & mesh_info)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  //typedef Tpetra::CrsMatrix<int,LO,GO> crs_type;
  typedef Tpetra::Map<LO,GO,panzer::TpetraNodeType> map_type;
  typedef Tpetra::Import<LO,GO,panzer::TpetraNodeType> import_type;
  //typedef Tpetra::MultiVector<double,LO,GO> mvec_type;
  //typedef Tpetra::MultiVector<GO,LO,GO> ordmvec_type;

  // Make sure the STK interface is valid
  TEUCHOS_ASSERT(mesh.isInitialized());

  // This is required by some of the STK stuff
  TEUCHOS_ASSERT(typeid(LO) == typeid(int));

  Teuchos::RCP<const Teuchos::Comm<int> > comm = mesh.getComm();

  TEUCHOS_FUNC_TIME_MONITOR_DIFF("panzer_stk::generateLocalMeshInfo",GenerateLocalMeshInfo);

  // This horrible line of code is required since the connection manager only takes rcps of a mesh
  RCP<const panzer_stk::STK_Interface> mesh_rcp = Teuchos::rcpFromRef(mesh);
  // We're allowed to do this since the connection manager only exists in this scope... even though it is also an RCP...

  // extract topology handle
  RCP<panzer::ConnManager<LO,GO> > conn_rcp = rcp(new panzer_stk::STKConnManager<GO>(mesh_rcp));
  panzer::ConnManager<LO,GO> & conn = *conn_rcp;

  // build cell to node map
  Kokkos::View<GO**> owned_cell_to_nodes;
  buildCellToNodes(conn, owned_cell_to_nodes);

  // build the local to global cell ID map
  ///////////////////////////////////////////////////////////
  Kokkos::View<GO*> owned_cells;
  buildCellGlobalIDs(conn, owned_cells);

  // get neighboring cells
  ///////////////////////////////////////////////////////////
  Kokkos::View<const GO*> ghstd_cells = buildGhostedCellOneRing<GO>(comm,owned_cells,owned_cell_to_nodes);

  // build cell maps
  /////////////////////////////////////////////////////////////////////

  RCP<map_type> owned_cell_map = rcp(new map_type(-1,owned_cells,  0,comm));
  RCP<map_type> ghstd_cell_map = rcp(new map_type(-1,ghstd_cells,0,comm));

  // build importer: cell importer, owned to ghstd
  RCP<import_type> cellimport_own2ghst = rcp(new import_type(owned_cell_map,ghstd_cell_map));

  // read all the vertices associated with these elements, get ghstd contributions
  /////////////////////////////////////////////////////////////////////
  std::vector<std::size_t> localCells(owned_cells.extent(0),-1);
  for(size_t i=0;i<localCells.size();i++){
    localCells[i] = i;
  }

  // TODO: This all needs to be rewritten for when element blocks have different cell topologies
  std::vector<std::string> element_block_names;
  mesh.getElementBlockNames(element_block_names);

  const shards::CellTopology & cell_topology = *(mesh.getCellTopology(element_block_names[0]));

  const int space_dim = cell_topology.getDimension();
  const int vertices_per_cell = cell_topology.getVertexCount();
  const int faces_per_cell = cell_topology.getSubcellCount(space_dim-1);

  // Kokkos::View<double***> owned_vertices("owned_vertices",localCells.size(),vertices_per_cell,space_dim);
  Kokkos::DynRankView<double,PHX::Device> owned_vertices("owned_vertices",localCells.size(),vertices_per_cell,space_dim);
  mesh.getElementVerticesNoResize(localCells,owned_vertices);

  // this builds a ghstd vertex array
  Kokkos::DynRankView<double,PHX::Device> ghstd_vertices = buildGhostedVertices(*cellimport_own2ghst,owned_vertices);

  // build edge to cell neighbor mapping
  //////////////////////////////////////////////////////////////////

  std::unordered_map<GO,int> global_to_local;
  global_to_local[-1] = -1; // this is the "no neighbor" flag
  for(size_t i=0;i<owned_cells.extent(0);i++)
    global_to_local[owned_cells(i)] = i;
  for(size_t i=0;i<ghstd_cells.extent(0);i++)
    global_to_local[ghstd_cells(i)] = i+Teuchos::as<int>(owned_cells.extent(0));

  // this class comes from Mini-PIC and Matt B
  RCP<panzer::FaceToElement<LO,GO> > faceToElement = rcp(new panzer::FaceToElement<LO,GO>());
  faceToElement->initialize(conn);
  auto elems_by_face = faceToElement->getFaceToElementsMap();
  auto face_to_lidx  = faceToElement->getFaceToCellLocalIdxMap();

  const int num_owned_cells =owned_cells.extent(0);
  const int num_ghstd_cells =ghstd_cells.extent(0);

//  print_view("elems_by_face",elems_by_face);

  // We also need to consider faces that connect to cells that do not exist, but are needed for boundary conditions
  // We dub them virtual cell since there should be no geometry associated with them, or topology really
  // They exist only for datastorage so that they are consistent with 'real' cells from an algorithm perspective

  // Each virtual face (face linked to a '-1' cell) requires a virtual cell (i.e. turn the '-1' into a virtual cell)
  // Virtual cells are those that do not exist but are connected to an owned cell
  // Note - in the future, ghosted cells will also need to connect to virtual cells at boundary conditions, but for the moment we will ignore this.

  // Iterate over all faces and identify the faces connected to a potential virtual cell
  std::vector<int> all_boundary_faces;
  const int num_faces = elems_by_face.extent(0);
  for(int face=0;face<num_faces;++face){
    if(elems_by_face(face,0) < 0 or elems_by_face(face,1) < 0){
      all_boundary_faces.push_back(face);
    }
  }
  const LO num_virtual_cells = all_boundary_faces.size();

  // total cells and faces include owned, ghosted, and virtual
  const LO num_real_cells = num_owned_cells + num_ghstd_cells;
  const LO num_total_cells = num_real_cells + num_virtual_cells;
  const LO num_total_faces = elems_by_face.extent(0);

  // Create some global indexes associated with the virtual cells
  // Note: We are assuming that virtual cells belong to ranks and are not 'shared' - this will change later on
  Kokkos::View<GO*> virtual_cells = Kokkos::View<GO*>("virtual_cells",num_virtual_cells);
  {
    PANZER_FUNC_TIME_MONITOR_DIFF("Initial global index creation",InitialGlobalIndexCreation);

    const int num_ranks = comm->getSize();
    const int rank = comm->getRank();

    std::vector<GO> owned_cell_distribution(num_ranks,0);
    {
      std::vector<GO> my_owned_cell_distribution(num_ranks,0);
      my_owned_cell_distribution[rank] = num_owned_cells;

      Teuchos::reduceAll(*comm,Teuchos::REDUCE_SUM, num_ranks, my_owned_cell_distribution.data(),owned_cell_distribution.data());
    }

    std::vector<GO> virtual_cell_distribution(num_ranks,0);
    {
      std::vector<GO> my_virtual_cell_distribution(num_ranks,0);
      my_virtual_cell_distribution[rank] = num_virtual_cells;

      Teuchos::reduceAll(*comm,Teuchos::REDUCE_SUM, num_ranks, my_virtual_cell_distribution.data(),virtual_cell_distribution.data());
    }

    GO num_global_real_cells=0;
    for(int i=0;i<num_ranks;++i){
      num_global_real_cells+=owned_cell_distribution[i];
    }

    GO global_virtual_start_idx = num_global_real_cells;
    for(int i=0;i<rank;++i){
      global_virtual_start_idx += virtual_cell_distribution[i];
    }

    for(int i=0;i<num_virtual_cells;++i){
      virtual_cells(i) = global_virtual_start_idx + GO(i);
    }

  }

  // Lookup cells connected to a face
  Kokkos::View<LO*[2]> face_to_cells = Kokkos::View<LO*[2]>("face_to_cells",num_total_faces);

  // Lookup local face indexes given cell and left/right state (0/1)
  Kokkos::View<LO*[2]> face_to_localidx = Kokkos::View<LO*[2]>("face_to_localidx",num_total_faces);

  // Lookup face index given a cell and local face index
  Kokkos::View<LO**> cell_to_face = Kokkos::View<LO**>("cell_to_face",num_total_cells,faces_per_cell);

  // initialize with negative one cells that are not associated with a face
  Kokkos::deep_copy(cell_to_face,-1);

  // Transfer information from 'faceToElement' datasets to local arrays
  {
    PANZER_FUNC_TIME_MONITOR_DIFF("Transer faceToElement to local",TransferFaceToElementLocal);

    int virtual_cell_index = num_real_cells;
    for(size_t f=0;f<elems_by_face.extent(0);f++) {

      const GO global_c0 = elems_by_face(f,0);
      const GO global_c1 = elems_by_face(f,1);

      // make sure that no bonus cells get in here
      TEUCHOS_ASSERT(global_to_local.find(global_c0)!=global_to_local.end());
      TEUCHOS_ASSERT(global_to_local.find(global_c1)!=global_to_local.end());

      auto c0 = global_to_local[global_c0];
      auto lidx0 = face_to_lidx(f,0);
      auto c1 = global_to_local[global_c1];
      auto lidx1 = face_to_lidx(f,1);

      // Test for virtual cells

      // Left cell
      if(c0 < 0){
        // Virtual cell - create it!
        c0 = virtual_cell_index++;

        // We need the subcell_index to line up between real and virtual cell
        // This way the face has the same geometry... though the face normal
        // will point in the wrong direction
        lidx0 = lidx1;
      }
      cell_to_face(c0,lidx0) = f;


      // Right cell
      if(c1 < 0){
        // Virtual cell - create it!
        c1 = virtual_cell_index++;

        // We need the subcell_index to line up between real and virtual cell
        // This way the face has the same geometry... though the face normal
        // will point in the wrong direction
        lidx1 = lidx0;
      }
      cell_to_face(c1,lidx1) = f;

      // Faces point from low cell index to high cell index
      if(c0<c1){
        face_to_cells(f,0) = c0;
        face_to_localidx(f,0) = lidx0;
        face_to_cells(f,1) = c1;
        face_to_localidx(f,1) = lidx1;
      } else {
        face_to_cells(f,0) = c1;
        face_to_localidx(f,0) = lidx1;
        face_to_cells(f,1) = c0;
        face_to_localidx(f,1) = lidx0;
      }

      // We should avoid having two virtual cells linked together.
      TEUCHOS_ASSERT(c0<num_real_cells or c1<num_real_cells)

    }
  }

  // at this point all the data structures have been built, so now we can "do" DG.
  // There are two of everything, an "owned" data structured corresponding to "owned"
  // cells. And a "ghstd" data structure corresponding to ghosted cells
  ////////////////////////////////////////////////////////////////////////////////////
  {
    PANZER_FUNC_TIME_MONITOR_DIFF("Assign Indices",AssignIndices);
    mesh_info.cell_to_faces           = cell_to_face;
    mesh_info.face_to_cells           = face_to_cells;      // faces
    mesh_info.face_to_lidx            = face_to_localidx;

    mesh_info.num_owned_cells = owned_cells.extent(0);
    mesh_info.num_ghstd_cells = ghstd_cells.extent(0);
    mesh_info.num_virtual_cells = virtual_cells.extent(0);

    mesh_info.global_cells = Kokkos::View<GO*>("global_cell_indices",num_total_cells);
    mesh_info.local_cells = Kokkos::View<LO*>("local_cell_indices",num_total_cells);

    for(int i=0;i<num_owned_cells;++i){
      mesh_info.global_cells(i) = owned_cells(i);
      mesh_info.local_cells(i) = i;
    }

    for(int i=0;i<num_ghstd_cells;++i){
      mesh_info.global_cells(i+num_owned_cells) = ghstd_cells(i);
      mesh_info.local_cells(i+num_owned_cells) = i+num_owned_cells;
    }

    for(int i=0;i<num_virtual_cells;++i){
      mesh_info.global_cells(i+num_real_cells) = virtual_cells(i);
      mesh_info.local_cells(i+num_real_cells) = i+num_real_cells;
    }

    mesh_info.cell_vertices = Kokkos::View<double***, PHX::Device>("cell_vertices",num_total_cells,vertices_per_cell,space_dim);

    // Initialize coordinates to zero
    Kokkos::deep_copy(mesh_info.cell_vertices, 0.);

    for(int i=0;i<num_owned_cells;++i){
      for(int j=0;j<vertices_per_cell;++j){
        for(int k=0;k<space_dim;++k){
          mesh_info.cell_vertices(i,j,k) = owned_vertices(i,j,k);
        }
      }
    }

    for(int i=0;i<num_ghstd_cells;++i){
      for(int j=0;j<vertices_per_cell;++j){
        for(int k=0;k<space_dim;++k){
          mesh_info.cell_vertices(i+num_owned_cells,j,k) = ghstd_vertices(i,j,k);
        }
      }
    }

    // This will backfire at some point, but we're going to make the virtual cell have the same geometry as the cell it interfaces with
    // This way we can define a virtual cell geometry without extruding the face outside of the domain
    {
      PANZER_FUNC_TIME_MONITOR_DIFF("Assign geometry traits",AssignGeometryTraits);
      for(int i=0;i<num_virtual_cells;++i){

        const LO virtual_cell = i+num_real_cells;
        bool exists = false;
        for(int local_face=0; local_face<faces_per_cell; ++local_face){
          const LO face = cell_to_face(virtual_cell, local_face);
          if(face >= 0){
            exists = true;
            const LO other_side = (face_to_cells(face, 0) == virtual_cell) ? 1 : 0;
            const LO real_cell = face_to_cells(face,other_side);
            TEUCHOS_ASSERT(real_cell < num_real_cells);
            for(int j=0;j<vertices_per_cell;++j){
              for(int k=0;k<space_dim;++k){
                mesh_info.cell_vertices(virtual_cell,j,k) = mesh_info.cell_vertices(real_cell,j,k);
              }
            }
            break;
          }
        }
        TEUCHOS_TEST_FOR_EXCEPT_MSG(!exists, "panzer_stk::generateLocalMeshInfo : Virtual cell is not linked to real cell");
      }
    }
  }

  // Setup element blocks and sidesets
  std::vector<std::string> sideset_names;
  mesh.getSidesetNames(sideset_names);

  for(const std::string & element_block_name : element_block_names){
    PANZER_FUNC_TIME_MONITOR_DIFF("Set up setupLocalMeshBlockInfo",SetupLocalMeshBlockInfo);
    panzer::LocalMeshBlockInfo<LO,GO> & block_info = mesh_info.element_blocks[element_block_name];
    setupLocalMeshBlockInfo(mesh, conn, mesh_info, element_block_name, block_info);

    // Setup sidesets
    for(const std::string & sideset_name : sideset_names){
      PANZER_FUNC_TIME_MONITOR_DIFF("Setup LocalMeshSidesetInfo",SetupLocalMeshSidesetInfo);
      panzer::LocalMeshSidesetInfo<LO,GO> & sideset_info = mesh_info.sidesets[element_block_name][sideset_name];
      setupLocalMeshSidesetInfo(mesh, conn, mesh_info, element_block_name, sideset_name, sideset_info);
    }

  }

}

}

// Explicit instantiation
template
void
panzer_stk::generateLocalMeshInfo<int,int>(const panzer_stk::STK_Interface & mesh,
                                           panzer::LocalMeshInfo<int,int> & mesh_info);

#ifndef PANZER_ORDINAL64_IS_INT
template
void
panzer_stk::generateLocalMeshInfo<int,panzer::Ordinal64>(const panzer_stk::STK_Interface & mesh,
                                           panzer::LocalMeshInfo<int,panzer::Ordinal64> & mesh_info);
#endif
