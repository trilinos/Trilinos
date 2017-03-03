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

#include "Panzer_STK_LocalMeshUtilities.hpp"

#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_FaceToElement.hpp"
#include "Panzer_ConnManager.hpp"
#include "Panzer_STKConnManager.hpp"

#include "Panzer_FieldPattern.hpp"
#include "Panzer_NodalFieldPattern.hpp"
#include "Panzer_EdgeFieldPattern.hpp"
#include "Panzer_FaceFieldPattern.hpp"
#include "Panzer_ElemFieldPattern.hpp"

#include "Teuchos_Assert.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"

#include <vector>

namespace panzer_stk
{

namespace tools
{

/** Build a Kokkos array of all the global cell IDs from a connection manager.
  * Note that this is mapping between local IDs to global IDs.
  */
template <typename LO,typename GO>
void
buildCellGlobalIDs(panzer::ConnManager<LO,GO> & conn, Kokkos::View<GO*> & globals)
{
  // extract topologies, and build global connectivity...currently assuming only one topology
  std::vector<shards::CellTopology> elementBlockTopologies;
  conn.getElementBlockTopologies(elementBlockTopologies);

  const shards::CellTopology & topology = elementBlockTopologies[0];
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
  panzer::NodalFieldPattern pattern(elementBlockTopologies[0]);
  conn.buildConnectivity(pattern);

  // calculate total number of local cells
  std::vector<std::string> block_ids;
  conn.getElementBlockIds(block_ids);

  // compute total cells and maximum nodes
  std::size_t totalCells=0, maxNodes=0;
  for (std::size_t which_blk=0;which_blk<block_ids.size();which_blk++) {
    // get the elem to face mapping
    const std::vector<LO> & localIDs = conn.getElementBlock(block_ids[which_blk]);
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
Teuchos::RCP<const Tpetra::Map<LO,GO> >
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

  typedef Tpetra::Map<LO,GO> map_type;

  // get locally unique global ids
  std::set<GO> global_nodes;
  for(unsigned int i=0;i<cells_to_nodes.dimension(0);i++)
    for(unsigned int j=0;j<cells_to_nodes.dimension(1);j++)
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
Teuchos::RCP<Tpetra::CrsMatrix<LO,LO,GO> >
buildNodeToCellMatrix(const Teuchos::RCP<const Teuchos::Comm<int> > & comm,
                      Kokkos::View<const GO*> owned_cells,
                      Kokkos::View<const GO**> owned_cells_to_nodes)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef Tpetra::Map<LO,GO> map_type;
  typedef Tpetra::CrsMatrix<LO,LO,GO> crs_type;
  typedef Tpetra::Import<LO,GO> import_type;

  TEUCHOS_ASSERT(owned_cells.dimension(0)==owned_cells_to_nodes.dimension(0));

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

    // The matrix is indexed by (global cell, global node) = local node
    cell_to_node = rcp(new crs_type(cell_map,0));
    cell_to_node->resumeFill();

    // fill in the cell to node matrix
    const unsigned int num_local_cells = owned_cells_to_nodes.dimension(0);
    const unsigned int num_nodes_per_cell = owned_cells_to_nodes.dimension(1);
    std::vector<LO> local_node_indexes(num_nodes_per_cell);
    std::vector<GO> global_node_indexes(num_nodes_per_cell);
    for(unsigned int i=0;i<num_local_cells;i++) {
      const GO global_cell_index = owned_cells(i);
  //    std::vector<LO> vals(cells_to_nodes.dimension(1));
  //    std::vector<GO> cols(cells_to_nodes.dimension(1));
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
    // Create an object designed to transpose the (global cell, global node) matrix to give
    // a (global node, global cell) matrix
    Tpetra::RowMatrixTransposer<LO,LO,GO> transposer(cell_to_node);

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
  typedef Tpetra::CrsMatrix<int,int,GO> crs_type;

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
  for(size_t i=0;i<cells.dimension(0);i++) {
    unique_cells.insert(cells(i));
  }

  // The set of ghost cells that share a global node with an owned cell
  std::set<GO> ghstd_cells_set;

  // Get a list of global node indexes associated with the cells owned by this process
//  auto node_map = node_to_cell->getRangeMap()->getMyGlobalIndices();
  auto node_map = node_to_cell->getMap()->getMyGlobalIndices();

  // Iterate through the global node indexes associated with this process
  for(size_t i=0;i<node_map.dimension(0);i++) {
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
buildGhostedVertices(const Tpetra::Import<int,GO> & importer,
                     Kokkos::DynRankView<const double,PHX::Device> owned_vertices)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef Tpetra::MultiVector<double,int,GO> mvec_type;
  typedef typename mvec_type::dual_view_type dual_view_type;

  size_t owned_cell_cnt = importer.getSourceMap()->getNodeNumElements();
  size_t ghstd_cell_cnt = importer.getTargetMap()->getNodeNumElements();
  int vertices_per_cell = owned_vertices.dimension(1);
  int space_dim         = owned_vertices.dimension(2);

  TEUCHOS_ASSERT(owned_vertices.dimension(0)==owned_cell_cnt);

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

}


template <typename LO, typename GO>
panzer::LocalMeshInfo<LO,GO>
generateLocalMeshInfo(const panzer_stk::STK_Interface & mesh,
                      const std::string & element_block_name)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef Tpetra::CrsMatrix<int,LO,GO> crs_type;
  typedef Tpetra::Map<LO,GO> map_type;
  typedef Tpetra::Import<LO,GO> import_type;
  typedef Tpetra::MultiVector<double,LO,GO> mvec_type;
  typedef Tpetra::MultiVector<GO,LO,GO> ordmvec_type;

  // This is required by some of the STK stuff
  TEUCHOS_ASSERT(typeid(LO) == typeid(int));

  Teuchos::RCP<const Teuchos::Comm<int> > comm = mesh.getComm();

  const shards::CellTopology & cell_topology = *(mesh.getCellTopology(element_block_name));

  const int space_dim = cell_topology.getDimension();
  const int vertices_per_cell = cell_topology.getVertexCount();
  const int faces_per_cell = cell_topology.getSubcellCount(space_dim-1);

  TEUCHOS_FUNC_TIME_MONITOR("panzer_stk::generateLocalMeshInfo");

  // This horrible line of code is required since the connection manager only takes rcps of a mesh
  RCP<const panzer_stk::STK_Interface> mesh_rcp = Teuchos::rcpFromRef(mesh);
  // We're allowed to do this since the connection manager only exists in this scope... even though it is also an RCP...

  // extract topology handle
  RCP<panzer::ConnManager<LO,GO> > conn_rcp = rcp(new panzer_stk::STKConnManager<GO>(mesh_rcp));
  panzer::ConnManager<LO,GO> & conn = *conn_rcp;

  // build cell to node map
  Kokkos::View<GO**> owned_cell_to_nodes;
  tools::buildCellToNodes(conn, owned_cell_to_nodes);

  // build the local to global cell ID map
  ///////////////////////////////////////////////////////////
  Kokkos::View<GO*> owned_cells;
  tools::buildCellGlobalIDs(conn, owned_cells);

  // get neighboring cells
  ///////////////////////////////////////////////////////////
  Kokkos::View<const GO*> ghstd_cells = tools::buildGhostedCellOneRing<GO>(comm,owned_cells,owned_cell_to_nodes);

  // build cell maps
  /////////////////////////////////////////////////////////////////////

  RCP<map_type> owned_cell_map   = rcp(new map_type(-1,owned_cells,  0,comm));
  RCP<map_type> ghstd_cell_map = rcp(new map_type(-1,ghstd_cells,0,comm));

  // build importer: cell importer, owned to ghstd
  RCP<import_type> cellimport_own2ghst = rcp(new import_type(owned_cell_map,ghstd_cell_map));

  // read all the vertices associated with these elements, get ghstd contributions
  /////////////////////////////////////////////////////////////////////
  std::vector<std::size_t> localCells(owned_cells.dimension(0),-1);
  for(size_t i=0;i<localCells.size();i++){
    localCells[i] = i;
  }

  // Kokkos::View<double***> owned_vertices("owned_vertices",localCells.size(),vertices_per_cell,space_dim);
  Kokkos::DynRankView<double,PHX::Device> owned_vertices("owned_vertices",localCells.size(),vertices_per_cell,space_dim);
  mesh.getElementVerticesNoResize(localCells,owned_vertices);

  // this builds a ghstd vertex array
  Kokkos::DynRankView<double,PHX::Device> ghstd_vertices = tools::buildGhostedVertices(*cellimport_own2ghst,owned_vertices);

  // build edge to cell neighbor mapping
  //////////////////////////////////////////////////////////////////

  std::unordered_map<GO,int> global_to_local;
  global_to_local[-1] = -1; // this is the "no neighbor" flag
  for(size_t i=0;i<owned_cells.dimension(0);i++)
    global_to_local[owned_cells(i)] = i;
  for(size_t i=0;i<ghstd_cells.dimension(0);i++)
    global_to_local[ghstd_cells(i)] = i+Teuchos::as<int>(owned_cells.dimension(0));

  // this class comes from Mini-PIC and Matt B
  RCP<panzer::FaceToElement<LO,GO> > faceToElement = rcp(new panzer::FaceToElement<LO,GO>());
  faceToElement->initialize(conn);
  auto elems_by_face = faceToElement->getFaceToElementsMap();
  auto face_to_lidx  = faceToElement->getFaceToCellLocalIdxMap();

  const int num_owned_cells =owned_cells.dimension(0);
  const int num_ghstd_cells =ghstd_cells.dimension(0);

//  print_view("elems_by_face",elems_by_face);

  // We also need to consider faces that connect to cells that do not exist, but are needed for boundary conditions
  // We dub them virtual cell since there should be no geometry associated with them, or topology really
  // They exist only for datastorage so that they are consistent with 'real' cells from an algorithm perspective

  // Each virtual face (face linked to a '-1' cell) requires a virtual cell (i.e. turn the '-1' into a virtual cell)
  // Virtual cells are those that do not exist but are connected to an owned cell
  // Note - in the future, ghosted cells will also need to connect to virtual cells at boundary conditions, but for the moment we will ignore this.

  // Iterate over all faces and identify the faces connected to a potential virtual cell
  std::vector<int> all_boundary_faces;
  const int num_faces = elems_by_face.dimension(0);
  for(int face=0;face<num_faces;++face){
    if(elems_by_face(face,0) < 0 or elems_by_face(face,1) < 0){
      all_boundary_faces.push_back(face);
    }
  }
  const LO num_virtual_cells = all_boundary_faces.size();

  // total cells and faces include owned, ghosted, and virtual
  const LO num_real_cells = num_owned_cells + num_ghstd_cells;
  const LO num_total_cells = num_real_cells + num_virtual_cells;
  const LO num_total_faces = elems_by_face.dimension(0);

  // Local indexes associated with virtual cells
  Kokkos::View<LO*> virtual_cells = Kokkos::View<LO*>("virtual_cells",num_virtual_cells);

  // Fill virtual cells - for now it is just an offset
  for(LO i=0; i<num_virtual_cells;++i){
    virtual_cells(i) = i + num_real_cells;
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
    int virtual_cell_index = num_real_cells;
    for(size_t f=0;f<elems_by_face.dimension(0);f++) {

//      printf("face %i has cells %i and %i\n",f,elems_by_face(f,0),elems_by_face(f,1));

      const GO global_c0 = elems_by_face(f,0);
      const GO global_c1 = elems_by_face(f,1);

      // make sure that no bonus cells get in here
      TEUCHOS_ASSERT(global_to_local.find(global_c0)!=global_to_local.end());
      TEUCHOS_ASSERT(global_to_local.find(global_c1)!=global_to_local.end());

      // Left cell
      auto c0 = global_to_local[elems_by_face(f,0)];
      auto lidx0 = face_to_lidx(f,0);
      if(c0 < 0){
        // Virtual cell - create it!
        c0 = virtual_cell_index++;

        // Virtual cells always connect to the real domain through face 0
        lidx0 = 0;
      }
      cell_to_face(c0,lidx0) = f;


      // Right cell
      auto c1 = global_to_local[elems_by_face(f,1)];
      auto lidx1 = face_to_lidx(f,1);
      if(c1 < 0){
        // Virtual cell - create it!
        c1 = virtual_cell_index++;

        // Virtual cells always connect to the real domain through face 0
        lidx1 = 0;
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

  panzer::LocalMeshInfo<LO,GO> localMeshInfo;

  localMeshInfo.element_block_name      = element_block_name;
  localMeshInfo.sideset_name            = "";
  localMeshInfo.cell_to_face            = cell_to_face;
  localMeshInfo.face_to_cells           = face_to_cells;      // faces
  localMeshInfo.face_to_lidx            = face_to_localidx;
  localMeshInfo.owned_cells             = owned_cells;        // cells
  localMeshInfo.ghstd_cells             = ghstd_cells;
  localMeshInfo.virtual_cells           = virtual_cells;
  localMeshInfo.owned_vertices          = owned_vertices;     // vertices
  localMeshInfo.ghstd_vertices          = ghstd_vertices;
  localMeshInfo.cell_topology           = mesh.getCellTopology(element_block_name);

  return localMeshInfo;
}





template <typename LO, typename GO>
panzer::LocalMeshInfo<LO,GO>
generateLocalSidesetInfo(const panzer_stk::STK_Interface & mesh,
                      const std::string & element_block_name,
                      const std::string & sideset_name)
{

  // This isn't setup yet
  TEUCHOS_ASSERT(false);


  // Get element block, local chunk of mesh, and bulk data
//  const stk::mesh::Part & element_block = *(mesh.getElementBlockPart(description.getElementBlock()));
//  const stk::mesh::Part & local_mesh = *(mesh.getOwnedPart());
//  const stk::mesh::BulkData & bulk_data = *(mesh.getBulkData());


//    const GO face_subcell_dimension = mesh.getCellTopology(description.getElementBlock())->getDimension()-1;
//
//    const int requested_chunk_size = description.getWorksetSize();
//
//
//    Kokkos::DynRankView<double,PHX::Device> data_allocation;
//


//
//    // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//    const bool is_domain_boundary = true;
//
//    //===================================================
//    // Sideset
//
//    // This is a sideset, therefore we must... do something...
//    const std::string & sideset_name = description.getSideset();
//
//    std::vector<stk::mesh::Entity> side_entities;
//
//    // Grab the side entities associated with this sideset on the mesh
//    // Note: Throws exception if element block or sideset doesn't exist
//    try{
//
//      mesh.getMySides(sideset_name,description.getElementBlock(),side_entities);
//
//    } catch(STK_Interface::SidesetException & e) {
//       std::stringstream ss;
//       std::vector<std::string> sideset_names;
//       mesh.getSidesetNames(sideset_names);
//
//       // build an error message
//       ss << e.what() << "\nChoose one of:\n";
//       for(const auto & optional_sideset_name : sideset_names){
//          ss << "\"" << optional_sideset_name << "\"\n";
//       }
//
//       TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(true,std::logic_error,ss.str());
//
//    } catch(STK_Interface::ElementBlockException & e) {
//       std::stringstream ss;
//       std::vector<std::string> element_block_names;
//       mesh.getElementBlockNames(element_block_names);
//
//       // build an error message
//       ss << e.what() << "\nChoose one of:\n";
//       for(const auto & optional_element_block_name : element_block_names){
//          ss << "\"" << optional_element_block_name << "\"\n";
//       }
//
//       TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(true,std::logic_error,ss.str());
//
//    } catch(std::logic_error & e) {
//       std::stringstream ss;
//       ss << e.what() << "\nUnrecognized logic error.\n";
//
//       TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(true,std::logic_error,ss.str());
//
//    }
//
//    // We now have a list of sideset entities, lets unwrap them and create some chunks!
//
//    // Lets go all out and grab everything associated with that side
//
//    std::map<std::pair<int,int>, std::vector<GO> > local_cell_indexes_by_subcell;
//
//    // List of all elements - not really needed
//    std::vector<stk::mesh::Entity> elements;
//
//    // List of local subcell indexes indexed by element:
//    // For example: a Tet (element) would have
//    //  - 4 triangular faces (subcell_index 0-3, subcell_dimension=2)
//    //  - 6 edges (subcell_index 0-5, subcell_dimension=1)
//    //  - 4 vertices (subcell_index 0-3, subcell_dimension=0)
//    // Another example: a Line (element) would have
//    //  - 2 vertices (subcell_index 0-1, subcell_dimension=0)
//    std::vector<GO> subcell_indexes;
//    std::vector<GO> subcell_dimensions;
//
//    // Grab everything
//    panzer_stk::workset_utils::getSideElementCascade(mesh, description.getElementBlock(),side_entities, subcell_dimensions,subcell_indexes,elements);
//
//    // build local cell_ids, mapped by local side id
//    for(GO element_index=0; element_index<elements.size(); ++element_index) {
//      const GO subcell_dimension = subcell_dimensions[element_index];
//      const GO subcell_index = subcell_indexes[element_index];
//      const GO element_local_index = mesh.elementLocalId(elements[element_index]);
//
//      // Add subcell to map
//      local_cell_indexes_by_subcell[std::pair<int,int>(subcell_dimension, subcell_index)].push_back(element_local_index);
//    }
//
//    // We now have a list of all local cells 'touching' the side set, as well as the faces, edges, and vertices that define the 'touch'
//
//    // For our purposes we only really want the cells and faces associated with the side
//
//    std::vector<GO> local_cell_indexes;
//    {
//      // We use a set to avoid duplicates
//      std::set<GO> local_cell_indexes_set;
//      for(const auto & cell_subcell_pair : local_cell_indexes_by_subcell){
//        const std::pair<GO,GO> & subcell_definition = cell_subcell_pair.first;
//        if(subcell_definition.first == face_subcell_dimension){
//          const std::vector<GO> & local_cell_indexes_for_subcell = cell_subcell_pair.second;
//          for(const GO & local_cell_index : local_cell_indexes_for_subcell){
//            local_cell_indexes_set.insert(local_cell_index);
//          }
//        }
//      }
//
//      for(const GO & cell_index : local_cell_indexes_set){
//        local_cell_indexes.push_back(cell_index);
//      }
//    }
//
//
//    const GO num_total_cells = local_cell_indexes.size();
//    const GO num_left_over = num_total_cells % requested_chunk_size;
//    const GO num_chunks = num_total_cells / requested_chunk_size + ((num_left_over > 0)?1:0);
//
//    // We create our chunks
//    for(GO chunk_index=0; chunk_index<num_chunks; ++chunk_index){
//      const GO cell_index_start = chunk_index * requested_chunk_size;
//      const GO cell_index_end = std::min(num_total_cells, cell_index_start + requested_chunk_size);
//      const GO num_real_cells = cell_index_end - cell_index_start;
//      const GO num_virtual_cells = (is_domain_boundary)?num_real_cells:0;
//      const GO num_cells = num_real_cells + num_virtual_cells;
//
//      chunks.push_back(panzer::LocalMeshChunk<LO,GO>());
//      panzer::LocalMeshChunk<LO,GO> & chunk = chunks.back();
//
//      chunk.num_cells = num_cells;
//      chunk.subcell_dimension = face_subcell_dimension;
//      chunk.element_block_name = description.getElementBlock();
//      chunk.sideset_name = sideset_name;
//      chunk.local_cell_global_indexes.assign(local_cell_indexes.begin()+cell_index_start,local_cell_indexes.begin()+cell_index_end);
//
//      // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//      // Side has either virtual or ghost cells
//      chunk.virtual_cell_global_indexes;
//
//      // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//      // Side has either virtual or ghost cells
//      chunk.ghost_cell_global_indexes;
//
//      // Must take virtual cells into account
//      {
//        mesh.getElementVertices(chunk.local_cell_global_indexes, description.getElementBlock(), data_allocation);
//
//        const GO num_cells_with_vertices = data_allocation.dimension_0();
//        const GO num_vertices_per_cell = data_allocation.dimension_1();
//        const GO num_dims = data_allocation.dimension_2();
//
//        TEUCHOS_ASSERT(num_real_cells == num_cells_with_vertices);
//
//        chunk.cell_vertices = Kokkos::View<double***>("cell_vertices", num_cells,num_vertices_per_cell,num_dims);
//        Kokkos::deep_copy(chunk.cell_vertices,0.);
//
//        for(GO i=0;i<num_cells_with_vertices;++i){
//          for(GO j=0;j<num_vertices_per_cell;++j){
//            for(GO k=0;k<num_dims;++k){
//              chunk.cell_vertices(i,j,k) = data_allocation(i,j,k);
//            }
//          }
//        }
//      }
//
//      // Each face in this side set connects a real cell to a ghost or virtual cell
//
//      const std::vector<GO> & other_side_cells = (is_domain_boundary) ? chunk.virtual_cell_global_indexes : chunk.ghost_cell_global_indexes;
//
//      chunk.face_to_cells = Kokkos::View<LO*[2],PHX::Device>("face_to_cells", num_real_cells);
//      for(GO i=0; i<num_real_cells; ++i){
//
//      }
//
//
//
//
//
//    }





























  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef Tpetra::CrsMatrix<int,LO,GO> crs_type;
  typedef Tpetra::Map<LO,GO> map_type;
  typedef Tpetra::Import<LO,GO> import_type;
  typedef Tpetra::MultiVector<double,LO,GO> mvec_type;
  typedef Tpetra::MultiVector<GO,LO,GO> ordmvec_type;

  // This is required by some of the STK stuff
  TEUCHOS_ASSERT(typeid(LO) == typeid(int));

  Teuchos::RCP<const Teuchos::Comm<int> > comm = mesh.getComm();

  const shards::CellTopology & cell_topology = *(mesh.getCellTopology(element_block_name));

  const int space_dim = cell_topology.getDimension();
  const int vertices_per_cell = cell_topology.getVertexCount();
  const int faces_per_cell = cell_topology.getSubcellCount(space_dim-1);

  TEUCHOS_FUNC_TIME_MONITOR("panzer_stk::generateLocalMeshInfo");

  // This horrible line of code is required since the connection manager only takes rcps of a mesh
  RCP<const panzer_stk::STK_Interface> mesh_rcp = Teuchos::rcpFromRef(mesh);
  // We're allowed to do this since the connection manager only exists in this scope... even though it is also an RCP...

  // extract topology handle
  RCP<panzer::ConnManager<LO,GO> > conn_rcp = rcp(new panzer_stk::STKConnManager<GO>(mesh_rcp));
  panzer::ConnManager<LO,GO> & conn = *conn_rcp;

  // build cell to node map
  Kokkos::View<GO**> owned_cell_to_nodes;
  tools::buildCellToNodes(conn, owned_cell_to_nodes);

  // build the local to global cell ID map
  ///////////////////////////////////////////////////////////
  Kokkos::View<GO*> owned_cells;
  tools::buildCellGlobalIDs(conn, owned_cells);

  // get neighboring cells
  ///////////////////////////////////////////////////////////
  Kokkos::View<const GO*> ghstd_cells = tools::buildGhostedCellOneRing<GO>(comm,owned_cells,owned_cell_to_nodes);

  // build cell maps
  /////////////////////////////////////////////////////////////////////

  RCP<map_type> owned_cell_map   = rcp(new map_type(-1,owned_cells,  0,comm));
  RCP<map_type> ghstd_cell_map = rcp(new map_type(-1,ghstd_cells,0,comm));

  // build importer: cell importer, owned to ghstd
  RCP<import_type> cellimport_own2ghst = rcp(new import_type(owned_cell_map,ghstd_cell_map));

  // read all the vertices associated with these elements, get ghstd contributions
  /////////////////////////////////////////////////////////////////////
  std::vector<std::size_t> localCells(owned_cells.dimension(0),-1);
  for(size_t i=0;i<localCells.size();i++){
    localCells[i] = i;
  }

  // Kokkos::View<double***> owned_vertices("owned_vertices",localCells.size(),vertices_per_cell,space_dim);
  Kokkos::DynRankView<double,PHX::Device> owned_vertices("owned_vertices",localCells.size(),vertices_per_cell,space_dim);
  mesh.getElementVerticesNoResize(localCells,owned_vertices);

  // this builds a ghstd vertex array
  Kokkos::DynRankView<double,PHX::Device> ghstd_vertices = tools::buildGhostedVertices(*cellimport_own2ghst,owned_vertices);

  // build edge to cell neighbor mapping
  //////////////////////////////////////////////////////////////////

  std::unordered_map<GO,int> global_to_local;
  global_to_local[-1] = -1; // this is the "no neighbor" flag
  for(size_t i=0;i<owned_cells.dimension(0);i++)
    global_to_local[owned_cells(i)] = i;
  for(size_t i=0;i<ghstd_cells.dimension(0);i++)
    global_to_local[ghstd_cells(i)] = i+Teuchos::as<int>(owned_cells.dimension(0));

  // this class comes from Mini-PIC and Matt B
  RCP<panzer::FaceToElement<LO,GO> > faceToElement = rcp(new panzer::FaceToElement<LO,GO>());
  faceToElement->initialize(conn);
  auto elems_by_face = faceToElement->getFaceToElementsMap();
  auto face_to_lidx  = faceToElement->getFaceToCellLocalIdxMap();

  const int num_owned_cells =owned_cells.dimension(0);
  const int num_ghstd_cells =ghstd_cells.dimension(0);

//  print_view("elems_by_face",elems_by_face);

  // We also need to consider faces that connect to cells that do not exist, but are needed for boundary conditions
  // We dub them virtual cell since there should be no geometry associated with them, or topology really
  // They exist only for datastorage so that they are consistent with 'real' cells from an algorithm perspective

  // Each virtual face (face linked to a '-1' cell) requires a virtual cell (i.e. turn the '-1' into a virtual cell)
  // Virtual cells are those that do not exist but are connected to an owned cell
  // Note - in the future, ghosted cells will also need to connect to virtual cells at boundary conditions, but for the moment we will ignore this.

  // Iterate over all faces and identify the faces connected to a potential virtual cell
  std::vector<int> all_boundary_faces;
  const int num_faces = elems_by_face.dimension(0);
  for(int face=0;face<num_faces;++face){
    if(elems_by_face(face,0) < 0 or elems_by_face(face,1) < 0){
      all_boundary_faces.push_back(face);
    }
  }
  const LO num_virtual_cells = all_boundary_faces.size();

  // total cells and faces include owned, ghosted, and virtual
  const LO num_real_cells = num_owned_cells + num_ghstd_cells;
  const LO num_total_cells = num_real_cells + num_virtual_cells;
  const LO num_total_faces = elems_by_face.dimension(0);

  // Local indexes associated with virtual cells
  Kokkos::View<LO*> virtual_cells = Kokkos::View<LO*>("virtual_cells",num_virtual_cells);

  // Fill virtual cells - for now it is just an offset
  for(LO i=0; i<num_virtual_cells;++i){
    virtual_cells(i) = i + num_real_cells;
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
    int virtual_cell_index = num_real_cells;
    for(size_t f=0;f<elems_by_face.dimension(0);f++) {

//      printf("face %i has cells %i and %i\n",f,elems_by_face(f,0),elems_by_face(f,1));

      const GO global_c0 = elems_by_face(f,0);
      const GO global_c1 = elems_by_face(f,1);

      // make sure that no bonus cells get in here
      TEUCHOS_ASSERT(global_to_local.find(global_c0)!=global_to_local.end());
      TEUCHOS_ASSERT(global_to_local.find(global_c1)!=global_to_local.end());

      // Left cell
      auto c0 = global_to_local[elems_by_face(f,0)];
      auto lidx0 = face_to_lidx(f,0);
      if(c0 < 0){
        // Virtual cell - create it!
        c0 = virtual_cell_index++;

        // Virtual cells always connect to the real domain through face 0
        lidx0 = 0;
      }
      cell_to_face(c0,lidx0) = f;


      // Right cell
      auto c1 = global_to_local[elems_by_face(f,1)];
      auto lidx1 = face_to_lidx(f,1);
      if(c1 < 0){
        // Virtual cell - create it!
        c1 = virtual_cell_index++;

        // Virtual cells always connect to the real domain through face 0
        lidx1 = 0;
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

  panzer::LocalMeshInfo<LO,GO> localMeshInfo;

  localMeshInfo.element_block_name      = element_block_name;
  localMeshInfo.sideset_name            = sideset_name;
  localMeshInfo.cell_to_face            = cell_to_face;
  localMeshInfo.face_to_cells           = face_to_cells;      // faces
  localMeshInfo.face_to_lidx            = face_to_localidx;
  localMeshInfo.owned_cells             = owned_cells;        // cells
  localMeshInfo.ghstd_cells             = ghstd_cells;
  localMeshInfo.virtual_cells           = virtual_cells;
  localMeshInfo.owned_vertices          = owned_vertices;     // vertices
  localMeshInfo.ghstd_vertices          = ghstd_vertices;
  localMeshInfo.cell_topology           = mesh.getCellTopology(element_block_name);

  return localMeshInfo;

}





}











// Explicit instantiation
template
panzer::LocalMeshInfo<int,int>
panzer_stk::generateLocalMeshInfo<int,int>(const panzer_stk::STK_Interface & mesh,
                      const std::string & element_block_name);

#ifndef PANZER_ORDINAL64_IS_INT
template
panzer::LocalMeshInfo<int,panzer::Ordinal64>
panzer_stk::generateLocalMeshInfo<int,panzer::Ordinal64>(const panzer_stk::STK_Interface & mesh,
                      const std::string & element_block_name);
#endif

template
panzer::LocalMeshInfo<int,int>
panzer_stk::generateLocalSidesetInfo<int,int>(const panzer_stk::STK_Interface & mesh,
                      const std::string & element_block_name,
                      const std::string & sideset_name);

#ifndef PANZER_ORDINAL64_IS_INT
template
panzer::LocalMeshInfo<int,panzer::Ordinal64>
panzer_stk::generateLocalSidesetInfo<int,panzer::Ordinal64>(const panzer_stk::STK_Interface & mesh,
                      const std::string & element_block_name,
                      const std::string & sideset_name);
#endif

