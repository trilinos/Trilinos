// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

#include "Panzer_FieldPattern.hpp"
#include "Panzer_NodalFieldPattern.hpp"
#include "Panzer_EdgeFieldPattern.hpp"
#include "Panzer_FaceFieldPattern.hpp"
#include "Panzer_ElemFieldPattern.hpp"

#include "Panzer_ConnManager.hpp"

#include "Phalanx_KokkosDeviceTypes.hpp"

#include "Teuchos_Assert.hpp"
#include "Teuchos_OrdinalTraits.hpp"

#include "Tpetra_Import.hpp"

#include <string>
#include <map>
#include <vector>
#include <unordered_set>

namespace panzer_stk
{

// No external access
namespace
{

/** This method takes a cell importer (owned to ghstd) and communicates nodes
  * of the ghstd elements.
  */
Kokkos::DynRankView<double,PHX::Device>
buildGhostedNodes(const Tpetra::Import<int,panzer::GlobalOrdinal,panzer::TpetraNodeType> & importer,
                  Kokkos::DynRankView<const double,PHX::Device> owned_nodes)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef Tpetra::MultiVector<double,int,panzer::GlobalOrdinal,panzer::TpetraNodeType> mvec_type;

  size_t owned_cell_cnt = importer.getSourceMap()->getLocalNumElements();
  size_t ghstd_cell_cnt = importer.getTargetMap()->getLocalNumElements();
  int nodes_per_cell    = owned_nodes.extent(1);
  int space_dim         = owned_nodes.extent(2);

  TEUCHOS_ASSERT(owned_nodes.extent(0)==owned_cell_cnt);

  // build node multivector
  RCP<mvec_type> owned_nodes_mv = rcp(new mvec_type(importer.getSourceMap(),nodes_per_cell*space_dim));
  RCP<mvec_type> ghstd_nodes_mv = rcp(new mvec_type(importer.getTargetMap(),nodes_per_cell*space_dim));

  {
    auto owned_nodes_view = owned_nodes_mv->getLocalViewDevice(Tpetra::Access::OverwriteAll);
    Kokkos::parallel_for(owned_cell_cnt, KOKKOS_LAMBDA (size_t i) {
      int l = 0;
      for(int j=0;j<nodes_per_cell;j++)
        for(int k=0;k<space_dim;k++,l++)
          owned_nodes_view(i,l) = owned_nodes(i,j,k);
      });
  }

  // communicate ghstd nodes
  ghstd_nodes_mv->doImport(*owned_nodes_mv,importer,Tpetra::INSERT);

  // copy multivector into ghstd node structure
  Kokkos::DynRankView<double,PHX::Device> ghstd_nodes("ghstd_nodes",ghstd_cell_cnt,nodes_per_cell,space_dim);
  {
    auto ghstd_nodes_view = ghstd_nodes_mv->getLocalViewDevice(Tpetra::Access::ReadOnly);
    Kokkos::parallel_for(ghstd_cell_cnt, KOKKOS_LAMBDA (size_t i) {
      int l = 0;
      for(int j=0;j<nodes_per_cell;j++)
        for(int k=0;k<space_dim;k++,l++)
          ghstd_nodes(i,j,k) = ghstd_nodes_view(i,l);
      } );
    Kokkos::fence();
  }

  return ghstd_nodes;
} // end build ghstd nodes
void
setupLocalMeshBlockInfo(const panzer_stk::STK_Interface & mesh,
                        panzer::ConnManager & conn,
                        const panzer::LocalMeshInfo & mesh_info,
                        const std::string & element_block_name,
                        panzer::LocalMeshBlockInfo & block_info)
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

  std::vector<panzer::LocalOrdinal> owned_block_cells;
  auto local_cells_h = Kokkos::create_mirror_view(mesh_info.local_cells);
  Kokkos::deep_copy(local_cells_h, mesh_info.local_cells);
  for(int parent_owned_cell=0;parent_owned_cell<num_parent_owned_cells;++parent_owned_cell){
    const panzer::LocalOrdinal local_cell = local_cells_h(parent_owned_cell);
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
    panzer::partitioning_utilities::setupSubLocalMeshInfo(mesh_info, owned_block_cells, block_info);
  }
}


void
setupLocalMeshSidesetInfo(const panzer_stk::STK_Interface & mesh,
                          panzer::ConnManager& /* conn */,
                          const panzer::LocalMeshInfo & mesh_info,
                          const std::string & element_block_name,
                          const std::string & sideset_name,
                          panzer::LocalMeshSidesetInfo & sideset_info)
{

  // We use these typedefs to make the algorithm slightly more clear
  using LocalOrdinal = panzer::LocalOrdinal;
  using ParentOrdinal = int;
  using SubcellOrdinal = short;

  // This function identifies all cells in mesh_info that belong to element_block_name
  // and creates a block_info from it.

  // This is a list of all entities that make up the 'side'
  std::vector<stk::mesh::Entity> side_entities;

  // Grab the side entities associated with this sideset on the mesh
  // Note: Throws exception if element block or sideset doesn't exist
  try{

    mesh.getAllSides(sideset_name, element_block_name, side_entities);

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
  std::map<ParentOrdinal,std::vector<SubcellOrdinal> > owned_parent_cell_to_subcell_indexes;
  {

    // This is the subcell dimension we are trying to line up on the sideset
    const size_t face_subcell_dimension = static_cast<size_t>(mesh.getCellTopology(element_block_name)->getDimension()-1);

    // List of local subcell indexes indexed by element:
    // For example: a Tet (element) would have
    //  - 4 triangular faces (subcell_index 0-3, subcell_dimension=2)
    //  - 6 edges (subcell_index 0-5, subcell_dimension=1)
    //  - 4 nodes (subcell_index 0-3, subcell_dimension=0)
    // Another example: a Line (element) would have
    //  - 2 nodes (subcell_index 0-1, subcell_dimension=0)
    // The nodes coincide with the element vertices for these first order examples
    std::vector<stk::mesh::Entity> elements;
    std::vector<size_t> subcell_indexes;
    std::vector<size_t> subcell_dimensions;
    panzer_stk::workset_utils::getSideElementCascade(mesh, element_block_name, side_entities, subcell_dimensions, subcell_indexes, elements);
    const size_t num_elements = subcell_dimensions.size();

    // We need a fast lookup for maping local indexes to parent indexes
    std::unordered_map<LocalOrdinal,ParentOrdinal> local_to_parent_lookup;
    auto local_cells_h = Kokkos::create_mirror_view(mesh_info.local_cells);
    Kokkos::deep_copy(local_cells_h, mesh_info.local_cells);
    for(ParentOrdinal parent_index=0; parent_index<static_cast<ParentOrdinal>(mesh_info.local_cells.extent(0)); ++parent_index)
      local_to_parent_lookup[local_cells_h(parent_index)] = parent_index;

    // Add the subcell indexes for each element on the sideset
    // TODO: There is a lookup call here to map from local indexes to parent indexes which slows things down. Maybe there is a way around that
    for(size_t element_index=0; element_index<num_elements; ++element_index) {
      const size_t subcell_dimension = subcell_dimensions[element_index];

      // Add subcell to map
      if(subcell_dimension == face_subcell_dimension){
        const SubcellOrdinal subcell_index = static_cast<SubcellOrdinal>(subcell_indexes[element_index]);
        const LocalOrdinal element_local_index = static_cast<LocalOrdinal>(mesh.elementLocalId(elements[element_index]));

        // Look up the parent cell index using the local cell index
        const auto itr = local_to_parent_lookup.find(element_local_index);
        TEUCHOS_ASSERT(itr!= local_to_parent_lookup.end());
        const ParentOrdinal element_parent_index = itr->second;

        owned_parent_cell_to_subcell_indexes[element_parent_index].push_back(subcell_index);
      }
    }
  }

  // We now know the mapping of parent cell indexes to subcell indexes touching the sideset

  const panzer::LocalOrdinal num_owned_cells = owned_parent_cell_to_subcell_indexes.size();

  sideset_info.element_block_name = element_block_name;
  sideset_info.sideset_name = sideset_name;
  sideset_info.cell_topology = mesh.getCellTopology(element_block_name);

  sideset_info.num_owned_cells = num_owned_cells;

  struct face_t{
    face_t(const ParentOrdinal c0,
           const ParentOrdinal c1,
           const SubcellOrdinal sc0,
           const SubcellOrdinal sc1)
    {
      cell_0=c0;
      cell_1=c1;
      subcell_index_0=sc0;
      subcell_index_1=sc1;
    }
    ParentOrdinal cell_0;
    ParentOrdinal cell_1;
    SubcellOrdinal subcell_index_0;
    SubcellOrdinal subcell_index_1;
  };


  // Figure out how many cells on the other side of the sideset are ghost or virtual
  std::unordered_set<panzer::LocalOrdinal> owned_parent_cells_set, ghstd_parent_cells_set, virtual_parent_cells_set;
  std::vector<face_t> faces;
  {
    auto cell_to_faces_h = Kokkos::create_mirror_view(mesh_info.cell_to_faces);
    auto face_to_cells_h = Kokkos::create_mirror_view(mesh_info.face_to_cells);
    auto face_to_lidx_h = Kokkos::create_mirror_view(mesh_info.face_to_lidx);
    Kokkos::deep_copy(cell_to_faces_h, mesh_info.cell_to_faces);
    Kokkos::deep_copy(face_to_cells_h, mesh_info.face_to_cells);
    Kokkos::deep_copy(face_to_lidx_h, mesh_info.face_to_lidx);
    panzer::LocalOrdinal parent_virtual_cell_offset = mesh_info.num_owned_cells + mesh_info.num_ghstd_cells;
    for(const auto & local_cell_index_pair : owned_parent_cell_to_subcell_indexes){

      const ParentOrdinal parent_cell = local_cell_index_pair.first;
      const auto & subcell_indexes = local_cell_index_pair.second;

      owned_parent_cells_set.insert(parent_cell);

      for(const SubcellOrdinal & subcell_index : subcell_indexes){

        const LocalOrdinal face = cell_to_faces_h(parent_cell, subcell_index);
        const LocalOrdinal face_other_side = (face_to_cells_h(face,0) == parent_cell) ? 1 : 0;

        TEUCHOS_ASSERT(subcell_index == face_to_lidx_h(face, 1-face_other_side));

        const ParentOrdinal other_side_cell = face_to_cells_h(face, face_other_side);
        const SubcellOrdinal other_side_subcell_index = face_to_lidx_h(face, face_other_side);

        faces.push_back(face_t(parent_cell, other_side_cell, subcell_index, other_side_subcell_index));

        if(other_side_cell >= parent_virtual_cell_offset){
          virtual_parent_cells_set.insert(other_side_cell);
        } else {
          ghstd_parent_cells_set.insert(other_side_cell);
        }
      }
    }
  }

  std::vector<ParentOrdinal> all_cells;
  all_cells.insert(all_cells.end(),owned_parent_cells_set.begin(),owned_parent_cells_set.end());
  all_cells.insert(all_cells.end(),ghstd_parent_cells_set.begin(),ghstd_parent_cells_set.end());
  all_cells.insert(all_cells.end(),virtual_parent_cells_set.begin(),virtual_parent_cells_set.end());

  sideset_info.num_ghstd_cells = ghstd_parent_cells_set.size();
  sideset_info.num_virtual_cells = virtual_parent_cells_set.size();

  const LocalOrdinal num_real_cells = sideset_info.num_owned_cells + sideset_info.num_ghstd_cells;
  const LocalOrdinal num_total_cells = num_real_cells + sideset_info.num_virtual_cells;
  const LocalOrdinal num_nodes_per_cell = mesh_info.cell_nodes.extent(1);
  const LocalOrdinal num_dims = mesh_info.cell_nodes.extent(2);

  // Copy local indexes, global indexes, and cell nodes to sideset info
  {
    sideset_info.global_cells = PHX::View<panzer::GlobalOrdinal*>("global_cells", num_total_cells);
    sideset_info.local_cells = PHX::View<LocalOrdinal*>("local_cells", num_total_cells);
    sideset_info.cell_nodes = PHX::View<double***>("cell_nodes", num_total_cells, num_nodes_per_cell, num_dims);
    Kokkos::deep_copy(sideset_info.cell_nodes,0.);

    typename PHX::View<ParentOrdinal*>::HostMirror all_cells_h("all_cells_h", num_total_cells);
    PHX::View<ParentOrdinal*> all_cells_d("all_cells_d", num_total_cells);
    for(LocalOrdinal i=0; i<num_total_cells; ++i)
      all_cells_h(i) = all_cells[i];
    Kokkos::deep_copy(all_cells_d, all_cells_h);
    Kokkos::parallel_for(num_total_cells, KOKKOS_LAMBDA (LocalOrdinal i) {
      const ParentOrdinal parent_cell = all_cells_d(i);
      sideset_info.local_cells(i) = mesh_info.local_cells(parent_cell);
      sideset_info.global_cells(i) = mesh_info.global_cells(parent_cell);
      for(LocalOrdinal j=0; j<num_nodes_per_cell; ++j)
	for(LocalOrdinal k=0; k<num_dims; ++k)
	  sideset_info.cell_nodes(i,j,k) = mesh_info.cell_nodes(parent_cell,j,k);
      });
  }

  // Now we have to set the connectivity for the faces.

  const LocalOrdinal num_faces = faces.size();
  const LocalOrdinal num_faces_per_cell = mesh_info.cell_to_faces.extent(1);

  sideset_info.face_to_cells = PHX::View<LocalOrdinal*[2]>("face_to_cells", num_faces);
  sideset_info.face_to_lidx = PHX::View<LocalOrdinal*[2]>("face_to_lidx", num_faces);
  sideset_info.cell_to_faces = PHX::View<LocalOrdinal**>("cell_to_faces", num_total_cells, num_faces_per_cell);
  auto cell_to_faces_h = Kokkos::create_mirror_view(sideset_info.cell_to_faces);
  auto face_to_cells_h = Kokkos::create_mirror_view(sideset_info.face_to_cells);
  auto face_to_lidx_h = Kokkos::create_mirror_view(sideset_info.face_to_lidx);

  // Default the system with invalid cell index - this will be most of the entries
  Kokkos::deep_copy(cell_to_faces_h, -1);

  // Lookup for sideset cell index given parent cell index
  std::unordered_map<ParentOrdinal,ParentOrdinal> parent_to_sideset_lookup;
  for(unsigned int i=0; i<all_cells.size(); ++i)
    parent_to_sideset_lookup[all_cells[i]] = i;

  for(int face_index=0;face_index<num_faces;++face_index){
    const face_t & face = faces[face_index];
    const ParentOrdinal & parent_cell_0 = face.cell_0;
    const ParentOrdinal & parent_cell_1 = face.cell_1;

    // Convert the parent cell index into a sideset cell index
    const auto itr_0 = parent_to_sideset_lookup.find(parent_cell_0);
    TEUCHOS_ASSERT(itr_0 != parent_to_sideset_lookup.end());
    const ParentOrdinal sideset_cell_0 = itr_0->second;

    const auto itr_1 = parent_to_sideset_lookup.find(parent_cell_1);
    TEUCHOS_ASSERT(itr_1 != parent_to_sideset_lookup.end());
    const ParentOrdinal sideset_cell_1 = itr_1->second;

//    const ParentOrdinal sideset_cell_0 = std::distance(all_cells.begin(), std::find(all_cells.begin(), all_cells.end(), parent_cell_0));
//    const ParentOrdinal sideset_cell_1 = std::distance(all_cells.begin(), std::find(all_cells.begin()+num_owned_cells, all_cells.end(), parent_cell_1));

    face_to_cells_h(face_index,0) = sideset_cell_0;
    face_to_cells_h(face_index,1) = sideset_cell_1;

    face_to_lidx_h(face_index,0) = face.subcell_index_0;
    face_to_lidx_h(face_index,1) = face.subcell_index_1;

    cell_to_faces_h(sideset_cell_0,face.subcell_index_0) = face_index;
    cell_to_faces_h(sideset_cell_1,face.subcell_index_1) = face_index;

  }
  Kokkos::deep_copy(sideset_info.cell_to_faces, cell_to_faces_h);
  Kokkos::deep_copy(sideset_info.face_to_cells, face_to_cells_h);
  Kokkos::deep_copy(sideset_info.face_to_lidx,  face_to_lidx_h );

}

} // namespace

Teuchos::RCP<panzer::LocalMeshInfo>
generateLocalMeshInfo(const panzer_stk::STK_Interface & mesh)
{
  TEUCHOS_FUNC_TIME_MONITOR_DIFF("panzer_stk::generateLocalMeshInfo",GenerateLocalMeshInfo);

  using Teuchos::RCP;
  using Teuchos::rcp;

  //typedef Tpetra::CrsMatrix<int,panzer::LocalOrdinal,panzer::GlobalOrdinal> crs_type;
  typedef Tpetra::Map<panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType> map_type;
  typedef Tpetra::Import<panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType> import_type;
  //typedef Tpetra::MultiVector<double,panzer::LocalOrdinal,panzer::GlobalOrdinal> mvec_type;
  //typedef Tpetra::MultiVector<panzer::GlobalOrdinal,panzer::LocalOrdinal,panzer::GlobalOrdinal> ordmvec_type;

  auto mesh_info_rcp = Teuchos::rcp(new panzer::LocalMeshInfo);
  auto & mesh_info = *mesh_info_rcp;

  // Make sure the STK interface is valid
  TEUCHOS_ASSERT(mesh.isInitialized());

  // This is required by some of the STK stuff
  TEUCHOS_ASSERT(typeid(panzer::LocalOrdinal) == typeid(int));

  Teuchos::RCP<const Teuchos::Comm<int> > comm = mesh.getComm();

  // This horrible line of code is required since the connection manager only takes rcps of a mesh
  RCP<const panzer_stk::STK_Interface> mesh_rcp = Teuchos::rcpFromRef(mesh);
  // We're allowed to do this since the connection manager only exists in this scope... even though it is also an RCP...

  // extract topology handle
  RCP<panzer::ConnManager> conn_rcp = rcp(new panzer_stk::STKConnManager(mesh_rcp));
  panzer::ConnManager & conn = *conn_rcp;

  PHX::View<panzer::GlobalOrdinal*> owned_cells, ghost_cells, virtual_cells;
  panzer::fillLocalCellIDs(comm, conn, owned_cells, ghost_cells, virtual_cells);

  // build cell maps
  /////////////////////////////////////////////////////////////////////

  RCP<map_type> owned_cell_map = rcp(new map_type(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),owned_cells,0,comm));
  RCP<map_type> ghstd_cell_map = rcp(new map_type(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),ghost_cells,0,comm));

  // build importer: cell importer, owned to ghstd
  RCP<import_type> cellimport_own2ghst = rcp(new import_type(owned_cell_map,ghstd_cell_map));

  // read all the nodes associated with these elements, get ghstd contributions
  /////////////////////////////////////////////////////////////////////

  // TODO: This all needs to be rewritten for when element blocks have different cell topologies
  std::vector<std::string> element_block_names;
  mesh.getElementBlockNames(element_block_names);

  const shards::CellTopology & cell_topology = *(mesh.getCellTopology(element_block_names[0]));

  const int space_dim = cell_topology.getDimension();
  const int nodes_per_cell = cell_topology.getNodeCount();
  const int faces_per_cell = cell_topology.getSubcellCount(space_dim-1);

  Kokkos::DynRankView<double,PHX::Device> owned_nodes("owned_nodes",owned_cells.extent(0),nodes_per_cell,space_dim);
  {
    std::vector<std::size_t> localCells(owned_cells.extent(0),Teuchos::OrdinalTraits<std::size_t>::invalid());
    for(size_t i=0;i<localCells.size();i++)
      localCells[i] = i;
    mesh.getElementNodesNoResize(localCells,owned_nodes);
  }

  // this builds a ghstd node array
  Kokkos::DynRankView<double,PHX::Device> ghstd_nodes = buildGhostedNodes(*cellimport_own2ghst,owned_nodes);

  // build edge to cell neighbor mapping
  //////////////////////////////////////////////////////////////////

  auto owned_cells_h = Kokkos::create_mirror_view(owned_cells);
  auto ghost_cells_h = Kokkos::create_mirror_view(ghost_cells);
  Kokkos::deep_copy(owned_cells_h, owned_cells);
  Kokkos::deep_copy(ghost_cells_h, ghost_cells);
  std::unordered_map<panzer::GlobalOrdinal,int> global_to_local;
  global_to_local[-1] = -1; // this is the "no neighbor" flag
  for(size_t i=0;i<owned_cells.extent(0);i++)
    global_to_local[owned_cells_h(i)] = i;
  for(size_t i=0;i<ghost_cells.extent(0);i++)
    global_to_local[ghost_cells_h(i)] = i+Teuchos::as<int>(owned_cells.extent(0));

  // this class comes from Mini-PIC and Matt B
  RCP<panzer::FaceToElement<panzer::LocalOrdinal,panzer::GlobalOrdinal> > faceToElement = rcp(new panzer::FaceToElement<panzer::LocalOrdinal,panzer::GlobalOrdinal>());
  faceToElement->initialize(conn, comm);
  auto elems_by_face = faceToElement->getFaceToElementsMap();
  auto face_to_lidx  = faceToElement->getFaceToCellLocalIdxMap();

  // We also need to consider faces that connect to cells that do not exist, but are needed for boundary conditions
  // We dub them virtual cell since there should be no geometry associated with them, or topology really
  // They exist only for datastorage so that they are consistent with 'real' cells from an algorithm perspective

  // Each virtual face (face linked to a '-1' cell) requires a virtual cell (i.e. turn the '-1' into a virtual cell)
  // Virtual cells are those that do not exist but are connected to an owned cell
  // Note - in the future, ghosted cells will also need to connect to virtual cells at boundary conditions, but for the moment we will ignore this.

  // Iterate over all faces and identify the faces connected to a potential virtual cell

  const panzer::LocalOrdinal num_owned_cells = owned_cells.extent(0);
  const panzer::LocalOrdinal num_ghstd_cells = ghost_cells.extent(0);
  const panzer::LocalOrdinal num_virtual_cells = virtual_cells.extent(0);

  // total cells and faces include owned, ghosted, and virtual
  const panzer::LocalOrdinal num_real_cells = num_owned_cells + num_ghstd_cells;
  const panzer::LocalOrdinal num_total_cells = num_real_cells + num_virtual_cells;
  const panzer::LocalOrdinal num_total_faces = elems_by_face.extent(0);

  // Lookup cells connected to a face
  PHX::View<panzer::LocalOrdinal*[2]> face_to_cells("face_to_cells",num_total_faces);

  // Lookup local face indexes given cell and left/right state (0/1)
  PHX::View<panzer::LocalOrdinal*[2]> face_to_localidx("face_to_localidx",num_total_faces);

  // Lookup face index given a cell and local face index
  PHX::View<panzer::LocalOrdinal**> cell_to_face("cell_to_face",num_total_cells,faces_per_cell);

  // Transfer information from 'faceToElement' datasets to local arrays
  {
    PANZER_FUNC_TIME_MONITOR_DIFF("Transer faceToElement to local",TransferFaceToElementLocal);

    int virtual_cell_index = num_real_cells;
    auto elems_by_face_h = Kokkos::create_mirror_view(elems_by_face);
    auto face_to_lidx_h = Kokkos::create_mirror_view(face_to_lidx);
    auto face_to_cells_h = Kokkos::create_mirror_view(face_to_cells);
    auto face_to_localidx_h = Kokkos::create_mirror_view(face_to_localidx);
    auto cell_to_face_h = Kokkos::create_mirror_view(cell_to_face);
    Kokkos::deep_copy(elems_by_face_h, elems_by_face);
    Kokkos::deep_copy(face_to_lidx_h, face_to_lidx);
    // initialize with negative one cells that are not associated with a face
    Kokkos::deep_copy(cell_to_face_h, -1);
    for(size_t f=0;f<elems_by_face.extent(0);f++) {

      const panzer::GlobalOrdinal global_c0 = elems_by_face_h(f,0);
      const panzer::GlobalOrdinal global_c1 = elems_by_face_h(f,1);

      // make sure that no bonus cells get in here
      TEUCHOS_ASSERT(global_to_local.find(global_c0)!=global_to_local.end());
      TEUCHOS_ASSERT(global_to_local.find(global_c1)!=global_to_local.end());

      auto c0 = global_to_local[global_c0];
      auto lidx0 = face_to_lidx_h(f,0);
      auto c1 = global_to_local[global_c1];
      auto lidx1 = face_to_lidx_h(f,1);

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
      cell_to_face_h(c0,lidx0) = f;


      // Right cell
      if(c1 < 0){
        // Virtual cell - create it!
        c1 = virtual_cell_index++;

        // We need the subcell_index to line up between real and virtual cell
        // This way the face has the same geometry... though the face normal
        // will point in the wrong direction
        lidx1 = lidx0;
      }
      cell_to_face_h(c1,lidx1) = f;

      // Faces point from low cell index to high cell index
      if(c0<c1){
        face_to_cells_h(f,0) = c0;
        face_to_localidx_h(f,0) = lidx0;
        face_to_cells_h(f,1) = c1;
        face_to_localidx_h(f,1) = lidx1;
      } else {
        face_to_cells_h(f,0) = c1;
        face_to_localidx_h(f,0) = lidx1;
        face_to_cells_h(f,1) = c0;
        face_to_localidx_h(f,1) = lidx0;
      }

      // We should avoid having two virtual cells linked together.
      TEUCHOS_ASSERT(c0<num_real_cells or c1<num_real_cells)

    }
    Kokkos::deep_copy(face_to_cells,  face_to_cells_h);
    Kokkos::deep_copy(face_to_localidx, face_to_localidx_h);
    Kokkos::deep_copy(cell_to_face,   cell_to_face_h);
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
    mesh_info.subcell_dimension       = space_dim;
    mesh_info.subcell_index           = -1;
    mesh_info.has_connectivity        = true;

    mesh_info.num_owned_cells = owned_cells.extent(0);
    mesh_info.num_ghstd_cells = ghost_cells.extent(0);
    mesh_info.num_virtual_cells = virtual_cells.extent(0);

    mesh_info.global_cells = PHX::View<panzer::GlobalOrdinal*>("global_cell_indices",num_total_cells);
    mesh_info.local_cells = PHX::View<panzer::LocalOrdinal*>("local_cell_indices",num_total_cells);

    Kokkos::parallel_for(num_owned_cells,KOKKOS_LAMBDA (int i) {
      mesh_info.global_cells(i) = owned_cells(i);
      mesh_info.local_cells(i) = i;
      });

    Kokkos::parallel_for(num_ghstd_cells,KOKKOS_LAMBDA (int i) {
      mesh_info.global_cells(i+num_owned_cells) = ghost_cells(i);
      mesh_info.local_cells(i+num_owned_cells) = i+num_owned_cells;
      });

    Kokkos::parallel_for(num_virtual_cells,KOKKOS_LAMBDA (int i) {
      mesh_info.global_cells(i+num_real_cells) = virtual_cells(i);
      mesh_info.local_cells(i+num_real_cells) = i+num_real_cells;
      });

    mesh_info.cell_nodes = PHX::View<double***>("cell_nodes",num_total_cells,nodes_per_cell,space_dim);

    // Initialize coordinates to zero
    Kokkos::deep_copy(mesh_info.cell_nodes, 0.);

    Kokkos::parallel_for(num_owned_cells,KOKKOS_LAMBDA (int i) {
      for(int j=0;j<nodes_per_cell;++j){
        for(int k=0;k<space_dim;++k){
          mesh_info.cell_nodes(i,j,k) = owned_nodes(i,j,k);
        }
      }
      });

    Kokkos::parallel_for(num_ghstd_cells,KOKKOS_LAMBDA (int i) {
      for(int j=0;j<nodes_per_cell;++j){
        for(int k=0;k<space_dim;++k){
          mesh_info.cell_nodes(i+num_owned_cells,j,k) = ghstd_nodes(i,j,k);
        }
      }
      });

    // This will backfire at some point, but we're going to make the virtual cell have the same geometry as the cell it interfaces with
    // This way we can define a virtual cell geometry without extruding the face outside of the domain
    // TODO BWR Certainly, this is an issue for curved meshes
    {
      PANZER_FUNC_TIME_MONITOR_DIFF("Assign geometry traits",AssignGeometryTraits);
      Kokkos::parallel_for(num_virtual_cells,KOKKOS_LAMBDA (int i) {
        const panzer::LocalOrdinal virtual_cell = i+num_real_cells;
        for(int local_face=0; local_face<faces_per_cell; ++local_face){
          const panzer::LocalOrdinal face = cell_to_face(virtual_cell, local_face);
          if(face >= 0){
            const panzer::LocalOrdinal other_side = (face_to_cells(face, 0) == virtual_cell) ? 1 : 0;
            const panzer::LocalOrdinal real_cell = face_to_cells(face,other_side);
            for(int j=0;j<nodes_per_cell;++j){
              for(int k=0;k<space_dim;++k){
                mesh_info.cell_nodes(virtual_cell,j,k) = mesh_info.cell_nodes(real_cell,j,k);
              }
            }
            break;
          }
        }
	});
      
    }
  }

  // Setup element blocks and sidesets
  std::vector<std::string> sideset_names;
  mesh.getSidesetNames(sideset_names);

  for(const std::string & element_block_name : element_block_names){
    PANZER_FUNC_TIME_MONITOR_DIFF("Set up setupLocalMeshBlockInfo",SetupLocalMeshBlockInfo);
    panzer::LocalMeshBlockInfo & block_info = mesh_info.element_blocks[element_block_name];
    setupLocalMeshBlockInfo(mesh, conn, mesh_info, element_block_name, block_info);
    block_info.subcell_dimension = space_dim;
    block_info.subcell_index = -1;
    block_info.has_connectivity = true;

    // Setup sidesets
    for(const std::string & sideset_name : sideset_names){
      PANZER_FUNC_TIME_MONITOR_DIFF("Setup LocalMeshSidesetInfo",SetupLocalMeshSidesetInfo);
      panzer::LocalMeshSidesetInfo & sideset_info = mesh_info.sidesets[element_block_name][sideset_name];
      setupLocalMeshSidesetInfo(mesh, conn, mesh_info, element_block_name, sideset_name, sideset_info);
      sideset_info.subcell_dimension = space_dim;
      sideset_info.subcell_index = -1;
      sideset_info.has_connectivity = true;
    }

  }

  return mesh_info_rcp;

}

}
