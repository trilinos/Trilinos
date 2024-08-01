// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_Relations.hpp"
#include "Panzer_NodalFieldPattern.hpp"
#include "Panzer_EdgeFieldPattern.hpp"
#include "Panzer_FaceFieldPattern.hpp"
#include "Panzer_ElemFieldPattern.hpp"
#include "Teuchos_DefaultMpiComm.hpp"

namespace panzer {

FaceToElems::FaceToElems(Teuchos::RCP<panzer::ConnManager> conn) :
  conn_(conn){

  std::vector<std::string> block_ids;
  conn_->getElementBlockIds(block_ids);
  num_blocks_ = block_ids.size();
  TEUCHOS_ASSERT(num_blocks_ > 0);

  conn_->getElementBlockTopologies(element_block_topologies_);
  dimension_ = element_block_topologies_[0].getDimension();

  for (int i=0;i<num_blocks_; ++i) {
    TEUCHOS_ASSERT(
      static_cast<int>(element_block_topologies_[i].getDimension()) ==
      dimension_);
    // Currently this only works if all the topologies are the same due to field pattern stuff
    TEUCHOS_ASSERT(element_block_topologies_[i] ==element_block_topologies_[0] );
  }

  std::set<GlobalOrdinal> face_gids;


  Teuchos::RCP<const Map> elem_map,face_map;

  Teuchos::RCP<const Teuchos::Comm<int>> comm(new Teuchos::MpiComm< int>(MPI_COMM_WORLD));


  if ( dimension_ == 1 ) {
    panzer::EdgeFieldPattern edge_pattern(element_block_topologies_[0]);
    conn_->buildConnectivity(edge_pattern);
  } else if ( dimension_ == 2 ){
    panzer::FaceFieldPattern face_pattern(element_block_topologies_[0]);
    conn_->buildConnectivity(face_pattern);
  } else {
    panzer::ElemFieldPattern elem_pattern(element_block_topologies_[0]);
    conn_->buildConnectivity(elem_pattern);
  }

  std::vector<GlobalOrdinal> elem_gids;
  {
    total_elements_ = 0;
    for (int iblk=0; iblk< num_blocks_; ++iblk){
      auto foo = conn_->getElementBlock(block_ids[iblk]);
      total_elements_ += conn_->getElementBlock(block_ids[iblk]).size();
    }
    elem_gids.resize(total_elements_);
    for (int iblk=0; iblk< num_blocks_; ++iblk){
      const std::vector<LocalOrdinal> &localIDs = conn_->getElementBlock(block_ids[iblk]);
      for (unsigned id=0;id<localIDs.size(); ++id) {
        int n_conn = conn_->getConnectivitySize(localIDs[id]);
        const auto * connectivity = conn_->getConnectivity(localIDs[id]);
        TEUCHOS_ASSERT(n_conn==1);
        elem_gids[localIDs[id]] = connectivity[0];
      }
    }
  }
  elem_map = Teuchos::RCP<Map>( new Map(Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(), &elem_gids[0], elem_gids.size(), 0, comm ));

  // Now build one level down
  if ( dimension_ == 1 ) {
    panzer::NodalFieldPattern node_pattern(element_block_topologies_[0]);
    conn_->buildConnectivity(node_pattern);
  } else if ( dimension_ == 2 ){
    panzer::EdgeFieldPattern edge_pattern(element_block_topologies_[0]);
    conn_->buildConnectivity(edge_pattern);
  } else {
    panzer::FaceFieldPattern face_pattern(element_block_topologies_[0]);
    conn_->buildConnectivity(face_pattern);
  }



  total_elements_ = 0;
  for (int iblk=0; iblk< num_blocks_; ++iblk){
    auto foo = conn_->getElementBlock(block_ids[iblk]);
    total_elements_ += conn_->getElementBlock(block_ids[iblk]).size();
  }
  elem_to_face_.resize(total_elements_);

  for (int iblk=0; iblk< num_blocks_; ++iblk){
    const std::vector<LocalOrdinal> &localIDs = conn_->getElementBlock(block_ids[iblk]);
    for (unsigned id=0;id<localIDs.size(); ++id) {
      int n_conn = conn_->getConnectivitySize(localIDs[id]);
      const auto * connectivity = conn_->getConnectivity(localIDs[id]);
      elem_to_face_[id].resize(n_conn);
      for (int iconn=0;iconn<n_conn; ++iconn) {
        elem_to_face_[id][iconn] = connectivity[iconn];
        face_gids.insert(connectivity[iconn]);
      }
    }
  }


  std::vector<GlobalOrdinal> face_gids_vec;
  face_gids_vec.insert(face_gids_vec.begin(), face_gids.begin(), face_gids.end());

  face_map = Teuchos::RCP<const Map>(new Map( Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(), &face_gids_vec[0], face_gids_vec.size(), 0, comm ));


  Teuchos::RCP<const Map>  owned_face_map = Tpetra::createOneToOne(face_map);
  Graph graph(owned_face_map, 2);
  for (int ielem=0;ielem< static_cast<int>(elem_to_face_.size()); ++ielem)
    for (int iface=0; iface<static_cast<int>(elem_to_face_[ielem].size()); ++iface ) {
      GlobalOrdinal go=elem_map->getGlobalElement(ielem);
      graph.insertGlobalIndices(elem_to_face_[ielem][iface], 1, &go);
    }
  graph.globalAssemble();


  Import imp(owned_face_map, face_map);
  Graph graph_overlap(face_map, 2);
  graph_overlap.doImport(graph, imp, Tpetra::ADD);


  face_to_elem_ = PHX::View<GlobalOrdinal*[2]>("FaceToElems::face_to_elem_",face_map->getLocalNumElements());
  auto face_to_elem_h = Kokkos::create_mirror_view(face_to_elem_);
  num_boundary_faces_=0;
  for (int i(0); i < face_to_elem_.extent_int(0); ++i)
  {
    typename Graph::nonconst_global_inds_host_view_type indices("indices", 2);
    size_t num_ent;
    graph_overlap.getGlobalRowCopy(face_map->getGlobalElement(i), indices, num_ent);
    assert(num_ent == 2 || num_ent == 1);
    face_to_elem_h(i,0) = indices(0);
    if ( num_ent == 2)
      face_to_elem_h(i,1) = indices(1);
    else {
      face_to_elem_h(i,1) = -1;
      num_boundary_faces_++;
    }
  }
  Kokkos::deep_copy(face_to_elem_, face_to_elem_h);

  // Now we can get the nodal values
  std::vector<std::vector<int>> face_to_node;
  {
    panzer::NodalFieldPattern node_pattern(element_block_topologies_[0]);
    conn_->buildConnectivity(node_pattern);
    int num_faces = element_block_topologies_[0].getSubcellCount(element_block_topologies_[0].getDimension()-1);
    face_to_node.resize(num_faces);
    for (int i=0; i<num_faces; ++i) {
      int num_nodes = element_block_topologies_[0].getNodeCount(element_block_topologies_[0].getDimension()-1, i);
      face_to_node[i].resize(num_nodes);
      TEUCHOS_ASSERT(face_to_node[i].size() ==face_to_node[0].size() );
      for (int j=0;j<num_nodes; ++j)
        face_to_node[i][j] = element_block_topologies_[0].getNodeMap(element_block_topologies_[0].getDimension()-1, i, j);
    }

  }
  face_to_node_ = PHX::View<GlobalOrdinal**>("face_to_node", face_to_elem_.extent(0), face_to_node[0].size());
  auto face_to_node_h = Kokkos::create_mirror_view(face_to_node_);
  Kokkos::deep_copy(face_to_node_h, -1);
  for (int ielem=0;ielem< static_cast<int>(elem_to_face_.size()); ++ielem) {
    const auto * connectivity = conn_->getConnectivity(ielem);
    for (int iface=0; iface <static_cast<int>(elem_to_face_[ielem].size()); ++iface ) {
      for (int inode(0); inode < face_to_node_.extent_int(1); ++inode)
      {
        GlobalOrdinal g_face = elem_to_face_[ielem][iface];
        LocalOrdinal l_face = face_map->getLocalElement(g_face);
        face_to_node_h(l_face, inode) = connectivity[face_to_node[iface][inode]];
      }
    }
  }
  for (int i(0); i < face_to_node_.extent_int(0); ++i)
    for (int j(0); j < face_to_node_.extent_int(1); ++j)
      TEUCHOS_ASSERT(face_to_node_h(i,j) >=0);
  Kokkos::deep_copy(face_to_node_, face_to_node_h);
}

void FaceToElems::setNormals(Teuchos::RCP<std::vector<panzer::Workset> > worksets){

  // Now we can get the nodal values
  std::vector<std::vector<int>> face_to_node;
  {
    panzer::NodalFieldPattern node_pattern(element_block_topologies_[0]);
    conn_->buildConnectivity(node_pattern);
    int num_faces = element_block_topologies_[0].getSubcellCount(element_block_topologies_[0].getDimension()-1);
    face_to_node.resize(num_faces);
    for (int i=0; i<num_faces; ++i) {
      int num_nodes = element_block_topologies_[0].getNodeCount(element_block_topologies_[0].getDimension()-1, i);
      face_to_node[i].resize(num_nodes);
      TEUCHOS_ASSERT(face_to_node[i].size() ==face_to_node[0].size() );
      for (int j=0;j<num_nodes; ++j)
        face_to_node[i][j] = element_block_topologies_[0].getNodeMap(element_block_topologies_[0].getDimension()-1, i, j);
    }

  }
  face_normal_ = PHX::View<double***> ("FaceToElems::face_normal_", total_elements_, face_to_node.size(), dimension_ );
  face_centroid_ = PHX::View<double***> ("FaceToElems::face_centroid_", total_elements_, face_to_node.size(), dimension_ );
  auto face_normal_h = Kokkos::create_mirror_view(face_normal_);
  auto face_centroid_h = Kokkos::create_mirror_view(face_centroid_);

  int num_worksets = worksets->size();

  for (int nwkst=0; nwkst<num_worksets; ++nwkst){
    panzer::Workset &workset = (*worksets)[nwkst];
    auto coords = workset.cell_node_coordinates;
    auto coords_h = Kokkos::create_mirror_view(coords.get_static_view());
    Kokkos::deep_copy(coords_h, coords.get_static_view());
    int num_cells = workset.num_cells;
    // Compute the rough cell face centroid
    for (int c=0; c<num_cells; ++c) {
      for (int nface=0;nface <static_cast<int>(face_to_node.size()); ++nface) {
        std::vector<double> center(dimension_, 0.0);
        for (int nnode=0; nnode < static_cast<int>(face_to_node[nface].size()); ++nnode)
          for (int idim=0;idim<dimension_; ++idim)
            center[idim] += coords_h(c,face_to_node[nface][nnode], idim);
        for (int idim=0;idim<dimension_; ++idim)
          face_centroid_h(workset.cell_local_ids[c], nface, idim) =  center[idim]/face_to_node[nface].size();

      }
    }
    // Now lets compute the face normals
    typename PHX::View<double**>::HostMirror edges("temp::Edges",40, dimension_);  // overkill on 40 size
    for (int c=0; c<num_cells; ++c) {

      std::vector<double> center(3,0.);
      for (int nface=0;nface <static_cast<int>(face_to_node.size()); ++nface)
        for (int idim=0;idim<dimension_; ++idim)
          center[idim] += face_centroid_h(workset.cell_local_ids[c], nface, idim)/face_to_node.size();

      for (int nface=0;nface <static_cast<int>(face_to_node.size()); ++nface) {
        std::vector<double> normal(3,0);

        // Create centroid to node edges.
        for (int nnode=0; nnode < static_cast<int>(face_to_node[nface].size()); ++nnode) {
          for (int idim=0;idim<dimension_; ++idim)
            edges(nnode,idim) = coords_h(c,face_to_node[nface][nnode], idim) - face_centroid_h(workset.cell_local_ids[c], nface, idim);
        }

        std::vector<double> approx_normal(3);
        for (int idim=0;idim<dimension_; ++idim)
          approx_normal[idim] = center[idim]-face_centroid_h(workset.cell_local_ids[c], nface, idim);

        if ( dimension_ == 1) {
          normal[0] = 1;
        } else if (dimension_ == 2) {
          normal[0] = +edges(0,1);
          normal[1] = -edges(0,0);
        } else {
          // One can do a cross product on these to get a normal direction
          for (int nnode=0; nnode < static_cast<int>(face_to_node[nface].size()); ++nnode) {
            int n1=nnode, n2=(nnode+1)%face_to_node[nface].size();
            normal[0] += edges(n1,1)*edges(n2,2) - edges(n1,2)*edges(n2,1);
            normal[1] += edges(n1,2)*edges(n2,0) - edges(n1,0)*edges(n2,2);
            normal[2] += edges(n1,0)*edges(n2,1) - edges(n1,1)*edges(n2,0);
          }
        }
        // Now we can normalize the normal
        double normal_norm = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
        double sign = normal[0]*approx_normal[0]+normal[1]*approx_normal[1]+normal[2]*approx_normal[2];
        if ( sign < 0 )
          normal_norm *= -1;
        for (int idim=0;idim<dimension_; ++idim)
          face_normal_h(workset.cell_local_ids[c], nface, idim) = normal[idim]/normal_norm;

      }
    }
  }
  Kokkos::deep_copy(face_centroid_, face_centroid_h);
  Kokkos::deep_copy(face_normal_, face_normal_h);
}
void FaceToElems::getNormal(LocalOrdinal ielem, int iface, std::vector<double> &normal){
  TEUCHOS_ASSERT(ielem < face_normal_.extent_int(0));
  TEUCHOS_ASSERT(iface < face_normal_.extent_int(1));
  normal.resize(dimension_);
  auto face_normal_h = Kokkos::create_mirror_view(face_normal_);
  Kokkos::deep_copy(face_normal_h, face_normal_);

  for (int idim=0;idim<dimension_; ++idim)
    normal[idim] = face_normal_h(ielem, iface, idim);
}


}// panzer namespace
