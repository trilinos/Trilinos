/*
#@HEADER
# ************************************************************************
#
#                          Moertel FE Package
#                 Copyright (2015) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Glen Hansen (gahanse@sandia.gov)
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */


#include <gtest/gtest.h>

#include "mrk_default_kokkos_device_type.hpp"
#include "mrk_data_types.hpp"
#include "mrk_api_classes.hpp"
#include "mrk_interface_impl.hpp"
#include "mrk_manager_impl.hpp"
#include "mrk_manager_tester.hpp"
#include "morkon_unit_test_utils_TRI3.hpp"

//
// Use default_kokkos_device_t (from mrk_default_kokkos_device_type.hpp) for
// now.  When we start using devices that require Kokkos::initialize(..) and
// Kokkos::finalize(), we probably want to use TEST_F, which automatically calls
// the SetUpTestCase() and TearDownTestCase() static member functions in a
// test fixture class.  See, e.g., kokkos/core/unit_test/TestOpenMP.cpp.


TEST(morkon,just_check_if_it_compiles)
{
  using namespace morkon_exp;
  typedef Morkon_Manager<default_kokkos_device_t, 3, MRK_QUAD4>    default_manager_3d_t;
  typedef Teuchos::RCP< default_manager_3d_t >        default_manager_3d_ptr;
  typedef Interface<default_kokkos_device_t, 3, MRK_QUAD4>     default_interface_3d_t;
  typedef Teuchos::RCP< default_interface_3d_t >  default_interface_3d_ptr;
  default_manager_3d_ptr manager = default_manager_3d_t::MakeInstance(0, FACET_NORMAL_PROJECTION, 0);
  default_interface_3d_ptr interface = manager->create_interface(0,0);
  manager->commit_interfaces();
  Tpetra::CrsMatrix<> *dummy_D = 0;
  Tpetra::CrsMatrix<> *dummy_M = 0;
  manager->mortar_integrate(dummy_D, dummy_M);
}


TEST(morkon, interface_host_side_adapter_stores_data)
{
  using namespace morkon_exp;
  typedef Morkon_Manager<default_kokkos_device_t, 3, MRK_TRI3>  manager_3d_t;
  typedef Teuchos::RCP< manager_3d_t >                        manager_3d_ptr;
  typedef Interface<default_kokkos_device_t, 3, MRK_TRI3>     interface_3d_t;
  typedef Teuchos::RCP< interface_3d_t >                    interface_3d_ptr;

  manager_3d_ptr manager = manager_3d_t::MakeInstance(0, FACET_NORMAL_PROJECTION, 0);
  interface_3d_ptr interface = manager->create_interface(0,0);

  Mrk_2x2_aligned_TriangleInterfaceFixture tris_2x2(interface);

  const interface_3d_t::host_side_adapter_t *non_mortar_side_hsa =
        interface->get_HostSideAdapter(InterfaceBase::NON_MORTAR_SIDE);
  EXPECT_NE(static_cast<interface_3d_t::host_side_adapter_t *>(0), non_mortar_side_hsa);

  const interface_3d_t::host_side_adapter_t *mortar_side_hsa =
        interface->get_HostSideAdapter(InterfaceBase::MORTAR_SIDE);
  EXPECT_NE(static_cast<interface_3d_t::host_side_adapter_t *>(0), mortar_side_hsa);

  size_t node_i = 0;
  for (auto node_entry : non_mortar_side_hsa->m_nodes)
  {
    EXPECT_EQ(tris_2x2.node_gids[node_i], node_entry.first);
    EXPECT_EQ(tris_2x2.node_gids[node_i], node_entry.second.m_id);
    EXPECT_EQ(InterfaceBase::NON_MORTAR_SIDE, node_entry.second.m_side);
    for (size_t j = 0; j < tris_2x2.NodesPerFace; ++j) {
      EXPECT_EQ(tris_2x2.node_coords[node_i][j], node_entry.second.m_coords[j]);
    }
    ++node_i;
  }
  for (auto node_entry : mortar_side_hsa->m_nodes)
  {
    EXPECT_EQ(tris_2x2.node_gids[node_i], node_entry.first);
    EXPECT_EQ(tris_2x2.node_gids[node_i], node_entry.second.m_id);
    EXPECT_EQ(InterfaceBase::MORTAR_SIDE, node_entry.second.m_side);
    for (size_t j = 0; j < tris_2x2.NodesPerFace; ++j) {
      EXPECT_EQ(tris_2x2.node_coords[node_i][j], node_entry.second.m_coords[j]);
    }
    ++node_i;
  }
  const size_t num_nodes = tris_2x2.NumNodes;
  EXPECT_EQ(num_nodes, node_i);

  size_t face_i = 0;
  for (auto face_entry : non_mortar_side_hsa->m_faces)
  {
    EXPECT_EQ(tris_2x2.face_gids[face_i], face_entry.first);
    EXPECT_EQ(tris_2x2.node_gids[face_i], face_entry.second.m_id);
    EXPECT_EQ(InterfaceBase::NON_MORTAR_SIDE, face_entry.second.m_side);
    for (size_t j = 0; j < tris_2x2.NodesPerFace; ++j) {
      EXPECT_EQ(tris_2x2.face_node_gids[face_i][j], face_entry.second.m_nodes[j]);
    }
    ++face_i;
  }
  for (auto face_entry : mortar_side_hsa->m_faces)
  {
    EXPECT_EQ(tris_2x2.face_gids[face_i], face_entry.first);
    EXPECT_EQ(tris_2x2.node_gids[face_i], face_entry.second.m_id);
    EXPECT_EQ(InterfaceBase::MORTAR_SIDE, face_entry.second.m_side);
    for (size_t j = 0; j < tris_2x2.NodesPerFace; ++j) {
      EXPECT_EQ(tris_2x2.face_node_gids[face_i][j], face_entry.second.m_nodes[j]);
    }
    ++face_i;
  }
  const size_t num_faces = tris_2x2.NumFaces;
  EXPECT_EQ(num_faces, face_i);
}


TEST(morkon, manager_commit_interfaces_one_interface)
{
  using namespace morkon_exp;
  typedef Morkon_Manager_Tester<default_kokkos_device_t, 3, MRK_TRI3>  manager_3d_t;
  typedef Teuchos::RCP< manager_3d_t >                        manager_3d_ptr;
  typedef Interface<default_kokkos_device_t, 3, MRK_TRI3>     interface_3d_t;
  typedef Teuchos::RCP< interface_3d_t >                    interface_3d_ptr;

  manager_3d_ptr manager = manager_3d_t::MakeInstance(0, FACET_NORMAL_PROJECTION, 0);

  const int interface_id = 17;
  interface_3d_ptr interface = manager->create_interface(interface_id, 0);

  Mrk_2x2_aligned_TriangleInterfaceFixture tris_2x2(interface);

  EXPECT_EQ(true, manager->commit_interfaces());

  manager->get_node_global_ids();
  manager->get_node_coords();
  manager->get_predicted_node_coords();

  const size_t interface_num_nodes = tris_2x2.NumNodes;
  EXPECT_EQ(interface_num_nodes, manager->hm_node_global_ids.dimension_0());
  for (size_t i = 0; i < interface_num_nodes; ++i)
  {
    EXPECT_EQ(tris_2x2.node_gids[i], manager->hm_node_global_ids(i));
    for (size_t j = 0; j < tris_2x2.NodesPerFace; ++j) {
      EXPECT_EQ(tris_2x2.node_coords[i][j], manager->hm_node_coords(i, j));
      EXPECT_EQ(tris_2x2.node_coords[i][j], manager->hm_predicted_node_coords(i, j));
    }
  }

  manager->get_face_global_ids();
  manager->get_face_to_interface_and_side();
  manager->get_face_to_num_nodes();
  manager->get_face_to_nodes();

  const size_t interface_num_faces = tris_2x2.NumFaces;
  const size_t interface_num_nodes_per_face = tris_2x2.NodesPerFace;
  EXPECT_EQ(interface_num_faces, manager->hm_face_global_ids.dimension_0());
  for (size_t i = 0; i < interface_num_faces; ++i)
  {
    EXPECT_EQ(tris_2x2.face_gids[i], manager->hm_face_global_ids(i));
    EXPECT_EQ(interface_id, manager->hm_face_to_interface_and_side(i, 0));
    EXPECT_EQ(tris_2x2.face_sides[i], manager->hm_face_to_interface_and_side(i, 1));
    EXPECT_EQ(interface_num_nodes_per_face, manager->hm_face_to_num_nodes(i));

    for (size_t j = 0; j < tris_2x2.NodesPerFace; ++j) {
      EXPECT_EQ(tris_2x2.face_node_gids[i][j],
                manager->hm_node_global_ids(manager->hm_face_to_nodes(i, j)));
    }
  }
}


TEST(morkon,compute_normals_single_tri) {
  using namespace morkon_exp;

  const int DIM(3);

  //make an empty mesh
  Mrk_SurfaceMesh<default_kokkos_device_t,DIM> theMesh;
  Mrk_SurfaceMesh<default_kokkos_device_t,DIM>::face_to_nodes_t &faceToNodes = theMesh.m_face_to_nodes;
  //resize to store one face
  Kokkos::resize(faceToNodes,1);
  //insert node connectivity data for one triangular face
  const int faceNumber(0);
  faceToNodes(faceNumber,0) = 0;
  faceToNodes(faceNumber,1) = 1;
  faceToNodes(faceNumber,2) = 2;

  //set up node-to-face connectivity matrix
  Mrk_SurfaceMesh<default_kokkos_device_t,DIM>::node_to_faces_t &nodeToFaces = theMesh.m_node_to_faces;
  const std::string label("node to face connectivity");
  const int nNodes(3);
  const int nFaces(1);
  const int nEntries(3);
  int val[] = {faceNumber,faceNumber,faceNumber};
  int rows[] = {0,1,2,3};
  int cols[] = {0,0,0};
  nodeToFaces.import( label, nNodes, nFaces, nEntries, val, rows, cols );

  //create empty field data
  Mrk_Fields<default_kokkos_device_t,DIM> fieldData;
  //extract nodal coordinates
  Mrk_Fields<default_kokkos_device_t,DIM>::points_t &nodeCoords = fieldData.m_node_coords;
  //resize to store three nodes
  Kokkos::resize(nodeCoords,3);
  //insert nodal coordinates
  int nodeNumber(0);
  nodeCoords(nodeNumber,0) = 0.0;
  nodeCoords(nodeNumber,1) = 0.0;
  nodeCoords(nodeNumber,2) = 5.0;
  nodeNumber = 1;
  nodeCoords(nodeNumber,0) = 1.0;
  nodeCoords(nodeNumber,1) = 0.0;
  nodeCoords(nodeNumber,2) = 5.0;
  nodeNumber = 2;
  nodeCoords(nodeNumber,0) = 0.0;
  nodeCoords(nodeNumber,1) = 1.0;
  nodeCoords(nodeNumber,2) = 5.0;
  //extract face normals
  Mrk_Fields<default_kokkos_device_t,DIM>::normals_t &faceNormals = fieldData.m_face_normals;
  //resize to store one face normal
  Kokkos::resize(faceNormals,1);
  //extract node normals
  Mrk_Fields<default_kokkos_device_t,DIM>::normals_t &nodeNormals = fieldData.m_node_normals;
  //resize to store three node normals
  Kokkos::resize(nodeNormals,3);

  //compute face normals
  compute_face_normals<default_kokkos_device_t,DIM,MRK_TRI3> cfn(theMesh,fieldData);
  const int faceIndex(0);
  cfn(faceIndex);
  EXPECT_DOUBLE_EQ( 0.0,faceNormals(faceIndex,0) );
  EXPECT_DOUBLE_EQ( 0.0,faceNormals(faceIndex,1) );
  EXPECT_DOUBLE_EQ( 1.0,faceNormals(faceIndex,2) );

  //compute node normals
  compute_node_normals_from_faces<default_kokkos_device_t,DIM> cnn(theMesh,fieldData);
  for (int nodeIndex=0; nodeIndex<3; ++nodeIndex) {
    cnn(nodeIndex);
    EXPECT_DOUBLE_EQ( 0.0,nodeNormals(nodeIndex,0) );
    EXPECT_DOUBLE_EQ( 0.0,nodeNormals(nodeIndex,1) );
    EXPECT_DOUBLE_EQ( 1.0,nodeNormals(nodeIndex,2) );
  }

}//end TEST(morkon,compute_normals_single_tri) 


TEST(morkon, manager_compute_face_normals_2x2)
{
  using namespace morkon_exp;
  typedef Morkon_Manager_Tester<default_kokkos_device_t, 3, MRK_TRI3>  manager_3d_t;
  typedef Teuchos::RCP< manager_3d_t >                               manager_3d_ptr;
  typedef Interface<default_kokkos_device_t, 3, MRK_TRI3>            interface_3d_t;
  typedef Teuchos::RCP< interface_3d_t >                           interface_3d_ptr;

  manager_3d_ptr manager = manager_3d_t::MakeInstance(0, FACET_NORMAL_PROJECTION, 0);

  const int interface_id = 17;
  interface_3d_ptr interface = manager->create_interface(interface_id, 0);

  Mrk_2x2_offset_TriangleInterfaceFixture tris_2x2(interface);

  EXPECT_EQ(true, manager->commit_interfaces());

  EXPECT_EQ(true, manager->compute_normals());
  manager->get_face_normals();

  for (size_t face_i = 0; face_i < tris_2x2.NumFaces; ++face_i)
  {
    EXPECT_DOUBLE_EQ(0.0, manager->hm_face_normals(face_i, 0));
    EXPECT_DOUBLE_EQ(0.0, manager->hm_face_normals(face_i, 1));
  }

  EXPECT_DOUBLE_EQ(1.0, manager->hm_face_normals(0, 2));
  EXPECT_DOUBLE_EQ(1.0, manager->hm_face_normals(1, 2));
  EXPECT_DOUBLE_EQ(-1.0, manager->hm_face_normals(2, 2));
  EXPECT_DOUBLE_EQ(-1.0, manager->hm_face_normals(3, 2));
}

TEST(morkon, manager_find_possible_contact_face_pairs_2x2)
{
  using namespace morkon_exp;
  typedef Morkon_Manager_Tester<default_kokkos_device_t, 3, MRK_TRI3>  manager_3d_t;
  typedef Teuchos::RCP< manager_3d_t >                               manager_3d_ptr;
  typedef Interface<default_kokkos_device_t, 3, MRK_TRI3>            interface_3d_t;
  typedef Teuchos::RCP< interface_3d_t >                           interface_3d_ptr;
  typedef MorkonCommonlyUsed<default_kokkos_device_t, 3>              morkon_common;
  typedef typename morkon_common::coarse_search_results_t    coarse_search_results_t;

  manager_3d_ptr manager = manager_3d_t::MakeInstance(0, FACET_NORMAL_PROJECTION, 0);

  const int interface_id = 17;
  interface_3d_ptr interface = manager->create_interface(interface_id, 0);

  Mrk_2x2_offset_TriangleInterfaceFixture tris_2x2(interface);

  EXPECT_EQ(true, manager->commit_interfaces());
  EXPECT_EQ(true, manager->compute_normals());

  coarse_search_results_t search_results= manager->find_possible_contact_face_pairs();
  coarse_search_results_t::HostMirror search_results_host = Kokkos::create_mirror_view(search_results);
  Kokkos::deep_copy(search_results_host, search_results);

  EXPECT_EQ(4, search_results_host.dimension_0());

  EXPECT_EQ(0, search_results_host(0,0));
  EXPECT_EQ(2, search_results_host(0,1));
  EXPECT_EQ(0, search_results_host(1,0));
  EXPECT_EQ(3, search_results_host(1,1));

  EXPECT_EQ(1, search_results_host(2,0));
  EXPECT_EQ(2, search_results_host(2,1));
  EXPECT_EQ(1, search_results_host(3,0));
  EXPECT_EQ(3, search_results_host(3,1));
}

TEST(morkon, manager_compute_face_normals_2x4)
{
  using namespace morkon_exp;
  typedef Morkon_Manager_Tester<default_kokkos_device_t, 3, MRK_TRI3>  manager_3d_t;
  typedef Teuchos::RCP< manager_3d_t >                               manager_3d_ptr;
  typedef Interface<default_kokkos_device_t, 3, MRK_TRI3>            interface_3d_t;
  typedef Teuchos::RCP< interface_3d_t >                           interface_3d_ptr;

  manager_3d_ptr manager = manager_3d_t::MakeInstance(0, FACET_NORMAL_PROJECTION, 0);

  const int interface_id = 17;
  interface_3d_ptr interface = manager->create_interface(interface_id, 0);

  Mrk_2x4_offset_TriangleInterfaceFixture tris_2x4(interface);

  EXPECT_EQ(true, manager->commit_interfaces());

  EXPECT_EQ(true, manager->compute_normals());
  manager->get_face_normals();

  for (size_t face_i = 0; face_i < tris_2x4.NumFaces; ++face_i)
  {
    EXPECT_DOUBLE_EQ(0.0, manager->hm_face_normals(face_i, 0));
    EXPECT_DOUBLE_EQ(0.0, manager->hm_face_normals(face_i, 1));
  }

  EXPECT_DOUBLE_EQ(1.0, manager->hm_face_normals(0, 2));
  EXPECT_DOUBLE_EQ(1.0, manager->hm_face_normals(1, 2));
  EXPECT_DOUBLE_EQ(1.0, manager->hm_face_normals(2, 2));
  EXPECT_DOUBLE_EQ(1.0, manager->hm_face_normals(3, 2));
  EXPECT_DOUBLE_EQ(-1.0, manager->hm_face_normals(4, 2));
  EXPECT_DOUBLE_EQ(-1.0, manager->hm_face_normals(5, 2));
  EXPECT_DOUBLE_EQ(-1.0, manager->hm_face_normals(6, 2));
  EXPECT_DOUBLE_EQ(-1.0, manager->hm_face_normals(7, 2));
}

TEST(morkon, manager_find_possible_contact_face_pairs_2x4)
{
  using namespace morkon_exp;
  typedef Morkon_Manager_Tester<default_kokkos_device_t, 3, MRK_TRI3>  manager_3d_t;
  typedef Teuchos::RCP< manager_3d_t >                               manager_3d_ptr;
  typedef Interface<default_kokkos_device_t, 3, MRK_TRI3>            interface_3d_t;
  typedef Teuchos::RCP< interface_3d_t >                           interface_3d_ptr;
  typedef MorkonCommonlyUsed<default_kokkos_device_t, 3>              morkon_common;
  typedef typename morkon_common::coarse_search_results_t    coarse_search_results_t;

  manager_3d_ptr manager = manager_3d_t::MakeInstance(0, FACET_NORMAL_PROJECTION, 0);

  const int interface_id = 17;
  interface_3d_ptr interface = manager->create_interface(interface_id, 0);

  Mrk_2x4_offset_TriangleInterfaceFixture tris_2x4(interface);

  EXPECT_EQ(true, manager->commit_interfaces());
  EXPECT_EQ(true, manager->compute_normals());

  coarse_search_results_t search_results= manager->find_possible_contact_face_pairs();
  coarse_search_results_t::HostMirror search_results_host = Kokkos::create_mirror_view(search_results);
  Kokkos::deep_copy(search_results_host, search_results);

  EXPECT_EQ(4, search_results_host.dimension_0());

  EXPECT_EQ(2, search_results_host(0,0));
  EXPECT_EQ(4, search_results_host(0,1));
  EXPECT_EQ(2, search_results_host(1,0));
  EXPECT_EQ(5, search_results_host(1,1));

  EXPECT_EQ(3, search_results_host(2,0));
  EXPECT_EQ(4, search_results_host(2,1));
  EXPECT_EQ(3, search_results_host(3,0));
  EXPECT_EQ(5, search_results_host(3,1));
}

int main( int argc, char *argv[] ) {
  ::testing::InitGoogleTest(&argc,argv);
  return RUN_ALL_TESTS();
}
