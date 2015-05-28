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

TEST(morkon,just_check_if_it_compiles) {
  using namespace morkon_exp;
  typedef Morkon_Manager<default_kokkos_device_t, 3, MRK_QUAD4>    default_manager_3d_t;
  typedef Teuchos::RCP< default_manager_3d_t >        default_manager_3d_ptr;
  typedef Interface<default_kokkos_device_t, 3, MRK_QUAD4>     default_interface_3d_t;
  typedef Teuchos::RCP< default_interface_3d_t >  default_interface_3d_ptr;
  default_manager_3d_ptr manager_0 = default_manager_3d_t::MakeInstance(0, 0);
  default_interface_3d_ptr interface_0 = manager_0->create_interface(0,0);
  manager_0->commit_interfaces();
  Tpetra::CrsMatrix<> *dummy_D = 0;
  Tpetra::CrsMatrix<> *dummy_M = 0;
  manager_0->mortar_integrate(dummy_D, dummy_M);
}


TEST(morkon,compute_normals_single_tri) {
  using namespace morkon_exp;

  const int DIM(3);

  //make an empty mesh
  Mrk_SkinOnlyMesh<default_kokkos_device_t,DIM> theMesh;
  //extract face connectivity data
  Mrk_SkinOnlyMesh<default_kokkos_device_t,DIM>::face_connectivity_data_t &faceData = theMesh.m_face_data;
  //extract face to node data
  FaceConnectivityData<default_kokkos_device_t>::face_to_nodes_t &faceToNodes = faceData.m_face_to_nodes;
  //resize to store one face
  Kokkos::resize(faceToNodes,1);
  //insert node connectivity data for one triangular face
  const int faceNumber(0);
  faceToNodes(faceNumber,0) = 0;
  faceToNodes(faceNumber,1) = 1;
  faceToNodes(faceNumber,2) = 2;

  //set up node-to-face connectivity matrix
  Mrk_SkinOnlyMesh<default_kokkos_device_t,DIM>::node_connectivity_data_t &nodeData = theMesh.m_node_data;
  NodeConnectivityData<default_kokkos_device_t,DIM>::node_to_faces_t &nodeToFaces = nodeData.m_node_to_faces;
  const std::string label("node to face connectivity");
  const int nNodes(3);
  const int nFaces(1);
  const int nEntries(3);
  int val[] = {faceNumber,faceNumber,faceNumber};
  int rows[] = {0,1,2};
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

int main( int argc, char *argv[] ) {
  ::testing::InitGoogleTest(&argc,argv);
  return RUN_ALL_TESTS();
}
