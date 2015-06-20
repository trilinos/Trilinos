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

#ifndef MORKON_EXP_FUNCTORS_H
#define MORKON_EXP_FUNCTORS_H

#include <mrk_data_types.hpp>

namespace morkon_exp {

template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE>
struct compute_face_normals
{

};

template <typename DeviceType>
struct compute_face_normals<DeviceType, 2,MRK_LINE2>
{
  typedef typename DeviceType::execution_space        execution_space;

  typedef Mrk_SkinOnlyMesh<DeviceType, 2>                      mesh_t;
  typedef typename mesh_t::face_connectivity_data_t     face_data_t;
  typedef typename face_data_t::face_to_nodes_t      face_to_nodes_t;

  typedef Mrk_Fields<DeviceType, 2>          fields_t;
  typedef typename fields_t::points_mrat  points_mrat;
  typedef typename fields_t::normals_t      normals_t;

  face_to_nodes_t    m_face_nodes;
  points_mrat     m_node_coords;
  normals_t       m_face_normals;

  compute_face_normals(mesh_t mesh, fields_t fields)
    : m_face_nodes(mesh.m_face_data.m_face_to_nodes)
    , m_node_coords(fields.m_node_coords)
    , m_face_normals(fields.m_face_normals)
  {
    assert(m_face_nodes.dimension_0() == m_face_normals.dimension_0());

    Kokkos::parallel_for(m_face_nodes.dimension_0(), *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (unsigned face_i) const
  {
    //this function is hard-wired for two-node segments.

    const unsigned int DIM(2);

    local_idx_t head_i, tail_i;
    tail_i = m_face_nodes(face_i, 0);
    head_i = m_face_nodes(face_i, 1);

    double tail[DIM], head[DIM], ext[DIM];
    tail[0] = m_node_coords(tail_i, 0);
    tail[1] = m_node_coords(tail_i, 1);
    head[0] = m_node_coords(head_i, 0);
    head[1] = m_node_coords(head_i, 1);

    // RHR
    ext[0] = head[1] - tail[1];
    ext[1] = tail[0] - head[0];

    // Normalize
    double len = sqrt(ext[0] * ext[0] + ext[1] * ext[1]);
    m_face_normals(face_i, 0) = ext[0] / len;
    m_face_normals(face_i, 1) = ext[1] / len;
  }
};


template <typename DeviceType>
struct compute_face_normals<DeviceType, 3,MRK_TRI3>
{
  typedef typename DeviceType::execution_space             execution_space;

  typedef Mrk_SkinOnlyMesh<DeviceType, 3>                           mesh_t;
  typedef typename mesh_t::face_connectivity_data_t          face_data_t;
  typedef typename face_data_t::face_to_num_nodes_t   face_to_num_nodes_t;
  typedef typename face_data_t::face_to_nodes_t           face_to_nodes_t;

  typedef Mrk_Fields<DeviceType, 3>          fields_t;
  typedef typename fields_t::points_mrat  points_mrat;
  typedef typename fields_t::normals_t      normals_t;

  face_to_num_nodes_t m_face_num_nodes;   // In case we want to generalize or do robust direction.
  face_to_nodes_t         m_face_nodes;
  points_mrat          m_node_coords;
  normals_t            m_face_normals;

  compute_face_normals(mesh_t mesh, fields_t fields)
    : m_face_num_nodes(mesh.m_face_data.m_face_to_num_nodes)
    , m_face_nodes(mesh.m_face_data.m_face_to_nodes)
    , m_node_coords(fields.m_node_coords)
    , m_face_normals(fields.m_face_normals)
  {
    assert(m_face_nodes.dimension_0() == m_face_normals.dimension_0());

    Kokkos::parallel_for(m_face_nodes.dimension_0(), *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (unsigned face_i) const
  {
    //this function is hard-wired for three-node triangles
    const unsigned int DIM(3); 
    const int NUM_NODES_PER_ELEMENT(3);
    
    local_idx_t node[NUM_NODES_PER_ELEMENT];
    for (unsigned j = 0; j < NUM_NODES_PER_ELEMENT; ++j)
      node[j] = m_face_nodes(face_i, j);

    //a vector has three components in three dimensions.
    double pts[NUM_NODES_PER_ELEMENT][DIM];
    for (unsigned j = 0; j < NUM_NODES_PER_ELEMENT; ++j) {
      for (unsigned k = 0; k < DIM; ++k)
        pts[j][k] = m_node_coords(node[j], k);
    }//end for for (unsigned j = 0; j < NUM_NODES_PER_ELEMENT; ++j) {

    //a vector has three components in three dimensions.
    double vecs[2][DIM];
    for (unsigned j = 0; j < 2; ++j) {
      for (unsigned k = 0; k < DIM; ++k)
        vecs[j][k] = pts[j+1][k] - pts[j][k];
    }//end for (unsigned j = 0; j < 2; ++j) {

    // WE NEED A 3D MATVEC & GEOMETRY HEADER LIBRARY THAT VECTORIZES, ETC.,

    // Cross product
    //a vector has three components in three dimensions.
    double ext[DIM];
    ext[0] = vecs[0][1] * vecs[1][2] - vecs[0][2] * vecs[1][1];
    ext[1] = vecs[0][2] * vecs[1][0] - vecs[0][0] * vecs[1][2];
    ext[2] = vecs[0][0] * vecs[1][1] - vecs[0][1] * vecs[1][0];

    // Normalize
    double len = sqrt(ext[0] * ext[0] + ext[1] * ext[1] + ext[2] * ext[2]);
    m_face_normals(face_i, 0) = ext[0] / len;
    m_face_normals(face_i, 1) = ext[1] / len;
    m_face_normals(face_i, 2) = ext[2] / len;
  }
};

template <typename DeviceType>
struct compute_face_normals<DeviceType, 3,MRK_QUAD4>
{
  typedef typename DeviceType::execution_space             execution_space;

  typedef Mrk_SkinOnlyMesh<DeviceType, 3>                           mesh_t;
  typedef typename mesh_t::face_connectivity_data_t          face_data_t;
  typedef typename face_data_t::face_to_num_nodes_t   face_to_num_nodes_t;
  typedef typename face_data_t::face_to_nodes_t           face_to_nodes_t;

  typedef Mrk_Fields<DeviceType, 3>          fields_t;
  typedef typename fields_t::points_mrat  points_mrat;
  typedef typename fields_t::normals_t      normals_t;

  face_to_num_nodes_t m_face_num_nodes;   // In case we want to generalize or do robust direction.
  face_to_nodes_t         m_face_nodes;
  points_mrat          m_node_coords;
  normals_t            m_face_normals;

  compute_face_normals(mesh_t mesh, fields_t fields)
    : m_face_num_nodes(mesh.m_face_data.m_face_to_num_nodes)
    , m_face_nodes(mesh.m_face_data.m_face_to_nodes)
    , m_node_coords(fields.m_node_coords)
    , m_face_normals(fields.m_face_normals)
  {
    assert(m_face_nodes.dimension_0() == m_face_normals.dimension_0());

    Kokkos::parallel_for(m_face_nodes.dimension_0(), *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (unsigned face_i) const
  {
    //this function is hard-wired for four-node quads
    const unsigned int DIM(3); 
    const int NUM_NODES_PER_ELEMENT(4);
    
    local_idx_t node[NUM_NODES_PER_ELEMENT];
    for (unsigned j = 0; j < NUM_NODES_PER_ELEMENT; ++j)
      node[j] = m_face_nodes(face_i, j);

    //a vector has three components in three dimensions.
    double pts[NUM_NODES_PER_ELEMENT][DIM];
    for (unsigned j = 0; j < NUM_NODES_PER_ELEMENT; ++j) {
      for (unsigned k = 0; k < DIM; ++k)
        pts[j][k] = m_node_coords(node[j], k);
    }//end for for (unsigned j = 0; j < NUM_NODES_PER_ELEMENT; ++j) {

    //a vector has three components in three dimensions.
    double vecs[2][DIM];
    for (unsigned k = 0; k < DIM; ++k) {
      vecs[0][k] = pts[2][k] - pts[0][k];
      vecs[1][k] = pts[3][k] - pts[1][k];
    }

    // WE NEED A 3D MATVEC & GEOMETRY HEADER LIBRARY THAT VECTORIZES, ETC.,

    // Cross product
    //a vector has three components in three dimensions.
    double ext[DIM];
    ext[0] = vecs[0][1] * vecs[1][2] - vecs[0][2] * vecs[1][1];
    ext[1] = vecs[0][2] * vecs[1][0] - vecs[0][0] * vecs[1][2];
    ext[2] = vecs[0][0] * vecs[1][1] - vecs[0][1] * vecs[1][0];

    // Normalize
    double len = sqrt(ext[0] * ext[0] + ext[1] * ext[1] + ext[2] * ext[2]);
    m_face_normals(face_i, 0) = ext[0] / len;
    m_face_normals(face_i, 1) = ext[1] / len;
    m_face_normals(face_i, 2) = ext[2] / len;
  }
};


template <typename DeviceType, unsigned int DIM >
struct compute_node_normals_from_faces
{
  typedef typename DeviceType::execution_space             execution_space;

  typedef Mrk_SkinOnlyMesh<DeviceType, DIM>                              mesh_t;
  typedef typename mesh_t::node_connectivity_data_t         node_connectivity_t;
  typedef typename node_connectivity_t::node_to_faces_t       node_to_faces_t;

  typedef Mrk_Fields<DeviceType, 3>            fields_t;
  typedef typename fields_t::normals_mrat  normals_mrat;
  typedef typename fields_t::normals_t        normals_t;

  node_to_faces_t     m_node_to_faces;
  normals_mrat    m_face_normals;
  normals_t          m_node_normals;

  compute_node_normals_from_faces(mesh_t skin_mesh, fields_t fields)
    : m_node_to_faces(skin_mesh.m_node_data.m_node_to_faces)
    , m_face_normals(fields.m_face_normals)
    , m_node_normals(fields.m_node_normals)
  {
    const int num_nodes = m_node_to_faces.numRows();
    const int num_node_normals =  m_node_normals.dimension_0();
    assert(num_nodes==num_node_normals); 
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (unsigned node_i) const
  {
    const int i_faces_begin = m_node_to_faces.graph.row_map(node_i);
    const int i_faces_end   = m_node_to_faces.graph.row_map(node_i + 1);
    const int num_faces = i_faces_end - i_faces_begin;

    double nml[DIM];
    for (int i = 0; i < DIM; ++i) 
      nml[i] = 0;

    // Sum the face normals
    for (int j = i_faces_begin; j < i_faces_end; ++j) {
      local_idx_t face_j = m_node_to_faces.graph.entries(j);
      for (int k = 0; k < DIM; ++k)
        nml[k] += m_face_normals(face_j, k);
    }//end for (int j = i_faces_begin; j < i_faces_end; ++j)

    // Average.
    for (int k=0; k < DIM; ++k)
      m_node_normals(node_i, k) = nml[k] /  num_faces;
  }

};

} // namespace morton_exp

#endif
