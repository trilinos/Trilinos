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

#ifndef MORKON_EXP_DATA_TYPES_H
#define MORKON_EXP_DATA_TYPES_H

#include <Kokkos_View.hpp>
#include <Kokkos_DualView.hpp>
#include <Kokkos_CrsMatrix.hpp>

#include <mrk_int_types.hpp>

// TO DO:
//     - make downward connectivity data Kokkos::DualViews for now
//     - member variables for Mrk_SkinOnlyMesh
//     - generalize to support higher-order elements



namespace morkon_exp {

enum MorkonFaceType { MRK_LINE2=0,
		      MRK_TRI3,
		      MRK_QUAD4 };

template <typename DeviceType, int DIM, int ORDER = 1>
struct FaceType
{
  unsigned m_num_nodes;
};


template  <typename DeviceType, int DIM>
struct NodeConnectivityData
{
};

template  <typename DeviceType>
struct NodeConnectivityData<DeviceType, 2>
{
  typedef typename DeviceType::execution_space execution_space;

  typedef Kokkos::CrsMatrix<local_idx_t, local_idx_t, execution_space>  Node2EdgesGraphType; // face id, node ordinal
  typedef Node2EdgesGraphType                              node_to_faces_t;

  node_to_faces_t m_node_to_faces;
};

template  <typename DeviceType>
struct NodeConnectivityData<DeviceType, 3>
{
  typedef typename DeviceType::execution_space execution_space;

  typedef Kokkos::CrsMatrix<local_idx_t, local_idx_t, execution_space >  Node2EdgesGraphType; // edge id, node ordinal
  typedef Kokkos::CrsMatrix<local_idx_t, local_idx_t, execution_space >  Node2FacesGraphType; // face id, node ordinal
  typedef Node2FacesGraphType                               node_to_faces_t;

  node_to_faces_t  m_node_to_faces;
  Node2FacesGraphType    m_node_to_edges;  // Will this be needed?
};


template <typename DeviceType, int DIM>
struct EdgeConnectivityData
{
};

template <typename DeviceType>
struct EdgeConnectivityData<DeviceType, 2>
{
  enum NodeInds { TAIL = 0, HEAD = 1 };
  enum FaceInds { FACE_ID = 0, EDGE_ORDINAL = 1 };

  typedef typename DeviceType::execution_space execution_space;

  typedef Kokkos::View<local_idx_t *[2], execution_space >  Edge2NodesDataType;  // head, tail
  typedef Edge2NodesDataType                       face_to_nodes_t;

  face_to_nodes_t m_face_to_nodes;
};

// For now, bi-linear edges in 2D.
template <typename DeviceType>
struct EdgeConnectivityData<DeviceType, 3>
{
  enum NodeInds            { TAIL = 0, HEAD = 1 };
  enum FaceInds            { FACE_ID = 0, EDGE_ORDINAL = 1 };
  enum EdgeConnectPolarity { POSITIVE_POLARITY = 0, NEGATIVE_POLARITY = 1 };

  typedef typename DeviceType::execution_space execution_space;

  typedef Kokkos::View<local_idx_t *[2], execution_space >          Edge2NodesDataType;  // head, tail
  typedef Kokkos::View<local_idx_t *[2], execution_space >  Edge2FacesConnectivityType;  // ordinal, polarity

  typedef Edge2FacesConnectivityType  edge_2_faces_connectivity_t;
  typedef Edge2NodesDataType                      face_to_nodes_t;

  face_to_nodes_t m_face_to_nodes;
};


// For now, bi-linear triangular or quad faces in 3D.
template  <typename DeviceType>
struct FaceConnectivityData
{
  typedef typename DeviceType::execution_space execution_space;

  typedef Kokkos::View<local_idx_t *, execution_space>              FaceNumSidesDataType;
  typedef Kokkos::View<local_idx_t *[4], execution_space>             Face2NodesDataType;
  typedef Kokkos::View<local_idx_t *[4], execution_space>             Face2EdgesDataType;
  typedef Kokkos::View<polarity_t *, DeviceType>              Face2EdgesPolartiyDataType;

  typedef FaceNumSidesDataType                           face_to_num_nodes_t;
  typedef Face2NodesDataType                                 face_to_nodes_t;

  FaceNumSidesDataType m_face_to_num_nodes;
  face_to_nodes_t       m_face_to_nodes;
};


template  <typename DeviceType, unsigned int DIM>
struct Mrk_SkinOnlyMesh
{
};

template  <typename DeviceType>
struct Mrk_SkinOnlyMesh<DeviceType, 2>
{
  typedef typename DeviceType::execution_space execution_space;

  // Will need one of these each for nodes and edges, respectively.
  typedef Kokkos::View<global_idx_t *, execution_space> local_to_global_idx_t;

  typedef NodeConnectivityData<DeviceType, 2>      node_connectivity_data_t;
  typedef EdgeConnectivityData<DeviceType, 2>   face_connectivity_data_t;

  node_connectivity_data_t        m_node_data;
  face_connectivity_data_t  m_face_data;
};

template  <typename DeviceType>
struct Mrk_SkinOnlyMesh<DeviceType, 3>
{
  typedef typename DeviceType::execution_space execution_space;

  typedef Kokkos::View<global_idx_t *, execution_space> local_to_global_idx_t;

  typedef NodeConnectivityData<DeviceType, 3>     node_connectivity_data_t;
  typedef EdgeConnectivityData<DeviceType, 3>     edge_connectivity_data_t;
  typedef FaceConnectivityData<DeviceType>     face_connectivity_data_t;

  node_connectivity_data_t       m_node_data;
  edge_connectivity_data_t       m_edge_data;
  face_connectivity_data_t m_face_data;

};

template  <typename DeviceType, unsigned int DIM>
struct Mrk_Fields
{
  typedef typename DeviceType::execution_space execution_space;

  typedef Kokkos::View<double *[DIM], execution_space>                                 points_t;
  typedef Kokkos::View<double *[DIM], execution_space, Kokkos::MemoryRandomAccess>  points_mrat;

  typedef points_t       normals_t;
  typedef points_mrat normals_mrat;

  points_t       m_node_coords;
  normals_t     m_node_normals;
  normals_t  m_face_normals;
};

template  <typename DeviceType, unsigned int DIM>
struct Mrk_MortarPallets
{
  typedef typename DeviceType::execution_space execution_space;

  typedef Kokkos::View<local_idx_t *[2], execution_space>       faces_t;

  typedef Kokkos::View<double *[DIM], execution_space>            points_t;
  typedef points_t                                               normals_t;

  typedef Kokkos::View<double *[DIM][DIM], execution_space>     vertices_t;
  typedef Kokkos::View<double *[2][DIM][DIM - 1], execution_space> shape_fn_parms_t;

  // Each pallet is associated with a mortar-side face and a non-mortar-side face.
  faces_t  m_generating_faces;

  // Each pallet has a projection (mortar) plane.
  normals_t         m_plane_normals;
  points_t        m_plane_witnesses;

  // Each pallet has DIM vertices that lie on that plane.
  vertices_t             m_vertices;

  // Each pallet has shape function parameters for the mortar and the non-mortar side
  // points that project to the pallet vertices.
  // Are they needed?  I.e., are they helpful for computing the integration points
  // of each pallet on the faces?
  shape_fn_parms_t  m_vertex_sf_parms;
};

}

#endif
