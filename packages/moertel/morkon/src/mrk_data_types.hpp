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
//     - member variables for Mrk_SurfaceMesh
//     - generalize to support higher-order elements


namespace morkon_exp {

enum MorkonFaceType { MRK_LINE2=0,
                      MRK_TRI3,
                      MRK_QUAD4 };

template <unsigned int MORKON_FACE_TYPE>
struct TopoConsts { };

template <> struct TopoConsts<MRK_TRI3>
{
    enum {SPATIAL_DIM = 3u,
          NODES_PER_FACE = 3u};
};

template <> struct TopoConsts<MRK_QUAD4>
{
    enum {SPATIAL_DIM = 3u,
          NODES_PER_FACE = 4u};
};

template  <typename DeviceType, unsigned int DIM>
struct MorkonCommonlyUsed
{
  typedef typename DeviceType::execution_space                            execution_space;
  typedef Kokkos::View<local_idx_t *[2], execution_space>         coarse_search_results_t;
  typedef Kokkos::View<local_idx_t *[2], execution_space>   const_coarse_search_results_t;
  typedef Kokkos::DualView<typename coarse_search_results_t::value_type *[2],
                           typename coarse_search_results_t::array_layout,
                           typename coarse_search_results_t::execution_space>  coarse_search_results_dvt;
};


template <typename DeviceType, int DIM, int ORDER = 1>
struct FaceType
{
  unsigned m_num_nodes;
};


template  <typename DeviceType, unsigned int DIM>
struct Mrk_SurfaceMesh
{
};

template  <typename DeviceType>
struct Mrk_SurfaceMesh<DeviceType, 3>
{
  typedef typename DeviceType::execution_space                                execution_space;
  typedef Kokkos::View<global_idx_t *, execution_space>                 local_to_global_idx_t;
  typedef Kokkos::View<local_idx_t *, execution_space>                   local_to_local_idx_t;
  typedef Kokkos::CrsMatrix<local_idx_t, local_idx_t, execution_space >       node_to_faces_t;
  typedef Kokkos::View<local_idx_t *, execution_space>                    face_to_num_nodes_t;
  typedef Kokkos::View<local_idx_t *[4], execution_space>                     face_to_nodes_t;
  typedef Kokkos::View<const local_idx_t *[4], execution_space,
                              Kokkos::MemoryRandomAccess>                  face_to_nodes_mrat;

  node_to_faces_t          m_node_to_faces;
  face_to_num_nodes_t  m_face_to_num_nodes;
  face_to_nodes_t          m_face_to_nodes;

  typedef Kokkos::DualView<typename local_to_global_idx_t::value_type *,
                           typename local_to_global_idx_t::array_layout,
                           typename local_to_global_idx_t::execution_space>  local_to_global_idx_dvt;

  typedef Kokkos::DualView<typename face_to_num_nodes_t::value_type *,
                           typename face_to_num_nodes_t::array_layout,
                           typename face_to_num_nodes_t::execution_space>      face_to_num_nodes_dvt;

  typedef Kokkos::DualView<typename face_to_nodes_t::value_type *[4],
                           typename face_to_nodes_t::array_layout,
                           typename face_to_num_nodes_t::execution_space>          face_to_nodes_dvt;
};

template  <typename DeviceType, unsigned int DIM>
struct Mrk_Fields
{
  typedef typename DeviceType::execution_space execution_space;

  typedef Kokkos::View<double *[DIM], execution_space>              points_t;
  typedef Kokkos::View<const double *[DIM], execution_space,
                       Kokkos::MemoryRandomAccess>               points_mrat;
  typedef typename points_t::HostMirror                           points_hmt;
  typedef points_t                                                 normals_t;
  typedef points_mrat                                           normals_mrat;

  points_t            m_node_coords;
  normals_t          m_node_normals;
  normals_t          m_face_normals;

  points_t  m_predicted_node_coords;

  typedef Kokkos::DualView<typename points_t::value_type *[DIM],
                           typename points_t::array_layout,
                           typename points_t::execution_space>    points_dvt;
};

template  <typename DeviceType, unsigned int DIM>
struct Mrk_MortarPallets
{
  typedef typename DeviceType::execution_space                      execution_space;
  typedef Kokkos::View<local_idx_t *[2], execution_space>        generating_faces_t;
  typedef Kokkos::View<double *[DIM], execution_space>                     points_t;
  typedef points_t                                                        normals_t;
  typedef Kokkos::View<double *[DIM][DIM], execution_space>              vertices_t;

  //  Each index quadruple is (pallet #, node # on pallet, interface side, psi or eta)
  typedef Kokkos::View<double *[DIM][2][DIM - 1], execution_space> shape_fn_parms_t;

  // Each pallet is associated with a mortar-side face and a non-mortar-side face.
  generating_faces_t  m_generating_faces;

  // Each pallet has a projection (mortar) plane.
  normals_t              m_plane_normals;
  points_t             m_plane_witnesses;

  // Each pallet has DIM vertices that lie on that plane.
  vertices_t                  m_vertices;

  // Each pallet has shape function parameters for the mortar and the non-mortar side
  // points that project to the pallet vertices.
  shape_fn_parms_t     m_vertex_sf_parms;

  void resize(int newsize)
  {
    Kokkos::resize(m_generating_faces, newsize);
    Kokkos::resize(m_plane_normals, newsize);
    Kokkos::resize(m_plane_witnesses, newsize);
    Kokkos::resize(m_vertices, newsize);
    Kokkos::resize(m_vertex_sf_parms, newsize);
  }
};

}

#endif
