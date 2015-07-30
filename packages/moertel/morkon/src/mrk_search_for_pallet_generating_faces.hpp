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

#ifndef MORKON_EXP_SEARCH_FOR_PALLET_GENERATING_FACES_H
#define MORKON_EXP_SEARCH_FOR_PALLET_GENERATING_FACES_H

#include <mrk_data_types.hpp>

namespace morkon_exp {

template <typename DeviceType, typename ScalarType>
struct AxisAlignedBB
{
    enum AABBIndices { X_MIN = 0, Y_MIN, Z_MIN, X_MAX, Y_MAX, Z_MAX, NUM_AABB_INDICES};

    typedef ScalarType                                        scalar_type;
    typedef typename DeviceType::execution_space          execution_space;
    typedef Kokkos::View<ScalarType *[NUM_AABB_INDICES]>          boxes_t;
};

template <typename DeviceType, unsigned int DIM>
struct search_for_pallet_generating_faces
{
    typedef typename DeviceType::execution_space                          execution_space;
    typedef Mrk_SurfaceMesh<DeviceType, DIM>                               surface_mesh_t;
    typedef typename surface_mesh_t::face_to_num_nodes_t              face_to_num_nodes_t;
    typedef typename surface_mesh_t::face_to_nodes_t                      face_to_nodes_t;
    typedef Mrk_Fields<DeviceType, DIM>                                          fields_t;
    typedef typename fields_t::points_t                                          points_t;
    typedef typename fields_t::points_mrat                                    points_mrat;
    typedef Kokkos::View<local_idx_t *[2], execution_space>  face_to_interface_and_side_t;
    typedef Kokkos::View<local_idx_t *[2], execution_space>      contact_search_results_t;
    typedef Kokkos::View<local_idx_t *, execution_space>                      faces_ids_t;

    typedef typename AxisAlignedBB<DeviceType, float>::boxes_t           bounding_boxes_t;

    face_to_num_nodes_t                     m_face_to_num_nodes;
    face_to_nodes_t                             m_face_to_nodes;
    points_mrat                                   m_node_coords;
    points_mrat                         m_predicted_node_coords;
    face_to_interface_and_side_t   m_face_to_interface_and_side;
    contact_search_results_t                   m_search_results;

    const float                                 m_boxes_epsilon;

    search_for_pallet_generating_faces(surface_mesh_t surface_mesh,
                                       points_t node_coords,
                                       points_t predicted_node_coords,
                                       face_to_interface_and_side_t face_to_interface_and_side,
                                       double epsilon,
                                       contact_search_results_t search_results)
        : m_face_to_num_nodes(surface_mesh.m_face_to_num_nodes)
        , m_face_to_nodes(surface_mesh.m_face_to_nodes)
        , m_node_coords(node_coords)
        , m_predicted_node_coords(predicted_node_coords)
        , m_face_to_interface_and_side(face_to_interface_and_side)
        , m_search_results(search_results)
        , m_boxes_epsilon(static_cast<float>(epsilon))
    {
        faces_ids_t  non_mortarside_faces;
        faces_ids_t      mortarside_faces;
        // Divide up the faces into non-mortarside and mortarside sets.

        bounding_boxes_t non_mortarside_boxes;
        bounding_boxes_t     mortarside_boxes;
        // Construct bounding boxes for the sets of spaces, using both the node_coords
        // and predicted node_coords and also the m_boxes_epsilon.

        // Do a (brute-force, for now) search for pairs with bounding boxes that overlap
        // Extra points if you disallow face pairs whose normals are incompatible with
        // contact.
    }

};


}

#endif
