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

#include <limits>
#include <mrk_data_types.hpp>

namespace morkon_exp {

template <typename DeviceType, typename ScalarType>
struct AxisAlignedBB
{
    enum AABBIndices { X_MIN = 0, Y_MIN, Z_MIN, X_MAX, Y_MAX, Z_MAX, NUM_AABB_INDICES};

    typedef ScalarType                                        scalar_type;
    typedef typename DeviceType::execution_space          execution_space;
    typedef Kokkos::View<ScalarType *[NUM_AABB_INDICES]>          boxes_t;

    constexpr static const scalar_type min_scalar = std::numeric_limits<scalar_type>::min();
    constexpr static const scalar_type max_scalar = std::numeric_limits<scalar_type>::max();

};

template <typename DeviceType, unsigned int DIM>
struct search_for_pallet_generating_faces
{
    typedef typename DeviceType::execution_space                          execution_space;
    typedef typename DeviceType::memory_space::size_type                        size_type;
    typedef Mrk_SurfaceMesh<DeviceType, DIM>                               surface_mesh_t;
    typedef typename surface_mesh_t::face_to_num_nodes_t              face_to_num_nodes_t;
    typedef typename surface_mesh_t::face_to_nodes_t                      face_to_nodes_t;
    typedef Mrk_Fields<DeviceType, DIM>                                          fields_t;
    typedef typename fields_t::points_t                                          points_t;
    typedef typename fields_t::points_mrat                                    points_mrat;
    typedef Kokkos::View<local_idx_t *[2], execution_space>  face_to_interface_and_side_t;
    typedef Kokkos::View<const local_idx_t *[2], execution_space>  face_to_interface_and_side_ct;
    typedef Kokkos::View<local_idx_t *[2], execution_space>          face_ids_on_a_side_t;
    typedef Kokkos::View<int *, execution_space>                               ints_vec_t;
    typedef Kokkos::View<int, execution_space>                                 int_view_t;
    typedef Kokkos::DualView<int, execution_space>                         int_dualview_t;

    typedef MorkonCommonlyUsed<DeviceType, DIM>                           morkon_common_t;
    typedef typename morkon_common_t::contact_search_results_t   contact_search_results_t;

    typedef typename AxisAlignedBB<DeviceType, float>::boxes_t           bounding_boxes_t;

    face_to_num_nodes_t                     m_face_to_num_nodes;
    face_to_nodes_t                             m_face_to_nodes;
    points_mrat                                   m_node_coords;
    points_mrat                         m_predicted_node_coords;
    face_to_interface_and_side_ct  m_face_to_interface_and_side;

    face_ids_on_a_side_t                 m_non_mortarside_faces;
    face_ids_on_a_side_t                     m_mortarside_faces;

    contact_search_results_t                   m_search_results;

    const float                                 m_boxes_epsilon;

    struct count_mortarside_faces
    {
      typedef int                                                                value_type;

      ints_vec_t                m_offsets;
      const size_type        m_last_index;
      int_view_t            m_total_found;
      face_to_interface_and_side_ct  m_face_to_interface_and_side;

      count_mortarside_faces(search_for_pallet_generating_faces &parent,
                             ints_vec_t offsets_out,
                             int_view_t total_found)
        : m_offsets(offsets_out)
        , m_last_index(m_offsets.dimension_0() == 0 ? 0 : m_offsets.dimension_0() - 1)
        , m_total_found(total_found)
        , m_face_to_interface_and_side(parent.m_face_to_interface_and_side)
      {
        assert(m_offsets.dimension_0() == parent.m_face_to_interface_and_side.dimension_0());
        size_type output_size = m_face_to_interface_and_side.dimension_0();
        if (output_size > 0)
          Kokkos::parallel_scan(m_face_to_interface_and_side.dimension_0(), *this);
      }

      KOKKOS_INLINE_FUNCTION
      void operator() (const size_type& face_id, value_type& offset, const bool& final_pass)  const
      {
        if(final_pass) {
          m_offsets(face_id) = offset;
        }
        offset += m_face_to_interface_and_side(face_id, 1);
        if (final_pass && face_id == m_last_index) {
          m_total_found = offset;
        }
      }
    };

    struct separate_into_sides
    {
      face_to_interface_and_side_ct  m_face_to_interface_and_side;
      size_type                                     m_total_faces;
      ints_vec_t                             m_mortarside_offsets;
      int_dualview_t                             m_num_mortarside;
      face_ids_on_a_side_t                 m_non_mortarside_faces;
      face_ids_on_a_side_t                     m_mortarside_faces;

      separate_into_sides(search_for_pallet_generating_faces &parent)
        : m_face_to_interface_and_side(parent.m_face_to_interface_and_side)
        , m_total_faces(m_face_to_interface_and_side.dimension_0())
        , m_num_mortarside("num_mortarside")
      {
        Kokkos::resize(m_mortarside_offsets, m_total_faces);
        m_num_mortarside. template modify<typename int_dualview_t::t_dev>();
        count_mortarside_faces(parent, m_mortarside_offsets, m_num_mortarside.d_view);
        m_num_mortarside. template sync<typename int_dualview_t::t_host>();
        int num_mortarside_faces = m_num_mortarside.h_view();

        Kokkos::resize(m_non_mortarside_faces, m_total_faces - num_mortarside_faces);
        Kokkos::resize(m_mortarside_faces, num_mortarside_faces);
        Kokkos::parallel_for(m_total_faces, *this);

        parent.m_non_mortarside_faces = m_non_mortarside_faces;
        parent.m_mortarside_faces = m_mortarside_faces;
      }

      KOKKOS_INLINE_FUNCTION
      void operator() (const size_type& face_id) const
      {
        int mortarside_offset = m_mortarside_offsets(face_id);
        int interface_id = m_face_to_interface_and_side(face_id, 0);
        if (m_face_to_interface_and_side(face_id, 1) == InterfaceBase::NON_MORTAR_SIDE)
        {
          int non_mortarside_offset = face_id - mortarside_offset;
          m_non_mortarside_faces(non_mortarside_offset, 0) = face_id;
          m_non_mortarside_faces(non_mortarside_offset, 1) = interface_id;
        }
        else
        {
          m_mortarside_faces(mortarside_offset, 0) = face_id;
          m_mortarside_faces(mortarside_offset, 1) = interface_id;
        }
      }
    };

    struct construct_bounding_boxes
    {
      face_ids_on_a_side_t       m_face_interface_pairs;
      KOKKOS_INLINE_FUNCTION
      void operator() (const unsigned &idx) const
      {
        typename bounding_boxes_t::scalar_type min_corner[3] =  { bounding_boxes_t::max_scalar,
                                                                  bounding_boxes_t::max_scalar,
                                                                  bounding_boxes_t::max_scalar};

        typename bounding_boxes_t::scalar_type max_corner[3] =  { bounding_boxes_t::min_scalar,
                                                                  bounding_boxes_t::min_scalar,
                                                                  bounding_boxes_t::min_scalar};

        // YOU ARE HERE
      }
    };


    struct filter_to_sides_tag {};
    struct construct__bounding_boxes_tag {};

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
        std::cout << "In search_for_pallet_generating_faces(..)" << std::endl;

        separate_into_sides(*this);

        bounding_boxes_t non_mortarside_aabbs, mortarside_aabbs;
        Kokkos::resize(non_mortarside_aabbs, m_non_mortarside_faces.dimension_0());
        Kokkos::resize(mortarside_aabbs, m_mortarside_faces.dimension_0());


        // construct_bounding_boxes(*this, non_mortarside_faces, non_mortarside_aabbs);
        // construct_bounding_boxes(*this, mortarside_faces, mortarside_aabbs);

        // Do a (brute-force, for now) search for pairs with bounding boxes that overlap
        // Extra points if you disallow face pairs whose normals are incompatible with
        // contact.
        //
        // Parallel_for to count.  Parallel_scan to compute offsets.  Parallel_for to fill.

    }

};


}

#endif
