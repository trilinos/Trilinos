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
    typedef typename DeviceType::memory_space::size_type                        size_type;
    typedef Mrk_SurfaceMesh<DeviceType, DIM>                               surface_mesh_t;
    typedef typename surface_mesh_t::face_to_num_nodes_t              face_to_num_nodes_t;
    typedef typename surface_mesh_t::face_to_nodes_t                      face_to_nodes_t;
    typedef Mrk_Fields<DeviceType, DIM>                                          fields_t;
    typedef typename fields_t::points_t                                          points_t;
    typedef typename fields_t::points_mrat                                    points_mrat;
    typedef Kokkos::View<local_idx_t *[2], execution_space>  face_to_interface_and_side_t;
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
    face_to_interface_and_side_t   m_face_to_interface_and_side;
    contact_search_results_t                   m_search_results;

    const float                                 m_boxes_epsilon;

    struct count_mortarside_faces
    {
      typedef Kokkos::DualView<int, execution_space>                        total_found_dvt;
      typedef int                                                                value_type;

      ints_vec_t                m_offsets;
      const size_type        m_last_index;
      int_view_t            m_total_found;
      face_to_interface_and_side_t   m_face_to_interface_and_side;

      count_mortarside_faces(search_for_pallet_generating_faces &parent,
                             ints_vec_t offsets_out,
                             int_view_t total_found)
        : m_offsets(offsets_out)
        , m_last_index(parent.m_face_to_interface_and_side.dimension_0() == 0
                       ? 0 : parent.m_face_to_interface_and_side.dimension_0() - 1)
        , m_total_found(total_found)
        , m_face_to_interface_and_side(parent.m_face_to_interface_and_side)
      {
        size_type output_size = m_face_to_interface_and_side.dimension_0();
        std::cout << "In count_mortarside_faces(.) " << output_size << " faces in" << std::endl;

        if (output_size > 0)
        {
            Kokkos::resize(m_offsets, output_size);
            Kokkos::parallel_scan(m_face_to_interface_and_side.dimension_0(), *this);
        }
        std::cout << "count_mortarside_faces found " << /* m_total_found() <<  */ " faces" << std::endl;
      }

      KOKKOS_INLINE_FUNCTION
      void operator() (const size_type& idx, int& offset, const bool& final_pass)  const {
        if(final_pass) {
          m_offsets(idx) = offset;
        }
        offset += m_face_to_interface_and_side(idx, 1);
        if (final_pass && idx == m_last_index) {
          m_total_found = offset;
        }
      }
    };

    struct separate_into_sides
    {
      face_ids_on_a_side_t m_non_mortarside_faces, m_mortarside_faces;
      int_view_t        m_total_faces;
      int_dualview_t m_num_mortarside;

      separate_into_sides(search_for_pallet_generating_faces &parent,
                          face_ids_on_a_side_t non_mortarside_faces,
                          face_ids_on_a_side_t mortarside_faces)
        : m_non_mortarside_faces(non_mortarside_faces)
        , m_mortarside_faces(mortarside_faces)
        , m_total_faces()
        , m_num_mortarside()
      {
        std::cout << "In separate_into_sides(..)" << std::endl;

        ints_vec_t mortarside_offsets;
        m_num_mortarside. template modify<typename int_dualview_t::t_dev>();
        count_mortarside_faces(parent, mortarside_offsets, m_num_mortarside.d_view);

        std::cout << "Back from counting mortarside faces " << std::endl;

        m_num_mortarside. template sync<typename int_dualview_t::t_host>();
        int num_mortarside_faces = m_num_mortarside.h_view();
        int total_faces = parent.m_face_to_interface_and_side.dimension_0();
        Kokkos::deep_copy(m_total_faces, total_faces);
        Kokkos::resize(m_non_mortarside_faces, total_faces - num_mortarside_faces);
        Kokkos::resize(m_mortarside_faces, num_mortarside_faces);
      }

      KOKKOS_INLINE_FUNCTION
      void operator() (const size_type& idx) {
        int mortarside_offset = m_mortarside_offsets(idx);
        int interface_id = m_face_to_interface_and_side(idx, 0);
        if (m_face_to_interface_and_side(idx, 1) == InterfaceBase::NON_MORTAR_SIDE)
        {
          int non_mortarside_offset = m_total_faces - mortarside_offset;
          m_non_mortarside_faces(non_mortarside_offset, 0) = idx;
          m_non_mortarside_faces(non_mortarside_offset, 1) = interface_id;
        }
        else
        {
          m_mortarside_faces(mortarside_offset, 0) = idx;
          m_mortarside_faces(mortarside_offset, 1) = interface_id;
        }
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
        face_ids_on_a_side_t  non_mortarside_faces;
        face_ids_on_a_side_t      mortarside_faces;
        separate_into_sides(*this, non_mortarside_faces, mortarside_faces);

        bounding_boxes_t non_mortarside_boxes;
        bounding_boxes_t     mortarside_boxes;
        // Construct bounding boxes for the sets of spaces, using both the node_coords
        // and predicted node_coords and also the m_boxes_epsilon.
        //
        // Resize the these two views.  Then use two parallel_fors to fill.

        // Do a (brute-force, for now) search for pairs with bounding boxes that overlap
        // Extra points if you disallow face pairs whose normals are incompatible with
        // contact.
        //
        // Parallel_for to count.  Parallel_scan to compute offsets.  Parallel_for to fill.

    }

};


}

#endif
