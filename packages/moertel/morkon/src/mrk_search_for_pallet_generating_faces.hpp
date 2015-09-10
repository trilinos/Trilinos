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
    enum AABBIndices { X_MIN = 0, X_MAX, Y_MIN, Y_MAX, Z_MIN, Z_MAX, NUM_AABB_INDICES};

    KOKKOS_INLINE_FUNCTION
    static local_idx_t min_val_idx_for_dim(local_idx_t dim) { return 2*dim; }

    KOKKOS_INLINE_FUNCTION
    static local_idx_t max_val_idx_for_dim(local_idx_t dim) { return 2*dim + 1; }

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
    typedef typename surface_mesh_t::face_to_nodes_mrat                face_to_nodes_mrat;
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
    typedef typename morkon_common_t::coarse_search_results_t     coarse_search_results_t;

    typedef float                                                        aabb_scalar_type;
    typedef AxisAlignedBB<DeviceType, aabb_scalar_type>                            aabb_t;
    typedef typename aabb_t::boxes_t                                     bounding_boxes_t;

    face_to_num_nodes_t                     m_face_to_num_nodes;
    face_to_nodes_t                             m_face_to_nodes;
    points_mrat                                   m_node_coords;
    points_mrat                         m_predicted_node_coords;
    face_to_interface_and_side_ct  m_face_to_interface_and_side;

    face_ids_on_a_side_t                 m_non_mortarside_faces;
    face_ids_on_a_side_t                     m_mortarside_faces;

    coarse_search_results_t                    m_search_results;

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
        const size_type output_size = m_face_to_interface_and_side.dimension_0();
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
          const int non_mortarside_offset = face_id - mortarside_offset;
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
      face_ids_on_a_side_t     m_face_interface_pairs;
      face_to_num_nodes_t         m_face_to_num_nodes;
      face_to_nodes_mrat              m_face_to_nodes;
      points_mrat                       m_node_coords;
      points_mrat             m_predicted_node_coords;
      const aabb_scalar_type                m_epsilon;
      bounding_boxes_t                    m_aabbs_out;

      construct_bounding_boxes(const search_for_pallet_generating_faces &parent,
                               face_ids_on_a_side_t face_interface_pairs,
                               bounding_boxes_t aabbs_out)
      : m_face_interface_pairs(face_interface_pairs),
        m_face_to_num_nodes(parent.m_face_to_num_nodes),
        m_face_to_nodes(parent.m_face_to_nodes),
        m_node_coords(parent.m_node_coords),
        m_predicted_node_coords(parent.m_predicted_node_coords),
        m_epsilon(parent.m_boxes_epsilon),
        m_aabbs_out(aabbs_out)
      {
        Kokkos::parallel_for(m_face_interface_pairs.dimension_0(), *this);
      }

      KOKKOS_INLINE_FUNCTION
      void operator() (const unsigned &idx) const
      {
        aabb_scalar_type min_corner[3] =  { aabb_t::max_scalar, aabb_t::max_scalar, aabb_t::max_scalar };
        aabb_scalar_type max_corner[3] =  { aabb_t::min_scalar, aabb_t::min_scalar, aabb_t::min_scalar };

        const local_idx_t face_id   = m_face_interface_pairs(idx, 0);
        const local_idx_t num_nodes = m_face_to_num_nodes(idx);

        for (local_idx_t node_i = 0; node_i < num_nodes; ++node_i)
        {
          const local_idx_t node_id = m_face_to_nodes(face_id, node_i);
          for (unsigned dim_j = 0; dim_j < DIM; ++dim_j)
          {
            const aabb_scalar_type val_j = m_node_coords(node_id, dim_j);
            const aabb_scalar_type pred_val_j = m_node_coords(node_id, dim_j);
            const aabb_scalar_type lb_val_j = (val_j < pred_val_j ? val_j : pred_val_j);
            const aabb_scalar_type ub_val_j = (val_j > pred_val_j ? val_j : pred_val_j);
            min_corner[dim_j] = (lb_val_j < min_corner[dim_j] ? lb_val_j : min_corner[dim_j]);
            max_corner[dim_j] = (ub_val_j > max_corner[dim_j] ? ub_val_j : max_corner[dim_j]);
          }
        }

        for (unsigned dim_j= 0; dim_j < DIM; ++dim_j)
        {
          m_aabbs_out(idx, aabb_t::min_val_idx_for_dim(dim_j)) = min_corner[dim_j] - m_epsilon;
          m_aabbs_out(idx, aabb_t::max_val_idx_for_dim(dim_j)) = max_corner[dim_j] + m_epsilon;
        }
      }
    };

    struct brute_force_search
    {
      typedef typename ints_vec_t::value_type     value_type;

      enum IntersectsFunctionMode {COUNT, FILL};

      const int m_num_nonmortar_faces;
      const int m_last_scan_index;
      face_ids_on_a_side_t  m_non_mortarside_face_interface_pairs;
      face_ids_on_a_side_t      m_mortarside_face_interface_pairs;
      bounding_boxes_t                     m_non_mortarside_aabbs;
      bounding_boxes_t                         m_mortarside_aabbs;
      coarse_search_results_t                    m_search_results;

      ints_vec_t                                         m_counts;
      ints_vec_t                                        m_offsets;
      int_dualview_t                                m_total_found;

      KOKKOS_INLINE_FUNCTION
      void find_intersections (const int &idx, const IntersectsFunctionMode mode) const
      {
        const local_idx_t non_mortar_face_idx = idx;
        const local_idx_t fill_offset = (mode == FILL? m_offsets[idx] : 0);
        const local_idx_t face_A = m_non_mortarside_face_interface_pairs(non_mortar_face_idx, 0);

        local_idx_t aabb_min_idx[DIM], aabb_max_idx[DIM];
        aabb_scalar_type min_corner_A[DIM], max_corner_A[DIM];
        for (unsigned i = 0; i < DIM; ++i)
        {
          aabb_min_idx[i] = aabb_t::min_val_idx_for_dim(i);
          aabb_max_idx[i] = aabb_t::max_val_idx_for_dim(i);

          min_corner_A[i] = m_non_mortarside_aabbs(non_mortar_face_idx, aabb_min_idx[i]);
          max_corner_A[i] = m_non_mortarside_aabbs(non_mortar_face_idx, aabb_max_idx[i]);
        }
        const local_idx_t interface_A = m_non_mortarside_face_interface_pairs(non_mortar_face_idx, 1);

        int count = 0;
        for (local_idx_t mortar_face_idx = 0; mortar_face_idx < m_num_nonmortar_faces; ++mortar_face_idx)
        {
          const local_idx_t interface_B = m_mortarside_face_interface_pairs(mortar_face_idx, 1);

          if (interface_A != interface_B)
            continue;

          local_idx_t face_B = m_mortarside_face_interface_pairs(mortar_face_idx, 0);
          bool intersects_along_all = true;;
          for (unsigned i = 0; i < DIM; ++i)
          {
            const aabb_scalar_type min_A = min_corner_A[i];
            const aabb_scalar_type max_A = max_corner_A[i];
            const aabb_scalar_type min_B = m_mortarside_aabbs(mortar_face_idx, aabb_min_idx[i]);
            const aabb_scalar_type max_B = m_mortarside_aabbs(mortar_face_idx, aabb_max_idx[i]);
            intersects_along_all &= ((min_A <= min_B) & (min_B <= max_A)) | ((min_A <= max_B) & (max_B <= max_A));
          }
          if (intersects_along_all)
          {
            if (mode == FILL)
            {
              m_search_results(fill_offset + count, 0) = face_A;
              m_search_results(fill_offset + count, 1) = face_B;
            }
            ++count;
          }
        }
        if (mode == COUNT)
          m_counts(idx) = count;
      }

      struct counts_tag { };
      struct offsets_tag { };
      struct fill_tag { };

      KOKKOS_INLINE_FUNCTION
      void operator() (const counts_tag &tag, const int &i) const {
          find_intersections(i, COUNT);
      }

      KOKKOS_INLINE_FUNCTION
      void operator() (const offsets_tag& tag, const int& i, value_type& offset, const bool& final)  const {
        if(final) {
          m_offsets(i) = offset;
        }
        offset+=m_counts(i);
        if (final && i == m_last_scan_index) {
          m_total_found.d_view() = offset;
        }
      }

      KOKKOS_INLINE_FUNCTION
      void operator() (const fill_tag &tag, const int &i) const {
          find_intersections(i, FILL);
      }

      brute_force_search(search_for_pallet_generating_faces &parent,
                         bounding_boxes_t non_mortarside_aabbs,
                         bounding_boxes_t mortarside_aabbs)
        : m_num_nonmortar_faces(parent.m_non_mortarside_faces.dimension_0()),
          m_last_scan_index(m_num_nonmortar_faces - 1),
          m_non_mortarside_face_interface_pairs(parent.m_non_mortarside_faces),
          m_mortarside_face_interface_pairs(parent.m_mortarside_faces),
          m_non_mortarside_aabbs(non_mortarside_aabbs),
          m_mortarside_aabbs(mortarside_aabbs),
          m_total_found("total_found")
      {
        Kokkos::resize(m_counts, m_num_nonmortar_faces);
        Kokkos::resize(m_offsets, m_num_nonmortar_faces);

        Kokkos::parallel_for(Kokkos::RangePolicy<execution_space, counts_tag>(0, m_num_nonmortar_faces), *this);

        m_total_found. template modify<typename int_dualview_t::t_dev>();
        Kokkos::parallel_scan(Kokkos::RangePolicy<execution_space, offsets_tag>(0, m_num_nonmortar_faces), *this);
        m_total_found. template sync<typename int_dualview_t::t_host>();

        Kokkos::resize(parent.m_search_results, m_total_found.h_view());
        m_search_results = parent.m_search_results;
        Kokkos::parallel_for(Kokkos::RangePolicy<execution_space, fill_tag>(0, m_num_nonmortar_faces), *this);
      }
    };

    search_for_pallet_generating_faces(const surface_mesh_t surface_mesh,
                                       points_t node_coords,
                                       points_t predicted_node_coords,
                                       face_to_interface_and_side_t face_to_interface_and_side,
                                       double epsilon)
        : m_face_to_num_nodes(surface_mesh.m_face_to_num_nodes)
        , m_face_to_nodes(surface_mesh.m_face_to_nodes)
        , m_node_coords(node_coords)
        , m_predicted_node_coords(predicted_node_coords)
        , m_face_to_interface_and_side(face_to_interface_and_side)
        , m_boxes_epsilon(static_cast<float>(epsilon))
    {
        separate_into_sides(*this);

        bounding_boxes_t non_mortarside_aabbs, mortarside_aabbs;
        Kokkos::resize(non_mortarside_aabbs, m_non_mortarside_faces.dimension_0());
        Kokkos::resize(mortarside_aabbs, m_mortarside_faces.dimension_0());

        construct_bounding_boxes(*this, m_non_mortarside_faces, non_mortarside_aabbs);
        construct_bounding_boxes(*this, m_mortarside_faces, mortarside_aabbs);

        brute_force_search(*this, non_mortarside_aabbs, mortarside_aabbs);
    }

};

}

#endif
