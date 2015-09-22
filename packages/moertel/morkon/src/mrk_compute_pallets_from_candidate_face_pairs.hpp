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

#ifndef MORKON_EXP_COMPUTE_PALLETS_FROM_CANDIDATE_FACE_PAIRS_H
#define MORKON_EXP_COMPUTE_PALLETS_FROM_CANDIDATE_FACE_PAIRS_H

#include <limits>
#include <mrk_data_types.hpp>

namespace morkon_exp {

template <typename DeviceType, unsigned int DIM>
struct compute_pallets_from_candidate_face_pairs
{
  enum ComputePalletsFunctionMode {COUNT, FILL};

  typedef typename DeviceType::execution_space                          execution_space;
  typedef typename DeviceType::memory_space::size_type                        size_type;

  typedef Kokkos::View<int *, execution_space>                               ints_vec_t;
  typedef Kokkos::DualView<int, execution_space>                         int_dualview_t;
  typedef typename ints_vec_t::value_type                                    value_type;

  typedef Mrk_SurfaceMesh<DeviceType, DIM>                               surface_mesh_t;
  typedef typename surface_mesh_t::face_to_num_nodes_t              face_to_num_nodes_t;
  typedef typename surface_mesh_t::face_to_nodes_t                      face_to_nodes_t;
  typedef typename surface_mesh_t::face_to_nodes_mrat                face_to_nodes_mrat;
  typedef Mrk_Fields<DeviceType, DIM>                                          fields_t;
  typedef typename fields_t::points_t                                          points_t;
  typedef typename fields_t::points_mrat                                    points_mrat;
  typedef MorkonCommonlyUsed<DeviceType, DIM>                           morkon_common_t;
  typedef typename morkon_common_t::coarse_search_results_t        coarse_search_results_t;
  typedef typename morkon_common_t::const_coarse_search_results_t  const_coarse_search_results_t;

  typedef Mrk_MortarPallets<DeviceType, DIM>                           mortar_pallets_t;

  const int                   m_num_input_pairs;
  const int                   m_last_scan_index;
  const_coarse_search_results_t   m_input_pairs;
  ints_vec_t                           m_counts;
  ints_vec_t                          m_offsets;
  int_dualview_t                  m_total_found;

  mortar_pallets_t             m_result_pallets;


  KOKKOS_INLINE_FUNCTION
  bool quick_overlap_check(const int &idx) const
  {
    return false;
  }

  KOKKOS_INLINE_FUNCTION
  void project_onto_nonmortarside(const int &idx) const
  {

  }

  KOKKOS_INLINE_FUNCTION
  void compute_halfspace_membership_tables(const int &idx) const
  {

  }

  KOKKOS_INLINE_FUNCTION
  void compute_pallet_generators(const int &idx) const
  {

  }

  KOKKOS_INLINE_FUNCTION
  void flesh_out_pallet_data(const int &idx) const
  {

  }

  KOKKOS_INLINE_FUNCTION
  void compute_pallets (const int &idx, const ComputePalletsFunctionMode mode) const
  {
    if (!quick_overlap_check(idx))
      return;

    if (mode == COUNT)
    {
      project_onto_nonmortarside(idx);
      compute_halfspace_membership_tables(idx);
      compute_pallet_generators(idx);
    }
    else
    {
      flesh_out_pallet_data(idx);
    }
  }

  struct counts_tag {};
  struct offsets_tag{};
  struct fill_tag{};

  KOKKOS_INLINE_FUNCTION
  void operator() (const counts_tag &tag, const int &i) const {
      compute_pallets(i, COUNT);
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
      compute_pallets(i, FILL);
  }

  compute_pallets_from_candidate_face_pairs(const surface_mesh_t &mesh,
                                            const fields_t &fields,
                                            coarse_search_results_t coarse_search_results)
    : m_num_input_pairs(coarse_search_results.dimension_0()),
      m_last_scan_index(m_num_input_pairs - 1),
      m_input_pairs(coarse_search_results),
      m_total_found("total_found")
  {
    Kokkos::resize(m_counts, m_num_input_pairs);
    Kokkos::resize(m_offsets, m_num_input_pairs);

    Kokkos::parallel_for(Kokkos::RangePolicy<execution_space, counts_tag>(0, m_num_input_pairs ), *this);

    m_total_found. template modify<typename int_dualview_t::t_dev>();
    Kokkos::parallel_scan(Kokkos::RangePolicy<execution_space, offsets_tag>(0, m_num_input_pairs ), *this);
    m_total_found. template sync<typename int_dualview_t::t_host>();

    m_result_pallets.resize(m_total_found.h_view());
    Kokkos::parallel_for(Kokkos::RangePolicy<execution_space, fill_tag>(0, m_num_input_pairs ), *this);
  }
};

}

#endif

