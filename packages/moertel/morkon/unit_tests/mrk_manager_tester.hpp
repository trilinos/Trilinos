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

#ifndef MORKON_MORKON_MANAGER_TESTER_HPP
#define MORKON_MORKON_MANAGER_TESTER_HPP

#include <mrk_api_classes.hpp>
#include <mrk_compute_normals.hpp>

namespace morkon_exp {

// Morkon_Manager subclass that enables public access to everything.
//
template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
class Morkon_Manager_Tester : public Morkon_Manager<DeviceType, DIM, FACE_TYPE>
{
public:

  typedef Morkon_Manager<DeviceType, DIM, FACE_TYPE>                 morkon_manager;

  typedef typename morkon_manager::interfaces_map_t                interfaces_map_t;
  typedef typename morkon_manager::local_to_global_idx_t      local_to_global_idx_t;
  typedef typename morkon_manager::surface_mesh_t                    surface_mesh_t;
  typedef typename morkon_manager::fields_t                                fields_t;
  typedef typename morkon_manager::on_boundary_table_t          on_boundary_table_t;
  typedef typename morkon_manager::node_support_sets_t          node_support_sets_t;

  typedef typename morkon_manager::local_to_global_idx_hmt  local_to_global_idx_hmt;
  typedef typename morkon_manager::face_to_num_nodes_hmt      face_to_num_nodes_hmt;
  typedef typename morkon_manager::face_to_nodes_hmt              face_to_nodes_hmt;
  typedef typename morkon_manager::points_hmt                            points_hmt;
  typedef typename morkon_manager::on_boundary_table_hmt      on_boundary_table_hmt;

  typedef typename morkon_manager::face_to_interface_and_side_t      face_to_interface_and_side_t;
  typedef typename morkon_manager::face_to_interface_and_side_hmt  face_to_interface_and_side_hmt;

  typedef typename morkon_manager::coarse_search_results_t  contact_search_results_t;
  typedef typename morkon_manager::mortar_pallets_t                  mortar_pallets_t;

  local_to_global_idx_hmt                    hm_node_global_ids;
  local_to_global_idx_hmt                    hm_face_global_ids;
  face_to_interface_and_side_hmt  hm_face_to_interface_and_side;
  face_to_num_nodes_hmt                    hm_face_to_num_nodes;
  face_to_nodes_hmt                            hm_face_to_nodes;
  points_hmt                                     hm_node_coords;
  points_hmt                           hm_predicted_node_coords;
  points_hmt                                    hm_face_normals;

  static Teuchos::RCP< Morkon_Manager_Tester<DeviceType, DIM, FACE_TYPE> >
    MakeInstance(MPI_Comm mpi_comm,  FaceProjectionMethod projection_method, int printlevel);

  MPI_Comm  get_mpi_comm() const     { return morkon_manager::m_mpi_comm; }
  int       get_printlevel() const { return morkon_manager::m_printlevel; }

  template <typename ViewType>
  void CopyToHostMirror(typename ViewType::HostMirror &dest, ViewType &src)
  {
    DeviceType().fence();
    if (dest.shape() != src.shape())
      dest = Kokkos::create_mirror_view(src);

    Kokkos::deep_copy(dest, src);
  }

  void get_node_global_ids() { CopyToHostMirror(hm_node_global_ids, morkon_manager::m_node_global_ids); }
  void get_face_global_ids() { CopyToHostMirror(hm_face_global_ids, morkon_manager::m_face_global_ids); }
  void get_face_to_interface_and_side() {
    CopyToHostMirror(hm_face_to_interface_and_side, morkon_manager::m_face_to_interface_and_side);
  }
  void get_face_to_num_nodes() {
    CopyToHostMirror(hm_face_to_num_nodes, morkon_manager::m_surface_mesh.m_face_to_num_nodes);
  }
  void get_face_to_nodes() {
    CopyToHostMirror(hm_face_to_nodes, morkon_manager::m_surface_mesh.m_face_to_nodes);
  }
  void get_node_coords() {
    CopyToHostMirror(hm_node_coords, morkon_manager::m_fields.m_node_coords);
  }
  void get_predicted_node_coords() {
    CopyToHostMirror(hm_predicted_node_coords, morkon_manager::m_fields.m_predicted_node_coords);
  }
  void get_face_normals() {
    CopyToHostMirror(hm_face_normals, morkon_manager::m_fields.m_face_normals);
  }

  Morkon_Manager_Tester(MPI_Comm mpi_comm, FaceProjectionMethod projection_method, int printlevel)
      : morkon_manager(mpi_comm, projection_method, printlevel) {}

  bool compute_normals() { return morkon_manager::compute_normals(); }

  typename morkon_manager::coarse_search_results_t find_possible_contact_face_pairs() {
    return morkon_manager::find_possible_contact_face_pairs();
  }

  typename morkon_manager::mortar_pallets_t
  compute_contact_pallets(typename morkon_manager::coarse_search_results_t coarse_search_results) {
      return morkon_manager::compute_contact_pallets(coarse_search_results);
  }
};

template  <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE>
Teuchos::RCP< Morkon_Manager_Tester<DeviceType, DIM, FACE_TYPE> >
Morkon_Manager_Tester<DeviceType, DIM, FACE_TYPE>::MakeInstance(MPI_Comm mpi_comm,
                                                                FaceProjectionMethod projection_method,
                                                                int printlevel)
{
  typedef Morkon_Manager_Tester<DeviceType, DIM, FACE_TYPE> morkon_manager_t;

  return Teuchos::RCP<morkon_manager_t>(new morkon_manager_t(mpi_comm, projection_method, printlevel) );
}

}

#endif
