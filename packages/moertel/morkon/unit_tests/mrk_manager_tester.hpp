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

  typedef typename morkon_manager::local_to_global_idx_dvt  local_to_global_idx_dvt;
  typedef typename morkon_manager::face_to_num_nodes_dvt      face_to_num_nodes_dvt;
  typedef typename morkon_manager::face_to_nodes_dvt              face_to_nodes_dvt;
  typedef typename morkon_manager::points_dvt                            points_dvt;
  typedef typename morkon_manager::on_boundary_table_dvt      on_boundary_table_dvt;

  typedef typename morkon_manager::face_to_interface_and_side_t      face_to_interface_and_side_t;
  typedef typename morkon_manager::face_to_interface_and_side_dvt  face_to_interface_and_side_dvt;

  typedef typename morkon_manager::contact_search_results_t  contact_search_results_t;
  typedef typename morkon_manager::mortar_pallets_t                  mortar_pallets_t;

  static Teuchos::RCP< Morkon_Manager_Tester<DeviceType, DIM, FACE_TYPE> > MakeInstance(MPI_Comm mpi_comm, int printlevel);

  MPI_Comm  get_mpi_comm() const     { return morkon_manager::m_mpi_comm; }
  int       get_printlevel() const { return morkon_manager::m_printlevel; }

  Teuchos::RCP<Tpetra::Map<> > get_problem_map_ptr() const { return morkon_manager::m_problem_map; }

  const interfaces_map_t             & get_interfaces_ref() const { return morkon_manager::m_interfaces; }
  const local_to_global_idx_t        & get_node_global_ids_ref() const { return morkon_manager::m_node_global_ids; }
  const local_to_global_idx_t        & get_face_global_ids_ref() const { return morkon_manager::m_face_global_ids; }
  const surface_mesh_t               & get_surface_mesh_ref() const { return  morkon_manager::m_surface_mesh; }
  const fields_t                     & get_fields_ref() const { return morkon_manager::m_fields; }
  const face_to_interface_and_side_t & get_face_to_interface_and_side_ref() const {
      return morkon_manager::m_face_to_interface_and_side;
  }

  on_boundary_table_t          & get_on_boundary_table_ref() const { return morkon_manager::m_is_ifc_boundary_node; }
  node_support_sets_t          & get_node_support_sets_ref() const { return morkon_manager::m_node_support_sets; }

  Morkon_Manager_Tester(MPI_Comm mpi_comm, int printlevel)
      : morkon_manager(mpi_comm, printlevel) {}

  bool internalize_interfaces() { return morkon_manager::internalize_interfaces(); }

  bool migrate_to_device(local_to_global_idx_dvt node_to_global_id,
                         local_to_global_idx_dvt face_to_global_id,
                         face_to_interface_and_side_dvt face_to_interface_and_side,
                         face_to_num_nodes_dvt face_to_num_nodes,
                         face_to_nodes_dvt face_to_nodes,
                         points_dvt node_coords,
                         on_boundary_table_dvt is_node_on_boundary)
  {
      return morkon_manager::migrate_to_device(node_to_global_id,
                                               face_to_global_id,
                                               face_to_interface_and_side,
                                               face_to_num_nodes,
                                               face_to_nodes,
                                               node_coords,
                                               is_node_on_boundary);
  }

  bool compute_face_and_node_normals() { return morkon_manager::compute_face_and_node_normals(); }

  bool find_possible_contact_face_pairs(contact_search_results_t coarse_search_results) {
      return morkon_manager::find_possible_contact_face_pairs(coarse_search_results);
  }

  bool compute_boundary_node_support_sets(contact_search_results_t coarse_search_results) {
      return morkon_manager::compute_boundary_node_support_sets(coarse_search_results);
  }

  bool compute_contact_pallets(contact_search_results_t coarse_search_results,
                               mortar_pallets_t &resulting_pallets) {
      return morkon_manager::compute_contact_pallets(coarse_search_results, resulting_pallets);
  }

  bool integrate_pallets_into_onrank_D(mortar_pallets_t pallets_to_integrate_on) {
      return morkon_manager::integrate_pallets_into_onrank_D(pallets_to_integrate_on);
  }

  bool integrate_pallets_into_onrank_M(mortar_pallets_t pallets_to_integrate_on) {
      return morkon_manager::integrate_pallets_into_onrank_M(pallets_to_integrate_on);
  }
};

template  <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE>
Teuchos::RCP< Morkon_Manager_Tester<DeviceType, DIM, FACE_TYPE> >
Morkon_Manager_Tester<DeviceType, DIM, FACE_TYPE>::MakeInstance(MPI_Comm mpi_comm, int printlevel)
{
  typedef Morkon_Manager_Tester<DeviceType, DIM, FACE_TYPE> morkon_manager_t;

  return Teuchos::RCP<morkon_manager_t>(new morkon_manager_t(mpi_comm, printlevel) );
}


}

#endif
