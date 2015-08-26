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

#ifndef MORKON_EXP_API_MANAGER_IMPL_H
#define MORKON_EXP_API_MANAGER_IMPL_H

#include <mrk_api_classes.hpp>
#include <mrk_compute_normals.hpp>
#include <mrk_search_for_pallet_generating_faces.hpp>
#include <mrk_interface_impl.hpp>
#include <mrk_interface_host_side_adapter.hpp>

namespace morkon_exp {

template  <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE>
Teuchos::RCP< Morkon_Manager<DeviceType, DIM, FACE_TYPE> >
Morkon_Manager<DeviceType, DIM, FACE_TYPE>::MakeInstance(MPI_Comm mpi_comm, int printlevel)
{
  typedef Morkon_Manager<DeviceType, DIM, FACE_TYPE> morkon_manager_t;

  return Teuchos::RCP<morkon_manager_t>(new morkon_manager_t(mpi_comm, printlevel) );
}

template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE>
bool Morkon_Manager<DeviceType, DIM, FACE_TYPE>::set_problem_map(Tpetra::Map<> *gp_map)
{
  m_problem_map = Teuchos::rcp(new Tpetra::Map<>(*gp_map));
  return true;
}

template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
typename Morkon_Manager<DeviceType, DIM, FACE_TYPE>::interface_ptr
Morkon_Manager<DeviceType, DIM, FACE_TYPE>::create_interface(int id, int printlevel)
{
  if (m_interfaces.find(id) != m_interfaces.end())
  {
    return interface_ptr(0);
  }

  interface_ptr new_interface = interface_ptr(new interface_t(this));
  m_interfaces[id] = new_interface;

  return new_interface;
}


// Convert interface information into mesh structure if needed.  Handle ghosting if needed.
template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
bool Morkon_Manager<DeviceType, DIM, FACE_TYPE>::commit_interfaces()
{
  for (typename interfaces_map_t::iterator ifcs_i =  m_interfaces.begin(); ifcs_i != m_interfaces.end(); ++ifcs_i)
  {
    Interface<DeviceType,DIM,FACE_TYPE> &interface = *ifcs_i->second;

    // No more changes allowed to this interface, even if some other interface
    // cannot be internalized and this function fails!
    interface.m_committed = true;
  }

  return internalize_interfaces();
}


template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
bool Morkon_Manager<DeviceType, DIM, FACE_TYPE>::mortar_integrate(Tpetra::CrsMatrix<> *D_to_overwrite, Tpetra::CrsMatrix<> *M_to_overwrite)
{
  typedef Morkon_Manager<DeviceType, DIM, FACE_TYPE> mgr_t;

  // Using the internal SurfaceMesh, populate
  //   - m_fields.m_node_normals
  //   - m_fields.m_face_normals
  if (!compute_face_and_node_normals())
  {
    return false;
  }

  // Generate vector of (nms_face_id, ms_face_id, interface_id) triples.
  contact_search_results_t  coarse_contacts;
  if (!find_possible_contact_face_pairs(coarse_contacts))
  {
    return false;
  }

  if (!compute_boundary_node_support_sets(coarse_contacts))
  {
    return false;
  }

  // Will our integration scheme require node_support_sets the way the legacy version does?

  mortar_pallets_t pallets_for_integration;
  if (!compute_contact_pallets(coarse_contacts, pallets_for_integration))
  {
    return false;
  }

  if (!integrate_pallets_into_onrank_D(pallets_for_integration))
  {
    return false;
  }

  if (!integrate_pallets_into_onrank_M(pallets_for_integration))
  {
    return false;
  }

  return true;
}


template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
int Morkon_Manager<DeviceType, DIM, FACE_TYPE>::globalize_LM_DOFs()
{
  std::cout << "Need to write globalize_LM_DOFs()" << std::endl;
  return 0;
}

template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
bool Morkon_Manager<DeviceType, DIM, FACE_TYPE>::build_sys_M_and_D(Tpetra::CrsMatrix<> *D_to_overwrite, Tpetra::CrsMatrix<> *M_to_overwrite)
{
  std::cout << "Need to write build_sys_M_and_D(..)" << std::endl;
  return false;
}


template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
Morkon_Manager<DeviceType, DIM, FACE_TYPE>::Morkon_Manager(MPI_Comm mpi_comm, int printlevel)
    : m_mpi_comm(mpi_comm)
    , m_printlevel(printlevel)
{
}


template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
bool Morkon_Manager<DeviceType, DIM, FACE_TYPE>::internalize_interfaces()
{
  // Count up the numbers of nodes, faces, and interfaces.  Fill in the following
  // on the host side, then migrate the data to the device.

  local_to_global_idx_dvt                  node_to_global_id;
  local_to_global_idx_dvt                  face_to_global_id;
  face_to_interface_and_side_dvt  face_to_interface_and_side;
  face_to_num_nodes_dvt                    face_to_num_nodes;
  face_to_nodes_dvt                            face_to_nodes;
  points_dvt                                     node_coords;
  points_dvt                           predicted_node_coords;
  on_boundary_table_dvt                  is_node_on_boundary;

  for (typename interfaces_map_t::iterator ifcs_i = m_interfaces.begin(); ifcs_i != m_interfaces.end(); ++ifcs_i)
  {
    Interface<DeviceType,DIM,FACE_TYPE> &interface = *ifcs_i->second;

    for (unsigned hsa_i = 0; hsa_i < interface.m_hs_adapters.size(); ++hsa_i)
    {
      Interface_HostSideAdapter<DIM> *adapter_rp = interface.m_hs_adapters[hsa_i];

      if (!adapter_rp)
        continue;

      Interface_HostSideAdapter<DIM> &adapter = *adapter_rp;

      // Convert the adapter's version of the side to one in terms of the internal face_ids in the
      // surface_mesh, inserting nodes and faces as needed.
      // WRITE ME!
      std::cout << "Let's see what faces are on this adapter!" << std::endl;
      adapter.m_faces.begin();
    }
  }

  // Now that the data is ready to move to the device side, commit it to there.
  return migrate_to_device(node_to_global_id, face_to_global_id,
                           face_to_interface_and_side, face_to_num_nodes, face_to_nodes,
                           node_coords, predicted_node_coords, is_node_on_boundary);
}

template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
bool
Morkon_Manager<DeviceType, DIM, FACE_TYPE>::migrate_to_device(
                                        local_to_global_idx_dvt node_to_global_id,
                                        local_to_global_idx_dvt face_to_global_id,
                                        face_to_interface_and_side_dvt face_to_interface_and_side,
                                        face_to_num_nodes_dvt face_to_num_nodes,
                                        face_to_nodes_dvt face_to_nodes,
                                        points_dvt node_coords,
                                        points_dvt predicted_node_coords,
                                        on_boundary_table_dvt is_node_on_boundary)
{
    std::cout << "In migrate_to_device()!"  << std::endl;

    face_to_global_id.template modify<typename local_to_global_idx_dvt::t_host>();
    node_to_global_id.template modify<typename local_to_global_idx_dvt::t_host>();
    face_to_interface_and_side.template modify<typename face_to_interface_and_side_dvt::t_host>();
    face_to_num_nodes.template modify<typename face_to_num_nodes_dvt::t_host>();
    face_to_nodes.template modify<typename face_to_nodes_dvt::t_host>();
    node_coords.template modify<typename points_dvt::t_host>();
    predicted_node_coords.template modify<typename points_dvt::t_host>();
    is_node_on_boundary.template modify<typename on_boundary_table_dvt::t_host>();

    face_to_global_id.template sync<typename local_to_global_idx_dvt::t_dev>();
    node_to_global_id.template sync<typename local_to_global_idx_dvt::t_dev>();
    face_to_interface_and_side.template sync<typename face_to_interface_and_side_dvt::t_dev>();
    face_to_num_nodes.template sync<typename face_to_num_nodes_dvt::t_dev>();
    face_to_nodes.template sync<typename face_to_nodes_dvt::t_dev>();
    node_coords.template sync<typename points_dvt::t_dev>();
    predicted_node_coords.template sync<typename points_dvt::t_dev>();
    is_node_on_boundary.template sync<typename on_boundary_table_dvt::t_dev>();

    m_node_global_ids                  = node_to_global_id.d_view;
    m_face_global_ids                  = face_to_global_id.d_view;
    m_face_to_interface_and_side       = face_to_interface_and_side.d_view;
    m_surface_mesh.m_face_to_num_nodes = face_to_num_nodes.d_view;
    m_surface_mesh.m_face_to_nodes     = face_to_nodes.d_view;
    m_fields.m_node_coords             = node_coords.d_view;
    m_fields.m_predicted_node_coords   = predicted_node_coords.d_view;
    m_is_ifc_boundary_node             = is_node_on_boundary.h_view;

    // TO DO: compute upward connectivities on the surface mesh, if needed.

    return true;
}

template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE>
bool Morkon_Manager<DeviceType, DIM, FACE_TYPE>::compute_face_and_node_normals()
{
  // We can make this function provide a useful return value having the implementations
  // do a parallel_reduce with a num_errs reduction variable as argument.

  compute_face_normals<DeviceType, DIM, FACE_TYPE>(m_surface_mesh, m_fields);

  compute_node_normals_from_faces<DeviceType, DIM >(m_surface_mesh, m_fields);

  return true;
}

template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
bool Morkon_Manager<DeviceType, DIM, FACE_TYPE>::find_possible_contact_face_pairs(contact_search_results_t search_results)
{
  const double bounding_boxes_epsilon = 0.001;

  search_for_pallet_generating_faces<DeviceType, DIM>(m_surface_mesh,
                                                      m_fields.m_node_coords,
                                                      m_fields.m_predicted_node_coords,
                                                      m_face_to_interface_and_side,
                                                      bounding_boxes_epsilon,
                                                      search_results);
  return false;
}

template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
bool 
Morkon_Manager<DeviceType, DIM, FACE_TYPE>::compute_boundary_node_support_sets(contact_search_results_t course_search_results)
{
  std::cout << "Need to write compute_boundary_node_support_sets()" << std::endl;
  return false;
}

template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
bool Morkon_Manager<DeviceType, DIM, FACE_TYPE>::compute_contact_pallets(contact_search_results_t course_search_results,
                                                                         mortar_pallets_t &resulting_pallets)
{
  // In the Serial prototype and the Cuda version, we can use atomic fetch and adds to allocate space for the
  // pallets resulting from pair of faces.
  //
  // In the OpenMP friendly version, the counting pass functor does the minimal work possible to figure out
  // how many pallets for each pair, and the compute-and-fill pass does the complete version of the work.

  std::cout << "Need to write compute_contact_pallets()" << std::endl;
  return false;
}

template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
bool Morkon_Manager<DeviceType, DIM, FACE_TYPE>::
integrate_pallets_into_onrank_D(mortar_pallets_t pallets_to_integrate_on)
{
  std::cout << "Need to write integrate_pallets_into_onrank_D()" << std::endl;
  return false;
}

template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
bool Morkon_Manager<DeviceType, DIM, FACE_TYPE>::
integrate_pallets_into_onrank_M(mortar_pallets_t pallets_to_integrate_on)
{
  std::cout << "Need to write integrate_pallets_into_onrank_M()" << std::endl;
  return false;
}


} // namespace morkon_exp


#endif
