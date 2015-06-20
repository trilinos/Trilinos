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
#include <mrk_manager_functors.hpp>
#include <mrk_interface_impl.hpp>

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

    if (interface.m_distributed)
    {
      return false;
    }
    interface.m_committed = true;
  }

  // Generate Kokkos::CrsMatrixOfKVecs that maps internal face_ids (e.g., face_ids) to sorted vectors of
  // (interface_id, side) pairs.
  //
  // As needed, generates and populates
  //   - m_non_dense_node_ids
  //   - m_skin_mesh
  //   - m_fields.m_node_coords
  //   - m_face_ifc_side, the sparce matrix that maps from face_id to its (interface_id, side) pair(s).
  return internalize_interfaces();
}


template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
bool
Morkon_Manager<DeviceType, DIM, FACE_TYPE>::declare_all_interfaces(face_interface_mat_t faces_in_ifcs, 
                                                        skin_only_mesh_t dense_idx_mesh,
                                                        points_t node_coords,
                                                        local_to_global_idx_t non_dense_node_ids,
                                                        on_boundary_table_t boundary_node_table)
{
  m_skin_mesh            =      dense_idx_mesh;
  m_face_ifc_side_mat     =        faces_in_ifcs;
  m_fields.m_node_coords =         node_coords;
  m_non_dense_node_ids   =  non_dense_node_ids;
  m_is_ifc_boundary_node = boundary_node_table;
}

template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
bool Morkon_Manager<DeviceType, DIM, FACE_TYPE>::mortar_integrate(Tpetra::CrsMatrix<> *D_to_overwrite, Tpetra::CrsMatrix<> *M_to_overwrite)
{
  typedef Morkon_Manager<DeviceType, DIM, FACE_TYPE> mgr_t;

  // Using the internal SkinOnlyMesh, populate
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

  node_support_sets_t node_support_sets;
  if (!compute_boundary_node_support_sets(coarse_contacts, node_support_sets))
  {
    return false;
  }

  mortar_pallets_t pallets_for_integration;
  if (!compute_contact_pallets(pallets_for_integration))
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
  : m_mpi_comm(mpi_comm), m_printlevel(printlevel)
{
}

template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
bool Morkon_Manager<DeviceType, DIM, FACE_TYPE>::internalize_interfaces()
{
  // Count up the numbers of nodes, faces, and interfaces.  Want to be able to go
  // from a (face_id, face_id) pair to a (non_mrtr_sd_face_id, mrtr_sd_face_id, ifc_id) triple.
  // This could be done by constructing a crs_repn from face_id to (ifc_id, side) here, with
  // the ifc_ids sorted for each face.  When the face-face searches are done, traverse the
  // lists for the pair to find the matching ifc_id.

  // Thus, gill in the following:
  face_interface_mat_t         faces_in_ifcs;
  skin_only_mesh_t              dense_idx_mesh;
  points_t                         node_coords;
  local_to_global_idx_t     non_dense_node_ids;
  on_boundary_table_t     is_ifc_boundary_node;

  for (typename interfaces_map_t::iterator ifcs_i = m_interfaces.begin(); ifcs_i != m_interfaces.end(); ++ifcs_i)
  {
    Interface<DeviceType,DIM,FACE_TYPE> &interface = *ifcs_i->second;

    for (int hsa_i = 0; hsa_i < interface.m_hs_adapters.size(); ++hsa_i)
    {
      Interface_HostSideAdapter<DIM> &adapter = *interface.m_hs_adapters[hsa_i];

      // Convert the adapter's version of the side to one in terms of the internal face_ids in the
      // skin_only_mesh, inserting nodes and faces as needed.
      // WRITE ME!
    }
    interface.m_committed = true;
  }

  // Now call declare_all_interfaces(..) since we have the data on the Device side we need for
  // its arguments.
  bool ok = declare_all_interfaces(faces_in_ifcs, dense_idx_mesh, node_coords, non_dense_node_ids,
                                   is_ifc_boundary_node);
  return false;
}

template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE>
bool Morkon_Manager<DeviceType, DIM, FACE_TYPE>::compute_face_and_node_normals()
{
  // We can make this function provide a useful return value having the implementations
  // do a parallel_reduce with a num_errs reduction variable as argument.

  compute_face_normals<DeviceType, DIM, FACE_TYPE>(m_skin_mesh, m_fields);

  compute_node_normals_from_faces<DeviceType, DIM >(m_skin_mesh, m_fields);

  return true;
}

template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
bool Morkon_Manager<DeviceType, DIM, FACE_TYPE>::find_possible_contact_face_pairs(contact_search_results_t &)
{
  // Implement with a functor over the faces.
  std::cout << "Need to write :find_possible_contact_face_pairs()" << std::endl;
  return false;
}

template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
bool 
Morkon_Manager<DeviceType, DIM, FACE_TYPE>::compute_boundary_node_support_sets(contact_search_results_t course_search_results,
                                                                    node_support_sets_t &support_sets)
{
  std::cout << "Need to write compute_boundary_node_support_sets()" << std::endl;
  return false;
}

template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
bool Morkon_Manager<DeviceType, DIM, FACE_TYPE>::compute_contact_pallets(mortar_pallets_t &resulting_pallets)
{
  std::cout << "Need to write compute_contact_pallets()" << std::endl;
  return false;
}

template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
bool Morkon_Manager<DeviceType, DIM, FACE_TYPE>::
integrate_pallets_into_onrank_D(mortar_pallets_t pallets_to_integrate_on /* additional arg(s)? */ )
{
  std::cout << "Need to write integrate_pallets_into_onrank_D()" << std::endl;
  return false;
}

template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
bool Morkon_Manager<DeviceType, DIM, FACE_TYPE>::
integrate_pallets_into_onrank_M(mortar_pallets_t pallets_to_integrate_on /* additional arg(s)? */ )
{
  std::cout << "Need to write integrate_pallets_into_onrank_M()" << std::endl;
  return false;
}


} // namespace morkon_exp


#endif
