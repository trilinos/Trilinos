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

#ifndef MORKON_EXP_API_CLASSES_H
#define MORKON_EXP_API_CLASSES_H

#include <cstdint>

#include <Teuchos_RCP.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <mrk_data_types.hpp>

namespace morkon_exp {

template <typename DeviceType, unsigned int DIM = 3, MorkonFaceType = MRK_TRI3 >
class Interface;

template <typename DeviceType, unsigned int DIM = 3, MorkonFaceType = MRK_TRI3 >
class Morkon_Manager;

template <unsigned int DIM = 3 >
struct Interface_HostSideAdapter;


template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
class Interface : public InterfaceBase
{
  friend  class Morkon_Manager<DeviceType, DIM, FACE_TYPE> ;

  typedef typename DeviceType::execution_space         execution_space;
  typedef Kokkos::View<local_idx_t *, execution_space>     faces_ids_t;
  typedef Kokkos::View<local_idx_t *, execution_space>  faces_ids_dv_t;

public:

  // For pulling data in from the host space.
  bool hsa_add_node(SideEnum which_side, global_idx_t gbl_node_id, const double coords[]);
  bool hsa_add_face(SideEnum which_side, global_idx_t gbl_face_id, int num_nodes, const global_idx_t gbl_node_id[]);

  // No more changes via public API after this.
  bool commited() const { return m_committed; }

private:

  Interface(Morkon_Manager<DeviceType, DIM, FACE_TYPE> *manager);

  Morkon_Manager<DeviceType, DIM, FACE_TYPE>   *m_manager;
  bool                             m_committed;

  std::vector<Interface_HostSideAdapter<DIM> *> m_hs_adapters;
};


template <typename DeviceType, unsigned int DIM, MorkonFaceType FACE_TYPE >
class Morkon_Manager
{
  typedef typename DeviceType::execution_space  execution_space;
  typedef Interface<DeviceType, DIM, FACE_TYPE>     interface_t;
  typedef Teuchos::RCP<interface_t>               interface_ptr;
  typedef std::map<int, interface_ptr>         interfaces_map_t;

  typedef Mrk_SurfaceMesh<DeviceType, DIM>                                   surface_mesh_t;
  typedef typename surface_mesh_t::local_to_global_idx_t              local_to_global_idx_t;
  typedef typename local_to_global_idx_t::HostMirror                local_to_global_idx_hmt;
  typedef typename surface_mesh_t::face_connectivity_t                  face_connectivity_t;
  typedef typename face_connectivity_t::face_to_num_nodes_t             face_to_num_nodes_t;
  typedef typename face_to_num_nodes_t::HostMirror                    face_to_num_nodes_hmt;
  typedef typename face_connectivity_t::face_to_nodes_t                     face_to_nodes_t;
  typedef typename face_to_nodes_t::HostMirror                            face_to_nodes_hmt;

  typedef Mrk_Fields<DeviceType, DIM>                                              fields_t;
  typedef typename fields_t::points_t                                              points_t;
  typedef typename points_t::HostMirror                                          points_hmt;

  typedef Kokkos::View<local_idx_t *[2], execution_space>      face_to_interface_and_side_t;
  typedef typename face_to_interface_and_side_t::HostMirror  face_to_interface_and_side_hmt;
  typedef Kokkos::View<local_idx_t *[3], execution_space>          contact_search_results_t;
  typedef Kokkos::View<bool *, execution_space>                         on_boundary_table_t;
  typedef typename on_boundary_table_t::HostMirror                    on_boundary_table_hmt;
  typedef Kokkos::CrsMatrix<local_idx_t, local_idx_t, DeviceType>       node_support_sets_t;

  typedef Mrk_MortarPallets<DeviceType, DIM>                               mortar_pallets_t;

public:

  static Teuchos::RCP< Morkon_Manager<DeviceType, DIM, FACE_TYPE> > MakeInstance(MPI_Comm mpi_comm, int printlevel);

  bool set_problem_map(Tpetra::Map<> *gp_map);

  // For creating and building Interfaces serially.
  interface_ptr create_interface(int id, int printlevel);

  // Convert serially-built Interfaces information into mesh structure on device.
  // Handle ghosting if needed in future?
  bool commit_interfaces();

  bool mortar_integrate(Tpetra::CrsMatrix<> *D_to_overwrite, Tpetra::CrsMatrix<> *M_to_overwrite);

  // On each MPI rank, returns # of LM dofs assigned by by it.
  int globalize_LM_DOFs();

  bool build_sys_M_and_D(Tpetra::CrsMatrix<> *D_to_overwrite, Tpetra::CrsMatrix<> *M_to_overwrite);

protected:

  // Set in constructor.
  MPI_Comm    m_mpi_comm;
  int       m_printlevel;

  // Input set/manipulated by the application, on the host side for now.
  Teuchos::RCP<Tpetra::Map<> >                 m_problem_map;
  interfaces_map_t                              m_interfaces;

  // On the Device, nodes and faces use local ids.
  local_to_global_idx_t                    m_node_global_ids;
  local_to_global_idx_t                   m_faces_global_ids;
  surface_mesh_t                              m_surface_mesh;
  fields_t                                          m_fields;
  face_to_interface_and_side_t  m_face_to_interface_and_side;  // Might be able to just use separate views for mortar-side face_id

  on_boundary_table_t                 m_is_ifc_boundary_node;  // Is node_id on an interface boundary?
  node_support_sets_t                    m_node_support_sets;

  Morkon_Manager(MPI_Comm mpi_comm, int printlevel);

  // Consider changing the following into free functions in the file that contains
  // the implementation of Morkon_Manager::mortar_integrate().

  bool internalize_interfaces();
  bool migrate_to_device(face_to_interface_and_side_hmt face_to_interface_and_side,
                         face_to_num_nodes_hmt face_to_num_nodes,
                         face_to_nodes_hmt face_to_nodes,
                         points_hmt node_coords,
                         on_boundary_table_hmt is_node_on_boundary);

  bool compute_face_and_node_normals();
  bool find_possible_contact_face_pairs(contact_search_results_t &);
  bool compute_boundary_node_support_sets(contact_search_results_t course_search_results,
                                          node_support_sets_t &support_sets);
  bool compute_contact_pallets(mortar_pallets_t &resulting_pallets);

  // Note that the non-mortar-side integration points needed in computing D are also
  // needed to compute M.  Store and re-use, or re-compute?
  bool integrate_pallets_into_onrank_D(mortar_pallets_t pallets_to_integrate_on,
                                       const node_support_sets_t &support_sets);
  bool integrate_pallets_into_onrank_M(mortar_pallets_t pallets_to_integrate_on,
                                       const node_support_sets_t &support_sets);

};


}

#endif

