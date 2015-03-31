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

template <typename DeviceType, unsigned int DIM = 3 >
class Interface;

template <typename DeviceType, unsigned int DIM = 3 >
class Morkon_Manager;

template <unsigned int DIM = 3 >
struct Interface_HostSideAdapter;


template <typename DeviceType, unsigned int DIM >
class Interface : public InterfaceBase
{
  friend  class Morkon_Manager<DeviceType, DIM> ;

  typedef typename DeviceType::execution_space  execution_space;
  typedef Kokkos::View<local_idx_t *, execution_space>     faces_ids_t;
  typedef Kokkos::View<local_idx_t *, execution_space>  faces_ids_dv_t;

public:

  // For when faces and other mesh data are already in the Kokkos::execution_space.
  bool define_side(SideEnum which_side, faces_ids_t faces_on_side);  // Global ids or local ids?!

  // For pulling data in from the host space.
  bool hsa_add_node(SideEnum which_side, global_idx_t gbl_node_id, const double coords[]);
  bool hsa_add_segment(SideEnum which_side, global_idx_t gbl_seg_id, int num_nodes, const global_idx_t glb_nids[]);

  // No more changes via public API after this.
  bool commited() const { return m_committed; }

  // Declare that this is a multi-MPI-rank Interface.  Must be called consistently on all ranks that know about this interface,
  bool set_distributed();

private:

  Interface(Morkon_Manager<DeviceType, DIM> *manager);

  Morkon_Manager<DeviceType, DIM>   *m_manager;
  bool                             m_committed;
  bool                           m_distributed;
  std::vector<faces_ids_t>             m_sides;

  std::vector<Interface_HostSideAdapter<DIM> *> m_hs_adapters;

};


template <typename DeviceType, unsigned int DIM >
class Morkon_Manager
{
  typedef typename DeviceType::execution_space  execution_space;
  typedef Interface<DeviceType, DIM>                interface_t;
  typedef Teuchos::RCP<interface_t>               interface_ptr;
  typedef std::map<int, interface_ptr>         interfaces_map_t;

  typedef Kokkos::View<local_idx_t *[2], execution_space>         faces2interface_t;
  typedef Kokkos::DualView<local_idx_t *[2], execution_space>  faces2interface_dv_t;

  typedef Mrk_SkinOnlyMesh<DeviceType, DIM>    skin_only_mesh_t;
  typedef skin_only_mesh_t                local_to_global_idx_t;
  typedef Mrk_Fields<DeviceType, DIM>                  fields_t;
  typedef typename fields_t::points_t                  points_t;

  typedef Mrk_MortarPallets<DeviceType, DIM>   mortar_pallets_t;

  typedef Kokkos::CrsMatrix<local_idx_t, local_idx_t, DeviceType>         segment_interface_mat_t;
  typedef Kokkos::View<local_idx_t *[3], execution_space>  contact_search_results_t;
  typedef Kokkos::CrsMatrix<bool, local_idx_t, DeviceType>                    on_boundary_table_t;
  typedef Kokkos::CrsMatrix<local_idx_t, local_idx_t, DeviceType>             node_support_sets_t;

public:

  static Teuchos::RCP< Morkon_Manager<DeviceType, DIM> > MakeInstance(MPI_Comm mpi_comm, int printlevel);

  bool set_problem_map(Tpetra::Map<> *gp_map);

  // For creating and building Interfaces serially.
  interface_ptr create_interface(int id, int printlevel);

  // Convert serially-built Interfaces information into mesh structure if needed.
  // Handle ghosting if needed in future?
  bool commit_interfaces();

  // When data is already on device; called at end of commit_interfaces().
  bool declare_all_interfaces(segment_interface_mat_t segs_in_ifcs, 
                              skin_only_mesh_t dense_idx_mesh,
                              points_t node_coords,
                              local_to_global_idx_t non_dense_node_ids,
                              on_boundary_table_t boundary_node_table);

  bool mortar_integrate(Tpetra::CrsMatrix<> *D_to_overwrite, Tpetra::CrsMatrix<> *M_to_overwrite);

  // On each MPI rank, returns # of LM dofs assigned by by it.
  int globalize_LM_DOFs();

  bool build_sys_M_and_D(Tpetra::CrsMatrix<> *D_to_overwrite, Tpetra::CrsMatrix<> *M_to_overwrite);

private:

  MPI_Comm    m_mpi_comm;
  int       m_printlevel;

  Teuchos::RCP<Tpetra::Map<> >  m_problem_map;
  interfaces_map_t               m_interfaces;
  skin_only_mesh_t                m_skin_mesh;
  segment_interface_mat_t  m_seg_ifc_side_mat;
  fields_t                           m_fields;

  local_to_global_idx_t  m_non_dense_node_ids;

  on_boundary_table_t  m_is_ifc_boundary_node;  // Is node_id on an interface boundary?

  Morkon_Manager(MPI_Comm mpi_comm, int printlevel);

  // Consider changing the following into free functions in the file that contains
  // the implementation of Morkon_Manager::mortar_integrate().

  bool internalize_interfaces();
  bool compute_face_and_node_normals();
  bool find_possible_contact_face_pairs(contact_search_results_t &);
  bool compute_boundary_node_support_sets(contact_search_results_t course_search_results,
                                          node_support_sets_t &support_sets);
  bool compute_contact_pallets(mortar_pallets_t &resulting_pallets);

  // Note that the non-mortar-side integration points needed in computing D are also
  // needed to compute M.  Store and re-use, or re-compute?
  bool integrate_pallets_into_onrank_D(mortar_pallets_t pallets_to_integrate_on /* additional arg(s)? */ );
  bool integrate_pallets_into_onrank_M(mortar_pallets_t pallets_to_integrate_on /* additional arg(s)? */ );

};


}

#endif

