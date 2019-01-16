// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef PANZER_WORKSET_UTILITIES_HPP
#define PANZER_WORKSET_UTILITIES_HPP

#include "Panzer_Workset.hpp"
#include "Panzer_WorksetDescriptor.hpp"

#include "Teuchos_RCP.hpp"

#include <vector>

namespace panzer
{
class WorksetDescriptor;
struct LocalMeshInfo;
}

namespace panzer
{

/**
 * Build worksets for a partitioned mesh
 *
 * \param[in] mesh_info Mesh info object
 * \param[in] description Description of workset
 *
 * \returns vector of worksets for the corresponding element block.
 */
Teuchos::RCP<std::vector<panzer::Workset> >  
buildWorksets(const panzer::LocalMeshInfo & mesh_info,
              const panzer::WorksetDescriptor & description,
              const Teuchos::RCP<const OrientationsInterface> & orientations = Teuchos::null);

/**
 * Legacy call for building a BC workset directly from cell indexes and coordinates
 *
 * \param[in] element_block Element block of workset
 * \param[in] sideset Sideset of workset (if no sideset exists use "")
 * \param[in] cell_topology Cell topology of element block
 * \param[in] local_cell_ids Cell indexes for local mesh
 * \param[in] local_side_ids Cell subcell indexes defining side (i.e. face defining boundary)
 * \param[in] vertex_coordinates Cell vertices
 * \param[in] force_side_assembly Force worksets to construct side assembled integrators
 * \param[in] workset_size Size of worksets to create (default is all cells for each subcell index)
 *
 * \return List of worksets defining boundary
 */
Teuchos::RCP<std::vector<Workset> >
buildGroupedSubcellWorksets(const std::string & element_block,
                            const std::string & sideset,
                            const Teuchos::RCP<const shards::CellTopology> & cell_topology,
                            const std::vector<std::size_t>& local_cell_ids,
                            const std::vector<std::size_t>& local_side_ids,
                            const Kokkos::DynRankView<double,PHX::Device> & vertex_coordinates,
                            const bool force_side_assembly,
                            const int workset_size = WorksetSizeType::ALL_ELEMENTS,
                            const Teuchos::RCP<const OrientationsInterface> & orientations = Teuchos::null);

/**
 * Legacy call for building inter-block side assembled worksets directly from cell indexes and coordinates
 *
 * \note There is no issue setting element_block_a = element_block_b, or having sidesets be empty strings
 *
 * \param[in] element_block_a Element block a of boundary
 * \param[in] sideset_a Sideset of a of boundary wrt block a
 * \param[in] cell_topology Cell topology of element block a
 * \param[in] local_cell_ids_a Cell indexes for local mesh in block a
 * \param[in] local_side_ids_a Cell subcell indexes defining boundary (i.e. face defining boundary) in block a
 * \param[in] vertex_coordinates_a Cell vertices in block a
 * \param[in] element_block_b Element block b of boundary
 * \param[in] sideset_b Sideset of a of boundary wrt block b
 * \param[in] cell_topology Cell topology of element block b
 * \param[in] local_cell_ids_b Cell indexes for local mesh in block b
 * \param[in] local_side_ids_b Cell subcell indexes defining boundary (i.e. face defining boundary) in block b
 * \param[in] vertex_coordinates_b Cell vertices in block b
 * \param[in] workset_size Size of worksets to create (default is all cells for each subcell index)
 *
 * \return List of worksets defining boundary
 */
Teuchos::RCP<std::vector<Workset> >
buildGroupedSubcellWorksets(const std::string & element_block_a,
                            const std::string & sideset_a,
                            const Teuchos::RCP<const shards::CellTopology> & cell_topology_a,
                            const std::vector<std::size_t>& local_cell_ids_a,
                            const std::vector<std::size_t>& local_side_ids_a,
                            const Kokkos::DynRankView<double,PHX::Device> & vertex_coordinates_a,
                            const std::string & element_block_b,
                            const std::string & sideset_b,
                            const Teuchos::RCP<const shards::CellTopology> & cell_topology_b,
                            const std::vector<std::size_t>& local_cell_ids_b,
                            const std::vector<std::size_t>& local_side_ids_b,
                            const Kokkos::DynRankView<double,PHX::Device> & vertex_coordinates_b,
                            const int workset_size = WorksetSizeType::ALL_ELEMENTS,
                            const Teuchos::RCP<const OrientationsInterface> & orientations = Teuchos::null);

}

#endif
