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

#ifndef PANZER_LOCAL_MESH_INFO_HPP
#define PANZER_LOCAL_MESH_INFO_HPP

#include "Kokkos_View.hpp"
#include "PanzerCore_config.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Shards_CellTopology.hpp"
#include "Teuchos_RCP.hpp"
#include <string>

namespace panzer
{

  /** Base class for LocalMeshInfo structures */
  struct LocalMeshInfoBase
  {
    panzer::LocalOrdinal num_owned_cells;
    panzer::LocalOrdinal num_ghstd_cells;
    panzer::LocalOrdinal num_virtual_cells;

    // Global cell indexes -> [owned] then [ghosted] then [virtual]
    Kokkos::View<panzer::GlobalOrdinal*> global_cells;

    // These are the cell indexes in the LocalMeshInfo class
    Kokkos::View<panzer::LocalOrdinal*> local_cells;

    // Vertices
    Kokkos::View<double***,PHX::Device> cell_vertices;

    // Face to neighbors
    Kokkos::View<panzer::LocalOrdinal*[2]> face_to_cells;
    Kokkos::View<panzer::LocalOrdinal*[2]> face_to_lidx;
    Kokkos::View<panzer::LocalOrdinal**> cell_to_faces;
  };

  /** Partition of LocalMeshInfo, used for generating worksets */
  struct LocalMeshPartition : public LocalMeshInfoBase
  {
    std::string element_block_name;
    Teuchos::RCP<const shards::CellTopology> cell_topology;

    // In case this is a sideset
    std::string sideset_name;
  };

  /** Portion of LocalMeshInfo associated with sidesets
   *
   * Used to represent a sideset found on the local process
   *
   */
  struct LocalMeshSidesetInfo : public LocalMeshInfoBase
  {
    std::string sideset_name;

    std::string element_block_name;

    // Cell topology associated with element_block_name
    Teuchos::RCP<const shards::CellTopology> cell_topology;
  };

  /** Portion of LocalMeshInfo associated with element block
   *
   * Used to represent an element block found on the local process
   *
   */
  struct LocalMeshBlockInfo : public LocalMeshInfoBase
  {
    std::string element_block_name;

    Teuchos::RCP<const shards::CellTopology> cell_topology;
  };

  /** Entire mesh found on a local process */
  struct LocalMeshInfo : public LocalMeshInfoBase
  {
    // Element block -> block info
    std::map<std::string, LocalMeshBlockInfo> element_blocks;

    // Element block, sideset -> sideset info
    std::map<std::string, std::map<std::string,LocalMeshSidesetInfo>> sidesets;
  };

}

#endif
