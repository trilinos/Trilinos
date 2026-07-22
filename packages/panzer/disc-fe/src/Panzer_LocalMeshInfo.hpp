// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_LOCAL_MESH_INFO_HPP
#define PANZER_LOCAL_MESH_INFO_HPP

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

    // For side support
    int subcell_index;
    int subcell_dimension;

    // Global cell indexes -> [owned] then [ghosted] then [virtual]
    PHX::View<panzer::GlobalOrdinal*> global_cells;

    // These are the cell indexes in the LocalMeshInfo class
    PHX::View<panzer::LocalOrdinal*> local_cells;

    // Nodes
    PHX::View<double***> cell_nodes;

    // Face to neighbors
    bool has_connectivity;
    PHX::View<panzer::LocalOrdinal*[2]> face_to_cells;
    PHX::View<panzer::LocalOrdinal*[2]> face_to_lidx;
    PHX::View<panzer::LocalOrdinal**> cell_to_faces;
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
