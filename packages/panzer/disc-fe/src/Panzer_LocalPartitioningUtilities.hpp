// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_LOCAL_PARTITIONING_UTILITIES_HPP
#define PANZER_LOCAL_PARTITIONING_UTILITIES_HPP

#include "Panzer_LocalMeshInfo.hpp"

#include <vector>

namespace panzer
{

class ConnManager;
class WorksetDescriptor;

/**
 * \brief Get the owned, ghost and virtual global cell ids for this process
 *
 * This is designed to support a single nodally connected halo of ghost cells.
 * Virtual cells, however, are defined by face connectivity since they are only used in discontinuous discretizations (e.g. FV and DG).
 *
 * \note This is a collective call over processors
 * \note Local cell indexing is evaluated as [owned][ghost][virtual]
 * \note Global ids for virtual cells exist at the end of the global indexing scheme, and can safely be ignored
 *
 * \param[in] comm Comm for global sync calls
 * \param[in] conn Connectivity manager
 * \param[out] owned_cells Will be filled with owned cell global indexes. These are cells owned by this process.
 * \param[out] ghost_cells Will be filled with ghost cell global indexes. These are halo cells owned by other processors.
 * \param[out] virtual_cells Will be filled with virtual cell global indexes. These are non-existant cells that represent the boundary of the mesh.
 */
void
fillLocalCellIDs(const Teuchos::RCP<const Teuchos::Comm<int>> & comm,
                 panzer::ConnManager & conn,
                 PHX::View<panzer::GlobalOrdinal*> & owned_cells,
                 PHX::View<panzer::GlobalOrdinal*> & ghost_cells,
                 PHX::View<panzer::GlobalOrdinal*> & virtual_cells);

/** Create a set of partitions given a descriptor containing:
 * 1) Volume worksets
 *  - Element Block Name
 *  - Workset Size
 * 2) Sideset worksets
 *  - Element Block Name
 *  - Sideset Name
 *  - Workset Size
 *
 * \param[in] mesh_info Reference to fully constructed mesh_info
 * \param[in] description Workset descriptor defining area to partition
 * \param[out] partitions Set of local mesh partitions for given region of mesh_info
 *
 */
void
generateLocalMeshPartitions(const panzer::LocalMeshInfo & mesh_info,
                            const panzer::WorksetDescriptor & description,
                            std::vector<panzer::LocalMeshPartition> & partitions);

namespace partitioning_utilities
{

/** Create a LocalMeshInfoBase from a parent LocalMeshInfoBase given a set of cell indexes
 *
 * \param[in] parent_info Reference to fully constructed LocalMeshInfoBase
 * \param[in] owned_parent_cells Vector of indexes (in parent's indexing scheme) for child to own
 * \param[out] child_info Child which will be generated
 *
 */
void
setupSubLocalMeshInfo(const panzer::LocalMeshInfoBase & parent_info,
                      const std::vector<panzer::LocalOrdinal> & owned_parent_cells,
                      panzer::LocalMeshInfoBase & child_info);
}

}

#endif
