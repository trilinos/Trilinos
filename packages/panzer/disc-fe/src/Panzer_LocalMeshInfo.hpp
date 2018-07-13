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
#include "Kokkos_DynRankView.hpp"

#include "Phalanx_KokkosDeviceTypes.hpp"

#include "Shards_CellTopology.hpp"

#include "Teuchos_RCP.hpp"

#include <string>

namespace panzer
{

/** Base class for LocalMeshInfo structures
 *
 * TODO: Replace with Connectivity manager
 */
template <typename LO, typename GO>
struct LocalMeshInfoBase
{

  LO num_owned_cells;
  LO num_ghstd_cells;
  LO num_virtual_cells;

  // Global cell indexes -> [owned] then [ghosted] then [virtual]
  Kokkos::View<GO*> global_cells;

  // These are the cell indexes in the LocalMeshInfo class
  Kokkos::View<LO*> local_cells;

  // Vertices
  Kokkos::View<double***,PHX::Device> cell_vertices;

  // Face to neighbors
  Kokkos::View<LO*[2]> face_to_cells;
  Kokkos::View<LO*[2]> face_to_lidx;
  Kokkos::View<LO**> cell_to_faces;

};

/** Partition of LocalMeshInfo
 *
 * Used for generating worksets
 *
 */
template <typename LO, typename GO>
struct LocalMeshPartition:
    public LocalMeshInfoBase<LO,GO>
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
template <typename LO, typename GO>
struct LocalMeshSidesetInfo:
    public LocalMeshInfoBase<LO,GO>
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
template <typename LO, typename GO>
struct LocalMeshBlockInfo:
    public LocalMeshInfoBase<LO,GO>
{
  std::string element_block_name;

  Teuchos::RCP<const shards::CellTopology> cell_topology;

};

/** Entire mesh found on a local process
 *
 */
template <typename LO, typename GO>
struct LocalMeshInfo:
    public LocalMeshInfoBase<LO,GO>
{

  // Element block -> block info
  std::map<std::string, LocalMeshBlockInfo<LO,GO> > element_blocks;

  // Element block, sideset -> sideset info
  std::map<std::string, std::map<std::string, LocalMeshSidesetInfo<LO,GO> > > sidesets;

};

}

#endif
