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

#ifndef PANZER_LOCAL_MESH_CHUNK_HPP
#define PANZER_LOCAL_MESH_CHUNK_HPP

#include "Kokkos_View.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Shards_CellTopology.hpp"
#include "Teuchos_RCP.hpp"
#include <vector>

namespace panzer
{

template <typename LO, typename GO>
struct LocalMeshChunk
{

  int num_cells;
  std::vector<GO> owned_cell_global_indexes;
  std::vector<GO> ghost_cell_global_indexes;
  std::vector<LO> virtual_cell_local_indexes;

  std::string element_block_name;
  std::string sideset_name;

//  // given cell index, gives starting index into local_face_indexes
//  std::vector<int> subcell_indexes_adj;
//
//  // List of subcell indexes for a given
//  std::vector<int> subcell_indexes;

  Kokkos::View<double***,PHX::Device> cell_vertices;

  Teuchos::RCP<const shards::CellTopology> cell_topology;

  int num_faces;
  Kokkos::View<const LO*[2],PHX::Device> face_to_cells;
  Kokkos::View<const LO**,PHX::Device> cell_to_faces;
  Kokkos::View<const LO*[2],PHX::Device> face_to_local_faces;

};

}

#endif
