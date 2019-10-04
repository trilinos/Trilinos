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

#ifndef PANZER_LOCAL_PARTITIONING_UTILITIES_HPP
#define PANZER_LOCAL_PARTITIONING_UTILITIES_HPP

#include "Panzer_LocalMeshInfo.hpp"

#include <vector>

namespace panzer
{

class WorksetDescriptor;

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
