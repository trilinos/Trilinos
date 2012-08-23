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


#ifndef CHARON_WORKSET_BUILDER_DECL_HPP
#define CHARON_WORKSET_BUILDER_DECL_HPP

#include <vector>
#include <map>
#include "Teuchos_RCP.hpp"

namespace shards {
  class CellTopology;
}

namespace panzer {
  
  struct Workset;
  class MeshData;
  class BoundaryCondition;
  class InputPhysicsBlock;
  class PhysicsBlock;
  class BC;

  template<typename ArrayT>
  Teuchos::RCP<std::vector<panzer::Workset> > 
  buildWorksets(const std::string& block_id,
                const Teuchos::RCP<const shards::CellTopology> & blockTopo,
		const std::vector<std::size_t>& local_cell_ids,
		const ArrayT& vertex_coordinates, 
		const panzer::InputPhysicsBlock& ipb,
		std::size_t workset_size,
		int base_cell_dimension);

  template<typename ArrayT>
  Teuchos::RCP<std::vector<panzer::Workset> > 
  buildWorksets(const panzer::PhysicsBlock & physBlk,
		const std::vector<std::size_t>& local_cell_ids,
		const ArrayT& vertex_coordinates, 
		std::size_t workset_size);
  
  template<typename ArrayT>
  Teuchos::RCP<std::map<unsigned,panzer::Workset> >
  buildBCWorkset(const panzer::BC& bc,
                 const Teuchos::RCP<const shards::CellTopology> & blockTopo,
		 const std::vector<std::size_t>& local_cell_ids,
		 const std::vector<std::size_t>& local_side_ids,
		 const ArrayT& vertex_coordinates, 
		 const panzer::InputPhysicsBlock& ipb,
		 unsigned base_cell_dim);

}

#endif
