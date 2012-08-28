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

#ifndef PANZER_PARAMETER_LIST_OBJECT_BUILDERS_HPP
#define PANZER_PARAMETER_LIST_OBJECT_BUILDERS_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include <map>
#include <vector>
#include <string>

#include "Shards_CellTopology.hpp"

namespace panzer {

  class InputPhysicsBlock;
  class BC;
  class EquationSetFactory;
  class PhysicsBlock;
  class GlobalData;

  void buildInputPhysicsBlocks(std::map<std::string,panzer::InputPhysicsBlock>& ipb,
			       const Teuchos::ParameterList& p);
  
  void buildBCs(std::vector<panzer::BC>& bcs, 
		const Teuchos::ParameterList& p);

  void buildBlockIdToPhysicsIdMap(std::map<std::string,std::string>& b_to_p, 
				  const Teuchos::ParameterList& p);

  /** This builds the physics block objects. In particular a map from the 
    * element block Id to the physics block is constructed.
    *
    * \param[in] block_ids_to_physics_ids A mapping from element block IDs to
    *                                     physics IDs 
    * \param[in] block_ids_to_cell_topo A mapping from element block IDs to
    *                                   their cell topology
    * \param[in] physics_id_to_input_physics_blocks This takes the physics IDs and
    *                                               maps to an input physics block.
    *                                               Essentially this is used to construct
    *                                               the physics block object
    * \param[in] base_cell_dimension What is the dimension of the volume element
    * \param[in] workset_size Size of volume worksets used in Phalanx assembly
    * \param[in] eqset_factory User defined factory for building equation sets
    * \param[in] build_transient_support Tells the equation sets to build evaluators 
    *                                    for a transient analysis
    * \param[in,out] physicsBlock A vector of pointers to the physics blocks
    */
  void buildPhysicsBlocks(const std::map<std::string,std::string>& block_ids_to_physics_ids,
                          const std::map<std::string,Teuchos::RCP<const shards::CellTopology> > & block_id_to_cell_topo,
                          const std::map<std::string,panzer::InputPhysicsBlock>& physics_id_to_input_physics_blocks,
                          const int base_cell_dimension, 
                          const std::size_t workset_size,
                          const panzer::EquationSetFactory & eqset_factory,
			  const Teuchos::RCP<panzer::GlobalData>& global_data,
	                  const bool build_transient_support,
                          std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks);

}

#endif
