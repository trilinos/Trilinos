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

#include "Panzer_ParameterList_ObjectBuilders.hpp"
#include "Panzer_String_Utilities.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_InputPhysicsBlock.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_EquationSet_Factory.hpp"
#include "Teuchos_TestForException.hpp"
#include <sstream>

namespace panzer {
  
  void buildInputPhysicsBlocks(std::map<std::string,panzer::InputPhysicsBlock>& ipb,
			       const Teuchos::ParameterList& p)
  {
    using std::string;
    using std::vector;
    using Teuchos::ParameterList;

    // check for deprecated input file item
    TEUCHOS_TEST_FOR_EXCEPTION(p.isParameter("Physics Blocks"),std::logic_error,
					     "Error - the use of the \"Physics Blocks\" parameter has been deprecated.  Please remove this from your xml input file.");

    // Loop over physics block sublists
    typedef ParameterList::ConstIterator pl_it;
    for (pl_it pb = p.begin(); pb != p.end(); ++pb) {

      std::string physics_block_id = pb->first;

      TEUCHOS_TEST_FOR_EXCEPTION( !(pb->second.isList()), std::logic_error,
				  "Error - All entries in the \"Physics Block\" sublist must be sublists!" );

      ParameterList& pb_sublist = pb->second.getValue(&pb_sublist);

      // check for deprecated input file item
      TEUCHOS_TEST_FOR_EXCEPTION(pb_sublist.isParameter("Number of Equation Sets"),std::logic_error,
				 "Error - the use of the \"Number of Equation Sets\"  parameter in the physics block named \"" << physics_block_id << "\" has been deprecated.  Please remove this from your physics block sublists in you xml input file.");

      // set the physics block id in the physics block
      ipb[physics_block_id].physics_block_id = physics_block_id;

      // Loop over each equation set in the physics block
      for (pl_it eq_set = pb_sublist.begin(); eq_set != pb_sublist.end(); ++eq_set) {

	TEUCHOS_TEST_FOR_EXCEPTION( !(eq_set->second.isList()), std::logic_error,
				    "Error - All entries in the physics block sublist \"" << physics_block_id 
				    << "\" sublist must be equaiton set sublists!" );

	ParameterList& eq_set_sublist = eq_set->second.getValue(&eq_set_sublist);
	panzer::InputEquationSet ies(eq_set_sublist);

	ipb[physics_block_id].eq_sets.push_back(ies);
      }

    }

  }
  
  void buildBCs(std::vector<panzer::BC>& bcs, 
		const Teuchos::ParameterList& p)
  {
    using Teuchos::ParameterList;

    bcs.clear();

    // Check for non-backward compatible change
    TEUCHOS_TEST_FOR_EXCEPTION(p.isParameter("Number of Boundary Conditions"),
			       std::logic_error,
			       "Error - the parameter \"Number of Boundary Conditions\" is no longer valid for the boundary condition sublist.  Please remove this from your input file!");
     
    std::size_t bc_index = 0;
    for (ParameterList::ConstIterator bc_pl=p.begin(); bc_pl != p.end(); ++bc_pl,++bc_index) {
      TEUCHOS_TEST_FOR_EXCEPTION( !(bc_pl->second.isList()), std::logic_error,
				  "Error - All objects in the boundary condition sublist must be BC sublists!" );
      ParameterList& sublist = bc_pl->second.getValue(&sublist);

      panzer::BC bc(bc_index,sublist);
      bcs.push_back(bc);
    }

  }

  void buildBlockIdToPhysicsIdMap(std::map<std::string,std::string>& b_to_p,
				  const Teuchos::ParameterList& p)
  {
    for (Teuchos::ParameterList::ConstIterator entry = p.begin();
	 entry != p.end(); ++entry) {
      std::string dummy_type;
      b_to_p[entry->first] = entry->second.getValue(&dummy_type);
    }
  }

  void buildPhysicsBlocks(const std::map<std::string,std::string>& block_ids_to_physics_ids,
                          const std::map<std::string,Teuchos::RCP<const shards::CellTopology> > & block_ids_to_cell_topo,
                          const std::map<std::string,panzer::InputPhysicsBlock>& physics_id_to_input_physics_blocks,
                          const int base_cell_dimension,
                          const std::size_t workset_size,
                          const panzer::EquationSetFactory & eqset_factory,
			  const Teuchos::RCP<panzer::GlobalData>& global_data,
                          const bool build_transient_support,
                          std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks)
  {
     using Teuchos::RCP;
     using Teuchos::rcp;
  
     // loop over all block id physics id pairs
     std::map<std::string,std::string>::const_iterator itr;
     for (itr = block_ids_to_physics_ids.begin(); itr!=block_ids_to_physics_ids.end();++itr) {
        std::string element_block_id = itr->first;
        std::string physics_block_id = itr->second;

        std::map<std::string,Teuchos::RCP<const shards::CellTopology> >::const_iterator ct_itr
              = block_ids_to_cell_topo.find(element_block_id);
        TEUCHOS_TEST_FOR_EXCEPTION(ct_itr==block_ids_to_cell_topo.end(),
  	  		           std::runtime_error,
  	          		   "Falied to find CellTopology for element block id: \""
  		          	   << element_block_id << "\"!");
        Teuchos::RCP<const shards::CellTopology> cellTopo = ct_itr->second; 
   
        const panzer::CellData volume_cell_data(workset_size, base_cell_dimension,cellTopo);
        
        // find InputPhysicsBlock that corresponds to a paricular block ID
        std::map<std::string,panzer::InputPhysicsBlock>::const_iterator ipb_it = 
              physics_id_to_input_physics_blocks.find(physics_block_id);
  
        // sanity check: passes only if there is a paricular physics ID
        TEUCHOS_TEST_FOR_EXCEPTION(ipb_it == physics_id_to_input_physics_blocks.end(),
  			 std::runtime_error,
  			 "Falied to find InputPhysicsBlock for physics id: "
  			 << physics_block_id << "!");
  
        const panzer::InputPhysicsBlock& ipb = ipb_it->second;
        RCP<panzer::PhysicsBlock> pb = 
	  rcp(new panzer::PhysicsBlock(ipb, element_block_id, volume_cell_data, eqset_factory, global_data, build_transient_support));
        physicsBlocks.push_back(pb);
     }
  }
}
