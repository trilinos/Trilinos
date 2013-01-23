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
#include "Teuchos_ParameterList.hpp"
#include "Panzer_BC.hpp"
#include "Teuchos_TestForException.hpp"

namespace panzer {

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

}
