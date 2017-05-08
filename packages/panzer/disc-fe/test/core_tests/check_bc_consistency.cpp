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

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_LocalTestingHelpers.hpp"
#include "Panzer_CheckBCConsistency.hpp"
#include "Panzer_BC.hpp"

namespace panzer_test {

  TEUCHOS_UNIT_TEST(bc, check_bc_consistency)
  {
    std::vector<std::string> element_block_names{"liquid","solid","gas"};
    std::vector<std::string> sideset_names{"left","right"};

    panzer::BC bc1(0,panzer::BCT_Dirichlet,"left","liquid","MyEquationSetName","MyStrategy");
    panzer::BC bc2(1,panzer::BCT_Neumann,"right","solid","MyEquationSetName","MyStrategy");
    panzer::BC bc_wrong_sideset_name(2,panzer::BCT_Neumann,"bogus","liquid","MyEquationSetName","MyStrategy");
    panzer::BC bc_wrong_element_block_name(3,panzer::BCT_Neumann,"left","bogus","MyEquationSetName","MyStrategy");

    std::vector<panzer::BC> bcs;
    bcs.push_back(bc1);
    bcs.push_back(bc2);

    panzer::checkBCConsistency(element_block_names,
                               sideset_names,
                               bcs);


    bcs.push_back(bc_wrong_sideset_name);
    TEST_THROW(panzer::checkBCConsistency(element_block_names,
                                          sideset_names,
                                          bcs),
               std::runtime_error);

    bcs.clear();
    bcs.push_back(bc1);
    bcs.push_back(bc2);
    bcs.push_back(bc_wrong_element_block_name);
    TEST_THROW(panzer::checkBCConsistency(element_block_names,
                                          sideset_names,
                                          bcs),
               std::runtime_error);
    

  }

}
