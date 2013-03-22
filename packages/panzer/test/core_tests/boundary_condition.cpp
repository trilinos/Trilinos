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

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_ParameterList.hpp>

#include "Panzer_BC.hpp"
#include <iostream>
#include <sstream>
#include <map>

namespace panzer {

  TEUCHOS_UNIT_TEST(bc, nonmember_ctor)
  {
    Teuchos::ParameterList bc_params;

    std::vector<panzer::BC> bcs;
    Teuchos::ParameterList& bc_0 = bc_params.sublist("BC 0");
    bc_0.set("Type", "Dirichlet");
    bc_0.set("Sideset ID", "4");
    bc_0.set("Element Block ID", "fluid");
    bc_0.set("Equation Set Name", "UX");
    bc_0.set("Strategy", "Constant");
    bc_0.sublist("Data").set("Value",1.0);
    Teuchos::ParameterList& bc_1 = bc_params.sublist("BC 1");
    bc_1.set("Type", "Neumann");
    bc_1.set("Sideset ID", "4");
    bc_1.set("Element Block ID", "fluid");
    bc_1.set("Equation Set Name", "UX");
    bc_1.set("Strategy", "Constant");
    bc_1.sublist("Data").set("Value",1.0);

    panzer::buildBCs(bcs, bc_params);

    TEST_EQUALITY(bcs.size(), 2);
    TEST_EQUALITY(bcs[0].bcID(), 0);
    TEST_EQUALITY(bcs[1].bcID(), 1);
    TEST_EQUALITY(bcs[0].bcType(), panzer::BCT_Dirichlet);
    TEST_EQUALITY(bcs[1].bcType(), panzer::BCT_Neumann);
  }


  TEUCHOS_UNIT_TEST(bc, neumann_no_param_list)
  {

    std::size_t bc_id = 0;
    panzer::BCType neumann = BCT_Dirichlet;
    std::string sideset_id = "4";
    std::string element_block_id = "fluid";
    std::string dof_name = "UX";
    std::string strategy = "Constant";
    double value = 5.0;
    Teuchos::ParameterList p;
    p.set("Value",value);
    panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name, 
		  strategy, p);

    TEST_EQUALITY(bc.bcID(), bc_id);
    TEST_EQUALITY(bc.bcType(), neumann);
    TEST_EQUALITY(bc.sidesetID(), sideset_id);
    TEST_EQUALITY(bc.elementBlockID(), element_block_id);
    TEST_EQUALITY(bc.equationSetName(), dof_name);

    std::stringstream s;
    s << bc << std::endl;
  }

  TEUCHOS_UNIT_TEST(bc, dirichlet_with_param_list)
  {
    std::size_t bc_id = 0;
    panzer::BCType dirichlet = BCT_Dirichlet;
    std::string sideset_id = "4";
    std::string element_block_id = "fluid";
    std::string dof_name = "UX";
    std::string strategy = "Constant";
    double value = 5.0;
    Teuchos::ParameterList p;
    p.set("Value",value);
    panzer::BC bc(bc_id, dirichlet, sideset_id, element_block_id, dof_name, 
		  strategy, p);

    TEST_EQUALITY(bc.bcID(), bc_id);
    TEST_EQUALITY(bc.bcType(), dirichlet);
    TEST_EQUALITY(bc.sidesetID(), sideset_id);
    TEST_EQUALITY(bc.elementBlockID(), element_block_id);
    TEST_EQUALITY(bc.equationSetName(), dof_name);

    std::stringstream s;
    s << bc << std::endl;
  }

  TEUCHOS_UNIT_TEST(bc, dirichlet_complete_param_list)
  {
    Teuchos::ParameterList p;
    p.set("Type", "Dirichlet");
    p.set("Sideset ID", "4");
    p.set("Element Block ID", "fluid");
    p.set("Equation Set Name", "UX");
    p.set("Strategy", "Constant");
    p.sublist("Data").set("Value",1.0);

    panzer::BC bc(0,p);

    TEST_EQUALITY(bc.bcID(), 0);
    TEST_EQUALITY(bc.bcType(), BCT_Dirichlet);
    TEST_EQUALITY(bc.sidesetID(), "4");
    TEST_EQUALITY(bc.elementBlockID(), "fluid");
    TEST_EQUALITY(bc.equationSetName(), "UX");

    std::stringstream s;
    s << bc << std::endl;
  }

  TEUCHOS_UNIT_TEST(bc, map_comparitor)
  {
    using panzer::BC;

    BC bc1(0,BCT_Dirichlet,"3","fluid","VELOCITY","Constant");
    BC bc2(1,BCT_Dirichlet,"3","fluid","VELOCITY","Constant");
    BC bc3(2,BCT_Dirichlet,"3","fluid","VELOCITY","Constant");

    std::map<BC,int,panzer::LessBC> my_bcs;

    my_bcs[bc1] = 8;
    my_bcs[bc2] = 2;
    my_bcs[bc3] = 11;

    TEST_EQUALITY(my_bcs[bc1], 8);
    TEST_EQUALITY(my_bcs[bc2], 2);
    TEST_EQUALITY(my_bcs[bc3], 11);

    TEST_INEQUALITY(my_bcs[bc3], 4);    
  }
}
