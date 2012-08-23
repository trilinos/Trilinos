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

#include "Panzer_InputEquationSet.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(input_equation_set, basic)
  {
    panzer::InputEquationSet ies_2;
    ies_2.name = "Continuity";
    ies_2.basis = "Q1";
    ies_2.integration_order = 1;
    ies_2.model_id = "fluid";
    ies_2.prefix = "ION_";
    ies_2.params.set("Tau_C", "Codina");

    Teuchos::ParameterList my_list;
    my_list.set("Name", "Continuity");
    my_list.set("Basis", "Q1");
    my_list.set("Integration Order", 1);
    my_list.set("Model ID", "fluid");
    my_list.set("Prefix", "ION_");
    my_list.sublist("Options").set("Tau_C", "Codina");
    
    panzer::InputEquationSet ies_3(my_list);
    
    TEST_EQUALITY(ies_2.name, ies_3.name);
    TEST_EQUALITY(ies_2.basis, ies_3.basis);
    TEST_EQUALITY(ies_2.integration_order, ies_3.integration_order);
    TEST_EQUALITY(ies_2.model_id, ies_3.model_id);
    TEST_EQUALITY(ies_2.prefix, ies_3.prefix);
    TEST_EQUALITY(ies_2.params, ies_3.params);
  }
}
