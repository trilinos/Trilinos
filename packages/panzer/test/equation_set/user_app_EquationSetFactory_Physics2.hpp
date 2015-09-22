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


#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_EquationSet_Factory_Defines.hpp"
#include "Panzer_CellData.hpp"

// Add my equation sets here
#include "user_app_EquationSet_Energy.hpp"

namespace user_app {

  /* comment out - already declared in physics 1
  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER("Energy 1", user_app::EquationSet_Energy,
					EquationSet_Energy)
  */
  
  class MyFactory2 : public panzer::EquationSetFactory {

    bool m_throw_on_failure;

  public:

    MyFactory2(bool throw_on_failure) : m_throw_on_failure(throw_on_failure) {}

    Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> >
    buildEquationSet(const Teuchos::RCP<Teuchos::ParameterList>& params,
		     const int& default_integration_order,
		     const panzer::CellData& cell_data,
		     const Teuchos::RCP<panzer::GlobalData>& global_data,
		     const bool build_transient_support) const
    {
      Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > eq_set= 
	Teuchos::rcp(new panzer::EquationSet_TemplateManager<panzer::Traits>);
      
      bool found = false;
      
      PANZER_BUILD_EQSET_OBJECTS("Energy 2", my_app::EquationSet_Energy,
				 EquationSet_Energy)
      
      if (!found && m_throw_on_failure) {
	std::string msg = "Error - the \"Equation Set\" with \"Type\"= \"" + params->get<std::string>("Type") +
	  "\" is not a valid equation set identifier. Please supply the correct factory.\n";
	TEUCHOS_TEST_FOR_EXCEPTION(!found && m_throw_on_failure, std::logic_error, msg);
      }
      
      if (!found)
	return Teuchos::null;

      return eq_set;
    }
    
  };

}

