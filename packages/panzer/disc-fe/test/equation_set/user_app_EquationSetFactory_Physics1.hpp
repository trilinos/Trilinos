// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_EquationSet_Factory_Defines.hpp"
#include "Panzer_CellData.hpp"

// Add my equation sets here
#include "user_app_EquationSet_Energy.hpp"

namespace user_app {

  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(user_app::EquationSet_Energy,
					EquationSet_Energy)

  class MyFactory1 : public panzer::EquationSetFactory {

    bool m_throw_on_failure;

  public:

    MyFactory1(bool throw_on_failure) : m_throw_on_failure(throw_on_failure) {}

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
      
      PANZER_BUILD_EQSET_OBJECTS("Energy 1", EquationSet_Energy)
      
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

