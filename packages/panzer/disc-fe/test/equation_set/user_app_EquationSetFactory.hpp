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
#include "user_app_EquationSet_MeshCoords.hpp"

namespace user_app {

  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(user_app::EquationSet_Energy,
					EquationSet_Energy)
  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(user_app::EquationSet_MeshCoords,
					EquationSet_MeshCoords)

  class MyFactory : public panzer::EquationSetFactory {

  public:

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
      
      PANZER_BUILD_EQSET_OBJECTS("Energy",     EquationSet_Energy)
      PANZER_BUILD_EQSET_OBJECTS("MeshCoords", EquationSet_MeshCoords)
      
      if (!found) {
	std::string msg = "Error - the \"Equation Set\" with \"Type\" = \"" + params->get<std::string>("Type") +
	  "\" is not a valid equation set identifier. Please supply the correct factory.\n";
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
      }
      
      return eq_set;
    }
    
  };

}
