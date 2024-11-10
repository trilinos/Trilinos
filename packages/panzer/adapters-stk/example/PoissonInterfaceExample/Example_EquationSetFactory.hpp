// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __PoissonExample_EquationSetFactory_hpp__
#define __PoissonExample_EquationSetFactory_hpp__

#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_EquationSet_Factory_Defines.hpp"
#include "Panzer_CellData.hpp"

// Add my equation sets here
#include "Example_PoissonEquationSet.hpp"

namespace Example {

// A macro that defines a class to make construction of the equation sets easier
//   - The equation set is constructed over a list of automatic differention types
PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(PoissonEquationSet, PoissonEquationSet)

// A user written factory that creates each equation set.  The key member here
// is buildEquationSet
class EquationSetFactory : public panzer::EquationSetFactory {
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
         
      bool found = false; // this is used by PANZER_BUILD_EQSET_OBJECTS
         
      // macro checks if(ies.name=="Poisson") then an EquationSet_Energy object is constructed
      PANZER_BUILD_EQSET_OBJECTS("Poisson", PoissonEquationSet)
         
      // make sure your equation set has been found
      if(!found) {
	std::string msg = "Error - the \"Equation Set\" called \"" + params->get<std::string>("Type") +
                           "\" is not a valid equation set identifier. Please supply the correct factory.\n";
         TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
      }
         
      return eq_set;
   }
    
};

}

#endif
