// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_EquationSet_Factory_Composite.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_GlobalData.hpp"
#include <sstream>

namespace panzer {
  
  EquationSet_FactoryComposite::EquationSet_FactoryComposite(const std::vector<Teuchos::RCP<panzer::EquationSetFactory> >& factories) :
    m_factories(factories)
  { }
  
  Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> >
  EquationSet_FactoryComposite::buildEquationSet(const Teuchos::RCP<Teuchos::ParameterList>& input_plist,
						 const int& default_integration_rule,
						 const panzer::CellData& cell_data,
						 const Teuchos::RCP<panzer::GlobalData>& global_data,
						 const bool build_transient_support) const
  {
    Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > eq_set;
    
    
    for (std::vector<Teuchos::RCP<panzer::EquationSetFactory> >::const_iterator factory = m_factories.begin();
	 factory != m_factories.end(); ++factory) {
      eq_set = (*factory)->buildEquationSet(input_plist,default_integration_rule,cell_data,global_data,build_transient_support);

      if (nonnull(eq_set))
	break;
    }
    
    std::ostringstream os;
    os << "The equation set factory failed to build and equation set for the following input parameter list.  Please correct the input list:\n";
    input_plist->print(os);
    TEUCHOS_TEST_FOR_EXCEPTION(is_null(eq_set), std::logic_error,os.str());
    
    return eq_set;
  }
  
}

