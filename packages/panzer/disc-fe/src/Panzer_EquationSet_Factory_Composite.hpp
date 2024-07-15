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
#include "Teuchos_RCP.hpp"

namespace Teuchos {
  class ParameterList;
}

namespace panzer {

  class CellData;
  struct GlobalData;

  class EquationSet_FactoryComposite : public panzer::EquationSetFactory {

    std::vector<Teuchos::RCP<panzer::EquationSetFactory> > m_factories;

  public:
    
    EquationSet_FactoryComposite(const std::vector<Teuchos::RCP<panzer::EquationSetFactory> >& factories);

    Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> >
    buildEquationSet(const Teuchos::RCP<Teuchos::ParameterList>& input_plist,
		     const int& default_integration_rule,
		     const panzer::CellData& cell_data,
		     const Teuchos::RCP<panzer::GlobalData>& global_data,
		     const bool build_transient_support) const;
    
  };

}

