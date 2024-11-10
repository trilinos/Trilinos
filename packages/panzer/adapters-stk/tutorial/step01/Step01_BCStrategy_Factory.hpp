// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Step01_BCStrategy_Factory_hpp__
#define __Step01_BCStrategy_Factory_hpp__

#include "Teuchos_RCP.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_BCStrategy_TemplateManager.hpp"
#include "Panzer_BCStrategy_Factory.hpp"
#include "Panzer_BCStrategy_Factory_Defines.hpp"
#include "Panzer_GlobalData.hpp"

namespace user_app {
  
class BCStrategyFactory : public panzer::BCStrategyFactory {
public:

  Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> >
  buildBCStrategy(const panzer::BC& /* bc */, const Teuchos::RCP<panzer::GlobalData>& /* global_data */) const
  {
    Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> > bcs_tm = 
        Teuchos::rcp(new panzer::BCStrategy_TemplateManager<panzer::Traits>);
      
    return bcs_tm;
  }

};
  
}

#endif
