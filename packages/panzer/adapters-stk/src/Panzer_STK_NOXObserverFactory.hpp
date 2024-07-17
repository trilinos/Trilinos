// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_STK_NOX_OBSERVER_FACTORY_HPP
#define PANZER_STK_NOX_OBSERVER_FACTORY_HPP

#include "NOX_Abstract_PrePostOperator.H"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"

#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_LinearObjFactory.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_Utilities.hpp"

namespace panzer_stk {

  class NOXObserverFactory {

  public:
    
    virtual ~NOXObserverFactory() {}

    virtual Teuchos::RCP<NOX::Abstract::PrePostOperator>
    buildNOXObserver(const Teuchos::RCP<panzer_stk::STK_Interface>& mesh,
		     const Teuchos::RCP<const panzer::GlobalIndexer>& dof_manager,
		     const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> >& lof) const = 0;
  };

}

#endif
