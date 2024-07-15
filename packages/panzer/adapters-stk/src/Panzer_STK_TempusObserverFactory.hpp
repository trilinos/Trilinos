// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_STK_TEMPUS_OBSERVER_FACTORY_HPP
#define PANZER_STK_TEMPUS_OBSERVER_FACTORY_HPP

#include "Tempus_IntegratorObserver.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_LinearObjFactory.hpp"

#include "Panzer_STK_Utilities.hpp"

namespace panzer_stk {

  class TempusObserverFactory {

  public:

    virtual ~TempusObserverFactory() {}

    //! Use the NOX observer as well?
    virtual bool useNOXObserver() const = 0;

    virtual Teuchos::RCP<Tempus::IntegratorObserver<double> >
    buildTempusObserver(const Teuchos::RCP<panzer_stk::STK_Interface>& mesh,
                        const Teuchos::RCP<const panzer::GlobalIndexer> & dof_manager,
                        const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> >& lof) const = 0;
  };

}

#endif
