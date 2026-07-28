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

  /** \brief Abstract factory for building a NOX::Observer that has
    * access to the STK mesh, DOF manager, and linear object factory,
    * e.g. for writing solution snapshots to disk between nonlinear
    * iterations.
    */
  class NOXObserverFactory {

  public:

    virtual ~NOXObserverFactory() {}

    /** \brief Builds a NOX::Observer configured with the given mesh, DOF manager, and linear object factory.
      *
      * \param[in] mesh the STK mesh the observer may write solution data to.
      * \param[in] dof_manager maps solution vector entries to mesh fields.
      * \param[in] lof linear object factory used to interpret the solution vector passed to the observer.
      */
    virtual Teuchos::RCP<NOX::Abstract::PrePostOperator>
    buildNOXObserver(const Teuchos::RCP<panzer_stk::STK_Interface>& mesh,
		     const Teuchos::RCP<const panzer::GlobalIndexer>& dof_manager,
		     const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> >& lof) const = 0;
  };

}

#endif
