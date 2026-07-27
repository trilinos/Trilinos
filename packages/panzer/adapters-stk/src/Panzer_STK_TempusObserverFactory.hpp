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

  /** \brief Abstract factory for building a Tempus integrator observer
    * that has access to the STK mesh, DOF manager, and linear object
    * factory, e.g. for writing solution snapshots to disk between time
    * steps.
    */
  class TempusObserverFactory {

  public:

    virtual ~TempusObserverFactory() {}

    /// \brief Returns true if a NOX observer should also be attached alongside the Tempus observer.
    virtual bool useNOXObserver() const = 0;

    /** \brief Builds a Tempus integrator observer configured with the given mesh, DOF manager, and linear object factory.
      *
      * \param[in] mesh the STK mesh the observer may write solution data to.
      * \param[in] dof_manager maps solution vector entries to mesh fields.
      * \param[in] lof linear object factory used to interpret the solution vector passed to the observer.
      */
    virtual Teuchos::RCP<Tempus::IntegratorObserver<double> >
    buildTempusObserver(const Teuchos::RCP<panzer_stk::STK_Interface>& mesh,
                        const Teuchos::RCP<const panzer::GlobalIndexer> & dof_manager,
                        const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> >& lof) const = 0;
  };

}

#endif
