// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef USER_APP_TEMPUS_OBSERVER_FACTORY_HPP
#define USER_APP_TEMPUS_OBSERVER_FACTORY_HPP

#include "PanzerAdaptersSTK_config.hpp"

#ifdef PANZER_HAVE_TEMPUS

#include "Panzer_STK_TempusObserverFactory.hpp"
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_Traits.hpp"

#include "Tempus_IntegratorObserverComposite.hpp"

// Individual Observers
#include "user_app_TempusObserver_WriteToExodus.hpp"

namespace user_app {

  class TempusObserverFactory : public panzer_stk::TempusObserverFactory {

  public:
    TempusObserverFactory(const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > & stkIOResponseLibrary,
                           const Teuchos::RCP<panzer::WorksetContainer> wkstContainer)
       : stkIOResponseLibrary_(stkIOResponseLibrary)
       , wkstContainer_(wkstContainer)
    {}

    bool useNOXObserver() const { return false; }
    
    Teuchos::RCP<Tempus::IntegratorObserver<double> >
    buildTempusObserver(const Teuchos::RCP<panzer_stk::STK_Interface>& mesh,
			 const Teuchos::RCP<const panzer::GlobalIndexer> & dof_manager,
			 const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> >& lof) const
    {
      Teuchos::RCP<Tempus::IntegratorObserverComposite<double> > composite_observer =
        Teuchos::rcp(new Tempus::IntegratorObserverComposite<double>);

      {
        Teuchos::RCP<user_app::TempusObserver_WriteToExodus> observer 
            = Teuchos::rcp(new user_app::TempusObserver_WriteToExodus(mesh,dof_manager,lof,stkIOResponseLibrary_));
        composite_observer->addObserver(observer);
      }

      return composite_observer;
    }

  private:
    //! Store STK IO response library...be careful, it will be modified externally
    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > stkIOResponseLibrary_;

    Teuchos::RCP<panzer::WorksetContainer> wkstContainer_;
  };

}

#endif // PANZER_HAVE_TEMPUS
#endif
