#ifndef USER_APP_RYTHMOS_OBSERVER_FACTORY_HPP
#define USER_APP_RYTHMOS_OBSERVER_FACTORY_HPP

#include "Panzer_STK_RythmosObserverFactory.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"

// Individual Observers
#include "user_app_RythmosObserver_EpetraToExodus.hpp"

namespace user_app {

  class RythmosObserverFactory_Epetra : public panzer_stk::RythmosObserverFactory {

  public:
    bool useNOXObserver() const { return false; }
    
    Teuchos::RCP<Rythmos::IntegrationObserverBase<double> >
    buildRythmosObserver(const Teuchos::RCP<panzer_stk::STK_Interface>& mesh,
			 const Teuchos::RCP<panzer::UniqueGlobalIndexer<int,int> >& dof_manager,
			 const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> >& lof) const
    {
      Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> > ep_lof
         = Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjFactory<panzer::Traits,int> >(lof);

      Teuchos::RCP<user_app::RythmosObserver_EpetraToExodus> observer = Teuchos::rcp(new user_app::RythmosObserver_EpetraToExodus(mesh,dof_manager,ep_lof));
      return observer;
    }

  };

}

#endif
