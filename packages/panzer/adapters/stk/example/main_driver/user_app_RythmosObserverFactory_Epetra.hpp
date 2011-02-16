#ifndef USER_APP_RYTHMOS_OBSERVER_FACTORY_HPP
#define USER_APP_RYTHMOS_OBSERVER_FACTORY_HPP

#include "Panzer_STK_RythmosObserverFactory_Epetra.hpp"
#include "user_app_RythmosObserver_Epetra.hpp"

namespace user_app {

  class RythmosObserverFactory_Epetra : public panzer_stk::RythmosObserverFactory_Epetra {

  public:
    
    Teuchos::RCP<Rythmos::IntegrationObserverBase<double> >
    buildRythmosObserver(const Teuchos::RCP<panzer_stk::STK_Interface>& mesh,
			 const RCP<panzer::UniqueGlobalIndexer<int,int> >& dof_manager,
			 const Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> >& lof) const
    {
      Teuchos::RCP<user_app::RythmosObserver_Epetra> observer = Teuchos::rcp(new user_app::RythmosObserver_Epetra(mesh,dof_manager,lof));
      return observer;
    }

  };

}

#endif
