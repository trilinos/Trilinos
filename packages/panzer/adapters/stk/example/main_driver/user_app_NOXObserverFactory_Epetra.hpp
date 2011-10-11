#ifndef USER_APP_NOX_OBSERVER_FACTORY_HPP
#define USER_APP_NOX_OBSERVER_FACTORY_HPP

#include "Panzer_STK_NOXObserverFactory.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "NOX_PrePostOperator_Vector.H"

// Individual Observers
#include "user_app_NOXObserver_EpetraToExodus.hpp"

namespace user_app {
  
  class NOXObserverFactory_Epetra : public panzer_stk::NOXObserverFactory {
    
  public:
    
    Teuchos::RCP<NOX::Abstract::PrePostOperator>
    buildNOXObserver(const Teuchos::RCP<panzer_stk::STK_Interface>& mesh,
		     const RCP<panzer::UniqueGlobalIndexer<int,int> >& dof_manager,
		     const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> >& lof,
                     const Teuchos::ParameterList & pl) const
    {
      Teuchos::RCP<NOX::PrePostOperatorVector> observer = 
	Teuchos::rcp(new NOX::PrePostOperatorVector);

      Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> > ep_lof
         = Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjFactory<panzer::Traits,int> >(lof);

      // Always register the exodus writer to output solution
      {
	Teuchos::RCP<NOX::Abstract::PrePostOperator> solution_writer = 
	  Teuchos::rcp(new user_app::NOXObserver_EpetraToExodus(mesh,dof_manager,ep_lof));
	observer->pushBack(solution_writer);
      }

      return observer;
    }

  };

}

#endif
