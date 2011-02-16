#ifndef PANZER_STK_RYTHMOS_OBSERVER_FACTORY_HPP
#define PANZER_STK_RYTHMOS_OBSERVER_FACTORY_HPP

#include "Rythmos_IntegrationObserverBase.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"

#include "Panzer_STK_Utilities.hpp"

namespace panzer_stk {

  class RythmosObserverFactory_Epetra {

  public:
    
    virtual Teuchos::RCP<Rythmos::IntegrationObserverBase<double> >
    buildRythmosObserver(const Teuchos::RCP<panzer_stk::STK_Interface>& mesh,
			 const RCP<panzer::UniqueGlobalIndexer<int,int> >& dof_manager,
			 const Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> >& lof) const = 0;
  };

}

#endif
