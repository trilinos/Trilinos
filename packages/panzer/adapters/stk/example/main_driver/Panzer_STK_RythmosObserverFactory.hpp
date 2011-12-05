#ifndef PANZER_STK_RYTHMOS_OBSERVER_FACTORY_HPP
#define PANZER_STK_RYTHMOS_OBSERVER_FACTORY_HPP

#include "Rythmos_IntegrationObserverBase.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_LinearObjFactory.hpp"

#include "Panzer_STK_Utilities.hpp"

namespace panzer_stk {

  class RythmosObserverFactory {

  public:

    virtual ~RythmosObserverFactory() {}

    virtual Teuchos::RCP<Rythmos::IntegrationObserverBase<double> >
    buildRythmosObserver(const Teuchos::RCP<panzer_stk::STK_Interface>& mesh,
			 const Teuchos::RCP<panzer::UniqueGlobalIndexer<int,int> >& dof_manager,
			 const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> >& lof) const = 0;
  };

}

#endif
