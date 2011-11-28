#include "Panzer_config.hpp"

#include "Panzer_ResponseFunctional_Aggregator.hpp"

namespace panzer {

template < >
Teuchos::RCP<ResponseAggregatorBase<panzer::Traits> > ResponseFunctional_Aggregator_Builder::
build<panzer::Traits::Residual>() const
{ 
   Teuchos::RCP<ResponseAggregatorBase<panzer::Traits> > respAgg =
      Teuchos::rcp(new ResponseFunctional_Aggregator<panzer::Traits::Residual,panzer::Traits>); 
   respAgg->setLinearObjFactory(getLinearObjFactory());
   respAgg->setGlobalIndexer(getGlobalIndexer());
   return respAgg;
}

}

#ifdef PANZER_EXPLICIT_TEMPLATE_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Panzer_ResponseFunctional_AggregatorT.hpp"

template class panzer::ResponseFunctional_Aggregator<panzer::Traits::Residual,panzer::Traits>;

#endif
