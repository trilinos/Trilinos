#include "Panzer_config.hpp"

#include "Panzer_ResponseAggregator_IPCoordinates.hpp"

namespace panzer {

template < >
Teuchos::RCP<ResponseAggregatorBase<panzer::Traits> > ResponseAggregator_IPCoordinates_Builder::
build<panzer::Traits::Residual>() const
{ 
   Teuchos::RCP<ResponseAggregatorBase<panzer::Traits> > respAgg =
      Teuchos::rcp(new ResponseAggregator_IPCoordinates<panzer::Traits::Residual,panzer::Traits>); 
   respAgg->setLinearObjFactory(getLinearObjFactory());
   respAgg->setGlobalIndexer(getGlobalIndexer());
   return respAgg;
}

}

#ifdef PANZER_EXPLICIT_TEMPLATE_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Panzer_ResponseAggregator_IPCoordinatesT.hpp"

template class panzer::ResponseAggregator_IPCoordinates<panzer::Traits::Residual,panzer::Traits>;

#endif
