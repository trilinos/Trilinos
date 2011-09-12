#include "Panzer_ResponseFunctional_Aggregator.hpp"

namespace panzer {

template < >
Teuchos::RCP<ResponseAggregatorBase<panzer::Traits> > ResponseFunctional_Aggregator_Builder::
build<panzer::Traits::Residual>() const
{ 
   return Teuchos::rcp(new ResponseFunctional_Aggregator<panzer::Traits::Residual,panzer::Traits>); 
}

}
