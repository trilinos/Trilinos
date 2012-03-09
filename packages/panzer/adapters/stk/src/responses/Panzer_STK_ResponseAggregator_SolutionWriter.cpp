#include "Panzer_config.hpp"
#include "Panzer_STK_config.hpp"

#include "Panzer_STK_ResponseAggregator_SolutionWriter_decl.hpp"

namespace panzer_stk {

template < >
Teuchos::RCP<panzer::ResponseAggregatorBase<panzer::Traits> > 
ResponseAggregator_SolutionWriter_Builder::
build<panzer::Traits::Residual>() const
{ 
   Teuchos::RCP<panzer::ResponseAggregatorBase<panzer::Traits> > respAgg =
      Teuchos::rcp(new ResponseAggregator_SolutionWriter<panzer::Traits::Residual,panzer::Traits>(mesh_)); 
   respAgg->setLinearObjFactory(getLinearObjFactory());
   respAgg->setGlobalIndexer(getGlobalIndexer());
   return respAgg;
}

#ifdef HAVE_STOKHOS
template < >
Teuchos::RCP<panzer::ResponseAggregatorBase<panzer::Traits> > 
ResponseAggregator_SolutionWriter_Builder::
build<panzer::Traits::SGResidual>() const
{ 
   Teuchos::RCP<panzer::ResponseAggregatorBase<panzer::Traits> > respAgg =
      Teuchos::rcp(new ResponseAggregator_SolutionWriter<panzer::Traits::SGResidual,panzer::Traits>(mesh_)); 
   respAgg->setLinearObjFactory(getLinearObjFactory());
   respAgg->setGlobalIndexer(getGlobalIndexer());
   return respAgg;
}
#endif

}

#ifdef HAVE_PANZER_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Panzer_STK_ResponseAggregator_SolutionWriter_impl.hpp"

template class panzer_stk::ResponseAggregator_SolutionWriter<panzer::Traits::Residual,panzer::Traits>;
#ifdef HAVE_STOKHOS
template class panzer_stk::ResponseAggregator_SolutionWriter<panzer::Traits::SGResidual,panzer::Traits>;
#endif

#endif
