#include "Panzer_config.hpp"
#include "Panzer_STK_config.hpp"

#ifdef HAVE_PANZER_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Panzer_STK_ResponseAggregator_SolutionWriter_decl.hpp"
#include "Panzer_STK_ResponseAggregator_SolutionWriter_impl.hpp"

template class panzer_stk::ResponseAggregator_SolutionWriter<panzer::Traits::Residual,panzer::Traits>;
#ifdef HAVE_STOKHOS
template class panzer_stk::ResponseAggregator_SolutionWriter<panzer::Traits::SGResidual,panzer::Traits>;
#endif

#endif
