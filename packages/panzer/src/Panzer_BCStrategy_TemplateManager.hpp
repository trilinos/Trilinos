#ifndef PANZER_BCSTRATEGY_TEMPLATE_MANAGER_HPP
#define PANZER_BCSTRATEGY_TEMPLATE_MANAGER_HPP

#include "Phalanx_TemplateManager.hpp"

#include "Panzer_Traits.hpp"
#include "Panzer_BCStrategy_Base.hpp"
#include "Panzer_BCStrategy.hpp"

#include "boost/mpl/placeholders.hpp"
using namespace boost::mpl::placeholders;

namespace panzer {

  template<typename Traits>
    class BCStrategy_TemplateManager : 
    public PHX::TemplateManager<typename Traits::EvalTypes,
				panzer::BCStrategyBase,
                                panzer::BCStrategy<_> > {

  public:

    BCStrategy_TemplateManager() {}

    ~BCStrategy_TemplateManager() {}

  };

} 

#endif 
