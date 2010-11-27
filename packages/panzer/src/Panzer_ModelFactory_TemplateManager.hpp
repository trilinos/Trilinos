#ifndef PANZER_MODEL_FACTORY_TEMPLATE_MANAGER_HPP
#define PANZER_MODEL_FACTORY_TEMPLATE_MANAGER_HPP

#include "Phalanx_TemplateManager.hpp"

#include "Panzer_Traits.hpp"
#include "Panzer_Base.hpp"
#include "Panzer_ModelFactory.hpp"

#include "boost/mpl/placeholders.hpp"
using namespace boost::mpl::placeholders;

namespace panzer {

  template<typename Traits>
  class ModelFactory_TemplateManager : 
    public PHX::TemplateManager<typename Traits::EvalTypes,
				panzer::Base,
                                panzer::ModelFactory<_> > {

  public:

    ModelFactory_TemplateManager() {}

    ~ModelFactory_TemplateManager() {}

  };

} 

#endif 
