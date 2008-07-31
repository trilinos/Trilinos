#ifndef PHX_SCALAR_CONTAINER_TEMPLATE_MANAGER_HPP
#define PHX_SCALAR_CONTAINER_TEMPLATE_MANAGER_HPP

#include "Phalanx_TemplateManager.hpp"
#include "Phalanx_EvaluationContainer.hpp"

#include "boost/mpl/placeholders.hpp"
using namespace boost::mpl::placeholders;

namespace PHX {

  template<typename Traits>
  class ScalarContainer_TemplateManager : 
    public PHX::TemplateManager<typename Traits::ScalarTypes,
				PHX::ScalarContainerBase<Traits>,
				PHX::ScalarContainer<_,Traits> > {

  public:

    ScalarContainer_TemplateManager() {}

    ~ScalarContainer_TemplateManager() {}

  };

} 

#endif 
