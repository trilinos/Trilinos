#ifndef PHX_SCALAR_CONTAINER_TEMPLATE_BUILDER_HPP
#define PHX_SCALAR_CONTAINER_TEMPLATE_BUILDER_HPP

#include "Phalanx_EvaluationContainer.hpp"

namespace PHX {

  template <typename Traits>
  class ScalarContainer_TemplateBuilder {

  public:
    
    ScalarContainer_TemplateBuilder() {}

    template <typename ScalarT>
    Teuchos::RCP<PHX::ScalarContainerBase<Traits> > build() const {
      return Teuchos::rcp( new ScalarContainer<ScalarT, Traits>);
    }
    
  };
  
}

#endif 
