#ifndef PHX_SCALAR_CONTAINER_TEMPLATE_BUILDER_HPP
#define PHX_SCALAR_CONTAINER_TEMPLATE_BUILDER_HPP

#include "Phalanx_EvaluationContainer.hpp"

namespace PHX {

  template <typename Traits>
  class EvaluationContainer_TemplateBuilder {

  public:
    
    EvaluationContainer_TemplateBuilder() {}

    template <typename ScalarT>
    Teuchos::RCP<PHX::EvaluationContainerBase<Traits> > build() const {
      return Teuchos::rcp( new EvaluationContainer<ScalarT, Traits>);
    }
    
  };
  
}

#endif 
