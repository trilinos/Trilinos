#ifndef PHX_FIELDEVALUATOR_TEMPLATE_BUILDER_HPP
#define PHX_FIELDEVALUATOR_TEMPLATE_BUILDER_HPP

#include <string>
#include "boost/mpl/apply.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

namespace PHX {

  template<typename Traits, typename ObjectT>
  class FieldEvaluator_TemplateBuilder {
    
    Teuchos::RCP<Teuchos::ParameterList> p;

  public:
    
    FieldEvaluator_TemplateBuilder(const Teuchos::RCP<Teuchos::ParameterList>& param) :
      p(param) {}

    template <typename ScalarT>
    Teuchos::RCP< PHX::FieldEvaluatorBase<Traits> > build() const {
      typedef typename boost::mpl::apply<ObjectT,ScalarT>::type type;
      return Teuchos::rcp( static_cast< PHX::FieldEvaluatorBase<Traits>* > (new type(*p)) );
    }
    
  };
  
}

#endif 
