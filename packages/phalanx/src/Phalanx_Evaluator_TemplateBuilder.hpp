#ifndef PHX_FIELDEVALUATOR_TEMPLATE_BUILDER_HPP
#define PHX_FIELDEVALUATOR_TEMPLATE_BUILDER_HPP

#include <string>
#include "boost/mpl/apply.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

namespace PHX {

  template<typename Traits, typename ObjectT>
  class Evaluator_TemplateBuilder {
    
    Teuchos::RCP<Teuchos::ParameterList> p;

  public:
    
    Evaluator_TemplateBuilder(const Teuchos::RCP<Teuchos::ParameterList>& param) :
      p(param) {}

    template <typename ScalarT>
    Teuchos::RCP< PHX::EvaluatorBase<Traits> > build() const {
      typedef typename boost::mpl::apply<ObjectT,ScalarT>::type type;
      return Teuchos::rcp( static_cast< PHX::EvaluatorBase<Traits>* > (new type(*p)) );
    }
    
  };
  
}

#endif 
