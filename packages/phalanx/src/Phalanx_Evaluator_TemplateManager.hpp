#ifndef PHX_FIELDEVALUATOR_TEMPLATE_MANAGER_HPP
#define PHX_FIELDEVALUATOR_TEMPLATE_MANAGER_HPP

#include "Phalanx_TemplateManager.hpp"
#include "Phalanx_Evaluator_Derived.hpp"

#include "boost/mpl/placeholders.hpp"
using namespace boost::mpl::placeholders;

namespace PHX {
  
  template<typename Traits>
  class Evaluator_TemplateManager : 
    public PHX::TemplateManager<typename Traits::EvalTypes,
				PHX::EvaluatorBase<Traits>,
				PHX::EvaluatorDerived<_,Traits> > {
    
  public:
    
    Evaluator_TemplateManager() {}
    
    ~Evaluator_TemplateManager() {}
    
  };

} 

#endif 
