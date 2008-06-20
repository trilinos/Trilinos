#ifndef PHX_FIELDEVALUATOR_TEMPLATE_MANAGER_HPP
#define PHX_FIELDEVALUATOR_TEMPLATE_MANAGER_HPP

#include "Phalanx_TemplateManager.hpp"
#include "Phalanx_FieldEvaluator_Derived.hpp"

#include "boost/mpl/placeholders.hpp"
using namespace boost::mpl::placeholders;

namespace PHX {
  
  template<typename Traits>
  class FieldEvaluator_TemplateManager : 
    public PHX::TemplateManager<typename Traits::ScalarTypes,
				PHX::FieldEvaluatorBase<Traits>,
				PHX::FieldEvaluatorDerived<_,Traits> > {
    
  public:
    
    FieldEvaluator_TemplateManager() {}
    
    ~FieldEvaluator_TemplateManager() {}
    
  };

} 

#endif 
