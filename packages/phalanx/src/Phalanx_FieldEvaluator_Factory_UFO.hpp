#ifndef PHX_FIELD_EVALUATOR_FACTORY_UFO_HPP
#define PHX_FIELD_EVALUATOR_FACTORY_UFO_HPP

#include <map>
#include <vector>
#include "Sacado_mpl_at.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Phalanx_FieldEvaluator_TemplateManager.hpp"
#include "Phalanx_FieldEvaluator_TemplateBuilder.hpp"

namespace PHX {

  //! Unary Function Object (UFO) - helper class required for mpl::for_each<>.
  template<typename Traits, typename FactoryTraits>
  struct UFO {
    
    int object_type;
    Teuchos::RCP<Teuchos::ParameterList> params;
    Teuchos::RCP< FieldEvaluator_TemplateManager<Traits> >& tm;
    bool& found_object;
    
    UFO(int v, const Teuchos::RCP<Teuchos::ParameterList>& p, 
	Teuchos::RCP< FieldEvaluator_TemplateManager<Traits> >& t,
	bool& found) 
      : object_type(v), params(p), tm(t), found_object(found) {}
    
    template<typename T> void operator()(T t)
    {
      if (object_type == t) {
	typedef typename Sacado::mpl::at<typename FactoryTraits::FieldEvaluatorTypes, T::value >::type type;
	PHX::FieldEvaluator_TemplateBuilder<Traits, type> builder(params);
	tm->buildObjects(builder);
	found_object = true;
      }
    }
    
  };
  
} 

#endif 
