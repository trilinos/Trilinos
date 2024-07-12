// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_FIELD_EVALUATOR_FACTORY_UFO_HPP
#define PHX_FIELD_EVALUATOR_FACTORY_UFO_HPP

#include <map>
#include <vector>
#include "Sacado_mpl_at.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Phalanx_Evaluator_TemplateManager.hpp"
#include "Phalanx_Evaluator_TemplateBuilder.hpp"

namespace PHX {

  //! Unary Function Object (UFO) - helper class required for mpl::for_each<>.
  template<typename Traits, typename FactoryTraits>
  struct UFO {
    
    int object_type;
    Teuchos::RCP<Teuchos::ParameterList> params;
    Teuchos::RCP< Evaluator_TemplateManager<Traits> >& tm;
    bool& found_object;
    
    UFO(int v, const Teuchos::RCP<Teuchos::ParameterList>& p, 
	Teuchos::RCP< Evaluator_TemplateManager<Traits> >& t,
	bool& found) 
      : object_type(v), params(p), tm(t), found_object(found) {}
    
    template<typename T> void operator()(T t) const
    {
      if (object_type == t) {
 	typedef typename Sacado::mpl::at<typename FactoryTraits::EvaluatorTypes, T::value>::type type;
	PHX::Evaluator_TemplateBuilder<Traits, type> builder(params);
	tm->buildObjects(builder);
	found_object = true;
      }
    }
    
  };
  
} 

#endif 
