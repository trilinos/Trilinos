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
#include "Phalanx_TypeList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Phalanx_Evaluator_TemplateManager.hpp"
#include "Phalanx_Evaluator_TemplateBuilder.hpp"

namespace PHX {

  //! Unary Function Object (UFO) - helper class required for PHX::mpl::for_each_index<>.
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
    
    template<std::size_t I> void operator()() const
    {
      if (object_type == static_cast<int>(I)) {
 	typedef PHX::mpl::at_t<typename FactoryTraits::EvaluatorTypes, I> type;
	PHX::Evaluator_TemplateBuilder<Traits, type> builder(params);
	tm->buildObjects(builder);
	found_object = true;
      }
    }
    
  };
  
} 

#endif 
