// @HEADER
// ************************************************************************
// 
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
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
    
    template<typename T> void operator()(T t)
    {
      if (object_type == t) {
	typedef typename Sacado::mpl::at<typename FactoryTraits::EvaluatorTypes, T::value >::type type;
	PHX::Evaluator_TemplateBuilder<Traits, type> builder(params);
	tm->buildObjects(builder);
	found_object = true;
      }
    }
    
  };
  
} 

#endif 
