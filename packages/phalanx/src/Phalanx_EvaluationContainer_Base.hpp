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

#ifndef PHX_EVALUATION_CONTAINER_BASE_HPP
#define PHX_EVALUATION_CONTAINER_BASE_HPP

#include <cstddef>
#include <string>
#include <map>
#include "Phalanx_Evaluator_Manager.hpp"

namespace PHX {

  template<typename Traits> class FieldManager;

  template<typename Traits>
  class EvaluationContainerBase {

  public:

    EvaluationContainerBase();

    virtual ~EvaluationContainerBase();

    virtual void requireField(const PHX::FieldTag& v);

    virtual void 
    registerEvaluator(const Teuchos::RCP<PHX::Evaluator<Traits> >& p);

    virtual void postRegistrationSetup(typename Traits::SetupData d,
				       PHX::FieldManager<Traits>& vm) = 0;

    virtual void evaluateFields(typename Traits::EvalData d) = 0;

    virtual void preEvaluate(typename Traits::PreEvalData d) = 0;

    virtual void postEvaluate(typename Traits::PostEvalData d) = 0;

    virtual void print(std::ostream& os) const = 0;
    
  protected:
    
    PHX::EvaluatorManager<Traits> vp_manager_;

  };

  template<typename Traits>
  std::ostream& operator<<(std::ostream& os, 
			   const PHX::EvaluationContainerBase<Traits>& sc);
  
}

#include "Phalanx_EvaluationContainer_Base_Def.hpp"

#endif 
