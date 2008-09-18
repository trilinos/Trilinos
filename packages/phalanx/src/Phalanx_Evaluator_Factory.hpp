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

#ifndef PHX_FIELD_EVALUATOR_FACTORY_HPP
#define PHX_FIELD_EVALUATOR_FACTORY_HPP

#include <map>
#include <vector>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Phalanx_Evaluator_TemplateManager.hpp"

namespace PHX {

  template<typename Traits, typename FactoryTraits>
  class EvaluatorFactory {
    
  public:
    typename Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator_TemplateManager<Traits> > > > 
    buildEvaluators(const std::map<std::string, 
			 Teuchos::RCP<Teuchos::ParameterList> >& data);
    
  };


  /*! \brief Nonmember helper function for registering field evaluators for all scalar types that are built with template managers.

  \relates PHX::EvaluatorFactory

  */
  template<typename Traits>
  void registerEvaluators(const Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator_TemplateManager<Traits> > > >& t, PHX::FieldManager<Traits>& fm);

} 

#include "Phalanx_Evaluator_Factory_Def.hpp"

#endif 
