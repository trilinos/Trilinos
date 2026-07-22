// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
