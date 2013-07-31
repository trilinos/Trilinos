// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef USER_APP_STK_CLOSURE_MODEL_FACTORY_T_HPP
#define USER_APP_STK_CLOSURE_MODEL_FACTORY_T_HPP

#include <iostream>
#include <sstream>
#include "Panzer_Constant.hpp"

#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_TypeNameTraits.hpp"

// User application evaluators for this factory
#include "user_app_ConstantModel.hpp"
#include "TestEvaluators.hpp"

// ********************************************************************
// ********************************************************************
template<typename EvalT>
Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > 
user_app::STKModelFactory<EvalT>::
buildClosureModels(const std::string& model_id,
		   const Teuchos::ParameterList& models,  
		   const panzer::FieldLayoutLibrary& fl,
		   const Teuchos::RCP<panzer::IntegrationRule>& ir,
		   const Teuchos::ParameterList& default_params,
		   const Teuchos::ParameterList& user_data,
		   const Teuchos::RCP<panzer::GlobalData>& global_data,
		   PHX::FieldManager<panzer::Traits>& fm) const
{

  using std::string;
  using std::vector;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using PHX::Evaluator;

  RCP< vector< RCP<Evaluator<panzer::Traits> > > > evaluators = 
    rcp(new vector< RCP<Evaluator<panzer::Traits> > > );

  const ParameterList& my_models = models.sublist(model_id);

  for (ParameterList::ConstIterator model_it = my_models.begin(); 
       model_it != my_models.end(); ++model_it) {
    const std::string key = model_it->first;
    const Teuchos::ParameterEntry& entry = model_it->second;
    const ParameterList& plist = Teuchos::getValue<Teuchos::ParameterList>(entry);
    if (plist.isType<double>("Value")) {
      Teuchos::ParameterList input;
      input.set("Name", key);
      input.set("Value", plist.get<double>("Value"));
      input.set("Data Layout", ir->dl_scalar);
      RCP< PHX::Evaluator<panzer::Traits> > e = 
        rcp(new panzer::Constant<EvalT,panzer::Traits>(input));
      evaluators->push_back(e);
    }
  }

  {
     RCP<PHX::Evaluator<panzer::Traits> > e
        = rcp(new panzer::TestEvaluator<EvalT,panzer::Traits>(user_data));

     evaluators->push_back(e);
  }


  return evaluators;
}

#endif
