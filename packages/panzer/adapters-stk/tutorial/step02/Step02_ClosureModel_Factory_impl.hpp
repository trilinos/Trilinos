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

#ifndef __Step02_ClosureModel_Factory_impl_hpp__
#define __Step02_ClosureModel_Factory_impl_hpp__

#include <iostream>
#include <sstream>
#include <typeinfo>

#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Phalanx_FieldTag_Tag.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"

#include "Panzer_Constant.hpp"
#include "Panzer_Integrator_Scalar.hpp"

#include "Step02_LinearFunction.hpp"

// ********************************************************************
// ********************************************************************
template<typename EvalT>
Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > 
user_app::ClosureModelFactory<EvalT>::
buildClosureModels(const std::string& model_id,
    const Teuchos::ParameterList& models,
    const panzer::FieldLayoutLibrary& fl,
    const Teuchos::RCP<panzer::IntegrationRule>& ir,
    const Teuchos::ParameterList& default_params,
    const Teuchos::ParameterList& user_data,
    const Teuchos::RCP<panzer::GlobalData>& global_data,
    PHX::FieldManager<panzer::Traits>& fm) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  RCP<std::vector< RCP<PHX::Evaluator<panzer::Traits> > > > evaluators = 
      rcp(new std::vector< RCP<PHX::Evaluator<panzer::Traits> > > );

  if (!models.isSublist(model_id)) {
    models.print(std::cout);
    std::stringstream msg;
    msg << "Falied to find requested model, \"" << model_id 
        << "\", for equation set:\n" << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(!models.isSublist(model_id), std::logic_error, msg.str());
  }

  std::vector<Teuchos::RCP<const panzer::PureBasis> > bases;
  fl.uniqueBases(bases);

  // for a model id get the parameter list that defines all the models
  const ParameterList& my_models = models.sublist(model_id);

  // loop over all the models
  for (ParameterList::ConstIterator model_it = my_models.begin(); 
      model_it != my_models.end(); ++model_it) {

    bool found = false;

    // extract the model information, the "key" and the parameter list that
    // defines the model
    const std::string key = model_it->first;
    ParameterList input;
    const Teuchos::ParameterEntry& entry = model_it->second;
    const ParameterList& plist = Teuchos::getValue<Teuchos::ParameterList>(entry);

    // add in a constant evaluator, using the PL key as the name of the field
    if (plist.isType<double>("Value")) {
      // add constant evaluator for each Integration Point (IP)
      {
        input.set("Name", key);
        input.set("Value", plist.get<double>("Value"));
        input.set("Data Layout", ir->dl_scalar);
        RCP<PHX::Evaluator<panzer::Traits> > e =
              rcp(new panzer::Constant<EvalT,panzer::Traits>(input));
        evaluators->push_back(e);
      }

      // add constant evaluator for each basis
      for (std::vector<Teuchos::RCP<const panzer::PureBasis> >::const_iterator basis_itr = bases.begin();
        basis_itr != bases.end(); ++basis_itr) {
        input.set("Name", key);
        input.set("Value", plist.get<double>("Value"));
        Teuchos::RCP<const panzer::BasisIRLayout> basis = basisIRLayout(*basis_itr,*ir);
        input.set("Data Layout", basis->functional);
        RCP<PHX::Evaluator<panzer::Traits> > e =
            rcp(new panzer::Constant<EvalT,panzer::Traits>(input));
        evaluators->push_back(e);
      }
      found = true;
    }
    else if (plist.isType<std::string>("Type")) {
      std::string evaluator_type = plist.get<std::string>("Type");

      // a linear fucntion across the domain
      if(evaluator_type == "Linear Function" ) {
        double acoeff = plist.get<double>("ACoeff"); 
        double bcoeff = plist.get<double>("BCoeff"); 

        RCP<PHX::Evaluator<panzer::Traits> > e =
            rcp(new user_app::LinearFunction<EvalT,panzer::Traits>(key,acoeff,bcoeff,*ir));
        evaluators->push_back(e);

        found = true;
      }
    }

    // fail if you can't find one of the models
    if (!found) {
      std::stringstream msg;
      msg << "ClosureModelFactory failed to build evaluator for key \"" << key 
          << "\"\nin model \"" << model_id
          << "\".  Please correct the type or add support to the \nfactory." <<std::endl;
      TEUCHOS_TEST_FOR_EXCEPTION(!found, std::logic_error, msg.str());
    }

  } // end loop over models in parameter list

  return evaluators;
}

#endif
