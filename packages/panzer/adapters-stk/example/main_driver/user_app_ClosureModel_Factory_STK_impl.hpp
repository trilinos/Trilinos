// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef USER_APP_CLOSURE_MODEL_FACTORY_STK_T_HPP
#define USER_APP_CLOSURE_MODEL_FACTORY_STK_T_HPP

#include <iostream>
#include <sstream>
#include <typeinfo>
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Integrator_Scalar.hpp"
#include "Panzer_String_Utilities.hpp"

#include "Phalanx_FieldTag_Tag.hpp"

#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

// Evaluators for this factory
#include "Panzer_STK_GatherExodusCellDataToIP.hpp"

// ********************************************************************
// ********************************************************************
template<typename EvalT>
Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > >
user_app::MyModelFactorySTK<EvalT>::
buildClosureModels(const std::string& model_id,
                   const Teuchos::ParameterList& models,
                   const panzer::FieldLayoutLibrary& fl,
                   const Teuchos::RCP<panzer::IntegrationRule>& ir,
                   const Teuchos::ParameterList& /* default_params */,
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

  if (!models.isSublist(model_id)) {
    models.print(std::cout);
    std::stringstream msg;
    msg << "Falied to find requested model, \"" << model_id
        << "\", for equation set:\n" << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(!models.isSublist(model_id), std::logic_error, msg.str());
  }

  std::vector<Teuchos::RCP<const panzer::PureBasis> > bases;
  fl.uniqueBases(bases);

  const ParameterList& my_models = models.sublist(model_id);

  for (ParameterList::ConstIterator model_it = my_models.begin();
       model_it != my_models.end(); ++model_it) {

    bool found = false;

    const std::string key = model_it->first;
    const Teuchos::ParameterEntry& entry = model_it->second;
    const ParameterList& plist = Teuchos::getValue<Teuchos::ParameterList>(entry);

    if (plist.isType<std::string>("Type")) {

      if (plist.get<std::string>("Type") == "Exodus Cell Data") {

        Teuchos::ParameterList validParams("Valid Params");
        auto plist_copy = plist;
        {
          validParams.set("Type","Exodus Cell Data");
          validParams.set("Field Names","");
          validParams.set("Exodus Names","");
          plist_copy.validateParametersAndSetDefaults(validParams);
        }

        RCP<std::vector<std::string> > fieldNames = Teuchos::make_rcp<std::vector<std::string>>();
        panzer::StringTokenizer(*fieldNames,plist_copy.get<std::string>("Field Names"));
        TEUCHOS_ASSERT(fieldNames->size() > 0);

        RCP<std::vector<std::string> > exodusNames = Teuchos::make_rcp<std::vector<std::string>>();
        panzer::StringTokenizer(*exodusNames,plist_copy.get<std::string>("Exodus Names"));
        TEUCHOS_ASSERT(exodusNames->size() > 0);

        auto e = Teuchos::make_rcp<panzer_stk::GatherExodusCellDataToIP<EvalT,panzer::Traits>>
                   (this->getMesh(),
                    *fieldNames,
                    *exodusNames,
                    ir);

          evaluators->push_back(e);
        found = true;
      }
    }

    if (!found && this->m_throw_if_model_not_found) {
      std::stringstream msg;
      msg << "ClosureModelFactory failed to build evaluator for key \"" << key
          << "\"\nin model \"" << model_id
          << "\".  Please correct the type or add support to the \nfactory." <<std::endl;
      TEUCHOS_TEST_FOR_EXCEPTION(!found, std::logic_error, msg.str());
    }
  }

  return evaluators;
}

#endif
