// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef USER_APP_CLOSURE_MODEL_FACTORY_T_HPP
#define USER_APP_CLOSURE_MODEL_FACTORY_T_HPP

#include <iostream>
#include <sstream>
#include <typeinfo>
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Integrator_Scalar.hpp"

#include "Phalanx_FieldTag_Tag.hpp"

#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_TypeNameTraits.hpp"

// User application evaluators for this factory
#include "user_app_ConstantModel.hpp"
#include "user_app_TSquaredModel.hpp"

#include "Panzer_Parameter.hpp"
#include "Panzer_GlobalStatistics.hpp"
#include "Panzer_CoordinatesEvaluator.hpp"
#include "Panzer_Constant.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_DOF.hpp"
#include "Panzer_FieldSpy.hpp"
#include "Panzer_Product.hpp"
#include "Panzer_String_Utilities.hpp"

// ********************************************************************
// ********************************************************************
template<typename EvalT>
Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > 
user_app::MyModelFactory<EvalT>::
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

        if (plist.get<std::string>("Type") == "Parameter") {
          TEUCHOS_ASSERT(!plist.isParameter("Value"));
          // Defaults for backward compatibility
          std::string parameter_name = key;
          std::string field_name = key;
          if (plist.isType<std::string>("Parameter Name"))
            parameter_name = plist.get<std::string>("Parameter Name");
          if (plist.isType<std::string>("Field Name"))
            field_name = plist.get<std::string>("Field Name");

          { // at IP
            RCP< Evaluator<panzer::Traits> > e =
                rcp(new panzer::Parameter<EvalT,panzer::Traits>(parameter_name,field_name,ir->dl_scalar,*global_data->pl));
            evaluators->push_back(e);
          }

          for (std::vector<Teuchos::RCP<const panzer::PureBasis> >::const_iterator basis_itr = bases.begin();
              basis_itr != bases.end(); ++basis_itr) { // at BASIS
            Teuchos::RCP<const panzer::BasisIRLayout> basis = basisIRLayout(*basis_itr,*ir);
            RCP< Evaluator<panzer::Traits> > e =
                rcp(new panzer::Parameter<EvalT,panzer::Traits>(parameter_name,field_name,basis->functional,*global_data->pl));
            evaluators->push_back(e);
          }

          found = true;

          continue;
        }
        else if (plist.get<std::string>("Type") == "Distributed Parameter") {
          // sanity check
          TEUCHOS_ASSERT(distr_param_lof!=Teuchos::null);

          // build a nodal basis
          Teuchos::RCP<const panzer::PureBasis> nodal_basis
          = Teuchos::rcp(new panzer::PureBasis("HGrad",1,bases[0]->numCells(),
              bases[0]->getCellTopology()));

          {
            Teuchos::RCP<std::vector<std::string> > dof_names = Teuchos::rcp(new std::vector<std::string>);
            dof_names->push_back(key);

            ParameterList p("Gather");
            p.set("Basis", nodal_basis);
            p.set("DOF Names", dof_names);
            p.set("Indexer Names", dof_names);
            p.set("Sensitivities Name", key);
            p.set("First Sensitivities Available", true);
            p.set("Gather Seed Index", 0);
            p.set("Global Data Key", key);

            RCP< PHX::Evaluator<panzer::Traits> > e = distr_param_lof->buildGatherDomain<EvalT>(p);

            evaluators->push_back(e);
          }

          {
            ParameterList p;
            p.set("Name", key);
            p.set("Basis", basisIRLayout(nodal_basis,*ir));
            p.set("IR", ir);

            RCP< PHX::Evaluator<panzer::Traits> > e = rcp(new panzer::DOF<EvalT,panzer::Traits>(p));

            evaluators->push_back(e);
          }

          found = true;
        }
        else if (plist.get<std::string>("Type") == "Product") {
          RCP<std::vector<std::string> > valuesNames = rcp(new std::vector<std::string>); 
          // Tokenize a string, put tokens in a vector
          panzer::StringTokenizer(*valuesNames,plist.get<std::string>("Term Names"));
          double scaling = 1.0; 
          if(plist.isType<double>("Scaling")) 
            scaling = plist.get<double>("Scaling");
  
          {
            Teuchos::ParameterList input;
            input.set("Scaling", scaling);
            input.set("Product Name", key);
            input.set("Values Names", valuesNames);
            input.set("Data Layout", ir->dl_scalar);
  
            RCP< panzer::Product<EvalT,panzer::Traits> > e = 
              rcp(new panzer::Product<EvalT,panzer::Traits>(input));
            evaluators->push_back(e);
          }
  
          found = true;
        }
        else if (plist.get<std::string>("Type") == "T Squared Source") {
          // Required parameters
          plist.get<double>("Multiplier");
          plist.get<std::string>("Source Name");
          plist.get<std::string>("Target Name");
          {
            Teuchos::ParameterList input = plist;
            input.set("Data Layout", ir->dl_scalar);
            RCP< user_app::TSquaredModel<EvalT,panzer::Traits> > e = 
              rcp(new user_app::TSquaredModel<EvalT,panzer::Traits>(input));
            evaluators->push_back(e);
          }
          found = true;
        }
      }
      else if (plist.isType<double>("Value")) {
        { // at IP
          ParameterList input;
          input.set("Name", key);
          input.set("Value", plist.get<double>("Value"));
          input.set("Data Layout", ir->dl_scalar);
          RCP< Evaluator<panzer::Traits> > e =
              rcp(new panzer::Constant<EvalT,panzer::Traits>(input));
          evaluators->push_back(e);
        }
        // at BASIS
        for (std::vector<Teuchos::RCP<const panzer::PureBasis> >::const_iterator basis_itr = bases.begin();
            basis_itr != bases.end(); ++basis_itr) {
          ParameterList input;
          input.set("Name", key);
          input.set("Value", plist.get<double>("Value"));
          Teuchos::RCP<const panzer::BasisIRLayout> basis = basisIRLayout(*basis_itr,*ir);
          input.set("Data Layout", basis->functional);
          RCP< Evaluator<panzer::Traits> > e =
              rcp(new panzer::Constant<EvalT,panzer::Traits>(input));
          evaluators->push_back(e);
        }
        found = true;
      }

    if (plist.isType<std::string>("Value")) {

      const std::string value = plist.get<std::string>("Value");

      if (key == "Global Statistics") {
        if (typeid(EvalT) == typeid(panzer::Traits::Residual)) {
          ParameterList input;
          input.set("Comm", user_data.get<Teuchos::RCP<const Teuchos::Comm<int> > >("Comm"));
          input.set("Names", value);
          input.set("IR", ir);
          input.set("Global Data", global_data);
          RCP< panzer::GlobalStatistics<EvalT,panzer::Traits> > e =
              rcp(new panzer::GlobalStatistics<EvalT,panzer::Traits>(input));
          evaluators->push_back(e);

          // Require certain fields be evaluated
          fm.template requireField<EvalT>(e->getRequiredFieldTag());
        }
        found = true;
      }
      else if(value=="Field Spy") {
        const std::string & source = plist.get<std::string>("Source Field");

        RCP<panzer::FieldSpy<EvalT,panzer::Traits> > e =
            rcp(new panzer::FieldSpy<EvalT,panzer::Traits>(source,ir->dl_scalar));
        evaluators->push_back(e);

        fm.template requireField<EvalT>(e->getRequiredFieldTag());

        found = true;
      }
      else if(value=="Field Spy Basis") {
        const std::string & source = plist.get<std::string>("Source Field");

        RCP<panzer::FieldSpy<EvalT,panzer::Traits> > e =
            rcp(new panzer::FieldSpy<EvalT,panzer::Traits>(source,bases[0]->functional));
        evaluators->push_back(e);

        fm.template requireField<EvalT>(e->getRequiredFieldTag());

        found = true;
      }

    }

    if (key == "Volume Integral") {

      {
        ParameterList input;
        input.set("Name", "Unit Value");
        input.set("Value", 1.0);
        input.set("Data Layout", ir->dl_scalar);
        RCP< Evaluator<panzer::Traits> > e =
            rcp(new panzer::Constant<EvalT,panzer::Traits>(input));
        evaluators->push_back(e);
      }

      {
        ParameterList input;
        input.set("Integral Name", "Volume_Integral");
        input.set("Integrand Name", "Unit Value");
        input.set("IR", ir);

        RCP< Evaluator<panzer::Traits> > e =
            rcp(new panzer::Integrator_Scalar<EvalT,panzer::Traits>(input));
        evaluators->push_back(e);
      }

      found = true;
    }

    if (key == "Coordinates") {
      std::string dim_str[3] = {"X","Y","Z"};
      panzer::CellData cell_data(ir->workset_size,ir->topology);
      panzer::PureBasis basis("HGrad",1,cell_data);

      for(int i=0;i<basis.dimension();i++) {
        ParameterList input;
        input.set("Field Name", "COORD"+dim_str[i]);
        input.set("Data Layout", basis.functional);
        input.set("Dimension", i);

        RCP< Evaluator<panzer::Traits> > e = 
            rcp(new panzer::CoordinatesEvaluator<EvalT,panzer::Traits>(input));
        evaluators->push_back(e);
      }

      found = true;
    }


    if (!found) {
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
