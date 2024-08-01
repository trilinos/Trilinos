// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Example_ClosureModelFactoryT_hpp__
#define __Example_ClosureModelFactoryT_hpp__

#include <iostream>
#include <sstream>
#include <typeinfo>
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Integrator_Scalar.hpp"
#include "Phalanx_FieldTag_Tag.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_TypeNameTraits.hpp"

#include "Panzer_Product.hpp"
#include "Panzer_DotProduct.hpp"
#include "Panzer_Sum.hpp"

#include "Example_SimpleSource.hpp"
#include "Example_SimpleSolution.hpp"

// ********************************************************************
// ********************************************************************
template<typename EvalT>
Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > 
Example::ModelFactory<EvalT>::
buildClosureModels(const std::string& model_id,
		   const Teuchos::ParameterList& models, 
		   const panzer::FieldLayoutLibrary& fl,
		   const Teuchos::RCP<panzer::IntegrationRule>& ir,
		   const Teuchos::ParameterList& /* default_params */,
		   const Teuchos::ParameterList& /* user_data */,
		   const Teuchos::RCP<panzer::GlobalData>& /* global_data */,
		   PHX::FieldManager<panzer::Traits>& /* fm */) const
{
  using std::string;
  using std::vector;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using PHX::Evaluator;

  RCP< vector< RCP<Evaluator<panzer::Traits> > > > evaluators = 
    rcp(new vector< RCP<Evaluator<panzer::Traits> > > );

  const bool curvilinear = models.get<bool>("Curvilinear");

  if (!models.isSublist(model_id)) {
    models.print(std::cout);
    std::stringstream msg;
    msg << "Failed to find requested model, \"" << model_id 
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
    ParameterList input;
    const Teuchos::ParameterEntry& entry = model_it->second;
    const ParameterList& plist = Teuchos::getValue<Teuchos::ParameterList>(entry);

    if (plist.isType<double>("Value")) {
      { // at IP
	input.set("Name", key);
	input.set("Value", plist.get<double>("Value"));
	input.set("Data Layout", ir->dl_scalar);
	RCP< Evaluator<panzer::Traits> > e = 
	  rcp(new panzer::Constant<EvalT,panzer::Traits>(input));
	evaluators->push_back(e);
      }
      
      for (std::vector<Teuchos::RCP<const panzer::PureBasis> >::const_iterator basis_itr = bases.begin();
	   basis_itr != bases.end(); ++basis_itr) { // at BASIS
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

    if (plist.isType<std::string>("Type")) {
      std::string type = plist.get<std::string>("Type");
      if(type=="SIMPLE SOURCE") {
        RCP< Evaluator<panzer::Traits> > e = 
          rcp(new Example::SimpleSource<EvalT,panzer::Traits>(key,*ir,curvilinear));
        evaluators->push_back(e);

        found = true;
      }
      else if(type=="TEMPERATURE_EXACT") {
        RCP< Evaluator<panzer::Traits> > e = 
          rcp(new Example::SimpleSolution<EvalT,panzer::Traits>(key,*ir,curvilinear));
        evaluators->push_back(e);

        found = true;
      }
      else if(type=="L2 ERROR_CALC") {
        {
          std::vector<std::string> values(2);
          values[0] = plist.get<std::string>("Field A");
          values[1] = plist.get<std::string>("Field B");
  
          std::vector<double> scalars(2); 
          scalars[0] = 1.0; 
          scalars[1] = -1.0;
  
          Teuchos::ParameterList p;
          p.set("Sum Name",key+"_DIFF"); // Name of sum
          p.set<RCP<std::vector<std::string> > >("Values Names",Teuchos::rcpFromRef(values));
          p.set<RCP<const std::vector<double> > >("Scalars",Teuchos::rcpFromRef(scalars));
          p.set("Data Layout",ir->dl_scalar);
  
          RCP< Evaluator<panzer::Traits> > e = 
                   rcp(new panzer::Sum<EvalT,panzer::Traits>(p));
  
          evaluators->push_back(e);
        }

        {
          Teuchos::RCP<const panzer::PointRule> pr = ir;
          std::vector<std::string> values(2);
          values[0] = key+"_DIFF";
          values[1] = key+"_DIFF";
  
          Teuchos::ParameterList p;
          p.set("Product Name",key);
          p.set("Values Names",Teuchos::rcpFromRef(values));
          p.set("Data Layout",ir->dl_scalar);
  
          RCP< Evaluator<panzer::Traits> > e = 
                   rcp(new panzer::Product<EvalT,panzer::Traits>(p));
  
          evaluators->push_back(e);
        }

        found = true;
      }
      else if(type=="H1 ERROR_CALC") {
        // Compute L2 contribution
        {
          std::vector<std::string> values(2);
          values[0] = plist.get<std::string>("Field A");
          values[1] = plist.get<std::string>("Field B");
  
          std::vector<double> scalars(2); 
          scalars[0] = 1.0; 
          scalars[1] = -1.0;
  
          Teuchos::ParameterList p;
          p.set("Sum Name",key+"_H1_L2DIFF"); // Name of sum
          p.set<RCP<std::vector<std::string> > >("Values Names",Teuchos::rcpFromRef(values));
          p.set<RCP<const std::vector<double> > >("Scalars",Teuchos::rcpFromRef(scalars));
          p.set("Data Layout",ir->dl_scalar);
  
          RCP< Evaluator<panzer::Traits> > e = 
                   rcp(new panzer::Sum<EvalT,panzer::Traits>(p));
  
          evaluators->push_back(e);
        }

        {
          std::vector<std::string> values(2);
          values[0] = "GRAD_"+plist.get<std::string>("Field A");
          values[1] = "GRAD_"+plist.get<std::string>("Field B");
  
          std::vector<double> scalars(2); 
          scalars[0] = 1.0; 
          scalars[1] = -1.0;
  
          Teuchos::ParameterList p;
          p.set("Sum Name","GRAD_"+key+"_DIFF"); // Name of sum
          p.set<RCP<std::vector<std::string> > >("Values Names",Teuchos::rcpFromRef(values));
          p.set<RCP<const std::vector<double> > >("Scalars",Teuchos::rcpFromRef(scalars));
          p.set("Data Layout",ir->dl_vector);
  
          RCP< Evaluator<panzer::Traits> > e = 
                   rcp(new panzer::Sum<EvalT,panzer::Traits>(p));
  
          evaluators->push_back(e);
        }

        {
          Teuchos::RCP<const panzer::PointRule> pr = ir;
          std::vector<std::string> values(2);
          values[0] = "GRAD_"+key+"_DIFF";
          values[1] = "GRAD_"+key+"_DIFF";
  
          Teuchos::ParameterList p;
          p.set("Product Name",key+"_H1Semi");
          p.set("Values Names",Teuchos::rcpFromRef(values));
          p.set("Data Layout",ir->dl_vector);
  
          RCP< Evaluator<panzer::Traits> > e = 
             panzer::buildEvaluator_DotProduct<EvalT,panzer::Traits>(key,*ir,values[0],values[1]);
  
          evaluators->push_back(e);
        }

        {
          Teuchos::RCP<const panzer::PointRule> pr = ir;
          std::vector<std::string> values(2);
          values[0] = key+"_H1_L2DIFF";
          values[1] = key+"_H1_L2DIFF";
  
          Teuchos::ParameterList p;
          p.set("Product Name",key+"_L2");
          p.set("Values Names",Teuchos::rcpFromRef(values));
          p.set("Data Layout",ir->dl_scalar);
  
          RCP< Evaluator<panzer::Traits> > e = 
                   rcp(new panzer::Product<EvalT,panzer::Traits>(p));
  
          evaluators->push_back(e);
        }

        {
          std::vector<std::string> values(2);
          values[0] = key+"_L2";
          values[1] = key+"_H1Semi";
  
          Teuchos::ParameterList p;
          p.set("Sum Name",key); // Name of sum
          p.set<RCP<std::vector<std::string> > >("Values Names",Teuchos::rcpFromRef(values));
          p.set("Data Layout",ir->dl_scalar);
  
          RCP< Evaluator<panzer::Traits> > e = 
                   rcp(new panzer::Sum<EvalT,panzer::Traits>(p));
  
          evaluators->push_back(e);
        }

        found = true;
      }
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
