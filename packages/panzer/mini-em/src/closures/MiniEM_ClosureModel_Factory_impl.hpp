#ifndef __MiniEM_ClosureModelFactoryT_hpp__
#define __MiniEM_ClosureModelFactoryT_hpp__

#include <iostream>
#include <sstream>
#include <typeinfo>
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Integrator_Scalar.hpp"
#include "Panzer_DotProduct.hpp"
#include "Panzer_Sum.hpp"
#include "Phalanx_FieldTag_Tag.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_TypeNameTraits.hpp"

#include "MiniEM_GaussianPulse.hpp"

// ********************************************************************
// ********************************************************************
template<typename EvalT>
Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > 
mini_em::ModelFactory<EvalT>::
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
      if(type=="GAUSSIAN PULSE") {
        double dt = plist.get<double>("dt");
	RCP< Evaluator<panzer::Traits> > e = 
	  rcp(new mini_em::GaussianPulse<EvalT,panzer::Traits>(key,*ir,fl,dt));
	evaluators->push_back(e);

        found = true;
      }
      if(type=="ELECTROMAGNETIC ENERGY") {
        double mu = plist.get<double>("mu");
        double epsilon = plist.get<double>("epsilon");

        // compute ||E||^2
        {
          Teuchos::ParameterList input;
          input.set("Result Name", "E_SQUARED");
          input.set<Teuchos::RCP<const panzer::PointRule> >("Point Rule", ir);
          input.set("Vector A Name", "E_edge");
          input.set("Vector B Name", "E_edge");

          RCP< Evaluator<panzer::Traits> > e = 
	    rcp(new panzer::DotProduct<EvalT,panzer::Traits>(input));
	  evaluators->push_back(e);
        }

        // compute ||B||^2
        {
          Teuchos::ParameterList input;
          input.set("Result Name", "B_SQUARED");
          input.set<Teuchos::RCP<const panzer::PointRule> >("Point Rule", ir);
          input.set("Vector A Name", "B_face");
          input.set("Vector B Name", "B_face");

          RCP< Evaluator<panzer::Traits> > e = 
	    rcp(new panzer::DotProduct<EvalT,panzer::Traits>(input));
	  evaluators->push_back(e);
        }

        // compute epsilon/2*||E||^2 + 1/(2*mu)*||B||^2
        {
          RCP<std::vector<double> > coeffs = rcp(new std::vector<double>);
          coeffs->push_back(0.5*epsilon);
          coeffs->push_back(0.5/mu);
  
          RCP<std::vector<std::string> > valuesNames = rcp(new std::vector<std::string>);
          valuesNames->push_back("E_SQUARED");
          valuesNames->push_back("B_SQUARED");

          Teuchos::ParameterList input;
          input.set("Sum Name","EM_ENERGY");
          input.set("Values Names",valuesNames);
          input.set("Data Layout",ir->dl_scalar);
          input.set<RCP<const std::vector<double> > >("Scalars", coeffs);
        
          RCP< Evaluator<panzer::Traits> > e = 
	    rcp(new panzer::Sum<EvalT,panzer::Traits>(input));
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
