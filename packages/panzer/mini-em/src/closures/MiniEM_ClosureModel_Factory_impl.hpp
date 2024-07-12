// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __MiniEM_ClosureModelFactoryT_hpp__
#define __MiniEM_ClosureModelFactoryT_hpp__

#include <iostream>
#include <sstream>
#include <typeinfo>
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Integrator_Scalar.hpp"
#include "Panzer_Product.hpp"
#include "Panzer_DotProduct.hpp"
#include "Panzer_Sum.hpp"
#include "Phalanx_FieldTag_Tag.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_TypeNameTraits.hpp"

#include "MiniEM_GaussianPulse.hpp"
#include "MiniEM_RandomForcing.hpp"
#include "MiniEM_MaxwellAnalyticForcing.hpp"
#include "MiniEM_MaxwellAnalyticSolution.hpp"
#include "MiniEM_DarcyAnalyticForcing.hpp"
#include "MiniEM_DarcyAnalyticSolution.hpp"
#include "MiniEM_PiecewiseConstant.hpp"
#include "MiniEM_TensorConductivity.hpp"
#include "MiniEM_VariableTensorConductivity.hpp"

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
    const Teuchos::ParameterEntry& entry = model_it->second;
    const ParameterList& plist = Teuchos::getValue<Teuchos::ParameterList>(entry);

    if (plist.isType<double>("Value")) {
      { // at IP
        ParameterList input;
	input.set("Name", key);
	input.set("Value", plist.get<double>("Value"));
	input.set("Data Layout", ir->dl_scalar);
	RCP< Evaluator<panzer::Traits> > e = 
	  rcp(new panzer::Constant<EvalT,panzer::Traits>(input));
	evaluators->push_back(e);
      }
      
      for (std::vector<Teuchos::RCP<const panzer::PureBasis> >::const_iterator basis_itr = bases.begin();
	   basis_itr != bases.end(); ++basis_itr) { // at BASIS
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

    if (plist.isType<std::string>("Type")) {
      std::string type = plist.get<std::string>("Type");
      if(type=="GAUSSIAN PULSE") {
        double dt = plist.get<double>("dt");
	RCP< Evaluator<panzer::Traits> > e = 
	  rcp(new mini_em::GaussianPulse<EvalT,panzer::Traits>(key,*ir,fl,dt));
	evaluators->push_back(e);

        found = true;
      }
      if(type=="RANDOM") {
        unsigned int seed = plist.get<unsigned int>("seed");
        double min = plist.get<double>("range min");
        double max = plist.get<double>("range max");
        std::string basisName = plist.get<std::string>("DoF Name");
	RCP< Evaluator<panzer::Traits> > e =
	  rcp(new mini_em::RandomForcing<EvalT,panzer::Traits>(key,*ir,fl,seed,min,max,basisName));
	evaluators->push_back(e);

        found = true;
      }
      if(type=="DARCY ANALYTIC FORCING") {
        double kappa = plist.get<double>("kappa");
	RCP<Evaluator<panzer::Traits> > e =
	  rcp(new mini_em::DarcyAnalyticForcing<EvalT,panzer::Traits>(key,*ir,fl, kappa));
	evaluators->push_back(e);

        found = true;
      }
      if(type=="DARCY ANALYTIC SOLUTION") {
	RCP<Evaluator<panzer::Traits> > e =
	  rcp(new mini_em::DarcyAnalyticSolution<EvalT,panzer::Traits>(key,*ir,fl));
	evaluators->push_back(e);

        found = true;
      }
      if(type=="MAXWELL ANALYTIC FORCING") {
        double epsilon = plist.get<double>("epsilon");
        double timeScale = plist.get<double>("Time scale");
	RCP<Evaluator<panzer::Traits> > e =
	  rcp(new mini_em::MaxwellAnalyticForcing<EvalT,panzer::Traits>(key,*ir,fl,epsilon,timeScale));
	evaluators->push_back(e);

        found = true;
      }
      if(type=="MAXWELL ANALYTIC SOLUTION") {
        double timeScale = plist.get<double>("Time scale");
	RCP<Evaluator<panzer::Traits> > e =
	  rcp(new mini_em::MaxwellAnalyticSolution<EvalT,panzer::Traits>(key,*ir,fl,timeScale));
	evaluators->push_back(e);

        found = true;
      }
      if(type=="PIECEWISE CONSTANT") {
        double value0 = plist.get<double>("value0");
	double value1 = plist.get<double>("value1");
        double xl = plist.get<double>("xl");
	double xr = plist.get<double>("xr");
	double yl = plist.get<double>("yl");
	double yr = plist.get<double>("yr");
	double zl = plist.get<double>("zl");
	double zr = plist.get<double>("zr");
        std::string DoF = plist.get<std::string>("DoF Name");
	RCP< Evaluator<panzer::Traits> > e =
	  rcp(new mini_em::PiecewiseConstant<EvalT,panzer::Traits>(key,*ir,fl,value0,value1,xl,xr,yl,yr,zl,zr,DoF));
	evaluators->push_back(e);

        found = true;
      }
      if(type=="TENSOR CONDUCTIVITY") {
        double sigma = plist.get<double>("sigma");
        double betax = plist.get<double>("betax");
        double betay = plist.get<double>("betay");
        double betaz = plist.get<double>("betaz");
        std::string DoF = plist.get<std::string>("DoF Name");
	RCP< Evaluator<panzer::Traits> > e =
	  rcp(new mini_em::TensorConductivity<EvalT,panzer::Traits>(key,*ir,fl,sigma,betax,betay,betaz,DoF));
	evaluators->push_back(e);

        found = true;
      }
      if(type=="VARIABLE TENSOR CONDUCTIVITY") {
        double sigma0 = plist.get<double>("sigma0");
        double betax0 = plist.get<double>("betax0");
        double betay0 = plist.get<double>("betay0");
        double betaz0 = plist.get<double>("betaz0");
        double sigma1 = plist.get<double>("sigma1");
        double betax1 = plist.get<double>("betax1");
        double betay1 = plist.get<double>("betay1");
        double betaz1 = plist.get<double>("betaz1");
        double sigma2 = plist.get<double>("sigma2");
        double betax2 = plist.get<double>("betax2");
        double betay2 = plist.get<double>("betay2");
        double betaz2 = plist.get<double>("betaz2");
        std::string DoF = plist.get<std::string>("DoF Name");
	RCP< Evaluator<panzer::Traits> > e =
	  rcp(new mini_em::VariableTensorConductivity<EvalT,panzer::Traits>(key,*ir,fl,sigma0,sigma1,sigma2,
                                                                            betax0,betay0,betaz0,
                                                                            betax1,betay1,betaz1,
                                                                            betax2,betay2,betaz2,
                                                                            DoF));
	evaluators->push_back(e);

        found = true;
      }
      if(type=="ELECTROMAGNETIC ENERGY") {
        // compute (E, epsilon*E)
        {
          Teuchos::ParameterList input;
          input.set("Result Name", "E_SQUARED");
          input.set<Teuchos::RCP<const panzer::PointRule> >("Point Rule", ir);
          input.set("Vector A Name", "E_edge");
          input.set("Vector B Name", "E_edge");
          input.set("Field Multiplier", "epsilon");

          RCP< Evaluator<panzer::Traits> > e = 
	    rcp(new panzer::DotProduct<EvalT,panzer::Traits>(input));
	  evaluators->push_back(e);
        }

        // compute (B, 1/mu * B)
        {
          if (ir->spatial_dimension == 3) {
            Teuchos::ParameterList input;
            input.set("Result Name", "B_SQUARED");
            input.set<Teuchos::RCP<const panzer::PointRule> >("Point Rule", ir);
            input.set("Vector A Name", "B_face");
            input.set("Vector B Name", "B_face");
            input.set("Field Multiplier", "1/mu");

            RCP< Evaluator<panzer::Traits> > e =
              rcp(new panzer::DotProduct<EvalT,panzer::Traits>(input));
            evaluators->push_back(e);
          } else if (ir->spatial_dimension == 2) {
            Teuchos::ParameterList input;
            input.set("Product Name", "B_SQUARED");
            RCP<std::vector<std::string> > valuesNames = rcp(new std::vector<std::string>);
            valuesNames->push_back("B_face");
            valuesNames->push_back("B_face");
            input.set("Values Names",valuesNames);
            input.set("Data Layout",ir->dl_scalar);
            input.set("Field Multiplier", "1/mu");

            RCP< Evaluator<panzer::Traits> > e =
              rcp(new panzer::Product<EvalT,panzer::Traits>(input));
            evaluators->push_back(e);
          }
        }

        // compute 1/2*(E, epsilon * E) + 1/2*(B, 1/mu * B)
        {
          RCP<std::vector<double> > coeffs = rcp(new std::vector<double>);
          coeffs->push_back(0.5);
          coeffs->push_back(0.5);
  
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

        {
          RCP<std::vector<std::string> > valuesNames = rcp(new std::vector<std::string>);
          valuesNames->push_back("EM_ENERGY");
          valuesNames->push_back("1/dt");
          valuesNames->push_back("1/dt");

          Teuchos::ParameterList input;
          input.set("Product Name","EM_ENERGY/dt^2");
          input.set("Values Names",valuesNames);
          input.set("Data Layout",ir->dl_scalar);

          RCP< Evaluator<panzer::Traits> > e =
	    rcp(new panzer::Product<EvalT,panzer::Traits>(input));
	  evaluators->push_back(e);
        }
 
        found = true;
      }
      if(type=="NORM") {
        // compute ||u||^2
        {
          Teuchos::ParameterList input;
          input.set("Product Name",key);
          RCP<std::vector<std::string> > valuesNames = rcp(new std::vector<std::string>);
          valuesNames->push_back(plist.get<std::string>("Field"));
          valuesNames->push_back(plist.get<std::string>("Field"));
          input.set("Values Names",valuesNames);
          input.set("Data Layout",ir->dl_scalar);

          RCP< Evaluator<panzer::Traits> > e =
            rcp(new panzer::Product<EvalT,panzer::Traits>(input));
          evaluators->push_back(e);
        }

        found = true;
      }
      if(type=="ERROR") {
        // compute ||E-E_ex||^2
        bool vectorial = false;
        if (plist.isType<bool>("Vectorial"))
          vectorial = plist.get<bool>("Vectorial");
        const std::string diffName = "DIFFERENCE_" + plist.get<std::string>("Field") + "_" + plist.get<std::string>("Exact Field");
        if (vectorial) {
          {
            RCP<std::vector<double> > coeffs = rcp(new std::vector<double>);
            coeffs->push_back(1);
            coeffs->push_back(-1);

            RCP<std::vector<std::string> > valuesNames = rcp(new std::vector<std::string>);
            valuesNames->push_back(plist.get<std::string>("Field"));
            valuesNames->push_back(plist.get<std::string>("Exact Field"));

            Teuchos::ParameterList input;
            input.set("Sum Name",diffName);
            input.set("Values Names",valuesNames);
            input.set("Data Layout",ir->dl_vector);
            input.set<RCP<const std::vector<double> > >("Scalars", coeffs);

            RCP< Evaluator<panzer::Traits> > e =
              rcp(new panzer::Sum<EvalT,panzer::Traits>(input));
            evaluators->push_back(e);
          }
          {
            Teuchos::ParameterList input;
            input.set("Result Name", key);
            input.set<Teuchos::RCP<const panzer::PointRule> >("Point Rule", ir);
            input.set("Vector A Name", diffName);
            input.set("Vector B Name", diffName);

            RCP< Evaluator<panzer::Traits> > e =
              rcp(new panzer::DotProduct<EvalT,panzer::Traits>(input));
            evaluators->push_back(e);
          }
        } else {
          {
            RCP<std::vector<double> > coeffs = rcp(new std::vector<double>);
            coeffs->push_back(1);
            coeffs->push_back(-1);

            RCP<std::vector<std::string> > valuesNames = rcp(new std::vector<std::string>);
            valuesNames->push_back(plist.get<std::string>("Field"));
            valuesNames->push_back(plist.get<std::string>("Exact Field"));

            Teuchos::ParameterList input;
            input.set("Sum Name",diffName);
            input.set("Values Names",valuesNames);
            input.set("Data Layout",ir->dl_scalar);
            input.set<RCP<const std::vector<double> > >("Scalars", coeffs);

            RCP< Evaluator<panzer::Traits> > e =
              rcp(new panzer::Sum<EvalT,panzer::Traits>(input));
            evaluators->push_back(e);
          }
          {
            Teuchos::ParameterList input;
            input.set("Product Name",key);
            RCP<std::vector<std::string> > valuesNames = rcp(new std::vector<std::string>);
            valuesNames->push_back(diffName);
            valuesNames->push_back(diffName);
            input.set("Values Names",valuesNames);
            input.set("Data Layout",ir->dl_scalar);

            RCP< Evaluator<panzer::Traits> > e =
              rcp(new panzer::Product<EvalT,panzer::Traits>(input));
            evaluators->push_back(e);
          }
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
