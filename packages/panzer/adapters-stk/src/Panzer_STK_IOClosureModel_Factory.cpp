// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_STK_IOClosureModel_Factory.hpp"
#include "Panzer_STK_ScatterCellAvgQuantity.hpp"
#include "Panzer_STK_ScatterCellAvgVector.hpp"
#include "Panzer_STK_ScatterCellQuantity.hpp"
#include "Panzer_STK_ScatterFields.hpp"

template< >
Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > 
panzer_stk::IOClosureModelFactory<panzer::Traits::Residual>::
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
  using PHX::Evaluator;

  // build user evaluators
  RCP< std::vector< RCP<Evaluator<panzer::Traits> > > > user_evals = 
    userCMF_->buildClosureModels(model_id,models,fl,ir,default_params,user_data,global_data,fm);

  // add user evaluators to evaluator list
  RCP< std::vector< RCP<Evaluator<panzer::Traits> > > > evaluators = 
    rcp(new std::vector< RCP<Evaluator<panzer::Traits> > > );

  // extract element block id
  std::string block_id = default_params.get<std::string>("Block ID");

  if(!blockIdEvaluated_[block_id]) {
     typedef std::map<std::string,std::vector<std::string> > BlockIdToFields;

     int worksetsize = ir->dl_scalar->extent(0);

     // see if there's a scaling parameter object
     Teuchos::RCP<std::map<std::string,double>> varScaleFactors;
     if (user_data.isParameter("Variable Scale Factors Map"))
     {
       varScaleFactors = user_data.get<Teuchos::RCP<std::map<std::string,double>>>("Variable Scale Factors Map");
     }

     // if a requested field is found then add in cell avg quantity evaluator
     BlockIdToFields::const_iterator cellAvgItr = blockIdToCellAvgFields_.find(block_id);
     if(cellAvgItr!=blockIdToCellAvgFields_.end() ) {
        Teuchos::RCP<std::vector<std::string> > fieldNames = Teuchos::rcp(new std::vector<std::string>(cellAvgItr->second));
   
        // setup averge cell fields
        Teuchos::ParameterList pl;
        pl.set("Mesh",mesh_);
        pl.set("IR",ir);
        pl.set("Field Names",fieldNames);
        pl.set("Scatter Name", block_id+"_Cell_Avg_Fields");
        pl.set("Variable Scale Factors Map", varScaleFactors);
        Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval
            = Teuchos::rcp(new panzer_stk::ScatterCellAvgQuantity<panzer::Traits::Residual,panzer::Traits>(pl));
        fm.registerEvaluator<panzer::Traits::Residual>(eval);
        fm.requireField<panzer::Traits::Residual>(*eval->evaluatedFields()[0]);
   
        evaluators->push_back(eval);
   
        blockIdEvaluated_[block_id] = true;
     } 

     // if a requested field is found then add in cell avg vector evaluator
     BlockIdToFields::const_iterator cellAvgVecItr = blockIdToCellAvgVectors_.find(block_id);
     if(cellAvgVecItr != blockIdToCellAvgVectors_.end() ) {
        Teuchos::RCP<std::vector<std::string> > fieldNames = Teuchos::rcp(new std::vector<std::string>(cellAvgVecItr->second));
   
        // setup cell average vectors
        Teuchos::ParameterList pl;
        pl.set("Mesh",mesh_);
        pl.set("IR",ir);
        pl.set("Field Names",fieldNames);
        pl.set("Scatter Name", block_id+"_Cell_Avg_Vectors");
        pl.set("Variable Scale Factors Map", varScaleFactors);
        Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval
            = Teuchos::rcp(new panzer_stk::ScatterCellAvgVector<panzer::Traits::Residual,panzer::Traits>(pl));
        fm.registerEvaluator<panzer::Traits::Residual>(eval);
        fm.requireField<panzer::Traits::Residual>(*eval->evaluatedFields()[0]);
   
        evaluators->push_back(eval);
   
        blockIdEvaluated_[block_id] = true;
     } 

     // if a requested field is found then add in cell quantity evaluator
     BlockIdToFields::const_iterator cellItr = blockIdToCellFields_.find(block_id);
     if(cellItr!=blockIdToCellFields_.end() ) {
        Teuchos::RCP<std::vector<std::string> > fieldNames = Teuchos::rcp(new std::vector<std::string>(cellItr->second));
   
        // setup averge cell fields
        Teuchos::ParameterList pl;
        pl.set("Mesh",mesh_);
        pl.set("Workset Size",worksetsize);
        pl.set("Field Names",fieldNames);
        pl.set("Scatter Name", block_id+"_Cell_Fields");
        pl.set("Variable Scale Factors Map", varScaleFactors);
        Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval
            = Teuchos::rcp(new panzer_stk::ScatterCellQuantity<panzer::Traits::Residual,panzer::Traits>(pl));
        fm.registerEvaluator<panzer::Traits::Residual>(eval);
        fm.requireField<panzer::Traits::Residual>(*eval->evaluatedFields()[0]);
   
        evaluators->push_back(eval);
   
        blockIdEvaluated_[block_id] = true;
     } 

     // if a requested field is found then add in cell quantity evaluator
     BlockIdToFields::const_iterator nodalItr = blockIdToNodalFields_.find(block_id);
     if(nodalItr!=blockIdToNodalFields_.end() ) {
        Teuchos::RCP<std::vector<std::string> > fieldNames = Teuchos::rcp(new std::vector<std::string>(nodalItr->second));

        Teuchos::RCP<const panzer::PureBasis> basis = Teuchos::rcp(new panzer::PureBasis("HGrad",1,ir->workset_size,ir->topology));

        // setup scaling factors as accepted by ScatterFields, if present
        std::vector<double> scale_factors_(fieldNames->size(),1.0);
        if (!varScaleFactors.is_null()) {
          for(size_t f=0; f < fieldNames->size(); ++f) {
            std::map<std::string,double>::const_iterator f2s_itr = varScaleFactors->find((*fieldNames)[f]);
            if(f2s_itr != varScaleFactors->end()) {
              scale_factors_[f] = f2s_itr->second;
            }
          }
        }

        // setup scatter nodal fields
        Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval
          = Teuchos::rcp(new panzer_stk::ScatterFields<panzer::Traits::Residual,panzer::Traits>(block_id+"Nodal_Fields",mesh_,basis,*fieldNames,scale_factors_));
        fm.registerEvaluator<panzer::Traits::Residual>(eval);
        fm.requireField<panzer::Traits::Residual>(*eval->evaluatedFields()[0]);
   
        evaluators->push_back(eval);
   
        blockIdEvaluated_[block_id] = true;
     } 
  }

  evaluators->insert(evaluators->end(),user_evals->begin(),user_evals->end()); 

  return evaluators;
}

#ifdef HAVE_PANZER_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Panzer_STK_IOClosureModel_Factory_decl.hpp"
#include "Panzer_STK_IOClosureModel_Factory_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_ONE_T(panzer_stk::IOClosureModelFactory)

#endif
