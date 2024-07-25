// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_STK_IOCLOSURE_MODEL_FACTORY_IMPL_HPP
#define PANZER_STK_IOCLOSURE_MODEL_FACTORY_IMPL_HPP

#include "Panzer_String_Utilities.hpp"

#include "Panzer_STK_ScatterCellAvgQuantity.hpp"
#include "Panzer_STK_ScatterCellQuantity.hpp"

// ********************************************************************
// ********************************************************************
template<typename EvalT>
panzer_stk::IOClosureModelFactory<EvalT>::
IOClosureModelFactory(const Teuchos::RCP<const panzer::ClosureModelFactory<EvalT> > userCMF,
                      const Teuchos::RCP<STK_Interface> & mesh,
                      const Teuchos::ParameterList & outputList)
   : mesh_(mesh), userCMF_(userCMF)
{
   parseOutputList(outputList.sublist("Cell Average Quantities"),blockIdToCellAvgFields_);
   parseOutputList(outputList.sublist("Cell Average Vectors"),blockIdToCellAvgVectors_);
   parseOutputList(outputList.sublist("Cell Quantities"),blockIdToCellFields_);
   parseOutputList(outputList.sublist("Nodal Quantities"),blockIdToNodalFields_);
}

// ********************************************************************
// ********************************************************************
template<typename EvalT>
panzer_stk::IOClosureModelFactory<EvalT>::
IOClosureModelFactory(const Teuchos::RCP<const panzer::ClosureModelFactory<EvalT> > userCMF,
                      const Teuchos::RCP<STK_Interface> & mesh,
                      const std::map<std::string,std::vector<std::string> > & nodalFields,
                      const std::map<std::string,std::vector<std::string> > & cellFields)
   : mesh_(mesh), userCMF_(userCMF)
{
   blockIdToNodalFields_ = nodalFields;
   blockIdToCellFields_ = cellFields;

   typedef std::map<std::string,std::vector<std::string> >::const_iterator const_iterator;

   for(const_iterator itr=nodalFields.begin();itr!=nodalFields.end();++itr) 
      blockIdEvaluated_[itr->first] = false;
   for(const_iterator itr=cellFields.begin();itr!=cellFields.end();++itr) 
      blockIdEvaluated_[itr->first] = false;
}

// ********************************************************************
// ********************************************************************
template<typename EvalT>
void panzer_stk::IOClosureModelFactory<EvalT>::
parseOutputList(const Teuchos::ParameterList & pl,
                std::map<std::string,std::vector<std::string> > & blockIdToFields) const
{
   for(Teuchos::ParameterList::ConstIterator itr=pl.begin();
       itr!=pl.end();++itr) {
      const std::string & blockId = itr->first;
      const std::string & fields = Teuchos::any_cast<std::string>(itr->second.getAny());
      std::vector<std::string> & tokens = blockIdToFields[blockId];
 
      // break up comma seperated fields
      panzer::StringTokenizer(tokens,fields,",",true);

      blockIdEvaluated_[blockId] = false; // initialize this value
   }
}

// ********************************************************************
// ********************************************************************
template<typename EvalT>
Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > 
panzer_stk::IOClosureModelFactory<EvalT>::
buildClosureModels(const std::string& model_id,
		   const Teuchos::ParameterList& models, 
		   const panzer::FieldLayoutLibrary& fl,
		   const Teuchos::RCP<panzer::IntegrationRule>& ir,
		   const Teuchos::ParameterList& default_params, 
		   const Teuchos::ParameterList& user_data,
		   const Teuchos::RCP<panzer::GlobalData>& global_data,
		   PHX::FieldManager<panzer::Traits>& fm) const
{
  // Note that the Residual version of this is in the cpp file!!!!

  return userCMF_->buildClosureModels(model_id,models,fl,ir,default_params,user_data,global_data,fm);
}

// ********************************************************************
// ********************************************************************

#endif
