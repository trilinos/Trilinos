#ifndef PANZER_STK_IOCLOSURE_MODEL_FACTORY_T_HPP
#define PANZER_STK_IOCLOSURE_MODEL_FACTORY_T_HPP

#include "Panzer_String_Utilities.hpp"

#include "Panzer_STK_ScatterCellAvgQuantity.hpp"

// ********************************************************************
// ********************************************************************
template<typename EvalT>
panzer_stk::IOClosureModelFactory<EvalT>::
IOClosureModelFactory(const Teuchos::RCP<const panzer::ClosureModelFactory<EvalT> > userCMF,
                      const Teuchos::RCP<STK_Interface> & mesh,
                      const Teuchos::ParameterList & outputList)
   : mesh_(mesh), userCMF_(userCMF)
{
   parseOutputList(outputList.sublist("Cell Average Quantities"));
}

// ********************************************************************
// ********************************************************************
template<typename EvalT>
void panzer_stk::IOClosureModelFactory<EvalT>::
parseOutputList(const Teuchos::ParameterList & pl)
{
   for(Teuchos::ParameterList::ConstIterator itr=pl.begin();
       itr!=pl.end();++itr) {
      const std::string & blockId = itr->first;
      const std::string & fields = Teuchos::any_cast<std::string>(itr->second.getAny());
      std::vector<std::string> & tokens = blockIdToFields_[blockId];
 
      // break up comma seperated fields
      panzer::StringTokenizer(tokens,fields);
   }
}

// ********************************************************************
// ********************************************************************
template<typename EvalT>
Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > 
panzer_stk::IOClosureModelFactory<EvalT>::
buildClosureModels(const std::string& model_id,
		   const panzer::InputEquationSet& set,
		   const Teuchos::ParameterList& models, 
		   const Teuchos::ParameterList& default_params, 
		   const Teuchos::ParameterList& user_data,
		   PHX::FieldManager<panzer::Traits>& fm) const
{
  return userCMF_->buildClosureModels(model_id,set,models,default_params,user_data,fm);
}

// ********************************************************************
// ********************************************************************

#endif
