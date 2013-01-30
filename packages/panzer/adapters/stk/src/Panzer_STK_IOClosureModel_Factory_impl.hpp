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
   parseOutputList(outputList.sublist("Cell Quantities"),blockIdToCellFields_);
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
  return userCMF_->buildClosureModels(model_id,models,fl,ir,default_params,user_data,global_data,fm);
}

// ********************************************************************
// ********************************************************************

#endif
