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

#ifndef __Panzer_ResponseScatterEvaluator_impl_hpp__
#define __Panzer_ResponseScatterEvaluator_impl_hpp__

#include "Intrepid_FunctionSpaceTools.hpp"
#include "Panzer_String_Utilities.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <boost/io/ios_state.hpp>
#include <iomanip>

namespace panzer {

//**********************************************************************
template <typename EvalT,typename Traits,typename AggregatorT>
ResponseScatterEvaluator<EvalT,Traits,AggregatorT>
::ResponseScatterEvaluator(const std::string & name,
                           const Teuchos::RCP<panzer::ResponseData<Traits> > & data,
                           const Teuchos::RCP<const AggregatorT> & aggregator,
                           const std::vector<std::string> & responseNames,
                           int worksetSize)
{
  using Teuchos::RCP;

  // read in some information from response parameter list
  responseData_ = data;
  responseAggregator_ = aggregator;

  // build dummy tag to register with field manager
  responseDummyTag_ =
    Teuchos::rcp(new PHX::Tag<ScalarT>("Response Scatter: " + name,Teuchos::rcp(new PHX::MDALayout<Dummy>(0))));

  // add evaluated dummy filed and dependent response fields
  this->addEvaluatedField(*responseDummyTag_);

  Teuchos::RCP<PHX::DataLayout> dl_cell = Teuchos::rcp(new PHX::MDALayout<Cell>(worksetSize));
  for (std::vector<std::string>::const_iterator s=responseNames.begin(); s!=responseNames.end(); ++s) {
    PHX::MDField<ScalarT,Cell> field(*s,dl_cell);
    responseFields_.push_back(field); // record field in evaluator
    this->addDependentField(field);   // regiseter this as a dependent field
  }
   
  // add dependent fields to evaluator

  std::string n = "Response Scatter: " + name;
  this->setName(n);
}

//**********************************************************************
template <typename EvalT,typename Traits,typename AggregatorT>
void ResponseScatterEvaluator<EvalT,Traits,AggregatorT>::
postRegistrationSetup(typename Traits::SetupData sd,PHX::FieldManager<Traits>& fm)
{
   for(typename std::vector<PHX::MDField<ScalarT,Cell> >::iterator field = responseFields_.begin();
      field != responseFields_.end(); ++field)
     this->utils.setFieldData(*field,fm);
}

//**********************************************************************
template <typename EvalT,typename Traits,typename AggregatorT>
void ResponseScatterEvaluator<EvalT,Traits,AggregatorT>::
evaluateFields(typename Traits::EvalData workset)
{ 
  if (workset.num_cells == 0)
    return;

  // do some work
  responseAggregator_->template evaluateFields<PHX::MDField<ScalarT,Cell> > (workset,*responseData_,responseFields_);
}

} // namespace panzer

#endif
