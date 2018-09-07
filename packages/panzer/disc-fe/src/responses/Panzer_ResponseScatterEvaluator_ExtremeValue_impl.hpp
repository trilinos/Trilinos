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

#ifndef PANZER_RESPONSE_SCATTER_EVALUATOR_EXTREMEVALUE_IMPL_HPP
#define PANZER_RESPONSE_SCATTER_EVALUATOR_EXTREMEVALUE_IMPL_HPP

#include <iostream>
#include <string>

#include "PanzerDiscFE_config.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Panzer_ResponseBase.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"

#include "Thyra_SpmdVectorBase.hpp"
#include "Teuchos_ArrayRCP.hpp"

namespace panzer {

/** This class handles responses with values aggregated
  * on each finite element cell.
  */
template<typename EvalT, typename Traits>
ResponseScatterEvaluator_ExtremeValue<EvalT,Traits>::
ResponseScatterEvaluator_ExtremeValue(const std::string & name,
                                    const CellData & cd,
                                    bool useMax,
                                    const Teuchos::RCP<ExtremeValueScatterBase> & extremeValueScatter)
  : responseName_(name)
  , scatterObj_(extremeValueScatter)
  , useMax_(useMax)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  std::string dummyName = ResponseBase::buildLookupName(name) + " dummy target";

  // build dummy target tag
  RCP<PHX::DataLayout> dl_dummy = rcp(new PHX::MDALayout<panzer::Dummy>(0));
  scatterHolder_ = rcp(new PHX::Tag<ScalarT>(dummyName,dl_dummy));
  this->addEvaluatedField(*scatterHolder_);

  // build dendent field
  RCP<PHX::DataLayout> dl_cell = rcp(new PHX::MDALayout<panzer::Cell>(cd.numCells()));
  cellExtremeValue_ = PHX::MDField<const ScalarT,panzer::Cell>(name,dl_cell);
  this->addDependentField(cellExtremeValue_);

  std::string n = "Extreme Value Response Scatter: " + name;
  this->setName(n);
}

template<typename EvalT, typename Traits>
ResponseScatterEvaluator_ExtremeValue<EvalT,Traits>::
ResponseScatterEvaluator_ExtremeValue(const std::string & integrandName,
                                    const std::string & responseName,
                                    const CellData & cd,
                                    bool useMax,
                                    const Teuchos::RCP<ExtremeValueScatterBase> & extremeValueScatter)
  : responseName_(responseName)
  , scatterObj_(extremeValueScatter)
  , useMax_(useMax)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  std::string dummyName = ResponseBase::buildLookupName(responseName) + " dummy target";

  // build dummy target tag
  RCP<PHX::DataLayout> dl_dummy = rcp(new PHX::MDALayout<panzer::Dummy>(0));
  scatterHolder_ = rcp(new PHX::Tag<ScalarT>(dummyName,dl_dummy));
  this->addEvaluatedField(*scatterHolder_);

  // build dendent field
  RCP<PHX::DataLayout> dl_cell = rcp(new PHX::MDALayout<panzer::Cell>(cd.numCells()));
  cellExtremeValue_ = PHX::MDField<const ScalarT,panzer::Cell>(integrandName,dl_cell);
  this->addDependentField(cellExtremeValue_);

  std::string n = "Extreme Value Response Scatter: " + responseName;
  this->setName(n);
}

template<typename EvalT, typename Traits>
void ResponseScatterEvaluator_ExtremeValue<EvalT,Traits>::
preEvaluate(typename Traits::PreEvalData d)
{
  // extract linear object container
  responseObj_ = Teuchos::rcp_dynamic_cast<Response_ExtremeValue<EvalT> >(
                                   d.gedc->getDataObject(ResponseBase::buildLookupName(responseName_)),true);
}


template<typename EvalT, typename Traits>
void ResponseScatterEvaluator_ExtremeValue<EvalT,Traits>::
evaluateFields(typename Traits::EvalData d)
{
  for(index_t i=0;i<d.num_cells;i++) {
    if(useMax_)
      responseObj_->value = (responseObj_->value < cellExtremeValue_(i)) ? cellExtremeValue_(i) : responseObj_->value;
    else
      responseObj_->value = (responseObj_->value > cellExtremeValue_(i)) ? cellExtremeValue_(i) : responseObj_->value;
  }
}

template < >
void ResponseScatterEvaluator_ExtremeValue<panzer::Traits::Jacobian,panzer::Traits>::
evaluateFields(panzer::Traits::EvalData d)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::SpmdVectorBase;

  TEUCHOS_ASSERT(scatterObj_!=Teuchos::null);

  // grab local data for inputing
  Teuchos::ArrayRCP<double> local_dgdx;
  RCP<SpmdVectorBase<double> > dgdx = rcp_dynamic_cast<SpmdVectorBase<double> >(responseObj_->getGhostedVector());
  dgdx->getNonconstLocalData(ptrFromRef(local_dgdx));
  TEUCHOS_ASSERT(!local_dgdx.is_null());

  scatterObj_->scatterDerivative(cellExtremeValue_,d,this->wda,local_dgdx);
}

}

#endif
