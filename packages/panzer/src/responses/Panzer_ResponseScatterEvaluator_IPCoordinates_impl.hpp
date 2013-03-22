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

#ifndef PANZER_RESPONSE_SCATTER_EVALUATOR_IPCoordinates_IMPL_HPP
#define PANZER_RESPONSE_SCATTER_EVALUATOR_IPCoordinates_IMPL_HPP

#include <iostream>
#include <string>

#include "Panzer_config.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Panzer_ResponseBase.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_Workset_Utilities.hpp"

namespace panzer {

/** This class handles responses with values aggregated
  * on each finite element cell.
  */
template<typename EvalT, typename Traits>
ResponseScatterEvaluator_IPCoordinates<EvalT,Traits>::
ResponseScatterEvaluator_IPCoordinates(const std::string & name,
                                       int ir_order)
  : responseName_(name), ir_order_(ir_order)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  std::string dummyName = ResponseBase::buildLookupName(name) + " dummy target";

  // build dummy target tag
  RCP<PHX::DataLayout> dl_dummy = rcp(new PHX::MDALayout<panzer::Dummy>(0));
  scatterHolder_ = rcp(new PHX::Tag<ScalarT>(dummyName,dl_dummy));
  this->addEvaluatedField(*scatterHolder_);

  std::string n = "IPCoordinates Response Scatter: " + name;
  this->setName(n);
}

template<typename EvalT, typename Traits>
void ResponseScatterEvaluator_IPCoordinates<EvalT,Traits>::
preEvaluate(typename Traits::PreEvalData d)
{
  // extract response object
  responseObj_ = Teuchos::rcp_dynamic_cast<Response_IPCoordinates<EvalT> >(
                                   d.getDataObject(ResponseBase::buildLookupName(responseName_)),true);
}


template<typename EvalT, typename Traits>
void ResponseScatterEvaluator_IPCoordinates<EvalT,Traits>::
postRegistrationSetup(typename Traits::SetupData sd,
                      PHX::FieldManager<Traits>& fm)
{
  ir_index_ = panzer::getIntegrationRuleIndex(ir_order_,(*sd.worksets_)[0]);
}

template<typename EvalT, typename Traits>
void ResponseScatterEvaluator_IPCoordinates<EvalT,Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  Intrepid::FieldContainer<double>& workset_coords = (workset.int_rules[ir_index_])->ip_coordinates;

  if (tmpCoords_.size() != Teuchos::as<std::size_t>(workset_coords.dimension(2))) {
    tmpCoords_.resize(workset_coords.dimension(2));
    for(std::size_t dim=0;dim<tmpCoords_.size();dim++)
      tmpCoords_[dim].clear();
  }

  // This ordering is for the DataTransferKit.  It blocks all x
  // coordinates for a set of points, then all y coordinates and if
  // required all z coordinates.
  for (int dim = 0; dim < workset_coords.dimension(2); ++dim)
    for (std::size_t cell = 0; cell < workset.num_cells; ++cell)
      for (int ip = 0; ip < workset_coords.dimension(1); ++ip)
        tmpCoords_[dim].push_back(workset_coords(static_cast<int>(cell),ip,dim));
}

//**********************************************************************
template<typename EvalT, typename Traits>
void ResponseScatterEvaluator_IPCoordinates<EvalT,Traits>::
postEvaluate(typename Traits::PostEvalData data)
{
  std::vector<panzer::Traits::Residual::ScalarT> & coords = *responseObj_->getNonconstCoords();
  coords.clear();

  for (std::size_t dim = 0; dim < tmpCoords_.size(); ++dim) {
    for (typename std::vector<ScalarT>::const_iterator x=tmpCoords_[dim].begin(); x != tmpCoords_[dim].end(); ++ x)
      coords.push_back(Sacado::ScalarValue<ScalarT>::eval(*x));
  }

  tmpCoords_.clear();
}

}

#endif
