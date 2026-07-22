// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_RESPONSE_SCATTER_EVALUATOR_IPCoordinates_IMPL_HPP
#define PANZER_RESPONSE_SCATTER_EVALUATOR_IPCoordinates_IMPL_HPP

#include <iostream>
#include <string>

#include "PanzerDiscFE_config.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Panzer_ResponseBase.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"

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
                                   d.gedc->getDataObject(ResponseBase::buildLookupName(responseName_)),true);
}


template<typename EvalT, typename Traits>
void ResponseScatterEvaluator_IPCoordinates<EvalT,Traits>::
postRegistrationSetup(typename Traits::SetupData sd,
                      PHX::FieldManager<Traits>& /* fm */)
{
  ir_index_ = panzer::getIntegrationRuleIndex(ir_order_,(*sd.worksets_)[0], this->wda);
}

template<typename EvalT, typename Traits>
void ResponseScatterEvaluator_IPCoordinates<EvalT,Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  // Kokkos::DynRankView<double,PHX::Device>& workset_coords = (this->wda(workset).int_rules[ir_index_])->ip_coordinates;
  IntegrationValues2<double> & iv = *this->wda(workset).int_rules[ir_index_];

  if (tmpCoords_.size() != Teuchos::as<std::size_t>(iv.ip_coordinates.extent(2))) {
    tmpCoords_.resize(iv.ip_coordinates.extent(2));
    for(std::size_t dim=0;dim<tmpCoords_.size();dim++)
      tmpCoords_[dim].clear();
  }

  auto ip_coordinates_h = Kokkos::create_mirror_view(PHX::as_view(iv.ip_coordinates));
  Kokkos::deep_copy(ip_coordinates_h, PHX::as_view(iv.ip_coordinates));

  // This ordering is for the DataTransferKit.  It blocks all x
  // coordinates for a set of points, then all y coordinates and if
  // required all z coordinates.
  for (int dim = 0; dim < iv.ip_coordinates.extent_int(2); ++dim)
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
      for (int ip = 0; ip < iv.ip_coordinates.extent_int(1); ++ip)
        tmpCoords_[dim].push_back(ip_coordinates_h(static_cast<int>(cell),ip,dim));
}

//**********************************************************************
template<typename EvalT, typename Traits>
void ResponseScatterEvaluator_IPCoordinates<EvalT,Traits>::
postEvaluate(typename Traits::PostEvalData /* data */)
{
  std::vector<panzer::Traits::Residual::ScalarT> & coords = *responseObj_->getNonconstCoords();
  coords.clear();

  for (std::size_t dim = 0; dim < tmpCoords_.size(); ++dim) {
    for (typename std::vector<ScalarT>::const_iterator x=tmpCoords_[dim].begin(); x != tmpCoords_[dim].end(); ++ x)
      coords.push_back(Sacado::scalarValue(*x));
  }

  tmpCoords_.clear();
}

}

#endif
