// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_RESPONSE_SCATTER_EVALUATOR_IPCoordinates_HPP
#define PANZER_RESPONSE_SCATTER_EVALUATOR_IPCoordinates_HPP

#include <iostream>
#include <string>

#include "PanzerDiscFE_config.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_Response_IPCoordinates.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {

/** This class handles responses with values aggregated
  * on each finite element cell.
  */
template<typename EvalT, typename Traits>
class ResponseScatterEvaluator_IPCoordinates : public panzer::EvaluatorWithBaseImpl<Traits>,
                                               public PHX::EvaluatorDerived<EvalT, Traits>  { 
public:

  //! A constructor with concrete arguments instead of a parameter list.
  ResponseScatterEvaluator_IPCoordinates(const std::string & name,int ir_order);

  void postRegistrationSetup(typename Traits::SetupData d,
                             PHX::FieldManager<Traits>& fm);

  void evaluateFields(typename Traits::EvalData d);

  void preEvaluate(typename Traits::PreEvalData d);
  void postEvaluate(typename Traits::PostEvalData d);

private:
  typedef typename EvalT::ScalarT ScalarT;

  std::string responseName_;
  int ir_order_;
  int ir_index_;

  Teuchos::RCP<Response_IPCoordinates<EvalT> > responseObj_;
  std::vector<std::vector<ScalarT> > tmpCoords_;

  Teuchos::RCP<PHX::FieldTag> scatterHolder_; // dummy target
};

}

#include "Panzer_ResponseScatterEvaluator_IPCoordinates_impl.hpp"

#endif
