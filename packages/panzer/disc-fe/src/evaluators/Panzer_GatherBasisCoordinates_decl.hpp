// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_GATHER_BASIS_COORDINATES_DECL_HPP
#define PANZER_EVALUATOR_GATHER_BASIS_COORDINATES_DECL_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CloneableEvaluator.hpp"

#include "Panzer_PureBasis.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {

/** \brief Gathers coordinates for the basis function from the
    workset and stores them in the field manager.
*/
template<typename EvalT, typename TRAITS>
class GatherBasisCoordinates
  : public panzer::EvaluatorWithBaseImpl<TRAITS>,
    public PHX::EvaluatorDerived<EvalT, TRAITS> {
public:
  
  GatherBasisCoordinates(const panzer::PureBasis & basis);
  
  void postRegistrationSetup(typename TRAITS::SetupData d,
			     PHX::FieldManager<TRAITS>& vm);
  
  void evaluateFields(typename TRAITS::EvalData d);

  static std::string fieldName(const std::string & basisName);

private:
  typedef typename EvalT::ScalarT ScalarT;

  std::string basisName_;
  std::vector<std::string>::size_type basisIndex_; 
  PHX::MDField<ScalarT,Cell,BASIS,Dim> basisCoordinates_;

  GatherBasisCoordinates();
};

}

// **************************************************************
#endif
