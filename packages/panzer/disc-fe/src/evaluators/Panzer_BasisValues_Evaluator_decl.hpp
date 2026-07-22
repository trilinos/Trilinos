// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_BASIS_VALUES_EVALUATOR_DECL_HPP
#define PANZER_BASIS_VALUES_EVALUATOR_DECL_HPP

#include <string>
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "PanzerDiscFE_config.hpp"
#include "Panzer_PointValues2.hpp"
#include "Panzer_BasisValues2.hpp"
#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {
    
//! Interpolates basis DOF values to IP DOF values
template<typename EvalT, typename Traits>
class BasisValues_Evaluator
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    BasisValues_Evaluator(
      const Teuchos::ParameterList& p);

    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;
 
  Teuchos::RCP<const panzer::PureBasis> basis;
  
  // is anything other than ScalarT really needed here?
  Teuchos::RCP<BasisValues2<double> > basisValues;
  PointValues2<double> pointValues;
  PointValues2<const double> constPointValues;

  Teuchos::RCP<const std::vector<Intrepid2::Orientation> > orientations;

  bool derivativesRequired_;
 
  //! Initialization method to unify the constructors.
  void initialize(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                  const Teuchos::RCP<const panzer::PureBasis> & basis,
                  bool derivativesRequired);

public:
  BasisValues_Evaluator(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                        const Teuchos::RCP<const panzer::PureBasis> & basis);

  BasisValues_Evaluator(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                        const Teuchos::RCP<const panzer::PureBasis> & basis,
                        bool derivativesRequired);

}; // end of class BasisValues_Evaluator


}

#endif
