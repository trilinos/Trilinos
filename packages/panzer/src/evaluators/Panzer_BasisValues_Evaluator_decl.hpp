#ifndef PANZER_BASIS_VALUES_EVALUATOR_DECL_HPP
#define PANZER_BASIS_VALUES_EVALUATOR_DECL_HPP

#include <string>
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_Field.hpp"

#include "Panzer_config.hpp"
#include "Panzer_PointValues.hpp"
#include "Panzer_BasisValues.hpp"

namespace panzer {
    
//! Interpolates basis DOF values to IP DOF values
PHX_EVALUATOR_CLASS(BasisValues_Evaluator)
 
  Teuchos::RCP<const panzer::PureBasis> basis;
  
  // is anything other than ScalarT really needed here?
  BasisValues<ScalarT,PHX::MDField<ScalarT> > basisValues;
  PointValues<ScalarT,PHX::MDField<ScalarT> > pointValues;
 
  //! Initialization method to unify the constructors.
  void initialize(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                  const Teuchos::RCP<const panzer::PureBasis> & basis);

public:
  BasisValues_Evaluator(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                        const Teuchos::RCP<const panzer::PureBasis> & basis);

PHX_EVALUATOR_CLASS_END

}

#endif
