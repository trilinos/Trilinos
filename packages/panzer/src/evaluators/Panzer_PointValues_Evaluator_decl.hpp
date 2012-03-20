#ifndef PANZER_POINT_VALUES_EVALUATOR_DECL_HPP
#define PANZER_POINT_VALUES_EVALUATOR_DECL_HPP

#include <string>
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_Field.hpp"

#include "Panzer_PointValues.hpp"

namespace panzer {
    
//! Interpolates basis DOF values to IP DOF values
PHX_EVALUATOR_CLASS(PointValues_Evaluator)

  // is anything other than ScalarT really needed here?
  PointValues<ScalarT,PHX::MDField<ScalarT> > pointValues;
 
  Intrepid::FieldContainer<double> refPointArray;

  //! Initialization method to unify the constructors.
  void initialize(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                  const Teuchos::RCP<const Intrepid::FieldContainer<double> > & userArray);

public:
  PointValues_Evaluator(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                        const Teuchos::RCP<const Intrepid::FieldContainer<double> > & userArray);

PHX_EVALUATOR_CLASS_END

}

#endif
