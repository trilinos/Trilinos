#ifndef PANZER_EVALUATOR_DOF_DECL_HPP
#define PANZER_EVALUATOR_DOF_DECL_HPP

#include <string>
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_Field.hpp"

#include "Panzer_config.hpp"
#include "Panzer_BasisValues.hpp"

namespace panzer {
    
//! Interpolates basis DOF values to IP DOF values
PHX_EVALUATOR_CLASS(DOF)
  
  PHX::MDField<ScalarT,Cell,Point> dof_basis;
  PHX::MDField<ScalarT> dof_ip;

  std::string basis_name;
  std::size_t basis_index;

  PHX::MDField<ScalarT,Cell,BASIS> dof_orientation;
  bool requires_orientation;

PHX_EVALUATOR_CLASS_END

//! Interpolates basis DOF values to IP DOF values
PHX_EVALUATOR_CLASS(DOF_PointValues)
  
  PHX::MDField<ScalarT,Cell,Point> dof_basis;
  PHX::MDField<ScalarT> dof_ip;

  PHX::MDField<ScalarT,Cell,BASIS> dof_orientation;
  bool requires_orientation;

  Teuchos::RCP<const PureBasis> basis;
  BasisValues<ScalarT,PHX::MDField<ScalarT> > basisValues;

PHX_EVALUATOR_CLASS_END

}

#endif
