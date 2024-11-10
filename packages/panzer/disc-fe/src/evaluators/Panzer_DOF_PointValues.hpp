// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_DOF_POINTVALUES_HPP
#define PANZER_EVALUATOR_DOF_POINTVALUES_HPP

#include <string>
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "PanzerDiscFE_config.hpp"
#include "Panzer_BasisValues2.hpp"
#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {
    
//! Interpolates basis DOF values to IP DOF Curl values
template<typename EvalT, typename TRAITS>                   
class DOF_PointValues : public panzer::EvaluatorWithBaseImpl<TRAITS>,      
                        public PHX::EvaluatorDerived<EvalT, TRAITS>  {   
public:

  DOF_PointValues(const Teuchos::ParameterList& p);

  void postRegistrationSetup(typename TRAITS::SetupData d,
                             PHX::FieldManager<TRAITS>& fm);

  void evaluateFields(typename TRAITS::EvalData d);

private:

  typedef typename EvalT::ScalarT ScalarT;
  
  PHX::MDField<const ScalarT,Cell,Point> dof_basis;
  PHX::MDField<ScalarT,Cell,Point> dof_ip_scalar;
  PHX::MDField<ScalarT,Cell,Point,Dim> dof_ip_vector;

  bool is_vector_basis;

  Teuchos::RCP<const PureBasis> basis;
  Teuchos::RCP<BasisValues2<double> > basisValues;
};

/** Interpolates basis DOF values to IP DOF Curl values (specialization for the jacobian)
  * Allows short cut for simple jacobian to dof structure.
  */
template<typename TRAITS>                   
class DOF_PointValues<typename TRAITS::Jacobian,TRAITS> : 
                        public panzer::EvaluatorWithBaseImpl<TRAITS>,      
                        public PHX::EvaluatorDerived<typename TRAITS::Jacobian, TRAITS>  {   
public:

  DOF_PointValues(const Teuchos::ParameterList& p);

  void postRegistrationSetup(typename TRAITS::SetupData d,
                             PHX::FieldManager<TRAITS>& fm);

  void evaluateFields(typename TRAITS::EvalData d);

private:

  typedef panzer::Traits::Jacobian::ScalarT ScalarT;
  
  PHX::MDField<const ScalarT,Cell,Point> dof_basis;
  PHX::MDField<ScalarT,Cell,Point> dof_ip_scalar;
  PHX::MDField<ScalarT,Cell,Point,Dim> dof_ip_vector;

  bool is_vector_basis;

  bool accelerate_jacobian;
  PHX::View<int*> offsets_array;

  Teuchos::RCP<const PureBasis> basis;
  Teuchos::RCP<BasisValues2<double> > basisValues;
};

}

#endif
