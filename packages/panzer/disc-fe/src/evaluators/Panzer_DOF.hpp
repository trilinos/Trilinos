// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_DOF_DECL_HPP
#define PANZER_EVALUATOR_DOF_DECL_HPP

#include <string>
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "PanzerDiscFE_config.hpp"
#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {
    
//! Interpolates basis DOF values to IP DOF values
template<typename EvalT, typename TRAITS>                   
class DOF : public panzer::EvaluatorWithBaseImpl<TRAITS>,      
            public PHX::EvaluatorDerived<EvalT, TRAITS>  {   
public:

  DOF(const Teuchos::ParameterList& p);

  /** \brief Ctor
    *
    * \param[in] input Tag that corresponds to the input DOF field (sized according to bd)
    * \param[in] output Tag that corresponds to the output field (sized according the id and if bd corresponds to a vector basis)
    * \param[in] bd Basis being used
    * \param[in] id Integration rule used
    */
  DOF(const PHX::FieldTag & input,
      const PHX::FieldTag & output,
      const panzer::BasisDescriptor & bd,
      const panzer::IntegrationDescriptor & id);

  void postRegistrationSetup(typename TRAITS::SetupData d,
                             PHX::FieldManager<TRAITS>& fm);

  void evaluateFields(typename TRAITS::EvalData d);

private:

  typedef typename EvalT::ScalarT ScalarT;

  bool use_descriptors_;
  panzer::BasisDescriptor bd_;
  panzer::IntegrationDescriptor id_;
  
  PHX::MDField<const ScalarT,Cell,Point> dof_basis;

  PHX::MDField<ScalarT,Cell,Point> dof_ip_scalar;
  PHX::MDField<ScalarT,Cell,Point,Dim> dof_ip_vector;

  std::string basis_name;
  std::size_t basis_index;

  bool is_vector_basis;
};

/** Interpolates basis DOF values to IP DOF Curl values (specialization for the jacobian)
  * Allows short cut for simple jacobian to dof structure.
  */
template<typename TRAITS>                   
class DOF<typename TRAITS::Jacobian,TRAITS> : 
            public panzer::EvaluatorWithBaseImpl<TRAITS>,      
            public PHX::EvaluatorDerived<typename TRAITS::Jacobian, TRAITS>  {   
public:

  DOF(const Teuchos::ParameterList& p);

  DOF(const PHX::FieldTag & input,
      const PHX::FieldTag & output,
      const panzer::BasisDescriptor & bd,
      const panzer::IntegrationDescriptor & id);

  void postRegistrationSetup(typename TRAITS::SetupData d,
                             PHX::FieldManager<TRAITS>& fm);

  void preEvaluate(typename TRAITS::PreEvalData d);

  void evaluateFields(typename TRAITS::EvalData d);

private:

  typedef panzer::Traits::Jacobian::ScalarT ScalarT;

  bool use_descriptors_;
  panzer::BasisDescriptor bd_;
  panzer::IntegrationDescriptor id_;
  
  PHX::MDField<const ScalarT,Cell,Point> dof_basis;

  PHX::MDField<ScalarT,Cell,Point> dof_ip_scalar;
  PHX::MDField<ScalarT,Cell,Point,Dim> dof_ip_vector;

  std::string basis_name;
  std::size_t basis_index;

  bool accelerate_jacobian_enabled;
  bool accelerate_jacobian;
  PHX::View<int*> offsets_array;
  std::string sensitivities_name; // This sets which gather operations have sensitivities
                                  // and thus which DOF operations can use accelerated jacobians

  bool is_vector_basis;
};

}

#endif
