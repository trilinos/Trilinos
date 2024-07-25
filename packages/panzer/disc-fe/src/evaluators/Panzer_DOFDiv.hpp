// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_DOF_DIV_DECL_HPP
#define PANZER_EVALUATOR_DOF_DIV_DECL_HPP

#include "Panzer_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

namespace panzer {

//! Interpolates basis DOF values to IP DOF Div values
template<typename EvalT, typename TRAITS>
class DOFDiv : public panzer::EvaluatorWithBaseImpl<TRAITS>,
	       public PHX::EvaluatorDerived<EvalT, TRAITS>  {
public:

  DOFDiv(const Teuchos::ParameterList& p);

  /** \brief Ctor
    *
    * \param[in] input Tag that corresponds to the input DOF field (sized according to bd)
    * \param[in] output Tag that corresponds to the output field (sized according the id)
    * \param[in] bd Basis being used
    * \param[in] id Integration rule used
    */
  DOFDiv(const PHX::FieldTag & input,
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

  PHX::MDField<const ScalarT,Cell,Point> dof_value;
  PHX::MDField<ScalarT,Cell,IP> dof_div;

  std::string basis_name;
  std::size_t basis_index;
};

// Specitialization for the Jacobian
template<typename TRAITS>
class DOFDiv<panzer::Traits::Jacobian,TRAITS> :
                public panzer::EvaluatorWithBaseImpl<TRAITS>,
                public PHX::EvaluatorDerived<panzer::Traits::Jacobian, TRAITS>  {
public:

  DOFDiv(const Teuchos::ParameterList& p);

  DOFDiv(const PHX::FieldTag & input,
         const PHX::FieldTag & output,
         const panzer::BasisDescriptor & bd,
         const panzer::IntegrationDescriptor & id);

  void postRegistrationSetup(typename TRAITS::SetupData d,
                             PHX::FieldManager<TRAITS>& fm);

  void evaluateFields(typename TRAITS::EvalData d);

private:

  typedef panzer::Traits::Jacobian::ScalarT ScalarT;

  bool use_descriptors_;
  panzer::BasisDescriptor bd_;
  panzer::IntegrationDescriptor id_;

  PHX::MDField<const ScalarT,Cell,Point> dof_value;
  PHX::MDField<ScalarT,Cell,IP> dof_div;

  std::string basis_name;
  std::size_t basis_index;

  bool accelerate_jacobian;
  std::vector<int> offsets;
};

}

#endif
