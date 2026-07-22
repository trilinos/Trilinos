// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_DOF_GRADIENT_DECL_HPP
#define PANZER_EVALUATOR_DOF_GRADIENT_DECL_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {
    
//! Interpolates basis DOF values to IP DOF Gradient values
template <typename EvalT, typename TRAITS>                   
class DOFGradient : public panzer::EvaluatorWithBaseImpl<TRAITS>,      
                    public PHX::EvaluatorDerived<EvalT, TRAITS>  {   
public:

  DOFGradient(const Teuchos::ParameterList& p);

  /** \brief Ctor
    *
    * \param[in] input Tag that corresponds to the input DOF field (sized according to bd)
    * \param[in] output Tag that corresponds to the output field (sized according the id, and the dimension)
    * \param[in] bd Basis being used
    * \param[in] id Integration rule used
    */
  DOFGradient(const PHX::FieldTag & input,
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

  // <cell,point>
  PHX::MDField<const ScalarT,Cell,Point> dof_value;
  // <cell,point,dim>
  PHX::MDField<ScalarT> dof_gradient;

  std::string basis_name;
  std::size_t basis_index;
};

}

#endif
