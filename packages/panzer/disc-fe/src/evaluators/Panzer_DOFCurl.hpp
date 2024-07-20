// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_DOF_CURL_DECL_HPP
#define PANZER_EVALUATOR_DOF_CURL_DECL_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {
    
//! Interpolates basis DOF values to IP DOF Curl values
template<typename EvalT, typename TRAITS>                   
class DOFCurl : public panzer::EvaluatorWithBaseImpl<TRAITS>,      
                public PHX::EvaluatorDerived<EvalT, TRAITS>  {   
public:

  DOFCurl(const Teuchos::ParameterList& p);

  /** \brief Ctor
    *
    * \param[in] input Tag that corresponds to the input DOF field (sized according to bd)
    * \param[in] output Tag that corresponds to the output field (sized according to bd and the basis_dim)
    * \param[in] bd Basis being differentiated
    * \param[in] id Integration rule used
    * \param[in] basis_dim Spatial dimension, used to choose the rank of the output field
    */
  DOFCurl(const PHX::FieldTag & input,
          const PHX::FieldTag & output,
          const panzer::BasisDescriptor & bd,
          const panzer::IntegrationDescriptor & id,
          int basis_dim);

  void postRegistrationSetup(typename TRAITS::SetupData d,
                             PHX::FieldManager<TRAITS>& fm);

  void evaluateFields(typename TRAITS::EvalData d);

private:

  typedef typename EvalT::ScalarT ScalarT;

  bool use_descriptors_;
  panzer::BasisDescriptor bd_;
  panzer::IntegrationDescriptor id_;

  PHX::MDField<const ScalarT,Cell,Point> dof_value;

  PHX::MDField<ScalarT,Cell,Point> dof_curl_scalar;
  PHX::MDField<ScalarT,Cell,Point,Dim> dof_curl_vector;

  std::string basis_name;
  std::size_t basis_index;
  int basis_dimension;
};

// Specitialization for the Jacobian
template<typename TRAITS>                   
class DOFCurl<typename TRAITS::Jacobian,TRAITS> : 
                public panzer::EvaluatorWithBaseImpl<TRAITS>,      
                public PHX::EvaluatorDerived<typename TRAITS::Jacobian, TRAITS>  {   
public:

  DOFCurl(const Teuchos::ParameterList& p);

  DOFCurl(const PHX::FieldTag & input,
          const PHX::FieldTag & output,
          const panzer::BasisDescriptor & bd,
          const panzer::IntegrationDescriptor & id,
          int basis_dim);

  void postRegistrationSetup(typename TRAITS::SetupData d,
                             PHX::FieldManager<TRAITS>& fm);

  void evaluateFields(typename TRAITS::EvalData d);

private:

  typedef panzer::Traits::Jacobian::ScalarT ScalarT;

  bool use_descriptors_;
  panzer::BasisDescriptor bd_;
  panzer::IntegrationDescriptor id_;

  PHX::MDField<const ScalarT,Cell,Point> dof_value;

  PHX::MDField<ScalarT,Cell,Point> dof_curl_scalar;
  PHX::MDField<ScalarT,Cell,Point,Dim> dof_curl_vector;

  PHX::View<const int*> offsets_array;
  std::vector<int> offsets;

  std::string basis_name;
  std::size_t basis_index;
  int basis_dimension;

  bool accelerate_jacobian;
};



}

#endif
