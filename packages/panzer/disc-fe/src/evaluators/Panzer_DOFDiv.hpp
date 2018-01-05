// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_DOF_DIV_DECL_HPP
#define PANZER_EVALUATOR_DOF_DIV_DECL_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Panzer_Evaluator_WithBaseImpl.hpp"

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
  DOFDiv(const PHX::Tag<typename EvalT::ScalarT> & input,
         const PHX::Tag<typename EvalT::ScalarT> & output,
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

  DOFDiv(const PHX::Tag<panzer::Traits::Jacobian::ScalarT> & input,
         const PHX::Tag<panzer::Traits::Jacobian::ScalarT> & output,
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
