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
  PHX::MDField<const double, BASIS, IP,    void, void, void, void, void, void>
    constBasisRefScalar_;
  PHX::MDField<const double, Cell,  BASIS, IP,   void, void, void, void, void>
    constBasisScalar_;
  PHX::MDField<const double, BASIS, IP,    Dim,  void, void, void, void, void>
    constBasisRefVector_;
  PHX::MDField<const double, Cell,  BASIS, IP,   Dim,  void, void, void, void>
    constBasisVector_;
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
  Kokkos::View<int*,PHX::Device> offsets_array;

  Teuchos::RCP<const PureBasis> basis;
  Teuchos::RCP<BasisValues2<double> > basisValues;
  PHX::MDField<const double, BASIS, IP,    void, void, void, void, void, void>
    constBasisRefScalar_;
  PHX::MDField<const double, Cell,  BASIS, IP,   void, void, void, void, void>
    constBasisScalar_;
  PHX::MDField<const double, BASIS, IP,    Dim,  void, void, void, void, void>
    constBasisRefVector_;
  PHX::MDField<const double, Cell,  BASIS, IP,   Dim,  void, void, void, void>
    constBasisVector_;
};

}

#endif
