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

#ifndef PANZER_CONVECTION_T_HPP
#define PANZER_CONVECTION_T_HPP

#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"

namespace user_app {

//**********************************************************************
template<typename EvalT, typename Traits>
Convection<EvalT, Traits>::
Convection(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;

  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");

  RCP<DataLayout> scalar = ir->dl_scalar;
  RCP<DataLayout> vector = ir->dl_vector;
  
  conv = MDField<ScalarT,Cell,Point>(p.get<string>("Operator Name"),scalar);

  a = MDField<const ScalarT,Cell,Point,Dim>(p.get<string>("A Name"),vector);

  grad_x = 
    MDField<const ScalarT,Cell,Point,Dim>(p.get<string>("Gradient Name"),vector);

  multiplier = p.get<double>("Multiplier");

  this->addEvaluatedField(conv);

  this->addDependentField(a);
  this->addDependentField(grad_x);
  
  std::string n = "Convection: " + conv.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
Convection<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ 
  typedef typename PHX::MDField<ScalarT,Cell,Point>::size_type size_type;
  
  for (panzer::index_t cell = 0; cell < workset.num_cells; ++cell) {    
    for (size_type point = 0; point < conv.extent(1); ++point) {
      
      conv(cell,point) = 0.0;
	
      for (size_type dim = 0; dim < a.extent(2); ++dim)
	conv(cell,point) += a(cell,point,dim) * grad_x(cell,point,dim);
      
      conv(cell,point) *= multiplier;
      
    }
  }

}

//**********************************************************************

}

#endif
