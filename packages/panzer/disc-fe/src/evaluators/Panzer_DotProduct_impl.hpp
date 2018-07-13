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

#ifndef PANZER_EVALUATOR_DotProduct_IMPL_HPP
#define PANZER_EVALUATOR_DotProduct_IMPL_HPP

#include <string>

#include "Panzer_PointRule.hpp"
#include "Panzer_Workset_Utilities.hpp"

#include "Teuchos_RCP.hpp"

namespace panzer {

template <typename EvalT,typename TraitsT>
Teuchos::RCP<DotProduct<EvalT,TraitsT> > 
buildEvaluator_DotProduct(const std::string & resultName,
                          const panzer::PointRule & pr,
                          const std::string & vecA,
                          const std::string & vecB,
                          double multiplier,
                          const std::string & fieldMultiplier)
{
  Teuchos::ParameterList pl;
  pl.set("Result Name",resultName);
  pl.set("Point Rule",Teuchos::rcpFromRef(pr));
  pl.set("Vector A Name",vecA);
  pl.set("Vector B Name",vecB);
  pl.set("Multiplier",multiplier);
  pl.set("Field Multiplier",fieldMultiplier);

  return Teuchos::rcp(new DotProduct<EvalT,TraitsT>(pl));
}

//**********************************************************************
template<typename EvalT, typename Traits>
DotProduct<EvalT, Traits>::
DotProduct(
  const Teuchos::ParameterList& p)
  : multiplier_field_on(false)
{
  std::string result_name = p.get<std::string>("Result Name");
  std::string vec_a_name = p.get<std::string>("Vector A Name");
  std::string vec_b_name = p.get<std::string>("Vector B Name");

  std::string multiplier_name = "";
  if(p.isType<std::string>("Field Multiplier"))
    multiplier_name = p.get<std::string>("Field Multiplier");

  multiplier_value = 1.0;
  if(p.isType<double>("Multiplier"))
    multiplier_value = p.get<double>("Multiplier");
  
  const Teuchos::RCP<const panzer::PointRule> pr = 
    p.get< Teuchos::RCP<const panzer::PointRule> >("Point Rule");

  vec_a_dot_vec_b = PHX::MDField<ScalarT>(result_name, pr->dl_scalar);
  vec_a = PHX::MDField<const ScalarT>(vec_a_name, pr->dl_vector);
  vec_b = PHX::MDField<const ScalarT>(vec_b_name, pr->dl_vector);

  if(multiplier_name!="") {
    multiplier_field = PHX::MDField<const ScalarT>(multiplier_name,pr->dl_scalar);
    multiplier_field_on = true;
    this->addDependentField(multiplier_field);
  }

  this->addEvaluatedField(vec_a_dot_vec_b);
  this->addDependentField(vec_a);
  this->addDependentField(vec_b);
 
  std::string n = "DotProduct: " + result_name + " = " + vec_a_name + " . " + vec_b_name;
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
DotProduct<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData  /* sd */,
  PHX::FieldManager<Traits>&  fm)
{
  this->utils.setFieldData(vec_a_dot_vec_b,fm);
  this->utils.setFieldData(vec_a,fm);
  this->utils.setFieldData(vec_b,fm);

  if(multiplier_field_on)
    this->utils.setFieldData(multiplier_field,fm);

  num_pts = vec_a.extent(1);
  num_dim = vec_a.extent(2);

  TEUCHOS_ASSERT(vec_a.extent(1) == vec_b.extent(1));
  TEUCHOS_ASSERT(vec_a.extent(2) == vec_b.extent(2));
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
DotProduct<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ 
  for (index_t cell = 0; cell < workset.num_cells; ++cell) {
    for (int p = 0; p < num_pts; ++p) {
      vec_a_dot_vec_b(cell,p) = ScalarT(0.0);
      for (int dim = 0; dim < num_dim; ++dim)
	vec_a_dot_vec_b(cell,p) += vec_a(cell,p,dim) * vec_b(cell,p,dim); 

      if(multiplier_field_on)
        vec_a_dot_vec_b(cell,p) *= multiplier_value*multiplier_field(cell,p);
      else
        vec_a_dot_vec_b(cell,p) *= multiplier_value;
    }
  }
}

//**********************************************************************

}

#endif
