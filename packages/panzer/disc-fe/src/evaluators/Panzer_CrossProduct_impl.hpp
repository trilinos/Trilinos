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

#ifndef PANZER_EVALUATOR_CrossProduct_IMPL_HPP
#define PANZER_EVALUATOR_CrossProduct_IMPL_HPP

#include <string>

#include "Panzer_PointRule.hpp"
#include "Panzer_Workset_Utilities.hpp"

#include "Teuchos_RCP.hpp"

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
CrossProduct<EvalT, Traits>::
CrossProduct(
  const Teuchos::ParameterList& p)
{
  std::string result_name = p.get<std::string>("Result Name");
  std::string vec_a_name = p.get<std::string>("Vector A Name");
  std::string vec_b_name = p.get<std::string>("Vector B Name");
  
  const Teuchos::RCP<const panzer::PointRule> pr = 
    p.get< Teuchos::RCP<const panzer::PointRule> >("Point Rule");

  // use a scalar field only if dimension is 2D
  useScalarField = (pr->spatial_dimension==2);

  if(!useScalarField)
    vec_a_cross_vec_b = PHX::MDField<ScalarT>(result_name, pr->dl_vector);
  else 
    vec_a_cross_vec_b = PHX::MDField<ScalarT>(result_name, pr->dl_scalar);

  vec_a = PHX::MDField<const ScalarT>(vec_a_name, pr->dl_vector);
  vec_b = PHX::MDField<const ScalarT>(vec_b_name, pr->dl_vector);

  this->addEvaluatedField(vec_a_cross_vec_b);
  this->addDependentField(vec_a);
  this->addDependentField(vec_b);
 
  std::string n = "CrossProduct: " + result_name + " = " + vec_a_name + " . " + vec_b_name;
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
CrossProduct<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData  /* sd */,
  PHX::FieldManager<Traits>&  fm)
{
  this->utils.setFieldData(vec_a_cross_vec_b,fm);
  this->utils.setFieldData(vec_a,fm);
  this->utils.setFieldData(vec_b,fm);

  num_pts = vec_a.dimension(1);
  num_dim = vec_a.dimension(2);

  TEUCHOS_ASSERT(vec_a.dimension(1) == vec_b.dimension(1));
  TEUCHOS_ASSERT(vec_a.dimension(2) == vec_b.dimension(2));
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
CrossProduct<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ 
  if(useScalarField) {
    for (index_t cell = 0; cell < workset.num_cells; ++cell) {
      for (int p = 0; p < num_pts; ++p) {
        vec_a_cross_vec_b(cell,p) = vec_a(cell,p,0)*vec_b(cell,p,1)-vec_a(cell,p,1)*vec_b(cell,p,0);
      }
    }
  }
  else {
    for (index_t cell = 0; cell < workset.num_cells; ++cell) {
      for (int p = 0; p < num_pts; ++p) {
        vec_a_cross_vec_b(cell,p,0) =   vec_a(cell,p,1)*vec_b(cell,p,2)-vec_a(cell,p,2)*vec_b(cell,p,1);
        vec_a_cross_vec_b(cell,p,1) = -(vec_a(cell,p,0)*vec_b(cell,p,2)-vec_a(cell,p,2)*vec_b(cell,p,0));
        vec_a_cross_vec_b(cell,p,2) =   vec_a(cell,p,0)*vec_b(cell,p,1)-vec_a(cell,p,1)*vec_b(cell,p,0);
      }
    }
  }
}

//**********************************************************************

}

#endif
