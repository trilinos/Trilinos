// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  typename Traits::SetupData /* sd */,
  PHX::FieldManager<Traits>& /* fm */)
{
  num_pts = vec_a.extent(1);
  num_dim = vec_a.extent(2);

  TEUCHOS_ASSERT(vec_a.extent(1) == vec_b.extent(1));
  TEUCHOS_ASSERT(vec_a.extent(2) == vec_b.extent(2));
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
