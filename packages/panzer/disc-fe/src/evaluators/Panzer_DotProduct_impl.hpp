// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
DotProduct<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ 

  auto vec_a_v = vec_a.get_static_view();
  auto vec_b_v = vec_b.get_static_view();
  auto vec_a_dot_vec_b_v = vec_a_dot_vec_b.get_static_view();
  auto multiplier_field_v = multiplier_field.get_static_view();

  int l_num_pts = num_pts, l_num_dim = num_dim;
  auto l_multiplier_field_on = multiplier_field_on;
  auto l_multiplier_value = multiplier_value;
  
  Kokkos::parallel_for (workset.num_cells, KOKKOS_LAMBDA (const int cell) {
      for (int p = 0; p < l_num_pts; ++p) {
	vec_a_dot_vec_b_v(cell,p) = ScalarT(0.0);
	for (int dim = 0; dim < l_num_dim; ++dim)
	  vec_a_dot_vec_b_v(cell,p) += vec_a_v(cell,p,dim) * vec_b_v(cell,p,dim); 
	
	if(l_multiplier_field_on)
	  vec_a_dot_vec_b_v(cell,p) *= l_multiplier_value*multiplier_field_v(cell,p);
	else
	  vec_a_dot_vec_b_v(cell,p) *= l_multiplier_value;
      }
    });
  Kokkos::fence();
}

//**********************************************************************

}

#endif
