// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_GATHER_INTEGRATION_COORDINATES_IMPL_HPP
#define PANZER_GATHER_INTEGRATION_COORDINATES_IMPL_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset_Utilities.hpp"

#include "Teuchos_FancyOStream.hpp"

template<typename EvalT,typename TRAITS>
std::string 
panzer::GatherIntegrationCoordinates<EvalT, TRAITS>::
fieldName(int degree)
{
   std::stringstream ss; 
   ss << "IR_" << degree << " IntegrationCoordinates";
   return ss.str();
}

template<typename EvalT,typename TRAITS>
panzer::GatherIntegrationCoordinates<EvalT, TRAITS>::
GatherIntegrationCoordinates(const panzer::IntegrationRule & quad)
{ 
  quadDegree_ = quad.cubature_degree;

  quadCoordinates_ = PHX::MDField<ScalarT,Cell,Point,Dim>(fieldName(quadDegree_),quad.dl_vector);

  this->addEvaluatedField(quadCoordinates_);
  this->addUnsharedField(Teuchos::rcp_const_cast<PHX::FieldTag>(quadCoordinates_.fieldTagPtr()));

  this->setName("GatherIntegrationCoordinates: "+fieldName(quadDegree_));
}

// **********************************************************************
template<typename EvalT,typename TRAITS>
void panzer::GatherIntegrationCoordinates<EvalT, TRAITS>::
postRegistrationSetup(typename TRAITS::SetupData sd, 
		      PHX::FieldManager<TRAITS>& /* fm */)
{
  quadIndex_ = panzer::getIntegrationRuleIndex(quadDegree_, (*sd.worksets_)[0], this->wda);
  // make sure to zero out derivative array as we only set the value for AD types in the loop below
  Kokkos::deep_copy(quadCoordinates_.get_static_view(),0.0);
}

// **********************************************************************
template<typename EvalT,typename TRAITS> 
void panzer::GatherIntegrationCoordinates<EvalT, TRAITS>::
evaluateFields(typename TRAITS::EvalData workset)
{ 
  // const Kokkos::DynRankView<double,PHX::Device> & quadCoords = this->wda(workset).int_rules[quadIndex_]->ip_coordinates;  
  const IntegrationValues2<double> & iv = *this->wda(workset).int_rules[quadIndex_];
  auto s_ip_coordinates = iv.ip_coordinates.get_static_view();
  auto d_quadCoordinates = quadCoordinates_.get_static_view();

  // just copy the array
  Kokkos::MDRangePolicy<PHX::Device,Kokkos::Rank<3>> policy({0,0,0},{int(workset.num_cells),s_ip_coordinates.extent_int(1),s_ip_coordinates.extent_int(2)});
  Kokkos::parallel_for("GatherIntegrationCoords", policy, KOKKOS_LAMBDA (const int i, const int j, const int k) {
    auto s_ip_coordinates_tmp = s_ip_coordinates;
    auto d_quadCoordinates_tmp = d_quadCoordinates;
    if constexpr(Sacado::IsADType<typename EvalT::ScalarT>::value) {
      d_quadCoordinates_tmp(i,j,k).val() = s_ip_coordinates_tmp(i,j,k);
    }
    else {
      d_quadCoordinates_tmp(i,j,k) = s_ip_coordinates_tmp(i,j,k);
    }
  });
}

// **********************************************************************
#endif
