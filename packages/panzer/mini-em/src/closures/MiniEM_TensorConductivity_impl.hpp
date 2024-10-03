// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MINIEM_TENSORCONDUCTIVITY_IMPL_HPP
#define MINIEM_TENSORCONDUCTIVITY_IMPL_HPP

#include <cmath>

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_GatherBasisCoordinates.hpp"

namespace mini_em {

//**********************************************************************
template <typename EvalT,typename Traits>
TensorConductivity<EvalT,Traits>::TensorConductivity(const std::string & name,
                                                     const panzer::IntegrationRule & ir,
                                                     const panzer::FieldLayoutLibrary & fl,
                                                     const double & sigma_,
                                                     const double & betax_,
                                                     const double & betay_,
                                                     const double & betaz_,
                                                     const std::string& DoF_)
{
  using Teuchos::RCP;

  Teuchos::RCP<PHX::DataLayout> data_layout = ir.dl_tensor;
  ir_degree = ir.cubature_degree;
  ir_dim = ir.spatial_dimension;

  conductivity = PHX::MDField<ScalarT,Cell,Point,Dim,Dim>(name, data_layout);
  this->addEvaluatedField(conductivity);

  betax = betax_;
  betay = betay_;
  betaz = betaz_;
  sigma = sigma_ / (1.0 + betax*betax + betay*betay + betaz*betaz);

  Teuchos::RCP<const panzer::PureBasis> basis = fl.lookupBasis(DoF_);
  const std::string coordName = panzer::GatherBasisCoordinates<EvalT,Traits>::fieldName(basis->name());
  coords = PHX::MDField<const ScalarT,Cell,Point,Dim>(coordName, basis->coordinates);
  this->addDependentField(coords);

  std::string n = "TensorConductivity";
  this->setName(n);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void TensorConductivity<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{
  using panzer::index_t;
  if (ir_dim == 3) {
    auto temp_conductivity = conductivity;
    auto temp_sigma = sigma;
    auto temp_betax = betax;
    auto temp_betay = betay;
    auto temp_betaz = betaz;
    Kokkos::MDRangePolicy<PHX::exec_space,Kokkos::Rank<2>> policy({0,0},{workset.num_cells,conductivity.extent_int(1)});
    Kokkos::parallel_for("panzer:TensorConductivity 3D",policy,KOKKOS_LAMBDA (const int cell,const int point) {
        // const ScalarT& x = coords(cell,point,0);
        // const ScalarT& y = coords(cell,point,1);
        // const ScalarT& z = coords(cell,point,2);
        temp_conductivity(cell,point,0,0) = temp_sigma * (1.0 + temp_betax*temp_betax);
        temp_conductivity(cell,point,0,1) = temp_sigma * (      temp_betax*temp_betay - temp_betaz);
        temp_conductivity(cell,point,0,2) = temp_sigma * (      temp_betax*temp_betaz + temp_betay);

        temp_conductivity(cell,point,1,0) = temp_sigma * (      temp_betay*temp_betax + temp_betaz);
        temp_conductivity(cell,point,1,1) = temp_sigma * (1.0 + temp_betay*temp_betay);
        temp_conductivity(cell,point,1,2) = temp_sigma * (      temp_betay*temp_betaz - temp_betax);

        temp_conductivity(cell,point,2,0) = temp_sigma * (      temp_betaz*temp_betax - temp_betay);
        temp_conductivity(cell,point,2,1) = temp_sigma * (      temp_betaz*temp_betay + temp_betax);
        temp_conductivity(cell,point,2,2) = temp_sigma * (1.0 + temp_betaz*temp_betaz);
      });
  } else {
    auto temp_conductivity = conductivity;
    auto temp_sigma = sigma;
    Kokkos::MDRangePolicy<PHX::exec_space,Kokkos::Rank<2>> policy({0,0},{workset.num_cells,conductivity.extent_int(1)});
    Kokkos::parallel_for("panzer:TensorConductivity 2D",policy,KOKKOS_LAMBDA (const int cell,const int point) {
        // const ScalarT& x = coords(cell,point,0);
        // const ScalarT& y = coords(cell,point,1);
        temp_conductivity(cell,point,0,0) = temp_sigma;
        temp_conductivity(cell,point,0,1) = 0.;

        temp_conductivity(cell,point,1,0) = 0.;
        temp_conductivity(cell,point,1,1) = temp_sigma;
      });
  }
}

//**********************************************************************
}

#endif
