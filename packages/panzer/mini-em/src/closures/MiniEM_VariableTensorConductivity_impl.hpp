// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MINIEM_VARIABLETENSORCONDUCTIVITY_IMPL_HPP
#define MINIEM_VARIABLETENSORCONDUCTIVITY_IMPL_HPP

#include <cmath>

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_GatherBasisCoordinates.hpp"

namespace mini_em {

//**********************************************************************
template <typename EvalT,typename Traits>
VariableTensorConductivity<EvalT,Traits>::VariableTensorConductivity(const std::string & name,
                                                                     const panzer::IntegrationRule & ir,
                                                                     const panzer::FieldLayoutLibrary & fl,
                                                                     const double & sigma0_,
                                                                     const double & sigma1_,
                                                                     const double & sigma2_,
                                                                     const double & betax0_,
                                                                     const double & betay0_,
                                                                     const double & betaz0_,
                                                                     const double & betax1_,
                                                                     const double & betay1_,
                                                                     const double & betaz1_,
                                                                     const double & betax2_,
                                                                     const double & betay2_,
                                                                     const double & betaz2_,
                                                                     const std::string& DoF_)
{
  using Teuchos::RCP;

  Teuchos::RCP<PHX::DataLayout> data_layout = ir.dl_tensor;
  ir_degree = ir.cubature_degree;
  ir_dim = ir.spatial_dimension;

  conductivity = PHX::MDField<ScalarT,Cell,Point,Dim,Dim>(name, data_layout);
  this->addEvaluatedField(conductivity);

  betax0 = betax0_;
  betay0 = betay0_;
  betaz0 = betaz0_;
  sigma0 = sigma0_ / (1.0 + betax0*betax0 + betay0*betay0 + betaz0*betaz0);

  betax1 = betax1_;
  betay1 = betay1_;
  betaz1 = betaz1_;
  sigma1 = sigma1_ / (1.0 + betax1*betax1 + betay1*betay1 + betaz1*betaz1);

  betax2 = betax2_;
  betay2 = betay2_;
  betaz2 = betaz2_;
  sigma2 = sigma2_ / (1.0 + betax2*betax2 + betay2*betay2 + betaz2*betaz2);

  Teuchos::RCP<const panzer::PureBasis> basis = fl.lookupBasis(DoF_);
  const std::string coordName = panzer::GatherBasisCoordinates<EvalT,Traits>::fieldName(basis->name());
  coords = PHX::MDField<const ScalarT,Cell,Point,Dim>(coordName, basis->coordinates);
  this->addDependentField(coords);

  std::string n = "VariableTensorConductivity";
  this->setName(n);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void VariableTensorConductivity<EvalT,Traits>::postRegistrationSetup(typename Traits::SetupData sd,
                                                                     PHX::FieldManager<Traits>& /* fm */)
{
  ir_index = panzer::getIntegrationRuleIndex(ir_degree,(*sd.worksets_)[0], this->wda);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void VariableTensorConductivity<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{
  using panzer::index_t;

  const double xl0 = 0.4;
  const double xr0 = 0.6;
  const double yl0 = 0.4;
  const double yr0 = 0.6;
  const double zl0 = 0.4;
  const double zr0 = 0.6;

  const double xl1 = 0.2;
  const double xr1 = 0.8;
  const double yl1 = 0.2;
  const double yr1 = 0.8;
  const double zl1 = 0.2;
  const double zr1 = 0.8;

  if (ir_dim == 3) {
    for (index_t cell = 0; cell < workset.num_cells; ++cell) {
      for (int point = 0; point < conductivity.extent_int(1); ++point) {

        auto x = workset.int_rules[ir_index]->ip_coordinates(cell,point,0);
        auto y = workset.int_rules[ir_index]->ip_coordinates(cell,point,1);
        auto z = workset.int_rules[ir_index]->ip_coordinates(cell,point,2);

        if ((xl0<=x) && (x<=xr0) &&
            (yl0<=y) && (y<=yr0) &&
            (zl0<=z) && (z<=zr0)) {

          conductivity(cell,point,0,0) = sigma0 * (1.0 + betax0*betax0);
          conductivity(cell,point,0,1) = sigma0 * (      betax0*betay0 - betaz0);
          conductivity(cell,point,0,2) = sigma0 * (      betax0*betaz0 + betay0);

          conductivity(cell,point,1,0) = sigma0 * (      betay0*betax0 + betaz0);
          conductivity(cell,point,1,1) = sigma0 * (1.0 + betay0*betay0);
          conductivity(cell,point,1,2) = sigma0 * (      betay0*betaz0 - betax0);

          conductivity(cell,point,2,0) = sigma0 * (      betaz0*betax0 - betay0);
          conductivity(cell,point,2,1) = sigma0 * (      betaz0*betay0 + betax0);
          conductivity(cell,point,2,2) = sigma0 * (1.0 + betaz0*betaz0);

        } else if ((xl1<=x) && (x<=xr1) &&
                   (yl1<=y) && (y<=yr1) &&
                   (zl1<=z) && (z<=zr1)) {

          conductivity(cell,point,0,0) = sigma1 * (1.0 + betax1*betax1);
          conductivity(cell,point,0,1) = sigma1 * (      betax1*betay1 - betaz1);
          conductivity(cell,point,0,2) = sigma1 * (      betax1*betaz1 + betay1);

          conductivity(cell,point,1,0) = sigma1 * (      betay1*betax1 + betaz1);
          conductivity(cell,point,1,1) = sigma1 * (1.0 + betay1*betay1);
          conductivity(cell,point,1,2) = sigma1 * (      betay1*betaz1 - betax1);

          conductivity(cell,point,2,0) = sigma1 * (      betaz1*betax1 - betay1);
          conductivity(cell,point,2,1) = sigma1 * (      betaz1*betay1 + betax1);
          conductivity(cell,point,2,2) = sigma1 * (1.0 + betaz1*betaz1);

        } else {

          conductivity(cell,point,0,0) = sigma2 * (1.0 + betax2*betax2);
          conductivity(cell,point,0,1) = sigma2 * (      betax2*betay2 - betaz2);
          conductivity(cell,point,0,2) = sigma2 * (      betax2*betaz2 + betay2);

          conductivity(cell,point,1,0) = sigma2 * (      betay2*betax2 + betaz2);
          conductivity(cell,point,1,1) = sigma2 * (1.0 + betay2*betay2);
          conductivity(cell,point,1,2) = sigma2 * (      betay2*betaz2 - betax2);

          conductivity(cell,point,2,0) = sigma2 * (      betaz2*betax2 - betay2);
          conductivity(cell,point,2,1) = sigma2 * (      betaz2*betay2 + betax2);
          conductivity(cell,point,2,2) = sigma2 * (1.0 + betaz2*betaz2);

        }
      }
    }
  } else {
    for (index_t cell = 0; cell < workset.num_cells; ++cell) {
      for (int point = 0; point < conductivity.extent_int(1); ++point) {

        auto x = workset.int_rules[ir_index]->ip_coordinates(cell,point,0);
        auto y = workset.int_rules[ir_index]->ip_coordinates(cell,point,1);

        if ((xl0<=x) && (x<=xr0) &&
            (yl0<=y) && (y<=yr0)) {

          conductivity(cell,point,0,0) = sigma0;
          conductivity(cell,point,0,1) = 0.;

          conductivity(cell,point,1,0) = 0.;
          conductivity(cell,point,1,1) = sigma0;

        } else if ((xl1<=x) && (x<=xr1) &&
                   (yl1<=y) && (y<=yr1)) {

          conductivity(cell,point,0,0) = sigma1;
          conductivity(cell,point,0,1) = 0.;

          conductivity(cell,point,1,0) = 0.;
          conductivity(cell,point,1,1) = sigma1;

        } else {

          conductivity(cell,point,0,0) = sigma2;
          conductivity(cell,point,0,1) = 0.;

          conductivity(cell,point,1,0) = 0.;
          conductivity(cell,point,1,1) = sigma2;

        }
      }
    }
  }
}

//**********************************************************************
}

#endif
