// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MINIEM_PIECEWISECONSTANT_IMPL_HPP
#define MINIEM_PIECEWISECONSTANT_IMPL_HPP

#include <cmath>

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_GatherBasisCoordinates.hpp"

namespace mini_em {

//**********************************************************************
template <typename EvalT,typename Traits>
PiecewiseConstant<EvalT,Traits>::PiecewiseConstant(const std::string & name,
                                                   const panzer::IntegrationRule & ir,
                                                   const panzer::FieldLayoutLibrary & fl,
                                                   const double value0,
                                                   const double value1,
                                                   const double xl,
                                                   const double xr,
                                                   const double yl,
                                                   const double yr,
                                                   const double zl,
                                                   const double zr,
                                                   const std::string& DoF_)
{
  using Teuchos::RCP;

  Teuchos::RCP<PHX::DataLayout> data_layout = ir.dl_scalar;
  ir_degree = ir.cubature_degree;
  ir_dim = ir.spatial_dimension;

  values = PHX::MDField<ScalarT,Cell,Point>(name, data_layout);
  this->addEvaluatedField(values);

  value0_ = value0;
  value1_ = value1;
  xl_ = xl;
  xr_ = xr;
  yl_ = yl;
  yr_ = yr;
  zl_ = zl;
  zr_ = zr;

  Teuchos::RCP<const panzer::PureBasis> basis = fl.lookupBasis(DoF_);
  const std::string coordName = panzer::GatherBasisCoordinates<EvalT,Traits>::fieldName(basis->name());
  coords = PHX::MDField<const ScalarT,Cell,Point,Dim>(coordName, basis->coordinates);
  this->addDependentField(coords);

  std::string n = "PiecewiseConstant";
  this->setName(n);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void PiecewiseConstant<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{
  using panzer::index_t;

  auto tmp_values = values.get_static_view();
  auto tmp_coords = coords.get_static_view();
  const double value0 = value0_;
  const double value1 = value1_;
  const double xl = xl_;
  const double xr = xr_;
  const double yl = yl_;
  const double yr = yr_;
  const double zl = zl_;
  const double zr = zr_;

  if (ir_dim == 3) {
    Kokkos::MDRangePolicy<PHX::exec_space,Kokkos::Rank<2>> policy({0,0},{workset.num_cells,values.extent_int(1)});
    Kokkos::parallel_for("panzer:PiecewiseConstant 3D",policy,KOKKOS_LAMBDA (const int cell,const int point) {

        auto x = tmp_coords(cell,point,0);
        auto y = tmp_coords(cell,point,1);
        auto z = tmp_coords(cell,point,2);

        if ((xl<=x) && (x<=xr) &&
            (yl<=y) && (y<=yr) &&
            (zl<=z) && (z<=zr))
          tmp_values(cell,point) = value0;
        else
          tmp_values(cell,point) = value1;
      });
  } else {
    Kokkos::MDRangePolicy<PHX::exec_space,Kokkos::Rank<2>> policy({0,0},{workset.num_cells,values.extent_int(1)});
    Kokkos::parallel_for("panzer:PiecewiseConstant 2D",policy,KOKKOS_LAMBDA (const int cell,const int point) {

        auto x = tmp_coords(cell,point,0);
        auto y = tmp_coords(cell,point,1);

        if ((xl<=x) && (x<=xr) &&
            (yl<=y) && (y<=yr))
          tmp_values(cell,point) = value0;
        else
          tmp_values(cell,point) = value1;
      });
  }
}

//**********************************************************************
}

#endif
