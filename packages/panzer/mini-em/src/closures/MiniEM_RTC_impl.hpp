// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MINIEM_RTC_CONDUCTIVITY_IMPL_HPP
#define MINIEM_RTC_CONDUCTIVITY_IMPL_HPP

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_GatherBasisCoordinates.hpp"

namespace mini_em {

//**********************************************************************
template <typename EvalT,typename Traits>
RTC<EvalT,Traits>::RTC(const std::string & name,
                       const panzer::IntegrationRule & ir,
                       const panzer::FieldLayoutLibrary & fl,
                       const std::string &funBody,
                       const std::string& DoF_)
{
  using Teuchos::RCP;

  Teuchos::RCP<PHX::DataLayout> data_layout = ir.dl_scalar;
  ir_degree = ir.cubature_degree;
  ir_dim = ir.spatial_dimension;

  values = PHX::MDField<ScalarT,Cell,Point>(name, data_layout);
  this->addEvaluatedField(values);

  Teuchos::RCP<const panzer::PureBasis> basis = fl.lookupBasis(DoF_);
  const std::string coordName = panzer::GatherBasisCoordinates<EvalT,Traits>::fieldName(basis->name());
  coords = PHX::MDField<const ScalarT,Cell,Point,Dim>(coordName, basis->coordinates);
  this->addDependentField(coords);

  fun_ = PG_RuntimeCompiler::Function();
  if (ir_dim == 3) {
    fun_.addVar("double","x");
    fun_.addVar("double","y");
    fun_.addVar("double","z");
  } else {
    fun_.addVar("double","x");
    fun_.addVar("double","y");
  }
  fun_.addVar("double","value");
  fun_.addBody(funBody);

  std::string n = "RTC: "+name;
  this->setName(n);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void RTC<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{
  using panzer::index_t;

  double x, y, z, value;
  if (ir_dim == 3) {
    fun_.varAddrFill(0, &x);
    fun_.varAddrFill(1, &y);
    fun_.varAddrFill(2, &z);
    fun_.varAddrFill(3, &value);
    for (index_t cell = 0; cell < workset.num_cells; ++cell) {
      for (int point = 0; point < values.extent_int(1); ++point) {

        if constexpr(Sacado::IsADType<typename EvalT::ScalarT>::value) {
          x = coords(cell,point,0).val();
          y = coords(cell,point,1).val();
          z = coords(cell,point,2).val();
        } else {
          x = coords(cell,point,0);
          y = coords(cell,point,1);
          z = coords(cell,point,2);
        }

        fun_.execute();
        values(cell,point) = value;
      }
    }
  } else {
    fun_.varAddrFill(0, &x);
    fun_.varAddrFill(1, &y);
    fun_.varAddrFill(2, &value);
    for (index_t cell = 0; cell < workset.num_cells; ++cell) {
      for (int point = 0; point < values.extent_int(1); ++point) {

        if constexpr(Sacado::IsADType<typename EvalT::ScalarT>::value) {
          x = coords(cell,point,0).val();
          y = coords(cell,point,1).val();
        } else {
          x = coords(cell,point,0);
          y = coords(cell,point,1);
        }

        fun_.execute();

        values(cell,point) = value;
      }
    }
  }
}

//**********************************************************************
}

#endif
