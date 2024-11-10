// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_CellExtreme_impl_hpp__
#define __Panzer_CellExtreme_impl_hpp__

#include <limits>

#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
CellExtreme<EvalT, Traits>::
CellExtreme(
  const Teuchos::ParameterList& p) : quad_index(0)
{
  Teuchos::RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  // setup default value
  use_max = true;
  if(p.isType<bool>("Use Max")) 
    use_max = p.get<bool>("Use Max");

  Teuchos::RCP<panzer::IntegrationRule> ir = p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR");
  quad_order = ir->cubature_degree;

  Teuchos::RCP<PHX::DataLayout> dl_cell = Teuchos::rcp(new PHX::MDALayout<Cell>(ir->dl_scalar->extent(0)));
  extreme = PHX::MDField<ScalarT>( p.get<std::string>("Extreme Name"), dl_cell);
  scalar = PHX::MDField<const ScalarT,Cell,IP>( p.get<std::string>("Field Name"), ir->dl_scalar);

  this->addEvaluatedField(extreme);
  this->addDependentField(scalar);
    
  multiplier = 1.0;
  if(p.isType<double>("Multiplier"))
     multiplier = p.get<double>("Multiplier");

  if (p.isType<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers")) {
    const std::vector<std::string>& field_multiplier_names = 
      *(p.get<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers"));

    for (std::vector<std::string>::const_iterator name = 
	   field_multiplier_names.begin(); 
	 name != field_multiplier_names.end(); ++name) {
      PHX::MDField<const ScalarT,Cell,IP> tmp_field(*name, p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar);
      field_multipliers.push_back(tmp_field);
    }
  }

  for (typename std::vector<PHX::MDField<const ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
       field != field_multipliers.end(); ++field)
    this->addDependentField(*field);

  std::string n = "CellExtreme: " + extreme.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
CellExtreme<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  num_qp = scalar.extent(1);
  quad_index =  panzer::getIntegrationRuleIndex(quad_order,(*sd.worksets_)[0], this->wda);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
CellExtreme<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
    double mult = this->multiplier;
    std::size_t num_pt = this->num_qp;
    bool extreme_max = this->use_max;
    auto scalar_view = scalar.get_view();
    auto extreme_view = extreme.get_view();

    // compute field weighted with multipliers
    PHX::View<ScalarT**> mult_field("Multiplied Field", workset.num_cells, scalar.extent(1));

    // initialize to scalar value
    Kokkos::parallel_for("Initialize Muliplier Field", workset.num_cells, KOKKOS_LAMBDA( const int cell) {
      for (std::size_t qp = 0; qp < num_pt; ++qp) {
        mult_field(cell, qp) = mult * scalar_view(cell, qp);
      }
    });

    // multiply by field values
    for (std::size_t field_num = 0; field_num < field_multipliers.size(); ++field_num)
    {
      auto field = field_multipliers[field_num];
      Kokkos::parallel_for("CellExtreme: Multiply Fields", workset.num_cells, KOKKOS_LAMBDA( const int cell) {
        for (std::size_t qp = 0; qp < num_pt; ++qp) {
          mult_field(cell, qp) *= field(cell, qp);
        }
      });
    }

    // take extreme over points in each cell
    if (extreme_max) {
      Kokkos::parallel_for ("CellExtreme (max)", workset.num_cells, KOKKOS_LAMBDA( const int cell) {
        for (std::size_t qp = 0; qp < num_pt; ++qp) {
          auto current = mult_field(cell, qp);

          // take first point
          if (qp == 0)
            extreme_view(cell) = current;
          else
            extreme_view(cell) = extreme_view(cell)<current ? current : extreme_view(cell);
        }
      });
    } else {
      Kokkos::parallel_for ("CellExtreme (min)", workset.num_cells, KOKKOS_LAMBDA( const int cell) {
        for (std::size_t qp = 0; qp < num_pt; ++qp) {
          auto current = mult_field(cell, qp);
        
          // take first point
          if (qp == 0)
            extreme_view(cell) = current;
          else // use_min
            extreme_view(cell) = extreme_view(cell)>current ? current : extreme_view(cell);
        }
      });
    }
}

//**********************************************************************
template<typename EvalT, typename TRAITS>
Teuchos::RCP<Teuchos::ParameterList> 
CellExtreme<EvalT, TRAITS>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Extreme Name", "?");
  p->set<std::string>("Field Name", "?");

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);
  p->set<bool>("Use Max", true);
  p->set<double>("Multiplier", 1.0);

  Teuchos::RCP<const std::vector<std::string> > fms;
  p->set("Field Multipliers", fms);
  return p;
}

//**********************************************************************

}

#endif
