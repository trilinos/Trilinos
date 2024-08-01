// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_COORDINATESEVALUTOR_IMPL_HPP
#define PANZER_COORDINATESEVALUTOR_IMPL_HPP

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
CoordinatesEvaluator<EvalT, Traits>::
CoordinatesEvaluator(
  const Teuchos::ParameterList& p) :
  dimension(p.get<int>("Dimension")),
  coordinate( p.get<std::string>("Field Name"), 
	      p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") )
{
  this->addEvaluatedField(coordinate);
  
  std::string n = "CoordinatesEvaluator: " + coordinate.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
CoordinatesEvaluator<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData /* worksets */,
  PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(coordinate,fm);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
CoordinatesEvaluator<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData d)
{ 
  auto coords = this->wda(d).cell_node_coordinates.get_static_view();
  auto coordinate_v = coordinate.get_static_view();
  auto l_dimension = dimension;

  // copy coordinates directly into the field
  Kokkos::parallel_for(d.num_cells, KOKKOS_LAMBDA (int i) {
      for(int j=0;j<coords.extent_int(1);j++)
	coordinate_v(i,j) = coords(i,j,l_dimension);       
    });
  Kokkos::fence();
}

//**********************************************************************

}

#endif
