// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_CONSTANT_FLUX_IMPL_HPP
#define PANZER_CONSTANT_FLUX_IMPL_HPP

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
ConstantFlux<EvalT, Traits>::
ConstantFlux(
  const Teuchos::ParameterList& p) :
  flux( p.get<std::string>("Flux Field Name"), 
	p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") )
{
  const Teuchos::ParameterList& flux_values = p.sublist("Flux Values");

  for (Teuchos::ParameterList::ConstIterator i = flux_values.begin(); i != flux_values.end(); ++i)
    values.push_back(Teuchos::getValue<double>(i->second));

  this->addEvaluatedField(flux);
  
  std::string n = "ConstantFlux: " + flux.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
ConstantFlux<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData  /* worksets */,
  PHX::FieldManager<Traits>&  fm)
{
  using namespace PHX;
  this->utils.setFieldData(flux,fm);

  TEUCHOS_ASSERT(static_cast<std::size_t>(flux.extent(2)) == values.size());

  auto flux_v = flux.get_static_view();
  
  for (int dim = 0; dim < flux_v.extent_int(2); ++dim) {
    auto val = values[dim];
    Kokkos::parallel_for ("ConstantFlux", flux.extent_int(0), KOKKOS_LAMBDA( const int cell) {
	for (int ip = 0; ip < flux_v.extent_int(1); ++ip)
	  flux_v(cell,ip,dim) = val;
      });
  }
  Kokkos::fence();
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
ConstantFlux<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  /* d */)
{ }

//**********************************************************************

}

#endif
