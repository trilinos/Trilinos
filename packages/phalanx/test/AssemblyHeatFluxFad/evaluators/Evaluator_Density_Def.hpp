// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//**********************************************************************
template<typename EvalT, typename Traits> Density<EvalT, Traits>::
Density(const Teuchos::ParameterList& p) :
  density("Density", p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") ),
  temp("Temperature", p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") )
{
  this->addEvaluatedField(density);

  // Should normally use a const dependent field. We will use a
  // non-const to test the method addNonConstDependentField.
  // this->addDependentField(temp);
  this->addNonConstDependentField(temp);
  this->setName("Density");
}

//**********************************************************************
template<typename EvalT, typename Traits>
void Density<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData /* d */,
		      PHX::FieldManager<Traits>& /* vm */)
{
  cell_data_size = density.size() / density.dimension(0);
}

//**********************************************************************
template<typename EvalT, typename Traits>
KOKKOS_INLINE_FUNCTION
void Density<EvalT, Traits>:: operator () (const int i) const
{
  for (PHX::index_size_type ip=0; ip< static_cast<PHX::index_size_type>(density.extent(1)); ip++)
    density(i,ip) =  temp(i,ip) * temp(i,ip);
}

//**********************************************************************
// template<typename EvalT, typename Traits>
// KOKKOS_INLINE_FUNCTION
// void Density<EvalT, Traits>:: operator () (const DensityTag, const int i) const
// {
//   for (PHX::index_size_type ip=0; ip< static_cast<PHX::index_size_type>(density.extent(1)); ip++)
//     density(i,ip) =  temp(i,ip) * temp(i,ip);
// }

//**********************************************************************
// template<typename EvalT, typename Traits>
// KOKKOS_INLINE_FUNCTION
// void Density<EvalT, Traits>:: operator () (const DensityTag, typename Kokkos::TeamPolicy<>::member_type & team) const
// {
//   for (PHX::index_size_type ip=0; ip< static_cast<PHX::index_size_type>(density.extent(1)); ip++)
//     density(0,ip) =  temp(0,ip) * temp(0,ip);
// }

//*********************************************************************
template<typename EvalT, typename Traits>
void Density<EvalT, Traits>::evaluateFields(typename Traits::EvalData d)
{
  // typedef Kokkos::TeamPolicy<DensityTag> team_policy ;
  // team_policy policy(d.num_cells,2);
  // Kokkos::parallel_for(policy, *this);

  Kokkos::parallel_for(d.num_cells, *this);
}

//**********************************************************************
