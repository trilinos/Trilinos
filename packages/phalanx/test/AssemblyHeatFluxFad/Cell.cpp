// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Cell.hpp"

//**********************************************************************
MyCell::MyCell()
{ 
  m_phi = Kokkos::View<double**,PHX::Device>("phi",4,4);
  m_grad_phi = Kokkos::View<double***,PHX::Device>("grad_phi",4,4,3);

  // just some garbage values for unit testing
  // for (PHX::Device::size_type i=0; i < m_phi.dimension(0); ++i) {
  //   for (PHX::Device::size_type j=0; j < m_phi.dimension(1); ++j) {
  //     m_phi(i,j) = 0.25;
  //     for (PHX::Device::size_type k=0; k < m_phi.dimension(2); ++k) {
  //       m_grad_phi(i,j,k) = 0.25;
  //     }
  //   }
  // }
  Kokkos::deep_copy(m_phi,0.25);
  Kokkos::deep_copy(m_grad_phi,0.25);
}

//**********************************************************************
Kokkos::View<double**,PHX::Device>
MyCell::getNodeCoordinates()
{
  return m_coords;
}

//**********************************************************************
Kokkos::View<double**,PHX::Device>
MyCell::getBasisFunctions() 
{
  return m_phi;
}

//**********************************************************************
Kokkos::View<double***,PHX::Device>
MyCell::getBasisFunctionGradients()
{
  return m_grad_phi;
}

//**********************************************************************
std::size_t MyCell::localIndex()
{
  return local_index_;
}

//**********************************************************************
void MyCell::setLocalIndex(std::size_t index)
{
  local_index_ = index;
}

//**********************************************************************
