// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "CellData.hpp"

//**********************************************************************
CellData::CellData()
{ 
  m_coords = Kokkos::View<double***,PHX::Device>("coords",4,4,3);
  m_phi = Kokkos::View<double**,PHX::Device>("phi",4,4);
  m_grad_phi = Kokkos::View<double***,PHX::Device>("grad_phi",4,4,3);

  // just some garbage values for unit testing
  for (PHX::Device::size_type i=0; i < m_phi.extent(0); ++i) {
    for (PHX::Device::size_type j=0; j < m_phi.extent(1); ++j) {
      m_phi(i,j) = 0.25;
      for (PHX::Device::size_type k=0; k < m_phi.extent(2); ++k) {
	m_coords(i,j,k) = 0.25;
	m_grad_phi(i,j,k) = 0.25;
      }
    }
  }
}
//**********************************************************************
Kokkos::View<double***,PHX::Device> CellData::getNodeCoordinates()
{
  return m_coords;
}

//**********************************************************************
Kokkos::View<double**,PHX::Device> CellData::getBasisFunctions()
{
  return m_phi;
}

//**********************************************************************
Kokkos::View<double***,PHX::Device>
CellData::getBasisFunctionGradients()
{
  return m_grad_phi;
}

//**********************************************************************
