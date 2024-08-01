// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_EXAMPLE_CELL_DATA_HPP
#define PHX_EXAMPLE_CELL_DATA_HPP

#include <vector>
#include "Phalanx_config.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include <vector>

class CellData {
  
public:

  CellData();
  
  virtual ~CellData() {}
    
  Kokkos::View<double***,PHX::Device> getNodeCoordinates();
  
  Kokkos::View<double**,PHX::Device> getBasisFunctions();
  
  Kokkos::View<double***,PHX::Device> getBasisFunctionGradients();
  
private:
  
  Kokkos::View<double***,PHX::Device> m_coords;
  Kokkos::View<double**,PHX::Device> m_phi;
  Kokkos::View<double***,PHX::Device> m_grad_phi;
};

#endif
