// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_EXAMPLE_MYCELL_HPP
#define PHX_EXAMPLE_MYCELL_HPP

#include "Phalanx_config.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"

class MyCell {
  
public:

  MyCell();
  
  virtual ~MyCell() {}
  
  Kokkos::View<double**,PHX::Device> getNodeCoordinates();
  
  Kokkos::View<double**,PHX::Device> getBasisFunctions();
  
  Kokkos::View<double***,PHX::Device> getBasisFunctionGradients();
  
  std::size_t localIndex();

  void setLocalIndex(std::size_t index);

private:
  
  std::size_t local_index_;

  Kokkos::View<double**,PHX::Device> m_coords;
  Kokkos::View<double**,PHX::Device> m_phi;
  Kokkos::View<double***,PHX::Device> m_grad_phi;

};

#endif
