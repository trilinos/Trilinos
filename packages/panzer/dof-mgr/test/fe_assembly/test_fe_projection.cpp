// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
//
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_01.cpp
    \brief  Test for checking orientation tools for hexahedral elements.
    \author Created by Mauro Perego
*/

#include "Kokkos_Core.hpp"
#include "test_fe_projection.hpp"


int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiScope (&argc, &argv);
  Kokkos::ScopeGuard kokkosScope (argc, argv);
  
  const int r_val = Discretization::Example::feProjection<double,PHX::Device>(argc, argv);

  return r_val;
}

