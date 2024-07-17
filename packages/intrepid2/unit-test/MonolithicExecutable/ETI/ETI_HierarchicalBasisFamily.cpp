// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   ETI_HierarchicalBasisFamily.cpp
    \brief  Explicit Template Instantiation for HierarchicalBasisFamily.  Each definition here should be paired with a declaration in ETI.hpp.
    \author Created by N.V. Roberts.
 */

#include "Intrepid2_HierarchicalBasisFamily.hpp"

#include "Intrepid2_ConfigDefs.hpp"

#include <Kokkos_Core.hpp>

// TODO: figure out the best way to define these… (they are template aliases)
//template class Intrepid2::HierarchicalBasisFamily  <Kokkos::DefaultExecutionSpace,double,double>;
//template class Intrepid2::DGHierarchicalBasisFamily<Kokkos::DefaultExecutionSpace,double,double>;
//
//#ifdef HAVE_INTREPID2_SACADO
//using Sacado_Fad_DFadType = Sacado::Fad::DFad<double>; // Sacado type used in Intrepid2 tests
//template class Intrepid2::HierarchicalBasisFamily  <Kokkos::DefaultExecutionSpace, Sacado::Fad::DFad<double>, Sacado::Fad::DFad<double> >;
//template class Intrepid2::DGHierarchicalBasisFamily<Kokkos::DefaultExecutionSpace, Sacado::Fad::DFad<double>, Sacado::Fad::DFad<double> >;
//#endif
