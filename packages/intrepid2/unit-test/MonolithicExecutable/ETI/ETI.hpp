// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   ETI.hpp
    \brief  Declarations for explicit template instantiations.  Each of these should have a corresponding definition in an ETI*.cpp file in the same directory.
    \author Created by N.V. Roberts.
 */

#include "Intrepid2_Sacado.hpp" // Sacado includes, guarded by the appropriate preprocessor variable

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_CellGeometry.hpp"
#include "Intrepid2_HierarchicalBasisFamily.hpp"

#include <Kokkos_Core.hpp>

#ifndef Intrepid2_ETI_h
#define Intrepid2_ETI_h

//#include "ETI_FunctionSpaceTools_declare.hpp"

// CellGeometry - double, up to 3D
extern template class Intrepid2::CellGeometry<double,1,Kokkos::DefaultExecutionSpace>;
extern template class Intrepid2::CellGeometry<double,2,Kokkos::DefaultExecutionSpace>;
extern template class Intrepid2::CellGeometry<double,3,Kokkos::DefaultExecutionSpace>;

#ifdef HAVE_INTREPID2_SACADO
// CellGeometry - DFad, up to 3D
extern template class Intrepid2::CellGeometry<Sacado::Fad::DFad<double>,1,Kokkos::DefaultExecutionSpace>;
extern template class Intrepid2::CellGeometry<Sacado::Fad::DFad<double>,2,Kokkos::DefaultExecutionSpace>;
extern template class Intrepid2::CellGeometry<Sacado::Fad::DFad<double>,3,Kokkos::DefaultExecutionSpace>;
#endif

// TODO: figure out the best way to declare these… (they are template aliases)
//extern template class Intrepid2::HierarchicalBasisFamily  <Kokkos::DefaultExecutionSpace,double,double>;
//extern template class Intrepid2::DGHierarchicalBasisFamily<Kokkos::DefaultExecutionSpace,double,double>;
//
//#ifdef HAVE_INTREPID2_SACADO
//extern template class Intrepid2::HierarchicalBasisFamily  <Kokkos::DefaultExecutionSpace, Sacado::Fad::DFad<double>, Sacado::Fad::DFad<double> >;
//extern template class Intrepid2::DGHierarchicalBasisFamily<Kokkos::DefaultExecutionSpace, Sacado::Fad::DFad<double>, Sacado::Fad::DFad<double> >;
//#endif


#endif /* Intrepid2_ETI_h */
