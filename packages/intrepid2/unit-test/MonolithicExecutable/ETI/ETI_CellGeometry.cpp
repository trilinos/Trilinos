// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   ETI_CellGeometry.cpp
    \brief  Explicit Template Instantiation for CellGeometry.  Each definition here should be paired with a declaration in ETI.hpp.
    \author Created by N.V. Roberts.
 */

#include "Intrepid2_CellGeometry.hpp"

#include "Intrepid2_ConfigDefs.hpp"

#include <Kokkos_Core.hpp>

template class Intrepid2::CellGeometry<double,1,Kokkos::DefaultExecutionSpace>;
template class Intrepid2::CellGeometry<double,2,Kokkos::DefaultExecutionSpace>;
template class Intrepid2::CellGeometry<double,3,Kokkos::DefaultExecutionSpace>;

#ifdef HAVE_INTREPID2_SACADO
// CellGeometry - DFad, up to 3D
using Sacado_Fad_DFadType = Sacado::Fad::DFad<double>; // Sacado type used in Intrepid2 tests
template class Intrepid2::CellGeometry<Sacado_Fad_DFadType,1,Kokkos::DefaultExecutionSpace>;
template class Intrepid2::CellGeometry<Sacado_Fad_DFadType,2,Kokkos::DefaultExecutionSpace>;
template class Intrepid2::CellGeometry<Sacado_Fad_DFadType,3,Kokkos::DefaultExecutionSpace>;
#endif
