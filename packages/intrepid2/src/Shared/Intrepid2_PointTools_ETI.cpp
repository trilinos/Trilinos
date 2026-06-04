// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_PointTools_ETI.cpp
    \brief  \brief ETI instantiation for Intrepid2 PointTools class.
    \author Created by N.V. Roberts.
 */

#include "Intrepid2_PointTools.hpp"

INTREPID2_POINTTOOLS_INSTANT(double, Kokkos::DefaultHostExecutionSpace, );
INTREPID2_POINTTOOLS_INSTANT(double, Kokkos::DefaultHostExecutionSpace::device_type, );
INTREPID2_POINTTOOLS_INSTANT(double, Kokkos::LayoutRight COMMA Kokkos::DefaultHostExecutionSpace, );
INTREPID2_POINTTOOLS_INSTANT(double, Kokkos::LayoutRight COMMA Kokkos::DefaultHostExecutionSpace::device_type, );
INTREPID2_POINTTOOLS_INSTANT(double, Kokkos::LayoutRight COMMA Kokkos::HostSpace, );
