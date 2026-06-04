// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_CellTools_ETI.cpp
    \brief  \brief ETI instantiation for Intrepid2 CellTools class.
    \author Created by N.V. Roberts.
 */

#include "Intrepid2_CellTools.hpp"

INTREPID2_CELLTOOLS_INSTANT(Kokkos::DefaultExecutionSpace::device_type, double, double, );
