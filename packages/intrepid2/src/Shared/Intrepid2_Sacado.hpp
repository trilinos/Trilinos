// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_Sacado.hpp
    \brief  Header file to include all Sacado headers that are required if using Intrepid2 with Sacado types.
    \author Created by N. Roberts.
*/

#ifndef __INTREPID2_SACADO_HPP__
#define __INTREPID2_SACADO_HPP__

#include "Intrepid2_ConfigDefs.hpp"

#ifdef HAVE_INTREPID2_SACADO
#include <Sacado.hpp>
#include <Kokkos_DynRankView_Fad.hpp> // This Sacado header defines deep_copy, subview overloads
#endif

#endif
