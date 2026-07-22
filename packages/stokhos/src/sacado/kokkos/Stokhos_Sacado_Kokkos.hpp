// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_SACADO_KOKKOS_HPP
#define STOKHOS_SACADO_KOKKOS_HPP

#include "Stokhos_ConfigDefs.h"

#ifdef HAVE_STOKHOS_ENSEMBLE_SCALAR_TYPE
#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"
#endif

#ifdef HAVE_STOKHOS_PCE_SCALAR_TYPE
#include "Stokhos_Sacado_Kokkos_UQ_PCE.hpp"
#endif

#endif // STOKHOS_SACADO_KOKKOS_HPP
