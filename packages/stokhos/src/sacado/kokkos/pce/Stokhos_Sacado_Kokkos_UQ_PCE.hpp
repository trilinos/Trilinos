// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_SACADO_KOKKOS_UQ_PCE_HPP
#define STOKHOS_SACADO_KOKKOS_UQ_PCE_HPP

#include "Stokhos_ConfigDefs.h"

#if defined(Stokhos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The PCE scalar type is deprecated."
#endif
#endif

#include "Kokkos_Macros.hpp"

#include "Stokhos_Sacado_Kokkos_MathFunctions.hpp"

#include "Stokhos_KokkosTraits.hpp"
#include "Stokhos_StaticFixedStorage.hpp"
#include "Stokhos_StaticStorage.hpp"
#include "Stokhos_DynamicStorage.hpp"
#include "Stokhos_DynamicStridedStorage.hpp"
#include "Stokhos_DynamicThreadedStorage.hpp"
#include "Stokhos_ViewStorage.hpp"

// Don't include "Sacado_UQ_PCE.hpp" here, it is included in the file below
// in a special order for Kokkos overloads to be found correctly
#include "Kokkos_View_UQ_PCE.hpp"
#include "Kokkos_Atomic_UQ_PCE.hpp"

#include "Teuchos_SerialQRDenseSolver_UQ_PCE.hpp"
#include "Teuchos_LAPACK_UQ_PCE.hpp"

#endif // STOKHOS_SACADO_KOKKOS_MP_VECTOR_HPP
