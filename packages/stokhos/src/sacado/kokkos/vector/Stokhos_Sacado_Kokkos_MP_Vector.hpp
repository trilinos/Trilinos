// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_SACADO_KOKKOS_MP_VECTOR_HPP
#define STOKHOS_SACADO_KOKKOS_MP_VECTOR_HPP

#include "Stokhos_ConfigDefs.h"

#include "Kokkos_Macros.hpp"

#include "Stokhos_Sacado_Kokkos_MathFunctions.hpp"

#include "Stokhos_KokkosTraits.hpp"
#include "Stokhos_StaticFixedStorage.hpp"
#include "Stokhos_StaticStorage.hpp"
#include "Stokhos_DynamicStorage.hpp"
#include "Stokhos_DynamicStridedStorage.hpp"
#include "Stokhos_DynamicThreadedStorage.hpp"
#include "Stokhos_ViewStorage.hpp"

#include "Sacado_MP_ExpressionTraits.hpp"
#include "Sacado_MP_VectorTraits.hpp"
#include "Sacado_MP_Vector.hpp"
#include "Kokkos_View_MP_Vector.hpp"
#include "Kokkos_Atomic_MP_Vector.hpp"

#include "Teuchos_SerialQRDenseSolver_MP_Vector.hpp"
#include "Teuchos_BLAS_MP_Vector.hpp"
#include "Teuchos_LAPACK_MP_Vector.hpp"

#endif // STOKHOS_SACADO_KOKKOS_MP_VECTOR_HPP
