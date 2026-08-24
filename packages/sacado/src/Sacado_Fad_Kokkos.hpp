// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_KOKKOS_HPP
#define SACADO_FAD_KOKKOS_HPP

#include "Sacado_ConfigDefs.h"

#if defined(HAVE_SACADO_KOKKOS)

#if defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

#include "Kokkos_Core_fwd.hpp"
#include "Sacado_Fad_Kokkos_View_Support.hpp"

#else

// Bring in some definitions when the Fad specialization is not defined
// (these need to be moved out of this file)
#include "KokkosExp_View_Fad_Contiguous.hpp"

#endif // defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

#endif // defined(HAVE_SACADO_KOKKOS)

#endif /* #ifndef SACADO_FAD_KOKKOS_HPP */
