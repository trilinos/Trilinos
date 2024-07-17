// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_DISABLE_KOKKOS_CUDA_HPP
#define SACADO_DISABLE_KOKKOS_CUDA_HPP

//
// Include this file in any translation unit to disable the use of Sacado
// classes on Cuda.  Several Sacado classes (e.g., Sacado::Fad::GeneralFad)
// are setup to work with Kokkos, but don't work with Cuda with some choices
// of their template parameters.  However
// if Cuda is enabled then __device__ is added to the KOKKOS_*_FUNCTION macros
// which prevents these classes from compiling.  By including this file, the
// __device__ annotation will be removed allowing these classes to be compiled
// by NVCC for host code.
//

// Include definitions of KOKKOS_*_FUNCTION macros
#include "Sacado_ConfigDefs.h"

// Redefine KOKKOS_*_FUNCTION macros to not include __device__
#if defined(HAVE_SACADO_KOKKOS) &&  ( defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) )
// TODO double check me
#if 1
#undef SACADO_FUNCTION
#undef SACADO_INLINE_FUNCTION
#undef SACADO_FORCEINLINE_FUNCTION

#define SACADO_FUNCTION /* */
#define SACADO_INLINE_FUNCTION inline
#define SACADO_FORCEINLINE_FUNCTION  inline
#endif
#define SACADO_DISABLE_CUDA_IN_KOKKOS 1

#endif


#endif // SACADO_NO_KOKKOS_HPP
