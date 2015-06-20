// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef SACADO_DISABLE_KOKKOS_CUDA_HPP
#define SACADO_DISABLE_KOKKOS_CUDA_HPP

//
// Include this file in any translation unit to disable the use of Sacado
// classes on Cuda.  Several Sacado classes (e.g., Sacado::Fad::GeneralFad)
// are setup to work with Kokkos, but don't work with Cuda with some choices
// of their template parameters (e.g., Sacado::Fad::MemPoolStorage).  However
// if Cuda is enabled then __device__ is added to the KOKKOS_*_FUNCTION macros
// which prevents these classes from compiling.  By including this file, the
// __device__ annotation will be removed allowing these classes to be compiled
// by NVCC for host code.
//

// Include definitions of KOKKOS_*_FUNCTION macros
#include "Sacado_ConfigDefs.h"

// Redefine KOKKOS_*_FUNCTION macros to not include __device__
#if defined(HAVE_SACADO_KOKKOSCORE) && defined(KOKKOS_HAVE_CUDA)

#undef KOKKOS_FUNCTION
#undef KOKKOS_INLINE_FUNCTION
#undef KOKKOS_FORCEINLINE_FUNCTION

#define KOKKOS_FUNCTION /* */
#define KOKKOS_INLINE_FUNCTION inline
#define KOKKOS_FORCEINLINE_FUNCTION  inline

#endif


#endif // SACADO_DISABLE_KOKKOS_CUDA_HPP
