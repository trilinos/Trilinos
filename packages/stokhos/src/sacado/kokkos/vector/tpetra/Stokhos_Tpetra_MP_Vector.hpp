// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_TPETRA_MP_VECTOR_HPP
#define STOKHOS_TPETRA_MP_VECTOR_HPP

// This header file should be included whenever compiling any Tpetra
// code with Stokhos scalar types

// MP includes and specializations
#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"

// Kokkos includes
#include "Tpetra_ConfigDefs.hpp"
#include "Kokkos_Core.hpp"
#if defined(HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT)
#include "Kokkos_BufferMacros.hpp"
#include "KokkosCompat_ClassicNodeAPI_Wrapper.hpp"
#include "KokkosCompat_View.hpp"
#include "KokkosCompat_View_def.hpp"
#endif

// Kokkos-Linalg
#if defined(TPETRA_HAVE_KOKKOS_REFACTOR)
#include "Kokkos_ArithTraits_MP_Vector.hpp"
#include "Kokkos_InnerProductSpaceTraits_MP_Vector.hpp"
#include "Kokkos_MV_MP_Vector.hpp"

// Order may be important here because Cuda provides a slightly more
// specialized version of MV_Multiply() for a single column
#include "Kokkos_CrsMatrix_MP_Vector_Cuda.hpp"
#include "Kokkos_CrsMatrix_MP_Vector.hpp"

#include "Kokkos_TeuchosCommAdapters_MP_Vector.hpp"
#include "Tpetra_KokkosRefactor_Details_MultiVectorDistObjectKernels_MP_Vector.hpp"
#include "Tpetra_KokkosRefactor_Details_MultiVectorLocalDeepCopy_MP_Vector.hpp"
#include "Tpetra_Details_fill_MP_Vector.hpp"
#include "Kokkos_Random_MP_Vector.hpp"
#endif

namespace Stokhos {

/// \brief Trait class that determines (new) Kokkos execution space
///   type from Kokkos(Classic) Node type.
/// \tparam Node Kokkos(Classic) Node type.
///
/// The \c type typedef of this class gives the (new) Kokkos execution
/// space corresponding to the given (classic) Kokkos \c Node type.
template <typename Node>
struct DeviceForNode {
#if defined(KOKKOS_ENABLE_SERIAL)
  // Prefer the Kokkos::Serial execution space if it exists.
  typedef Kokkos::Serial type;
#else
  // If the Kokkos::Serial execution space does not exist, use the
  // default host execution space.  Kokkos::HostSpace (the host memory
  // space) always exists, and it always has an execution_space
  // typedef, which corresponds to the default host execution space.
  typedef Kokkos::HostSpace::execution_space type;
#endif // defined(KOKKOS_ENABLE_SERIAL)
};

#if defined(HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT)
/// \brief Partial specialization of DeviceForNode, for new Kokkos
///   "wrapper" Node types.
/// \tparam Device (New) Kokkos execution space type.
template <typename Device>
struct DeviceForNode< Kokkos::Compat::KokkosDeviceWrapperNode<Device> > {
  typedef Device type;
};
#endif

}

#endif // STOKHOS_TPETRA_MP_VECTOR_HPP
