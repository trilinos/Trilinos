// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include <Tpetra_KokkosCompat_ClassicNodeAPI_Wrapper.hpp>
#include "KokkosCompat_View.hpp"
#include "KokkosCompat_View_def.hpp"
#endif

// Kokkos-Linalg
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

#include "KokkosBlas_MP_Vector.hpp"

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
template <typename ExecSpace, typename MemSpace>
struct DeviceForNode< Tpetra::KokkosCompat::KokkosDeviceWrapperNode<ExecSpace, MemSpace> > {
  typedef typename Tpetra::KokkosCompat::KokkosDeviceWrapperNode<ExecSpace, MemSpace>::device_type type;
};
#endif

}

#endif // STOKHOS_TPETRA_MP_VECTOR_HPP
