// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_DeviceAssert.hpp
    \brief  Implementation of an assert that can safely be called from device code.
    \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_DeviceAssert_h
#define Intrepid2_DeviceAssert_h

#include <Kokkos_Core_fwd.hpp>
#include <impl/Kokkos_Error.hpp>

#include <cstdio>
#include <cassert>

namespace Intrepid2 {
#ifdef NDEBUG
#define device_assert( v ) ((void)0)
#else

  // assert that can reasonably be called in device code
  KOKKOS_INLINE_FUNCTION
  void device_assert(bool val) {
    if (!val)
      /// KK: this assert header does not need to be separated from utils.hpp
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      Kokkos::abort("ASSERT IN CUDA CALL, SHOULD ABORT\n");
#else
    assert(false);
#endif
  }
#endif

}


#endif /* Interpid2_DeviceAssert_h */
