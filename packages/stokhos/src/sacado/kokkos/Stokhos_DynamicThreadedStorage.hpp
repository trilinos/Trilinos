// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_DYNAMIC_THREADED_STORAGE_HPP
#define STOKHOS_DYNAMIC_THREADED_STORAGE_HPP

#include "Stokhos_DynArrayTraits.hpp"

#include "Sacado_Traits.hpp"
#include "Stokhos_KokkosTraits.hpp"
#include <sstream>

namespace Stokhos {

  //! Dynamically allocated storage class with striding
  template <typename ordinal_t, typename value_t, typename device_t>
  class DynamicThreadedStorage {};

}

#include "Stokhos_StorageHelpers.hpp"
STOKHOS_STORAGE_HELPER_STRINGNAME_DYNAMIC(DynamicThreadedStorage)

// No Host specialization

// Cuda specialization
#include "Kokkos_Core_fwd.hpp"
#include "Stokhos_DynamicThreadedStorage_cuda.hpp"

#endif // STOKHOS_DYNAMIC_STORAGE_HPP
