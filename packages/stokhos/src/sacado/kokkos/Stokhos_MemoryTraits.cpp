// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <new>

#include "Stokhos_MemoryTraits.hpp"

// When aligning memory fo the MIC architecture, we need to globally replace
// new/delete with aligned versions because of seg faults in the MueLu setup.
// This is true even though MP::Vector provides its own aligned new/delete
// overloads, as well as aligned STL allocators.  Somewhere in the MueLu setup
// code, a raw new/delete must be called.
//
// Note, this still doesn't resolve seg faults for MIC when using (Intel) MPI.

#if STOKHOS_ALIGN_MEMORY && defined(__MIC__)
void* operator new(std::size_t count) throw (std::bad_alloc) {
  return Stokhos::MemoryTraits< Kokkos::HostSpace >::alloc(count);
}
void* operator new[](std::size_t count) throw (std::bad_alloc) {
  return Stokhos::MemoryTraits< Kokkos::HostSpace >::alloc(count);
}
void operator delete(void *ptr) throw() {
  Stokhos::MemoryTraits< Kokkos::HostSpace >::free(ptr);
}
void operator delete[](void *ptr) throw() {
  Stokhos::MemoryTraits< Kokkos::HostSpace >::free(ptr);
}
#endif
