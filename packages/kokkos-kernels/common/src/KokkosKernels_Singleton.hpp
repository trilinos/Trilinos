//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOSKERNELS_SINGLETON_HPP
#define KOKKOSKERNELS_SINGLETON_HPP

#include <Kokkos_Core.hpp>

namespace KokkosKernels {
namespace Impl {

// Singleton structure for a default-constructible type.
// If initialized, the object will be destructed only during Kokkos::finalize.
// This is safe to use as a global/static variable type.
//   This is unlike std::unique_ptr, whose destructor frees the object but doesn't
//   set its internal pointer to null. This can cause a double-free error when Kokkos::finalize
//   tries to free the same object later.
template <typename T>
struct Singleton {
  // Get the underlying singleton object.
  // If it hasn't been constructed yet, lazily construct it.
  T& get() {
    if (!ptr) {
      ptr = new T();
      Kokkos::push_finalize_hook([this]() {
        delete this->ptr;
        this->ptr = nullptr;
      });
    }
    return *ptr;
  }

  bool is_initialized() const { return ptr != nullptr; }

  T* ptr = nullptr;
};

}  // namespace Impl
}  // namespace KokkosKernels

#endif
