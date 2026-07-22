// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOKERNELS_EAGER_INITIALIZE_HPP
#define KOKKOKERNELS_EAGER_INITIALIZE_HPP

namespace KokkosKernels {
// \brief Eagerly initialize handles for all enabled TPLs, as well
// as any other globally shared resources that would otherwise be lazily initialized.
//
// Eagerly initializing a TPL means that it doesn't have to be
// lazily initialized when first calling a kernel that uses it.
// For example, \c eager_initialize() will call \c cusparseCreate() upfront
// so that the first call to \c KokkosSparse::spmv doesn't have to.
// This can add a significant amount of apparent runtime to that first kernel
// call, even though the added time isn't really spent in the kernel.
//
// Calling this before using any kernels/TPLs is optional.
// This function is idempotent (any calls after the first have no effect).
//
// \pre \c Kokkos::initialize() has been called.
void eager_initialize();
}  // namespace KokkosKernels

#endif
