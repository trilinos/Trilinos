/*
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
*/

#ifndef KOKKOSSPARSE_IMPL_GMRES_SPEC_HPP_
#define KOKKOSSPARSE_IMPL_GMRES_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_BsrMatrix.hpp"
#include "KokkosKernels_Handle.hpp"

// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosSparse_gmres_impl.hpp>
#endif

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class KernelHandle, class AT, class AO, class AD, class AM, class AS, class BType, class XType>
struct gmres_eti_spec_avail {
  enum : bool { value = false };
};

}  // namespace Impl
}  // namespace KokkosSparse

#define KOKKOSSPARSE_GMRES_ETI_SPEC_AVAIL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE,  \
                                          MEM_SPACE_TYPE)                                                        \
  template <>                                                                                                    \
  struct gmres_eti_spec_avail<                                                                                   \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      const SCALAR_TYPE, const ORDINAL_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                    \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE,                                                \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> > > {                          \
    enum : bool { value = true };                                                                                \
  };

// Include the actual specialization declarations
#include <KokkosSparse_gmres_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosSparse_gmres_eti_spec_avail.hpp>

namespace KokkosSparse {
namespace Impl {

// Unification layer
/// \brief Implementation of KokkosSparse::gmres

template <class KernelHandle, class AT, class AO, class AD, class AM, class AS, class BType, class XType,
          bool tpl_spec_avail = gmres_tpl_spec_avail<KernelHandle, AT, AO, AD, AM, AS, BType, XType>::value,
          bool eti_spec_avail = gmres_eti_spec_avail<KernelHandle, AT, AO, AD, AM, AS, BType, XType>::value>
struct GMRES {
  using AMatrix  = CrsMatrix<AT, AO, AD, AM, AS>;
  using BAMatrix = KokkosSparse::Experimental::BsrMatrix<AT, AO, AD, AM, AS>;
  static void gmres(KernelHandle *handle, const AMatrix &A, const BType &B, XType &X,
                    KokkosSparse::Experimental::Preconditioner<AMatrix> *precond = nullptr);

  static void gmres(KernelHandle *handle, const BAMatrix &A, const BType &B, XType &X,
                    KokkosSparse::Experimental::Preconditioner<BAMatrix> *precond = nullptr);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of gmres
// Unification layer
template <class KernelHandle, class AT, class AO, class AD, class AM, class AS, class BType, class XType>
struct GMRES<KernelHandle, AT, AO, AD, AM, AS, BType, XType, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  using AMatrix = CrsMatrix<AT, AO, AD, AM, AS>;
  static void gmres(KernelHandle *handle, const AMatrix &A, const BType &B, XType &X,
                    KokkosSparse::Experimental::Preconditioner<AMatrix> *precond = nullptr) {
    auto gmres_handle = handle->get_gmres_handle();
    using Gmres       = Experimental::GmresWrap<typename std::remove_pointer<decltype(gmres_handle)>::type>;

    Gmres::gmres(*gmres_handle, A, B, X, precond);
  }

  using BAMatrix = KokkosSparse::Experimental::BsrMatrix<AT, AO, AD, AM, AS>;
  static void gmres(KernelHandle *handle, const BAMatrix &A, const BType &B, XType &X,
                    KokkosSparse::Experimental::Preconditioner<BAMatrix> *precond = nullptr) {
    auto gmres_handle = handle->get_gmres_handle();
    using Gmres       = Experimental::GmresWrap<typename std::remove_pointer<decltype(gmres_handle)>::type>;

    Gmres::gmres(*gmres_handle, A, B, X, precond);
  }
};

#endif
}  // namespace Impl
}  // namespace KokkosSparse

//
// Macro for declaration of full specialization of
// This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSSPARSE_GMRES_ETI_SPEC_DECL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE,   \
                                         MEM_SPACE_TYPE)                                                         \
  extern template struct GMRES<                                                                                  \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      const SCALAR_TYPE, const ORDINAL_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                    \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE,                                                \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      false, true>;

#define KOKKOSSPARSE_GMRES_ETI_SPEC_INST(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE,   \
                                         MEM_SPACE_TYPE)                                                         \
  template struct GMRES<                                                                                         \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      const SCALAR_TYPE, const ORDINAL_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                    \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE,                                                \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >,                             \
      false, true>;

#include <KokkosSparse_gmres_tpl_spec_decl.hpp>

#endif
