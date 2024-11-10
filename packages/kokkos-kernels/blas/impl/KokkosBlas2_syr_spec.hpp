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

#ifndef KOKKOSBLAS2_SYR_SPEC_HPP_
#define KOKKOSBLAS2_SYR_SPEC_HPP_

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosBlas2_syr_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class EXEC_SPACE, class XMV, class ZMV>
struct syr_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::SYR. This is NOT for users!!! All the declarations of full
// specializations go in this header file. We may spread out definitions (see
// _INST macro below) across one or more .cpp files.
//
#define KOKKOSBLAS2_SYR_ETI_SPEC_AVAIL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                          \
  template <>                                                                                          \
  struct syr_eti_spec_avail<EXEC_SPACE,                                                                \
                            Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                    \
                            Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,      \
                                         Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                 \
    enum : bool { value = true };                                                                      \
  };

// Include the actual specialization declarations
#include <KokkosBlas2_syr_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas2_syr_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

//
// syr
//

// Implementation of KokkosBlas::syr.
template <class ExecutionSpace, class XViewType, class AViewType,
          bool tpl_spec_avail = syr_tpl_spec_avail<ExecutionSpace, XViewType, AViewType>::value,
          bool eti_spec_avail = syr_eti_spec_avail<ExecutionSpace, XViewType, AViewType>::value>
struct SYR {
  static void syr(const ExecutionSpace& space, const char trans[], const char uplo[],
                  const typename AViewType::const_value_type& alpha, const XViewType& x, const AViewType& A)
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
  {
    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY ? "KokkosBlas::syr[ETI]"
                                                                     : "KokkosBlas::syr[noETI]");

    typedef typename AViewType::size_type size_type;
    const size_type numRows = A.extent(0);
    const size_type numCols = A.extent(1);

    bool justTranspose = (trans[0] == 'T') || (trans[0] == 't');
    bool justUp        = (uplo[0] == 'U') || (uplo[0] == 'u');

    // Prefer int as the index type, but use a larsyr type if needed.
    if ((numRows < static_cast<size_type>(INT_MAX)) && (numCols < static_cast<size_type>(INT_MAX))) {
      if (justTranspose) {
        if (justUp) {
          generalSyrImpl<ExecutionSpace, XViewType, AViewType, int, true, true>(space, alpha, x, A);
        } else {
          generalSyrImpl<ExecutionSpace, XViewType, AViewType, int, true, false>(space, alpha, x, A);
        }
      } else {
        if (justUp) {
          generalSyrImpl<ExecutionSpace, XViewType, AViewType, int, false, true>(space, alpha, x, A);
        } else {
          generalSyrImpl<ExecutionSpace, XViewType, AViewType, int, false, false>(space, alpha, x, A);
        }
      }
    } else {
      if (justTranspose) {
        if (justUp) {
          generalSyrImpl<ExecutionSpace, XViewType, AViewType, int64_t, true, true>(space, alpha, x, A);
        } else {
          generalSyrImpl<ExecutionSpace, XViewType, AViewType, int64_t, true, false>(space, alpha, x, A);
        }
      } else {
        if (justUp) {
          generalSyrImpl<ExecutionSpace, XViewType, AViewType, int64_t, false, true>(space, alpha, x, A);
        } else {
          generalSyrImpl<ExecutionSpace, XViewType, AViewType, int64_t, false, false>(space, alpha, x, A);
        }
      }
    }

    Kokkos::Profiling::popRegion();
  }
#else
      ;
#endif  // if !defined(KOKKOSKERNELS_ETI_ONLY) ||
        // KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
};

}  // namespace Impl
}  // namespace KokkosBlas

//
// Macro for declaration of full specialization of KokkosBlas::Impl::SYR.
// This is NOT for users!!!
// All the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or more .cpp
// files.
//
#define KOKKOSBLAS2_SYR_ETI_SPEC_DECL(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                           \
  extern template struct SYR<                                                                                          \
      EXEC_SPACE,                                                                                                      \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      false, true>;

#define KOKKOSBLAS2_SYR_ETI_SPEC_INST(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                                           \
  template struct SYR<                                                                                                 \
      EXEC_SPACE,                                                                                                      \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      false, true>;

#include <KokkosBlas2_syr_tpl_spec_decl.hpp>

#endif  // KOKKOSBLAS2_SYR_SPEC_HPP_
