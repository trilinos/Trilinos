/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#ifndef KOKKOSBLAS_TRTRI_SPEC_HPP_
#define KOKKOSBLAS_TRTRI_SPEC_HPP_

#include "KokkosKernels_config.h"
#include "Kokkos_Core.hpp"

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include <KokkosBlas_trtri_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class RVIT, class AVIT>
struct trtri_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

//
// This Macros provides the ETI specialization of trtri, currently not
// available.
//
#define KOKKOSBLAS_TRTRI_ETI_SPEC_AVAIL(SCALAR, LAYOUTA, EXEC_SPACE,         \
                                        MEM_SPACE)                           \
  template <>                                                                \
  struct trtri_eti_spec_avail<                                               \
      Kokkos::View<int, LAYOUTA, Kokkos::HostSpace,                          \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                \
      Kokkos::View<SCALAR**, LAYOUTA, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {             \
    enum : bool { value = true };                                            \
  };

// Include the actual specialization declarations
#include <KokkosBlas_trtri_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosBlas_trtri_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

//
// trtri
//

// Unification layer
template <class RVIT, class AVIT,
          bool tpl_spec_avail = trtri_tpl_spec_avail<RVIT, AVIT>::value,
          bool eti_spec_avail = trtri_eti_spec_avail<RVIT, AVIT>::value>
struct TRTRI {
  static void trtri(const RVIT& R, const char uplo[], const char diag[],
                    const AVIT& A);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
template <class RVIT, class AVIT>
struct TRTRI<RVIT, AVIT, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void trtri(const RVIT& R, const char uplo[], const char diag[],
                    const AVIT& A) {
    static_assert(Kokkos::is_view<AVIT>::value, "AVIT must be a Kokkos::View.");
    static_assert(static_cast<int>(AVIT::rank) == 2, "AVIT must have rank 2.");

    Kokkos::Profiling::pushRegion(KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
                                      ? "KokkosBlas::trtri[ETI]"
                                      : "KokkosBlas::trtri[noETI]");

    typename AVIT::HostMirror host_A = Kokkos::create_mirror_view(A);
    typename RVIT::HostMirror host_R = Kokkos::create_mirror_view(R);

    Kokkos::deep_copy(host_A, A);

    SerialTrtri_Invoke<typename RVIT::HostMirror, typename AVIT::HostMirror>(
        R, uplo, diag, host_A);

    Kokkos::deep_copy(A, host_A);

    Kokkos::Profiling::popRegion();
  }
};
#endif  //! defined(KOKKOSKERNELS_ETI_ONLY) ||
        //! KOKKOSKERNELS_IMPL_COMPILE_LIBRARY

}  // namespace Impl
}  // namespace KokkosBlas

//
// These Macros are only included when we are not compiling libkokkoskernels but
// are auto generating files. These macros provide the explicit instantiation
// declaration and definition of TRTRI, potentially reducing user code size. The
// "extern template" skips the implicit instatiation step ensuring that the
// callers code uses this explicit instantiation definition of TRTRI.
//
#define KOKKOSBLAS_TRTRI_ETI_SPEC_DECL(SCALAR, LAYOUTA, EXEC_SPACE, MEM_SPACE) \
  extern template struct TRTRI<                                                \
      Kokkos::View<int, LAYOUTA, Kokkos::HostSpace,                            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<SCALAR**, LAYOUTA, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      false, true>;

#define KOKKOSBLAS_TRTRI_ETI_SPEC_INST(SCALAR, LAYOUTA, EXEC_SPACE, MEM_SPACE) \
  template struct TRTRI<                                                       \
      Kokkos::View<int, LAYOUTA, Kokkos::HostSpace,                            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<SCALAR**, LAYOUTA, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      false, true>;

#include <KokkosBlas_trtri_tpl_spec_decl.hpp>
#include <generated_specializations_hpp/KokkosBlas_trtri_eti_spec_decl.hpp>

#endif  // KOKKOSBLAS_TRTRI_SPEC_HPP_
