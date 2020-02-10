/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#ifndef KOKKOSBLAS_IMPL_GESV_SPEC_HPP_
#define KOKKOSBLAS_IMPL_GESV_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>

// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY 
#include <KokkosBlas_gesv_impl.hpp>
#endif

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template<class AVT, class BVT>
struct gesv_eti_spec_avail {
  enum : bool { value = false };
};
}
}

//
// Macro for declaration of full specialization availability
// KokkosBlas::Impl::GESV.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _INST macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS_GESV_ETI_SPEC_AVAIL( SCALAR_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE) \
    template<> \
    struct gesv_eti_spec_avail< \
                  Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                  Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged> > > \
    { enum : bool { value = true }; };

// Include the actual specialization declarations
#include<KokkosBlas_gesv_tpl_spec_avail.hpp>
#include<generated_specializations_hpp/KokkosBlas_gesv_eti_spec_avail.hpp>

namespace KokkosBlas {
namespace Impl {

// Unification layer
/// \brief Implementation of KokkosBlas::gesv.

template<class AMatrix,
         class BXMV,
         class IPIVV,
         bool tpl_spec_avail = gesv_tpl_spec_avail<AMatrix, BXMV>::value,
         bool eti_spec_avail = gesv_eti_spec_avail<AMatrix, BXMV>::value
        >
struct GESV{
  static void
  gesv (AMatrix& A,
        BXMV& B,
        IPIVV& IPIV);
};


#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of gesv for multi vectors.
// Unification layer
template<class AMatrix,
         class BXMV,
         class IPIVV>
struct GESV<AMatrix, BXMV, IPIVV, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY>{
  static void
  gesv (const AMatrix& A,
        const BXMV& B,
        const IPIVV& IPIV)
  {
   //NOTE: Might add the implementation of KokkosBlas::gesv later
  }
};

#endif
}// namespace Impl
}// namespace KokkosBlas

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::GESV.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSBLAS_GESV_ETI_SPEC_DECL( SCALAR_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE ) \
    extern template struct  \
    GESV<             Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE,  \
                                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE,  \
                                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                      Kokkos::View<int *, LAYOUT_TYPE,  \
                                   Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                      false, true >; \

#define KOKKOSBLAS_GESV_ETI_SPEC_INST( SCALAR_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE) \
    template struct  \
    GESV<             Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE,  \
                                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE,  \
                                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                      Kokkos::View<int *, LAYOUT_TYPE,  \
                                   Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                      false, true > ;

#include<KokkosBlas_gesv_tpl_spec_decl.hpp>
#include<generated_specializations_hpp/KokkosBlas_gesv_eti_spec_decl.hpp>


#endif // KOKKOSBLAS_IMPL_GESV_SPEC_HPP_
