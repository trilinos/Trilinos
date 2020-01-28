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
#ifndef KOKKOSSPARSE_IMPL_SPTRSV_SYMBOLIC_SPEC_HPP_
#define KOKKOSSPARSE_IMPL_SPTRSV_SYMBOLIC_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosKernels_Handle.hpp"

// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY 
#include <KokkosSparse_sptrsv_symbolic_impl.hpp>
#endif

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template<class KernelHandle,
         class RowMapType,
         class EntriesType>
struct sptrsv_symbolic_eti_spec_avail {
  enum : bool { value = false };
};

}
}


#define KOKKOSSPARSE_SPTRSV_SYMBOLIC_ETI_SPEC_AVAIL( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE) \
    template<> \
    struct sptrsv_symbolic_eti_spec_avail< \
                  KokkosKernels::Experimental::KokkosKernelsHandle<\
                               const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE,  \
                               EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE> , \
                  Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > > \
    { enum : bool { value = true }; };

// Include the actual specialization declarations
#include<KokkosSparse_sptrsv_symbolic_tpl_spec_avail.hpp>
#include<generated_specializations_hpp/KokkosSparse_sptrsv_symbolic_eti_spec_avail.hpp>

namespace KokkosSparse {
namespace Impl {

// Unification layer
/// \brief Implementation of KokkosSparse::sptrsv_symbolic

template<class KernelHandle,
         class RowMapType,
         class EntriesType,
         bool tpl_spec_avail =
             sptrsv_symbolic_tpl_spec_avail< KernelHandle,
                                  RowMapType,
                                  EntriesType>::value,
         bool eti_spec_avail =
             sptrsv_symbolic_eti_spec_avail< KernelHandle,
                                  RowMapType,
                                  EntriesType>::value >
struct SPTRSV_SYMBOLIC{
  static void
  sptrsv_symbolic (KernelHandle *handle,
                   const RowMapType row_map,
                   const EntriesType entries);
};


#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of sptrsv_symbolic
// Unification layer
template<class KernelHandle,
         class RowMapType,
         class EntriesType>
struct SPTRSV_SYMBOLIC<KernelHandle, RowMapType, EntriesType, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY>{
  static void
  sptrsv_symbolic (KernelHandle *handle,
                   const RowMapType row_map,
                   const EntriesType entries)
  {
    auto sptrsv_handle = handle->get_sptrsv_handle();
    auto nrows = row_map.extent(0)-1;
    sptrsv_handle->new_init_handle(nrows);

    if ( sptrsv_handle->is_lower_tri() ) {
      Experimental::lower_tri_symbolic(*sptrsv_handle, row_map, entries);
      sptrsv_handle->set_symbolic_complete();
    }
    else {
      Experimental::upper_tri_symbolic(*sptrsv_handle, row_map, entries);
      sptrsv_handle->set_symbolic_complete();
    }
  }

};
#endif
} // namespace Impl
} // namespace KokkosSparse

//
// Macro for declaration of full specialization of
// This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSSPARSE_SPTRSV_SYMBOLIC_ETI_SPEC_DECL( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE ) \
    extern template struct  \
    SPTRSV_SYMBOLIC<\
                  KokkosKernels::Experimental::KokkosKernelsHandle<\
                               const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE,  \
                               EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE> , \
                  Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                      false, true >; \

#define KOKKOSSPARSE_SPTRSV_SYMBOLIC_ETI_SPEC_INST( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE) \
    template struct  \
    SPTRSV_SYMBOLIC<\
                  KokkosKernels::Experimental::KokkosKernelsHandle<\
                               const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE,  \
                               EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE> , \
                  Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                      false, true > ;

#include<KokkosSparse_sptrsv_symbolic_tpl_spec_decl.hpp>
#include<generated_specializations_hpp/KokkosSparse_sptrsv_symbolic_eti_spec_decl.hpp>


#endif
