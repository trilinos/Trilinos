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
#ifndef KOKKOSSPARSE_IMPL_SPILUK_NUMERIC_SPEC_HPP_
#define KOKKOSSPARSE_IMPL_SPILUK_NUMERIC_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosKernels_Handle.hpp"

// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY 
#include <KokkosSparse_spiluk_numeric_impl.hpp>
#include <KokkosSparse_spiluk_symbolic_impl.hpp>
#endif

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template<class KernelHandle,
         class ARowMapType,
         class AEntriesType,
         class AValuesType,
         class LRowMapType,
         class LEntriesType,
         class LValuesType,
         class URowMapType,
         class UEntriesType,
         class UValuesType>
struct spiluk_numeric_eti_spec_avail {
  enum : bool { value = false };
};

}
}


#define KOKKOSSPARSE_SPILUK_NUMERIC_ETI_SPEC_AVAIL( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE) \
    template<> \
    struct spiluk_numeric_eti_spec_avail< \
                  KokkosKernels::Experimental::KokkosKernelsHandle<\
                               const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE,  \
                               EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE> , \
                  Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > > \
    { enum : bool { value = true }; };

// Include the actual specialization declarations
#include<KokkosSparse_spiluk_numeric_tpl_spec_avail.hpp>
#include<generated_specializations_hpp/KokkosSparse_spiluk_numeric_eti_spec_avail.hpp>

namespace KokkosSparse {
namespace Impl {

// Unification layer
/// \brief Implementation of KokkosSparse::spiluk_numeric

template<class KernelHandle,
         class ARowMapType,
         class AEntriesType,
         class AValuesType,
         class LRowMapType,
         class LEntriesType,
         class LValuesType,
         class URowMapType,
         class UEntriesType,
         class UValuesType,
         bool tpl_spec_avail =
              spiluk_numeric_tpl_spec_avail<KernelHandle,
                                            ARowMapType,
                                            AEntriesType,
                                            AValuesType,
                                            LRowMapType,
                                            LEntriesType,
                                            LValuesType,
                                            URowMapType,
                                            UEntriesType,
                                            UValuesType>::value,
         bool eti_spec_avail =
              spiluk_numeric_eti_spec_avail<KernelHandle,
                                            ARowMapType,
                                            AEntriesType,
                                            AValuesType,
                                            LRowMapType,
                                            LEntriesType,
                                            LValuesType,
                                            URowMapType,
                                            UEntriesType,
                                            UValuesType>::value>
struct SPILUK_NUMERIC{
  static void
  spiluk_numeric (KernelHandle *handle,
                  const typename KernelHandle::const_nnz_lno_t &fill_lev,
                  const ARowMapType&  A_row_map,
                  const AEntriesType& A_entries,
                  const AValuesType&  A_values,
                        LRowMapType&  L_row_map,
                        LEntriesType& L_entries,
                        LValuesType&  L_values,
                        URowMapType&  U_row_map,
                        UEntriesType& U_entries,
                        UValuesType&  U_values);
};


#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of spiluk_numeric
// Unification layer
template<class KernelHandle,
         class ARowMapType,
         class AEntriesType,
         class AValuesType,
         class LRowMapType,
         class LEntriesType,
         class LValuesType,
         class URowMapType,
         class UEntriesType,
         class UValuesType>
struct SPILUK_NUMERIC<KernelHandle, ARowMapType, AEntriesType, AValuesType, LRowMapType, LEntriesType, LValuesType, URowMapType, UEntriesType, UValuesType, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY>{
  static void
  spiluk_numeric (KernelHandle *handle,
                  const typename KernelHandle::const_nnz_lno_t &fill_lev,
                  const ARowMapType&  A_row_map,
                  const AEntriesType& A_entries,
                  const AValuesType&  A_values,
                        LRowMapType&  L_row_map,
                        LEntriesType& L_entries,
                        LValuesType&  L_values,
                        URowMapType&  U_row_map,
                        UEntriesType& U_entries,
                        UValuesType&  U_values)
  {
    // Call specific algorithm type
    auto spiluk_handle = handle->get_spiluk_handle();

    if ( spiluk_handle->is_symbolic_complete() == false ) {
      Experimental::iluk_symbolic(*spiluk_handle, fill_lev, A_row_map, A_entries, L_row_map, L_entries, U_row_map, U_entries);
    }
    Experimental::iluk_numeric(*spiluk_handle, A_row_map, A_entries, A_values,
                                               L_row_map, L_entries, L_values, 
                                               U_row_map, U_entries, U_values);
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
#define KOKKOSSPARSE_SPILUK_NUMERIC_ETI_SPEC_DECL( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE ) \
    extern template struct  \
    SPILUK_NUMERIC<\
                  KokkosKernels::Experimental::KokkosKernelsHandle<\
                               const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE,  \
                               EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE> , \
                  Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, false, true >;

#define KOKKOSSPARSE_SPILUK_NUMERIC_ETI_SPEC_INST( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE) \
    template struct  \
    SPILUK_NUMERIC<\
                  KokkosKernels::Experimental::KokkosKernelsHandle<\
                               const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE,  \
                               EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE> , \
                  Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, false, true >;

#include<KokkosSparse_spiluk_numeric_tpl_spec_decl.hpp>
#include<generated_specializations_hpp/KokkosSparse_spiluk_numeric_eti_spec_decl.hpp>


#endif
