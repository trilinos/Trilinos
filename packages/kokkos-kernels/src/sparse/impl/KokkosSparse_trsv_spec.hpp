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
#ifndef KOKKOSSPARSE_IMPL_TRSV_SPEC_HPP_
#define KOKKOSSPARSE_IMPL_TRSV_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include "KokkosSparse_CrsMatrix.hpp"

// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY 
#include <KokkosSparse_trsv_impl.hpp>
#endif

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
struct trsv_eti_spec_avail {
  enum : bool { value = false };
};

}
}


#define KOKKOSSPARSE_TRSV_ETI_SPEC_AVAIL( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE) \
    template<> \
    struct trsv_eti_spec_avail< \
                  KokkosSparse::CrsMatrix< const SCALAR_TYPE, const ORDINAL_TYPE, \
                                           Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                           Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
                                           const OFFSET_TYPE>, \
                  Kokkos::View<const SCALAR_TYPE **, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                  Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE,  \
                               Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged> > > \
    { enum : bool { value = true }; };

// Include the actual specialization declarations
#include<KokkosSparse_trsv_tpl_spec_avail.hpp>
#include<generated_specializations_hpp/KokkosSparse_trsv_eti_spec_avail.hpp>

namespace KokkosSparse {
namespace Impl {

// Unification layer
/// \brief Implementation of KokkosSparse::trsv (sparse matrix - dense
///   vector multiply) for single vectors (1-D Views).

template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType,
         bool tpl_spec_avail =
             trsv_tpl_spec_avail< CrsMatrixType,
                                   DomainMultiVectorType,
                                   RangeMultiVectorType>::value,
         bool eti_spec_avail =
             trsv_eti_spec_avail<   CrsMatrixType,
                                   DomainMultiVectorType,
                                   RangeMultiVectorType>::value >
struct TRSV{
  static void
  trsv (const char uplo[],
        const char trans[],
        const char diag[],
        const CrsMatrixType& A,
        DomainMultiVectorType B,
        RangeMultiVectorType X);
};


#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of trsv for multi vectors.
// Unification layer
template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
struct TRSV< CrsMatrixType, DomainMultiVectorType, RangeMultiVectorType, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY>{
  static void
  trsv (const char uplo[],
        const char trans[],
        const char diag[],
        const CrsMatrixType& A,
        DomainMultiVectorType B,
        RangeMultiVectorType X) // X is the output MV
  {
    if (trans[0] == 'N' || trans[0] == 'n') {       // no transpose
      if (uplo[0] == 'L' || uplo[0] == 'l') {   // lower triangular
        if (diag[0] == 'U' || diag[0] == 'u') {    // unit diagonal
          Sequential::lowerTriSolveCsrUnitDiag (X, A, B);
        } else {                               // non unit diagonal
          Sequential::lowerTriSolveCsr (X, A, B);
        }
      } else {                                  // upper triangular
        if (diag[0] == 'U' || diag[0] == 'u') {    // unit diagonal
          Sequential::upperTriSolveCsrUnitDiag (X, A, B);
        } else {                               // non unit diagonal
          Sequential::upperTriSolveCsr (X, A, B);
        }
      }
    }
    else if (trans[0] == 'T' || trans[0] == 't') {     // transpose
      if (uplo[0] == 'L' || uplo[0] == 'l') {   // lower triangular
        // Transposed lower tri CSR => upper tri CSC.
        if (diag[0] == 'U' || diag[0] == 'u') {    // unit diagonal
          Sequential::upperTriSolveCscUnitDiag (X, A, B);
        } else {                               // non unit diagonal
          Sequential::upperTriSolveCsc (X, A, B);
        }
      }
      else {                                    // upper triangular
        // Transposed upper tri CSR => lower tri CSC.
        if (diag[0] == 'U' || diag[0] == 'u') {    // unit diagonal
          Sequential::lowerTriSolveCscUnitDiag (X, A, B);
        } else {                               // non unit diagonal
          Sequential::lowerTriSolveCsc (X, A, B);
        }
      }
    }
    else if (trans[0] == 'C' || trans[0] == 'c') { // conj transpose
      if (uplo[0] == 'L' || uplo[0] == 'l') {    // lower triangular
        // Transposed lower tri CSR => upper tri CSC.
        if (diag[0] == 'U' || diag[0] == 'u') {     // unit diagonal
          Sequential::upperTriSolveCscUnitDiagConj (X, A, B);
        } else {                                // non unit diagonal
          Sequential::upperTriSolveCscConj (X, A, B);
        }
      }
      else {                                     // upper triangular
        // Transposed upper tri CSR => lower tri CSC.
        if (diag[0] == 'U' || diag[0] == 'u') {     // unit diagonal
          Sequential::lowerTriSolveCscUnitDiagConj (X, A, B);
        } else {                                // non unit diagonal
          Sequential::lowerTriSolveCscConj (X, A, B);
        }
      }
    }
  }
};

#endif
}
}

//
// Macro for declaration of full specialization of
// KokkosSparse::Impl::Dot for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSSPARSE_TRSV_ETI_SPEC_DECL( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE ) \
    extern template struct  \
    TRSV<             KokkosSparse::CrsMatrix< const SCALAR_TYPE, const ORDINAL_TYPE, \
                                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
                                   const OFFSET_TYPE>, \
                      Kokkos::View<const SCALAR_TYPE **, LAYOUT_TYPE,  \
                                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE,  \
                                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                      false, true >; \

#define KOKKOSSPARSE_TRSV_ETI_SPEC_INST( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE) \
    template struct  \
    TRSV< KokkosSparse::CrsMatrix< const SCALAR_TYPE, const ORDINAL_TYPE, \
                                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
                                   const OFFSET_TYPE>, \
                      Kokkos::View<const SCALAR_TYPE **, LAYOUT_TYPE,  \
                                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> >, \
                      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE,  \
                                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                      false, true > ;

#include<KokkosSparse_trsv_tpl_spec_decl.hpp>
#include<generated_specializations_hpp/KokkosSparse_trsv_eti_spec_decl.hpp>


#endif // KOKKOS_BLAS1_MV_IMPL_DOT_HPP_
