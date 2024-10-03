// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef stk_util_config_h
#define stk_util_config_h

#ifdef STK_BUILT_WITH_BJAM

#define STK_HAS_MPI
#define STK_HAS_ARBORX
#define STK_HAVE_BOOST
#define STK_HAVE_KOKKOS
#define STK_HAVE_INTREPID2
#define STK_HAVE_STKMESH
#define STK_16BIT_CONNECTIVITY_ORDINAL
#define STK_HAVE_STKIO
#define STK_HAVE_STKNGP_TEST
#define STK_HAS_SEACAS_IOSS
#define STK_HAS_SEACAS_EXODUS
#define STK_HAS_SEACAS_NEMESIS

#else
// This file gets created by cmake during a Trilinos build
// and will not be present in a sierra build using bjam or associated wrappers
#include "STK_Trilinos_config.h"

#ifndef STK_HAS_MPI

#ifndef MPI_Comm
#define MPI_Comm int
#endif

#ifndef MPI_COMM_NULL
#define MPI_COMM_NULL 0
#endif

#ifndef MPI_COMM_SELF
#define MPI_COMM_SELF 0
#endif

#endif // STK_HAS_MPI
#endif // STK_BUILT_WITH_BJAM

// GCC address sanitizer
#ifdef __SANITIZE_ADDRESS__
#  define STK_ASAN_IS_ON
#endif

// Clang address sanitizer
#if !defined(STK_ASAN_IS_ON) && defined(__has_feature)
#  if __has_feature(address_sanitizer)
#    define STK_ASAN_IS_ON
#  endif
#endif

//----------------------------------------------------------------------

// Use macro below to deprecate:
//   classes (class STK_DEPRECATED Class;), 
//   structs (struct STK_DEPRECATED Struct;), 
//   typedefs (STK_DEPRECATED typedef Type 1 Type2;, using Type1 STK_DEPRECATED = Type2;), 
//   variables (STK_DEPRECATED int variable;), 
//   non-static data members (union Union { STK_DEPRECATED int variable; }), 
//   functions (STK_DEPRECATED void function();),
//   inline functions (STK_DEPRECATED inline void function();),
//   namespaces (namespace STK_DEPRECATED stk { int variable; }),
//   enumeration (enum STK_DEPRECATED Enum{};),
//   enumerators (enum {Type 1 STK_DEPRECATED, Type2 DEPRECATED};), and
//   template specialization (template<> struct STK_DEPRECATED Struct<int>;).
//
// This is basically copied from the Trilinos version in Tribits to maintain some compatibility
/* Usage Example
 * #ifndef STK_HIDE_DEPRECATED_CODE // Delete after FILL_IN_DATE_TWO_SPRINTS_AFTER_END_OF_THIS_SPRINT_HERE
 * STK_DEPRECATED bool modification_end(impl::MeshModification::modification_optimization opt)
 * {
 *     if (impl::MeshModification::MOD_END_SORT == opt) {
 *         return m_meshModification.modification_end();
 *     } else {
 *         return m_meshModification.modification_end_with_compress();
 *     }
 * }
 * #endif // STK_HIDE_DEPRECATED_CODE
 */
#ifdef STK_SHOW_DEPRECATED_WARNINGS
#  ifndef STK_DEPRECATED
#    define STK_DEPRECATED [[deprecated]]
#  endif
#  ifndef STK_DEPRECATED_MSG
#    define STK_DEPRECATED_MSG(MSG) [[deprecated(#MSG)]]
#  endif
#else
#  ifndef STK_DEPRECATED
#    define STK_DEPRECATED
#    define STK_DEPRECATED_MSG(MSG)
#  endif
#endif

#include "stk_util/stk_kokkos_macros.h"

#endif /* stk_util_config_h */
