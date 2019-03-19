/*
 * Copyright (c) 2013, Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */

#ifndef stk_util_config_h
#define stk_util_config_h

#ifdef STK_BUILT_IN_SIERRA
#define STK_HAS_MPI
#else
// This file gets created by cmake during a Trilinos build
// and will not be present in a sierra build using bjam or associated wrappers
#include <stk_util/STK_Trilinos_config.h>
#ifdef HAVE_MPI
#define STK_HAS_MPI
#else

#ifndef MPI_Comm
#define MPI_Comm int
#endif

#ifndef MPI_COMM_NULL
#define MPI_COMM_NULL 0
#endif

#ifndef MPI_COMM_SELF
#define MPI_COMM_SELF 0
#endif

#endif
#endif

#define STK_PACKAGE stk
#define STK_HAS_SNL_EXODUSII

//----------------------------------------------------------------------

// Use macro below to deprecate functions (place at beginning of function or class method)
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
#    if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#      define STK_DEPRECATED  __attribute__((__deprecated__))
#    else
#      define STK_DEPRECATED
#    endif
#  endif
#  ifndef STK_DEPRECATED_MSG
#    if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 5))
#      define STK_DEPRECATED_MSG(MSG)  __attribute__((__deprecated__ (#MSG) ))
#    elif (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#      define STK_DEPRECATED_MSG(MSG)  __attribute__((__deprecated__))
#    else
#      define STK_DEPRECATED_MSG(MSG)
#    endif
#  endif
#else
#  ifndef STK_DEPRECATED
#    define STK_DEPRECATED
#    define STK_DEPRECATED_MSG(MSG)
#  endif
#endif

#include <stk_util/stk_kokkos_macros.h>

#endif /* stk_util_config_h */
