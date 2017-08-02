/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/


#ifndef KOKKOS_SPINWAIT_HPP
#define KOKKOS_SPINWAIT_HPP

#include <Kokkos_Macros.hpp>

#include <cstdint>

namespace Kokkos {
namespace Impl {

#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )

void spinwait_while_equal( volatile int32_t & flag , const int32_t value );
void spinwait_until_equal( volatile int32_t & flag , const int32_t value );

void spinwait_while_equal( volatile int64_t & flag , const int64_t value );
void spinwait_until_equal( volatile int64_t & flag , const int64_t value );

void yield_while_equal( volatile int32_t & flag , const int32_t value );
void yield_until_equal( volatile int32_t & flag , const int32_t value );

void yield_while_equal( volatile int64_t & flag , const int64_t value );
void yield_until_equal( volatile int64_t & flag , const int64_t value );

#else

KOKKOS_INLINE_FUNCTION
void spinwait_while_equal( volatile int32_t & , const int32_t ) {}
KOKKOS_INLINE_FUNCTION
void spinwait_until_equal( volatile int32_t & , const int32_t ) {}

KOKKOS_INLINE_FUNCTION
void spinwait_while_equal( volatile int64_t & , const int64_t ) {}
KOKKOS_INLINE_FUNCTION
void spinwait_until_equal( volatile int64_t & , const int64_t ) {}

KOKKOS_INLINE_FUNCTION
void yield_while_equal( volatile int32_t & , const int32_t ) {}
KOKKOS_INLINE_FUNCTION
void yield_until_equal( volatile int32_t & , const int32_t ) {}

KOKKOS_INLINE_FUNCTION
void yield_while_equal( volatile int64_t & , const int64_t ) {}
KOKKOS_INLINE_FUNCTION
void yield_until_equal( volatile int64_t & , const int64_t ) {}

#endif

} /* namespace Impl */
} /* namespace Kokkos */

#endif /* #ifndef KOKKOS_SPINWAIT_HPP */

