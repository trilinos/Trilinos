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

/// \file Test_Common_Test_All_Type_Combos.hpp

/**
 * KOKKOSKERNELS_EXECUTE_TEST should take (SCALAR, ORDINAL, OFFSET, DEVICE). All
 * these args are types.
 * #define NO_TEST_COMPLEX to skip testing of kokkos complex types
 */

#if !defined(KOKKOSKERNELS_EXECUTE_TEST)
#error Test_Common_Test_All_Type_Combos.hpp requires KOKKOSKERNELS_EXECUTE_TEST to be set
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && \
     !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))

// ETI is off, test all possible type combos

KOKKOSKERNELS_EXECUTE_TEST(double, int, int, TestExecSpace)
KOKKOSKERNELS_EXECUTE_TEST(double, int64_t, int, TestExecSpace)
KOKKOSKERNELS_EXECUTE_TEST(double, int, size_t, TestExecSpace)
KOKKOSKERNELS_EXECUTE_TEST(double, int64_t, size_t, TestExecSpace)
KOKKOSKERNELS_EXECUTE_TEST(float, int, int, TestExecSpace)
KOKKOSKERNELS_EXECUTE_TEST(float, int64_t, int, TestExecSpace)
KOKKOSKERNELS_EXECUTE_TEST(float, int, size_t, TestExecSpace)
KOKKOSKERNELS_EXECUTE_TEST(float, int64_t, size_t, TestExecSpace)

#if !defined(NO_TEST_COMPLEX)

KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_double, int, int, TestExecSpace)
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_double, int64_t, int, TestExecSpace)
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_double, int, size_t, TestExecSpace)
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_double, int64_t, size_t,
                           TestExecSpace)
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_float, int, int, TestExecSpace)
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_float, int64_t, int, TestExecSpace)
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_float, int, size_t, TestExecSpace)
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_float, int64_t, size_t, TestExecSpace)

#endif

#else

// ETI is on, only test instantiated type combos

#if (defined(KOKKOSKERNELS_INST_DOUBLE) &&      \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT))
KOKKOSKERNELS_EXECUTE_TEST(double, int, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_DOUBLE) &&          \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT))
KOKKOSKERNELS_EXECUTE_TEST(double, int64_t, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_DOUBLE) &&      \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) && \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T))
KOKKOSKERNELS_EXECUTE_TEST(double, int, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_DOUBLE) &&          \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T))
KOKKOSKERNELS_EXECUTE_TEST(double, int64_t, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) &&       \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT))
KOKKOSKERNELS_EXECUTE_TEST(float, int, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) &&           \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT))
KOKKOSKERNELS_EXECUTE_TEST(float, int64_t, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) &&       \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) && \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T))
KOKKOSKERNELS_EXECUTE_TEST(float, int, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) &&           \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T))
KOKKOSKERNELS_EXECUTE_TEST(float, int64_t, size_t, TestExecSpace)
#endif

#if !defined(NO_TEST_COMPLEX)

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&            \
     defined(KOKKOSKERNELS_INST_OFFSET_INT))
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_double, int, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) &&        \
     defined(KOKKOSKERNELS_INST_OFFSET_INT))
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_double, int64_t, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&            \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T))
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_double, int, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) &&        \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T))
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_double, int64_t, size_t,
                           TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&           \
     defined(KOKKOSKERNELS_INST_OFFSET_INT))
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_float, int, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) &&       \
     defined(KOKKOSKERNELS_INST_OFFSET_INT))
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_float, int64_t, int, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT) &&           \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T))
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_float, int, size_t, TestExecSpace)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && \
     defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) &&       \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T))
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_float, int64_t, size_t, TestExecSpace)
#endif

#endif  // !NO_TEST_COMPLEX

#endif  // ETI ON
