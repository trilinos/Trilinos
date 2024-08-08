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

/// \file Test_Common_Test_All_Type_Combos.hpp

/**
 * KOKKOSKERNELS_EXECUTE_TEST should take (SCALAR, ORDINAL, OFFSET, DEVICE). All
 * these args are types.
 * #define NO_TEST_COMPLEX to skip testing of kokkos complex types
 */

#if !defined(KOKKOSKERNELS_EXECUTE_TEST)
#error Test_Common_Test_All_Type_Combos.hpp requires KOKKOSKERNELS_EXECUTE_TEST to be set
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))

// ETI is off, test all possible type combos

KOKKOSKERNELS_EXECUTE_TEST(double, int, int, TestDevice)
KOKKOSKERNELS_EXECUTE_TEST(double, int64_t, int, TestDevice)
KOKKOSKERNELS_EXECUTE_TEST(double, int, size_t, TestDevice)
KOKKOSKERNELS_EXECUTE_TEST(double, int64_t, size_t, TestDevice)
KOKKOSKERNELS_EXECUTE_TEST(float, int, int, TestDevice)
KOKKOSKERNELS_EXECUTE_TEST(float, int64_t, int, TestDevice)
KOKKOSKERNELS_EXECUTE_TEST(float, int, size_t, TestDevice)
KOKKOSKERNELS_EXECUTE_TEST(float, int64_t, size_t, TestDevice)

#if !defined(NO_TEST_COMPLEX)

KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_double, int, int, TestDevice)
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_double, int64_t, int, TestDevice)
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_double, int, size_t, TestDevice)
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_double, int64_t, size_t, TestDevice)
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_float, int, int, TestDevice)
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_float, int64_t, int, TestDevice)
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_float, int, size_t, TestDevice)
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_float, int64_t, size_t, TestDevice)

#endif

#else

// ETI is on, only test instantiated type combos

#if (defined(KOKKOSKERNELS_INST_DOUBLE) && defined(KOKKOSKERNELS_INST_ORDINAL_INT) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT))
KOKKOSKERNELS_EXECUTE_TEST(double, int, int, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_DOUBLE) && defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT))
KOKKOSKERNELS_EXECUTE_TEST(double, int64_t, int, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_DOUBLE) && defined(KOKKOSKERNELS_INST_ORDINAL_INT) && \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T))
KOKKOSKERNELS_EXECUTE_TEST(double, int, size_t, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_DOUBLE) && defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T))
KOKKOSKERNELS_EXECUTE_TEST(double, int64_t, size_t, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) && defined(KOKKOSKERNELS_INST_ORDINAL_INT) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT))
KOKKOSKERNELS_EXECUTE_TEST(float, int, int, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) && defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT))
KOKKOSKERNELS_EXECUTE_TEST(float, int64_t, int, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) && defined(KOKKOSKERNELS_INST_ORDINAL_INT) && \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T))
KOKKOSKERNELS_EXECUTE_TEST(float, int, size_t, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT) && defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T))
KOKKOSKERNELS_EXECUTE_TEST(float, int64_t, size_t, TestDevice)
#endif

#if !defined(NO_TEST_COMPLEX)

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && defined(KOKKOSKERNELS_INST_ORDINAL_INT) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT))
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_double, int, int, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT))
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_double, int64_t, int, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && defined(KOKKOSKERNELS_INST_ORDINAL_INT) && \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T))
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_double, int, size_t, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) && defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T))
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_double, int64_t, size_t, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && defined(KOKKOSKERNELS_INST_ORDINAL_INT) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT))
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_float, int, int, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_INT))
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_float, int64_t, int, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && defined(KOKKOSKERNELS_INST_ORDINAL_INT) && \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T))
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_float, int, size_t, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) && defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && \
     defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T))
KOKKOSKERNELS_EXECUTE_TEST(kokkos_complex_float, int64_t, size_t, TestDevice)
#endif

#endif  // !NO_TEST_COMPLEX

#endif  // ETI ON
