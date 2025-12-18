// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, batched_scalar_serial_axpy_nt_float_float) { test_batched_axpy<TestDevice, float, float>(); }
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, batched_scalar_serial_axpy_nt_double_double) { test_batched_axpy<TestDevice, double, double>(); }
#endif
