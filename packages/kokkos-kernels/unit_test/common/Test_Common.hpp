#ifndef TEST_COMMON_HPP
#define TEST_COMMON_HPP

// FIXME_SYCL still some uses of the wrong namespace
#ifndef KOKKOS_ENABLE_SYCL
#include <Test_Common_ArithTraits.hpp>
#endif
// #include<Test_Common_float128.hpp>
#include <Test_Common_set_bit_count.hpp>
#include <Test_Common_Sorting.hpp>
#include <Test_Common_Transpose.hpp>
#include <Test_Common_IOUtils.hpp>
#include <Test_Common_Error.hpp>

#endif  // TEST_COMMON_HPP
