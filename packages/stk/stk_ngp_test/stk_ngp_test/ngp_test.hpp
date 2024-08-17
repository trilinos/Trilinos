#ifndef STK_NGP_TEST_NGP_TEST_HPP
#define STK_NGP_TEST_NGP_TEST_HPP

#include <gtest/gtest.h>
#include "GlobalReporter.hpp"
#include <Kokkos_Core.hpp>
#if KOKKOS_VERSION >= 40200
#include <Kokkos_Printf.hpp>
#endif
#include <stk_ngp_test/NgpTestDeviceMacros.hpp>

namespace ngp_testing {

struct NgpTestEnvironment {
  NgpTestEnvironment(int* argc, char** argv);
  ~NgpTestEnvironment();
  int run_all_tests();
  void finalize();
};

class Test : public ::testing::Test {
 protected:
  virtual void NGPSetUp() {}
  virtual void NGPTearDown() {}

 private:
  void SetUp() override final;
  void TearDown() override final;
  virtual void set_failure_in_teardown() const { ADD_FAILURE(); }
};

int get_max_failure_reports_per_test();
void set_max_failure_reports_per_test(const int n);

namespace internal {

NGP_TEST_FUNCTION
void add_failure(const char* condition, const char* location);
int report_failures();
void clear_failures();

template<typename T>
NGP_TEST_INLINE
bool expect_near(const T a, const T b, const T tolerance) {
  T diff = a - b;
  if(diff < 0) diff *= T(-1);
  return diff <= tolerance;
}

}
}


#define NGP_TEST(testCaseName, testName)                                \
  GTEST_TEST_(testCaseName, testName, ::ngp_testing::Test, ::testing::internal::GetTestTypeId())

#define NGP_TEST_F(testFixture, testName)                               \
  static_assert(std::is_base_of<::ngp_testing::Test, testFixture>::value, \
                "Test fixture must inherit from ::ngp_testing::Test");  \
  GTEST_TEST_(testFixture, testName, testFixture, ::testing::internal::GetTypeId<testFixture>())

#define NGP_TEST_STRINGIZE(x) #x
#define NUM_TO_STR(x) NGP_TEST_STRINGIZE(x)
#define LOCATION __FILE__ ":" NUM_TO_STR(__LINE__)

#ifdef __HIP_DEVICE_COMPILE__
//FIXME: unsupported indirect call to function on HIP-Clang
#define NGP_EXPECT_TRUE(cond)
#define NGP_ASSERT_TRUE(cond)

#else
#define NGP_EXPECT_TRUE(cond)                                   \
  do {                                                          \
    if (!(cond)) {                                              \
      ::ngp_testing::internal::add_failure(#cond, LOCATION);    \
    }                                                           \
  } while (false)

#define NGP_ASSERT_TRUE(cond)                                   \
  do {                                                          \
    if (!(cond)) {                                              \
      ::ngp_testing::internal::add_failure(#cond, LOCATION);    \
      return;                                                   \
    }                                                           \
  } while (false)
#endif

#define NGP_EXPECT_FALSE(cond) NGP_EXPECT_TRUE(!(cond))
#define NGP_ASSERT_FALSE(cond) NGP_ASSERT_TRUE(!(cond))

#define NGP_EXPECT_EQ(a, b) NGP_EXPECT_TRUE((a) == (b))
#define NGP_ASSERT_EQ(a, b) NGP_ASSERT_TRUE((a) == (b))

#define NGP_EXPECT_NE(a, b) NGP_EXPECT_TRUE((a) != (b))
#define NGP_ASSERT_NE(a, b) NGP_ASSERT_TRUE((a) != (b))

#define NGP_EXPECT_LT(a, b) NGP_EXPECT_TRUE((a) < (b))
#define NGP_ASSERT_LT(a, b) NGP_ASSERT_TRUE((a) < (b))

#define NGP_EXPECT_LE(a, b) NGP_EXPECT_TRUE((a) <= (b))
#define NGP_ASSERT_LE(a, b) NGP_ASSERT_TRUE((a) <= (b))

#define NGP_EXPECT_GT(a, b) NGP_EXPECT_TRUE((a) > (b))
#define NGP_ASSERT_GT(a, b) NGP_ASSERT_TRUE((a) > (b))

#define NGP_EXPECT_GE(a, b) NGP_EXPECT_TRUE((a) >= (b))
#define NGP_ASSERT_GE(a, b) NGP_ASSERT_TRUE((a) >= (b))

#ifdef __HIP_DEVICE_COMPILE__
//FIXME: unsupported indirect call to function on HIP-Clang
#define NGP_EXPECT_NEAR(a, b, tolerance)
#define NGP_ASSERT_NEAR(a, b, tolerance)

#else
#define NGP_EXPECT_NEAR(a, b, tolerance)                                \
  do {                                                                  \
    if (!::ngp_testing::internal::expect_near(a, b, tolerance)) {       \
      ::ngp_testing::internal::add_failure("|(" #a ") - (" #b ")| <= " #tolerance, LOCATION); \
    }                                                                   \
  } while (false)

#define NGP_ASSERT_NEAR(a, b, tolerance)                                \
  do {                                                                  \
    if (!::ngp_testing::internal::expect_near(a, b, tolerance)) {       \
      ::ngp_testing::internal::add_failure("|(" #a ") - (" #b ")| <= " #tolerance, LOCATION); \
      return;                                                           \
    }                                                                   \
  } while (false)
#endif

namespace ngp_testing {

inline
NgpTestEnvironment::NgpTestEnvironment(int* argc, char** argv) {
  Kokkos::initialize(*argc, argv);
  initialize_reporters();
  ::testing::InitGoogleTest(argc, argv);
  Kokkos::push_finalize_hook(finalize_reporters);
}

inline
NgpTestEnvironment::~NgpTestEnvironment() {
  finalize();
}

inline
void NgpTestEnvironment::finalize() {
  if(Kokkos::is_initialized()) Kokkos::finalize();
}

inline
int NgpTestEnvironment::run_all_tests() {
  return RUN_ALL_TESTS();
}

inline
void Test::SetUp() {
  internal::clear_failures();
  NGPSetUp();
}

inline
void Test::TearDown() {
  int numReports = internal::report_failures();
  if(numReports > 0) {
    set_failure_in_teardown();
  }
  NGPTearDown();
}

inline
int get_max_failure_reports_per_test() {
  return get_host_reporter()->get_capacity();
}

inline
void set_max_failure_reports_per_test(const int n) {
  get_host_reporter()->resize(n);
  get_device_reporter_on_host()->resize(n);
}


namespace internal {

NGP_TEST_INLINE
void add_failure(const char* condition, const char* location) {
  KOKKOS_IF_ON_HOST((
    get_host_reporter()->add_failure(condition, location);
  ))  

  KOKKOS_IF_ON_DEVICE((
    get_device_reporter()->add_failure(condition, location);
  ))  
}

inline
int report_failures() {
  int numFailures = get_host_reporter()->report_failures();
  numFailures += get_device_reporter_on_host()->report_failures();
  return numFailures;
}

inline
void clear_failures() {
  get_host_reporter()->clear();
  get_device_reporter_on_host()->clear();
}

}

}

#endif

