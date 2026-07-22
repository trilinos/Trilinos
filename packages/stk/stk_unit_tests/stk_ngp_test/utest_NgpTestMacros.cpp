#include <limits>
#include <stk_ngp_test/ngp_test.hpp>
#include <stk_util/ngp/NgpSpaces.hpp>

static constexpr double eps = std::numeric_limits<double>::epsilon();

class TestReports : public ::ngp_testing::Test {
 protected:
  int numExpectedReports;

  TestReports(int numReports) : numExpectedReports(numReports) {}

  void set_failure_in_teardown() const override {}
  void NGPTearDown() override {
    int numReports = ::ngp_testing::internal::report_failures();
    EXPECT_EQ(numExpectedReports, numReports);
  }
};

class TestReportHost : public TestReports {
 protected:
  TestReportHost() : TestReports(1) {}
};

class TestReportDevice : public TestReportHost {};

class TestReportHostAndDevice : public TestReports {
 protected:
  TestReportHostAndDevice() : TestReports(2) {}
};

class TestReportHostAndDeviceMany : public TestReports {
 protected:
  const int loopLen;
  TestReportHostAndDeviceMany() : TestReports(1024), loopLen(numExpectedReports / 2) {}
};

NGP_TEST_INLINE
void expect_true() {
  NGP_EXPECT_TRUE(0 == 1);
}

NGP_TEST_INLINE
void expect_false() {
  NGP_EXPECT_FALSE(1 == 1);
}

NGP_TEST_INLINE
void expect_eq() {
  NGP_EXPECT_EQ(0, 1);
}

NGP_TEST_INLINE
void expect_ne() {
  NGP_EXPECT_NE(1, 1);
}

NGP_TEST_INLINE
void expect_lt() {
  NGP_EXPECT_LT(1, 0);
}

NGP_TEST_INLINE
void expect_le() {
  NGP_EXPECT_LE(1, 0);
}

NGP_TEST_INLINE
void expect_gt() {
  NGP_EXPECT_GT(0, 1);
}

NGP_TEST_INLINE
void expect_ge() {
  NGP_EXPECT_GE(0, 1);
}

NGP_TEST_INLINE
void expect_near() {
  NGP_EXPECT_NEAR(-1.0 - 2*eps, -1.0, eps);
}

#define DEFINE_run_unary_assert(TEST) \
NGP_TEST_INLINE void run_assert_##TEST(int a, int b) {\
  NGP_ASSERT_##TEST(a == b); \
  NGP_EXPECT_##TEST(a == b); \
}

#define DEFINE_run_binary_assert(TEST) \
NGP_TEST_INLINE void run_assert_##TEST(int a, int b) {\
  NGP_ASSERT_##TEST(a, b); \
  NGP_EXPECT_##TEST(a, b); \
}

DEFINE_run_unary_assert(TRUE)
DEFINE_run_unary_assert(FALSE)
DEFINE_run_binary_assert(EQ)
DEFINE_run_binary_assert(NE)
DEFINE_run_binary_assert(LT)
DEFINE_run_binary_assert(LE)
DEFINE_run_binary_assert(GT)
DEFINE_run_binary_assert(GE)

NGP_TEST_INLINE
void assert_true() {
  run_assert_TRUE(0, 1);
}

NGP_TEST_INLINE
void assert_false() {
  run_assert_FALSE(1, 1);
}

NGP_TEST_INLINE
void assert_eq() {
  run_assert_EQ(1, 0);
}

NGP_TEST_INLINE
void assert_ne() {
  run_assert_NE(1, 1);
}

NGP_TEST_INLINE
void assert_lt() {
  run_assert_LT(1, 0);
}

NGP_TEST_INLINE
void assert_le() {
  run_assert_LE(1, 0);
}

NGP_TEST_INLINE
void assert_gt() {
  run_assert_GT(0, 1);
}

NGP_TEST_INLINE
void assert_ge() {
  run_assert_GE(0, 1);
}

NGP_TEST_INLINE
void run_assert_NEAR(const double a, const double b, const double tol) {
  NGP_ASSERT_NEAR(a, b, tol);
  NGP_EXPECT_NEAR(a, b, tol);
}

NGP_TEST_INLINE
void assert_near() {
  run_assert_NEAR(-1.0 - 2*eps, -1.0, eps);
}

template<class Func>
void execute(Func func, int howmany = 1) {
  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, howmany), func);
  Kokkos::fence();
}

#define DEFINE_DEVICE_FUNCTION(function)        \
  void device_##function(int howmany = 1) {     \
    execute(KOKKOS_LAMBDA(int /*i*/) {              \
        function();                             \
      }, howmany);                              \
  }

#define TEST_REPORT_ON_HOST(function)           \
  NGP_TEST_F(TestReportHost, function) {        \
    function();                                 \
  }

#define TEST_REPORT_ON_DEVICE(function)         \
  NGP_TEST_F(TestReportDevice, function) {      \
    device_##function();                        \
  }

#define TEST_REPORT_ON_HOST_AND_DEVICE(function)        \
  NGP_TEST_F(TestReportHostAndDevice, function) {       \
    function();                                         \
    device_##function();                                \
  }

#define TEST_REPORT_ON_HOST_AND_DEVICE_MANY(function)   \
  NGP_TEST_F(TestReportHostAndDeviceMany, function) {   \
    for(int i = 0; i < loopLen; ++i) {                  \
      function();                                       \
    }                                                   \
    device_##function(loopLen);                         \
  }

#define TEST_REPORTING(function)                \
  DEFINE_DEVICE_FUNCTION(function)              \
  TEST_REPORT_ON_HOST(function)                 \
  TEST_REPORT_ON_DEVICE(function)               \
  TEST_REPORT_ON_HOST_AND_DEVICE(function)      \
  TEST_REPORT_ON_HOST_AND_DEVICE_MANY(function)

TEST_REPORTING(expect_true)
TEST_REPORTING(assert_true)

TEST_REPORTING(expect_false)
TEST_REPORTING(assert_false)

TEST_REPORTING(expect_eq)
TEST_REPORTING(assert_eq)

TEST_REPORTING(expect_ne)
TEST_REPORTING(assert_ne)

TEST_REPORTING(expect_lt)
TEST_REPORTING(assert_lt)

TEST_REPORTING(expect_le)
TEST_REPORTING(assert_le)

TEST_REPORTING(expect_gt)
TEST_REPORTING(assert_gt)

TEST_REPORTING(expect_ge)
TEST_REPORTING(assert_ge)

TEST_REPORTING(expect_near)
TEST_REPORTING(assert_near)
