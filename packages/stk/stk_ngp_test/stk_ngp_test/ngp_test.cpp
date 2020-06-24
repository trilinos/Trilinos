#ifndef _ngp_test_cpp_
#define _ngp_test_cpp_

#include "GlobalReporter.hpp"
#include "Reporter.hpp"
#include <Kokkos_Core.hpp>
#include <stk_ngp_test/ngp_test.hpp>

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
  if(Kokkos::is_initialized()) Kokkos::finalize_all();
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
  return get_reporter()->get_capacity();
}

inline
void set_max_failure_reports_per_test(const int n) {
  get_reporter()->resize(n);
  get_device_reporter()->resize(n);
}


namespace internal {

NGP_TEST_INLINE
void add_failure(const char* condition, const char* location) {
  get_reporter()->add_failure(condition, location);
}

inline
int report_failures() {
  int numFailures = get_reporter()->report_failures();
  numFailures += get_device_reporter()->report_failures();
  return numFailures;
}

inline
void clear_failures() {
  get_reporter()->clear();
  get_device_reporter()->clear();
}

}

}

#endif

