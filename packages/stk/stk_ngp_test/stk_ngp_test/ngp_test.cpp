#include "GlobalReporter.hpp"
#include "Reporter.hpp"
#include <Kokkos_Core.hpp>
#include <stk_ngp_test/ngp_test.hpp>

namespace ngp_testing {

NgpTestEnvironment::NgpTestEnvironment(int* argc, char** argv) {
  Kokkos::initialize(*argc, argv);
  initialize_reporters();
  ::testing::InitGoogleTest(argc, argv);
  Kokkos::push_finalize_hook(finalize_reporters);
}

NgpTestEnvironment::~NgpTestEnvironment() {
  Kokkos::finalize_all();
}

int NgpTestEnvironment::run_all_tests() {
  return RUN_ALL_TESTS();
}

void Test::SetUp() {
  internal::clear_failures();
  NGPSetUp();
}

void Test::TearDown() {
  int numReports = internal::report_failures();
  if(numReports > 0) {
    set_failure_in_teardown();
  }
  NGPTearDown();
}

int get_max_failure_reports_per_test() {
  return get_reporter()->get_capacity();
}

void set_max_failure_reports_per_test(const int n) {
  get_reporter()->resize(n);
  get_device_reporter()->resize(n);
}


namespace internal {

NGP_TEST_FUNCTION
void add_failure(const char* condition, const char* location) {
  get_reporter()->add_failure(condition, location);
}

int report_failures() {
  int numFailures = get_reporter()->report_failures();
  numFailures += get_device_reporter()->report_failures();
  return numFailures;
}

void clear_failures() {
  get_reporter()->clear();
  get_device_reporter()->clear();
}

}

}
