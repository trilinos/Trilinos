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

/*
  Notes on performance_validate test

  This is intended to stress test the Kokkos_Performance.hpp Performance class.
  It runs through a series of alterations to mimic running the same test
  repeatedly with conditions changing. It validates all the features of the
  yaml archiver checking for proper pass/fails.

  Note this test deletes the yaml archive before running which is a special
  case for testing and not something normally done.

  See the update() method below for the sequence of events that are tested.
*/

#include "Kokkos_Performance.hpp"  // provides performance archiver

// a manager class to run the series of tests
class TestManager {
 public:
  TestManager()
      : is_error(false),                           // tracks if there was an unexpected error
        update_index(0),                           // utility for running through a sequence of tests
        archiveName("performance_validate.yaml"),  // the yaml file used
        testName("test1"),                         // the initial test name appearing in the archive
        hostName(""),                              // blank means auto detect - otherwise sets host name
        changeCompiler(""),                        // override option for simulating different machine config
        filename("somefilename"),                  // arbitrary configuration parameter
        mpi_ranks(1),                              // arbitrary configuration parameter
        teams(1),                                  // arbitrary configuration parameter
        threads(1),                                // arbitrary configuration parameter
        tolerance(0.1),                            // tolerance on time values
        time1(10.0),                               // arbitrary time value
        time2(13.3),                               // arbitrary time value
        niters(44),                                // arbitrary result value
        residual(0.001),                           // arbitrary result value
        bExtraParameters(false),                   // used to validate and test the update mechanism
        extra_time(22.3),                          // extra time we'll add to preexisiting data
        extra_result(40),                          // extra result we'll add to preexisting data
        bExactParameters(false),                   // add exact match parameters
        exact_int(5),                              // test for exact match int
        exact_string("somestring") {               // test for exact match string

    using KokkosKernels::Performance;

    // normally we would not delete the yaml since we want to compare but this
    // test is designed to run through a series of steps including starting
    // from nothing, so write out a blank node to start clean
    Performance::erase_archive(archiveName);

    // now loop all the tests in update until we run out of things to try
    while (update()) {
    }

    // print the archive for inspection
    Performance::print_archive(archiveName);
  }

  bool isError() const { return is_error; }

  int update() {
    using KokkosKernels::Performance;

    // some arbitrary delta used to shift values up and down during test
    const double small_time_change = 0.01;  // a small value not exceeding tolerance
    const double big_time_change   = 10.0;  // a big value exceeding tolerance
    const int integer_change       = 2;     // some non 0 value

    Performance::Result expected_result = Performance::Unknown;
    switch (++update_index) {
      case 1:
        std::cout << "Create a new machine entry with no prexisting archive." << std::endl;
        expected_result = Performance::NewMachine;
        break;
      case 2:
        std::cout << "Run again with same data - we expect it to pass." << std::endl;
        expected_result = Performance::Passed;
        break;
      case 3:
        std::cout << "Make a small time change not exceeding tolerance should "
                     "still pass."
                  << std::endl;
        time1 += small_time_change;
        expected_result = Performance::Passed;
        break;
      case 4:
        std::cout << "Make a big time change exceeding tolerance to trigger failure." << std::endl;
        time1 -= small_time_change;
        time1 += big_time_change;
        expected_result = Performance::Failed;
        break;
      case 5:
        std::cout << "Run with the original time again and should pass with a match." << std::endl;
        time1 -= big_time_change;
        expected_result = Performance::Passed;
        break;
      case 6:
        std::cout << "Change the host name to 'custom hostname 1' which "
                     "creates a new entry."
                  << std::endl;
        hostName        = "custom hostname 1";
        expected_result = Performance::NewMachine;
        break;
      case 7:
        std::cout << "Change the host name to 'custom hostname 2' which "
                     "creates another entry."
                  << std::endl;
        hostName        = "custom hostname 2";
        expected_result = Performance::NewMachine;
        break;
      case 8:
        std::cout << "Change back to default host name so we should match data." << std::endl;
        hostName        = "";
        expected_result = Performance::Passed;
        break;
      case 9:
        std::cout << "Change the test configuration to create a new config type." << std::endl;
        mpi_ranks += integer_change;
        expected_result = Performance::NewTestConfiguration;
        break;
      case 10:
        std::cout << "Go back to original configuration which we should match." << std::endl;
        mpi_ranks -= integer_change;
        expected_result = Performance::Passed;
        break;
      case 11:
        std::cout << "Create a new test name which will be a new entry." << std::endl;
        testName        = "test2";
        expected_result = Performance::NewTest;
        break;
      case 12:
        std::cout << "Change test name and the niters result value at same "
                     "time which creates a new entry."
                  << std::endl;
        testName = "test3";
        niters += integer_change;
        expected_result = Performance::NewTest;
        break;
      case 13:
        std::cout << "Still on test3 change niters back so we should fail - "
                     "doesn't match test3."
                  << std::endl;
        niters -= integer_change;
        expected_result = Performance::Failed;  // because still on test3 which had +5
        break;
      case 14:
        std::cout << "Changing back to test1 we should now pass because niters "
                     "matches test1."
                  << std::endl;
        testName        = "test1";
        expected_result = Performance::Passed;  // should work because niters matches
        break;
      case 15:
        std::cout << "Simulate a different compiler by overriding the machine "
                     "configuration to create a new entry. Also change test "
                     "config to mix it up."
                  << std::endl;
        changeCompiler = "SpecialCompiler";
        teams += integer_change;
        expected_result = Performance::NewConfiguration;
        break;
      case 16:
        std::cout << "Run again - should match our machine configuration." << std::endl;
        expected_result = Performance::Passed;
        break;
      case 17:
        std::cout << "Change test config back but stay on this machine "
                     "configuration - creates a new entry."
                  << std::endl;
        teams -= integer_change;
        expected_result = Performance::NewTestConfiguration;
        break;
      case 18:
        std::cout << "Go back to the original machine configuration - we "
                     "should succeed."
                  << std::endl;
        changeCompiler  = "";
        expected_result = Performance::Passed;
        break;
      case 19:
        std::cout << "Add an additional parameters to the test entry. This is "
                     "allowed."
                  << std::endl;
        bExtraParameters = true;
        expected_result  = Performance::UpdatedTest;
        break;
      case 20:
        std::cout << "Run again with extra parameters. Should be fine and pass." << std::endl;
        expected_result = Performance::Passed;
        break;
      case 21:
        std::cout << "Change an extra parameter - should fail." << std::endl;
        extra_time += big_time_change;
        expected_result = Performance::Failed;
        break;
      case 22:
        std::cout << "Try test2 without extra parameters and make sure it "
                     "still passes."
                  << std::endl;
        bExtraParameters = false;
        testName         = "test2";
        expected_result  = Performance::Passed;
        break;
      case 23:
        std::cout << "Back to test1 with extra parameters still off - this is "
                     "not allowed."
                  << std::endl;
        testName = "test1";
        extra_time -= big_time_change;  // restore the value so it's correct - we fail
                                        // here because we don't have the parameters, not
                                        // because they are wrong
        expected_result = Performance::Failed;
        break;
      case 24:
        std::cout << "Restore extra parameters and check it's ok." << std::endl;
        bExtraParameters = true;
        expected_result  = Performance::Passed;
        break;
      case 25:
        std::cout << "Turn on exact parameters." << std::endl;
        bExactParameters = true;
        expected_result  = Performance::UpdatedTest;
        break;
      case 26:
        std::cout << "Change exact string - should fail." << std::endl;
        exact_string    = "anotherstring";
        expected_result = Performance::Failed;
        break;
      case 27:
        std::cout << "Change exact string back. Should pass." << std::endl;
        exact_string    = "somestring";
        expected_result = Performance::Passed;
        break;
      case 28:
        std::cout << "Change exact int - should fail." << std::endl;
        exact_int += integer_change;
        expected_result = Performance::Failed;
        break;
      case 29:
        std::cout << "Change it back - should pass." << std::endl;
        exact_int -= integer_change;
        expected_result = Performance::Passed;
        break;
      default:
        std::cout << "All tests completed!" << std::endl << std::endl;
        return false;  // we finished all the tests
        break;
    }

    Performance archiver;
    // this step would normally not exist as we'd want the correct machine
    // config but to test we override a setting to simulate running a different
    // config
    if (changeCompiler != "") {
      archiver.set_machine_config("Compiler", changeCompiler);
    }

    // set up test config
    archiver.set_config("MPI Ranks", mpi_ranks);
    archiver.set_config("Teams", teams);
    archiver.set_config("Threads", threads);
    archiver.set_config("Filename", filename);

    // set up test results
    archiver.set_result("Time1", time1, tolerance);
    archiver.set_result("Time2", time2, tolerance);
    archiver.set_result("Iterations", niters, niters > 0 ? niters - 1 : 0, niters + 1);
    archiver.set_result("Residual", residual, tolerance);

    // add optionals if turned on
    if (bExtraParameters) {
      archiver.set_result("Extra time", extra_time, tolerance);
      archiver.set_result("Extra result", extra_result, extra_result - 3, extra_result + 3);
    }
    if (bExactParameters) {
      // if we don't set tolerance it will be assigned as an exact match
      // currently internally using the key character '!' to identify these
      archiver.set_result("Exact string", exact_string);
      archiver.set_result("Exact int", exact_int);
    }

    // update the archive
    Performance::Result result = archiver.run(archiveName, testName, hostName);

    // Get results string
    std::string result_string;
    switch (result) {
      case Performance::Passed: result_string = "Archiver Passed"; break;
      case Performance::Failed: result_string = "Archiver Failed"; break;
      case Performance::NewMachine: result_string = "Archiver Passed. Adding new machine entry."; break;
      case Performance::NewConfiguration: result_string = "Archiver Passed. Adding new machine configuration."; break;
      case Performance::NewTest: result_string = "Archiver Passed. Adding new test entry."; break;
      case Performance::NewTestConfiguration:
        result_string = "Archiver Passed. Adding new test entry configuration.";
        break;
      case Performance::UpdatedTest: result_string = "Archiver Passed. Updating test entry."; break;
      default: throw std::logic_error("Unexpected result code."); break;
    }

    // validate the returned result is as expected
    // note that a Fail event from the archiver may be the expected result
    if (result != expected_result) {
      std::cout << "  Invalid result from archive! " << result_string << std::endl << std::endl;
      is_error = true;  // we failed!
      return false;     // we can stop now - something went wrong
    } else {
      std::cout << "  Correct: " << result_string << std::endl << std::endl;
    }

    return true;
  }

 private:
  bool is_error;               // tracks if there was an unexpected error
  int update_index;            // utility for running through a sequence of tests
  std::string archiveName;     // the yaml file used
  std::string testName;        // the test name appearing in the archive
  std::string hostName;        // blank means auto detect - otherwise sets host name
  std::string changeCompiler;  // override option to simulate different machine config
  std::string filename;        // arbitrary configuration parameter
  int mpi_ranks;               // arbitrary configuration parameter
  int teams;                   // arbitrary configuration parameter
  int threads;                 // arbitrary configuration parameter
  double tolerance;            // tolerance on time values
  double time1;                // arbitrary time value
  double time2;                // arbitrary time value
  int niters;                  // arbitrary result value
  double residual;             // arbitrary result value
  bool bExtraParameters;       // used to validate and test the update mechanism
  double extra_time;           // extra time we'll add to preexisiting data
  int extra_result;            // extra result we'll add to preexisting data
  bool bExactParameters;       // used to turn on exact match parameters
  int exact_int;               // test for exact match int
  std::string exact_string;    // test for exact match string
};

int main(int argc, char *argv[]) {
  TestManager manager;  // loops through a series of changes and archive each
  bool success = !manager.isError();

  if (success) {
    std::cout << "End Result: TEST PASSED" << std::endl;
  } else {
    std::cout << "End Result: TEST FAILED" << std::endl;
  }

  return EXIT_SUCCESS;
}
