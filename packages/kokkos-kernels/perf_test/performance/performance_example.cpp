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
  Notes on performance_demo test

  This is intended to be a minimal example of using the new YAML archiver.
  The times and residuals are dummy values to mimic a real test.
  First time running the test it should create the new yaml archive with 1
  entry. Subsequent runs will validate the values and pass.

  To play around with this change the time1 value and run again to see it fail.
  Or see performance_validate which runs through all the things the archiver
  does.
*/

#include "Kokkos_Performance.hpp"  // provides performance archiver

bool run_example() {
  // Some tests are run and produce some times...
  double time1 = 10.0;
  double time2 = 13.3;

  // and they produce some results...
  double residual        = 0.001;
  int some_exact_counter = 22;

  // set up some user options
  std::string archiveName("performance_example.yaml");  //  name of the archive
  std::string testName = "performance_example";         // name of test
  std::string hostName;                                 // optional hostname - auto detected if blank
  double tolerance = 0.1;                               // for residual and times

  using KokkosKernels::Performance;

  // Create an archiver - steps are create, fill with members, then run
  Performance archiver;

  // Example of how to set customized machine config - to be developed
  // Change to generate new entries in the yaml under MachineConfiguration
  archiver.set_machine_config("Kokkos Config", "some node type");

  // Fill config
  archiver.set_config("MPI_Ranks", 2);
  archiver.set_config("Teams", 1);    // just arbitrary right now
  archiver.set_config("Threads", 1);  // just arbitrary right now
  archiver.set_config("Filename",
                      "somefilename");  // arbitrary - example of a string

  // Fill results
  archiver.set_result("Time1", time1, tolerance);
  archiver.set_result("Time2", time2, tolerance);
  archiver.set_result("Residual", residual, tolerance);
  archiver.set_result("Counter", some_exact_counter);  // must match exactly

  // run it
  Performance::Result result = archiver.run(archiveName, testName, hostName);

  // print the yaml file for inspection
  Performance::print_archive(archiveName);

  // Print results
  switch (result) {
    case Performance::Passed: std::cout << "Archiver Passed" << std::endl; break;
    case Performance::Failed: std::cout << "Archiver Failed" << std::endl; break;
    case Performance::NewMachine: std::cout << "Archiver Passed. Adding new machine entry." << std::endl; break;
    case Performance::NewConfiguration:
      std::cout << "Archiver Passed. Adding new machine configuration." << std::endl;
      break;
    case Performance::NewTest: std::cout << "Archiver Passed. Adding new test entry." << std::endl; break;
    case Performance::NewTestConfiguration:
      std::cout << "Archiver Passed. Adding new test entry configuration." << std::endl;
      break;
    case Performance::UpdatedTest: std::cout << "Archiver Passed. Updating test entry." << std::endl; break;
    default: throw std::logic_error("Unexpected result code."); break;
  }

  return (result != Performance::Failed);
}

int main(int argc, char *argv[]) {
  bool success = run_example();

  if (success) {
    std::cout << "End Result: TEST PASSED" << std::endl;
  } else {
    std::cout << "End Result: TEST FAILED" << std::endl;
  }

  return EXIT_SUCCESS;
}
