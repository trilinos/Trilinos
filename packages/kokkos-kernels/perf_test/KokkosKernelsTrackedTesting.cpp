// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
//
// Created by Poliakoff, David Zoeller on 4/26/21.
//
#include <common/RAJAPerfSuite.hpp>
#include <common/Executor.hpp>
#include "sparse/tracked_testing.hpp"
#include <iostream>
#include <Kokkos_Core.hpp>
// For RPS version of BLAS Level-1 Tests
#include "blas/blas1/tracked_testing.hpp"
#include "blas/blas2/tracked_testing.hpp"
#include "blas/blas3/tracked_testing.hpp"
int main(int argc, char* argv[]) {
  {
    // argument parsing for setting input data at runtime

    std::string inputDataPath;
    if (argc == 1) {
      //    print_help();
      std::cout << "Please provide input data directory: --input-data "
                   "/PATH/TO/KOKKOS-KERNELS/INPUT/DATA"
                << std::endl;
      return 0;
    }

    for (int i = 0; i < argc; i++) {
      // if((strcmp(argv[i],"-v")==0)) {numVecs=atoi(argv[++i]); continue;}
      if ((strcmp(argv[i], "--input-data") == 0)) {
        i++;

        if (i == argc) {
          std::cerr << "Must pass desired input data after '--input-data'";
          exit(1);
        }
        inputDataPath = std::string(argv[i]);
        continue;
      }
    }

    test::set_input_data_path(inputDataPath);

    // set up Executor
    rajaperf::Executor exec(0, argv);
    // rajaperf::Executor exec(argc, argv);
    rajaperf::RunParams run_params(0, argv);
    // Initialize Kokkos
    Kokkos::initialize(argc, argv);

    Kokkos::print_configuration(std::cout);

    // sparse , spmv
    test::sparse::build_executor(exec, argc, argv, run_params);

    // All BLAS tests (Dot, Team Dot)
    test::blas::build_blas_executor(exec, argc, argv, run_params);

    test::blas2::build_blas2_executor(exec, argc, argv, run_params);

    test::blas3::build_blas3_executor(exec, argc, argv, run_params);

    exec.setupSuite();

    // STEP 3: Report suite run summary
    //         (enable users to catch errors before entire suite is run)
    exec.reportRunSummary(std::cout);

    // STEP 4: Execute suite
    exec.runSuite();

    // STEP 5: Generate suite execution reports
    exec.outputRunData();
  }
  Kokkos::finalize();
  return 0;
}
