// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <stdlib.h>
#include <cstdio>
#include <random>

#include <Compadre_Config.h>
#include <Compadre_GMLS.hpp>
#include <Compadre_Evaluator.hpp>
#include <Compadre_PointCloudSearch.hpp>

#include "GMLS_Tutorial.hpp"

#ifdef COMPADRE_USE_MPI
#include <mpi.h>
#endif

#include <Kokkos_Timer.hpp>
#include <Kokkos_Core.hpp>

using namespace Compadre;

//! [Parse Command Line Arguments]

// called from command line
int main (int argc, char* args[]) {

#ifdef COMPADRE_USE_MPI
// initialized MPI (if avaialble) with command line arguments given
MPI_Init(&argc, &args);
#endif

// initializes Kokkos with command line arguments given
Kokkos::initialize(argc, args);

// becomes false if there is unwanted index in the filtered flags
bool all_passed = true;

// code block to reduce scope for all Kokkos View allocations
// otherwise, Views may be deallocating when we call Kokkos finalize() later
{
    // set the number of columns
    int num_cols = 50; // default 50 columns
    if (argc >= 3) {
        if (args[2] != NULL) {
            auto arg3toi = atoi(args[2]);
            if (isdigit(arg3toi)) {
                if (arg3toi > 0) {
                    num_cols = arg3toi;
                }
            }
        }
    }

    // set the number of flags
    int num_flags = 200; // default 200 flags
    if (argc >= 2) {
        if (args[1] != NULL) {
            auto arg2toi = atoi(args[1]);
            if (isdigit(arg2toi)) {
                if (arg2toi > 0) {
                    num_flags = arg2toi;
                }
            }
        }
    }
    //! [Parse Command Line Arguments]

    //! [Setting Up Data]
    Kokkos::Timer timer;
    Kokkos::Profiling::pushRegion("Setup Data");

    // create a 2D view of inputs
    Kokkos::View<int**, Kokkos::DefaultExecutionSpace> data_device("data", num_flags, num_cols);
    Kokkos::View<int**>::HostMirror data = Kokkos::create_mirror_view(data_device);

    // create a view of flags
    Kokkos::View<int*, Kokkos::DefaultExecutionSpace> flags_device("flags", num_flags);
    Kokkos::View<int*>::HostMirror flags = Kokkos::create_mirror_view(flags_device);

    //! [Setting Up Data]

    Kokkos::Profiling::popRegion();
    Kokkos::Profiling::pushRegion("Filter And Extract Data");

    //! [Filtering And Extracting Data]
    // create arbitrary data
    for (int i=0; i<num_flags; i++) {
        for (int j=0; j<num_cols; j++) {
            if ((i % 2) == 0) {
                data(i, j) = 1;
            } else {
                data(i, j) = 0;
            }
        }
    }
    // copy the data from host to device
    Kokkos::deep_copy(data_device, data);

    // create arbitrary flags
    int num_filtered_flags = 0; // number of filtered flags
    for (int i=0; i<num_flags; i++) {
        if ((i % 2) == 0) {
            flags(i) = 1;
            num_filtered_flags++;
        } else {
            flags(i) = 0;
        }
    }
    // copy the flags from host to device
    Kokkos::deep_copy(flags_device, flags);

    // Then call out the function to create view
    auto filtered_flags = filterViewByID<Kokkos::HostSpace>(flags_device, 1);
    auto extracted_data = Extract::extractViewByIndex<Kokkos::HostSpace>(data_device, filtered_flags);

    //! [Filtering Data]

    Kokkos::Profiling::popRegion();
    Kokkos::Profiling::pushRegion("Check Filtered And Extracted Data");

    //! [Checking Filtered And Extracted Data]

    if (filtered_flags.extent(0) != (size_t)num_filtered_flags) {
        all_passed = false;
        std::cout << "Failed - number of filtered flags not matched!" << filtered_flags.extent(0) << " " << num_filtered_flags << std::endl;
    }
    for (size_t i=0; i<filtered_flags.extent(0); i++) {
        if (filtered_flags(i) % 2 != 0) {
            all_passed = false;
            std::cout << "Failed - incorrect filtered flags " << filtered_flags(i) << std::endl;
        }
    }
    // All values inside extracted data should now be 1
    for (size_t i=0; i<extracted_data.extent(0); i++) {
        for (size_t j=0; j<extracted_data.extent(1); j++) {
            if (extracted_data(i, j) != 1) {
                all_passed = false;
                std::cout << "Failed - incorrect values in extracted view at index " << i << " " << j << " " << extracted_data(i, j) << std::endl;
            }
        }
    }

    //! [Checking Filtered And Extracted Data]
    // stop timing comparison loop
    Kokkos::Profiling::popRegion();
    //! [Finalize Program]

} // end of code block to reduce scope, causing Kokkos View de-allocations
// otherwise, Views may be deallocating when we call Kokkos finalize() later

// finalize Kokkos and MPI (if available)
Kokkos::finalize();
#ifdef COMPADRE_USE_MPI
MPI_Finalize();
#endif

// output to user that test passed or failed
if (all_passed) {
    fprintf(stdout, "Passed test \n");
    return 0;
} else {
    fprintf(stdout, "Failed test \n");
    return -1;
}

} // main

//! [Finalize Program]
