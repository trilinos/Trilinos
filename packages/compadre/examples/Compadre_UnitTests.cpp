// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#include <gtest/gtest.h>
#include <Compadre_KokkosParser.hpp>
#include "unittests/test_XYZ.hpp"
#include "unittests/test_NeighborLists.hpp"
#include "unittests/test_PointCloudSearch.hpp"
#include "unittests/test_LinearAlgebra.hpp"
#include "unittests/test_Targets.hpp"
#ifdef COMPADRE_USE_MPI
#include <mpi.h>
#endif

// this provides main(),
// but all other tests come from ./unittests/*.cpp
int main(int argc, char **argv) {

    // initializes MPI (if available) with command line arguments given
    #ifdef COMPADRE_USE_MPI
    MPI_Init(&argc, &argv);
    #endif
    // this MPI init is rerun when using the "threadsafe" test style below

    ::testing::InitGoogleTest(&argc, argv);

    // initializes kokkos
    Kokkos::initialize(argc, argv);

    // execute all tests
    ::testing::GTEST_FLAG(filter) = "*";
    int sig = RUN_ALL_TESTS();

    // finalize Kokkos and MPI (if available)
    Kokkos::finalize();

    // finialize MPI (if available)
    #ifdef COMPADRE_USE_MPI
    MPI_Finalize();
    #endif
    return sig; 
}
