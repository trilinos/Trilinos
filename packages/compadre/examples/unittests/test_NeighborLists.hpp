// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef TEST_NEIGHBORLISTS
#define TEST_NEIGHBORLISTS

#include "Compadre_NeighborLists.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace Compadre;

class NeighborListsTest: public ::testing::Test {
public:
    Kokkos::View<double**, host_execution_space> source_coords, target_coords;
    Kokkos::View<int**, host_execution_space> internal_2d_nl_view;

    NeighborListsTest( ) {
        // initialization
    }

    void SetUp( ) {

        internal_2d_nl_view = Kokkos::View<int**, host_execution_space>("2d_NL", 4, 5);

        // set number of neighbors
        internal_2d_nl_view(0,0)=1; 
        internal_2d_nl_view(1,0)=2; 
        internal_2d_nl_view(2,0)=3; 
        internal_2d_nl_view(3,0)=1; 
        // set neighbors
        internal_2d_nl_view(0,1)=0; 
        internal_2d_nl_view(1,1)=1; 
        internal_2d_nl_view(1,2)=2; 
        internal_2d_nl_view(2,1)=3; 
        internal_2d_nl_view(2,2)=4; 
        internal_2d_nl_view(2,3)=5; 
        internal_2d_nl_view(3,1)=6; 

    }

    void TearDown( ) {
        // after test completes
    }
};

TEST_F (NeighborListsTest, 2D_to_CompressedRow) {
    auto nl_f_2d = Convert2DToCompressedRowNeighborLists(internal_2d_nl_view);
    // what should be true
    ASSERT_EQ (4, nl_f_2d.getNumberOfTargets());
    ASSERT_EQ (1, nl_f_2d.getNumberOfNeighborsHost(0));
    ASSERT_EQ (2, nl_f_2d.getNumberOfNeighborsHost(1));
    ASSERT_EQ (3, nl_f_2d.getNumberOfNeighborsHost(2));
    ASSERT_EQ (1, nl_f_2d.getNumberOfNeighborsHost(3));
    ASSERT_EQ ((size_t)0, nl_f_2d.getRowOffsetHost(0));
    ASSERT_EQ ((size_t)1, nl_f_2d.getRowOffsetHost(1));
    ASSERT_EQ ((size_t)3, nl_f_2d.getRowOffsetHost(2));
    ASSERT_EQ ((size_t)6, nl_f_2d.getRowOffsetHost(3));
    ASSERT_EQ (0, nl_f_2d.getNeighborHost(0,0));
    ASSERT_EQ (1, nl_f_2d.getNeighborHost(1,0));
    ASSERT_EQ (2, nl_f_2d.getNeighborHost(1,1));
    ASSERT_EQ (3, nl_f_2d.getNeighborHost(2,0));
    ASSERT_EQ (4, nl_f_2d.getNeighborHost(2,1));
    ASSERT_EQ (5, nl_f_2d.getNeighborHost(2,2));
    ASSERT_EQ (6, nl_f_2d.getNeighborHost(3,0));
    ASSERT_EQ (3, nl_f_2d.getMaxNumNeighbors());
    ASSERT_EQ ((size_t)7, nl_f_2d.getTotalNeighborsOverAllListsHost());
}

#ifdef COMPADRE_EXTREME_DEBUG
TEST_F (NeighborListsTest, 2D_to_CompressedRow_EXTREME_DEBUG) {
    auto nl_f_2d = Convert2DToCompressedRowNeighborLists(internal_2d_nl_view);
    ASSERT_THROW(nl_f_2d.getNumberOfNeighborsHost(4), std::exception);
    ASSERT_THROW(nl_f_2d.getRowOffsetHost(4), std::exception);
    ASSERT_THROW(nl_f_2d.getNeighborHost(0,2), std::exception);
    ASSERT_THROW(nl_f_2d.getNeighborHost(4,0), std::exception);
}
#endif

#endif
