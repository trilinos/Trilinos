// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef TEST_POINTCLOUDSEARCH
#define TEST_POINTCLOUDSEARCH

#include "Compadre_PointCloudSearch.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace Compadre;

class PointCloudSearchTest: public ::testing::Test {
public:
    Kokkos::View<double**, host_execution_space> source_coords, target_coords;
    const int number_source_coords;
    const int number_target_coords;

    Kokkos::View<int**, host_execution_space> internal_2d_nl_view;

    PointCloudSearch<decltype(source_coords)> point_cloud_search;

    PointCloudSearchTest( ) : number_source_coords(5), number_target_coords(2),
                              point_cloud_search(decltype(point_cloud_search)(source_coords)) {
        // initialization
    }

    void SetUp( ) {

        //
        // 1D setup
        //

        // sets up reference neighbor list
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

        // coordinates of source sites
        source_coords = Kokkos::View<double**, host_execution_space>("source coordinates", 
                number_source_coords, 1);
        
        // coordinates of target sites
        target_coords = Kokkos::View<double**, host_execution_space>("target coordinates", 
                number_target_coords, 1);

        const double h_s = 1.0/(number_source_coords-1);
        // set up 1D problem
        for (int i=0; i<number_source_coords; ++i) {
            source_coords(i,0) = i*h_s;
        }
        target_coords(0,0) = 1.0/3.0;
        target_coords(1,0) = 2.0/3.0;

        point_cloud_search = CreatePointCloudSearch(source_coords, 1 /*dimension*/);

    }

    void TearDown( ) {
        // after test completes
    }
};

TEST_F (PointCloudSearchTest, 1D_Radius_Search) {
    // Empty views to be resized/filled
    Kokkos::View<int*, host_execution_space> neighbor_lists("neighbor lists", 0);
    Kokkos::View<int*, host_execution_space> number_of_neighbors_list("number of neighbor lists", 
            number_target_coords); 
    Kokkos::View<double*, host_execution_space> epsilon("h supports", 
            number_target_coords);

    // This dry run populates number_of_neighbors_list with neighborhood sizes
    size_t storage_size = 
            point_cloud_search.generateCRNeighborListsFromRadiusSearch(true /* dry run */,
                    target_coords, neighbor_lists, number_of_neighbors_list, epsilon, 0.2 /*radius*/);

    // Resize based on dry-run
    Kokkos::resize(neighbor_lists, storage_size);

    // Search again, now that we know that there is enough room to store the results
    point_cloud_search.generateCRNeighborListsFromRadiusSearch(false /* dry run */,
            target_coords, neighbor_lists, number_of_neighbors_list, epsilon, 0.2 /*radius*/);

    auto nla(CreateNeighborLists(neighbor_lists, number_of_neighbors_list));

    // 2 targets
    ASSERT_EQ(2, nla.getNumberOfTargets());
    // 2 neighbors for target 0
    ASSERT_EQ(2, nla.getNumberOfNeighborsHost(0));
    // 2 neighbors for target 1
    ASSERT_EQ(2, nla.getNumberOfNeighborsHost(1));
    // radius is 0.2
    ASSERT_DOUBLE_EQ(0.2, epsilon(0));
    ASSERT_DOUBLE_EQ(0.2, epsilon(1));
    // index of first neighbor for target 0 is 1
    ASSERT_EQ(1, nla.getNeighborHost(0,0));
    // index of second neighbor for target 0 is 2
    ASSERT_EQ(2, nla.getNeighborHost(0,1));
    // index of first neighbor for target 1 is 3
    ASSERT_EQ(3, nla.getNeighborHost(1,0));
    // index of second neighbor for target 1 is 2
    ASSERT_EQ(2, nla.getNeighborHost(1,1));
}

TEST_F (PointCloudSearchTest, 1D_Dynamic_Search) {
    // Empty views to be resized/filled
    Kokkos::View<int*, host_execution_space> neighbor_lists("neighbor lists", 0);
    Kokkos::View<int*, host_execution_space> number_of_neighbors_list("number of neighbor lists", 
            number_target_coords); 
    Kokkos::View<double*, host_execution_space> epsilon("h supports", 
            number_target_coords);

    double epsilon_multiplier = 1.5;
    // This dry run populates number_of_neighbors_list with neighborhood sizes
    size_t storage_size = 
            point_cloud_search.generateCRNeighborListsFromKNNSearch(true /*dry run*/, 
                    target_coords, neighbor_lists, 
                    number_of_neighbors_list, epsilon, 
                    3 /*min_neighbors*/, epsilon_multiplier);

    // Resize based on dry-run
    Kokkos::resize(neighbor_lists, storage_size);

    // Search again, now that we know that there is enough room to store the results
    point_cloud_search.generateCRNeighborListsFromKNNSearch(false /*dry run*/, 
                    target_coords, neighbor_lists, 
                    number_of_neighbors_list, epsilon, 
                    3 /*min_neighbors*/, epsilon_multiplier);

    auto nla(CreateNeighborLists(neighbor_lists, number_of_neighbors_list));

    // 2 targets
    ASSERT_EQ(2, nla.getNumberOfTargets());
    // 4 neighbors for target 0
    ASSERT_EQ(4, nla.getNumberOfNeighborsHost(0));
    // 4 neighbors for target 1
    ASSERT_EQ(4, nla.getNumberOfNeighborsHost(1));
    // radius is 0.5
    ASSERT_DOUBLE_EQ(0.5, epsilon(0));
    ASSERT_DOUBLE_EQ(0.5, epsilon(1));

    std::set<int> t0_neighbors, t1_neighbors;
    for (int j=0; j<4; ++j) {
        t0_neighbors.insert(nla.getNeighborHost(0,j));
        t1_neighbors.insert(nla.getNeighborHost(1,j));
    }

    // index of first neighbor for target 0 is 1
    // index of second neighbor for target 0 is 2
    // index of third neighbor for target 0 is 0
    // index of fourth neighbor for target 0 is 3
    // however, neighbor lists are not sorted (except for first entry)
    //ASSERT_EQ(1, nla.getNeighborHost(0,0));
    //ASSERT_EQ(2, nla.getNeighborHost(0,1));
    //ASSERT_EQ(0, nla.getNeighborHost(0,2));
    //ASSERT_EQ(3, nla.getNeighborHost(0,3));
    ASSERT_TRUE(t0_neighbors.find(1) != t0_neighbors.end());
    ASSERT_TRUE(t0_neighbors.find(2) != t0_neighbors.end());
    ASSERT_TRUE(t0_neighbors.find(0) != t0_neighbors.end());
    ASSERT_TRUE(t0_neighbors.find(3) != t0_neighbors.end());

    // index of first neighbor for target 1 is 3
    // index of second neighbor for target 1 is 2
    // index of third neighbor for target 1 is 4
    // index of fourth neighbor for target 1 is 1
    // however, neighbor lists are not sorted (except for first entry)
    //ASSERT_EQ(3, nla.getNeighborHost(1,0));
    //ASSERT_EQ(2, nla.getNeighborHost(1,1));
    //ASSERT_EQ(4, nla.getNeighborHost(1,2));
    //ASSERT_EQ(1, nla.getNeighborHost(1,3));
    ASSERT_TRUE(t1_neighbors.find(3) != t1_neighbors.end());
    ASSERT_TRUE(t1_neighbors.find(2) != t1_neighbors.end());
    ASSERT_TRUE(t1_neighbors.find(4) != t1_neighbors.end());
    ASSERT_TRUE(t1_neighbors.find(1) != t1_neighbors.end());
}

#endif
