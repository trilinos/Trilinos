// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#include <vector>
#include <map>
#include <random>

#include <Compadre_Config.h>
#include <Compadre_GMLS.hpp>
#include <Compadre_PointCloudSearch.hpp>

#ifdef COMPADRE_USE_MPI
#include <mpi.h>
#endif

template<typename T1, typename T2>
std::pair<T2,T1> swap_pair(const std::pair<T1,T2> &item)
{
    return std::pair<T2,T1>(item.second, item.first);
}

template<typename T1, typename T2>
std::multimap<T2,T1> invert_map(const std::map<T1,T2> &original_map)
{
    std::multimap<T2,T1> new_map;
    std::transform(original_map.begin(), original_map.end(), std::inserter(new_map, new_map.begin()), swap_pair<T1,T2>);
    return new_map;
}

using namespace Compadre;

int main (int argc, char* args[]) {

// initializes MPI (if available) with command line arguments given
#ifdef COMPADRE_USE_MPI
MPI_Init(&argc, &args);
#endif

// initializes Kokkos with command line arguments given
Kokkos::initialize(argc, args);

// becomes false if the computed solution not within the failure_threshold of the actual solution
bool all_passed = true;

{
    // seed random generator
    srand(1234);
    Kokkos::Timer timer;

    int search_type = 0; // 0 - radius, 1 - knn
    if (argc >= 5) {
        int arg5toi = atoi(args[4]);
        if (arg5toi >= 0) {
            search_type = arg5toi;
        }
    }
    bool do_radius_search = (search_type==0);
    bool do_knn_search = (search_type==1);
    printf("do_radius_search: %d\n", do_radius_search);
    printf("do_knn_search: %d\n", do_knn_search);

    // check if 4 arguments are given from the command line
    //  multiplier times h spacing for search
    double h_multiplier = 3.0; // dimension 3 by default
    if (argc >= 4) {
        double arg4tof = atof(args[3]);
        if (arg4tof > 0) {
            h_multiplier = arg4tof;
        }
    }

    // check if 3 arguments are given from the command line
    //  set the number of target sites where we will reconstruct the target functionals at
    int number_target_coords = 200; // 200 target sites by default
    if (argc >= 3) {
        int arg3toi = atoi(args[2]);
        if (arg3toi > 0) {
            number_target_coords = arg3toi;
        }
    }
    
    // check if 2 arguments are given from the command line
    //  set the number of dimensions for points and the search
    int dimension = 3; // 3D by default
    if (argc >= 2) {
        int arg2toi = atoi(args[1]);
        if (arg2toi > 0) {
            dimension = arg2toi;
        }
    }
    
    //// minimum neighbors for unisolvency is the same as the size of the polynomial basis 
    //const int min_neighbors = Compadre::GMLS::getNP(order, dimension);
    
    // approximate spacing of source sites
    double h_spacing = 1./std::pow(number_target_coords, 0.50 + 0.2*(dimension==2));
    int n_neg1_to_1 = 2*(1/h_spacing) + 1; // always odd
    
    // number of source coordinate sites that will fill a box of [-1,1]x[-1,1]x[-1,1] with a spacing approximately h
    int number_source_coords = (n_neg1_to_1+1)*(n_neg1_to_1+1);
    if (dimension==3) number_source_coords *= (n_neg1_to_1+1);
    printf("target coords: %d\n", number_target_coords);
    printf("source coords: %d\n", number_source_coords);
    
    // coordinates of source sites
    Kokkos::View<double**, Kokkos::DefaultHostExecutionSpace> source_coords("source coordinates", 
            number_source_coords, 3);
    
    // coordinates of target sites
    Kokkos::View<double**, Kokkos::DefaultHostExecutionSpace> target_coords("target coordinates", number_target_coords, 3);

    // fill source coordinates with a uniform grid
    int source_index = 0;
    double this_coord[3] = {0,0,0};
    for (int i=-n_neg1_to_1/2; i<n_neg1_to_1/2+1; ++i) {
        this_coord[0] = i*h_spacing;
        for (int j=-n_neg1_to_1/2; j<n_neg1_to_1/2+1; ++j) {
            this_coord[1] = j*h_spacing;
            for (int k=-n_neg1_to_1/2; k<n_neg1_to_1/2+1; ++k) {
                this_coord[2] = k*h_spacing;
                if (dimension==3) {
                    source_coords(source_index,0) = this_coord[0]; 
                    source_coords(source_index,1) = this_coord[1]; 
                    source_coords(source_index,2) = this_coord[2]; 
                    source_index++;
                }
            }
            if (dimension==2) {
                source_coords(source_index,0) = this_coord[0]; 
                source_coords(source_index,1) = this_coord[1]; 
                source_coords(source_index,2) = 0;
                source_index++;
            }
        }
        if (dimension==1) {
            source_coords(source_index,0) = this_coord[0]; 
            source_coords(source_index,1) = 0;
            source_coords(source_index,2) = 0;
            source_index++;
        }
    }
    
    // fill target coords somewhere inside of [-0.5,0.5]x[-0.5,0.5]x[-0.5,0.5]
    for(int i=0; i<number_target_coords; i++){
    
        // first, we get a uniformly random distributed direction
        double rand_dir[3] = {0,0,0};
    
        for (int j=0; j<dimension; ++j) {
            // rand_dir[j] is in [-0.5, 0.5]
            rand_dir[j] = ((double)rand() / (double) RAND_MAX) - 0.5;
        }
    
        // then we get a uniformly random radius
        for (int j=0; j<dimension; ++j) {
            target_coords(i,j) = rand_dir[j];
        }
    
    }

    if (do_radius_search) {
        timer.reset();
        { // 1D compressed row neighbor lists
            printf("\n1D compressed row neighbor lists:\n");
            // Point cloud construction for neighbor search
            // CreatePointCloudSearch constructs an object of type PointCloudSearch, but deduces the templates for you
            auto point_cloud_search(CreatePointCloudSearch(source_coords, dimension));

            Kokkos::View<int*, Kokkos::DefaultHostExecutionSpace> neighbor_lists("neighbor lists", 
                    number_target_coords); // first column is # of neighbors
            Kokkos::View<int*, Kokkos::DefaultHostExecutionSpace> number_neighbors_list("number of neighbors list", 
                    number_target_coords); // first column is # of neighbors
            
            // each target site has a window size
            Kokkos::View<double*, Kokkos::DefaultHostExecutionSpace> epsilon("h supports", number_target_coords);
            Kokkos::View<double*, Kokkos::DefaultHostExecutionSpace> random_vec("random vector", number_target_coords);

            Kokkos::parallel_for("random perturbation of search size", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>
                    (0,number_target_coords), [&](const int i) {
                // value in [-.25, .25]
                random_vec(i) = 0.5 * ((double)rand() / (double) RAND_MAX) - 0.5;
                epsilon(i) = (h_multiplier + random_vec(i)) * h_spacing;
            });
            Kokkos::fence();
            
            // query the point cloud to generate the neighbor lists using a radius search

            // start with a dry-run, but without enough room to store the results
            size_t total_num_neighbors = point_cloud_search.generateCRNeighborListsFromRadiusSearch(true /* dry run */,
                target_coords, neighbor_lists, number_neighbors_list, epsilon);
            printf("total num neighbors: %lu\n", total_num_neighbors);

            // resize neighbor lists to be large enough to hold the results
            Kokkos::resize(neighbor_lists, total_num_neighbors);

            // search again, now that we know that there is enough room to store the results
            point_cloud_search.generateCRNeighborListsFromRadiusSearch(false /* dry run */,
                target_coords, neighbor_lists, number_neighbors_list, epsilon);

            auto nla(CreateNeighborLists(neighbor_lists, number_neighbors_list));

            double radius_search_time = timer.seconds();
            printf("nanoflann search time: %f s\n", radius_search_time);

            // convert point cloud search to vector of maps
            timer.reset();
            std::vector<std::map<int, double> > point_cloud_neighbor_list(number_target_coords, std::map<int, double>());
            Kokkos::parallel_for("point search conversion", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>
                    (0,number_target_coords), [&](const int i) {
                auto &point_cloud_neighbor_list_i = const_cast<std::map<int, double>& >(point_cloud_neighbor_list[i]);
                for (int j=0; j<nla.getNumberOfNeighborsHost(i); ++j) {
                    point_cloud_neighbor_list_i[nla.getNeighborHost(i,j)] = 1.0;
                }
            });
            double point_cloud_convert_time = timer.seconds();
            printf("point cloud convert time: %f s\n", point_cloud_convert_time);
            Kokkos::fence();

            // loop over targets, finding all of their neighbors by brute force
            timer.reset();
            std::vector<std::map<int, double> > brute_force_neighbor_list(number_target_coords, std::map<int, double>());
            Kokkos::parallel_for("brute force radius search", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>
                    (0,number_target_coords), [&](const int i) {
                auto &brute_force_neighbor_list_i = const_cast<std::map<int, double>& >(brute_force_neighbor_list[i]);
                for (int j=0; j<number_source_coords; ++j) {
                    double dist = 0;
                    for (int k=0; k<dimension; ++k) {
                        dist += (target_coords(i,k)-source_coords(j,k))*(target_coords(i,k)-source_coords(j,k));
                    }
                    dist = std::sqrt(dist);

                    if (dist<(h_multiplier + random_vec(i))*h_spacing) {
                        brute_force_neighbor_list_i[j]=dist;
                    }
                }
            });
            double brute_force_search_time = timer.seconds();
            printf("brute force search time: %f s\n", brute_force_search_time);
            Kokkos::fence();

            timer.reset();
            // check that all neighbors in brute force list are in point cloud search list
            int num_passed = 0;
            Kokkos::parallel_reduce("check brute force search in point cloud search", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>
                    (0,number_target_coords), [&](const int i, int& t_num_passed) {
                bool all_found = true;
                for (auto it=brute_force_neighbor_list[i].begin(); it!=brute_force_neighbor_list[i].end(); ++it) {
                    if (point_cloud_neighbor_list[i].count(it->first)!=1) {
                        all_found = false;
                    }
                }
                if (all_found) t_num_passed++;
            }, Kokkos::Sum<int>(num_passed));
            Kokkos::fence();
            if (num_passed != number_target_coords) {
                printf("Brute force neighbor not found in point cloud list\n");
                all_passed = false;
            }

            num_passed = 0;
            // check that all neighbors in point cloud search list are in brute force list
            Kokkos::parallel_reduce("original in brute force search", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>
                    (0,number_target_coords), KOKKOS_LAMBDA(const int i, int& t_num_passed) {
                bool all_found = true;
                for (int j=0; j<nla.getNumberOfNeighborsHost(i); ++j) {
                    if (brute_force_neighbor_list[i].count(nla.getNeighborHost(i,j))!=1) {
                        all_found = false;
                    }
                }
                if (all_found) t_num_passed++;
            }, Kokkos::Sum<int>(num_passed));
            Kokkos::fence();
            if (num_passed != number_target_coords) {
                printf("Point cloud neighbor not found in brute force list\n");
                all_passed = false;
            }

            double check_time = timer.seconds();
            printf("check time: %f s\n", check_time);
        }
        timer.reset();
        { // 2D neighbor lists
            printf("\n2D neighbor lists:\n");
            // Point cloud construction for neighbor search
            // CreatePointCloudSearch constructs an object of type PointCloudSearch, but deduces the templates for you
            auto point_cloud_search(CreatePointCloudSearch(source_coords, dimension));

            Kokkos::View<int**, Kokkos::DefaultHostExecutionSpace> neighbor_lists("neighbor lists", 
                    number_target_coords, 2); // first column is # of neighbors
            
            // each target site has a window size
            Kokkos::View<double*, Kokkos::DefaultHostExecutionSpace> epsilon("h supports", number_target_coords);
            Kokkos::View<double*, Kokkos::DefaultHostExecutionSpace> random_vec("random vector", number_target_coords);

            Kokkos::parallel_for("random perturbation of search size", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>
                    (0,number_target_coords), [&](const int i) {
                // value in [-.25, .25]
                random_vec(i) = 0.5 * ((double)rand() / (double) RAND_MAX) - 0.5;
                epsilon(i) = (h_multiplier + random_vec(i)) * h_spacing;
            });
            Kokkos::fence();
            
            // query the point cloud to generate the neighbor lists using a radius search

            // start with a dry-run, but without enough room to store the results
            auto max_num_neighbors = point_cloud_search.generate2DNeighborListsFromRadiusSearch(true /* dry run */,
                target_coords, neighbor_lists, epsilon);
            printf("max num neighbors: %lu\n", max_num_neighbors);

            // resize neighbor lists to be large enough to hold the results
            neighbor_lists = Kokkos::View<int**, Kokkos::DefaultHostExecutionSpace>("neighbor lists", 
                number_target_coords, 1+max_num_neighbors); // first column is # of neighbors

            // search again, now that we know that there is enough room to store the results
            max_num_neighbors = point_cloud_search.generate2DNeighborListsFromRadiusSearch(false /* dry run */,
                target_coords, neighbor_lists, epsilon);

            double radius_search_time = timer.seconds();
            printf("nanoflann search time: %f s\n", radius_search_time);

            // convert point cloud search to vector of maps
            timer.reset();
            std::vector<std::map<int, double> > point_cloud_neighbor_list(number_target_coords, std::map<int, double>());
            Kokkos::parallel_for("point search conversion", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>
                    (0,number_target_coords), [&](const int i) {
                auto &point_cloud_neighbor_list_i = const_cast<std::map<int, double>& >(point_cloud_neighbor_list[i]);
                for (int j=1; j<=neighbor_lists(i,0); ++j) {
                    point_cloud_neighbor_list_i[neighbor_lists(i,j)] = 1.0;
                }
            });
            double point_cloud_convert_time = timer.seconds();
            printf("point cloud convert time: %f s\n", point_cloud_convert_time);
            Kokkos::fence();

            // loop over targets, finding all of their neighbors by brute force
            timer.reset();
            std::vector<std::map<int, double> > brute_force_neighbor_list(number_target_coords, std::map<int, double>());
            Kokkos::parallel_for("brute force radius search", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>
                    (0,number_target_coords), [&](const int i) {
                auto &brute_force_neighbor_list_i = const_cast<std::map<int, double>& >(brute_force_neighbor_list[i]);
                for (int j=0; j<number_source_coords; ++j) {
                    double dist = 0;
                    for (int k=0; k<dimension; ++k) {
                        dist += (target_coords(i,k)-source_coords(j,k))*(target_coords(i,k)-source_coords(j,k));
                    }
                    dist = std::sqrt(dist);

                    if (dist<(h_multiplier + random_vec(i))*h_spacing) {
                        brute_force_neighbor_list_i[j]=dist;
                    }
                }
            });
            double brute_force_search_time = timer.seconds();
            printf("brute force search time: %f s\n", brute_force_search_time);
            Kokkos::fence();

            timer.reset();
            // check that all neighbors in brute force list are in point cloud search list
            int num_passed = 0;
            Kokkos::parallel_reduce("check brute force search in point cloud search", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>
                    (0,number_target_coords), [&](const int i, int& t_num_passed) {
                bool all_found = true;
                for (auto it=brute_force_neighbor_list[i].begin(); it!=brute_force_neighbor_list[i].end(); ++it) {
                    if (point_cloud_neighbor_list[i].count(it->first)!=1) {
                        all_found = false;
                    }
                }
                if (all_found) t_num_passed++;
            }, Kokkos::Sum<int>(num_passed));
            Kokkos::fence();
            if (num_passed != number_target_coords) {
                printf("Brute force neighbor not found in point cloud list\n");
                all_passed = false;
            }

            num_passed = 0;
            // check that all neighbors in point cloud search list are in brute force list
            Kokkos::parallel_reduce("original in brute force search", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>
                    (0,number_target_coords), KOKKOS_LAMBDA(const int i, int& t_num_passed) {
                bool all_found = true;
                for (int j=1; j<=neighbor_lists(i,0); ++j) {
                    if (brute_force_neighbor_list[i].count(neighbor_lists(i,j))!=1) {
                        all_found = false;
                    }
                }
                if (all_found) t_num_passed++;
            }, Kokkos::Sum<int>(num_passed));
            Kokkos::fence();
            if (num_passed != number_target_coords) {
                printf("Point cloud neighbor not found in brute force list\n");
                all_passed = false;
            }

            double check_time = timer.seconds();
            printf("check time: %f s\n", check_time);
        }
        timer.reset();
        { // convert 2D neighbor lists to 1D compressed row neighbor lists
            printf("\n2D neighbor lists converted to 1D compressed row neighbor lists:\n");
            // Point cloud construction for neighbor search
            // CreatePointCloudSearch constructs an object of type PointCloudSearch, but deduces the templates for you
            auto point_cloud_search(CreatePointCloudSearch(source_coords, dimension));

            Kokkos::View<int**, Kokkos::DefaultHostExecutionSpace> neighbor_lists("neighbor lists", 
                    number_target_coords, 2); // first column is # of neighbors
            
            // each target site has a window size
            Kokkos::View<double*, Kokkos::DefaultHostExecutionSpace> epsilon("h supports", number_target_coords);
            Kokkos::View<double*, Kokkos::DefaultHostExecutionSpace> random_vec("random vector", number_target_coords);

            Kokkos::parallel_for("random perturbation of search size", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>
                    (0,number_target_coords), [&](const int i) {
                // value in [-.25, .25]
                random_vec(i) = 0.5 * ((double)rand() / (double) RAND_MAX) - 0.5;
                epsilon(i) = (h_multiplier + random_vec(i)) * h_spacing;
            });
            Kokkos::fence();
            
            // query the point cloud to generate the neighbor lists using a radius search

            // start with a dry-run, but without enough room to store the results
            auto max_num_neighbors = point_cloud_search.generate2DNeighborListsFromRadiusSearch(true /* dry run */,
                target_coords, neighbor_lists, epsilon);
            printf("max num neighbors: %lu\n", max_num_neighbors);

            // resize neighbor lists to be large enough to hold the results
            neighbor_lists = Kokkos::View<int**, Kokkos::DefaultHostExecutionSpace>("neighbor lists", 
                number_target_coords, 1+max_num_neighbors); // first column is # of neighbors

            // search again, now that we know that there is enough room to store the results
            max_num_neighbors = point_cloud_search.generate2DNeighborListsFromRadiusSearch(false /* dry run */,
                target_coords, neighbor_lists, epsilon);

            auto nla = Convert2DToCompressedRowNeighborLists(neighbor_lists);

            double radius_search_time = timer.seconds();
            printf("nanoflann search time: %f s\n", radius_search_time);

            // convert point cloud search to vector of maps
            timer.reset();
            std::vector<std::map<int, double> > point_cloud_neighbor_list(number_target_coords, std::map<int, double>());
            Kokkos::parallel_for("point search conversion", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>
                    (0,number_target_coords), [&](const int i) {
                auto &point_cloud_neighbor_list_i = const_cast<std::map<int, double>& >(point_cloud_neighbor_list[i]);
                for (int j=0; j<nla.getNumberOfNeighborsHost(i); ++j) {
                    point_cloud_neighbor_list_i[nla.getNeighborHost(i,j)] = 1.0;
                }
            });
            double point_cloud_convert_time = timer.seconds();
            printf("point cloud convert time: %f s\n", point_cloud_convert_time);
            Kokkos::fence();

            // loop over targets, finding all of their neighbors by brute force
            timer.reset();
            std::vector<std::map<int, double> > brute_force_neighbor_list(number_target_coords, std::map<int, double>());
            Kokkos::parallel_for("brute force radius search", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>
                    (0,number_target_coords), [&](const int i) {
                auto &brute_force_neighbor_list_i = const_cast<std::map<int, double>& >(brute_force_neighbor_list[i]);
                for (int j=0; j<number_source_coords; ++j) {
                    double dist = 0;
                    for (int k=0; k<dimension; ++k) {
                        dist += (target_coords(i,k)-source_coords(j,k))*(target_coords(i,k)-source_coords(j,k));
                    }
                    dist = std::sqrt(dist);

                    if (dist<(h_multiplier + random_vec(i))*h_spacing) {
                        brute_force_neighbor_list_i[j]=dist;
                    }
                }
            });
            double brute_force_search_time = timer.seconds();
            printf("brute force search time: %f s\n", brute_force_search_time);
            Kokkos::fence();

            timer.reset();
            // check that all neighbors in brute force list are in point cloud search list
            int num_passed = 0;
            Kokkos::parallel_reduce("check brute force search in point cloud search", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>
                    (0,number_target_coords), [&](const int i, int& t_num_passed) {
                bool all_found = true;
                for (auto it=brute_force_neighbor_list[i].begin(); it!=brute_force_neighbor_list[i].end(); ++it) {
                    if (point_cloud_neighbor_list[i].count(it->first)!=1) {
                        all_found = false;
                    }
                }
                if (all_found) t_num_passed++;
            }, Kokkos::Sum<int>(num_passed));
            Kokkos::fence();
            if (num_passed != number_target_coords) {
                printf("Brute force neighbor not found in point cloud list\n");
                all_passed = false;
            }

            num_passed = 0;
            // check that all neighbors in point cloud search list are in brute force list
            Kokkos::parallel_reduce("original in brute force search", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>
                    (0,number_target_coords), KOKKOS_LAMBDA(const int i, int& t_num_passed) {
                bool all_found = true;
                for (int j=0; j<nla.getNumberOfNeighborsHost(i); ++j) {
                    if (brute_force_neighbor_list[i].count(nla.getNeighborHost(i,j))!=1) {
                        all_found = false;
                    }
                }
                if (all_found) t_num_passed++;
            }, Kokkos::Sum<int>(num_passed));
            Kokkos::fence();
            if (num_passed != number_target_coords) {
                printf("Point cloud neighbor not found in brute force list\n");
                all_passed = false;
            }

            double check_time = timer.seconds();
            printf("check time: %f s\n", check_time);
        }
    }
    else if (do_knn_search) {

        // do knn search
        timer.reset();
        { // 1D compressed row neighbor lists
            printf("\n1D compressed row neighbor lists:\n");
            // Point cloud construction for neighbor search
            // CreatePointCloudSearch constructs an object of type PointCloudSearch, but deduces the templates for you
            auto point_cloud_search(CreatePointCloudSearch(source_coords, dimension));

            const int min_neighbors = Compadre::GMLS::getNP(3, dimension);
            Kokkos::View<int*, Kokkos::DefaultHostExecutionSpace> neighbor_lists("neighbor lists", 
                    number_target_coords); // first column is # of neighbors
            Kokkos::View<int*, Kokkos::DefaultHostExecutionSpace> number_neighbors_list("number of neighbors list", 
                    number_target_coords); // first column is # of neighbors
            
            // each target site has a window size
            Kokkos::View<double*, Kokkos::DefaultHostExecutionSpace> epsilon("h supports", number_target_coords);
            
            // query the point cloud to generate the neighbor lists using a KNN search

            // start with a dry-run, but without enough room to store the results
            auto total_num_neighbors = point_cloud_search.generateCRNeighborListsFromKNNSearch(true /* dry-run for sizes */,
                    target_coords, neighbor_lists, number_neighbors_list, epsilon, min_neighbors, 1.5 /* cutoff_multiplier */);
            printf("total num neighbors: %lu\n", total_num_neighbors);

            // resize with room to store results
            Kokkos::resize(neighbor_lists, total_num_neighbors);

            // real knn search with space to store
            point_cloud_search.generateCRNeighborListsFromKNNSearch(false /*not dry run*/, 
                    target_coords, neighbor_lists, number_neighbors_list, epsilon, min_neighbors, 1.5 /* cutoff_multiplier */);

            auto nla(CreateNeighborLists(neighbor_lists, number_neighbors_list));

            // convert point cloud search to vector of maps
            timer.reset();
            std::vector<std::map<int, double> > point_cloud_neighbor_list(number_target_coords, std::map<int, double>());
            Kokkos::parallel_for("point search conversion", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>
                    (0,number_target_coords), [&](const int i) {
                auto &point_cloud_neighbor_list_i = const_cast<std::map<int, double>& >(point_cloud_neighbor_list[i]);
                for (int j=0; j<nla.getNumberOfNeighborsHost(i); ++j) {
                    double dist = 0;
                    for (int k=0; k<dimension; ++k) {
                        dist += (target_coords(i,k)-source_coords(nla.getNeighborHost(i,j),k))*(target_coords(i,k)-source_coords(nla.getNeighborHost(i,j),k));
                    }
                    dist = std::sqrt(dist);
                    point_cloud_neighbor_list_i[nla.getNeighborHost(i,j)] = dist;
                }
            });
            double point_cloud_convert_time = timer.seconds();
            printf("point cloud convert time: %f s\n", point_cloud_convert_time);
            Kokkos::fence();

            double eps = 1e-14;
            // need to do a sort by dist, then verify epsilon(i) = 1.5 * knn_distance
            int epsilon_diff_from_knn_dist_multiplied = 0;
            Kokkos::parallel_reduce("check kth neighbor in point cloud", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>
                    (0,number_target_coords), [&](const int i, int& t_epsilon_diff_from_knn_dist_multiplied) {
                auto &point_cloud_neighbor_list_i = const_cast<std::map<int, double>& >(point_cloud_neighbor_list[i]);
                std::multimap<double, int> inverted_map = invert_map(point_cloud_neighbor_list_i);
                int j = 0;
                for (auto it=inverted_map.begin(); it!=inverted_map.end(); ++it) {
                    if (j==min_neighbors-1) {
                        if (std::abs(epsilon(i)/1.5 - it->first) > eps) {
                            t_epsilon_diff_from_knn_dist_multiplied++;
                        }
                        break;
                    }
                    j++;
                }
            }, Kokkos::Sum<int>(epsilon_diff_from_knn_dist_multiplied));
            if (epsilon_diff_from_knn_dist_multiplied > 0) {
                printf("1.5*kth_neighbor_distance differs from calculation in epsilon vector.\n");
                all_passed = false;
            }

            // loop over targets, finding knn by brute force
            // point cloud multiplied kth by 1.5, but we are searching by brute force to ensure that for the distance h,
            // we find at least as many neighbors as min_neighbors
            // (not for the 1.5x as many, but rather for 1.0x as many)
            timer.reset();
            int total_insufficient_neighbors = 0;
            std::vector<std::map<int, double> > brute_force_neighbor_list(number_target_coords, std::map<int, double>());
            Kokkos::parallel_reduce("brute force knn search", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>
                    (0,number_target_coords), [&](const int i, int& t_total_insufficient_neighbors) {
                auto &brute_force_neighbor_list_i = const_cast<std::map<int, double>& >(brute_force_neighbor_list[i]);
                for (int j=0; j<number_source_coords; ++j) {
                    double dist = 0;
                    for (int k=0; k<dimension; ++k) {
                        dist += (target_coords(i,k)-source_coords(j,k))*(target_coords(i,k)-source_coords(j,k));
                    }
                    dist = std::sqrt(dist);

                    // epsilon was set to 1.5 * distance to kth neighbor, but we need distance to kth neighbor
                    if (dist<(epsilon(i)/1.5 + eps)) {
                        brute_force_neighbor_list_i[j]=dist;
                    }
                }
                // verify k neighbors found
                if (brute_force_neighbor_list_i.size() < (size_t)min_neighbors) t_total_insufficient_neighbors++;
            }, Kokkos::Sum<int>(total_insufficient_neighbors));
            if (total_insufficient_neighbors > 0) {
                printf("less neighbors found using kth_neighbor_distance+eps than using knn search.\n");
                all_passed = false;
            }
            double brute_force_search_time = timer.seconds();
            printf("brute force search time: %f s\n", brute_force_search_time);
            Kokkos::fence();
        }
        timer.reset();
        { // 2D neighbor lists
            printf("\n2D neighbor lists:\n");
            // Point cloud construction for neighbor search
            // CreatePointCloudSearch constructs an object of type PointCloudSearch, but deduces the templates for you
            auto point_cloud_search(CreatePointCloudSearch(source_coords, dimension));

            const int min_neighbors = Compadre::GMLS::getNP(3, dimension);
            Kokkos::View<int**, Kokkos::DefaultHostExecutionSpace> neighbor_lists("neighbor lists", 
                    number_target_coords, 1+min_neighbors); // first column is # of neighbors
            
            // each target site has a window size
            Kokkos::View<double*, Kokkos::DefaultHostExecutionSpace> epsilon("h supports", number_target_coords);
            
            // query the point cloud to generate the neighbor lists using a KNN search
            // start with a dry-run, but without enough room to store the results
            auto max_num_neighbors = point_cloud_search.generate2DNeighborListsFromKNNSearch(true /* dry-run for sizes */,
                    target_coords, neighbor_lists, epsilon, min_neighbors, 1.5 /* cutoff_multiplier */);
            printf("max num neighbors: %lu\n", max_num_neighbors);

            // resize with room to store results
            Kokkos::resize(neighbor_lists, neighbor_lists.extent(0), 1+max_num_neighbors);

            // real knn search with space to store
            max_num_neighbors = point_cloud_search.generate2DNeighborListsFromKNNSearch(false /*not dry run*/, 
                    target_coords, neighbor_lists, epsilon, min_neighbors, 1.5 /* cutoff_multiplier */);

            // convert point cloud search to vector of maps
            timer.reset();
            std::vector<std::map<int, double> > point_cloud_neighbor_list(number_target_coords, std::map<int, double>());
            Kokkos::parallel_for("point search conversion", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>
                    (0,number_target_coords), [&](const int i) {
                auto &point_cloud_neighbor_list_i = const_cast<std::map<int, double>& >(point_cloud_neighbor_list[i]);
                for (int j=1; j<=neighbor_lists(i,0); ++j) {
                    double dist = 0;
                    for (int k=0; k<dimension; ++k) {
                        dist += (target_coords(i,k)-source_coords(neighbor_lists(i,j),k))*(target_coords(i,k)-source_coords(neighbor_lists(i,j),k));
                    }
                    dist = std::sqrt(dist);
                    point_cloud_neighbor_list_i[neighbor_lists(i,j)] = dist;
                }
            });
            double point_cloud_convert_time = timer.seconds();
            printf("point cloud convert time: %f s\n", point_cloud_convert_time);
            Kokkos::fence();

            double eps = 1e-14;
            // need to do a sort by dist, then verify epsilon(i) = 1.5 * knn_distance
            int epsilon_diff_from_knn_dist_multiplied = 0;
            Kokkos::parallel_reduce("check kth neighbor in point cloud", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>
                    (0,number_target_coords), [&](const int i, int& t_epsilon_diff_from_knn_dist_multiplied) {
                auto &point_cloud_neighbor_list_i = const_cast<std::map<int, double>& >(point_cloud_neighbor_list[i]);
                std::multimap<double, int> inverted_map = invert_map(point_cloud_neighbor_list_i);
                int j = 0;
                for (auto it=inverted_map.begin(); it!=inverted_map.end(); ++it) {
                    if (j==min_neighbors-1) {
                        if (std::abs(epsilon(i)/1.5 - it->first) > eps) {
                            t_epsilon_diff_from_knn_dist_multiplied++;
                        }
                        break;
                    }
                    j++;
                }
            }, Kokkos::Sum<int>(epsilon_diff_from_knn_dist_multiplied));
            if (epsilon_diff_from_knn_dist_multiplied > 0) {
                printf("1.5*kth_neighbor_distance differs from calculation in epsilon vector.\n");
                all_passed = false;
            }

            // loop over targets, finding knn by brute force
            // point cloud multiplied kth by 1.5, but we are searching by brute force to ensure that for the distance h,
            // we find at least as many neighbors as min_neighbors
            // (not for the 1.5x as many, but rather for 1.0x as many)
            timer.reset();
            int total_insufficient_neighbors = 0;
            std::vector<std::map<int, double> > brute_force_neighbor_list(number_target_coords, std::map<int, double>());
            Kokkos::parallel_reduce("brute force knn search", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>
                    (0,number_target_coords), [&](const int i, int& t_total_insufficient_neighbors) {
                auto &brute_force_neighbor_list_i = const_cast<std::map<int, double>& >(brute_force_neighbor_list[i]);
                for (int j=0; j<number_source_coords; ++j) {
                    double dist = 0;
                    for (int k=0; k<dimension; ++k) {
                        dist += (target_coords(i,k)-source_coords(j,k))*(target_coords(i,k)-source_coords(j,k));
                    }
                    dist = std::sqrt(dist);

                    // epsilon was set to 1.5 * distance to kth neighbor, but we need distance to kth neighbor
                    if (dist<(epsilon(i)/1.5 + eps)) {
                        brute_force_neighbor_list_i[j]=dist;
                    }
                }
                // verify k neighbors found
                if (brute_force_neighbor_list_i.size() < (size_t)min_neighbors) t_total_insufficient_neighbors++;
            }, Kokkos::Sum<int>(total_insufficient_neighbors));
            if (total_insufficient_neighbors > 0) {
                printf("less neighbors found using kth_neighbor_distance+eps than using knn search.\n");
                all_passed = false;
            }
            double brute_force_search_time = timer.seconds();
            printf("brute force search time: %f s\n", brute_force_search_time);
            Kokkos::fence();
        }
        timer.reset();
        { // convert 2D neighbor lists to 1D compressed row neighbor lists
            printf("\n2D neighbor lists converted to 1D compressed row neighbor lists:\n");
            // Point cloud construction for neighbor search
            // CreatePointCloudSearch constructs an object of type PointCloudSearch, but deduces the templates for you
            auto point_cloud_search(CreatePointCloudSearch(source_coords, dimension));

            const int min_neighbors = Compadre::GMLS::getNP(3, dimension);
            Kokkos::View<int**, Kokkos::DefaultHostExecutionSpace> neighbor_lists("neighbor lists", 
                    number_target_coords, 1+min_neighbors); // first column is # of neighbors
            
            // each target site has a window size
            Kokkos::View<double*, Kokkos::DefaultHostExecutionSpace> epsilon("h supports", number_target_coords);
            
            // query the point cloud to generate the neighbor lists using a KNN search
            // start with a dry-run, but without enough room to store the results
            auto max_num_neighbors = point_cloud_search.generate2DNeighborListsFromKNNSearch(true /* dry-run for sizes */,
                    target_coords, neighbor_lists, epsilon, min_neighbors, 1.5 /* cutoff_multiplier */);
            printf("max num neighbors: %lu\n", max_num_neighbors);

            // resize with room to store results
            Kokkos::resize(neighbor_lists, neighbor_lists.extent(0), 1+max_num_neighbors);

            // real knn search with space to store
            max_num_neighbors = point_cloud_search.generate2DNeighborListsFromKNNSearch(false /*not dry run*/, 
                    target_coords, neighbor_lists, epsilon, min_neighbors, 1.5 /* cutoff_multiplier */);

            auto nla = Convert2DToCompressedRowNeighborLists(neighbor_lists);

            // convert point cloud search to vector of maps
            timer.reset();
            std::vector<std::map<int, double> > point_cloud_neighbor_list(number_target_coords, std::map<int, double>());
            Kokkos::parallel_for("point search conversion", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>
                    (0,number_target_coords), [&](const int i) {
                auto &point_cloud_neighbor_list_i = const_cast<std::map<int, double>& >(point_cloud_neighbor_list[i]);
                for (int j=0; j<nla.getNumberOfNeighborsHost(i); ++j) {
                    double dist = 0;
                    for (int k=0; k<dimension; ++k) {
                        dist += (target_coords(i,k)-source_coords(nla.getNeighborHost(i,j),k))*(target_coords(i,k)-source_coords(nla.getNeighborHost(i,j),k));
                    }
                    dist = std::sqrt(dist);
                    point_cloud_neighbor_list_i[nla.getNeighborHost(i,j)] = dist;
                }
            });
            double point_cloud_convert_time = timer.seconds();
            printf("point cloud convert time: %f s\n", point_cloud_convert_time);
            Kokkos::fence();

            double eps = 1e-14;
            // need to do a sort by dist, then verify epsilon(i) = 1.5 * knn_distance
            int epsilon_diff_from_knn_dist_multiplied = 0;
            Kokkos::parallel_reduce("check kth neighbor in point cloud", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>
                    (0,number_target_coords), [&](const int i, int& t_epsilon_diff_from_knn_dist_multiplied) {
                auto &point_cloud_neighbor_list_i = const_cast<std::map<int, double>& >(point_cloud_neighbor_list[i]);
                std::multimap<double, int> inverted_map = invert_map(point_cloud_neighbor_list_i);
                int j = 0;
                for (auto it=inverted_map.begin(); it!=inverted_map.end(); ++it) {
                    if (j==min_neighbors-1) {
                        if (std::abs(epsilon(i)/1.5 - it->first) > eps) {
                            t_epsilon_diff_from_knn_dist_multiplied++;
                        }
                        break;
                    }
                    j++;
                }
            }, Kokkos::Sum<int>(epsilon_diff_from_knn_dist_multiplied));
            if (epsilon_diff_from_knn_dist_multiplied > 0) {
                printf("1.5*kth_neighbor_distance differs from calculation in epsilon vector.\n");
                all_passed = false;
            }

            // loop over targets, finding knn by brute force
            // point cloud multiplied kth by 1.5, but we are searching by brute force to ensure that for the distance h,
            // we find at least as many neighbors as min_neighbors
            // (not for the 1.5x as many, but rather for 1.0x as many)
            timer.reset();
            int total_insufficient_neighbors = 0;
            std::vector<std::map<int, double> > brute_force_neighbor_list(number_target_coords, std::map<int, double>());
            Kokkos::parallel_reduce("brute force knn search", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>
                    (0,number_target_coords), [&](const int i, int& t_total_insufficient_neighbors) {
                auto &brute_force_neighbor_list_i = const_cast<std::map<int, double>& >(brute_force_neighbor_list[i]);
                for (int j=0; j<number_source_coords; ++j) {
                    double dist = 0;
                    for (int k=0; k<dimension; ++k) {
                        dist += (target_coords(i,k)-source_coords(j,k))*(target_coords(i,k)-source_coords(j,k));
                    }
                    dist = std::sqrt(dist);

                    // epsilon was set to 1.5 * distance to kth neighbor, but we need distance to kth neighbor
                    if (dist<(epsilon(i)/1.5 + eps)) {
                        brute_force_neighbor_list_i[j]=dist;
                    }
                }
                // verify k neighbors found
                if (brute_force_neighbor_list_i.size() < (size_t)min_neighbors) t_total_insufficient_neighbors++;
            }, Kokkos::Sum<int>(total_insufficient_neighbors));
            if (total_insufficient_neighbors > 0) {
                printf("less neighbors found using kth_neighbor_distance+eps than using knn search.\n");
                all_passed = false;
            }
            double brute_force_search_time = timer.seconds();
            printf("brute force search time: %f s\n", brute_force_search_time);
            Kokkos::fence();
        }
    } else {
        all_passed = false;
        printf("do_radius_search and do_knn_search both false. Invalid combination.\n");
    }
}

// finalize Kokkos and MPI (if available)
Kokkos::finalize();
#ifdef COMPADRE_USE_MPI
MPI_Finalize();
#endif

// output to user that test passed or failed
if(all_passed) {
    fprintf(stdout, "Passed test \n");
    return 0;
} else {
    fprintf(stdout, "Failed test \n");
    return -1;
}

} // main
