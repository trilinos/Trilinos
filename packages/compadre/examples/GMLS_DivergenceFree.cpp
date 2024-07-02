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
#include "CommandLineProcessor.hpp"

#ifdef COMPADRE_USE_MPI
#include <mpi.h>
#endif

#include <Kokkos_Timer.hpp>
#include <Kokkos_Core.hpp>

using namespace Compadre;

//! [Parse Command Line Arguments]

// called from command line
int main (int argc, char* args[]) {

// initializes MPI (if available) with command line arguments given
#ifdef COMPADRE_USE_MPI
MPI_Init(&argc, &args);
#endif

// initializes Kokkos with command line arguments given
Kokkos::initialize(argc, args);

// becomes false if the computed solution not within the failure_threshold of the actual solution
bool all_passed = true;

// code block to reduce scope for all Kokkos View allocations
// otherwise, Views may be deallocating when we call Kokkos::finalize() later
{

    CommandLineProcessor clp(argc, args);
    auto order = clp.order;
    auto dimension = clp.dimension;
    auto number_target_coords = clp.number_target_coords;
    auto constraint_name = clp.constraint_name;
    auto solver_name = clp.solver_name;
    auto problem_name = clp.problem_name;

    // the functions we will be seeking to reconstruct are in the span of the basis
    // of the reconstruction space we choose for GMLS, so the error should be very small
    const double failure_tolerance = 1e-9;

    // minimum neighbors for unisolvency is the same as the size of the polynomial basis
    const int min_neighbors = Compadre::GMLS::getNP(order, dimension, DivergenceFreeVectorTaylorPolynomial);

    //! [Parse Command Line Arguments]
    Kokkos::Timer timer;
    Kokkos::Profiling::pushRegion("Setup Point Data");
    //! [Setting Up The Point Cloud]

    // approximate spacing of source sites
    double h_spacing = 0.05;
    int n_neg1_to_1 = 2*(1/h_spacing) + 1; // always odd

    // number of source coordinate sites that will fill a box of [-1,1]x[-1,1]x[-1,1] with a spacing approximately h
    const int number_source_coords = std::pow(n_neg1_to_1, dimension);

    // coordinates of source sites
    Kokkos::View<double**, Kokkos::DefaultExecutionSpace> source_coords_device("source coordinates",
            number_source_coords, 3);
    Kokkos::View<double**>::HostMirror source_coords = Kokkos::create_mirror_view(source_coords_device);

    // coordinates of target sites
    Kokkos::View<double**, Kokkos::DefaultExecutionSpace> target_coords_device ("target coordinates", number_target_coords, 3);
    Kokkos::View<double**>::HostMirror target_coords = Kokkos::create_mirror_view(target_coords_device);


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

    // Generate target points - these are random permutation from available source points
    // Note that this is assuming that the number of targets in this test will not exceed
    // the number of source coords, which is 41^3 = 68921
    // seed random number generator
    std::mt19937 rng(50);
    // generate random integers in [0..number_source_coords-1] (used to pick target sites)
    std::uniform_int_distribution<int> gen_num_neighbours(0, source_index);
    // fill target sites with random selections from source sites
    for (int i=0; i<number_target_coords; ++i) {
        const int source_site_to_copy = gen_num_neighbours(rng);
        for (int j=0; j<3; j++) {
            target_coords(i, j) = source_coords(source_site_to_copy, j);
        }
    }

    //! [Setting Up The Point Cloud]

    Kokkos::Profiling::popRegion();
    Kokkos::Profiling::pushRegion("Creating Data");

    //! [Creating The Data]


    // source coordinates need copied to device before using to construct sampling data
    Kokkos::deep_copy(source_coords_device, source_coords);

    // target coordinates copied next, because it is a convenient time to send them to device
    Kokkos::deep_copy(target_coords_device, target_coords);

    // need Kokkos View storing true solution
    Kokkos::View<double**, Kokkos::DefaultExecutionSpace> vector_sampling_span_basis_data_device("samples of true vector solution",
            source_coords_device.extent(0), dimension);
    Kokkos::View<double**, Kokkos::DefaultExecutionSpace> vector_sampling_single_polynomial_data_device("samples of true vector solution",
            source_coords_device.extent(0), dimension);

    Kokkos::parallel_for("Sampling Manufactured Solutions", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>
            (0,source_coords.extent(0)), KOKKOS_LAMBDA(const int i) {
        // coordinates of source site i
        double xval = source_coords_device(i,0);
        double yval = (dimension>1) ? source_coords_device(i,1) : 0;
        double zval = (dimension>2) ? source_coords_device(i,2) : 0;

        // data for targets with vector input
        for (int j=0; j<dimension; ++j) {
            vector_sampling_span_basis_data_device(i, j) = divfreeTestSolution_span_basis(xval, yval, zval, j, dimension, order);
            vector_sampling_single_polynomial_data_device(i, j) = divfreeTestSolution_single_polynomial(xval, yval, zval, j, dimension);
        }
    });


    //! [Creating The Data]

    Kokkos::Profiling::popRegion();
    Kokkos::Profiling::pushRegion("Neighbor Search");

    //! [Performing Neighbor Search]


    // Point cloud construction for neighbor search
    // CreatePointCloudSearch constructs an object of type PointCloudSearch, but deduces the templates for you
    auto point_cloud_search(CreatePointCloudSearch(source_coords, dimension));

    // each row is a neighbor list for a target site, with the first column of each row containing
    // the number of neighbors for that rows corresponding target site
    // for the default values in this test, the multiplier is suggested to be 2.2
    double epsilon_multiplier = 2.2;
    int estimated_upper_bound_number_neighbors =
        point_cloud_search.getEstimatedNumberNeighborsUpperBound(min_neighbors, dimension, epsilon_multiplier);

    Kokkos::View<int**, Kokkos::DefaultExecutionSpace> neighbor_lists_device("neighbor lists",
            number_target_coords, estimated_upper_bound_number_neighbors); // first column is # of neighbors
    Kokkos::View<int**>::HostMirror neighbor_lists = Kokkos::create_mirror_view(neighbor_lists_device);

    // each target site has a window size
    Kokkos::View<double*, Kokkos::DefaultExecutionSpace> epsilon_device("h supports", number_target_coords);
    Kokkos::View<double*>::HostMirror epsilon = Kokkos::create_mirror_view(epsilon_device);

    // query the point cloud to generate the neighbor lists using a kdtree to produce the n nearest neighbor
    // to each target site, adding (epsilon_multiplier-1)*100% to whatever the distance away the further neighbor used is from
    // each target to the view for epsilon
    point_cloud_search.generate2DNeighborListsFromKNNSearch(false /*not dry run*/, target_coords, neighbor_lists, 
            epsilon, min_neighbors, epsilon_multiplier);

    //! [Performing Neighbor Search]

    Kokkos::Profiling::popRegion();
    Kokkos::fence(); // let call to build neighbor lists complete before copying back to device
    timer.reset();

    //! [Setting Up The GMLS Object]


    // Copy data back to device (they were filled on the host)
    // We could have filled Kokkos Views with memory space on the host
    // and used these instead, and then the copying of data to the device
    // would be performed in the GMLS class
    Kokkos::deep_copy(neighbor_lists_device, neighbor_lists);
    Kokkos::deep_copy(epsilon_device, epsilon);

    // initialize an instance of the GMLS class
    // NULL in manifold order indicates non-manifold case
    // Vector basis to perform GMLS on divergence free basis
    GMLS vector_divfree_basis_gmls(DivergenceFreeVectorTaylorPolynomial,
                                   VectorPointSample,
                                   order, dimension,
                                   solver_name.c_str(), problem_name.c_str(), constraint_name.c_str(),
                                   0 /*manifold order*/);

    // pass in neighbor lists, source coordinates, target coordinates, and window sizes
    //
    // neighbor lists have the format:
    //      dimensions: (# number of target sites) X (# maximum number of neighbors for any given target + 1)
    //                  the first column contains the number of neighbors for that rows corresponding target index
    //
    // source coordinates have the format:
    //      dimensions: (# number of source sites) X (dimension)
    //                  entries in the neighbor lists (integers) correspond to rows of this 2D array
    //
    // target coordinates have the format:
    //      dimensions: (# number of target sites) X (dimension)
    //                  # of target sites is same as # of rows of neighbor lists
    //
    vector_divfree_basis_gmls.setProblemData(neighbor_lists_device, source_coords_device, target_coords_device, epsilon_device);

    // create a vector of target operations
    std::vector<TargetOperation> lro(4);
    lro[0] = VectorPointEvaluation;
    lro[1] = CurlOfVectorPointEvaluation;
    lro[2] = CurlCurlOfVectorPointEvaluation;
    lro[3] = GradientOfVectorPointEvaluation;

    // and then pass them to the GMLS class
    vector_divfree_basis_gmls.addTargets(lro);

    // sets the weighting kernel function from WeightingFunctionType
    vector_divfree_basis_gmls.setWeightingType(WeightingFunctionType::Power);

    // power to use in that weighting kernel function
    vector_divfree_basis_gmls.setWeightingParameter(2);

    // generate the alphas that to be combined with data for each target operation requested in lro
    vector_divfree_basis_gmls.generateAlphas(15 /* # batches */);

    //! [Setting Up The GMLS Object]

    double instantiation_time = timer.seconds();
    std::cout << "Took " << instantiation_time << "s to complete alphas generation." << std::endl;
    Kokkos::fence(); // let generateAlphas finish up before using alphas
    Kokkos::Profiling::pushRegion("Apply Alphas to Data");

    //! [Apply GMLS Alphas To Data]

    // it is important to note that if you expect to use the data as a 1D view, then you should use double*
    // however, if you know that the target operation will result in a 2D view (vector or matrix output),
    // then you should template with double** as this is something that can not be infered from the input data
    // or the target operator at compile time. Additionally, a template argument is required indicating either
    // Kokkos::HostSpace or Kokkos::DefaultExecutionSpace::memory_space()

    // The Evaluator class takes care of handling input data views as well as the output data views.
    // It uses information from the GMLS class to determine how many components are in the input
    // as well as output for any choice of target functionals and then performs the contactions
    // on the data using the alpha coefficients generated by the GMLS class, all on the device.
    Evaluator gmls_evaluator_vector(&vector_divfree_basis_gmls);

    auto output_vector_evaluation = gmls_evaluator_vector.applyAlphasToDataAllComponentsAllTargetSites<double**, Kokkos::HostSpace>
      (vector_sampling_span_basis_data_device, VectorPointEvaluation, VectorPointSample);
    auto output_curl_vector_evaluation = gmls_evaluator_vector.applyAlphasToDataAllComponentsAllTargetSites<double**, Kokkos::HostSpace>
      (vector_sampling_span_basis_data_device, CurlOfVectorPointEvaluation, VectorPointSample);
    auto output_curlcurl_vector_evaluation = gmls_evaluator_vector.applyAlphasToDataAllComponentsAllTargetSites<double**, Kokkos::HostSpace>
      (vector_sampling_single_polynomial_data_device, CurlCurlOfVectorPointEvaluation, VectorPointSample);
    auto output_gradient_vector_evaluation = gmls_evaluator_vector.applyAlphasToDataAllComponentsAllTargetSites<double**, Kokkos::HostSpace>
      (vector_sampling_span_basis_data_device, GradientOfVectorPointEvaluation, VectorPointSample);

    //! [Apply GMLS Alphas To Data]

    Kokkos::fence(); // let application of alphas to data finish before using results
    Kokkos::Profiling::popRegion();
    // times the Comparison in Kokkos
    Kokkos::Profiling::pushRegion("Comparison");

    //! [Check That Solutions Are Correct]

    // loop through the target sites
    for (int i=0; i<number_target_coords; i++) {
        // load vector components from output
        double GMLS_DivFree_VectorX = output_vector_evaluation(i, 0);
        double GMLS_DivFree_VectorY = (dimension>1) ? output_vector_evaluation(i, 1) : 0;
        double GMLS_DivFree_VectorZ = (dimension>2) ? output_vector_evaluation(i, 2) : 0;

        // load curl of vector components from output
        double GMLS_Curl_DivFree_VectorX = output_curl_vector_evaluation(i, 0);
        double GMLS_Curl_DivFree_VectorY = (dimension>2) ? output_curl_vector_evaluation(i, 1) : 0;
        double GMLS_Curl_DivFree_VectorZ = (dimension>2) ? output_curl_vector_evaluation(i, 2) : 0;

        // load curl curl of vector components from output
        double GMLS_CurlCurl_DivFree_VectorX = output_curlcurl_vector_evaluation(i, 0);
        double GMLS_CurlCurl_DivFree_VectorY = (dimension>1) ? output_curlcurl_vector_evaluation(i, 1) : 0;
        double GMLS_CurlCurl_DivFree_VectorZ = (dimension>2) ? output_curlcurl_vector_evaluation(i, 2) : 0;

        // load gradient of vector components from output
        double GMLS_Gradient_DivFree_VectorXX = output_gradient_vector_evaluation(i, 0);
        double GMLS_Gradient_DivFree_VectorXY = output_gradient_vector_evaluation(i, 1);
        double GMLS_Gradient_DivFree_VectorXZ = (dimension==3) ? output_gradient_vector_evaluation(i, 2) : 0.0;
        double GMLS_Gradient_DivFree_VectorYX = (dimension==3) ? output_gradient_vector_evaluation(i, 3) : output_gradient_vector_evaluation(i, 2);
        double GMLS_Gradient_DivFree_VectorYY = (dimension==3) ? output_gradient_vector_evaluation(i, 4) : output_gradient_vector_evaluation(i, 3);
        double GMLS_Gradient_DivFree_VectorYZ = (dimension==3) ? output_gradient_vector_evaluation(i, 5) : 0.0;
        double GMLS_Gradient_DivFree_VectorZX = (dimension==3) ? output_gradient_vector_evaluation(i, 6) : 0.0;
        double GMLS_Gradient_DivFree_VectorZY = (dimension==3) ? output_gradient_vector_evaluation(i, 7) : 0.0;
        double GMLS_Gradient_DivFree_VectorZZ = (dimension==3) ? output_gradient_vector_evaluation(i, 8) : 0.0;

        // target site i's coordinate
        double xval = target_coords(i,0);
        double yval = (dimension>1) ? target_coords(i,1) : 0;
        double zval = (dimension>2) ? target_coords(i,2) : 0;

        // evaluation of vector exact solutions
        double actual_vector[3] = {0, 0, 0};
        if (dimension>=1) {
            actual_vector[0] = divfreeTestSolution_span_basis(xval, yval, zval, 0, dimension, order);
            if (dimension>=2) {
                actual_vector[1] = divfreeTestSolution_span_basis(xval, yval, zval, 1, dimension, order);
                if (dimension==3) {
                    actual_vector[2] = divfreeTestSolution_span_basis(xval, yval, zval, 2, dimension, order);
                }
            }
        }

        // evaluation of curl of vector exact solutions
        double actual_curl_vector[3] = {0, 0, 0};
        if (dimension>=1) {
            actual_curl_vector[0] = curldivfreeTestSolution_span_basis(xval, yval, zval, 0, dimension, order);
            if (dimension==3) {
                actual_curl_vector[1] = curldivfreeTestSolution_span_basis(xval, yval, zval, 1, dimension, order);
                actual_curl_vector[2] = curldivfreeTestSolution_span_basis(xval, yval, zval, 2, dimension, order);
            }
        }

        // evaluation of curl of curl of vector exact solutions
        double actual_curlcurl_vector[3] = {0, 0, 0};
        if (dimension>=1) {
            actual_curlcurl_vector[0] = curlcurldivfreeTestSolution_single_polynomial(xval, yval, zval, 0, dimension);
            if (dimension>=2) {
                actual_curlcurl_vector[1] = curlcurldivfreeTestSolution_single_polynomial(xval, yval, zval, 1, dimension);
                if (dimension==3) {
                    actual_curlcurl_vector[2] = curlcurldivfreeTestSolution_single_polynomial(xval, yval, zval, 2, dimension);
                }
            }
        }

        double actual_gradient_vector[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        if (dimension==3) {
            for (int axes = 0; axes < 9; ++axes)
                actual_gradient_vector[axes] = gradientdivfreeTestSolution_span_basis(xval, yval, zval, axes/dimension, axes%dimension, dimension, order);
        }
        if (dimension==2) {
            for (int axes = 0; axes < 4; ++axes) 
                actual_gradient_vector[axes] = gradientdivfreeTestSolution_span_basis(xval, yval, zval, axes/dimension, axes%dimension, dimension, order);
        }

        // check vector evaluation
        if(std::abs(actual_vector[0] - GMLS_DivFree_VectorX) > failure_tolerance) {
            all_passed = false;
            std::cout << i << " Failed VectorX by: " << std::abs(actual_vector[0] - GMLS_DivFree_VectorX) << std::endl;
            if (dimension>1) {
                if(std::abs(actual_vector[1] - GMLS_DivFree_VectorY) > failure_tolerance) {
                    all_passed = false;
                    std::cout << i << " Failed VectorY by: " << std::abs(actual_vector[1] - GMLS_DivFree_VectorY) << std::endl;
                }
            }
            if (dimension>2) {
                if(std::abs(actual_vector[2] - GMLS_DivFree_VectorZ) > failure_tolerance) {
                    all_passed = false;
                    std::cout << i << " Failed VectorZ by: " << std::abs(actual_vector[2] - GMLS_DivFree_VectorZ) << std::endl;
                }
            }
        }

        // check curl of vector evaluation
        if (dimension==2) {
            if(std::abs(actual_curl_vector[0] - GMLS_Curl_DivFree_VectorX) > failure_tolerance) {
                all_passed = false;
                std::cout << i << " Failed curl VectorX by: " << std::abs(actual_curl_vector[2] - GMLS_Curl_DivFree_VectorX) << std::endl;
            }
        } else if (dimension>2) {
            if(std::abs(actual_curl_vector[0] - GMLS_Curl_DivFree_VectorX) > failure_tolerance) {
                all_passed = false;
                std::cout << i << " Failed curl VectorX by: " << std::abs(actual_curl_vector[0] - GMLS_Curl_DivFree_VectorX) << std::endl;
            }
            if(std::abs(actual_curl_vector[1] - GMLS_Curl_DivFree_VectorY) > failure_tolerance) {
                all_passed = false;
                std::cout << i << " Failed curl VectorY by: " << std::abs(actual_curl_vector[1] - GMLS_Curl_DivFree_VectorY) << std::endl;
            }
            if(std::abs(actual_curl_vector[2] - GMLS_Curl_DivFree_VectorZ) > failure_tolerance) {
                all_passed = false;
                std::cout << i << " Failed curl VectorZ by: " << std::abs(actual_curl_vector[2] - GMLS_Curl_DivFree_VectorZ) << std::endl;
            }
        }

        // check curlcurl curlcurl of vector evaluation
        if(std::abs(actual_curlcurl_vector[0] - GMLS_CurlCurl_DivFree_VectorX) > failure_tolerance) {
            all_passed = false;
            std::cout << i << " Failed curl curl VectorX by: " << std::abs(actual_curlcurl_vector[0] - GMLS_CurlCurl_DivFree_VectorX) << std::endl;
        }
        if (dimension>1) {
            if(std::abs(actual_curlcurl_vector[1] - GMLS_CurlCurl_DivFree_VectorY) > failure_tolerance) {
                all_passed = false;
                std::cout << i << " Failed curl curl VectorY by: " << std::abs(actual_curlcurl_vector[1] - GMLS_CurlCurl_DivFree_VectorY) << std::endl;
            }
        }
        if (dimension>2) {
            if(std::abs(actual_curlcurl_vector[2] - GMLS_CurlCurl_DivFree_VectorZ) > failure_tolerance) {
                all_passed = false;
                std::cout << i << " Failed curl curl VectorZ by: " << std::abs(actual_curlcurl_vector[2] - GMLS_CurlCurl_DivFree_VectorZ) << std::endl;
            }
        }

        // check gradient of vector evaluation
        if (dimension==3) {
            if (std::abs(actual_gradient_vector[0] - GMLS_Gradient_DivFree_VectorXX) > failure_tolerance) {
                all_passed = false;
                std::cout << i << " Failed gradient_x VectorX by: " << std::abs(actual_gradient_vector[0] - GMLS_Gradient_DivFree_VectorXX) << std::endl;
            }
            if (std::abs(actual_gradient_vector[1] - GMLS_Gradient_DivFree_VectorXY) > failure_tolerance) {
                all_passed = false;
                std::cout << i << " Failed gradient_y VectorX by: " << std::abs(actual_gradient_vector[1] - GMLS_Gradient_DivFree_VectorXY) << std::endl;
            }
            if (std::abs(actual_gradient_vector[2] - GMLS_Gradient_DivFree_VectorXZ) > failure_tolerance) {
                all_passed = false;
                std::cout << i << " Failed gradient_z VectorX by: " << std::abs(actual_gradient_vector[2] - GMLS_Gradient_DivFree_VectorXZ) << std::endl;
            }

            if (std::abs(actual_gradient_vector[3] - GMLS_Gradient_DivFree_VectorYX) > failure_tolerance) {
                all_passed = false;
                std::cout << i << " Failed gradient_x VectorY by: " << std::abs(actual_gradient_vector[3] - GMLS_Gradient_DivFree_VectorYX) << std::endl;
            }
            if (std::abs(actual_gradient_vector[4] - GMLS_Gradient_DivFree_VectorYY) > failure_tolerance) {
                all_passed = false;
                std::cout << i << " Failed gradient_y VectorY by: " << std::abs(actual_gradient_vector[4] - GMLS_Gradient_DivFree_VectorYY) << std::endl;
            }
            if (std::abs(actual_gradient_vector[5] - GMLS_Gradient_DivFree_VectorYZ) > failure_tolerance) {
                all_passed = false;
                std::cout << i << " Failed gradient_z VectorY by: " << std::abs(actual_gradient_vector[5] - GMLS_Gradient_DivFree_VectorYZ) << std::endl;
            }

            if (std::abs(actual_gradient_vector[6] - GMLS_Gradient_DivFree_VectorZX) > failure_tolerance) {
                all_passed = false;
                std::cout << i << " Failed gradient_x VectorZ by: " << std::abs(actual_gradient_vector[6] - GMLS_Gradient_DivFree_VectorZX) << std::endl;
            }
            if (std::abs(actual_gradient_vector[7] - GMLS_Gradient_DivFree_VectorZY) > failure_tolerance) {
                all_passed = false;
                std::cout << i << " Failed gradient_y VectorZ by: " << std::abs(actual_gradient_vector[7] - GMLS_Gradient_DivFree_VectorZY) << std::endl;
            }
            if (std::abs(actual_gradient_vector[8] - GMLS_Gradient_DivFree_VectorZZ) > failure_tolerance) {
                all_passed = false;
                std::cout << i << " Failed gradient_z VectorZ by: " << std::abs(actual_gradient_vector[8] - GMLS_Gradient_DivFree_VectorZZ) << std::endl;
            }
        }

        if (dimension==2) {
            if (std::abs(actual_gradient_vector[0] - GMLS_Gradient_DivFree_VectorXX) > failure_tolerance) {
                all_passed = false;
                std::cout << i << " Failed gradient_x VectorX by: " << std::abs(actual_gradient_vector[0] - GMLS_Gradient_DivFree_VectorXX) << std::endl;
            }
            if (std::abs(actual_gradient_vector[1] - GMLS_Gradient_DivFree_VectorXY) > failure_tolerance) {
                all_passed = false;
                std::cout << i << " Failed gradient_y VectorX by: " << std::abs(actual_gradient_vector[1] - GMLS_Gradient_DivFree_VectorXY) << std::endl;
            }

            if (std::abs(actual_gradient_vector[2] - GMLS_Gradient_DivFree_VectorYX) > failure_tolerance) {
                all_passed = false;
                std::cout << i << " Failed gradient_x VectorY by: " << std::abs(actual_gradient_vector[2] - GMLS_Gradient_DivFree_VectorYX) << std::endl;
            }
            if (std::abs(actual_gradient_vector[3] - GMLS_Gradient_DivFree_VectorYY) > failure_tolerance) {
                all_passed = false;
                std::cout << i << " Failed gradient_y VectorY by: " << (std::abs(actual_gradient_vector[3] - GMLS_Gradient_DivFree_VectorYY) > failure_tolerance) << std::endl;
            }
        }
    }

    //! [Check That Solutions Are Correct]
    // popRegion hidden from tutorial
    // stop timing comparison loop
    Kokkos::Profiling::popRegion();
    //! [Finalize Program]


} // end of code block to reduce scope, causing Kokkos View de-allocations
// otherwise, Views may be deallocating when we call Kokkos::finalize() later

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


//! [Finalize Program]
