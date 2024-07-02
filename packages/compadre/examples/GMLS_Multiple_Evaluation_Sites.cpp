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
    const int min_neighbors = Compadre::GMLS::getNP(order, dimension);
    
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

    // coordinates of additional evaluation sites
    Kokkos::View<double**, Kokkos::DefaultExecutionSpace> additional_target_coords_device ("additional target coordinates", 2*number_target_coords /* multiple evaluation sites for each target index */, 3);
    Kokkos::View<double**>::HostMirror additional_target_coords = Kokkos::create_mirror_view(additional_target_coords_device);

    // additional target site indices
    Kokkos::View<int**, Kokkos::DefaultExecutionSpace> additional_target_indices_device ("additional target indices", number_target_coords, 4 /* # of extra evaluation sites plus index for each */);
    Kokkos::View<int**>::HostMirror additional_target_indices = Kokkos::create_mirror_view(additional_target_indices_device);
    
    
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

    // generate coordinates to test multiple site evaluations
    // strategy is to have a variable number of evaluation sites per target site
    // so as to fully test the multi-site evaluation 
    int extra_evaluation_coordinates_count = 0;
    for(int i=0; i<number_target_coords; i++){

        // set list of indices for extra evaluations
        additional_target_indices(i,0) = (i%3)+1;
    
        // evaluation sites are same as target plus some perturbation
        for (int k=0; k<(i%3+1); ++k) {
            for (int j=0; j<dimension; ++j) {
                additional_target_coords(extra_evaluation_coordinates_count,j) = target_coords(i,j) + (j==0)*1e-3 + (j==1)*1e-2 + (j==1)*(-1e-1);
            }
            additional_target_indices(i,k+1) = extra_evaluation_coordinates_count;
            extra_evaluation_coordinates_count++;
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
    
    // additional evaluation coordinates copied next, because it is a convenient time to send them to device
    Kokkos::deep_copy(additional_target_coords_device, additional_target_coords);

    // additional evaluation indices copied next, because it is a convenient time to send them to device
    Kokkos::deep_copy(additional_target_indices_device, additional_target_indices);

    // need Kokkos View storing true solution
    Kokkos::View<double*, Kokkos::DefaultExecutionSpace> sampling_data_device("samples of true solution", 
            source_coords_device.extent(0));
    
    Kokkos::View<double**, Kokkos::DefaultExecutionSpace> gradient_sampling_data_device("samples of true gradient", 
            source_coords_device.extent(0), dimension);
    
    Kokkos::View<double**, Kokkos::DefaultExecutionSpace> divergence_sampling_data_device
            ("samples of true solution for divergence test", source_coords_device.extent(0), dimension);
    
    Kokkos::parallel_for("Sampling Manufactured Solutions", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>
            (0,source_coords.extent(0)), KOKKOS_LAMBDA(const int i) {
    
        // coordinates of source site i
        double xval = source_coords_device(i,0);
        double yval = (dimension>1) ? source_coords_device(i,1) : 0;
        double zval = (dimension>2) ? source_coords_device(i,2) : 0;
    
        // data for targets with scalar input
        sampling_data_device(i) = trueSolution(xval, yval, zval, order, dimension);
    
        // data for targets with vector input (divergence)
        double true_grad[3] = {0,0,0};
        trueGradient(true_grad, xval, yval,zval, order, dimension);
    
        for (int j=0; j<dimension; ++j) {
            gradient_sampling_data_device(i,j) = true_grad[j];
    
            // data for target with vector input (curl)
            divergence_sampling_data_device(i,j) = divergenceTestSamples(xval, yval, zval, j, dimension);
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
    double epsilon_multiplier = 1.5;
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
    GMLS my_GMLS(VectorOfScalarClonesTaylorPolynomial, VectorPointSample,
                 order, dimension,
                 solver_name.c_str(), problem_name.c_str(), constraint_name.c_str(),
                 2 /*manifold order*/);

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
    my_GMLS.setProblemData(neighbor_lists_device, source_coords_device, target_coords_device, epsilon_device);

    // set up additional sites to evaluate target operators
    my_GMLS.setAdditionalEvaluationSitesData(additional_target_indices_device, additional_target_coords_device);
    
    // create a vector of target operations
    std::vector<TargetOperation> lro(2);
    lro[0] = ScalarPointEvaluation;
    lro[1] = GradientOfScalarPointEvaluation;
    
    // and then pass them to the GMLS class
    my_GMLS.addTargets(lro);
    
    // sets the weighting kernel function from WeightingFunctionType
    my_GMLS.setWeightingType(WeightingFunctionType::Power);
    
    // power to use in that weighting kernel function
    my_GMLS.setWeightingParameter(2);
    
    // generate the alphas that to be combined with data for each target operation requested in lro
    my_GMLS.generateAlphas();
    
    
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
    Evaluator gmls_evaluator(&my_GMLS);
    
    auto output_value1 = gmls_evaluator.applyAlphasToDataAllComponentsAllTargetSites<double*, Kokkos::HostSpace>
            (sampling_data_device, ScalarPointEvaluation, PointSample, 
             true /*scalar_as_vector_if_needed*/, 1 /*evaluation site index*/);
    
    auto output_gradient1 = gmls_evaluator.applyAlphasToDataAllComponentsAllTargetSites<double**, Kokkos::HostSpace>
            (sampling_data_device, GradientOfScalarPointEvaluation, PointSample, 
             true /*scalar_as_vector_if_needed*/, 1 /*evaluation site index*/);

    auto output_value2 = gmls_evaluator.applyAlphasToDataAllComponentsAllTargetSites<double*, Kokkos::HostSpace>
            (sampling_data_device, ScalarPointEvaluation, PointSample, 
             true /*scalar_as_vector_if_needed*/, 2 /*evaluation site index*/);
    
    auto output_gradient2 = gmls_evaluator.applyAlphasToDataAllComponentsAllTargetSites<double**, Kokkos::HostSpace>
            (sampling_data_device, GradientOfScalarPointEvaluation, PointSample, 
             true /*scalar_as_vector_if_needed*/, 2 /*evaluation site index*/);

    auto output_value3 = gmls_evaluator.applyAlphasToDataAllComponentsAllTargetSites<double*, Kokkos::HostSpace>
            (sampling_data_device, ScalarPointEvaluation, PointSample, 
             true /*scalar_as_vector_if_needed*/, 3 /*evaluation site index*/);
    
    auto output_gradient3 = gmls_evaluator.applyAlphasToDataAllComponentsAllTargetSites<double**, Kokkos::HostSpace>
            (sampling_data_device, GradientOfScalarPointEvaluation, PointSample, 
             true /*scalar_as_vector_if_needed*/, 3 /*evaluation site index*/);
    
    //! [Apply GMLS Alphas To Data]
    
    Kokkos::fence(); // let application of alphas to data finish before using results
    Kokkos::Profiling::popRegion();
    // times the Comparison in Kokkos
    Kokkos::Profiling::pushRegion("Comparison");
    
    //! [Check That Solutions Are Correct]
    
    // load value from output
    double GMLS_value;
    // load partial x from gradient
    double GMLS_GradX;
    // load partial y from gradient
    double GMLS_GradY;
    // load partial z from gradient
    double GMLS_GradZ;
    
    // loop through the target sites
    extra_evaluation_coordinates_count = 0;
    for (int i=0; i<number_target_coords; i++) {

        for (int k=0; k<(i%3)+1; ++k) {
            if (k==0) {
                // load value from output
                GMLS_value = output_value1(i);
                // load partial x from gradient
                GMLS_GradX = output_gradient1(i,0);
                // load partial y from gradient
                GMLS_GradY = (dimension>1) ? output_gradient1(i,1) : 0;
                // load partial z from gradient
                GMLS_GradZ = (dimension>2) ? output_gradient1(i,2) : 0;
            } else if (k==1) {
                // load value from output
                GMLS_value = output_value2(i);
                // load partial x from gradient
                GMLS_GradX = output_gradient2(i,0);
                // load partial y from gradient
                GMLS_GradY = (dimension>1) ? output_gradient2(i,1) : 0;
                // load partial z from gradient
                GMLS_GradZ = (dimension>2) ? output_gradient2(i,2) : 0;
            } else if (k==2) {
                // load value from output
                GMLS_value = output_value3(i);
                // load partial x from gradient
                GMLS_GradX = output_gradient3(i,0);
                // load partial y from gradient
                GMLS_GradY = (dimension>1) ? output_gradient3(i,1) : 0;
                // load partial z from gradient
                GMLS_GradZ = (dimension>2) ? output_gradient3(i,2) : 0;
            }
    
            // target site i's coordinate
            double xval = additional_target_coords(extra_evaluation_coordinates_count,0);
            double yval = additional_target_coords(extra_evaluation_coordinates_count,1);
            double zval = additional_target_coords(extra_evaluation_coordinates_count,2);
    
            // evaluation of various exact solutions
            double actual_value = trueSolution(xval, yval, zval, order, dimension);
    
            double actual_Gradient[3] = {0,0,0}; // initialized for 3, but only filled up to dimension
            trueGradient(actual_Gradient, xval, yval, zval, order, dimension);
    
            // check actual function value
            if(GMLS_value!=GMLS_value || std::abs(actual_value - GMLS_value) > failure_tolerance) {
                all_passed = false;
                std::cout << i << " Failed Actual by: " << std::abs(actual_value - GMLS_value) << " for evaluation site: " << k << std::endl;
            }
    
            // check gradient
            if(std::abs(actual_Gradient[0] - GMLS_GradX) > failure_tolerance) {
                all_passed = false;
                std::cout << i << " Failed GradX by: " << std::abs(actual_Gradient[0] - GMLS_GradX) << " for evaluation site: " << k << std::endl;
                if (dimension>1) {
                    if(std::abs(actual_Gradient[1] - GMLS_GradY) > failure_tolerance) {
                        all_passed = false;
                        std::cout << i << " Failed GradY by: " << std::abs(actual_Gradient[1] - GMLS_GradY) << " for evaluation site: " << k << std::endl;
                    }
                }
                if (dimension>2) {
                    if(std::abs(actual_Gradient[2] - GMLS_GradZ) > failure_tolerance) {
                        all_passed = false;
                        std::cout << i << " Failed GradZ by: " << std::abs(actual_Gradient[2] - GMLS_GradZ) << " for evaluation site: " << k << std::endl;
                    }
                }
            }
            extra_evaluation_coordinates_count++;
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
