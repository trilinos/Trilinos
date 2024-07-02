// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
/*
 *
 * This examples tests the ability to reuse a GMLS class instance, 
 * changing the target and neighbor list for the new target.
 *
 */

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
    
    // Laplacian is a second order differential operator, which we expect to be slightly less accurate
    const double laplacian_failure_tolerance = 1e-9;
    
    // minimum neighbors for unisolvency is the same as the size of the polynomial basis 
    const int min_neighbors = Compadre::GMLS::getNP(order, dimension);
    
    //! [Parse Command Line Arguments]
    Kokkos::Timer timer;
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
    
    
    //! [Setting Up The Point Cloud]
    
    
    //! [Creating The Data]
    
    
    // source coordinates need copied to device before using to construct sampling data
    Kokkos::deep_copy(source_coords_device, source_coords);
    
    // target coordinates copied next, because it is a convenient time to send them to device
    Kokkos::deep_copy(target_coords_device, target_coords);
    
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
    
    
    //! [Setting Up The GMLS Object]
    
    // initialize an instance of the GMLS class 
    GMLS my_GMLS(VectorOfScalarClonesTaylorPolynomial, VectorPointSample,
                 order, dimension,
                 solver_name.c_str(), problem_name.c_str(), constraint_name.c_str(),
                 2 /*manifold order*/);
    
    // create a vector of target operations
    std::vector<TargetOperation> lro(5);
    lro[0] = ScalarPointEvaluation;
    lro[1] = LaplacianOfScalarPointEvaluation;
    lro[2] = GradientOfScalarPointEvaluation;
    lro[3] = DivergenceOfVectorPointEvaluation;
    lro[4] = CurlOfVectorPointEvaluation;
    
    // and then pass them to the GMLS class
    my_GMLS.addTargets(lro);
    
    // sets the weighting kernel function from WeightingFunctionType
    my_GMLS.setWeightingType(WeightingFunctionType::Power);
    
    // power to use in that weighting kernel function
    my_GMLS.setWeightingParameter(2);

    // set source sites once for all targets
    my_GMLS.setSourceSites(source_coords_device);

        
    // Point cloud construction for neighbor search
    // CreatePointCloudSearch constructs an object of type PointCloudSearch, but deduces the templates for you
    auto point_cloud_search(CreatePointCloudSearch(source_coords, dimension));

    // loop through the target sites
    for (int i=0; i<number_target_coords; i++) {
        timer.reset();
        // pass in neighbor lists, source coordinates, target coordinates, and window sizes
        //
        // single neighbor lists have the format:
        //      dimensions: (1 single neighbor list for one target) X (# maximum number of neighbors for any given target + 1)
        //                  the first column contains the number of neighbors for that rows corresponding target index
        //
        // source coordinates have the format:
        //      dimensions: (# number of source sites) X (dimension)
        //                  entries in the neighbor lists (integers) correspond to rows of this 2D array
        //
        // single target coordinates have the format:
        //      dimensions: (1 single target site) X (dimension)
        //                  # of target sites is same as # of rows of neighbor lists
        //
        
        
        // coordinates of single target sites
        Kokkos::View<double**, Kokkos::DefaultExecutionSpace> single_target_coords_device ("single target coordinates", 1, 3);
        Kokkos::View<double**>::HostMirror single_target_coords = Kokkos::create_mirror_view(single_target_coords_device);
        for (int j=0; j<3; ++j) {
            single_target_coords(0,j) = target_coords(i,j);
        }
        // target coordinates copied next, because it is a convenient time to send them to device
        Kokkos::deep_copy(single_target_coords_device, single_target_coords);
        Kokkos::fence();
    
        //! [Performing Neighbor Search]

        // each row is a neighbor list for a target site, with the first column of each row containing
        // the number of neighbors for that rows corresponding target site
        double epsilon_multiplier = 1.5;
        int estimated_upper_bound_number_neighbors = 
            point_cloud_search.getEstimatedNumberNeighborsUpperBound(min_neighbors, dimension, epsilon_multiplier);

        Kokkos::View<int**, Kokkos::DefaultExecutionSpace> single_neighbor_lists_device("neighbor lists", 
                1, estimated_upper_bound_number_neighbors); // first column is # of neighbors
        Kokkos::View<int**>::HostMirror single_neighbor_lists = Kokkos::create_mirror_view(single_neighbor_lists_device);
        
        // each target site has a window size
        Kokkos::View<double*, Kokkos::DefaultExecutionSpace> single_epsilon_device("h supports", 1);
        Kokkos::View<double*>::HostMirror single_epsilon = Kokkos::create_mirror_view(single_epsilon_device);
        
        // query the point cloud to generate the neighbor lists using a kdtree to produce the n nearest neighbor
        // to each target site, adding (epsilon_multiplier-1)*100% to whatever the distance away the further neighbor used is from
        // each target to the view for epsilon
        point_cloud_search.generate2DNeighborListsFromKNNSearch(false /*not dry run*/, single_target_coords, 
                single_neighbor_lists, single_epsilon, min_neighbors, epsilon_multiplier);
        
        //! [Performing Neighbor Search]
        
        
        // Copy data back to device (they were filled on the host)
        // We could have filled Kokkos Views with memory space on the host
        // and used these instead, and then the copying of data to the device
        // would be performed in the GMLS class
        Kokkos::deep_copy(single_neighbor_lists_device, single_neighbor_lists);
        Kokkos::deep_copy(single_epsilon_device, single_epsilon);
        Kokkos::fence(); // let call to build neighbor lists complete before copying back to device
        
        // Create temporary small views to hold just one target's information at a time
        my_GMLS.setNeighborLists(single_neighbor_lists_device);
        my_GMLS.setTargetSites(single_target_coords_device);
        my_GMLS.setWindowSizes(single_epsilon_device);
        
        
        // generate the alphas that to be combined with data for each target operation requested in lro
        my_GMLS.generateAlphas(1, true /* keep polynomial coefficients, only needed for a test later in this program */);
        
        
        //! [Setting Up The GMLS Object]
        
        double instantiation_time = timer.seconds();
        std::cout << "Took " << instantiation_time << "s to complete alphas generation." << std::endl;
        Kokkos::fence(); // let generateAlphas finish up before using alphas

        
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
        
        auto output_value = gmls_evaluator.applyAlphasToDataAllComponentsAllTargetSites<double*, Kokkos::HostSpace>
                (sampling_data_device, ScalarPointEvaluation);
        
        auto output_laplacian = gmls_evaluator.applyAlphasToDataAllComponentsAllTargetSites<double*, Kokkos::HostSpace>
                (sampling_data_device, LaplacianOfScalarPointEvaluation);
        
        auto output_gradient = gmls_evaluator.applyAlphasToDataAllComponentsAllTargetSites<double**, Kokkos::HostSpace>
                (sampling_data_device, GradientOfScalarPointEvaluation);
        
        auto output_divergence = gmls_evaluator.applyAlphasToDataAllComponentsAllTargetSites<double*, Kokkos::HostSpace>
                (gradient_sampling_data_device, DivergenceOfVectorPointEvaluation, VectorPointSample);
        
        auto output_curl = gmls_evaluator.applyAlphasToDataAllComponentsAllTargetSites<double**, Kokkos::HostSpace>
                (divergence_sampling_data_device, CurlOfVectorPointEvaluation);
        
        // retrieves polynomial coefficients instead of remapped field
        auto scalar_coefficients = gmls_evaluator.applyFullPolynomialCoefficientsBasisToDataAllComponents<double**, Kokkos::HostSpace>
                (sampling_data_device);
        
        //! [Apply GMLS Alphas To Data]
        
        Kokkos::fence(); // let application of alphas to data finish before using results
        // times the Comparison in Kokkos
        
        //! [Check That Solutions Are Correct]
    
        // load value from output
        double GMLS_value = output_value(0);
    
        // load laplacian from output
        double GMLS_Laplacian = output_laplacian(0);
    
        // load partial x from gradient
        // this is a test that the scalar_coefficients 2d array returned hold valid entries
        // scalar_coefficients(i,1)*1./epsilon(i) is equivalent to the target operation acting 
        // on the polynomials applied to the polynomial coefficients
        double GMLS_GradX = scalar_coefficients(0,1)*1./single_epsilon(0);
                            //output_gradient(i,0);
    
        // load partial y from gradient
        double GMLS_GradY = (dimension>1) ? output_gradient(0,1) : 0;
    
        // load partial z from gradient
        double GMLS_GradZ = (dimension>2) ? output_gradient(0,2) : 0;
    
        // load divergence from output
        double GMLS_Divergence = output_divergence(0);
    
        // load curl from output
        double GMLS_CurlX = (dimension>1) ? output_curl(0,0) : 0;
        double GMLS_CurlY = (dimension>1) ? output_curl(0,1) : 0;
        double GMLS_CurlZ = (dimension>2) ? output_curl(0,2) : 0;
    
    
        // target site i's coordinate
        double xval = target_coords(i,0);
        double yval = (dimension>1) ? target_coords(i,1) : 0;
        double zval = (dimension>2) ? target_coords(i,2) : 0;
    
        // evaluation of various exact solutions
        double actual_value = trueSolution(xval, yval, zval, order, dimension);
        double actual_Laplacian = trueLaplacian(xval, yval, zval, order, dimension);
    
        double actual_Gradient[3] = {0,0,0}; // initialized for 3, but only filled up to dimension
        trueGradient(actual_Gradient, xval, yval, zval, order, dimension);
    
        double actual_Divergence;
        actual_Divergence = trueLaplacian(xval, yval, zval, order, dimension);
    
        double actual_Curl[3] = {0,0,0}; // initialized for 3, but only filled up to dimension 
        // (and not at all for dimimension = 1)
        if (dimension>1) {
            actual_Curl[0] = curlTestSolution(xval, yval, zval, 0, dimension);
            actual_Curl[1] = curlTestSolution(xval, yval, zval, 1, dimension);
            if (dimension>2) {
                actual_Curl[2] = curlTestSolution(xval, yval, zval, 2, dimension);
            }
        }
    
        // check actual function value
        if(GMLS_value!=GMLS_value || std::abs(actual_value - GMLS_value) > failure_tolerance) {
            all_passed = false;
            std::cout << i << " Failed Actual by: " << std::abs(actual_value - GMLS_value) << std::endl;
        }
    
        // check Laplacian
        if(std::abs(actual_Laplacian - GMLS_Laplacian) > laplacian_failure_tolerance) {
            all_passed = false;
            std::cout << i <<" Failed Laplacian by: " << std::abs(actual_Laplacian - GMLS_Laplacian) << std::endl;
        }
    
        // check gradient
        if(std::abs(actual_Gradient[0] - GMLS_GradX) > failure_tolerance) {
            all_passed = false;
            std::cout << i << " Failed GradX by: " << std::abs(actual_Gradient[0] - GMLS_GradX) << std::endl;
            if (dimension>1) {
                if(std::abs(actual_Gradient[1] - GMLS_GradY) > failure_tolerance) {
                    all_passed = false;
                    std::cout << i << " Failed GradY by: " << std::abs(actual_Gradient[1] - GMLS_GradY) << std::endl;
                }
            }
            if (dimension>2) {
                if(std::abs(actual_Gradient[2] - GMLS_GradZ) > failure_tolerance) {
                    all_passed = false;
                    std::cout << i << " Failed GradZ by: " << std::abs(actual_Gradient[2] - GMLS_GradZ) << std::endl;
                }
            }
        }
    
        // check divergence
        if(std::abs(actual_Divergence - GMLS_Divergence) > failure_tolerance) {
            all_passed = false;
            std::cout << i << " Failed Divergence by: " << std::abs(actual_Divergence - GMLS_Divergence) << std::endl;
        }
    
        // check curl
        if (order > 2) { // reconstructed solution not in basis unless order greater than 2 used
            double tmp_diff = 0;
            if (dimension>1)
                tmp_diff += std::abs(actual_Curl[0] - GMLS_CurlX) + std::abs(actual_Curl[1] - GMLS_CurlY);
            if (dimension>2)
                tmp_diff += std::abs(actual_Curl[2] - GMLS_CurlZ);
            if(std::abs(tmp_diff) > failure_tolerance) {
                all_passed = false;
                std::cout << i << " Failed Curl by: " << std::abs(tmp_diff) << std::endl;
            }
        }
    }
    
    
    //! [Check That Solutions Are Correct] 
    // stop timing comparison loop
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
