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
    const double failure_tolerance = 9e-8;

    // minimum neighbors for unisolvency is the same as the size of the polynomial basis
    const int min_neighbors = Compadre::GMLS::getNP(order+1, dimension);

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

    // tangent bundle for each target sites
    Kokkos::View<double***, Kokkos::DefaultExecutionSpace> tangent_bundles_device ("tangent bundles", number_target_coords, 3, 3);
    Kokkos::View<double***>::HostMirror tangent_bundles = Kokkos::create_mirror_view(tangent_bundles_device);

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
        if (constraint_name == "NEUMANN_GRAD_SCALAR") { // create bundles of normal vectors
            if (dimension == 3) {
                tangent_bundles(i, 0, 0) = 0.0;
                tangent_bundles(i, 0, 1) = 0.0;
                tangent_bundles(i, 0, 2) = 0.0;
                tangent_bundles(i, 1, 0) = 0.0;
                tangent_bundles(i, 1, 1) = 0.0;
                tangent_bundles(i, 1, 2) = 0.0;
                tangent_bundles(i, 2, 0) = 1.0/(sqrt(3.0));
                tangent_bundles(i, 2, 1) = 1.0/(sqrt(3.0));
                tangent_bundles(i, 2, 2) = 1.0/(sqrt(3.0));
            }
            if (dimension == 2) {
                tangent_bundles(i, 0, 0) = 0.0;
                tangent_bundles(i, 0, 1) = 0.0;
                tangent_bundles(i, 1, 0) = 1.0/(sqrt(2.0));
                tangent_bundles(i, 1, 1) = 1.0/(sqrt(2.0));
            }
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
    Kokkos::View<double*, Kokkos::DefaultExecutionSpace> sampling_data_device("samples of true solution",
            source_coords_device.extent(0));

    Kokkos::parallel_for("Sampling Manufactured Solutions", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>
            (0,source_coords.extent(0)), KOKKOS_LAMBDA(const int i) {
        // coordinates of source site i
        double xval = source_coords_device(i,0);
        double yval = (dimension>1) ? source_coords_device(i,1) : 0;
        double zval = (dimension>2) ? source_coords_device(i,2) : 0;

        // data for targets with scalar input
        sampling_data_device(i) = trueSolution(xval, yval, zval, order, dimension);
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
    // First, analytica gradient on scalar polynomial basis
    GMLS scalar_basis_gmls(ScalarTaylorPolynomial,
                           StaggeredEdgeAnalyticGradientIntegralSample,
                           order, dimension,
                           solver_name.c_str(), problem_name.c_str(), constraint_name.c_str(),
                           0 /*manifold order*/);

    // Another class performing Gaussian quadrature integration on vector polynomial basis
    GMLS vector_basis_gmls(VectorTaylorPolynomial,
                           StaggeredEdgeIntegralSample,
                           StaggeredEdgeAnalyticGradientIntegralSample,
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
    scalar_basis_gmls.setProblemData(neighbor_lists_device, source_coords_device, target_coords_device, epsilon_device);
    vector_basis_gmls.setProblemData(neighbor_lists_device, source_coords_device, target_coords_device, epsilon_device);
    if (constraint_name == "NEUMANN_GRAD_SCALAR") {
        scalar_basis_gmls.setTangentBundle(tangent_bundles_device);
        vector_basis_gmls.setTangentBundle(tangent_bundles_device);
    }

    // create a vector of target operations
    std::vector<TargetOperation> lro(2);
    lro[0] = DivergenceOfVectorPointEvaluation;
    lro[1] = GradientOfScalarPointEvaluation;

    // and then pass them to the GMLS class
    scalar_basis_gmls.addTargets(lro);
    vector_basis_gmls.addTargets(lro);

    // sets the weighting kernel function from WeightingFunctionType
    scalar_basis_gmls.setWeightingType(WeightingFunctionType::Power);
    vector_basis_gmls.setWeightingType(WeightingFunctionType::Power);

    // power to use in that weighting kernel function
    scalar_basis_gmls.setWeightingParameter(2);
    vector_basis_gmls.setWeightingParameter(2);

    // setup quadrature for StaggeredEdgeIntegralSample
    vector_basis_gmls.setOrderOfQuadraturePoints(order);
    vector_basis_gmls.setDimensionOfQuadraturePoints(1);
    vector_basis_gmls.setQuadratureType("LINE");

    // generate the alphas that to be combined with data for each target operation requested in lro
    scalar_basis_gmls.generateAlphas();
    vector_basis_gmls.generateAlphas();

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
    Evaluator gmls_evaluator_scalar(&scalar_basis_gmls);
    Evaluator gmls_evaluator_vector(&vector_basis_gmls);

    auto output_divergence_scalar = gmls_evaluator_scalar.applyAlphasToDataAllComponentsAllTargetSites<double*, Kokkos::HostSpace>
      (sampling_data_device, DivergenceOfVectorPointEvaluation, StaggeredEdgeAnalyticGradientIntegralSample);

    auto output_divergence_vector = gmls_evaluator_vector.applyAlphasToDataAllComponentsAllTargetSites<double*, Kokkos::HostSpace>
      (sampling_data_device, DivergenceOfVectorPointEvaluation, StaggeredEdgeAnalyticGradientIntegralSample);

    // evaluating the gradient
    auto output_gradient_scalar = gmls_evaluator_scalar.applyAlphasToDataAllComponentsAllTargetSites<double**, Kokkos::HostSpace>
      (sampling_data_device, GradientOfScalarPointEvaluation, StaggeredEdgeAnalyticGradientIntegralSample);

    auto output_gradient_vector = gmls_evaluator_vector.applyAlphasToDataAllComponentsAllTargetSites<double**, Kokkos::HostSpace>
      (sampling_data_device, GradientOfScalarPointEvaluation, StaggeredEdgeAnalyticGradientIntegralSample);

    //! [Apply GMLS Alphas To Data]

    Kokkos::fence(); // let application of alphas to data finish before using results
    Kokkos::Profiling::popRegion();
    // times the Comparison in Kokkos
    Kokkos::Profiling::pushRegion("Comparison");

    //! [Check That Solutions Are Correct]

    // loop through the target sites
    for (int i=0; i<number_target_coords; i++) {

        // target site i's coordinate
        double xval = target_coords(i,0);
        double yval = (dimension>1) ? target_coords(i,1) : 0;
        double zval = (dimension>2) ? target_coords(i,2) : 0;

        // evaluation of various exact solutions
        double actual_Divergence;
        actual_Divergence = trueLaplacian(xval, yval, zval, order, dimension);
        double actual_Gradient[3] = {0,0,0}; // initialized for 3, but only filled up to dimension
        trueGradient(actual_Gradient, xval, yval, zval, order, dimension);

        // calculate correction factor
        double g = (dimension == 3) ? (1.0/sqrt(3.0))*(actual_Gradient[0] + actual_Gradient[1] + actual_Gradient[2]) 
                                        : (1.0/sqrt(2.0))*(actual_Gradient[0] + actual_Gradient[1]);
        // obtain number of neighbor for each target
        // in order to exploit the index where the value for the Lagrange multiplier is stored
        int num_neigh_i = neighbor_lists(i, 0);

        // load divergence from output
        double GMLS_Divergence_Scalar = output_divergence_scalar(i);
        double GMLS_Divergence_Vector = output_divergence_vector(i);

        // obtain adjusted value for divergence
        if (constraint_name == "NEUMANN_GRAD_SCALAR") {
            double b_i_scalar = scalar_basis_gmls.getSolutionSetHost()->getAlpha0TensorTo0Tensor(DivergenceOfVectorPointEvaluation, i, num_neigh_i);
            GMLS_Divergence_Scalar = GMLS_Divergence_Scalar + b_i_scalar*g;

            double b_i_vector = vector_basis_gmls.getSolutionSetHost()->getAlpha0TensorTo0Tensor(DivergenceOfVectorPointEvaluation, i, num_neigh_i);
            GMLS_Divergence_Vector = GMLS_Divergence_Vector + b_i_vector*g;
        }

        // check divergence - scalar basis
        if(std::abs(actual_Divergence - GMLS_Divergence_Scalar) > failure_tolerance) {
            all_passed = false;
            std::cout << i << " Failed Divergence on SCALAR basis by: " << std::abs(actual_Divergence - GMLS_Divergence_Scalar) << std::endl;
            std::cout << i << " GMLS " << GMLS_Divergence_Scalar << " adjusted " << GMLS_Divergence_Scalar << " actual " << actual_Divergence << std::endl;
        }

        // check divergence - vector basis
        if(std::abs(actual_Divergence - GMLS_Divergence_Vector) > failure_tolerance) {
            all_passed = false;
            std::cout << i << " Failed Divergence on VECTOR basis by: " << std::abs(actual_Divergence - GMLS_Divergence_Vector) << std::endl;
            std::cout << i << " GMLS " << GMLS_Divergence_Vector << " adjusted " << GMLS_Divergence_Vector << " actual " << actual_Divergence << std::endl;
        }

        // load gradient from output
        double Scalar_GMLS_GradX = (dimension>1) ? output_gradient_scalar(i,0) : 0;
        double Scalar_GMLS_GradY = (dimension>1) ? output_gradient_scalar(i,1) : 0;
        double Scalar_GMLS_GradZ = (dimension>2) ? output_gradient_scalar(i,2) : 0;

        double Vector_GMLS_GradX = (dimension>1) ? output_gradient_vector(i,0) : 0;
        double Vector_GMLS_GradY = (dimension>1) ? output_gradient_vector(i,1) : 0;
        double Vector_GMLS_GradZ = (dimension>2) ? output_gradient_vector(i,2) : 0;

        // Obtain adjusted value
        if (constraint_name == "NEUMANN_GRAD_SCALAR") {
            double bx_i_scalar = scalar_basis_gmls.getSolutionSetHost()->getAlpha0TensorTo1Tensor(GradientOfScalarPointEvaluation, i, 0, num_neigh_i);
            double by_i_scalar = scalar_basis_gmls.getSolutionSetHost()->getAlpha0TensorTo1Tensor(GradientOfScalarPointEvaluation, i, 1, num_neigh_i);
            double bz_i_scalar = scalar_basis_gmls.getSolutionSetHost()->getAlpha0TensorTo1Tensor(GradientOfScalarPointEvaluation, i, 2, num_neigh_i);
            Scalar_GMLS_GradX = Scalar_GMLS_GradX + bx_i_scalar*g;
            Scalar_GMLS_GradY = Scalar_GMLS_GradY + by_i_scalar*g;
            Scalar_GMLS_GradZ = Scalar_GMLS_GradZ + bz_i_scalar*g;

            double bx_i_vector = vector_basis_gmls.getSolutionSetHost()->getAlpha0TensorTo1Tensor(GradientOfScalarPointEvaluation, i, 0, num_neigh_i);
            double by_i_vector = vector_basis_gmls.getSolutionSetHost()->getAlpha0TensorTo1Tensor(GradientOfScalarPointEvaluation, i, 1, num_neigh_i);
            double bz_i_vector = vector_basis_gmls.getSolutionSetHost()->getAlpha0TensorTo1Tensor(GradientOfScalarPointEvaluation, i, 2, num_neigh_i);
            Vector_GMLS_GradX = Vector_GMLS_GradX + bx_i_vector*g;
            Vector_GMLS_GradY = Vector_GMLS_GradY + by_i_vector*g;
            Vector_GMLS_GradZ = Vector_GMLS_GradZ + bz_i_vector*g;
        }

        // check gradient - scalar gmls
        if(std::abs(actual_Gradient[0] - Scalar_GMLS_GradX) > failure_tolerance) {
            all_passed = false;
            std::cout << i << " Failed Scalar GradX by: " << std::abs(actual_Gradient[0] - Scalar_GMLS_GradX) << std::endl;
            if (dimension>1) {
                if(std::abs(actual_Gradient[1] - Scalar_GMLS_GradY) > failure_tolerance) {
                    all_passed = false;
                    std::cout << i << " Failed Scalar GradY by: " << std::abs(actual_Gradient[1] - Scalar_GMLS_GradY) << std::endl;
                }
            }
            if (dimension>2) {
                if(std::abs(actual_Gradient[2] - Scalar_GMLS_GradZ) > failure_tolerance) {
                    all_passed = false;
                    std::cout << i << " Failed Scalar GradZ by: " << std::abs(actual_Gradient[2] - Scalar_GMLS_GradZ) << std::endl;
                }
            }
        }

        // check gradient - vector gmls
        if(std::abs(actual_Gradient[0] - Vector_GMLS_GradX) > failure_tolerance) {
            all_passed = false;
            std::cout << i << " Failed Vector GradX by: " << std::abs(actual_Gradient[0] - Vector_GMLS_GradX) << std::endl;
            if (dimension>1) {
                if(std::abs(actual_Gradient[1] - Vector_GMLS_GradY) > failure_tolerance) {
                    all_passed = false;
                    std::cout << i << " Failed Vector GradY by: " << std::abs(actual_Gradient[1] - Vector_GMLS_GradY) << std::endl;
                }
            }
            if (dimension>2) {
                if(std::abs(actual_Gradient[2] - Vector_GMLS_GradZ) > failure_tolerance) {
                    all_passed = false;
                    std::cout << i << " Failed Vector GradZ by: " << std::abs(actual_Gradient[2] - Vector_GMLS_GradZ) << std::endl;
                }
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
