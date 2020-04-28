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

#include "GMLS_Manifold.hpp"

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

// code block to reduce scope for all Kokkos View allocations
// otherwise, Views may be deallocating when we call Kokkos::finalize() later
{
    // check if 8 arguments are given from the command line, the first being the program name
    //  constraint_type used in solving each GMLS problem:
    //      0 - No constraints used in solving each GMLS problem
    //      1 - Neumann Gradient Scalar used in solving each GMLS problem
    int constraint_type = 0; // No constraints by default
    if (argc >= 8) {
        int arg8toi = atoi(args[7]);
        if (arg8toi > 0) {
            constraint_type = arg8toi;
        }
    }

    // check if 7 arguments are given from the command line, the first being the program name
    // problem_type used in solving each GMLS problem:
    //      0 - Standard GMLS problem
    //      1 - Manifold GMLS problem
    int problem_type = 1; // Manifold for this example
    if (argc >= 7) {
        int arg7toi = atoi(args[6]);
        if (arg7toi > 0) {
            problem_type = arg7toi;
        }
    }

    // check if 6 arguments are given from the command line, the first being the program name
    //  solver_type used for factorization in solving each GMLS problem:
    //      0 - SVD used for factorization in solving each GMLS problem
    //      1 - QR  used for factorization in solving each GMLS problem
    //      2 - LU  used for factorization in solving each GMLS problem
    int solver_type = 1; // QR by default
    if (argc >= 6) {
        int arg6toi = atoi(args[5]);
        if (arg6toi >= 0) {
            solver_type = arg6toi;
        }
    }

    // check if 5 arguments are given from the command line, the first being the program name
    //  N_pts_on_sphere used to determine spatial resolution
    int N_pts_on_sphere = 1000; // 1000 points by default
    if (argc >= 5) {
        int arg5toi = atoi(args[4]);
        if (arg5toi > 0) {
            N_pts_on_sphere = arg5toi;
        }
    }
    
    // check if 4 arguments are given from the command line
    //  dimension for the coordinates and the solution
    int dimension = 3; // dimension 3 by default
    if (argc >= 4) {
        int arg4toi = atoi(args[3]);
        if (arg4toi > 0) {
            dimension = arg4toi;
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
    //  set the number of target sites where we will reconstruct the target functionals at
    int order = 3; // 3rd degree polynomial basis by default
    if (argc >= 2) {
        int arg2toi = atoi(args[1]);
        if (arg2toi > 0) {
            order = arg2toi;
        }
    }
    
    // minimum neighbors for unisolvency is the same as the size of the polynomial basis 
    // dimension has one subtracted because it is a D-1 manifold represented in D dimensions
    const int min_neighbors = Compadre::GMLS::getNP(order, dimension-1);
    
    //! [Parse Command Line Arguments]
    Kokkos::Timer timer;
    Kokkos::Profiling::pushRegion("Setup Point Data");
    //! [Setting Up The Point Cloud]
    

    // coordinates of source sites, bigger than needed then resized later
    Kokkos::View<double**, Kokkos::DefaultExecutionSpace> source_coords_device("source coordinates", 
            1.25*N_pts_on_sphere, 3);
    Kokkos::View<double**>::HostMirror source_coords = Kokkos::create_mirror_view(source_coords_device);

    double r = 1.0;

    // set number of source coordinates from what was calculated
    int number_source_coords;

    { // fill source coordinates from surface of a sphere with quasiuniform points
        // algorithm described at https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
        int N_count = 0;
        double a = 4*PI*r*r/N_pts_on_sphere;
        double d = std::sqrt(a);
        int M_theta = std::round(PI/d);
        double d_theta = PI/M_theta;
        double d_phi = a/d_theta;
        for (int i=0; i<M_theta; ++i) {
            double theta = PI*(i + 0.5)/M_theta;
            int M_phi = std::round(2*PI*std::sin(theta)/d_phi);
            for (int j=0; j<M_phi; ++j) {
                double phi = 2*PI*j/M_phi;
                source_coords(N_count, 0) = theta;
                source_coords(N_count, 1) = phi;
                N_count++;
            }
        }

        for (int i=0; i<N_count; ++i) {
            double theta = source_coords(i,0);
            double phi = source_coords(i,1);
            source_coords(i,0) = r*std::sin(theta)*std::cos(phi);
            source_coords(i,1) = r*std::sin(theta)*std::sin(phi);
            source_coords(i,2) = r*cos(theta);
            //printf("%f %f %f\n", source_coords(i,0), source_coords(i,1), source_coords(i,2));
        }
        number_source_coords = N_count;
    }

    // resize source_coords to the size actually needed
    Kokkos::resize(source_coords, number_source_coords, 3);
    Kokkos::resize(source_coords_device, number_source_coords, 3);

    // coordinates of target sites
    Kokkos::View<double**, Kokkos::DefaultExecutionSpace> target_coords_device("target coordinates", 
            number_target_coords, 3);
    Kokkos::View<double**>::HostMirror target_coords = Kokkos::create_mirror_view(target_coords_device);

    // seed random number generator
    std::mt19937 rng(50);
    
    // generate random integers in [0...number_source_coords-1] (used to pick target sites)
    std::uniform_int_distribution<int> gen_num_neighbors(0, number_source_coords-1); // uniform, unbiased

    // fill target sites with random selections from source sites
    for (int i=0; i<number_target_coords; ++i) {
        const int source_site_to_copy = gen_num_neighbors(rng);
        for (int j=0; j<3; ++j) {
            target_coords(i,j) = source_coords(source_site_to_copy,j);
        }
    }


    //! [Setting Up The Point Cloud]
    
    Kokkos::Profiling::popRegion();
    Kokkos::fence();
    Kokkos::Profiling::pushRegion("Creating Data");
    
    //! [Creating The Data]
    
    
    // source coordinates need copied to device before using to construct sampling data
    Kokkos::deep_copy(source_coords_device, source_coords);
    Kokkos::deep_copy(target_coords_device, target_coords);

    // ensure that source coordinates are sent to device before evaluating sampling data based on them
    Kokkos::fence(); 

    
    // need Kokkos View storing true solution (for samples)
    Kokkos::View<double*, Kokkos::DefaultExecutionSpace> sampling_data_device("samples of true solution", 
            source_coords_device.extent(0));

    // need Kokkos View storing true vector solution (for samples)
    Kokkos::View<double**, Kokkos::DefaultExecutionSpace> sampling_vector_data_device("samples of vector true solution", 
            source_coords_device.extent(0), 3);
    
    Kokkos::parallel_for("Sampling Manufactured Solutions", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>
            (0,source_coords.extent(0)), KOKKOS_LAMBDA(const int i) {
    
        // coordinates of source site i
        double xval = source_coords_device(i,0);
        double yval = (dimension>1) ? source_coords_device(i,1) : 0;
        double zval = (dimension>2) ? source_coords_device(i,2) : 0;
    
        // data for targets with scalar input
        sampling_data_device(i) = sphere_harmonic54(xval, yval, zval);
        //printf("%f\n", sampling_data_device(i));

        for (int j=0; j<3; ++j) {
            double gradient[3] = {0,0,0};
            gradient_sphereHarmonic54_ambient(gradient, xval, yval, zval);
            sampling_vector_data_device(i,j) = gradient[j];
        }
        //printf("%f %f %f\n", sampling_vector_data_device(i,0), sampling_vector_data_device(i,1), sampling_vector_data_device(i,2));
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
    point_cloud_search.generateNeighborListsFromKNNSearch(false /*not dry run*/, target_coords, neighbor_lists, 
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
    
    // solver name for passing into the GMLS class
    std::string solver_name;
    if (solver_type == 0) { // SVD
        solver_name = "SVD";
    } else if (solver_type == 1) { // QR
        solver_name = "QR";
    } else if (solver_type == 2) { // LU
        solver_name = "LU";
    }

    // problem name for passing into the GMLS class
    std::string problem_name;
    if (problem_type == 0) { // Standard
        problem_name = "STANDARD";
    } else if (problem_type == 1) { // Manifold
        problem_name = "MANIFOLD";
    }

    // boundary name for passing into the GMLS class
    std::string constraint_name;
    if (constraint_type == 0) { // No constraints
        constraint_name = "NO_CONSTRAINT";
    } else if (constraint_type == 1) { // Neumann Gradient Scalar
        constraint_name = "NEUMANN_GRAD_SCALAR";
    }
    
    // initialize an instance of the GMLS class
    GMLS my_GMLS_vector_1(ReconstructionSpace::VectorTaylorPolynomial,
                          StaggeredEdgeIntegralSample,
                          order, dimension,
                          solver_name.c_str(), problem_name.c_str(), constraint_name.c_str(),
                          order /*manifold order*/);
    
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
    my_GMLS_vector_1.setProblemData(neighbor_lists_device, source_coords_device, target_coords_device, epsilon_device);
    
    // create a vector of target operations
    std::vector<TargetOperation> lro_vector_1(1);
    lro_vector_1[0] = DivergenceOfVectorPointEvaluation;

    // and then pass them to the GMLS class
    my_GMLS_vector_1.addTargets(lro_vector_1);

    // sets the weighting kernel function from WeightingFunctionType for curvature
    my_GMLS_vector_1.setCurvatureWeightingType(WeightingFunctionType::Power);
    
    // power to use in the weighting kernel function for curvature coefficients
    my_GMLS_vector_1.setCurvatureWeightingPower(2);
    
    // sets the weighting kernel function from WeightingFunctionType
    my_GMLS_vector_1.setWeightingType(WeightingFunctionType::Power);
    
    // power to use in that weighting kernel function
    my_GMLS_vector_1.setWeightingPower(2);

    // setup quadrature for StaggeredEdgeIntegralSample
    my_GMLS_vector_1.setOrderOfQuadraturePoints(2);
    my_GMLS_vector_1.setDimensionOfQuadraturePoints(1);
    my_GMLS_vector_1.setQuadratureType("LINE");
    
    // generate the alphas that to be combined with data for each target operation requested in lro
    my_GMLS_vector_1.generateAlphas();
    
    // initialize another instance of the GMLS class
    GMLS my_GMLS_vector_2(ReconstructionSpace::VectorTaylorPolynomial, 
                          StaggeredEdgeIntegralSample,
                          StaggeredEdgeAnalyticGradientIntegralSample,
                          order, dimension,
                          solver_name.c_str(), problem_name.c_str(), constraint_name.c_str(),
                          order /*manifold order*/);

    my_GMLS_vector_2.setProblemData(neighbor_lists_device, source_coords_device, target_coords_device, epsilon_device);
    std::vector<TargetOperation> lro_vector_2(2);
    lro_vector_2[0] = ChainedStaggeredLaplacianOfScalarPointEvaluation;
    lro_vector_2[1] = DivergenceOfVectorPointEvaluation;
    //lro_vector_2[2] = GradientOfScalarPointEvaluation;
    my_GMLS_vector_2.addTargets(lro_vector_2);
    my_GMLS_vector_2.setCurvatureWeightingType(WeightingFunctionType::Power);
    my_GMLS_vector_2.setCurvatureWeightingPower(2);
    my_GMLS_vector_2.setWeightingType(WeightingFunctionType::Power);
    my_GMLS_vector_2.setWeightingPower(2);
    my_GMLS_vector_2.setOrderOfQuadraturePoints(2);
    my_GMLS_vector_2.setDimensionOfQuadraturePoints(1);
    my_GMLS_vector_2.setQuadratureType("LINE");
    my_GMLS_vector_2.generateAlphas();

    // initialize another instance of the GMLS class
    GMLS my_GMLS_scalar(ReconstructionSpace::ScalarTaylorPolynomial,
                        StaggeredEdgeAnalyticGradientIntegralSample,
                        order, dimension,
                        solver_name.c_str(), problem_name.c_str(), constraint_name.c_str(),
                        order /*manifold order*/);

    my_GMLS_scalar.setProblemData(neighbor_lists_device, source_coords_device, target_coords_device, epsilon_device);

    std::vector<TargetOperation> lro_scalar(1);
    lro_scalar[0] = ChainedStaggeredLaplacianOfScalarPointEvaluation;
    //lro_scalar[1] = GradientOfScalarPointEvaluation;
    my_GMLS_scalar.addTargets(lro_scalar);
    my_GMLS_scalar.setCurvatureWeightingType(WeightingFunctionType::Power);
    my_GMLS_scalar.setCurvatureWeightingPower(2);
    my_GMLS_scalar.setWeightingType(WeightingFunctionType::Power);
    my_GMLS_scalar.setWeightingPower(2);
    my_GMLS_scalar.generateAlphas();


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
    Evaluator vector_1_gmls_evaluator(&my_GMLS_vector_1);
    Evaluator vector_2_gmls_evaluator(&my_GMLS_vector_2);
    Evaluator scalar_gmls_evaluator(&my_GMLS_scalar);
    

    //auto output_gradient_vectorbasis = 
    //    vector_2_gmls_evaluator.applyAlphasToDataAllComponentsAllTargetSites<double**, Kokkos::HostSpace>
    //        (sampling_data_device, GradientOfScalarPointEvaluation, StaggeredEdgeAnalyticGradientIntegralSample);

    //auto output_gradient_scalarbasis = 
    //    scalar_gmls_evaluator.applyAlphasToDataAllComponentsAllTargetSites<double**, Kokkos::HostSpace>
    //        (sampling_data_device, GradientOfScalarPointEvaluation, StaggeredEdgeAnalyticGradientIntegralSample);

    auto output_divergence_vectorsamples = 
        vector_1_gmls_evaluator.applyAlphasToDataAllComponentsAllTargetSites<double*, Kokkos::HostSpace>
            (sampling_vector_data_device, DivergenceOfVectorPointEvaluation, StaggeredEdgeIntegralSample);

    auto output_divergence_scalarsamples = 
        vector_2_gmls_evaluator.applyAlphasToDataAllComponentsAllTargetSites<double*, Kokkos::HostSpace>
            (sampling_data_device, DivergenceOfVectorPointEvaluation, StaggeredEdgeAnalyticGradientIntegralSample);

    auto output_laplacian_vectorbasis = 
        vector_2_gmls_evaluator.applyAlphasToDataAllComponentsAllTargetSites<double*, Kokkos::HostSpace>
            (sampling_data_device, ChainedStaggeredLaplacianOfScalarPointEvaluation, StaggeredEdgeAnalyticGradientIntegralSample);

    auto output_laplacian_scalarbasis = 
        scalar_gmls_evaluator.applyAlphasToDataAllComponentsAllTargetSites<double*, Kokkos::HostSpace>
            (sampling_data_device, ChainedStaggeredLaplacianOfScalarPointEvaluation, StaggeredEdgeAnalyticGradientIntegralSample);


    //! [Apply GMLS Alphas To Data]
    
    Kokkos::fence(); // let application of alphas to data finish before using results
    Kokkos::Profiling::popRegion();
    // times the Comparison in Kokkos
    Kokkos::Profiling::pushRegion("Comparison");
    
    //! [Check That Solutions Are Correct]
    
    double laplacian_vectorbasis_error = 0;
    double laplacian_vectorbasis_norm = 0;

    double laplacian_scalarbasis_error = 0;
    double laplacian_scalarbasis_norm = 0;

    double gradient_vectorbasis_ambient_error = 0;
    double gradient_vectorbasis_ambient_norm = 0;
    
    double gradient_scalarbasis_ambient_error = 0;
    double gradient_scalarbasis_ambient_norm = 0;

    double divergence_vectorsamples_ambient_error = 0;
    double divergence_vectorsamples_ambient_norm = 0;

    double divergence_scalarsamples_ambient_error = 0;
    double divergence_scalarsamples_ambient_norm = 0;
    
    // loop through the target sites
    for (int i=0; i<number_target_coords; i++) {
    
        // target site i's coordinate
        double xval = target_coords(i,0);
        double yval = (dimension>1) ? target_coords(i,1) : 0;
        double zval = (dimension>2) ? target_coords(i,2) : 0;
    
        // evaluation of various exact solutions
        double actual_Laplacian = laplace_beltrami_sphere_harmonic54(xval, yval, zval);
        double actual_Gradient_ambient[3] = {0,0,0}; // initialized for 3, but only filled up to dimension
        gradient_sphereHarmonic54_ambient(actual_Gradient_ambient, xval, yval, zval);
    
        laplacian_vectorbasis_error += (output_laplacian_vectorbasis(i) - actual_Laplacian)*(output_laplacian_vectorbasis(i) - actual_Laplacian);
        laplacian_vectorbasis_norm += actual_Laplacian*actual_Laplacian;

        //printf("Error of %f, %f vs %f\n", (output_laplacian_scalarbasis(i) - actual_Laplacian), output_laplacian_scalarbasis(i), actual_Laplacian);
        laplacian_scalarbasis_error += (output_laplacian_scalarbasis(i) - actual_Laplacian)*(output_laplacian_scalarbasis(i) - actual_Laplacian);
        laplacian_scalarbasis_norm += actual_Laplacian*actual_Laplacian;

        //for (int j=0; j<dimension; ++j) {
        //    //printf("VectorBasis Error of %f, %f vs %f\n", (output_gradient_vectorbasis(i,j) - actual_Gradient_ambient[j]), output_gradient_vectorbasis(i,j), actual_Gradient_ambient[j]);
        //    gradient_vectorbasis_ambient_error += (output_gradient_vectorbasis(i,j) - actual_Gradient_ambient[j])*(output_gradient_vectorbasis(i,j) - actual_Gradient_ambient[j]);
        //    gradient_vectorbasis_ambient_norm += actual_Gradient_ambient[j]*actual_Gradient_ambient[j];
        //}

        //for (int j=0; j<dimension; ++j) {
        //    //printf("ScalarBasis Error of %f, %f vs %f\n", (output_gradient_scalarbasis(i,j) - actual_Gradient_ambient[j]), output_gradient_scalarbasis(i,j), actual_Gradient_ambient[j]);
        //    gradient_scalarbasis_ambient_error += (output_gradient_scalarbasis(i,j) - actual_Gradient_ambient[j])*(output_gradient_scalarbasis(i,j) - actual_Gradient_ambient[j]);
        //    gradient_scalarbasis_ambient_norm += actual_Gradient_ambient[j]*actual_Gradient_ambient[j];
        //}

        //printf("Error of %f, %f vs %f\n", (output_divergence(i) - actual_Laplacian), output_divergence(i), actual_Laplacian);
        divergence_vectorsamples_ambient_error += (output_divergence_vectorsamples(i) - actual_Laplacian)*(output_divergence_vectorsamples(i) - actual_Laplacian);
        divergence_vectorsamples_ambient_norm += actual_Laplacian*actual_Laplacian;

        divergence_scalarsamples_ambient_error += (output_divergence_scalarsamples(i) - actual_Laplacian)*(output_divergence_scalarsamples(i) - actual_Laplacian);
        divergence_scalarsamples_ambient_norm += actual_Laplacian*actual_Laplacian;

    }

    laplacian_vectorbasis_error /= number_target_coords;
    laplacian_vectorbasis_error = std::sqrt(laplacian_vectorbasis_error);
    laplacian_vectorbasis_norm /= number_target_coords;
    laplacian_vectorbasis_norm = std::sqrt(laplacian_vectorbasis_norm);

    laplacian_scalarbasis_error /= number_target_coords;
    laplacian_scalarbasis_error = std::sqrt(laplacian_scalarbasis_error);
    laplacian_scalarbasis_norm /= number_target_coords;
    laplacian_scalarbasis_norm = std::sqrt(laplacian_scalarbasis_norm);

    gradient_vectorbasis_ambient_error /= number_target_coords;
    gradient_vectorbasis_ambient_error = std::sqrt(gradient_vectorbasis_ambient_error);
    gradient_vectorbasis_ambient_norm /= number_target_coords;
    gradient_vectorbasis_ambient_norm = std::sqrt(gradient_vectorbasis_ambient_norm);
    
    gradient_scalarbasis_ambient_error /= number_target_coords;
    gradient_scalarbasis_ambient_error = std::sqrt(gradient_scalarbasis_ambient_error);
    gradient_scalarbasis_ambient_norm /= number_target_coords;
    gradient_scalarbasis_ambient_norm = std::sqrt(gradient_scalarbasis_ambient_norm);

    divergence_vectorsamples_ambient_error /= number_target_coords;
    divergence_vectorsamples_ambient_error = std::sqrt(divergence_vectorsamples_ambient_error);
    divergence_vectorsamples_ambient_norm /= number_target_coords;
    divergence_vectorsamples_ambient_norm = std::sqrt(divergence_vectorsamples_ambient_norm);

    divergence_scalarsamples_ambient_error /= number_target_coords;
    divergence_scalarsamples_ambient_error = std::sqrt(divergence_scalarsamples_ambient_error);
    divergence_scalarsamples_ambient_norm /= number_target_coords;
    divergence_scalarsamples_ambient_norm = std::sqrt(divergence_scalarsamples_ambient_norm);

    printf("Staggered Laplace-Beltrami (VectorBasis) Error: %g\n", laplacian_vectorbasis_error / laplacian_vectorbasis_norm);  
    printf("Staggered Laplace-Beltrami (ScalarBasis) Error: %g\n", laplacian_scalarbasis_error / laplacian_scalarbasis_norm);  
    printf("Surface Staggered Gradient (VectorBasis) Error: %g\n", gradient_vectorbasis_ambient_error / gradient_vectorbasis_ambient_norm);  
    printf("Surface Staggered Gradient (ScalarBasis) Error: %g\n", gradient_scalarbasis_ambient_error / gradient_scalarbasis_ambient_norm);  
    printf("Surface Staggered Divergence (VectorSamples) Error: %g\n", divergence_vectorsamples_ambient_error / divergence_vectorsamples_ambient_norm);  
    printf("Surface Staggered Divergence (ScalarSamples) Error: %g\n", divergence_scalarsamples_ambient_error / divergence_scalarsamples_ambient_norm);  
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

return 0;

} // main


//! [Finalize Program]
