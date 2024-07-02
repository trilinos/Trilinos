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

#include "GMLS_Manifold.hpp"
#include "CommandLineProcessor.hpp"

#ifdef COMPADRE_USE_MPI
#include <mpi.h>
#endif

#include <Kokkos_Timer.hpp>
#include <Kokkos_Core.hpp>

using namespace Compadre;
//! [Ambient to Local Back To Ambient Helper Function]
void AmbientLocalAmbient(XYZ& u, double* T_data, double* P_data) {
    //  reconstructions with vector output on a manifold are defined
    //  in their local chart (2D). Next, they are mapped to 3D using
    //  the tangent vector calculated at the target site.
    //
    //  We are interested in a mapping using the tangent vector at 
    //  additional evaluation site instead, so we map back to the 
    //  local chart using the tangent vector defined at the target 
    //  site (T), then mapping from the local chart to ambient space
    //  using the tangent vector calculated at the extra site (P).
    

    host_scratch_matrix_right_type T(T_data, 3, 3);
    host_scratch_matrix_right_type P(P_data, 3, 3);

    // first get back to local chart
    double local_vec[3] = {0,0};
    for (int j=0; j<3; ++j) {
        local_vec[0] += T(0,j) * u[j];
        local_vec[1] += T(1,j) * u[j];
    }
    // second go to ambient space using tangent for first neighbor
    for (int j=0; j<3; ++j) u[j] = 0;
    for (int j=0; j<3; ++j) {
        u[j] += P(0, j) * local_vec[0];
        u[j] += P(1, j) * local_vec[1];
    }
    
}

//! [Ambient to Local Back To Ambient Helper Function]

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

    CommandLineProcessor clp(argc, args);
    auto order = clp.order;
    auto dimension = clp.dimension;
    auto number_target_coords = clp.number_target_coords;
    auto constraint_name = clp.constraint_name;
    auto solver_name = clp.solver_name;
    auto problem_name = clp.problem_name;
    int N_pts_on_sphere = (clp.number_source_coords>=0) ? clp.number_source_coords : 1000;
    
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

    {
        // fill source coordinates from surface of a sphere with quasiuniform points
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
    Kokkos::View<double**, Kokkos::DefaultExecutionSpace> target_coords_device ("target coordinates", 
            number_target_coords, 3);
    Kokkos::View<double**>::HostMirror target_coords = Kokkos::create_mirror_view(target_coords_device);
    
    {
        bool enough_pts_found = false;
        // fill target coordinates from surface of a sphere with quasiuniform points
        // stop when enough points are found
        int N_count = 0;
        double a = 4.0*PI*r*r/((double)(5.0*number_target_coords)); // 5.0 is in case number_target_coords is set to something very small (like 1)
        double d = std::sqrt(a);
        int M_theta = std::round(PI/d);
        double d_theta = PI/((double)M_theta);
        double d_phi = a/d_theta;
        for (int i=0; i<M_theta; ++i) {
            double theta = PI*(i + 0.5)/M_theta;
            int M_phi = std::round(2*PI*std::sin(theta)/d_phi);
            for (int j=0; j<M_phi; ++j) {
                double phi = 2*PI*j/M_phi;
                target_coords(N_count, 0) = theta;
                target_coords(N_count, 1) = phi;
                N_count++;
                if (N_count == number_target_coords) {
                    enough_pts_found = true;
                    break;
                }
            }
            if (enough_pts_found) break;
        }

        for (int i=0; i<N_count; ++i) {
            double theta = target_coords(i,0);
            double phi = target_coords(i,1);
            target_coords(i,0) = r*std::sin(theta)*std::cos(phi);
            target_coords(i,1) = r*std::sin(theta)*std::sin(phi);
            target_coords(i,2) = r*cos(theta);
            //printf("%f %f %f\n", target_coords(i,0), target_coords(i,1), target_coords(i,2));
        }
    }

    
    //! [Setting Up The Point Cloud]
    
    Kokkos::Profiling::popRegion();
    Kokkos::fence();
    Kokkos::Profiling::pushRegion("Creating Data");
    
    //! [Creating The Data]
    
    
    // source coordinates need copied to device before using to construct sampling data
    Kokkos::deep_copy(source_coords_device, source_coords);
    
    // target coordinates copied next, because it is a convenient time to send them to device
    Kokkos::deep_copy(target_coords_device, target_coords);

    // ensure that source coordinates are sent to device before evaluating sampling data based on them
    Kokkos::fence(); 

    
    // need Kokkos View storing true solution (for samples)
    Kokkos::View<double*, Kokkos::DefaultExecutionSpace> sampling_data_device("samples of true solution", 
            source_coords_device.extent(0));
    Kokkos::View<double*, Kokkos::DefaultExecutionSpace> ones_data_device("samples of ones", 
            source_coords_device.extent(0));
    Kokkos::deep_copy(ones_data_device, 1.0);

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

        for (int j=0; j<3; ++j) {
            double gradient[3] = {0,0,0};
            gradient_sphereHarmonic54_ambient(gradient, xval, yval, zval);
            sampling_vector_data_device(i,j) = gradient[j];
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
    double epsilon_multiplier = 1.7;
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
    
    // initialize an instance of the GMLS class for problems with a scalar basis and 
    // traditional point sampling as the sampling functional
    GMLS my_GMLS_scalar(order, dimension,
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
    my_GMLS_scalar.setProblemData(neighbor_lists_device, source_coords_device, target_coords_device, epsilon_device);

    // set up additional sites to evaluate target operators
    // these sites will be the neighbors of the target site
    my_GMLS_scalar.setAdditionalEvaluationSitesData(neighbor_lists_device, source_coords_device);

    // set a reference outward normal direction, causing a surface orientation after
    // the GMLS instance computes an approximate tangent bundle
    // on a sphere, the ambient coordinates are the outward normal direction
    my_GMLS_scalar.setReferenceOutwardNormalDirection(target_coords_device, true /* use to orient surface */);
    
    // create a vector of target operations
    std::vector<TargetOperation> lro_scalar(3);
    lro_scalar[0] = ScalarPointEvaluation;
    lro_scalar[1] = GaussianCurvaturePointEvaluation;
    lro_scalar[2] = CurlOfVectorPointEvaluation;

    // and then pass them to the GMLS class
    my_GMLS_scalar.addTargets(lro_scalar);

    // sets the weighting kernel function from WeightingFunctionType for curvature
    my_GMLS_scalar.setCurvatureWeightingType(WeightingFunctionType::Power);
    
    // power to use in the weighting kernel function for curvature coefficients
    my_GMLS_scalar.setCurvatureWeightingParameter(2);
    
    // sets the weighting kernel function from WeightingFunctionType
    my_GMLS_scalar.setWeightingType(WeightingFunctionType::Power);
    
    // power to use in that weighting kernel function
    my_GMLS_scalar.setWeightingParameter(2);
    
    // generate the alphas that to be combined with data for each target operation requested in lro
    my_GMLS_scalar.generateAlphas();

    Kokkos::Profiling::pushRegion("Full Polynomial Basis GMLS Solution");
    // initialize another instance of the GMLS class for problems with a vector basis on a manifold and point 
    // evaluation of that vector as the sampling functional
    // VectorTaylorPolynomial indicates that the basis will be a polynomial with as many components as the
    // dimension of the manifold. This differs from another possibility, which follows this class.
    GMLS my_GMLS_vector(ReconstructionSpace::VectorTaylorPolynomial, VaryingManifoldVectorPointSample,
                          order, dimension,
                          solver_name.c_str(), problem_name.c_str(), constraint_name.c_str(),
                          order /*manifold order*/);

    my_GMLS_vector.setProblemData(neighbor_lists_device, source_coords_device, target_coords_device, epsilon_device);
    my_GMLS_vector.setAdditionalEvaluationSitesData(neighbor_lists_device, source_coords_device);
    std::vector<TargetOperation> lro_vector(2);
    lro_vector[0] = VectorPointEvaluation;
    //lro_vector[1] = DivergenceOfVectorPointEvaluation;
    my_GMLS_vector.addTargets(lro_vector);
    my_GMLS_vector.setCurvatureWeightingType(WeightingFunctionType::Power);
    my_GMLS_vector.setCurvatureWeightingParameter(2);
    my_GMLS_vector.setWeightingType(WeightingFunctionType::Power);
    my_GMLS_vector.setWeightingParameter(2);
    my_GMLS_vector.generateAlphas();
    Kokkos::Profiling::popRegion();

    Kokkos::Profiling::pushRegion("Scalar Polynomial Basis Repeated to Form a Vector GMLS Solution");
    // initialize another instance of the GMLS class for problems with a vector basis on a manifold and point 
    // evaluation of that vector as the sampling functional
    // VectorOfScalarClonesTaylorPolynomial indicates a scalar polynomial will be solved for, since
    // each component of the reconstructed vector are independent. This results in solving a smaller system
    // for each GMLS problem, and is the suggested way to do vector reconstructions when sampling functionals
    // acting on the basis would not create non-zero offdiagonal blocks.
    //
    // As an example, consider a 2D manifold in 3D ambient space. The GMLS problem is posed in the local chart,
    // meaning that the problem being solved looks like 
    //
    //  [P_0  0 | where P_0 has dimension #number of neighbors for a target X #dimension of a scalar basis
    //  | 0  P_1]
    //
    //  P_1 is similarly shaped, and for sampling functional that is a point evaluation, P_0 and P_1 are 
    //  identical and their degrees of freedom in this system are disjoint, allowing us to solve for the
    //  degrees of freedom for either block independently. Additionally, the will produce the exact
    //  same polynomial coefficients for the point sampling functional, therefore it makes sense to use
    //  VectorOfScalarClonesTaylorPolynomial.
    //
    //  In the print-out for this program, we include the timings and errors on this and VectorTaylorPolynomial
    //  in order to demonstrate that they produce exactly the same answer, but that one is much more efficient.
    //
    GMLS my_GMLS_vector_of_scalar_clones(ReconstructionSpace::VectorOfScalarClonesTaylorPolynomial, VaryingManifoldVectorPointSample,
                                         order, dimension,
                                         solver_name.c_str(), problem_name.c_str(), constraint_name.c_str(),
                                         order /*manifold order*/);

    my_GMLS_vector_of_scalar_clones.setProblemData(neighbor_lists_device, source_coords_device, target_coords_device, epsilon_device);
    my_GMLS_vector_of_scalar_clones.setAdditionalEvaluationSitesData(neighbor_lists_device, source_coords_device);
    std::vector<TargetOperation> lro_vector_of_scalar_clones(2);
    lro_vector_of_scalar_clones[0] = VectorPointEvaluation;
    //lro_vector_of_scalar_clones[1] = DivergenceOfVectorPointEvaluation;
    my_GMLS_vector_of_scalar_clones.addTargets(lro_vector_of_scalar_clones);
    my_GMLS_vector_of_scalar_clones.setCurvatureWeightingType(WeightingFunctionType::Power);
    my_GMLS_vector_of_scalar_clones.setCurvatureWeightingParameter(2);
    my_GMLS_vector_of_scalar_clones.setWeightingType(WeightingFunctionType::Power);
    my_GMLS_vector_of_scalar_clones.setWeightingParameter(2);
    my_GMLS_vector_of_scalar_clones.generateAlphas();
    Kokkos::Profiling::popRegion();


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
    Evaluator scalar_gmls_evaluator(&my_GMLS_scalar);
    Evaluator vector_gmls_evaluator(&my_GMLS_vector);
    Evaluator vector_gmls_evaluator_of_scalar_clones(&my_GMLS_vector_of_scalar_clones);
    
    auto output_value = scalar_gmls_evaluator.applyAlphasToDataAllComponentsAllTargetSites<double*, Kokkos::HostSpace>
            (sampling_data_device, ScalarPointEvaluation, PointSample, 
             true /*scalar_as_vector_if_needed*/, 1 /*evaluation site index*/);

    auto output_gaussian_curvature = scalar_gmls_evaluator.applyAlphasToDataAllComponentsAllTargetSites<double*, Kokkos::HostSpace>
            (ones_data_device, GaussianCurvaturePointEvaluation, PointSample, 
             true /*scalar_as_vector_if_needed*/, 1 /*evaluation site index*/);

    auto output_curl = scalar_gmls_evaluator.applyAlphasToDataAllComponentsAllTargetSites<double**, Kokkos::HostSpace>
            (sampling_data_device, CurlOfVectorPointEvaluation, PointSample, 
             true /*scalar_as_vector_if_needed*/, 1 /*evaluation site index*/);
    
    //auto output_laplacian = scalar_gmls_evaluator.applyAlphasToDataAllComponentsAllTargetSites<double*, Kokkos::HostSpace>
    //        (sampling_data_device, LaplacianOfScalarPointEvaluation, PointSample, 
    //         true /*scalar_as_vector_if_needed*/, 1 /*evaluation site index*/);

    //auto output_gradient = scalar_gmls_evaluator.applyAlphasToDataAllComponentsAllTargetSites<double**, Kokkos::HostSpace>
    //        (sampling_data_device, GradientOfScalarPointEvaluation, PointSample, 
    //         true /*scalar_as_vector_if_needed*/, 1 /*evaluation site index*/);

    auto output_vector = vector_gmls_evaluator.applyAlphasToDataAllComponentsAllTargetSites<double**, Kokkos::HostSpace>
            (sampling_vector_data_device, VectorPointEvaluation, VaryingManifoldVectorPointSample,
             true /*scalar_as_vector_if_needed*/, 1 /*evaluation site index*/);
    
    //auto output_divergence = vector_gmls_evaluator.applyAlphasToDataAllComponentsAllTargetSites<double*, Kokkos::HostSpace>
    //        (sampling_vector_data_device, DivergenceOfVectorPointEvaluation, ManifoldVectorPointSample,
    //         true /*scalar_as_vector_if_needed*/, 1 /*evaluation site index*/);

    auto output_vector_of_scalar_clones = 
        vector_gmls_evaluator_of_scalar_clones.applyAlphasToDataAllComponentsAllTargetSites<double**, 
        Kokkos::HostSpace>(sampling_vector_data_device, VectorPointEvaluation, VaryingManifoldVectorPointSample,
             true /*scalar_as_vector_if_needed*/, 1 /*evaluation site index*/);
    
    //auto output_divergence_of_scalar_clones = 
    //    vector_gmls_evaluator_of_scalar_clones.applyAlphasToDataAllComponentsAllTargetSites<double*, 
    //    Kokkos::HostSpace>(sampling_vector_data_device, DivergenceOfVectorPointEvaluation, ManifoldVectorPointSample,
    //         true /*scalar_as_vector_if_needed*/, 1 /*evaluation site index*/);


    //    Kokkos::fence(); // let application of alphas to data finish before using results
    //
    //// move gradient data to device so that it can be transformed into velocity
    //auto output_gradient_device_mirror = Kokkos::create_mirror(Kokkos::DefaultExecutionSpace::memory_space(), output_gradient);
    //Kokkos::deep_copy(output_gradient_device_mirror, output_gradient);
    //Kokkos::parallel_for("Create Velocity From Surface Gradient", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>
    //        (0,target_coords.extent(0)), KOKKOS_LAMBDA(const int i) {
    //
    //    // coordinates of target site i
    //    double xval = target_coords_device(i,0);
    //    double yval = (dimension>1) ? target_coords_device(i,1) : 0;
    //    double zval = (dimension>2) ? target_coords_device(i,2) : 0;

    //    double gradx = output_gradient_device_mirror(i,0);
    //    double grady = output_gradient_device_mirror(i,1);
    //    double gradz = output_gradient_device_mirror(i,2);
    //
    //    // overwrites gradient with velocity
    //    output_gradient_device_mirror(i,0) = (grady*zval - yval*gradz);
    //    output_gradient_device_mirror(i,1) = (-(gradx*zval - xval*gradz));
    //    output_gradient_device_mirror(i,2) = (gradx*yval - xval*grady);
    //
    //});
    //Kokkos::deep_copy(output_gradient, output_gradient_device_mirror);


    //! [Apply GMLS Alphas To Data]
    
    Kokkos::fence(); // let application of alphas to data finish before using results
    Kokkos::Profiling::popRegion();
    // times the Comparison in Kokkos
    Kokkos::Profiling::pushRegion("Comparison");
    
    //! [Check That Solutions Are Correct]
    
    double tangent_bundle_error = 0;
    double tangent_bundle_norm = 0;
    double values_error = 0;
    double values_norm = 0;
    double gc_error = 0;
    double gc_norm = 0;
    double curl_ambient_error = 0;
    double curl_ambient_norm = 0;
    //double laplacian_error = 0;
    //double laplacian_norm = 0;
    //double gradient_ambient_error = 0;
    //double gradient_ambient_norm = 0;
    double vector_ambient_error = 0;
    double vector_ambient_norm = 0;
    //double divergence_ambient_error = 0;
    //double divergence_ambient_norm = 0;
    double vector_of_scalar_clones_ambient_error = 0;
    double vector_of_scalar_clones_ambient_norm = 0;
    //double divergence_of_scalar_clones_ambient_error = 0;
    //double divergence_of_scalar_clones_ambient_norm = 0;
 
    // tangent vectors for each source coordinate are stored here
    auto d_prestencil_weights = my_GMLS_vector_of_scalar_clones.getPrestencilWeights();
    auto prestencil_weights = Kokkos::create_mirror_view(d_prestencil_weights);
    Kokkos::deep_copy(prestencil_weights, d_prestencil_weights);

    // tangent vector at target sites are stored here
    auto d_tangent_directions = *(my_GMLS_vector_of_scalar_clones.getTangentDirections());
    auto tangent_directions = Kokkos::create_mirror_view(d_tangent_directions);
    Kokkos::deep_copy(tangent_directions, d_tangent_directions);

    // loop through the target sites
    for (int i=0; i<number_target_coords; i++) {

        host_scratch_matrix_right_type T
                (tangent_directions.data() + TO_GLOBAL(i)*TO_GLOBAL(3)*TO_GLOBAL(3), 3, 3);
        XYZ u;

    
        // load value from output
        double GMLS_value = output_value(i);
        //printf("GMLS val: %f, %d\n", GMLS_value, i);
    
        double GMLS_gc = output_gaussian_curvature(i);
        //printf("GMLS gc: %f, %d\n", GMLS_gc, i);

    //    // load laplacian from output
    //    double GMLS_Laplacian = output_laplacian(i);
    
        // target site i's coordinate
        //double xval = target_coords(i,0);
        //double yval = (dimension>1) ? target_coords(i,1) : 0;
        //double zval = (dimension>2) ? target_coords(i,2) : 0;
        double xval = source_coords(neighbor_lists(i,1),0);
        double yval = (dimension>1) ? source_coords(neighbor_lists(i,1),1) : 0;
        double zval = (dimension>2) ? source_coords(neighbor_lists(i,1),2) : 0;
        double coord[3] = {xval, yval, zval};

        // get tangent vector and see if orthgonal to coordinate (it should be on a sphere)
        for (int j=0; j<dimension-1; ++j) {
            double tangent_inner_prod = 0;
            for (int k=0; k<std::min(dimension,3); ++k) {
                tangent_inner_prod += coord[k] * prestencil_weights(0, i, 0 /* local neighbor index */, j, k);
            }
            tangent_bundle_error += tangent_inner_prod * tangent_inner_prod;
        }
        double normal_inner_prod = 0;
        for (int k=0; k<dimension; ++k) {
            normal_inner_prod += coord[k] * my_GMLS_scalar.getTangentBundle(i, dimension-1, k);
        }
        // inner product could be plus or minus 1 (depends on tangent direction ordering)
        double normal_error_to_sum = (normal_inner_prod > 0) ? normal_inner_prod - 1 : normal_inner_prod + 1;
        tangent_bundle_error += normal_error_to_sum * normal_error_to_sum;
        tangent_bundle_norm += 1;
    
        // evaluation of various exact solutions
        double actual_value = sphere_harmonic54(xval, yval, zval);
        double actual_gc = 1.0; // Gaussian curvature is constant
    //    double actual_Laplacian = laplace_beltrami_sphere_harmonic54(xval, yval, zval);
        double actual_Gradient_ambient[3] = {0,0,0}; // initialized for 3, but only filled up to dimension
        gradient_sphereHarmonic54_ambient(actual_Gradient_ambient, xval, yval, zval);
    //    //velocity_sphereHarmonic54_ambient(actual_Gradient_ambient, xval, yval, zval);
        double actual_Curl_ambient[3] = {0,0,0};
        curl_sphere_harmonic54(actual_Curl_ambient, xval, yval, zval);

        values_error += (GMLS_value - actual_value)*(GMLS_value - actual_value);
        values_norm  += actual_value*actual_value;

        gc_error += (GMLS_gc - actual_gc)*(GMLS_gc - actual_gc);
        gc_norm  += 1;

        for (int j=0; j<dimension; ++j) u[j] = output_curl(i,j);
        AmbientLocalAmbient(u, T.data(), &prestencil_weights(0, i, 0 /* local neighbor index */, 0, 0));
        for (int j=0; j<dimension; ++j) output_curl(i,j) = u[j];

        for (int j=0; j<dimension; ++j) {
            curl_ambient_error += (output_curl(i,j) - actual_Curl_ambient[j])*(output_curl(i,j) - actual_Curl_ambient[j]);
            curl_ambient_norm += actual_Curl_ambient[j]*actual_Curl_ambient[j];
        }
    
    //    laplacian_error += (GMLS_Laplacian - actual_Laplacian)*(GMLS_Laplacian - actual_Laplacian);
    //    laplacian_norm += actual_Laplacian*actual_Laplacian;

    //    for (int j=0; j<dimension; ++j) {
    //        gradient_ambient_error += (output_gradient(i,j) - actual_Gradient_ambient[j])*(output_gradient(i,j) - actual_Gradient_ambient[j]);
    //        gradient_ambient_norm += actual_Gradient_ambient[j]*actual_Gradient_ambient[j];
    //    }

        for (int j=0; j<dimension; ++j) u[j] = output_vector(i,j);
        AmbientLocalAmbient(u, T.data(), &prestencil_weights(0, i, 0 /* local neighbor index */, 0, 0));
        for (int j=0; j<dimension; ++j) output_vector(i,j) = u[j];

        for (int j=0; j<dimension; ++j) {
            vector_ambient_error += (output_vector(i,j) - actual_Gradient_ambient[j])*(output_vector(i,j) - actual_Gradient_ambient[j]);
            vector_ambient_norm += actual_Gradient_ambient[j]*actual_Gradient_ambient[j];
        }

    //    divergence_ambient_error += (output_divergence(i) - actual_Laplacian)*(output_divergence(i) - actual_Laplacian);
    //    divergence_ambient_norm += actual_Laplacian*actual_Laplacian;
  
        
        for (int j=0; j<dimension; ++j) u[j] = output_vector_of_scalar_clones(i,j);
        AmbientLocalAmbient(u, T.data(), &prestencil_weights(0, i, 0 /* local neighbor index */, 0, 0));
        for (int j=0; j<dimension; ++j) output_vector_of_scalar_clones(i,j) = u[j];
        //// first get back to local chart
        //double local_vec[3] = {0,0};
        //for (int j=0; j<dimension; ++j) {
        //    local_vec[0] += T(0,j) * output_vector_of_scalar_clones(i,j);
        //    local_vec[1] += T(1,j) * output_vector_of_scalar_clones(i,j);
        //}
        //// second go to ambient space using tangent for first neighbor
        //for (int j=0; j<dimension; ++j) output_vector_of_scalar_clones(i,j) = 0;
        //for (int j=0; j<dimension; ++j) {
        //    output_vector_of_scalar_clones(i,j) += prestencil_weights(0, i, 0 /* local neighbor index */, 0, j) * local_vec[0];
        //    output_vector_of_scalar_clones(i,j) += prestencil_weights(0, i, 0 /* local neighbor index */, 1, j) * local_vec[1];
        //}


        for (int j=0; j<dimension; ++j) {
            vector_of_scalar_clones_ambient_error += (output_vector_of_scalar_clones(i,j) - actual_Gradient_ambient[j])*(output_vector_of_scalar_clones(i,j) - actual_Gradient_ambient[j]);
            vector_of_scalar_clones_ambient_norm += actual_Gradient_ambient[j]*actual_Gradient_ambient[j];
        }

    //    divergence_of_scalar_clones_ambient_error += (output_divergence_of_scalar_clones(i) - actual_Laplacian)*(output_divergence_of_scalar_clones(i) - actual_Laplacian);
    //    divergence_of_scalar_clones_ambient_norm += actual_Laplacian*actual_Laplacian;

    }

    tangent_bundle_error /= number_target_coords;
    tangent_bundle_error = std::sqrt(tangent_bundle_error);
    tangent_bundle_norm /= number_target_coords;
    tangent_bundle_norm = std::sqrt(tangent_bundle_norm);
    
    values_error /= number_target_coords;
    values_error = std::sqrt(values_error);
    values_norm /= number_target_coords;
    values_norm = std::sqrt(values_norm);

    gc_error /= number_target_coords;
    gc_error = std::sqrt(gc_error);
    gc_norm /= number_target_coords;
    gc_norm = std::sqrt(gc_norm);

    curl_ambient_error /= number_target_coords;
    curl_ambient_error = std::sqrt(curl_ambient_error);
    curl_ambient_norm /= number_target_coords;
    curl_ambient_norm = std::sqrt(curl_ambient_norm);

    //laplacian_error /= number_target_coords;
    //laplacian_error = std::sqrt(laplacian_error);
    //laplacian_norm /= number_target_coords;
    //laplacian_norm = std::sqrt(laplacian_norm);

    //gradient_ambient_error /= number_target_coords;
    //gradient_ambient_error = std::sqrt(gradient_ambient_error);
    //gradient_ambient_norm /= number_target_coords;
    //gradient_ambient_norm = std::sqrt(gradient_ambient_norm);
   
    vector_ambient_error /= number_target_coords;
    vector_ambient_error = std::sqrt(vector_ambient_error);
    vector_ambient_norm /= number_target_coords;
    vector_ambient_norm = std::sqrt(vector_ambient_norm);

    //divergence_ambient_error /= number_target_coords;
    //divergence_ambient_error = std::sqrt(divergence_ambient_error);
    //divergence_ambient_norm /= number_target_coords;
    //divergence_ambient_norm = std::sqrt(divergence_ambient_norm);

    vector_of_scalar_clones_ambient_error /= number_target_coords;
    vector_of_scalar_clones_ambient_error = std::sqrt(vector_of_scalar_clones_ambient_error);
    vector_of_scalar_clones_ambient_norm /= number_target_coords;
    vector_of_scalar_clones_ambient_norm = std::sqrt(vector_of_scalar_clones_ambient_norm);

    //divergence_of_scalar_clones_ambient_error /= number_target_coords;
    //divergence_of_scalar_clones_ambient_error = std::sqrt(divergence_of_scalar_clones_ambient_error);
    //divergence_of_scalar_clones_ambient_norm /= number_target_coords;
    //divergence_of_scalar_clones_ambient_norm = std::sqrt(divergence_of_scalar_clones_ambient_norm);

    printf("Tangent Bundle Error: %g\n", tangent_bundle_error / tangent_bundle_norm);  
    printf("Point Value Error: %g\n", values_error / values_norm);  
    printf("Gaussian Curvature Error: %g\n", gc_error / gc_norm);  
    printf("Surface Curl (Ambient) Error: %g\n", curl_ambient_error / curl_ambient_norm);  
    //printf("Laplace-Beltrami Error: %g\n", laplacian_error / laplacian_norm);  
    //printf("Surface Gradient (Ambient) Error: %g\n", gradient_ambient_error / gradient_ambient_norm);  
    printf("Surface Vector (VectorBasis) Error: %g\n", vector_ambient_error / vector_ambient_norm);  
    //printf("Surface Divergence (VectorBasis) Error: %g\n", divergence_ambient_error / divergence_ambient_norm);  
    printf("Surface Vector (ScalarClones) Error: %g\n", 
            vector_of_scalar_clones_ambient_error / vector_of_scalar_clones_ambient_norm);  
    //printf("Surface Divergence (ScalarClones) Error: %g\n", 
    //        divergence_of_scalar_clones_ambient_error / divergence_of_scalar_clones_ambient_norm);  
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
