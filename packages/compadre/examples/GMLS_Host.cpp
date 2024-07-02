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
#include "GMLS_Tutorial.hpp"
#include "CommandLineProcessor.hpp"

#ifdef COMPADRE_USE_MPI
#include <mpi.h>
#endif

#include <Kokkos_Timer.hpp>
#include <Kokkos_Core.hpp>

using namespace Compadre;

int main (int argc, char* args[])
{

#ifdef COMPADRE_USE_MPI
    MPI_Init(&argc, &args);
#endif
    
bool all_passed = true;

{

    CommandLineProcessor clp(argc, args);
    auto order = clp.order;
    auto dimension = clp.dimension;
    auto number_target_coords = clp.number_target_coords;
    auto constraint_name = clp.constraint_name;
    auto solver_name = clp.solver_name;
    auto problem_name = clp.problem_name;

    const double failure_tolerance = 1e-9;

    const int offset = 15;
    std::mt19937 rng(50);
    const int min_neighbors = 1*Compadre::GMLS::getNP(order);
    const int max_neighbors = 1*Compadre::GMLS::getNP(order)*1.15;
    std::cout << min_neighbors << " " << max_neighbors << std::endl;
    std::uniform_int_distribution<int> gen_num_neighbors(min_neighbors, max_neighbors); // uniform, unbiased


    Kokkos::initialize(argc, args);
    Kokkos::Timer timer;
    Kokkos::Profiling::pushRegion("Setup");


    const int N = 40000;
    std::uniform_int_distribution<int> gen_neighbor_number(offset, N); // 0 to 10 are junk (part of test)


    Kokkos::View<int**, Kokkos::HostSpace>    neighbor_lists("neighbor lists", number_target_coords, max_neighbors+1); // first column is # of neighbors
    Kokkos::View<double**, Kokkos::HostSpace> source_coords("neighbor coordinates", N, dimension);
    Kokkos::View<double*, Kokkos::HostSpace> epsilon("h supports", number_target_coords);

    for (int i=0; i<number_target_coords; i++) {
        epsilon(i) = 0.5;
    }

//    // fake coordinates not to be used
    for(int i = 0; i < offset; i++){
        for(int j = 0; j < dimension; j++){
            source_coords(i,j) = 0.1;
        }
    }

    // filling others with random coordinates
    for(int i = offset; i < N; i++){ //ignore first ten entries
        double randx = (2.0*(double)rand() / (double) RAND_MAX - 1.0)*epsilon(0)/2.0;
        double randy = (2.0*(double)rand() / (double) RAND_MAX - 1.0)*epsilon(0)/2.0;
        double randz = (2.0*(double)rand() / (double) RAND_MAX - 1.0)*epsilon(0)/2.0;
        source_coords(i,0) = randx;
        if (dimension>1) source_coords(i,1) = randy;
        if (dimension>2) source_coords(i,2) = randz;
    }

    const double target_epsilon = 0.1;
    // fill target coords
    Kokkos::View<double**, Kokkos::HostSpace> target_coords ("target coordinates", number_target_coords, dimension);
    for(int i = 0; i < number_target_coords; i++){ //ignore first ten entries
        double randx = (2.0*(double)rand() / (double) RAND_MAX - 1.0)*target_epsilon/2.0;
        double randy = (2.0*(double)rand() / (double) RAND_MAX - 1.0)*target_epsilon/2.0;
        double randz = (2.0*(double)rand() / (double) RAND_MAX - 1.0)*target_epsilon/2.0;
        target_coords(i,0) = randx;
        if (dimension>1) target_coords(i,1) = randy;
        if (dimension>2) target_coords(i,2) = randz;
    }

    // randomly fill neighbor lists
    for (int i=0; i<number_target_coords; i++) {
//        int r = gen_num_neighbors(rng);
//        assert(r<source_coords.extent(0)-offset);
        int r = max_neighbors;
        neighbor_lists(i,0) = r; // number of neighbors is random between max and min

        for(int j=0; j<r; j++){
            neighbor_lists(i,j+1) = offset + j + 1;
//            bool placed = false;
//            while (!placed) {
//                int ind = gen_neighbor_number(rng);
//                bool located = false;
//                for (int k=1; k<j+1; k++) {
//                    if (neighbor_lists(i,k) == ind) {
//                        located = true;
//                        break;
//                    }
//                }
//                if (!located) {
//                    neighbor_lists(i,j+1) = ind;
//                    placed = true;
//                } // neighbors can be from  10,...,N-1
//            }
        }
    }

    Kokkos::Profiling::popRegion();
    timer.reset();

    GMLS my_GMLS(order, dimension,
                 solver_name.c_str(), problem_name.c_str(), constraint_name.c_str(),
                 2 /*manifold order*/);
    my_GMLS.setProblemData(neighbor_lists, source_coords, target_coords, epsilon);
    my_GMLS.setWeightingParameter(10);

    std::vector<TargetOperation> lro(3);
    lro[0] = ScalarPointEvaluation;
    lro[1] = LaplacianOfScalarPointEvaluation;
    lro[2] = GradientOfScalarPointEvaluation;
    my_GMLS.addTargets(lro);
    // add two more targets individually to test addTargets(...) function
    my_GMLS.addTargets(DivergenceOfVectorPointEvaluation);
    my_GMLS.addTargets(CurlOfVectorPointEvaluation);
    my_GMLS.generateAlphas();

    double instantiation_time = timer.seconds();
    std::cout << "Took " << instantiation_time << "s to complete instantiation." << std::endl;

    Kokkos::Profiling::pushRegion("Creating Data");

    
    // need Kokkos View storing true solution
    Kokkos::View<double*, Kokkos::HostSpace> sampling_data("samples of true solution", source_coords.extent(0));
    Kokkos::View<double**, Kokkos::HostSpace> gradient_sampling_data("samples of true gradient", source_coords.extent(0), dimension);
    Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::HostSpace> divergence_sampling_data("samples of true solution for divergence test", source_coords.extent(0), dimension);
    Kokkos::parallel_for("Sampling Manufactured Solutions", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0,source_coords.extent(0)), KOKKOS_LAMBDA(const int i) {
        double xval = source_coords(i,0);
        double yval = (dimension>1) ? source_coords(i,1) : 0;
        double zval = (dimension>2) ? source_coords(i,2) : 0;
        sampling_data(i) = trueSolution(xval, yval, zval, order, dimension);
        double true_grad[3] = {0,0,0};
        trueGradient(true_grad, xval, yval,zval, order, dimension);
        for (int j=0; j<dimension; ++j) {
            divergence_sampling_data(i,j) = divergenceTestSamples(xval, yval, zval, j, dimension);
            gradient_sampling_data(i,j) = true_grad[j];
        }
    });
    Kokkos::Profiling::popRegion();

    Evaluator gmls_evaluator(&my_GMLS);

    for (int i=0; i<number_target_coords; i++) {

        Kokkos::Profiling::pushRegion("Apply Alphas to Data");

        double GMLS_value = gmls_evaluator.applyAlphasToDataSingleComponentSingleTargetSite(sampling_data, 0, ScalarPointEvaluation, i, 0, 0, 0, 0, 0);
        //for (int j = 0; j< neighbor_lists(i,0); j++){
        //    double xval = source_coords(neighbor_lists(i,j+1),0);
        //    double yval = (dimension>1) ? source_coords(neighbor_lists(i,j+1),1) : 0;
        //    double zval = (dimension>2) ? source_coords(neighbor_lists(i,j+1),2) : 0;
        //    GMLS_value += gmls_evaluator.getAlpha0TensorTo0Tensor(ScalarPointEvaluation, i, j)*trueSolution(xval, yval, zval, order, dimension);
        //}

        double GMLS_Laplacian = gmls_evaluator.applyAlphasToDataSingleComponentSingleTargetSite(sampling_data, 0, LaplacianOfScalarPointEvaluation, i, 0, 0, 0, 0, 0);
        //double GMLS_Laplacian = 0.0;
        //for (int j = 0; j< neighbor_lists(i,0); j++){
        //    double xval = source_coords(neighbor_lists(i,j+1),0);
        //    double yval = (dimension>1) ? source_coords(neighbor_lists(i,j+1),1) : 0;
        //    double zval = (dimension>2) ? source_coords(neighbor_lists(i,j+1),2) : 0;
        //    GMLS_Laplacian += gmls_evaluator.getAlpha0TensorTo0Tensor(LaplacianOfScalarPointEvaluation, i, j)*trueSolution(xval, yval, zval, order, dimension);
        //}

        double GMLS_GradX = gmls_evaluator.applyAlphasToDataSingleComponentSingleTargetSite(sampling_data, 0, GradientOfScalarPointEvaluation, i, 0, 0, 0, 0, 0);
        //double GMLS_GradX = 0.0;
        //for (int j = 0; j< neighbor_lists(i,0); j++){
        //    double xval = source_coords(neighbor_lists(i,j+1),0);
        //    double yval = (dimension>1) ? source_coords(neighbor_lists(i,j+1),1) : 0;
        //    double zval = (dimension>2) ? source_coords(neighbor_lists(i,j+1),2) : 0;
        //    GMLS_GradX += gmls_evaluator.getAlpha0TensorTo1Tensor(GradientOfScalarPointEvaluation, i, 0, j)*trueSolution(xval, yval, zval, order, dimension);
        //}

        double GMLS_GradY = (dimension>1) ? gmls_evaluator.applyAlphasToDataSingleComponentSingleTargetSite(sampling_data, 0, GradientOfScalarPointEvaluation, i, 0, 1, 0, 0, 0) : 0;
        //double GMLS_GradY = 0.0;
        //if (dimension>1) {
        //    for (int j = 0; j< neighbor_lists(i,0); j++){
        //        double xval = source_coords(neighbor_lists(i,j+1),0);
        //        double yval = source_coords(neighbor_lists(i,j+1),1);
        //        double zval = (dimension>2) ? source_coords(neighbor_lists(i,j+1),2) : 0;
        //        GMLS_GradY += gmls_evaluator.getAlpha0TensorTo1Tensor(GradientOfScalarPointEvaluation, i, 1, j)*trueSolution(xval, yval, zval, order, dimension);
        //    }
        //}

        double GMLS_GradZ = (dimension>2) ? gmls_evaluator.applyAlphasToDataSingleComponentSingleTargetSite(sampling_data, 0, GradientOfScalarPointEvaluation, i, 0, 2, 0, 0, 0) : 0;
        //double GMLS_GradZ = 0.0;
        //if (dimension>2) {
        //    for (int j = 0; j< neighbor_lists(i,0); j++){
        //        double xval = source_coords(neighbor_lists(i,j+1),0);
        //        double yval = source_coords(neighbor_lists(i,j+1),1);
        //        double zval = source_coords(neighbor_lists(i,j+1),2);
        //        GMLS_GradZ += gmls_evaluator.getAlpha0TensorTo1Tensor(GradientOfScalarPointEvaluation, i, 2, j)*trueSolution(xval, yval, zval, order, dimension);
        //    }
        //}

        double GMLS_Divergence  = gmls_evaluator.applyAlphasToDataSingleComponentSingleTargetSite(gradient_sampling_data, 0, DivergenceOfVectorPointEvaluation, i, 0, 0, 0, 0, 0);
        if (dimension>1) GMLS_Divergence += gmls_evaluator.applyAlphasToDataSingleComponentSingleTargetSite(gradient_sampling_data, 1, DivergenceOfVectorPointEvaluation, i, 0, 0, 0, 1, 0);
        if (dimension>2) GMLS_Divergence += gmls_evaluator.applyAlphasToDataSingleComponentSingleTargetSite(gradient_sampling_data, 2, DivergenceOfVectorPointEvaluation, i, 0, 0, 0, 2, 0);

        //double GMLS_Divergence  = gmls_evaluator.applyAlphasToDataSingleComponentSingleTargetSite(Kokkos::subview(gradient_sampling_data,Kokkos::ALL,0), 0, DivergenceOfVectorPointEvaluation, i, 0, 0, 0, 0);
        //if (dimension>1) GMLS_Divergence += gmls_evaluator.applyAlphasToDataSingleComponentSingleTargetSite(Kokkos::subview(gradient_sampling_data,Kokkos::ALL,1), 1, DivergenceOfVectorPointEvaluation, i, 0, 0, 1, 0);
        //if (dimension>2) GMLS_Divergence += gmls_evaluator.applyAlphasToDataSingleComponentSingleTargetSite(Kokkos::subview(gradient_sampling_data,Kokkos::ALL,2), 2, DivergenceOfVectorPointEvaluation, i, 0, 0, 2, 0);
        //double GMLS_Divergence = 0.0;
        //for (int j = 0; j< neighbor_lists(i,0); j++){
        //    double xval = source_coords(neighbor_lists(i,j+1),0);
        //    double yval = (dimension>1) ? source_coords(neighbor_lists(i,j+1),1) : 0;
        //    double zval = (dimension>2) ? source_coords(neighbor_lists(i,j+1),2) : 0;
        //    // TODO: use different functions for the vector components
        //    if (use_arbitrary_order_divergence) {
        //        GMLS_Divergence += gmls_evaluator.getAlpha1TensorTo0Tensor(DivergenceOfVectorPointEvaluation, i, j, 0)*trueSolution(xval, yval, zval, order, dimension);
        //        if (dimension>1) GMLS_Divergence += gmls_evaluator.getAlpha1TensorTo0Tensor(DivergenceOfVectorPointEvaluation, i, j, 1)*trueSolution(xval, yval, zval, order, dimension);
        //        if (dimension>2) GMLS_Divergence += gmls_evaluator.getAlpha1TensorTo0Tensor(DivergenceOfVectorPointEvaluation, i, j, 2)*trueSolution(xval, yval, zval, order, dimension);
        //    } else {
        //        for (int k=0; k<dimension; ++k) {
        //            GMLS_Divergence += gmls_evaluator.getAlpha1TensorTo0Tensor(DivergenceOfVectorPointEvaluation, i, j, k)*divergenceTestSamples(xval, yval, zval, k, dimension);
        //        }
        //    }
        //}

        double GMLS_CurlX = 0.0;
        double GMLS_CurlY = 0.0;
        double GMLS_CurlZ = 0.0;
        if (dimension>1) {
            for (int j=0; j<dimension; ++j) {
                GMLS_CurlX += gmls_evaluator.applyAlphasToDataSingleComponentSingleTargetSite(divergence_sampling_data, j, CurlOfVectorPointEvaluation, i, 0, 0, 0, j, 0);
                GMLS_CurlY += gmls_evaluator.applyAlphasToDataSingleComponentSingleTargetSite(divergence_sampling_data, j, CurlOfVectorPointEvaluation, i, 0, 1, 0, j, 0);
            }
        } 

        if (dimension>2) {
            for (int j=0; j<dimension; ++j) {
                GMLS_CurlZ += gmls_evaluator.applyAlphasToDataSingleComponentSingleTargetSite(divergence_sampling_data, j, CurlOfVectorPointEvaluation, i, 0, 2, 0, j, 0);
            }

        }

        Kokkos::Profiling::popRegion();
        //if (dimension>1) {
        //    for (int j = 0; j< neighbor_lists(i,0); j++){
        //        double xval = source_coords(neighbor_lists(i,j+1),0);
        //        double yval = source_coords(neighbor_lists(i,j+1),1);
        //        double zval = (dimension>2) ? source_coords(neighbor_lists(i,j+1),2) : 0;
        //        for (int k=0; k<dimension; ++k) {
        //            GMLS_CurlX += my_GMLS.getAlpha1TensorTo1Tensor(CurlOfVectorPointEvaluation, i, 0, j, k)*divergenceTestSamples(xval, yval, zval, k, dimension);
        //        }
        //    }

        //    for (int j = 0; j< neighbor_lists(i,0); j++){
        //        double xval = source_coords(neighbor_lists(i,j+1),0);
        //        double yval = source_coords(neighbor_lists(i,j+1),1);
        //        double zval = (dimension>2) ? source_coords(neighbor_lists(i,j+1),2) : 0;
        //        for (int k=0; k<dimension; ++k) {
        //            GMLS_CurlY += my_GMLS.getAlpha1TensorTo1Tensor(CurlOfVectorPointEvaluation, i, 1, j, k)*divergenceTestSamples(xval, yval, zval, k, dimension);
        //        }
        //    }
        //}

        //if (dimension>2) {
        //    for (int j = 0; j< neighbor_lists(i,0); j++){
        //        double xval = source_coords(neighbor_lists(i,j+1),0);
        //        double yval = source_coords(neighbor_lists(i,j+1),1);
        //        double zval = source_coords(neighbor_lists(i,j+1),2);
        //        for (int k=0; k<dimension; ++k) {
        //            GMLS_CurlZ += my_GMLS.getAlpha1TensorTo1Tensor(CurlOfVectorPointEvaluation, i, 2, j, k)*divergenceTestSamples(xval, yval, zval, k, dimension);
        //        }
        //    }
        //}
        //
        Kokkos::Profiling::pushRegion("Comparison");

        double xval = target_coords(i,0);
        double yval = (dimension>1) ? target_coords(i,1) : 0;
        double zval = (dimension>2) ? target_coords(i,2) : 0;

        double actual_value = trueSolution(xval, yval, zval, order, dimension);
        double actual_Laplacian = trueLaplacian(xval, yval, zval, order, dimension);
        double actual_Gradient[3] = {0,0,0};
        trueGradient(actual_Gradient, xval, yval, zval, order, dimension);
        double actual_Divergence;
        actual_Divergence = trueLaplacian(xval, yval, zval, order, dimension);

        double actual_CurlX = 0;
        double actual_CurlY = 0;
        double actual_CurlZ = 0;
        if (dimension>1) {
            actual_CurlX = curlTestSolution(xval, yval, zval, 0, dimension);
            actual_CurlY = curlTestSolution(xval, yval, zval, 1, dimension);
        }
        if (dimension>2) {
            actual_CurlZ = curlTestSolution(xval, yval, zval, 2, dimension);
        }

//        fprintf(stdout, "Reconstructed value: %f \n", GMLS_value);
//        fprintf(stdout, "Actual value: %f \n", actual_value);
//        fprintf(stdout, "Reconstructed Laplacian: %f \n", GMLS_Laplacian);
//        fprintf(stdout, "Actual Laplacian: %f \n", actual_Laplacian);

        if(GMLS_value!=GMLS_value || std::abs(actual_value - GMLS_value) > failure_tolerance) {
            all_passed = false;
            std::cout << "Failed Actual by: " << std::abs(actual_value - GMLS_value) << std::endl;
        }

        if(std::abs(actual_Laplacian - GMLS_Laplacian) > failure_tolerance) {
            all_passed = false;
            std::cout << "Failed Laplacian by: " << std::abs(actual_Laplacian - GMLS_Laplacian) << std::endl;
        }

        if(std::abs(actual_Gradient[0] - GMLS_GradX) > failure_tolerance) {
            all_passed = false;
            std::cout << "Failed GradX by: " << std::abs(actual_Gradient[0] - GMLS_GradX) << std::endl;
        }

        if (dimension>1) {
            if(std::abs(actual_Gradient[1] - GMLS_GradY) > failure_tolerance) {
                all_passed = false;
                std::cout << "Failed GradY by: " << std::abs(actual_Gradient[1] - GMLS_GradY) << std::endl;
            }
        }

        if (dimension>2) {
            if(std::abs(actual_Gradient[2] - GMLS_GradZ) > failure_tolerance) {
                all_passed = false;
                std::cout << "Failed GradZ by: " << std::abs(actual_Gradient[2] - GMLS_GradZ) << std::endl;
            }
        }

        if(std::abs(actual_Divergence - GMLS_Divergence) > failure_tolerance) {
            all_passed = false;
            std::cout << "Failed Divergence by: " << std::abs(actual_Divergence - GMLS_Divergence) << std::endl;
        }

        double tmp_diff = 0;
        if (dimension>1)
            tmp_diff += std::abs(actual_CurlX - GMLS_CurlX) + std::abs(actual_CurlY - GMLS_CurlY);
        if (dimension>2)
            tmp_diff += std::abs(actual_CurlZ - GMLS_CurlZ);
        if(std::abs(tmp_diff) > failure_tolerance) {
            all_passed = false;
            std::cout << "Failed Curl by: " << std::abs(tmp_diff) << std::endl;
        }
        Kokkos::Profiling::popRegion();
    }

}

    Kokkos::finalize();
#ifdef COMPADRE_USE_MPI
    MPI_Finalize();
#endif

if(all_passed) {
    fprintf(stdout, "Passed test \n");
    return 0;
} else {
    fprintf(stdout, "Failed test \n");
    return -1;
}
}
