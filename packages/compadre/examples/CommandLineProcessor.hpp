// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef COMMANDLINEPROCESSOR
#define COMMANDLINEPROCESSOR

struct CommandLineProcessor {

    int order, curvature_order, dimension, number_target_coords, number_source_coords, number_of_batches;
    std::string constraint_name, solver_name, problem_name;

    Compadre::WeightingFunctionType kernel_type, curvature_kernel_type;

    CommandLineProcessor(int argc, char* args[], const bool print=true) {

        order = 2;
        curvature_order = -1; // defaults to max(order, 2)
        dimension = 3;
        number_target_coords = 200; 
        number_source_coords = -1; 
        number_of_batches = 1; 
    
        constraint_name = "NO_CONSTRAINT"; // "NEUMANN_GRAD_SCALAR"
        solver_name = "QR"; // or CHOL
        problem_name = "STANDARD"; // MANIFOLD
        kernel_type = Compadre::GMLS::translateWeightingType("POWER");
        curvature_kernel_type = Compadre::GMLS::translateWeightingType("POWER");

        for (int i = 1; i < argc; ++i) {
            if (i + 1 < argc) { // not at end
                if (std::string(args[i]) == "--p") {
                   order = atoi(args[i+1]); 
                } else if (std::string(args[i]) == "--curvature-order") {
                   curvature_order = atoi(args[i+1]);
                } else if (std::string(args[i]) == "--d") {
                   dimension = atoi(args[i+1]); 
                } else if (std::string(args[i]) == "--nt") {
                   number_target_coords = atoi(args[i+1]); 
                } else if (std::string(args[i]) == "--ns") {
                   number_source_coords = atoi(args[i+1]); 
                } else if (std::string(args[i]) == "--nb") {
                   number_of_batches = atoi(args[i+1]); 
                } else if (std::string(args[i]) == "--solver") {
                   solver_name = std::string(args[i+1]); 
                } else if (std::string(args[i]) == "--problem") {
                   problem_name = std::string(args[i+1]); 
                } else if (std::string(args[i]) == "--constraint") {
                   constraint_name = std::string(args[i+1]); 
                } else if (std::string(args[i]) == "--kernel-type") {
                   kernel_type = Compadre::GMLS::translateWeightingType(std::string(args[i+1])); 
                } else if (std::string(args[i]) == "--curvature-kernel-type") {
                   curvature_kernel_type = Compadre::GMLS::translateWeightingType(std::string(args[i+1])); 
                }
            }
        }

        // default curvature order is at least 2: a linear (order 1) curvature fit gives a tangent plane
        // accurate only to O(h), with an error that varies irregularly with the neighbor arrangement
        if (curvature_order < 0) curvature_order = (order > 2) ? order : 2;

        if (print) {
            std::cout << "************************************************************************" << std::endl;
            std::cout << "order: " << order << ", curvature order: " << curvature_order << ", dimension: " << dimension << ", number_target_coords: " << number_target_coords << std::endl;
            std::cout << "solver: " << solver_name << ", problem: " << problem_name << ", constraint: " << constraint_name << std::endl;
            std::cout << "************************************************************************" << std::endl;

        }
    }
};

#endif
