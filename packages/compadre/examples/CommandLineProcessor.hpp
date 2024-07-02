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

    int order, dimension, number_target_coords, number_source_coords, number_of_batches;
    std::string constraint_name, solver_name, problem_name;

    CommandLineProcessor(int argc, char* args[], const bool print=true) {

        order = 2;
        dimension = 3;
        number_target_coords = 200; 
        number_source_coords = -1; 
        number_of_batches = 1; 
    
        constraint_name = "NO_CONSTRAINT"; // "NEUMANN_GRAD_SCALAR"
        solver_name = "QR"; // LU
        problem_name = "STANDARD"; // MANIFOLD

        for (int i = 1; i < argc; ++i) {
            if (i + 1 < argc) { // not at end
                if (std::string(args[i]) == "--p") {
                   order = atoi(args[i+1]); 
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
                }
            }
        }

        if (print) {
            std::cout << "************************************************************************" << std::endl;
            std::cout << "order: " << order << ", dimension: " << dimension << ", number_target_coords: " << number_target_coords << std::endl;
            std::cout << "solver: " << solver_name << ", problem: " << problem_name << ", constraint: " << constraint_name << std::endl;
            std::cout << "************************************************************************" << std::endl;

        }
    }
};

#endif
