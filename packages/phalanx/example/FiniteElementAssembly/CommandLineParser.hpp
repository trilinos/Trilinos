// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHALANX_COMMAND_LINE_PARSER_HPP
#define PHALANX_COMMAND_LINE_PARSER_HPP

#include <cstdlib>
#include <cstring>
#include "Teuchos_Assert.hpp"
#include "Kokkos_Core.hpp"

namespace phx_example {

  class CommandLineParser {
    int nx_;
    int ny_;
    int nz_;
    double lx_;
    double ly_;
    double lz_;
    int num_equations_;
    int workset_size_;
    bool do_residual_;
    bool do_jacobian_;
    bool print_residual_;
    bool print_jacobian_;
    bool print_to_file_;
    int team_size_;
    int vector_size_;
    bool do_graph_analysis_;

  public:

    CommandLineParser(int argc, char *argv[]) :
      nx_(4),
      ny_(2),
      nz_(2),
      lx_(1.0),
      ly_(1.0),
      lz_(1.0),
      num_equations_(2),
      workset_size_(3), // (num_elements % workset_size != 0) test partial workset
      do_residual_(true),
      do_jacobian_(true),
      print_residual_(false),
      print_jacobian_(false),
      print_to_file_(false),
      team_size_(1),
      vector_size_(1),
      do_graph_analysis_(true)
    {
      // Set better defaults for team/vector size based on arcitecture
#if defined(KOKKOS_ENABLE_CUDA)
      vector_size_ = 32;
      team_size_ = 256 / vector_size_;
#endif

      for (int i=1; i < argc; ++i) {
        // TEUCHOS_TEST_FOR_EXCEPTION(i+1 == argc,std::runtime_error,
        //                            "ERROR: the final value for the input parameter is missing!");
        if (std::strncmp(argv[i],"-nx",3)==0)
          nx_ = std::atoi(argv[++i]);
        else if (std::strncmp(argv[i],"-ny",3)==0)
          ny_ = std::atoi(argv[++i]);
        else if (std::strncmp(argv[i],"-nz",3)==0)
          nz_ = std::atoi(argv[++i]);
        else if (std::strncmp(argv[i],"-lx",3)==0)
          lx_ = std::atof(argv[++i]);
        else if (std::strncmp(argv[i],"-ly",3)==0)
          ly_ = std::atof(argv[++i]);
        else if (std::strncmp(argv[i],"-lz",3)==0)
          lz_ = std::atof(argv[++i]);
        else if (std::strncmp(argv[i],"-ne",3)==0)
          num_equations_ = std::atoi(argv[++i]);
        else if (std::strncmp(argv[i],"-ws",3)==0)
          workset_size_ = std::atoi(argv[++i]);
        else if (std::strncmp(argv[i],"-nr",3)==0)
          do_residual_ = false;
        else if (std::strncmp(argv[i],"-nj",3)==0)
          do_jacobian_ = false;
        else if (std::strncmp(argv[i],"-pr",3)==0)
          print_residual_ = true;
        else if (std::strncmp(argv[i],"-pj",3)==0)
          print_jacobian_ = true;
        else if (std::strncmp(argv[i],"-pa",3)==0) {
          print_residual_ = true;
	  print_jacobian_ = true;
	}
        else if (std::strncmp(argv[i],"-f",2)==0)
          print_to_file_ = true;
        else if (std::strncmp(argv[i],"-ts",3)==0)
          team_size_ = std::atoi(argv[++i]);
        else if (std::strncmp(argv[i],"-vs",3)==0)
          vector_size_ = std::atoi(argv[++i]);
        else if (std::strncmp(argv[i],"-nga",4)==0)
          do_graph_analysis_ = false;
        else if (std::strncmp(argv[i],"-h",2)==0) {
	  std::cout << "============================================" << std::endl;
	  std::cout << "usage: " << argv[0] << " [options]" << std::endl;
	  std::cout << "============================================" << std::endl;
	  std::cout << "-h           Print help message" << std::endl;
	  std::cout << "-nx <val>    Number of elements in the x direction [4]" << std::endl;
	  std::cout << "-ny <val>    Number of elements in the y direction [2]" << std::endl;
	  std::cout << "-nz <val>    Number of elements in the z direction [2]" << std::endl;
	  std::cout << "-lx <val>    Length of domain in the x direction [1.0]" << std::endl;
	  std::cout << "-ly <val>    Length of domain in the y direction [1.0]" << std::endl;
	  std::cout << "-lz <val>    Length of domain in the z direction [1.0]" << std::endl;
	  std::cout << "-ne <val>    Number of PDE equations per node [2]" << std::endl;
	  std::cout << "-ws <val>    Workset size [3]" << std::endl;
	  std::cout << "-nr          Disable the Residual computation" << std::endl;
	  std::cout << "-nj          Disable the Jacobian computation" << std::endl;
	  std::cout << "-pr          Print the Residual values" << std::endl;
	  std::cout << "-pj          Print the Jacobian values (also prints Residual values\n"
		    << "             that were evaluated during Jacobian computation)" << std::endl;
	  std::cout << "-pa          Print all evalaution type values" << std::endl;
	  std::cout << "-f           Print values to file instead of the screen" << std::endl;
	  std::cout << "-ts          Set the team size for Kokkos kernel launch" << std::endl;
	  std::cout << "-vs          Set the vector size for Kokkos kernel launch" << std::endl;
	  std::cout << "-nga         Disable the graph analysis" << std::endl;
	  std::cout << "============================================" << std::endl;
	  exit(0);
	}
        else {
          TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                                     "ERROR: the input parameter \"" << argv[i] << "\" is not valid!");
        }
      }
    }
    const int& nx() const {return nx_;}
    const int& ny() const {return ny_;}
    const int& nz() const {return nz_;}
    const double& lx() const {return lx_;}
    const double& ly() const {return ly_;}
    const double& lz() const {return lz_;}
    const int& numEquations() const {return num_equations_;}
    const int& worksetSize() const {return workset_size_;}
    const bool& doResidual() const {return do_residual_;}
    const bool& doJacobian() const {return do_jacobian_;}
    const bool& printResidual() const {return print_residual_;}
    const bool& printJacobian() const {return print_jacobian_;}
    const bool& printToFile() const {return print_to_file_;}
    const int& teamSize() const {return team_size_;}
    const int& vectorSize() const {return vector_size_;}
    const bool& doGraphAnalysis() const {return do_graph_analysis_;}
  };
}
#endif
