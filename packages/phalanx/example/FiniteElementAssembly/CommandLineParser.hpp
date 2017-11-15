//@HEADER
//@HEADER
#ifndef PHALANX_COMMAND_LINE_PARSER_HPP
#define PHALANX_COMMAND_LINE_PARSER_HPP

#include <cstdlib>
#include <cstring>
#include "Teuchos_Assert.hpp"

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
    bool print_residual_;
    bool print_jacobian_;
    bool print_to_file_;

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
      print_residual_(false),
      print_jacobian_(false),
      print_to_file_(false)
    {
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
    const bool& printResidual() const {return print_residual_;}
    const bool& printJacobian() const {return print_jacobian_;}
    const bool& printToFile() const {return print_to_file_;}
  };
}
#endif
