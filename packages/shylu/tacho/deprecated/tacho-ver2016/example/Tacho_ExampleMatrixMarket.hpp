#ifndef __TACHO_EXAMPLE_MATRIX_MARKET_HPP__
#define __TACHO_EXAMPLE_MATRIX_MARKET_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "Tacho_Util.hpp"
#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_CrsMatrixTools.hpp"
#include "Tacho_MatrixMarket.hpp"

namespace Tacho {
  
  template<typename DeviceSpaceType>
  int exampleMatrixMarket(const std::string file_input,
                          const bool verbose) {
    typedef typename
      Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;

    const bool detail = false;
    std::cout << "DeviceSpace::  "; DeviceSpaceType::print_configuration(std::cout, detail);
    std::cout << "HostSpace::    ";   HostSpaceType::print_configuration(std::cout, detail);
  
    typedef CrsMatrixBase<value_type,ordinal_type,size_type,HostSpaceType>   CrsMatrixBaseHostType;

    int r_val = 0;

    Kokkos::Impl::Timer timer;
    
    CrsMatrixBaseHostType AA("AA");
    timer.reset();
    {
      std::ifstream in;
      in.open(file_input);
      if (!in.good()) {
        std::cout << "Failed in open the file: " << file_input << std::endl;
        return -1;
      }
      MatrixMarket::read(AA, in);
    }
    double t_read = timer.seconds();
    
    timer.reset();
    {
      std::string file_output = "mm-test-output.mtx";
      std::ofstream out;
      out.open(file_output);
      if (!out.good()) {
        std::cout << "Failed in open the file: " << file_output << std::endl;
        return -1;
      }
      MatrixMarket::write(out, AA, "%% Test output");
    }
    double t_write = timer.seconds();

    {
      const auto prec = std::cout.precision();
      std::cout.precision(4);

      std::cout << std::scientific
                << "MatrixMarket:: dimension = " << AA.NumRows() << " x " << AA.NumCols() 
                << ", " << " nnz = " << AA.NumNonZeros() << ", "
                << "read = " << t_read << " [sec], "
                << "write = " << t_write << " [sec] "
                << std::endl;

      std::cout.unsetf(std::ios::scientific);
      std::cout.precision(prec);
    }
    
    if (verbose) {
      AA.showMe(std::cout) << std::endl;
    }

    CrsMatrixBaseHostType BB("BB");
    BB.createConfTo(AA);

    CrsMatrixTools::copy(BB, Uplo::Upper, 0, AA);
    if (verbose) {
      BB.setLabel("Copy::AA:Upper::0"); BB.showMe(std::cout) << std::endl;
    }

    CrsMatrixTools::copy(BB, Uplo::Upper, 1, AA);
    if (verbose) {
      BB.setLabel("Copy::AA:Upper::1"); BB.showMe(std::cout) << std::endl;    
    }

    CrsMatrixTools::copy(BB, Uplo::Lower, 0, AA);
    if (verbose) {
      BB.setLabel("Copy::AA:Lower::0"); BB.showMe(std::cout) << std::endl;    
    }

    CrsMatrixTools::copy(BB, Uplo::Lower, 1, AA);
    if (verbose) {
      BB.setLabel("Copy::AA:Lower::1"); BB.showMe(std::cout) << std::endl;    
    }

    return r_val;
  }
}

#endif
