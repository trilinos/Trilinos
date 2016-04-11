#ifndef __TACHO_EXAMPLE_SOLVER_HPP__
#define __TACHO_EXAMPLE_SOLVER_HPP__

#include "Tacho_Solver.hpp"

#ifdef HAVE_SHYLUTACHO_VTUNE
#include "ittnotify.h"
#endif

namespace Tacho {

#define TACHO_SOLVER_RUN(run, time) timer.reset(); run; const double time = timer.seconds();

  template<typename DeviceSpaceType>
  int exampleSolver(const std::string file_input,
                    const int prunecut,
                    const int max_concurrency,
                    const int nrhs,
                    const int mb,
                    const int nb,
                    const bool verbose) {
    typedef typename
      Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;
    
    const bool detail = false;
    std::cout << "DeviceSpace::  "; DeviceSpaceType::print_configuration(std::cout, detail);
    std::cout << "HostSpace::    ";   HostSpaceType::print_configuration(std::cout, detail);

    typedef Solver<value_type,ordinal_type,size_type,DeviceSpaceType> SolverType;

#ifdef HAVE_SHYLUTACHO_VTUNE
    __itt_pause();
#endif

    int r_val = 0;
    Kokkos::Impl::Timer timer;

    SolverType tacho("Tacho::CholSolver");
    tacho.setPolicy(max_concurrency);
    tacho.setBlocksize(mb, nb);

    ///
    /// Read from matrix market
    ///
    ///     input  - file
    ///     output - AA
    ///
    typename SolverType::CrsMatrixBaseHostType AA("AA");
    timer.reset();
    {
      std::ifstream in;
      in.open(file_input);
      if (!in.good()) {
        std::cout << "Failed in open the file: " << file_input << std::endl;
        return -1;
      }
      MatrixMarket::read(AA, in);
      
      typename SolverType::DenseMatrixBaseHostType BB("BB", AA.NumRows(), nrhs), XX("XX");
      XX.createConfTo(BB);
      
      const ordinal_type m = BB.NumRows();
      for (auto rhs=0;rhs<nrhs;++rhs) 
        for (auto i=0;i<m;++i) 
          BB.Value(i, rhs) = (rhs + 1);
      
      tacho.setProblem(AA, BB, XX);
    }
    const double t_input = timer.seconds();    
    std::cout << "Solver:: "
              << "input generation = " << t_input << " [sec] "
              << std::endl;

    if (verbose) {
      tacho.getA().showMe(std::cout) << std::endl;
      tacho.getB().showMe(std::cout) << std::endl;
      tacho.getX().showMe(std::cout) << std::endl;
    }


    ///
    /// Solver interface
    ///
    TACHO_SOLVER_RUN(tacho.reorder(prunecut), t_reorder);
    TACHO_SOLVER_RUN(tacho.analyze(), t_analyze);

#ifdef HAVE_SHYLUTACHO_VTUNE
    __itt_resume();
#endif
    TACHO_SOLVER_RUN(tacho.factorize(), t_factorize);
#ifdef HAVE_SHYLUTACHO_VTUNE
    __itt_pause();
#endif

    TACHO_SOLVER_RUN(tacho.solve(), t_solve);

    double norm, error;
    tacho.check(norm, error);
  
    ///
    /// Print out 
    ///
    {
      const auto prec = std::cout.precision();
      std::cout.precision(4);
      
      const auto AA = tacho.getA();
      const auto UU = tacho.getFlatU();
      const auto HU = tacho.getHierU();

      std::cout << std::scientific;
      std::cout << "Solver:: Given    matrix = " << AA.NumRows() << " x " << AA.NumCols() << ", nnz = " << AA.NumNonZeros() << std::endl;
      std::cout << "Solver:: Factored matrix = " << UU.NumRows() << " x " << UU.NumCols() << ", nnz = " << UU.NumNonZeros() << std::endl;
      std::cout << "Solver:: Hier     matrix = " << HU.NumRows() << " x " << HU.NumCols() << ", nnz = " << HU.NumNonZeros() << std::endl;
      
      std::cout << "Solver:: "
                << "input generation = " << t_input << " [sec], "
                << std::endl
                << "Solver:: "
                << "reorder          = " << t_reorder << " [sec] "
                << std::endl
                << "Solver:: "
                << "analyze          = " << t_analyze << " [sec] "
                << std::endl
                << "Solver:: "
                << "factorize        = " << t_factorize << " [sec] "
                << std::endl
                << "Solver:: "
                << "solve            = " << t_solve << " [sec] "
                << std::endl
                << "Solver:: "
                << "norm = " << norm << ", error = " << error << ", rel error = " << (error/norm) 
                << std::endl;

      std::cout << std::endl;
      
      std::cout.unsetf(std::ios::scientific);
      std::cout.precision(prec);
    }
    
    return r_val;
  }
}

#endif
