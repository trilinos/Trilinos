// @HEADER
// @HEADER

#ifndef NOX_SOLVER_STATS_HPP
#define NOX_SOLVER_STATS_HPP

#include "NOX_Common.H"
#include "NOX_LineSearch_Utils_Counters.H"

namespace NOX {

  //! Container for solver statistics
  struct SolverStats {

    SolverStats() :
      num_nonlinear_solves(0),
      num_nonlinear_iterations(0),
      num_total_nonlinear_iterations(0)
    {}
      
    // Resets counters. Does not reset counters that aggregate across multiple nonlinear solves.
    void resetForNewNonlinearSolve()
    {
      num_nonlinear_iterations = 0;
      line_search.reset();
      trust_region.reset();
    }

    //! Increases the number of nonlinear solves by one.
    void incrementNumNonlinearSolves()
    { ++num_nonlinear_solves; }

    void incrementNumNonlinearIterations()
    {
      ++num_nonlinear_iterations;
      ++num_total_nonlinear_iterations;
    }

    //! Total number of nonlinear solves attempted.
    int num_nonlinear_solves;
    
    //! Number of nonlinear iterations in the last nonlinear solve.
    int num_nonlinear_iterations;

    //! Total number of nonlinear iterations for all nonlinear solves.
    int num_total_nonlinear_iterations;

    //! Statistics for the linear solve.
    struct LinearSolveStats {
      bool last_linear_solve_converged;
      int num_iterations_last_linear_solve;
      double linear_residual_norm_last_line_solve;
      int num_iterations_all_linear_solves;      
    };

    //! Line search stats for the last nonlinear solve
    NOX::LineSearchCounters line_search;

    //! Container for trust region statistics
    struct TrustRegionStats {
      TrustRegionStats() :
        num_cauchy_steps(0),
        num_newton_steps(0),
        num_dogleg_steps(0),
        dogleg_fraction_cauchy_to_newton(0.0),
        dogleg_fraction_newton_length(0.0)
      {}

      void reset()
      {
        num_cauchy_steps = 0;
        num_newton_steps = 0;
        num_dogleg_steps = 0;
        dogleg_fraction_cauchy_to_newton = 0.0;
        dogleg_fraction_newton_length = 0.0;
      }
      
      int num_cauchy_steps; //! Number of pure Cauchy steps taken
      int num_newton_steps; //! Number of pure Newton steps taken
      int num_dogleg_steps; //! Number of dogleg steps taken
      int num_trust_region_inner_iterations;
      double dogleg_fraction_cauchy_to_newton;
      double dogleg_fraction_newton_length;
    };
    
    //! Trust Regions stats for the last nonlinear solve
    TrustRegionStats trust_region;

  };

}

#endif
