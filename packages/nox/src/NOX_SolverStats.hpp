// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_SOLVER_STATS_HPP
#define NOX_SOLVER_STATS_HPP

#include "NOX_Common.H"
#include "NOX_LineSearch_Utils_Counters.H"
#include <iosfwd>

namespace NOX {

  //! Container for solver statistics
  struct SolverStats {

    SolverStats() :
      numNonlinearSolves(0),
      numNonlinearIterations(0),
      numTotalNonlinearIterations(0),
      residual2Norm(0.0)
    {}
      
    // Resets counters for a new nonlinear solve. Does not reset counters that aggregate across multiple nonlinear solves.
    void reset()
    {
      numNonlinearIterations = 0;
      residual2Norm = 0.0;
      linearSolve.reset();
      lineSearch.reset();
      trustRegion.reset();
    }

    //! Increases the number of nonlinear solves by one.
    void incrementNumNonlinearSolves()
    { ++numNonlinearSolves; }

    void incrementNumNonlinearIterations()
    {
      ++numNonlinearIterations;
      ++numTotalNonlinearIterations;
    }

    //! Total number of nonlinear solves attempted.
    int numNonlinearSolves;
    
    //! Number of nonlinear iterations in the last nonlinear solve.
    int numNonlinearIterations;

    //! Total number of nonlinear iterations for all nonlinear solves.
    int numTotalNonlinearIterations;

    //! The 2-Norm of F at the end of the last nonlinear solve.
    double residual2Norm;

    //! Statistics for the linear solve.
    struct LinearSolveStats {
      LinearSolveStats() :
        lastLinearSolve_Converged(false),
        lastLinearSolve_NumIterations(0),
        lastLinearSolve_AchievedTolerance(0.0),
        lastLinearSolve_InitialResidualNorm(0.0),
        lastLinearSolve_FinalResidualNorm(0.0),
        lastNonlinearSolve_NumLinearSolves(0),
        lastNonlinearSolve_NumLinearIterations(0),
        allNonlinearSolves_NumLinearSolves(0),
        allNonlinearSolves_NumLinearIterations(0)
      {}

      void logLinearSolve(const bool converged,
                          const int numIterations,
                          const double achievedTolerance,
                          const double initialResidualNorm,
                          const double finalResidualNorm)
      {
        lastLinearSolve_Converged = converged;
        lastLinearSolve_NumIterations = numIterations;
        lastLinearSolve_AchievedTolerance = achievedTolerance;
        lastLinearSolve_InitialResidualNorm = initialResidualNorm;
        lastLinearSolve_FinalResidualNorm = finalResidualNorm;
        ++lastNonlinearSolve_NumLinearSolves;
        lastNonlinearSolve_NumLinearIterations += numIterations;
        ++allNonlinearSolves_NumLinearSolves;
        allNonlinearSolves_NumLinearIterations += numIterations;
      }

      // Don't reset the aggregate counters
      void reset() {
        lastLinearSolve_Converged = false;
        lastLinearSolve_NumIterations = 0;
        lastLinearSolve_AchievedTolerance = 0.0;
        lastLinearSolve_InitialResidualNorm = 0.0;
        lastLinearSolve_FinalResidualNorm = 0.0;
        lastNonlinearSolve_NumLinearSolves = 0;
        lastNonlinearSolve_NumLinearIterations = 0;
      }

      bool   lastLinearSolve_Converged;
      int    lastLinearSolve_NumIterations;
      double lastLinearSolve_AchievedTolerance;
      double lastLinearSolve_InitialResidualNorm;
      double lastLinearSolve_FinalResidualNorm;
      int    lastNonlinearSolve_NumLinearSolves;
      int    lastNonlinearSolve_NumLinearIterations;
      int    allNonlinearSolves_NumLinearSolves;
      int    allNonlinearSolves_NumLinearIterations;
    };

    LinearSolveStats linearSolve;

    //! Line search stats for the last nonlinear solve
    NOX::LineSearchCounters lineSearch;

    //! Container for trust region statistics
    struct TrustRegionStats {
      TrustRegionStats() :
        numCauchySteps(0),
        numNewtonSteps(0),
        numDoglegSteps(0),
        numTrustRegionInnerIterations(0),
        sumDoglegFractionCauchyToNewton(0.0),
        sumDoglegFractionNewtonLength(0.0)
      {}

      void reset()
      {
        numCauchySteps = 0;
        numNewtonSteps = 0;
        numDoglegSteps = 0;
        numTrustRegionInnerIterations = 0;
        sumDoglegFractionCauchyToNewton = 0.0;
        sumDoglegFractionNewtonLength = 0.0;
      }
      
      int numCauchySteps; //! Number of pure Cauchy steps taken
      int numNewtonSteps; //! Number of pure Newton steps taken
      int numDoglegSteps; //! Number of dogleg steps taken

      //! Number of inner iterations required to adjust the trust region radius.
      int numTrustRegionInnerIterations;

      /** \brief Holds the sum of the value of the fraction a dogleg
          step took between the Cauchy and Newton directions.  This is
          the \f$ \gamma \f$ variable in the standard dogleg algorithm
          and the \f$ \tau \f$ parameter in the inexact dogleg
          algorithm.  A value of 0.0 is a full step in the Cauchy
          direction and a value of 1.0 is a full step in the Newton
          direction. To get the average value for a nonlinear solve,
          divide this value by the number of dogleg steps.
       */
      double sumDoglegFractionCauchyToNewton;

      /** \brief Holds the sum of the values of the fraction a dogleg
          step took compared to the full Newton step.  The fractional
          value is computed as \f$ \mbox{frac} = \frac{\| d \|}{\|
          n\|} \f$. To get the average value for a nonlinear solve,
          divide this value by the number of dogleg steps.
      */
      double sumDoglegFractionNewtonLength;
    };
    
    //! Trust Regions stats for the last nonlinear solve
    TrustRegionStats trustRegion;

  };

}

#endif
