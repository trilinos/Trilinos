// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Belos_InnerSolverResult_hpp
#define __Belos_InnerSolverResult_hpp

#include <BelosTypes.hpp> // ReturnType
#include <map>
#include <string>

namespace Belos {

  /// \class InnerSolveResult
  /// \brief Represents the result of an inner solve.
  /// \author Mark Hoemmen
  ///
  /// \warning This is EXPERIMENTAL CODE.  DO NOT RELY ON THIS CODE.
  ///   The interface or implementation may change at any time.
  ///
  /// An "inner solve" is an invocation of an iterative method, itself
  /// used as the operator or preconditioner for another iterative
  /// method.  Inner solves may converge or not converge; they take
  /// some number of iterations to do either, and that number is less
  /// than or equal to some maximum number of iterations.  The inner
  /// solve result may represent the result of solving for a single
  /// right-hand side, or it may aggregate the results of solving for
  /// multiple right-hand side(s).
  class InnerSolveResult {
  public:
    /// \brief Constructor
    ///
    /// \param theResult [in] Result of the solve, in numerical terms
    ///   (Converged or Unconverged)
    ///
    /// \param theNumRestartCycles [in] Number of restart cycles
    ///   started (not necessarily completed, if convergence happened
    ///   before the end of the restart cycle).  Zero means the
    ///   initial guess was accurate enough to be judged Converged.
    ///
    /// \param theTotalNumIters [in] Total number of iterations,
    ///   summed over all restart cycles.
    ///
    /// \param theExtraData [in] Any more data that the inner solve
    ///   implementation would like to expose to the outside world.
    ///   The data are stored as a map from the string label, to a
    ///   double-precision floating-point value.
    InnerSolveResult (const ReturnType theResult,
		      const int theNumRestartCycles,
		      const int theTotalNumIters,
		      const std::map<std::string, double>& theExtraData);

    //! Constructor, with no "extra data" input argument
    InnerSolveResult (const ReturnType theResult,
		      const int theNumRestartCycles,
		      const int theTotalNumIters);

    /// \brief Did the inner solve converge?    
    ///
    /// The ReturnType enum currently only has two values: Converged
    /// and Unconverged.  This may or may not change in the future,
    /// which is why we express convergence using ReturnType rather
    /// than a Boolean value.
    ReturnType result() {
      return result_;
    }

    /// \brief Total number of iterations completed over all restarts.
    /// 
    /// This is the sum of all iterations over all restarts.
    int totalNumIters() const { 
      return totalNumIters_;
    }

    //! Total number of restart cycles.
    int numRestartCycles() const { 
      return numRestartCycles_;
    }

    /// \brief "Extra" data from the inner solve.
    ///
    /// The inner solve may choose to expose more data to the outside
    /// world than just the data above.  The data are stored as a map
    /// from the string label, to a double-precision floating-point
    /// value.
    const std::map<std::string, double>& extraData () const {
      return extraData_;
    }

  private:
    ReturnType result_;
    int numRestartCycles_, totalNumIters_;
    std::map<std::string, double> extraData_;

    static ReturnType 
    validatedReturnType (const ReturnType ret);

    static int 
    requireNonNegInt (const int k);
  };

} // namespace Belos


#endif // __Belos_InnerSolverResult_hpp
