//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef __Belos_GmresInnerSolver_hpp
#define __Belos_GmresInnerSolver_hpp

#include <BelosInnerSolver.hpp>
#include <BelosGmresSolMgr.hpp>

namespace Belos {

  template<class Scalar, class MV, class OP>
  class GmresInnerSolver : public InnerSolver<Scalar, MV, OP> {
  public:
    typedef InnerSolver<Scalar, MV, OP> base_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::magnitude_type magnitude_type;
    typedef typename base_type::multivector_type multivector_type;
    typedef typename base_type::operator_type operator_type;

    /// \brief Constructor.
    ///
    /// \param problem [in/out] The linear problem to solve.  Its
    ///   initial guess ("left-hand side") and its right-hand side
    ///   will be changed every time solve() is called.  The contents
    ///   of those multivectors won't be overwritten -- just the
    ///   pointers (RCPs) will be changed -- so if you want to keep
    ///   either of them, you need only save the RCP before calling
    ///   any of the solve() methods.
    ///
    /// \param params [in] Parameters for the solve.  This is 
    ///   a pointer to nonconst only because the SolutionManager
    ///   interface demands it (more or less...).
    ///
    /// \param debug [in] Whether or not to run the GMRES
    ///   implementation in debug mode (which will produce verbose
    ///   debugging output, and may also result in extra computation
    ///   in order to test certain invariants or display status).
    GmresInnerSolver (const Teuchos::RCP<LinearProblem<Scalar,MV,OP> >& problem,
		      const Teuchos::RCP<Teuchos::ParameterList>& params,
		      const bool debug = false)
      : solMgr_ (problem, params, debug)
    {}

    /// \brief Solve \f$AX=B\f$ for the given right-hand side(s) B.
    ///
    /// \param X [in/out] On input: The initial guess for the inner
    ///   solver, if the inner solver accepts an initial guess (it is
    ///   not required to do so).  On output: the approximate solution
    ///   to Ax=B as computed by the inner solver.  Whether or not the
    ///   solver accepts an initial guess, X must be allocated to hold
    ///   the output, and it must be in the correct vector space for
    ///   the solution vectors.
    ///
    /// \param B [in] Right-hand side(s) for which to solve
    ///
    /// \param convTol [in] "Convergence tolerance," the meaning of
    ///   which depends on the subclass
    /// \param maxItersPerRestart [in] Maximum number of iterations
    ///   per restart cycle in the inner solve.
    /// \param maxNumRestarts [in] Maximum number of restart cycle(s) 
    ///   in the inner solve.
    ///
    /// \return The result of the inner solve.  It is a single result,
    ///   aggregated over all right-hand side(s).
    InnerSolveResult
    solve (const Teuchos::RCP<MV>& X,
	   const Teuchos::RCP<const MV>& B,
	   const magnitude_type convTol,
	   const int maxItersPerRestart,
	   const int maxNumRestarts)
    {
      solMgr_.setRHS (B);
      solMgr_.setLHS (X);
      // FIXME (mfh 07 Mar 2011) 
      //
      // I have a vague feeling that the stopping criteria should be
      // changed _after_ setting the RHS, but I'm not sure if this is
      // true.
      solMgr_.changeStoppingCriteria (convTol, maxItersPerRestart, maxNumRestarts);
      return invokeSolver ();
    }

    /// \brief Solve \f$AX=B\f$ for the given right-hand side(s) B.
    ///
    /// This should do the same thing as the five-argument version of
    /// solve(), except it should pick reasonable defaults for the
    /// convergence tolerance, maximum number of iterations, and
    /// maximum number of restart cycles.
    /// 
    /// \param X [in/out] On input: The initial guess for the inner
    ///   solver, if the inner solver accepts an initial guess (it is
    ///   not required to do so).  On output: the approximate solution
    ///   to Ax=B as computed by the inner solver.  
    ///
    /// \param B [in] Right-hand side(s) for which to solve
    ///
    /// \return The result of the inner solve.  It is a single result,
    ///   aggregated over all right-hand side(s).
    ///
    InnerSolveResult
    solve (const Teuchos::RCP<MV>& X,
	   const Teuchos::RCP<const MV>& B)
    {
      // The solver manager retains the last configuration of
      // convergence tolerance, max number of iterations per restart
      // cycle, and max number of restart cycles. Just set the
      // left-hand side X and right-hand side B, and call solve().
      solMgr_.setLHS (X);
      solMgr_.setRHS (B);
      return invokeSolver ();
    }

  private:
    /// \brief GMRES solver manager (implementation of GMRES).
    /// 
    /// The solver manager configures GMRES on construction, and can
    /// be used to solve different linear systems with the same matrix
    /// and different right-hand sides.
    GmresSolMgr<Scalar, MV, OP> solMgr_;

    /// \brief Invoke the solver for the current linear system.
    ///
    /// This method assumes that the linear system (more precisely,
    /// the left-hand side X and the right-hand side B) have already
    /// been set.
    InnerSolveResult
    invokeSolver ()
    {
      using std::pair;
      using std::vector;

      ReturnType result = solMgr_.solve ();
      TEST_FOR_EXCEPTION(result != Converged && result != Unconverged,
			 std::logic_error,
			 "The solver manager returned an invalid ResultType "
			 "enum value " << result << ".  Valid values are "
			 "Converged=" << Converged << " and Unconverged=" 
			 << Unconverged << ".");
      //
      // Compute max (total number of iterations, number of restart cycles)
      // for all right-hand side(s).
      //
      const vector<pair<int, int> > iterInfo = solMgr_.totalNumIters();
      // If we could use lambdas, or if C++ had better iteration
      // primitives, we wouldn't need a for loop here.
      typedef vector<pair<int, int>> >::const_iterator iter_type;
      int totalNumIters = 0, numRestartCycles = 0;
      for (iter_type it = iterInfo.begin(); it != iterInfo.end(); ++it)
	{
	  numRestartCycles = std::max (it->first, numRestartCycles);
	  totalNumIters = std::max (it->second, totalNumIters);
	}
      return InnerSolveResult (result, numRestartCycles, totalNumIters);
    }

  };

} // namespace Belos

#endif // __Belos_GmresInnerSolver_hpp
