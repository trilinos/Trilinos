// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_PROGRESSIVEHEDGING_H
#define ROL_PROGRESSIVEHEDGING_H

#include "ROL_OptimizationSolver.hpp"
#include "ROL_ProgHedgeObjective.hpp"

/** @ingroup algo_group
    \class ROL::ProgressiveHedging
    \brief Provides the interface to solve a stochastic program using progressive hedging.

    ---
*/


namespace ROL {

template <class Real>
class ProgressiveHedging {
private:
  const Ptr<OptimizationProblem<Real>> input_;
  const Ptr<SampleGenerator<Real>> sampler_;
  Teuchos::ParameterList parlist_;
  Real penaltyParam_;
  Real update_;
  int  freq_;
  Real ztol_;
  int maxit_;
  bool print_;

  Ptr<ProgHedgeObjective<Real>>  ph_objective_;
  Ptr<OptimizationProblem<Real>> ph_problem_;
  Ptr<OptimizationSolver<Real>>  ph_solver_;
  Ptr<Vector<Real>> z_psum_, z_gsum_;
  std::vector<Ptr<Vector<Real>>> zvec_, wvec_;

public:
  ProgressiveHedging(const Ptr<OptimizationProblem<Real>> &input,
                     const Ptr<SampleGenerator<Real>> &sampler,
                     Teuchos::ParameterList &parlist)
    : input_(input), sampler_(sampler), parlist_(parlist) {
    // Get algorithmic parameters
    penaltyParam_ = parlist.sublist("SOL").sublist("Progressive Hedging").get("Initial Penalty Parameter",10.0);
    update_       = parlist.sublist("SOL").sublist("Progressive Hedging").get("Penalty Update Scale",10.0);
    freq_         = parlist.sublist("SOL").sublist("Progressive Hedging").get("Penalty Update Frequency",0);
    ztol_         = parlist.sublist("SOL").sublist("Progressive Hedging").get("Nonanticipativity Constraint Tolerance",1e-4);
    maxit_        = parlist.sublist("SOL").sublist("Progressive Hedging").get("Iteration Limit",100);
    print_        = parlist.sublist("SOL").sublist("Progressive Hedging").get("Print Subproblem Solve History",false);
    // Create progressive hedging objective function
    ph_objective_ = makePtr<ProgHedgeObjective<Real>>(input_->getObjective(),
                                                      *input_->getSolutionVector(),
                                                      penaltyParam_);
    // Build progressive hedging subproblems
    ph_problem_   = makePtr<OptimizationProblem<Real>>(ph_objective_,
                                                       input_->getSolutionVector(),
                                                       input_->getBoundConstraint(),
                                                       input_->getConstraint(),
                                                       input_->getMultiplierVector());
    // Build progressive hedging subproblem solver
    ph_solver_    = makePtr<OptimizationSolver<Real>>(*ph_problem_, parlist);
    // Initialize vector storage
    z_psum_       = ph_problem_->getSolutionVector()->clone();
    z_gsum_       = ph_problem_->getSolutionVector()->clone();
    z_gsum_->set(*ph_problem_->getSolutionVector());
    zvec_.resize(sampler_->numMySamples());
    wvec_.resize(sampler_->numMySamples());
    for (int i = 0; i < sampler_->numMySamples(); ++i) {
      zvec_[i] = z_psum_->clone();        zvec_[i]->zero();
      wvec_[i] = z_psum_->dual().clone(); wvec_[i]->zero();
    }
  }

  void run(std::ostream &outStream = std::cout) {
    const Real zero(0), one(1);
    Real znorm_p(0), znorm_g(0), znorm(ROL_INF<Real>());
    int iter(0), spiter(0), tspiter(0);
    // Output
    outStream << std::scientific << std::setprecision(6);
    outStream << std::endl << "Progressive Hedging"
              << std::endl << "  "
              << std::setw(8)  << std::left << "iter"
              << std::setw(15) << std::left << "||z-Ez||"
              << std::setw(15) << std::left << "penalty"
              << std::setw(8)  << std::left << "subiter"
              << std::endl;
    while (iter < maxit_ && znorm > ztol_) {
      z_psum_->zero(); spiter = 0;
      ph_problem_->getSolutionVector()->set(*z_gsum_);
      // Solve concurrent optimization problems
      for (int j = 0; j < sampler_->numMySamples(); ++j) {
        ph_objective_->setData(*z_gsum_,*wvec_[j],penaltyParam_);
        ph_problem_->getObjective()->setParameter(sampler_->getMyPoint(j));
        if (ph_problem_->getConstraint() != nullPtr) {
          ph_problem_->getConstraint()->setParameter(sampler_->getMyPoint(j));
        }
        if (print_) {
          ph_solver_->solve(outStream);
        }
        else {
          ph_solver_->solve();
        }
        zvec_[j]->set(*ph_problem_->getSolutionVector());
        spiter += ph_solver_->getAlgorithmState()->iter;
        ph_solver_->reset();
        z_psum_->axpy(sampler_->getMyWeight(j),*zvec_[j]);
      }
      // Aggregation
      z_gsum_->zero();
      sampler_->sumAll(*z_psum_,*z_gsum_);
      // Multiplier Update
      znorm_p = zero; znorm_g = zero;
      for (int j = 0; j < sampler_->numMySamples(); ++j) {
        wvec_[j]->axpy(penaltyParam_,*zvec_[j]);
        wvec_[j]->axpy(-penaltyParam_,*z_gsum_);
        z_psum_->set(*zvec_[j]);
        z_psum_->axpy(-one,*z_gsum_);
        znorm_p += z_psum_->dot(*z_psum_);
      }
      sampler_->sumAll(&znorm_p,&znorm_g,1);
      znorm = std::sqrt(znorm_g);
      iter++;
      tspiter += spiter;
      // Output
      outStream << "  "
                << std::setw(8)  << std::left << iter
                << std::setw(15) << std::left << znorm
                << std::setw(15) << std::left << penaltyParam_
                << std::setw(8)  << std::left << spiter
                << std::endl;
      // Update penalty parameter
      if (freq_ > 0 && iter%freq_ == 0) {
        penaltyParam_ *= update_;
      }
    }
    input_->getSolutionVector()->set(*z_gsum_);
    // Output reason for termination
    if (znorm < ztol_) {
      outStream << "Converged: Nonanticipativity constraint tolerance satisfied!" << std::endl;
    }
    else {
      if (iter > maxit_) {
        outStream << "Maximum number of iterations exceeded" << std::endl;
      }
    }
    outStream << "Total number of subproblem iterations per sample: "
              << tspiter << " / " << sampler_->numGlobalSamples()
              << " ~ " << static_cast<int>(std::ceil(static_cast<Real>(tspiter)/static_cast<Real>(sampler_->numGlobalSamples())))
              << std::endl;
  }

}; // class ProgressiveHedging

} // namespace ROL

#endif
