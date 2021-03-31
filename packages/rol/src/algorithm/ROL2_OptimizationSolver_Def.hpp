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

#pragma once
#ifndef ROL2_OPTIMIZATIONSOLVER_DEF_HPP
#define ROL2_OPTIMIZATIONSOLVER_DEF_HPP

namespace ROL {
namespace ROL2 {

template<typename Real>
OptimizationSolver<Real>::OptimizationSolver( const Ptr<OptimizationProblem<Real>> &opt,
                                                    ParameterList                           &parlist )
  : opt_(opt), problemType_(opt_->getProblemType()) {
  switch (problemType_) {
    case ProblemType::U:  algoU_ = TypeU::make_Algorithm<Real>(parlist); break;
    case ProblemType::B:  algoB_ = TypeB::make_Algorithm<Real>(parlist); break;
    case ProblemType::E:  algoE_ = TypeE::make_Algorithm<Real>(parlist); break;
    case ProblemType::G:  algoG_ = TypeG::make_Algorithm<Real>(parlist); break;
    case ProblemType::Last:
      ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
        "Error in OptimizationSolver::solve() : Unsupported problem type");
  }
}

template<typename Real>
std::vector<std::string> OptimizationSolver<Real>::getOutput() const {
  return output_;
}

template<typename Real>
int OptimizationSolver<Real>::solve(const Ptr<StatusTest<Real>> &status,
                                    bool combineStatus) {
  nullstream bhs;
  return solve(bhs,status,combineStatus);
}

template<typename Real>
int OptimizationSolver<Real>::solve( std::ostream &outStream,
                                     const Ptr<StatusTest<Real>> &status,
                                     bool combineStatus ) {
  switch (problemType_) {
    case ProblemType::U:
      if (status != nullPtr) algoU_->setStatusTest(status,combineStatus);
      output_ = algoU_->run(*opt_,outStream);
      break;
    case ProblemType::B:
      if (status != nullPtr) algoB_->setStatusTest(status,combineStatus);
      output_ = algoB_->run(*opt_,outStream);
      break;
    case ProblemType::E:
      if (status != nullPtr) algoE_->setStatusTest(status,combineStatus);
      output_ = algoE_->run(*opt_,outStream);
      break;
    case ProblemType::G:
      if (status != nullPtr) algoG_->setStatusTest(status,combineStatus);
      output_ = algoG_->run(*opt_,outStream);
      break;
    case ProblemType::Last:
      ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
        "Error in OptimizationSolver::solve() : Unsupported problem type");
  }
  // TODO: Interrogate AlgorithmState and StatusTest to generate a return code
  //       that indicates why the solver has stopped

  // Return an integer code
  return 0;
}

template<typename Real>
Ptr<const AlgorithmState<Real>>& OptimizationSolver<Real>::getAlgorithmState() const {
  switch (problemType_) {
    case ProblemType::U:  return algoU_->getState();
    case ProblemType::B:  return algoB_->getState();
    case ProblemType::E:  return algoE_->getState();
    case ProblemType::G:  return algoG_->getState();
    case ProblemType::Last:
    default:
      ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
        "Error in OptimizationSolver::getAlgorithmState() : Unsupported problem type");
  }
}

template<typename Real>
void OptimizationSolver<Real>::reset() {
  switch (problemType_) {
    case ProblemType::U:  algoU_->reset(); break;
    case ProblemType::B:  algoB_->reset(); break;
    case ProblemType::E:  algoE_->reset(); break;
    case ProblemType::G:  algoG_->reset(); break;
    case ProblemType::Last:
    default:
      ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
        "Error in OptimizationSolver::reset() : Unsupported problem type");
  }
}

} // namespace ROL2
} // namespace ROL

#endif // ROL2_OPTIMIZATIONSOLVER_DEF_HPP


