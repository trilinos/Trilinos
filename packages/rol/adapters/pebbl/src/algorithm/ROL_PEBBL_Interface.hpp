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

#ifndef ROL_PEBBL_INTERFACE_HPP
#define ROL_PEBBL_INTERFACE_HPP

//#include <acro_config.h>
#include <pebbl/bb/branching.h>
#include <pebbl/utilib/BitArray.h>
#include <pebbl/utilib/IntVector.h>
#include <pebbl/utilib/DoubleVector.h>
#include <pebbl/misc/chunkAlloc.h>

#include "ROL_NewOptimizationSolver.hpp"
#include "ROL_PEBBL_IntegerProblemFactory.hpp"
#include "ROL_PEBBL_BranchHelper.hpp"
#include "ROL_PEBBL_IntegerConstraint.hpp"

using namespace utilib;

namespace ROL {
namespace PEBBL {

template<class Real>
class BranchSub;

template<class Real>
class IntegerSolution : public pebbl::solution {
private:
  Ptr<Vector<Real>> solution_;

  class SnapInt : public Elementwise::UnaryFunction<Real> {
  private:
    const Real tol_;
  public:
    SnapInt(void) : tol_(1e-6) {}
    Real apply(const Real &x) const {
      Real fx = std::floor(x);
      Real cx = std::ceil(x);
      if (std::abs(fx - x) < tol_)      return fx;
      else if (std::abs(cx - x) < tol_) return cx;
      return x;
    }
  } snap;

  Ptr<Vector<Real>> getOptVector(const Ptr<Vector<Real>> &xs ) const {
    try {
      return dynamicPtrCast<PartitionedVector<Real>>(xs)->get(0);
    }
    catch (std::exception &e) {
      return xs;
    }
  }

  Ptr<Vector<Real>> getIntegerVector(const Ptr<Vector<Real>> &x) const {
    try {
      return dynamicPtrCast<MixedVector<Real>>(getOptVector(x))->getIntegerVariables();
    }
    catch (std::exception &e) {
      return getOptVector(x);
    }
  }

public:
  IntegerSolution(const Vector<Real> &x, const Real newvalue) {
    solution_ = x.clone();
    solution_->set(x);
    getIntegerVector(solution_)->applyUnary(snap);
    value = newvalue;
  }

  void copy(pebbl::solution* toCopy) {
    pebbl::solution::copy(toCopy);
    solution_.set(*dynamic_cast<IntegerSolution>(toCopy)->getVector());
  }

  const char* typeDescription() const { return "ROL/PEBBL integer solution"; }

  void print(std::ostream& s) { solution_->print(s); }

  //void printContents(std::ostream& /*s*/) { }

  const Ptr<const Vector<Real>> getVector(void) const { return solution_; }
};

template<class Real>
class Branching : public pebbl::branching {
protected:
  // OptimizationProblem encapsulates the following problem
  // min        obj(x)
  // subject to xl <= x <= xu
  //            econ(x) = 0         (Lagrange Multiplier: emul)
  //            cl <= icon(x) <= cu (Lagrange Multiplier: imul)
  const Ptr<IntegerProblemFactory<Real>> factory_;
  // Solver parameters
  const Ptr<ParameterList>               parlist_;
  // Vector specific branching helper
  const Ptr<BranchHelper<Real>>          bHelper_;

  const int                              verbosity_;
  const Ptr<std::ostream>                outStream_;

public:
  Branching(const Ptr<IntegerProblemFactory<Real>> &factory,
            const Ptr<ParameterList>               &parlist,
            const Ptr<BranchHelper<Real>>          &bHelper,
            int                                     verbosity = 0,
            const Ptr<std::ostream>                &outStream = nullPtr)
    : factory_(factory), parlist_(parlist), bHelper_(bHelper),
      verbosity_(verbosity), outStream_(outStream) {}

  virtual pebbl::branchSub* blankSub() {
    return new BranchSub<Real>(makePtrFromRef<Branching<Real>>(*this),verbosity_,outStream_);
  }

  void writeLoadLog(std::ostream &llFile, int proc = 0) {}

  const Ptr<IntegerProblemFactory<Real>> getProblemFactory(void) const {
    return factory_;
  }

  bool haveIncumbentHeuristic() { return true; }

  const Ptr<ParameterList> getSolverParameters(void) const { return parlist_; }

//  const Ptr<Vector<Real>> getIncumbent(void) const {
//    return problem_->getSolutionVector();
//  }

  const Ptr<BranchHelper<Real>> getBranchHelper(void) const { return bHelper_; }
};

template<class Real>
class BranchSub : public pebbl::branchSub {
protected:
  const Ptr<Branching<Real>>       branching_;
  const Ptr<BranchHelper<Real>>    bHelper_;
  std::map<int,Real>               fixed_;
  Ptr<IntegerTransformation<Real>> ptrans_;
  Ptr<IntegerProblem<Real>>        problem0_;
  Ptr<NewOptimizationSolver<Real>> solver_;
  Ptr<Vector<Real>>                solution_, rndSolution_;
  Ptr<Vector<Real>>                multiplier_;
  Ptr<Vector<Real>>                gradient_, dwa_;
  int                              nfrac_, index_;
  Real                             integralityMeasure_;
  bool                             hasConstraint_;

  const int                        verbosity_;
  const Ptr<std::ostream>          outStream_;

  class round : public Elementwise::UnaryFunction<Real> {
  public:
    Real apply(const Real &x) const { return std::round(x); }
  } rnd;

  Ptr<Vector<Real>> getOptVector(const Ptr<Vector<Real>> &xs ) const {
    try {
      return dynamicPtrCast<PartitionedVector<Real>>(xs)->get(0);
    }
    catch (std::exception &e) {
      return xs;
    }
  }

  Ptr<Vector<Real>> getIntegerVector(const Ptr<Vector<Real>> &x) const {
    try {
      return dynamicPtrCast<MixedVector<Real>>(getOptVector(x))->getIntegerVariables();
    }
    catch (std::exception &e) {
      return getOptVector(x);
    }
  }

public:
  BranchSub(const Ptr<Branching<Real>> &branching,
            int                         verbosity = 0,
            const Ptr<std::ostream>    &outStream = nullPtr)
    : branching_(branching),
      bHelper_(branching_->getBranchHelper()),
      nfrac_(-1), index_(-1), integralityMeasure_(-1),
      hasConstraint_(false),
      verbosity_(verbosity), outStream_(outStream) {
    problem0_   = branching_->getProblemFactory()->build();
    solution_   = problem0_->getPrimalOptimizationVector()->clone();
    rndSolution_= problem0_->getPrimalOptimizationVector()->clone();
    solution_->set(*problem0_->getPrimalOptimizationVector());
    rndSolution_->set(*problem0_->getPrimalOptimizationVector());
    rndSolution_->applyUnary(rnd);
    gradient_   = solution_->dual().clone();
    dwa_        = gradient_->clone();
    if ( problem0_->getMultiplierVector() != nullPtr ) {
      hasConstraint_ = true;
      multiplier_ = problem0_->getMultiplierVector()->clone();
    }
    ptrans_     = branching_->getBranchHelper()->createTransform();
  }

  BranchSub(const BranchSub &rpbs)
    : branching_(rpbs.branching_),
      bHelper_(rpbs.bHelper_),
      fixed_(rpbs.fixed_),
      nfrac_(-1), index_(-1), integralityMeasure_(-1), 
      hasConstraint_(rpbs.hasConstraint_),
      verbosity_(rpbs.verbosity_), outStream_(rpbs.outStream_) {
    branchSubAsChildOf(this);
    problem0_   = rpbs.branching_->getProblemFactory()->build();
    solution_   = rpbs.solution_->clone();
    solution_->set(*rpbs.solution_);
    rndSolution_ = rpbs.rndSolution_->clone();
    rndSolution_->set(*rpbs.rndSolution_);
    gradient_   = rpbs.gradient_->clone();
    gradient_->set(*rpbs.gradient_);
    dwa_        = rpbs.gradient_->clone();
    if (hasConstraint_) {
      multiplier_ = rpbs.multiplier_->clone();
      multiplier_->set(*rpbs.multiplier_);
    }
    ptrans_     = rpbs.branching_->getBranchHelper()->createTransform();
    bound       = rpbs.bound;
  }

  pebbl::branching* bGlobal() const { return getRawPtr<Branching<Real>>(branching_); }

  void setRootComputation() { fixed_.clear(); }

  void boundComputation(double* controlParam) {
    std::ios_base::fmtflags flags(outStream_->flags());
    if (verbosity_ > 0) {
      *outStream_ << std::scientific << std::setprecision(3);
      *outStream_ << std::endl << "In boundCompuation" << std::endl;
    }
    if (verbosity_ > 1) {
      *outStream_ << "  Fixed Components" << std::endl;
      *outStream_ << fixed_ << std::endl << std::endl;
    }
    // Get base optimization solver parameters
    const Ptr<ParameterList>
      parlist = branching_->getSolverParameters();
    // Set fixed constraint
    problem0_->edit();
    problem0_->deleteTransformation();
    if (!fixed_.empty()) {
      ptrans_->add(fixed_);
      problem0_->setTransformation(ptrans_);
      //problem0_->finalize(false,false);
      //problem0_->edit();
      //problem0_->deleteTransformation();
      //Ptr<IntegerConstraint<Real>> con = makePtr<IntegerConstraint<Real>>();
      //con->add(fixed_);
      //Ptr<Vector<Real>> mul = con->makeConstraintVector();
      //problem0_->removeLinearConstraint("Integer");
      //problem0_->addLinearConstraint("Integer",con,mul);
    }
    problem0_->setProjectionAlgorithm(*parlist);
    if (verbosity_ > 0) problem0_->finalize(false,true,*outStream_);
    else                problem0_->finalize(false,false);

    // Construct optimization solver
    solver_ = makePtr<NewOptimizationSolver<Real>>(problem0_,*parlist);
    // Solve optimization problem
    if (verbosity_ > 0) {
      if (verbosity_ > 3) problem0_->check(true,*outStream_);
      solver_->solve(*outStream_);
      if (verbosity_ > 1) *outStream_ << "  Problem ID:           " << id.serial << std::endl;
    }
    else {
      solver_->solve();
    }
    // Reset problem0 for use in other computations
    if (!fixed_.empty()) {
      problem0_->edit();
      problem0_->deleteTransformation();
    }
    // Get objective function value and status
    Real tol = static_cast<Real>(1e-10);
    Ptr<const AlgorithmState<Real>> state = solver_->getAlgorithmState();
    Real value = state->value;
    if (!fixed_.empty()) ptrans_->value(*solution_,*state->iterateVec,tol);
    else                 solution_->set(*state->iterateVec);
    if (hasConstraint_) multiplier_->set(*state->lagmultVec);
    EExitStatus statusFlag = state->statusFlag;
    ROL_TEST_FOR_EXCEPTION(statusFlag!=EXITSTATUS_CONVERGED,std::logic_error,
      ">>> ROL_PEBBL_Interface::boundComputation : Bound problem did not converge!");
    if (verbosity_ > 2) {
      *outStream_ << std::endl << "  ";
      solution_->print(*outStream_);
      *outStream_ << std::endl;
    }
    if (value > static_cast<Real>(bound)) {
      if (verbosity_ > 1) *outStream_ << "  Previous Bound:       " << bound << std::endl;
      bound = static_cast<double>(value);
      if (verbosity_ > 1) *outStream_ << "  New Bound:            " << bound << std::endl;
    }
    bHelper_->getNumFrac(nfrac_,integralityMeasure_,*solution_);
    if (verbosity_ > 1) {
      *outStream_ << "  Number of Fractional: " << nfrac_ << std::endl;
      *outStream_ << "  Integrality Measure:  " << integralityMeasure_ << std::endl;
      outStream_->flags(flags);
    }
    setState(bounded);
  }

  virtual pebbl::branchSub* makeChild(int whichChild = anyChild) {
    ROL_TEST_FOR_EXCEPTION(whichChild==anyChild,std::logic_error,
      ">>> ROL_PEBBL_BranchSub::makeChild: whichChild is equal to anyChild!");
    BranchSub<Real>* child = new BranchSub<Real>(*this);
    child->updateFixed(index_, (whichChild==0 ? static_cast<Real>(1) : static_cast<Real>(0)));
    return child;
  }

  bool candidateSolution() { return (nfrac_==0); }

  virtual void incumbentHeuristic() {
    Real tol(std::sqrt(ROL_EPSILON<Real>()));
    rndSolution_->set(*solution_);
    getIntegerVector(rndSolution_)->applyUnary(rnd);
    problem0_->getObjective()->update(*rndSolution_,UPDATE_TEMP);
    Real val = problem0_->getObjective()->value(*rndSolution_,tol);
    branching_->foundSolution(new IntegerSolution<Real>(*rndSolution_,val));
  }

  pebbl::solution* extractSolution() {
    Ptr<const AlgorithmState<Real>> state = solver_->getAlgorithmState();
    return new IntegerSolution<Real>(*solution_,state->value);
  }

  int splitComputation() {
    Real tol(std::sqrt(ROL_EPSILON<Real>()));
    problem0_->getObjective()->gradient(*gradient_,*solution_,tol);
    if (hasConstraint_) {
      problem0_->getConstraint()->applyAdjointJacobian(*dwa_,*multiplier_,*solution_,tol);
      gradient_->plus(*dwa_);
    }
    index_ = bHelper_->getIndex(*solution_, *gradient_);
    if (verbosity_ > 1) {
      std::ios_base::fmtflags flags(outStream_->flags());
      *outStream_ << std::scientific << std::setprecision(3);
      *outStream_ << std::endl << "In splitComputation" << std::endl;
      *outStream_ << "  Index:                " << index_ << std::endl;
      outStream_->flags(flags);
    }
    setState(separated);
    return 2;
  }

  void updateFixed(const int ind, const Real val) { 
    if (verbosity_ > 1) {
      std::ios_base::fmtflags flags(outStream_->flags());
      *outStream_ << std::scientific << std::setprecision(3);
      *outStream_ << std::endl << "In updateFixed:" << std::endl;
      *outStream_ << "  Index:                " << ind << std::endl;
      *outStream_ << "  Value:                " << val << std::endl;
      outStream_->flags(flags);
    }
    fixed_.insert(std::pair<int,Real>(ind,val));
  }

};

template<class Real>
void runSerial(const Ptr<Branching<Real>> &instance) {
  utilib::exception_mngr::set_stack_trace(false);
  utilib::exception_mngr::set_stack_trace(true);
  instance->reset();
  instance->solve();
}

} // namespace PEBBL
} // namespace ROL
#endif
