
#ifndef ROL_PEBBL_INTERFACE_HPP
#define ROL_PEBBL_INTERFACE_HPP

//#include <acro_config.h>
#include <pebbl/bb/branching.h>
#include <pebbl/utilib/BitArray.h>
#include <pebbl/utilib/IntVector.h>
#include <pebbl/utilib/DoubleVector.h>
#include <pebbl/misc/chunkAlloc.h>

#include "ROL_OptimizationSolver.hpp"
#include "ROL_OptimizationProblemFactory.hpp"
#include "ROL_BranchHelper_PEBBL.hpp"
#include "ROL_Constraint_PEBBL.hpp"
#include "ROL_TransformedObjective_PEBBL.hpp"
#include "ROL_TransformedConstraint_PEBBL.hpp"


using namespace utilib;

namespace ROL {

template<class Real>
class ROL_PEBBL_BranchSub;

template<class Real>
class ROL_PEBBL_Solution : public pebbl::solution {
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
      if (std::abs(fx - x) < tol_) {
        return fx;
      }
      else if (std::abs(cx - x) < tol_) {
        return cx;
      }
      return x;
    }
  } snap;

public:
  ROL_PEBBL_Solution(const Vector<Real> &x, const Real newvalue) {
    solution_ = x.clone();
    solution_->set(x);
    solution_->applyUnary(snap);
    value = newvalue;
  }

  void copy(pebbl::solution* toCopy) {
    pebbl::solution::copy(toCopy);
    solution_.set(*dynamic_cast<ROL_PEBBL_Solution>(toCopy)->getVector());
  }

  const char* typeDescription() const { return "ROL solution"; }

  void print(std::ostream& s) {
    solution_->print(s);
  }

  //void printContents(std::ostream& /*s*/) { }

  const Ptr<Vector<Real>> getVector(void) const {
    return solution_;
  }
};

template<class Real>
class ROL_PEBBL_Branching : public pebbl::branching {
protected:
  // OptimizationProblem encapsulates the following problem
  // min        obj(x)
  // subject to xl <= x <= xu
  //            econ(x) = 0         (Lagrange Multiplier: emul)
  //            cl <= icon(x) <= cu (Lagrange Multiplier: imul)
  const Ptr<OptimizationProblemFactory<Real>> factory_;
  // Solver parameters
  const Ptr<ParameterList>                    parlist_;
  // Vector specific branching helper
  const Ptr<BranchHelper_PEBBL<Real>>         bHelper_;

  const int verbosity_;
  const Ptr<std::ostream> outStream_;

public:
  ROL_PEBBL_Branching(const Ptr<OptimizationProblemFactory<Real>> &factory,
                      const Ptr<ParameterList>                    &parlist,
                      const Ptr<BranchHelper_PEBBL<Real>>         &bHelper,
                      const int                                    verbosity = 0,
                      const Ptr<std::ostream>                     &outStream = nullPtr)
    : factory_(factory), parlist_(parlist), bHelper_(bHelper),
      verbosity_(verbosity), outStream_(outStream) {}

  virtual pebbl::branchSub* blankSub() {
    return new ROL_PEBBL_BranchSub<Real>(makePtrFromRef<ROL_PEBBL_Branching<Real>>(*this),verbosity_,outStream_);
  }

  void writeLoadLog(std::ostream &llFile, int proc = 0) {}

  const Ptr<OptimizationProblemFactory<Real>> getProblemFactory(void) const {
    return factory_;
  }

  bool haveIncumbentHeuristic() { return true; }

  const Ptr<ParameterList> getSolverParameters(void) const {
    return parlist_;
  }

//  const Ptr<Vector<Real>> getIncumbent(void) const {
//    return problem_->getSolutionVector();
//  }

  const Ptr<BranchHelper_PEBBL<Real>> getBranchHelper(void) const {
    return bHelper_;
  }
};

template<class Real>
class ROL_PEBBL_BranchSub : public pebbl::branchSub {
protected:
  const Ptr<ROL_PEBBL_Branching<Real>> branching_;
  const Ptr<BranchHelper_PEBBL<Real>> bHelper_;
  std::map<int,Real> fixed_;
  Ptr<Transform_PEBBL<Real>>  ptrans_;
  Ptr<Objective<Real>> tobj_;
  Ptr<Constraint<Real>> tcon_;
  Ptr<OptimizationProblem<Real>> problem_, problem0_;
  Ptr<OptimizationSolver<Real>> solver_;
  Ptr<Vector<Real>> solution_, rndSolution_;
  Ptr<Vector<Real>> multiplier_;
  Ptr<Vector<Real>> gradient_, dwa_;
  int nfrac_, index_;
  Real integralityMeasure_;
  bool hasConstraint_;

  const int verbosity_;
  const Ptr<std::ostream> outStream_;

  class round : public Elementwise::UnaryFunction<Real> {
  public:
    Real apply(const Real &x) const {
      return std::round(x);
    }
  } rnd;

public:
  ROL_PEBBL_BranchSub(const Ptr<ROL_PEBBL_Branching<Real>> &branching,
                      const int verbosity = 0,
                      const Ptr<std::ostream> &outStream = nullPtr)
    : branching_(branching),
      bHelper_(branching_->getBranchHelper()),
      nfrac_(-1), index_(-1), integralityMeasure_(-1),
      hasConstraint_(false),
      verbosity_(verbosity), outStream_(outStream) {
    problem0_   = branching_->getProblemFactory()->build();
    solution_   = problem0_->getSolutionVector()->clone();
    solution_->set(*problem0_->getSolutionVector());
    rndSolution_= problem0_->getSolutionVector()->clone();
    rndSolution_->set(*problem0_->getSolutionVector());
    rndSolution_->applyUnary(rnd);
    gradient_   = solution_->dual().clone();
    dwa_        = gradient_->clone();
    if ( problem0_->getMultiplierVector() != nullPtr ) {
      hasConstraint_ = true;
      multiplier_ = problem0_->getMultiplierVector()->clone();
    }
    ptrans_     = branching_->getBranchHelper()->createTransform();
  }

  ROL_PEBBL_BranchSub(const ROL_PEBBL_BranchSub &rpbs)
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

  pebbl::branching* bGlobal() const {
    return getRawPtr<ROL_PEBBL_Branching<Real>>(branching_);
  }

  void setRootComputation() {
    fixed_.clear();
  }

  void boundComputation(double* controlParam) {
    std::ios_base::fmtflags flags(outStream_->flags());
    if (verbosity_ > 0) {
      *outStream_ << std::scientific << std::setprecision(3);
      *outStream_ << std::endl << "In boundCompuation" << std::endl;
    }
    // Get base optimization solver parameters
    const Ptr<ParameterList> parlist
      = branching_->getSolverParameters();
    // Set fixed constraint
    ptrans_->add(fixed_);
    // Set solution and multiplier vectors
    Real tol = static_cast<Real>(1e-10);
    ptrans_->value(*solution_,*solution_,tol);
    if (verbosity_ > 2) {
      *outStream_ << std::endl << "After apply transformation" << std::endl;
      solution_->print(*outStream_);
      *outStream_ << std::endl;
    }
    // Construct new optimization problem from problem0_
    tobj_ = makePtr<TransformedObjective_PEBBL<Real>>(
              problem0_->getObjective(),ptrans_);
    if (hasConstraint_) {
      multiplier_->set(*problem0_->getMultiplierVector());
      tcon_ = makePtr<TransformedConstraint_PEBBL<Real>>(
                problem0_->getConstraint(),ptrans_);
      problem_ = makePtr<OptimizationProblem<Real>>(tobj_,
                                                    solution_,
                                                    problem0_->getBoundConstraint(),
                                                    tcon_,
                                                    multiplier_);
    }
    else {
      problem_ = makePtr<OptimizationProblem<Real>>(tobj_,
                                                    solution_,
                                                    problem0_->getBoundConstraint());
    }
    // Construct optimization solver
    solver_ = makePtr<OptimizationSolver<Real>>(*problem_,*parlist);
    // Solve optimization problem
    if (verbosity_ > 0) {
      if (verbosity_ > 3) {
        problem_->check(*outStream_);
      }
      solver_->solve(*outStream_);
      ptrans_->value(*solution_,*problem_->getSolutionVector(),tol);
      if (verbosity_ > 2) {
        *outStream_ << std::endl << "  ";
        solution_->print(*outStream_);
        *outStream_ << std::endl;
      }
      if (verbosity_ > 1) {
        *outStream_ << "  Problem ID:           " << id.serial << std::endl;
      }
    }
    else {
      solver_->solve();
    }
    // Get objective function value and status
    Ptr<const AlgorithmState<Real>> state = solver_->getAlgorithmState();
    Real value = state->value;
    EExitStatus statusFlag = state->statusFlag;
    if (statusFlag != EXITSTATUS_CONVERGED) {
      throw Exception::NotImplemented(">>> ROL_PEBBL_Interface::boundComputation : Bound problem did not converge!");
    }
    if (value > static_cast<Real>(bound)) {
      if (verbosity_ > 1) {
        *outStream_ << "  Previous Bound:       " << bound << std::endl;
      }
      bound = static_cast<double>(value);
      if (verbosity_ > 1) {
        *outStream_ << "  New Bound:            " << bound << std::endl;
      }
    }
    bHelper_->getNumFrac(nfrac_,integralityMeasure_,*solution_);
    if (verbosity_ > 1) {
      *outStream_ << "  Number of Fractional: " << nfrac_ << std::endl;
      *outStream_ << "  Integrality Measure:  " << integralityMeasure_ << std::endl;
      outStream_->flags(flags);
    }
    setState(bounded);
  }

  pebbl::branchSub* makeChild(int whichChild = anyChild) {
    if (whichChild == anyChild) {
      throw Exception::NotImplemented(">>> ROL_PEBBL_BranchSub::makeChild: whichChild is equal to anyChild!");
    }
    ROL_PEBBL_BranchSub<Real>* child
      = new ROL_PEBBL_BranchSub<Real>(*this);
    child->updateFixed(index_,
      (whichChild==0 ? static_cast<Real>(1) : static_cast<Real>(0)));
    return child;
  }

  bool candidateSolution() {
    return (nfrac_==0);
  }

  virtual void incumbentHeuristic() {
    Real tol(std::sqrt(ROL_EPSILON<Real>()));
    rndSolution_->set(*solution_);
    rndSolution_->applyUnary(rnd);
    problem0_->getObjective()->update(*rndSolution_);
    Real val = problem0_->getObjective()->value(*rndSolution_,tol);
    branching_->foundSolution(new ROL_PEBBL_Solution<Real>(*rndSolution_,val));
  }

  pebbl::solution* extractSolution() {
    Ptr<const AlgorithmState<Real>> state = solver_->getAlgorithmState();
    return new ROL_PEBBL_Solution<Real>(*solution_,state->value);
  }

  int splitComputation() {
    Real tol(std::sqrt(ROL_EPSILON<Real>()));
    tobj_->gradient(*gradient_,*solution_,tol);
    if (hasConstraint_) {
      tcon_->applyAdjointJacobian(*dwa_,*multiplier_,*solution_,tol);
      gradient_->plus(*dwa_);
    }
    index_ = bHelper_->getIndex(*solution_, *gradient_);
    //index_ = bHelper_->getIndex(*solution_);
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
void runSerial(const Ptr<ROL_PEBBL_Branching<Real>> &instance) {
  utilib::exception_mngr::set_stack_trace(false);
  utilib::exception_mngr::set_stack_trace(true);
  instance->reset();
  instance->solve();
}


} // ROL namespace
#endif
