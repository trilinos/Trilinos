
#ifndef ROL_PEBBL_INTERFACE_HPP
#define ROL_PEBBL_INTERFACE_HPP

//#include <acro_config.h>
#include <pebbl/bb/branching.h>
#include <pebbl/utilib/CharString.h>
#include <pebbl/utilib/BasicArray.h>
#include <pebbl/utilib/BitArray.h>
#include <pebbl/utilib/IntVector.h>
#include <pebbl/utilib/DoubleVector.h>
#include <pebbl/utilib/_math.h>
#include <pebbl/utilib/ParameterSet.h>
#ifdef ACRO_HAVE_MPI
#include <pebbl/utilib/PackBuf.h>
#endif
#include <pebbl/misc/chunkAlloc.h>

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
    Real tol_;
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

  const char* typeDescription() const { return "ROL solution"; };

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
private:
  // OptimizationProblem encapsulates the following problem
  // min        obj(x)
  // subject to xl <= x <= xu
  //            econ(x) = 0         (Lagrange Multiplier: emul)
  //            cl <= icon(x) <= cu (Lagrange Multiplier: imul)
  const Ptr<OptimizationProblem<Real>> problem_;
  // Solver parameters
  const Ptr<ParameterList>             parlist_;
  // Vector specific branching helper
  const Ptr<BranchHelper_PEBBL<Real>>  bHelper_;

  const int verbosity_;
  const Ptr<std::ostream> outStream_;

public:
  ROL_PEBBL_Branching(const Ptr<OptimizationProblem<Real>> &problem,
                      const Ptr<ParameterList>             &parlist,
                      const Ptr<BranchHelper_PEBBL<Real>>  &bHelper,
                      const int verbosity = 0,
                      const Ptr<std::ostream> &outStream = nullPtr)
    : problem_(problem), parlist_(parlist), bHelper_(bHelper),
      verbosity_(verbosity), outStream_(outStream) {}

  pebbl::branchSub* blankSub() {
    return new ROL_PEBBL_BranchSub<Real>(makePtrFromRef<ROL_PEBBL_Branching<Real>>(*this),verbosity_,outStream_);
  }

  void writeLoadLog(std::ostream &llFile, int proc = 0) {}

  const Ptr<OptimizationProblem<Real>> getOptimizationProblem(void) const {
    return problem_;
  }

  const Ptr<ParameterList> getSolverParameters(void) const {
    return parlist_;
  }

  const Ptr<Vector<Real>> getIncumbent(void) const {
    return problem_->getSolutionVector();
  }

  const Ptr<BranchHelper_PEBBL<Real>> getBranchHelper(void) const {
    return bHelper_;
  }
};

template<class Real>
class ROL_PEBBL_BranchSub : public pebbl::branchSub {
private:
  const Ptr<ROL_PEBBL_Branching<Real>> branching_;
  const Ptr<BranchHelper_PEBBL<Real>> bHelper_;
  std::map<int,Real> fixed_;
  //Ptr<Constraint_PEBBL<Real>> fixedCon_;
  Ptr<Transform_PEBBL<Real>>  ptrans_;
  Ptr<Objective<Real>> tobj_;
  Ptr<Constraint<Real>> tcon_;
  Ptr<OptimizationProblem<Real>> problem_;
  Ptr<OptimizationSolver<Real>> solver_;
  Ptr<Vector<Real>> solution_;
  Ptr<Vector<Real>> multiplier_;
  int nfrac_, index_;
  Real integralityMeasure_;

  const int verbosity_;
  const Ptr<std::ostream> outStream_;

public:
  ROL_PEBBL_BranchSub(const Ptr<ROL_PEBBL_Branching<Real>> &branching,
                      const int verbosity = 0,
                      const Ptr<std::ostream> &outStream = nullPtr)
    : branching_(branching),
      bHelper_(branching_->getBranchHelper()),
      //fixedCon_(makePtr<Constraint_PEBBL<Real>>()),
      nfrac_(-1), index_(-1), integralityMeasure_(-1),
      verbosity_(verbosity), outStream_(outStream) {
    solution_   = branching_->getOptimizationProblem()->getSolutionVector()->clone();
    multiplier_ = branching_->getOptimizationProblem()->getMultiplierVector()->clone();
    ptrans_     = branching_->getBranchHelper()->createTransform();
  }

  ROL_PEBBL_BranchSub(const ROL_PEBBL_BranchSub &rpbs)
    : branching_(rpbs.branching_),
      bHelper_(rpbs.bHelper_),
      fixed_(rpbs.fixed_),
      //fixedCon_(makePtr<Constraint_PEBBL<Real>>()),
      nfrac_(-1), index_(-1), integralityMeasure_(-1), 
      verbosity_(rpbs.verbosity_), outStream_(rpbs.outStream_) {
    branchSubAsChildOf(this);
    solution_   = branching_->getOptimizationProblem()->getSolutionVector()->clone();
    multiplier_ = branching_->getOptimizationProblem()->getMultiplierVector()->clone();
    ptrans_     = branching_->getBranchHelper()->createTransform();
    bound = rpbs.bound;
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
    // Get base optimization problem 
    const Ptr<OptimizationProblem<Real>> problem0
      = branching_->getOptimizationProblem();
    // Get base optimization solver parameters
    const Ptr<ParameterList> parlist
      = branching_->getSolverParameters();
    // Set fixed constraint
    //fixedCon_->add(fixed_);
    ptrans_->add(fixed_);
    // Set solution and multiplier vectors
    //solution_->set(*problem0->getSolutionVector());
    Real tol = static_cast<Real>(1e-10);
    ptrans_->value(*solution_,*problem0->getSolutionVector(),tol);
    multiplier_->set(*problem0->getMultiplierVector());
//    // Build vector of equality constraints/multipliers
//    std::vector<Ptr<Constraint<Real>>> econ_vec(2,nullPtr);
//    econ_vec[0] = problem0->getConstraint();
//    econ_vec[1] = fixedCon_;
//    std::vector<Ptr<Vector<Real>>> emul_vec(2,nullPtr);
//    emul_vec[0] = multiplier_;
//    emul_vec[1] = fixedCon_->makeConstraintVector();
    // Construct new optimization problem from base
//    problem_
//      = makePtr<OptimizationProblem<Real>>(problem0->getObjective(),
//                                           solution_,
//                                           problem0->getBoundConstraint(),
//                                           econ_vec,
//                                           emul_vec);
    tobj_ = makePtr<TransformedObjective_PEBBL<Real>>(
              problem0->getObjective(),ptrans_);
    tcon_ = makePtr<TransformedConstraint_PEBBL<Real>>(
              problem0->getConstraint(),ptrans_);
    problem_ = makePtr<OptimizationProblem<Real>>(tobj_,
                                                  solution_,
                                                  problem0->getBoundConstraint(),
                                                  tcon_,
                                                  multiplier_);
    // Construct optimization solver
    solver_ = makePtr<OptimizationSolver<Real>>(*problem_,*parlist);
    // Solve optimization problem
    if (verbosity_ > 0) {
      if (verbosity_ > 3) {
        problem_->check(*outStream_);
      }
      solver_->solve(*outStream_);
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

  pebbl::solution* extractSolution() {
    Ptr<const AlgorithmState<Real>> state = solver_->getAlgorithmState();
    return new ROL_PEBBL_Solution<Real>(*solution_,state->value);
  }

  int splitComputation() {
    index_ = bHelper_->getIndex(*solution_);
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
