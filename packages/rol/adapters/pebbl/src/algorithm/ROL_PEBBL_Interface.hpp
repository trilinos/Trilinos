// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PEBBL_INTERFACE_HPP
#define ROL_PEBBL_INTERFACE_HPP

//#include <acro_config.h>
#include <pebbl/bb/branching.h>
#ifdef HAVE_MPI
#include <pebbl/pbb/parTeamBranching.h>
#endif
#include <pebbl/utilib/BitArray.h>
#include <pebbl/utilib/IntVector.h>
#include <pebbl/utilib/DoubleVector.h>
#include <pebbl/misc/chunkAlloc.h>

#include "ROL_Solver.hpp"
#include "ROL_PEBBL_IntegerProblemFactory.hpp"
#include "ROL_PEBBL_BranchHelper.hpp"
#include "ROL_PEBBL_IntegerConstraint.hpp"
#include "ROL_PEBBL_MixedVector.hpp"

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
      return dynamic_cast<PartitionedVector<Real>&>(*xs).get(0);
    }
    catch (std::exception &e) {
      return xs;
    }
  }

  Ptr<Vector<Real>> getIntegerVector(const Ptr<Vector<Real>> &x) const {
    try {
      return dynamic_cast<MixedVector<Real>&>(*getOptVector(x)).getIntegerVariables();
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

#ifdef HAVE_MPI
template<class Real>
class Branching : virtual public pebbl::parallelTeamBranching {
#else
template<class Real>
class Branching : virtual public pebbl::branching {
#endif
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

  Ptr<Vector<Real>>                      initGuess_;
#ifdef HAVE_MPI
  pebbl::parallelBranchSub* mySP_;
  //pebbl::mpiComm teamComm;
#endif

  Ptr<Vector<Real>> getOptVector( Vector<Real> &xs ) const {
    try {
      return dynamic_cast<PartitionedVector<Real>&>(xs).get(0);
    }
    catch (std::exception &e) {
      return makePtrFromRef(xs);
    }
  }

  Ptr<Vector<Real>> getIntegerVector(Vector<Real> &xs) const {
    try {
      return dynamic_cast<MixedVector<Real>&>(*getOptVector(xs)).getIntegerVariables();
    }
    catch (std::exception &e) {
      return getOptVector(xs);
    }
  }

public:
  Branching(const Ptr<IntegerProblemFactory<Real>> &factory,
            const Ptr<ParameterList>               &parlist,
            const Ptr<BranchHelper<Real>>          &bHelper,
            int                                     verbosity = 0,
            const Ptr<std::ostream>                &outStream = nullPtr)
    : factory_(factory), parlist_(parlist), bHelper_(bHelper),
      verbosity_(verbosity), outStream_(outStream) {}

  /****************************************************************************/
  /* BEGIN VIRTUAL FUNCTION DECLARATION                                       */
  /****************************************************************************/
  // Initialize application branching method
  virtual void initialize(void) {};
  // Create new serial subproblem
  virtual pebbl::branchSub* blankSub() {
    return new BranchSub<Real>(makePtrFromRef<Branching<Real>>(*this),verbosity_,outStream_);
  }
#ifdef HAVE_MPI
  // Create new parallel subproblem
  virtual pebbl::parallelBranchSub* blankParallelSub() {
    return new BranchSub<Real>(makePtrFromRef<Branching<Real>>(*this),verbosity_,outStream_);
  }
#endif
  /****************************************************************************/
  /* END VIRTUAL FUNCTION DECLARATIONS                                        */
  /****************************************************************************/


  /****************************************************************************/
  /* BEGIN pebbl::branching DECLARATIONS                                      */
  /****************************************************************************/
  void writeLoadLog(std::ostream &llFile, int proc = 0) {}
  const Ptr<IntegerProblemFactory<Real>> getProblemFactory(void) const {
    return factory_;
  }
  bool haveIncumbentHeuristic() {
#ifdef HAVE_MPI
    return false;
#else
    return true;
#endif
}
  /****************************************************************************/
  /* END pebbl::branching DECLARATIONS                                        */
  /****************************************************************************/


  /****************************************************************************/
  /* BEGIN pebbl::parallelBranching DECLARATIONS                              */
  /****************************************************************************/
#ifdef HAVE_MPI
  void setupSearchComm() {
    searchComm.free();
    teamComm.free();
    splitCommunicator();
    //std::cout << "END setupSearchComm " << teamComm.myRank() << "  " << searchComm.myRank() << "\n";
    factory_->setCommunicator(teamComm.myComm());
    initialize();
  }
  void pack(pebbl::PackBuffer & outBuffer) {}
  void unpack(pebbl::UnPackBuffer & outBuffer) {}
  int spPackSize() {
    return getIntegerVector(*factory_->build()->getPrimalOptimizationVector())->dimension() * (sizeof(Real)+sizeof(int)) + sizeof(int);
  }
  //// True if this processor is the head of a team
  //bool iAmHead() {
  //  return teamComm.myRank() == 0;
  //}
  //void reset(bool VBFlag=true) {
  //  if (iAmHead()) {
  //    parallelBranching::reset(VBFlag);
  //  }
  //  else {
  //    branching::reset(VBFlag);
  //  }
  //}
  //double search(){
  //  if(iAmHead()){
  //    double objVal = parallelSearchFramework(NULL);
  //    return objVal;
  //  }
  //  else {
  //    return nan("");
  //  }
  //}
  void teamOrganize() {
    ucout << "ROL::PEBBL::Branching::teamOrganize called" << std::endl;
    int one = 1;
    int countedSize = -1;
    teamComm.reduceCast(&one,&countedSize,1,MPI_INT,MPI_SUM);
    mySP_ = blankParallelSub();
    ucout << "teamComm current appears to have " << countedSize << " MPI ranks" << std::endl;
  }
  void printSerial(const string& descrip) {
    int serialNumber = -1;
    teamComm.broadcast(&serialNumber,1,MPI_INT,0);
    ucout << "Minion " << teamComm.myRank() << ' ' << descrip << " subproblem " 
          << serialNumber << std::endl;
   
  }
  void minionBound() {
    //std::cout << "In MinionBound: "
    //          << teamComm.myRank()   << "  " << teamComm.mySize()   << "  "
    //          << passedComm.myRank() << "  " << passedComm.mySize() << std::endl;
    double controlParam(0);
    mySP_->boundComputation(&controlParam);
  }
  void minionSplit() {}
  void minionMakeChild() {}
#endif
  /****************************************************************************/
  /* END pebbl::parallelBranching DECLARATIONS                                */
  /****************************************************************************/

  void setInitialGuess(const Vector<Real> &initGuess) {
    if (initGuess_ == nullPtr) initGuess_ = initGuess.clone();
    initGuess_->set(initGuess);
  }
  void getInitialGuess(Vector<Real> &vec) const {
    if (initGuess_ != nullPtr) vec.set(*initGuess_);
  }

  // Setup function from pebbl::branching
  bool setup(int& argc, char**& argv) {
    bool flag(false);
#ifdef HAVE_MPI
    flag = pebbl::parallelBranching::setup(argc,argv);
#else
    initialize();
    flag = pebbl::branching::setup(argc,argv);
#endif
    return flag;
  }


  /****************************************************************************/
  /* BEGIN ROL/PEBBL INTERFACE DECLARATIONS                                   */
  /****************************************************************************/
  const Ptr<ParameterList> getSolverParameters(void) const { return parlist_; }
//  const Ptr<Vector<Real>> getIncumbent(void) const {
//    return problem_->getSolutionVector();
//  }
  const Ptr<BranchHelper<Real>> getBranchHelper(void) const { return bHelper_; }
  /****************************************************************************/
  /* END ROL/PEBBL INTERFACE DECLARATIONS                                     */
  /****************************************************************************/
};

#ifdef HAVE_MPI
template<class Real>
class BranchSub : public pebbl::parallelBranchSub {
#else
template<class Real>
class BranchSub : public pebbl::branchSub {
#endif
protected:
  const Ptr<Branching<Real>>       branching_;
  const Ptr<BranchHelper<Real>>    bHelper_;
  std::map<int,Real>               fixed_;
  Ptr<IntegerTransformation<Real>> ptrans_;
  Ptr<IntegerProblem<Real>>        problem0_;
  Ptr<Solver<Real>>                solver_;
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
      return dynamic_cast<PartitionedVector<Real>&>(*xs).get(0);
    }
    catch (std::exception &e) {
      return xs;
    }
  }

  Ptr<Vector<Real>> getIntegerVector(const Ptr<Vector<Real>> &x) const {
    try {
      return dynamic_cast<MixedVector<Real>&>(*getOptVector(x)).getIntegerVariables();
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
#ifndef HAVE_MPI
    initialize();
#endif
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
    ptrans_     = rpbs.branching_->getBranchHelper()->createTransform();
    bound       = rpbs.bound;
#ifndef HAVE_MPI
    initialize();
    solution_->set(*rpbs.solution_);
    rndSolution_->set(*rpbs.rndSolution_);
    gradient_->set(*rpbs.gradient_);
    if (hasConstraint_) multiplier_->set(*rpbs.multiplier_);
#endif
  }

  void initialize(void) {
    problem0_    = branching_->getProblemFactory()->build();
    solution_    = problem0_->getPrimalOptimizationVector()->clone();
    rndSolution_ = solution_->clone();
    gradient_    = solution_->dual().clone();
    dwa_         = gradient_->clone();
    if ( problem0_->getMultiplierVector() != nullPtr ) {
      hasConstraint_ = true;
      multiplier_ = problem0_->getMultiplierVector()->clone();
    }
    solution_->set(*problem0_->getPrimalOptimizationVector());
    rndSolution_->set(*problem0_->getPrimalOptimizationVector());
    rndSolution_->applyUnary(rnd);
  }

  pebbl::branching* bGlobal() const { return getRawPtr<Branching<Real>>(branching_); }
#ifdef HAVE_MPI
  pebbl::parallelBranching* pGlobal() const { return getRawPtr<Branching<Real>>(branching_); }

  void pack(pebbl::PackBuffer & outBuffer) {
    outBuffer << (int)fixed_.size();
    for (auto it = fixed_.begin(); it != fixed_.end(); ++it) {
      outBuffer << *it;
    }
    //outBuffer << fixed_;
    appPack(outBuffer);
  }

  void unpack(pebbl::UnPackBuffer & inBuffer) {
    int size(0);
    inBuffer >> size;
    fixed_.clear();
    std::pair<int,Real> val;
    for (int i = 0; i < size; ++i) {
      inBuffer >> val;
      fixed_.insert(val);
    }
    //fixed_ << inBuffer;
    appUnpack(inBuffer);
  }

  virtual void appPack(pebbl::PackBuffer & outBuffer) {}
  virtual void appUnpack(pebbl::UnPackBuffer & inBuffer) {}
#endif

  void setRootComputation() { fixed_.clear(); }

  void boundComputation(double* controlParam) {
#ifdef HAVE_MPI
    branching_->alertBound();
    if (branching_->teamComm.myRank()==0) {
      pebbl::PackBuffer outBuffer;
      pack(outBuffer);
      int bufferSize = outBuffer.size();
      branching_->teamComm.broadcast(&bufferSize,
                                     1,
                                     MPI_INT,
                                     0);
      branching_->teamComm.broadcast((void*) outBuffer.buf(),
                                     outBuffer.size(),
                                     MPI_PACKED,
                                     0);
    }
    else {
      int bufferSize(0);
      branching_->teamComm.broadcast(&bufferSize,
                                     1,
                                     MPI_INT,
                                     0);
      pebbl::UnPackBuffer inBuffer(bufferSize);
      branching_->teamComm.broadcast((void*) inBuffer.buf(),
                                     bufferSize,
                                     MPI_PACKED,
                                     0);
      unpack(inBuffer);
    }
    initialize();
#endif
*outStream_ << "  Fixed Components " << std::endl << fixed_ << std::endl
            << branching_->teamComm.myRank()   << "  " << branching_->teamComm.mySize()   << "  "
            << branching_->passedComm.myRank() << "  " << branching_->passedComm.mySize() << std::endl;
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
//*outStream_ << "Before getIntialGuess" << std::endl;
//    branching_->getInitialGuess(*problem0_->getPrimalOptimizationVector());
//*outStream_ << "After getIntialGuess" << std:endl;
    problem0_->getPrimalOptimizationVector()->zero();

    // Construct optimization solver
    solver_ = makePtr<Solver<Real>>(problem0_,*parlist);
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
    branching_->setInitialGuess(*solution_);
    setState(bounded);
#ifdef HAVE_MPI
    if (!candidateSolution()) {
      Real IV = incumbentValue();
      branching_->teamComm.broadcast(&IV,
                                     1,
                                     MPI_DOUBLE,
                                     0);

      branching_->incumbentValue = IV;
      incumbentHeuristic();
    }
    //splitComputation();
    splitProblem();
#endif
  }

  virtual pebbl::branchSub* makeChild(int whichChild = anyChild) {
    ROL_TEST_FOR_EXCEPTION(whichChild==anyChild,std::logic_error,
      ">>> ROL_PEBBL_BranchSub::makeChild: whichChild is equal to anyChild!");
    BranchSub<Real>* child = new BranchSub<Real>(*this);
    child->updateFixed(index_, (whichChild==0 ? static_cast<Real>(1) : static_cast<Real>(0)));
    return child;
  }

#ifdef HAVE_MPI
  virtual pebbl::parallelBranchSub* makeParallelChild(int whichChild) {
    ROL_TEST_FOR_EXCEPTION(whichChild==anyChild,std::logic_error,
      ">>> ROL_PEBBL_BranchSub::makeParallelChild: whichChild is equal to anyChild!");
    BranchSub<Real>* child = new BranchSub<Real>(*this);
    child->updateFixed(index_, (whichChild==0 ? static_cast<Real>(1) : static_cast<Real>(0)));
*outStream_ << "In makeParallelChild: whichChild = " << whichChild << std::endl;
    return child;
  }
#endif

  bool candidateSolution() { return (nfrac_==0); }

  virtual void incumbentHeuristic() {
    Real tol(std::sqrt(ROL_EPSILON<Real>()));
    rndSolution_->set(*solution_);
    getIntegerVector(rndSolution_)->applyUnary(rnd);
    problem0_->getObjective()->update(*rndSolution_,UpdateType::Temp);
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
