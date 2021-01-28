#pragma once
#ifndef ROL2_TYPEU_ALGORITHM_DEF_HPP
#define ROL2_TYPEU_ALGORITHM_DEF_HPP

namespace ROL2 {
namespace TypeU {

//------------------------------------------------------------------------- 
// Algorithm::State method definitions

template<class Real>
void Algorithm<Real>::State::initialize( const Vector<Real>& x ) const {

  iterateVec_  = x.clone();
  minIterVec_  = x.clone();
  stepVec_     = x.clone();
  gradientVec_ = x.dual().clone();

  is_initialized_ = true;
}

template<class Real>
bool Algorithm<Real>::State::is_initialized() const {
  return is_initialized_;
}

template<class Real>
void Algorithm<Real>::State::reset() {
  ROL_TEST_FOR_EXCEPTION( !is_initialized(), std::logic_error, "Storage vector memory must be allocated" );

  iter_     = 0;
  minIter_  = 0;
  nfval_    = 0;
  ngrad_    = 0;
  value_    = ROL_INF<Real>;
  minValue_ = ROL_INF<Real>;
  gnorm_    = ROL_INF<Real>;
  snorm_    = ROL_INF<Real>;

  exitStatus_ = ExitStatus::Last;

  iterateVec_->zero();
  minIteVec_->zero();
  stepVec_->zero();
  gradientVec_->zero();
}

//------------------------------------------------------------------------- 
// Algorithm method definitions

template<class Real>
void Algorithm<Real>::setStatusTest( const Ptr<StatusTest<Real>& status, 
                                       bool appendStatus ) {
  if( appendStatus ) status_->add(status);
  else status_->reset();
}

template<typename Real>
void Algorithm<Real>::run( OptimizationProblem<Real>& problem,
                              std::ostream&              outStream ) {

  ROL_TEST_FOR_EXCEPTION( problem.getProblemType() != TYPE_U, std::logic_error, "Incompatible Problem Type" ) {

  auto output = run( *problem.getPrimalOptimizationVector(),
                     *problem.getDualOptimizationVector(),
                     *problem.getObjective(),
                     outStream);

    problem.finalizeIteration();
    return output;
  }
}

template<typename Real>
void Algorithm<Real>::run( Vector<Real>&    x,
                             Objective<Real>& obj,
                             std::ostream&    outStream ) {
  return run(x,x.dual(),obj,outStream);
}

template<typename Real>
void Algorithm<Real>::run( Vector<Real>&     x,
                             Objective<Real>&  obj,
                             Constraint<Real>& linear_con,
                             Vector<Real>&     linear_mul,
                             std::ostream&     outStream ) {
  return run(x,x.dual(),obj,linear_con,linear_mul,linear_mul.dual(),outStream);
}

template<typename Real>
void Algorithm<Real>::run( Vector<Real>&       x,
                             const Vector<Real>& g,
                             Objective<Real>&    obj,
                             Constraint<Real>&   linear_con,
                             Vector<Real>&       linear_mul,
                             const Vector<Real>& linear_c,
                             std::ostream&       outStream ) {
  auto xfeas = x.clone(); 
  xfeas->set(x);
  ReduceLinearConstraint<Real> rlc(makePtrFromRef(linear_con),xfeas,makePtrFromRef(linear_c));
  auto s = x.clone(); 
  s->zero();
  auto output = run(*s,g,*rlc.transform(makePtrFromRef(obj)),outStream);
  rlc.project(x,*s);
  x.plus(*rlc.getFeasibleVector());
  return output;
}

template<typename Real>
void Algorithm<Real>::writeHeader( std::ostream& os ) const {
  os << "  ";
  os << std::setw(6)  << std::left << "iter";
  os << std::setw(15) << std::left << "value";
  os << std::setw(15) << std::left << "gnorm";
  os << std::setw(15) << std::left << "snorm";
  os << std::setw(10) << std::left << "#fval";
  os << std::setw(10) << std::left << "#grad";
  os << std::endl;
  return hist.str();
}

template<typename Real>
void Algorithm<Real>::writeName( std::ostream& ) const {
  ROL_TEST_FOR_EXCEPTION(
  throw Exception::NotImplemented(">>> ROL::Algorithm::writeName() is not implemented!");
}

template<typename Real>
void Algorithm<Real>::writeOutput( std::ostream os, const bool print_header ) const {

  os << std::scientific << std::setprecision(6);

  if ( print_header ) writeHeader(os);

    os << "  ";
    os << std::setw(6)  << std::left << state_->iter;
    os << std::setw(15) << std::left << state_->value;
    os << std::setw(15) << std::left << state_->gnorm;

    if ( state_->iter == 0 ) {
      os << std::setw(15) << std::left << state_->snorm;
      os << std::setw(10) << std::left << state_->nfval;
      os << std::setw(10) << std::left << state_->ngrad;
    }

    os << std::endl;
}

template<typename Real>
const Algorithm<Real>::State& 
Algorithm<Real>::getState() const { return *state_; }

template<typename Real>
Algorithm<Real>::State& 
Algorithm<Real>::getState() { return *state_; }

template<typename Real>
const StatusTest<Real>&
const Algorithm<Real>::State& getState() const { return *status_; }

template<typename Real>
StatusTest<Real>&
Algorithm<Real>::State& getState() { return *status_; }


} // namespace TypeU
} // namespace ROL2

#endif //ROL2_TYPEU_ALGORITHM_DEF_HPP

