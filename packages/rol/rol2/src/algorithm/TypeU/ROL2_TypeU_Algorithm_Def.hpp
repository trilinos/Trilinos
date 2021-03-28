#pragma once
#ifndef ROL2_TYPEU_ALGORITHM_DEF_HPP
#define ROL2_TYPEU_ALGORITHM_DEF_HPP

namespace ROL2 {
namespace TypeU {

//------------------------------------------------------------------------- 
// Algorithm::State method definitions

template<class Real>
bool Algorithm<Real>::State::is_initialized() const {
  return ROL2::Algorithm<Real>::State::is_initialized() &&
         stepVec_ != nullPtr &&
         gradientVec_ != nullPtr;
}

template<class Real>
void Algorithm<Real>::State::reset() {
//  ROL_TEST_FOR_EXCEPTION( !is_initialized(), std::logic_error, "Storage vector memory must be allocated" );
  ROL2::Algorithm<Real>::State::reset();
  stepVec_->zero();
  gradientVec_->zero();
}

//------------------------------------------------------------------------- 
// Algorithm method definitions

template<class Real>
Algorithm<Real>::Algorithm() 
  : state_(makePtr<State>()),
    status_(makePtr<CombinedStatusTest<Real>>()) {}

template<class Real>
void Algorithm<Real>::initialize( const Vector<Real>& x,
                                  const Vector<Real>& g ) {
  ROL2::Algorithm<Real>::initialize(x);
  auto& state = getState();
  state.stepVec_     = x.clone();
  state.gradientVec_ = g.clone();
}


template<class Real>
void Algorithm<Real>::setStatusTest( const Ptr<StatusTest<Real>>& status, 
                                           bool                   appendStatus ) {
  if( appendStatus ) status_->add(status);
  else status_->reset();
}

//template<typename Real>
//void Algorithm<Real>::run( OptimizationProblem<Real>& problem,
//                              std::ostream&              os ) {
//
//  ROL_TEST_FOR_EXCEPTION( problem.getProblemType() != TYPE_U, std::logic_error, "Incompatible Problem Type" ) {
//
//  auto output = run( *problem.getPrimalOptimizationVector(),
//                     *problem.getDualOptimizationVector(),
//                     *problem.getObjective(),
//                     os);
//
//    problem.finalizeIteration();
//    return output;
//  }
//}

template<typename Real>
void Algorithm<Real>::run( Vector<Real>&    x,
                           Objective<Real>& obj,
                           std::ostream&    os ) {
  return run(x,x.dual(),obj,os);
}

//template<typename Real>
//void Algorithm<Real>::run( Vector<Real>&       x,
//                           Objective<Real>&  obj,
//                           Constraint<Real>& linear_con,
//                           Vector<Real>&     linear_mul,
//                           std::ostream&     os ) {
//  return run(x,x.dual(),obj,linear_con,linear_mul,linear_mul.dual(),os);
//}

//template<typename Real>
//void Algorithm<Real>::run(       Vector<Real>&     x,
//                           const Vector<Real>&     g,
//                                 Objective<Real>&  obj,
//                                 Constraint<Real>& linear_con,
//                                 Vector<Real>&     linear_mul,
//                           const Vector<Real>&     linear_c,
//                                 std::ostream&     os ) {
//  auto xfeas = x.clone(); 
//  xfeas->set(x);
//  ReduceLinearConstraint<Real> rlc(makePtrFromRef(linear_con),xfeas,makePtrFromRef(linear_c));
//  auto s = x.clone(); 
//  s->zero();
//  auto output = run(*s,g,*rlc.transform(makePtrFromRef(obj)),os);
//  rlc.project(x,*s);
//  x.plus(*rlc.getFeasibleVector());
//  return output;
//}
//
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
}

template<typename Real>
void Algorithm<Real>::writeName( std::ostream& ) const {
//  ROL_TEST_FOR_EXCEPTION(
//  throw Exception::NotImplemented(">>> ROL::Algorithm::writeName() is not implemented!");
}

template<typename Real>
void Algorithm<Real>::writeOutput( std::ostream& os, bool print_header ) const {

  os << std::scientific << std::setprecision(6);

  if ( print_header ) writeHeader(os);

  os << "  ";
  os << std::setw(6)  << std::left << state_->iter_;
  os << std::setw(15) << std::left << state_->value_;
  os << std::setw(15) << std::left << state_->gnorm_;

  if ( state_->iter_ == 0 ) {
    os << std::setw(15) << std::left << state_->snorm_;
    os << std::setw(10) << std::left << state_->nfval_;
    os << std::setw(10) << std::left << state_->ngrad_;
  }

  os << std::endl;
}

template<class Real>
void Algorithm<Real>::reset() {
  state_->reset();
}

} // namespace TypeU
} // namespace ROL2

#endif //ROL2_TYPEU_ALGORITHM_DEF_HPP

