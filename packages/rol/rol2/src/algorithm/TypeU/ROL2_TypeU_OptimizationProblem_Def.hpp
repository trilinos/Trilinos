#pragma once
#ifndef ROL2_TYPEU_OPTIMIZATIONPROBLEM_DEF_HPP
#define ROL2_TYPEU_OPTIMIZATIONPROBLEM_DEF_HPP

namespace ROL2 {
namespace TypeU {

template<class Real>
OptimizationProblem<Real>::
OptimizationProblem( const Ptr<Objective<Real>>& obj,
                     const Ptr<Vector<Real>>&    x ) 
  : obj_(obj), sol_(x) {} 

template<class Real>
const Ptr<Objective<Real>>
OptimizationProblem<Real>::getObjective() { return obj_; }

template<class Real>
const Ptr<Vector<Real>>
OptimizationProblem<Real>::getSolutionVector() { return sol_; }

template<class Real>
void OptimizationProblem<Real>::
checkObjective( OptimizationProblem<Real>::CheckData& data,
                Vector<Real>&                         x,
                Vector<Real>&                         u,
                Vector<Real>&                         v,
                std::ostream&                         os,
                int                                   numSteps,
                int                                   fd_order ) {
  if( !obj_.is_nullPtr() ) {
      os << std::endl << "Performing OptimizationProblem diagnostics."
                << std::endl << std::endl;
      os << "Checking objective function." << std::endl;
      data.checkGradient = obj_->checkGradient(x,v,true,os,numSteps,order);
      os << std::endl;
      data.checkHessVec  = obj_->checkHessVec(x,u,true,os,numSteps,order);
      os << std::endl;
      data.checkHessSym  = obj_->checkHessSym(x,u,v,true,os);
      os << std::endl;
  }
}

template<class Real>
void OptimizationProblem<Real>::
checkObjective( Vector<Real>& x,
                Vector<Real>& u,
                Vector<Real>& v,
                std::ostream& os,
                int           numSteps,
                int           fd_order ) {
  CheckData data;
  checkObjective(data,x,u,v,os,numSteps,fd_order);
}

template<class Real>
void OptimizationProblem<Real>::
checkSolutionVector( OptimizationProblem<Real>::CheckData& data,
                     Vector<Real>&                         x,
                     Vector<Real>&                         y,
                     Vector<Real>&                         v,
                     std::ostream&                         os ) {
  if( !obj_.is_nullPtr() ) {
    os << "\nPerforming OptimizationProblem diagnostics." << std::endl << std::endl;
    os << "Checking vector operations in optimization vector space X." << std::endl;
    data.checkSolutionVector = x.checkVector(y,u,true,os);
  }
}

template<class Real>
void OptimizationProblem<Real>::
checkSolutionVector( Vector<Real>& x,
                     Vector<Real>& y,
                     Vector<Real>& v,
                     std::ostream& os ) {
  CheckData data;
  if( !obj_.is_nullPtr() ) {
    os << "\nPerforming OptimizationProblem diagnostics." << std::endl << std::endl;
    os << "Checking vector operations in optimization vector space X." << std::endl;
    data.checkSolutionVector = x.checkVector(y,u,true,os);
  }
}




} // namespace TypeU
} // namespace ROL2

#endif //ROL2_TYPEU_OPTIMIZATIONPROBLEM_DEF_HPP

