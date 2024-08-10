// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SOLVER_DEF_HPP
#define ROL_SOLVER_DEF_HPP

namespace ROL {

template<typename Real>
Solver<Real>::Solver( const Ptr<Problem<Real>> &opt,
                      ParameterList            &parlist,
                      const Ptr<Secant<Real>>  &secant )
  : opt_(opt), problemType_(opt_->getProblemType()) {
  switch (problemType_) {
    case TYPE_U:  algoU_ = TypeU::AlgorithmFactory<Real>(parlist,secant); break;
    case TYPE_B:  algoB_ = TypeB::AlgorithmFactory<Real>(parlist,secant); break;
    case TYPE_E:  algoE_ = TypeE::AlgorithmFactory<Real>(parlist,secant); break;
    case TYPE_EB: algoG_ = TypeG::AlgorithmFactory<Real>(parlist,secant); break;
    case TYPE_LAST:
      ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
        "Error in Solver::solve() : Unsupported problem type");
  }
}

template<typename Real>
int Solver<Real>::solve( const Ptr<StatusTest<Real>> &status,
                         bool combineStatus) {
  nullstream bhs;
  return solve(bhs,status,combineStatus);
}

template<typename Real>
int Solver<Real>::solve( std::ostream &outStream,
                         const Ptr<StatusTest<Real>> &status,
                         bool combineStatus ) {
  switch (problemType_) {
    case TYPE_U:
      if (status != nullPtr) algoU_->setStatusTest(status,combineStatus);
      algoU_->run(*opt_,outStream);
      break;
    case TYPE_B:
      if (status != nullPtr) algoB_->setStatusTest(status,combineStatus);
      algoB_->run(*opt_,outStream);
      break;
    case TYPE_E:
      if (status != nullPtr) algoE_->setStatusTest(status,combineStatus);
      algoE_->run(*opt_,outStream);
      break;
    case TYPE_EB:
      if (status != nullPtr) algoG_->setStatusTest(status,combineStatus);
      algoG_->run(*opt_,outStream);
      break;
    case TYPE_LAST:
      ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
        "Error in Solver::solve() : Unsupported problem type");
  }
  // TODO: Interrogate AlgorithmState and StatusTest to generate a return code
  //       that indicates why the solver has stopped

  // Return an integer code
  return 0;
}

template<typename Real>
Ptr<const AlgorithmState<Real>> Solver<Real>::getAlgorithmState() const {
//Ptr<const AlgorithmState<Real>>& Solver<Real>::getAlgorithmState() const {
  switch (problemType_) {
    case TYPE_U:  return algoU_->getState();
    case TYPE_B:  return algoB_->getState();
    case TYPE_E:  return algoE_->getState();
    case TYPE_EB: return algoG_->getState();
    case TYPE_LAST:
    default:
      ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
        "Error in Solver::getAlgorithmState() : Unsupported problem type");
  }
}

template<typename Real>
void Solver<Real>::reset() {
  switch (problemType_) {
    case TYPE_U:  algoU_->reset(); break;
    case TYPE_B:  algoB_->reset(); break;
    case TYPE_E:  algoE_->reset(); break;
    case TYPE_EB: algoG_->reset(); break;
    case TYPE_LAST:
    default:
      ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
        "Error in Solver::reset() : Unsupported problem type");
  }
}

} // namespace ROL

#endif // ROL_SOLVER_DEF_HPP


