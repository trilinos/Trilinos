#pragma once
#ifndef ROL2_ALGORITHM_DEF_HPP
#define ROL2_ALGORITHM_DEF_HPP

namespace ROL2 {

template<class Real>
void Algorithm<Real>::State::reset() {
  iter_       = 0;
  minIter_    = 0;
  nfval_      = 0;
  ngrad_      = 0;
  value_      = ROL_INF<Real>;
  minValue_   = ROL_INF<Real>;
  gnorm_      = ROL_INF<Real>;
  snorm_      = ROL_INF<Real>;
  statusFlag_ = ExitStatus::Last;
  flag_                  = false;
  if( iterateVec_ != nullPtr ) iterateVec_->zero();
  if( minIterVec_ != nullPtr ) minIterVec_->zero();
} // ROL2::Algorithm<Real>::State::reset


template<class Real>
bool Algorithm<Real>::State::is_initialized() const  {
  return iterateVec_ != nullPtr &&
         minIterVec_ != nullPtr;
}

template<class Real>
void Algorithm<Real>::initialize( const Vector<Real>& x ) {
  auto& state = getState();
  if (state.iterateVec_ == nullPtr) {
    state.iterateVec_ = x.clone();
  }
  state.iterateVec_->set(x);
  if (state.minIterVec_ == nullPtr) {
    state.minIterVec_ = x.clone();
  }
  state.minIterVec_->set(x);
  state.minIter_ = state.iter_;
  state.minValue_ = state.value_;
} 

template<class Real>
void Algorithm<Real>::reset() {
  auto& state = getState();
  state.reset();
}

} // namesapce ROL2

#endif // ROL2_ALGORITHM_DEF_HPP

