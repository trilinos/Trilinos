// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SAMPLEDSCALAR_H
#define ROL_SAMPLEDSCALAR_H

namespace ROL {

template <class Real, class Key=std::vector<Real>>
class SampledScalar {
private:
  // Storage
  std::map<Key, int> indices_;
  std::vector<bool>  flags_;
  std::vector<Real>  scalars_;
  int maxIndex_;

  // Update flags
  bool updated_;

  void reset(const bool flag = true) {
    if ( flag ) {
      typename std::map<Key, int>::iterator it;
      for (it = indices_.begin(); it != indices_.end(); ++it) {
        flags_[it->second] = false;
      }
    }
  }

public:
  /** \brief Constructor.
  */
  SampledScalar(void)
    : maxIndex_(0), updated_(false) { 
    indices_.clear();
    flags_.clear();
    scalars_.clear();
  }

  /** \brief Update for SampledScalar storage.
  */
  void update(const bool flag = true) {
    updated_ = flag;
    reset(flag);
  }

  /** \brief Return vector corresponding to input parameter.
  */
  bool get(Real &x, const Key &param) {
    int count = indices_.count(param);
    bool flag = false;
    int index = maxIndex_;
    if (count) {
      typename std::map<Key, int>::iterator it = indices_.find(param);
      index = it->second;
      flag  = flags_[index];
      if (flag) {
        x = scalars_[index];
      }
    }
    else {
      indices_.insert(std::pair<Key, int>(param, index));
      flags_.push_back(false);
      scalars_.push_back(static_cast<Real>(0)); 
      maxIndex_++;
    }
    return flag;
  }

  /** \brief Set vector corresponding to input parameter.
  */
  void set(const Real &x, const Key &param) {
    int count = indices_.count(param);
    int index = maxIndex_;
    if (count) {
      typename std::map<Key, int>::iterator it = indices_.find(param);
      index = it->second;
      flags_[index] = true;
      scalars_[index] = x;
    }
    else {
      indices_.insert(std::pair<Key, int>(param, index));
      flags_.push_back(true);
      scalars_.push_back(x); 
      maxIndex_++;
    }
  }
}; // class SampledScalar

} // namespace ROL

#endif
