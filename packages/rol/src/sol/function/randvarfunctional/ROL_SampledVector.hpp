// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SAMPLEDVECTOR_H
#define ROL_SAMPLEDVECTOR_H

namespace ROL {

template <class Real, class Key=std::vector<Real>>
class SampledVector {
private:
  // Storage
  std::map<Key, int>                  indices_;
  std::vector<bool>                   flags_;
  std::vector<ROL::Ptr<Vector<Real>>> vectors_;
  int maxIndex_;

  // Update flags
  bool updated_;

  void reset(const bool flag = true) {
    if ( flag ) {
      flags_.assign(flags_.size(),false);
//      typename std::map<Key, int>::iterator it;
//      for (it = indices_.begin(); it != indices_.end(); ++it) {
//        flags_[it->second] = false;
//      }
    }
  }

public:
  /** \brief Constructor.
  */
  SampledVector(void)
    : maxIndex_(0), updated_(false) { 
    indices_.clear();
    flags_.clear();
    vectors_.clear();
  }

  /** \brief Update for SampledVector storage.
  */
  void update(const bool flag = true) {
    updated_ = flag;
    reset(flag);
  }

  /** \brief Return vector corresponding to input parameter.
  */
  bool get(Vector<Real> &x, const Key &param) {
    int count = indices_.count(param);
    bool flag = false;
    int index = maxIndex_;
    if (count) {
      typename std::map<Key, int>::iterator it = indices_.find(param);
      index = it->second;
      flag  = flags_[index];
      if (flag) {
        x.set(*vectors_[index]);
      }
    }
    else {
      indices_.insert(std::pair<Key, int>(param, index));
      flags_.push_back(false);
      vectors_.push_back(x.clone()); 
      maxIndex_++;
    }
    return flag;
  }

  /** \brief Set vector corresponding to input parameter.
  */
  void set(const Vector<Real> &x, const Key &param) {
    int count = indices_.count(param);
    int index = maxIndex_;
    if (count) {
      typename std::map<Key, int>::iterator it = indices_.find(param);
      index = it->second;
      flags_[index] = true;
      vectors_[index]->set(x);
    }
    else {
      indices_.insert(std::pair<Key, int>(param, index));
      flags_.push_back(true);
      vectors_.push_back(x.clone()); 
      vectors_[index]->set(x);
      maxIndex_++;
    }
  }
}; // class SampledVector

} // namespace ROL

#endif
