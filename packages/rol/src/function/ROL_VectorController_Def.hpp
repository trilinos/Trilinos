// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SIMCONTROLLER_DEF_H
#define ROL_SIMCONTROLLER_DEF_H

namespace ROL {

template <class Real, class Key>
VectorController<Real,Key>::VectorController(void)
  : maxIndex_(0), maxIndex_trial_(0), maxIndex_temp_(0),
    trial_(false), temp_(false), objUpdated_(false), conUpdated_(false) { 
  indices_.clear(); indices_trial_.clear(); indices_temp_.clear();
  flags_.clear();   flags_trial_.clear();   flags_temp_.clear();
  vectors_.clear(); vectors_trial_.clear(); vectors_temp_.clear();
}

template <class Real, class Key>
void VectorController<Real,Key>::reset(bool flag) {
  if ( flag ) {
    for (auto it = indices_.begin(); it != indices_.end(); ++it) {
      flags_[it->second] = false;
    }
  }
}

template <class Real, class Key>
void VectorController<Real,Key>::objectiveUpdate(bool flag) {
  temp_  = false;
  trial_ = false;
  if (!conUpdated_) {
    reset(flag);
  }
  objUpdated_ = true;
  if (conUpdated_ && objUpdated_) {
    objUpdated_ = false;
    conUpdated_ = false;
  }
}

template <class Real, class Key>
void VectorController<Real,Key>::constraintUpdate(bool flag) {
  temp_  = false;
  trial_ = false;
  if (!objUpdated_) {
    reset(flag);
  }
  conUpdated_ = true;
  if (conUpdated_ && objUpdated_) {
    objUpdated_ = false;
    conUpdated_ = false;
  }
}

template <class Real, class Key>
void VectorController<Real,Key>::objectiveUpdate(UpdateType type) {
  if (type == UpdateType::Temp) {
    temp_       = true;
    trial_      = false;
    objUpdated_ = false;
    conUpdated_ = false;
    resetTemp();
  }
  else {
    if (!conUpdated_) {
      switch(type) {
        case UpdateType::Initial: temp_ = false; trial_ = false; reset(true);  break;
        case UpdateType::Trial:   temp_ = false; trial_ = true;  resetTrial(); break;
        case UpdateType::Accept:  temp_ = false; trial_ = false; accept();     break;
        case UpdateType::Revert:  temp_ = false; trial_ = false;               break;
        case UpdateType::Temp:    temp_ = true;  trial_ = false; resetTemp();  break;
      }
    }
    objUpdated_ = true;
    if (conUpdated_ && objUpdated_) {
      objUpdated_ = false;
      conUpdated_ = false;
    }
  }
}

template <class Real, class Key>
void VectorController<Real,Key>::constraintUpdate(UpdateType type) {
  if (type == UpdateType::Temp) {
    temp_       = true;
    trial_      = false;
    objUpdated_ = false;
    conUpdated_ = false;
    resetTemp();
  }
  else {
    if (!objUpdated_) {
      switch(type) {
        case UpdateType::Initial: temp_ = false; trial_ = false; reset(true);  break;
        case UpdateType::Trial:   temp_ = false; trial_ = true;  resetTrial(); break;
        case UpdateType::Accept:  temp_ = false; trial_ = false; accept();     break;
        case UpdateType::Revert:  temp_ = false; trial_ = false;               break;
        case UpdateType::Temp:    temp_ = true;  trial_ = false; resetTemp();  break;
      }
    }
    conUpdated_ = true;
    if (conUpdated_ && objUpdated_) {
      objUpdated_ = false;
      conUpdated_ = false;
    }
  }
}

template <class Real, class Key>
bool VectorController<Real,Key>::isNull(const Key &param) const {
  if (!temp_) {
    if (trial_) {
      return isNull(param,indices_trial_);
    }
    else {
      return isNull(param,indices_);
    }
  }
  else {
    return isNull(param,indices_temp_);
  }
  return false;
}

template <class Real, class Key>
bool VectorController<Real,Key>::isComputed(const Key &param) const {
  if (!temp_) {
    if (trial_) {
      return isComputed(param,indices_trial_,flags_trial_);
    }
    else {
      return isComputed(param,indices_,flags_);
    }
  }
  else {
    return isComputed(param,indices_temp_,flags_temp_);
  }
  return false;
}

template <class Real, class Key>
void VectorController<Real,Key>::allocate(const Vector<Real> &x, const Key &param) {
  if (!temp_) {
    if (trial_) {
      allocate(x,param,indices_trial_,flags_trial_,vectors_trial_,maxIndex_trial_);
    }
    else {
      allocate(x,param,indices_,flags_,vectors_,maxIndex_);
    }
  }
  else {
    allocate(x,param,indices_temp_,flags_temp_,vectors_temp_,maxIndex_temp_);
  }
}

template <class Real, class Key>
const Ptr<const Vector<Real>> VectorController<Real,Key>::get(const Key &param) const {
  if (!temp_) {
    if (trial_) {
      return get(param,indices_trial_,flags_trial_,vectors_trial_,maxIndex_trial_);
    }
    else {
      return get(param,indices_,flags_,vectors_,maxIndex_);
    }
  }
  else {
    return get(param,indices_temp_,flags_temp_,vectors_temp_,maxIndex_temp_);
  }
  return nullPtr;
}

template <class Real, class Key>
const Ptr<Vector<Real>> VectorController<Real,Key>::set(const Key &param) {
  if (!temp_) {
    if (trial_) {
      return set(param,indices_trial_,flags_trial_,vectors_trial_,maxIndex_trial_);
    }
    else {
      return set(param,indices_,flags_,vectors_,maxIndex_);
    }
  }
  else {
    return set(param,indices_temp_,flags_temp_,vectors_temp_,maxIndex_temp_);
  }
  return nullPtr;
}

template <class Real, class Key>
bool VectorController<Real,Key>::get(Vector<Real> &x, const Key &param) {
  bool flag = false;
  if (!temp_) {
    if (trial_) {
      flag = get(x,param,indices_trial_,flags_trial_,vectors_trial_,maxIndex_trial_);
    }
    else {
      flag = get(x,param,indices_,flags_,vectors_,maxIndex_);
    }
  }
  else {
    flag = get(x,param,indices_temp_,flags_temp_,vectors_temp_,maxIndex_temp_);
  }
  return flag;
}

template <class Real, class Key>
void VectorController<Real,Key>::set(const Vector<Real> &x, const Key &param) {
  if (!temp_) {
    if (trial_) {
      set(x,param,indices_trial_,flags_trial_,vectors_trial_,maxIndex_trial_);
    }
    else {
      set(x,param,indices_,flags_,vectors_,maxIndex_);
    }
  }
  else {
    set(x,param,indices_temp_,flags_temp_,vectors_temp_,maxIndex_temp_);
  }
}

template <class Real, class Key>
void VectorController<Real,Key>::push(VectorController<Real,Key> &to) const {
  for (auto it = indices_.begin(); it != indices_.end(); ++it) {
    to.set(*vectors_[it->second],it->first);
  }
}

template <class Real, class Key>
void VectorController<Real,Key>::resetTrial(void) {
  for (auto it = indices_trial_.begin(); it != indices_trial_.end(); ++it) {
    flags_trial_[it->second] = false;
  }
}

template <class Real, class Key>
void VectorController<Real,Key>::resetTemp(void) {
  for (auto it = indices_temp_.begin(); it != indices_temp_.end(); ++it) {
    flags_temp_[it->second] = false;
  }
}

template <class Real, class Key>
bool VectorController<Real,Key>::isNull(const Key &param,
         const std::map<Key,int> &indices) const {
  return (indices.count(param)==0);
}

template <class Real, class Key>
bool VectorController<Real,Key>::isComputed(const Key &param,
         const std::map<Key,int> &indices, const std::vector<bool> &flags) const {
  if (indices.count(param)>0) {
    auto it = indices.find(param);
    return flags[it->second];
  }
  return false;
}

template <class Real, class Key>
void VectorController<Real,Key>::allocate(const Vector<Real> &x, const Key &param,
         std::map<Key,int> &indices, std::vector<bool> &flags,
         std::vector<Ptr<Vector<Real>>> &vectors, int &maxIndex) const {
  if (isNull(param,indices)) {
    int index = maxIndex;
    indices.insert(std::pair<Key,int>(param, index));
    flags.push_back(false);
    vectors.push_back(x.clone()); 
    maxIndex++;
  }
}

template <class Real, class Key>
const Ptr<const Vector<Real>> VectorController<Real,Key>::get(const Key &param,
         const std::map<Key,int> &indices, const std::vector<bool> &flags,
         const std::vector<Ptr<Vector<Real>>> &vectors, const int &maxIndex) const {
  int count = indices.count(param);
  int index = maxIndex;
  if (count) {
    auto it = indices.find(param);
    index = it->second;
    return vectors[index];
  }
  return nullPtr;
}

template <class Real, class Key>
const Ptr<Vector<Real>> VectorController<Real,Key>::set(const Key &param,
         std::map<Key,int> &indices, std::vector<bool> &flags,
         std::vector<Ptr<Vector<Real>>> &vectors, int &maxIndex) const {
  ROL_TEST_FOR_EXCEPTION(isNull(param,indices),std::logic_error,
    ">>> ROL::VectorController::set : Vector has not been allocated!");
  ROL_TEST_FOR_EXCEPTION(isComputed(param,indices,flags),std::logic_error,
    ">>> ROL::VectorController::set : Vector is already computed!");
  int count = indices.count(param);
  int index = maxIndex;
  if (count) {
    auto it = indices.find(param);
    index = it->second;
    flags[index] = true;
    return vectors[index];
  }
  return nullPtr;
}

template <class Real, class Key>
bool VectorController<Real,Key>::get(Vector<Real> &x, const Key &param,
         std::map<Key,int> &indices, std::vector<bool> &flags,
         std::vector<Ptr<Vector<Real>>> &vectors, int &maxIndex) const {
  int count = indices.count(param);
  bool flag = false;
  int index = maxIndex;
  if (count) {
    auto it = indices.find(param);
    index = it->second;
    flag  = flags[index];
    if (flag) {
      x.set(*vectors[index]);
    }
  }
  else {
    indices.insert(std::pair<Key,int>(param, index));
    flags.push_back(false);
    vectors.push_back(x.clone()); 
    maxIndex++;
  }
  return flag;
}

template <class Real, class Key>
void VectorController<Real,Key>::set(const Vector<Real> &x,const Key &param,
         std::map<Key,int> &indices, std::vector<bool> &flags,
         std::vector<Ptr<Vector<Real>>> &vectors, int &maxIndex) const {
  int count = indices.count(param);
  int index = maxIndex;
  if (count) {
    auto it = indices.find(param);
    index = it->second;
    flags[index] = true;
    vectors[index]->set(x);
  }
  else {
    indices.insert(std::pair<Key,int>(param, index));
    flags.push_back(true);
    vectors.push_back(x.clone()); 
    vectors[index]->set(x);
    maxIndex++;
  }
}

template <class Real, class Key>
void VectorController<Real,Key>::accept(void) {
  reset(true);
  for (auto it = indices_trial_.begin(); it != indices_trial_.end(); ++it) {
    set(*vectors_trial_[it->second],it->first,indices_,flags_,vectors_,maxIndex_);
  }
}

} // namespace ROL

#endif
