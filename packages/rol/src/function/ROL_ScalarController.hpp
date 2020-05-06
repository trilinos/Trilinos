// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER


#ifndef ROL_SCALARCONTROLLER_H
#define ROL_SCALARCONTROLLER_H

#include "ROL_UpdateTypes.hpp"

namespace ROL {

template <class Real, class Key=std::vector<Real>>
class ScalarController {
private:
  // Storage
  std::map<Key, int> indices_, indices_trial_, indices_temp_;
  std::vector<bool>  flags_, flags_trial_, flags_temp_;
  std::vector<Real>  scalars_, scalars_trial_, scalars_temp_;
  int maxIndex_, maxIndex_trial_, maxIndex_temp_;

  // Update flags
  bool objUpdated_, conUpdated_;
  bool temp_, trial_;

public:
  /** \brief Constructor.
  */
  ScalarController(void)
    : maxIndex_(0), maxIndex_trial_(0), maxIndex_temp_(0),
      objUpdated_(false), conUpdated_(false),
      temp_(false), trial_(false) { 
    indices_.clear(); indices_trial_.clear(); indices_temp_.clear();
    flags_.clear(); flags_trial_.clear(); flags_temp_.clear();
    scalars_.clear(); scalars_trial_.clear(); scalars_temp_.clear();
  }

  void reset(bool flag = true) {
    if ( flag ) {
      for (auto it = indices_.begin(); it != indices_.end(); ++it) {
        flags_[it->second] = false;
      }
    }
  }

  /** \brief Update for SampledScalar storage.
  */
  void objectiveUpdate(bool flag = true) {
    if (!conUpdated_) {
      reset(flag);
    }
    objUpdated_ = true;
    if (conUpdated_ && objUpdated_) {
      objUpdated_ = false;
      conUpdated_ = false;
    }
  }

  /** \brief Equality constraint update for SimController storage.
  */
  void constraintUpdate(bool flag = true) {
    if (!objUpdated_) {
      reset(flag);
    }
    conUpdated_ = true;
    if (conUpdated_ && objUpdated_) {
      objUpdated_ = false;
      conUpdated_ = false;
    }
  }

  void objectiveUpdate(EUpdateType type) {
    if (!conUpdated_) {
      switch(type) {
        case UPDATE_INITIAL: temp_ = false; trial_ = false; reset(true);  break;
        case UPDATE_TRIAL:   temp_ = false; trial_ = true;  resetTrial(); break;
        case UPDATE_ACCEPT:  temp_ = false; trial_ = false; accept();     break;
        case UPDATE_REVERT:  temp_ = false; trial_ = false;               break;
        case UPDATE_TEMP:    temp_ = true;  trial_ = false; resetTemp();  break;
      }
    }
    objUpdated_ = true;
    if (conUpdated_ && objUpdated_) {
      objUpdated_ = false;
      conUpdated_ = false;
    }
  }

  void constraintUpdate(EUpdateType type) {
    if (!objUpdated_) {
      switch(type) {
        case UPDATE_INITIAL: temp_ = false; trial_ = false; reset(true);  break;
        case UPDATE_TRIAL:   temp_ = false; trial_ = true;  resetTrial(); break;
        case UPDATE_ACCEPT:  temp_ = false; trial_ = false; accept();     break;
        case UPDATE_REVERT:  temp_ = false; trial_ = false;               break;
        case UPDATE_TEMP:    temp_ = true;  trial_ = true;  resetTemp();  break;
      }
    }
    conUpdated_ = true;
    if (conUpdated_ && objUpdated_) {
      objUpdated_ = false;
      conUpdated_ = false;
    }
  }

  /** \brief Return vector corresponding to input parameter.
  */
  bool get(Real &x, const Key &param) {
    bool flag = false;
    if (!temp_) {
      if (trial_) {
        flag = get(x,param,indices_trial_,flags_trial_,scalars_trial_,maxIndex_trial_);
      }
      else {
        flag = get(x,param,indices_,flags_,scalars_,maxIndex_);
      }
    }
    else {
      flag = get(x,param,indices_temp_,flags_temp_,scalars_temp_,maxIndex_temp_);
    }
    return flag;
  }

  /** \brief Set vector corresponding to input parameter.
  */
  void set(const Real &x, const Key &param) {
    if (!temp_) {
      if (trial_) {
        set(x,param,indices_trial_,flags_trial_,scalars_trial_,maxIndex_trial_);
      }
      else {
        set(x,param,indices_,flags_,scalars_,maxIndex_);
      }
    }
    else {
      set(x,param,indices_temp_,flags_temp_,scalars_temp_,maxIndex_temp_);
    }
  }

private:

  void resetTrial(void) {
    for (auto it = indices_trial_.begin(); it != indices_trial_.end(); ++it) {
      flags_trial_[it->second] = false;
    }
  }

  void resetTemp(void) {
    for (auto it = indices_temp_.begin(); it != indices_temp_.end(); ++it) {
      flags_temp_[it->second] = false;
    }
  }

  bool get(Real &x, const Key &param, std::map<Key,int> &indices,
           std::vector<bool> &flags, std::vector<Real> &scalars, int &maxIndex) {
    int count = indices.count(param);
    bool flag = false;
    int index = maxIndex;
    if (count) {
      auto it = indices.find(param);
      index = it->second;
      flag  = flags[index];
      if (flag) {
        x = scalars[index];
      }
    }
    else {
      indices.insert(std::pair<Key, int>(param, index));
      flags.push_back(false);
      scalars.push_back(static_cast<Real>(0)); 
      maxIndex++;
    }
    return flag;
  }

  void set(const Real &x, const Key &param, std::map<Key,int> &indices,
           std::vector<bool> &flags, std::vector<Real> &scalars, int &maxIndex) {
    int count = indices.count(param);
    int index = maxIndex;
    if (count) {
      auto it = indices.find(param);
      index = it->second;
      flags[index] = true;
      scalars[index] = x;
    }
    else {
      indices.insert(std::pair<Key, int>(param, index));
      flags.push_back(true);
      scalars.push_back(x); 
      maxIndex++;
    }
  }

  void accept(void) {
    reset(true);
    for (auto it = indices_trial_.begin(); it != indices_trial_.end(); ++it) {
      set(scalars_trial_[it->second],it->first,indices_,flags_,scalars_,maxIndex_);
    }
  }
}; // class SampledScalar

} // namespace ROL

#endif
