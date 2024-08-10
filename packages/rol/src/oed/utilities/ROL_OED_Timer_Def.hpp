// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_TIMER_DEF_HPP
#define ROL_OED_TIMER_DEF_HPP

namespace ROL {
namespace OED {

template<typename Real, typename Key>
Timer<Real,Key>::Timer(const std::string name) : name_(name) {
  reset();
}

template<typename Real, typename Key>
void Timer<Real,Key>::rename(std::string name) {
  name_ = name;
}

template<typename Real, typename Key>
void Timer<Real,Key>::reset() {
  timer_.clear(); time_.clear(); count_.clear();
}

template<typename Real, typename Key>
void Timer<Real,Key>::start(Key i) {
  if (timer_.count(i) == 0) {
    timer_.insert(std::pair<Key,std::clock_t>(i,std::clock()));
    time_.insert(std::pair<Key,Real>(i,static_cast<Real>(0)));
    count_.insert(std::pair<Key,int>(i,0));
  }
  else {
    timer_[i] = std::clock();
  }
}

template<typename Real, typename Key>
void Timer<Real,Key>::stop(Key i) {
  if (timer_.count(i) == 0) {
    throw Exception::NotImplemented(">>> OED::Timer::stop : Key does not exist!");
  }
  time_[i] += static_cast<Real>(std::clock()-timer_[i])/static_cast<Real>(CLOCKS_PER_SEC);
  count_[i]++;
}

template<typename Real, typename Key>
void Timer<Real,Key>::summarize(std::ostream &stream,
         const Ptr<BatchManager<Real>> &bman) const {
  Real nbatch(1);
  if (bman != nullPtr) {
    nbatch = static_cast<Real>(bman->numBatches());
    std::vector<Real> mv;
    for (typename std::map<Key,int>::iterator it = count_.begin(); it != count_.end(); ++it) {
      mv.push_back(static_cast<Real>(count_[it->first]));
      mv.push_back(time_[it->first]);
    }
    int size = mv.size(), cnt = 0;
    std::vector<Real> gv(size,0);
    bman->sumAll(&mv[0],&gv[0],size);
    for (typename std::map<Key,int>::iterator it = count_.begin(); it != count_.end(); ++it) {
      count_[it->first] = static_cast<int>(gv[cnt]);
      time_[it->first] = gv[cnt+1];
      cnt += 2;
    }
  }
  std::ios_base::fmtflags old(stream.flags());
  stream << std::setprecision(6);
  stream << "  " << name_ << std::endl;
  stream << std::setw(50) << std::right << "Ave. Time (s)"
         << std::setw(25) << std::right << "Ave. #Calls"
         << std::endl;
  for (typename std::map<Key,int>::iterator it = count_.begin(); it != count_.end(); ++it) {
    stream << std::setw(30) << std::right << it->first
           << std::setw(20) << std::right << time_[it->first]/nbatch
           << std::setw(25) << std::right << static_cast<Real>(count_[it->first])/nbatch
           << std::endl;
  }
  stream.flags(old);
}

} // End OED Namespace
} // End ROL Namespace

#endif
