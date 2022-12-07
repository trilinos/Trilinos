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
