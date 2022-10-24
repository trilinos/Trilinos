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

#ifndef ROL_OED_RADAMACHER_DEF_HPP
#define ROL_OED_RADAMACHER_DEF_HPP

namespace ROL {
namespace OED {

template<typename Real>
Radamacher<Real>::Radamacher(const Ptr<Vector<Real>> &theta, int size) {
  ProfiledClass<Real,std::string>::rename("OED::Radamacher");
  startTimer("Radamacher");
  Ptr<Vector<Real>> g = theta->dual().clone();
  std::vector<Real> param(1);
  for (int i = 0; i < size; ++i) {
    param[0] = static_cast<Real>(i);
    generate(*g);
    TraceSampler<Real>::setInStorage(*g,param);
  }
  stopTimer("Radamacher");
}

template<typename Real>
void Radamacher<Real>::generate(Vector<Real> &g) const {
  startTimer("generate");
  g.randomize();
  g.applyUnary(Elementwise::Round<Real>());
  g.scale(static_cast<Real>(2));
  g.applyUnary(Elementwise::Shift<Real>(static_cast<Real>(-1)));
  stopTimer("generate");
}

template<typename Real>
void Radamacher<Real>::get(Vector<Real> &F, const std::vector<Real> &param) {
  startTimer("get");
  bool isComputed = TraceSampler<Real>::getFromStorage(F,param);
  if (!isComputed) {
    generate(F);
    TraceSampler<Real>::setInStorage(F,param);
  }
  stopTimer("get");
}

} // End OED Namespace
} // End ROL Namespace

#endif
