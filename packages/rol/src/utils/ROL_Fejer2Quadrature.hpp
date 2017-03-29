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

#ifndef ROL_FEJER2QUADRATURE_HPP
#define ROL_FEJER2QUADRATURE_HPP

#include "ROL_Quadrature1D.hpp"

namespace ROL {

template<class Real>
class Fejer2Quadrature : public Quadrature1D<Real> {
private:
  const int nQuad_;

public:
  Fejer2Quadrature(const int nQuad) : nQuad_(nQuad) {
    std::vector<Real> pts, wts;
    buildQuadrature(pts,wts);
    Quadrature1D<Real>::set(pts,wts);
  }

  std::vector<std::vector<Real> > test(const bool printToStream = true,
                                       std::ostream &outStream = std::cout) const {
    const int deg = nQuad_ - (nQuad_%1==0 ? 1 : 0);
    std::vector<Real> tmp(4);
    std::vector<std::vector<Real> > out(deg+1,tmp);
    std::vector<Real> pts, wts;
    Quadrature1D<Real>::get(pts,wts);
    for (int i = 0; i < deg+1; ++i) {
      if (printToStream) {
        if (i==0) {
          outStream << std::right
                    << std::setw(20) << "Poly order"
                    << std::setw(20) << "integral"
                    << std::setw(20) << "quad approx"
                    << std::setw(20) << "abs error"
                    << std::endl;
        }
      }
      out[i][0] = static_cast<Real>(i);
      out[i][1] = ((i%2) ? static_cast<Real>(0)
                         : static_cast<Real>(2)/static_cast<Real>(i+1));
      for (int j = 0; j < nQuad_; ++j) {
        out[i][2] += wts[j]*std::pow(pts[j],out[i][0]);
      }
      out[i][3] = std::abs(out[i][2] - out[i][1]);
      if (printToStream) {
        outStream << std::scientific << std::setprecision(11) << std::right
                  << std::setw(20) << out[i][0]
                  << std::setw(20) << out[i][1]
                  << std::setw(20) << out[i][2]
                  << std::setw(20) << out[i][3]
                  << std::endl;
      }
    }
    return out;
  }

private:
  void buildQuadrature(std::vector<Real> &pts, std::vector<Real> &wts) const {
    pts.resize(nQuad_); wts.resize(nQuad_);
    const Real zero(0), one(1), two(2), four(4), pi(Teuchos::ScalarTraits<Real>::pi()), n = static_cast<Real>(nQuad_+1);
    const int n2 = static_cast<int>(0.5*n);
    Real theta(0), sum(0), arg(0), tmp(0);
    for (int i = 1; i < nQuad_+1; ++i) {
      theta = static_cast<Real>(i)*pi/n;
      pts[i-1] = std::cos(theta);
      tmp = zero;
      for (int j = 1; j < n2+1; ++j) {
        arg = two*static_cast<Real>(j)-one;
        tmp += std::sin(arg*theta)/arg;
      }
      wts[i-1] = four/n * std::sin(theta)*tmp;
      sum += wts[i-1];
    }
    if ( nQuad_%2 == 1 ) {
      pts[static_cast<int>((n-two)/two)] = zero;
    }
    for (int i = 0; i < nQuad_; ++i) {
      wts[i] *= two/sum;
    }
  }

};

}

#endif
