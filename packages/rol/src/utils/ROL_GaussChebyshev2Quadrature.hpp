// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_GAUSSCHEBYSHEV2QUADRATURE_HPP
#define ROL_GAUSSCHEBYSHEV2QUADRATURE_HPP

#include "ROL_Quadrature1D.hpp"
#include <cmath>

namespace ROL {

template<class Real>
class GaussChebyshev2Quadrature : public Quadrature1D<Real> {
private:
  const int nQuad_;

public:
  GaussChebyshev2Quadrature(const int nQuad) : nQuad_(nQuad) {
    // Check inputs
    std::vector<Real> pts, wts;
    buildQuadrature(pts,wts);
    Quadrature1D<Real>::set(pts,wts);
  }

  std::vector<std::vector<Real> > test(const bool printToStream = true,
                                       std::ostream &outStream = std::cout) const {
    const int deg = 2*nQuad_-1;
    const Real pi(ROL::ScalarTraits<Real>::pi()), two(2), one(1), four(4), half(0.5), C(4.0/ROL::ScalarTraits<Real>::pi());
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
      if ( i == 0 ) {
        out[i][1] = static_cast<Real>(2);
      }
      else {
        out[i][1] = ((i%2) ? static_cast<Real>(0)
                           : C*two*std::sqrt(pi)*std::tgamma(half*(out[i][0]+one))
                               /(four*std::tgamma(half*out[i][0]+two)));
      }
      for (int j = 0; j < nQuad_; ++j) {
        out[i][2] += wts[j]*std::pow(pts[j],out[i][0]);
      }
      out[i][3] = std::abs(out[i][2] - out[i][1]);
      if (printToStream) {
        outStream << std::fixed << std::setprecision(0) << std::right
                  << std::setw(20) << out[i][0]
                  << std::scientific << std::setprecision(11) << std::right
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
    Real sum(0), pi(ROL::ScalarTraits<Real>::pi()), two(2), one(1), n = static_cast<Real>(nQuad_);
    for (int i = 0; i < nQuad_; ++i) {
      pts[i] = std::cos(static_cast<Real>(i+1)*pi/(n+one));
      wts[i] = pi/(n+one) * std::pow(std::sin(static_cast<Real>(i+1)*pi/(n+one)),two);
      sum += wts[i];
    }
    for (int i = 0; i < nQuad_; ++i) {
      wts[i] *= two/sum;
    }
  }

};

}

#endif
