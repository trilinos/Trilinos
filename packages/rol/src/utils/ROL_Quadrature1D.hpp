// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_QUADRATURE1D_HPP
#define ROL_QUADRATURE1D_HPP

#include <vector>

namespace ROL {

template<class Real>
class Quadrature1D {
private:
  std::vector<Real> pts_;
  std::vector<Real> wts_;
  
protected:
  void set(const std::vector<Real> &pts, const std::vector<Real> &wts) {
    pts_.clear(); pts_.assign(pts.begin(),pts.end());
    wts_.clear(); wts_.assign(wts.begin(),wts.end());
  }

public:
  virtual ~Quadrature1D() {}
  Quadrature1D() {}

  void get(std::vector<Real> &pts, std::vector<Real> &wts) const {
    pts.clear(); pts.assign(pts_.begin(),pts_.end());
    wts.clear(); wts.assign(wts_.begin(),wts_.end());
  }

  virtual std::vector<std::vector<Real> > test(const bool printToStream = true,
                                               std::ostream &outStream = std::cout) const = 0;
};

}

#endif
