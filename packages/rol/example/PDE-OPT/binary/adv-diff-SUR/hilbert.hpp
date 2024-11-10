// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HILBERT_HPP
#define HILBERT_HPP

namespace hilbert {
  void rot(const int n, int &x, int &y, const int rx, const int ry) {
    int t(0);
    if (ry == 0) {
      if (rx == 1) {
        x = n - 1 - x;
        y = n - 1 - y;
      }
      t = x;
      x = y;
      y = t;
    }
  }

  void d2xy(const int m, const int d, int &x, int &y) {
    int n(0), rx(0), ry(0), t(d);
    n = std::pow(2,m);
    x = 0; y = 0;
    for (int s = 1; s < n; s *= 2) {
      rx = (1 & (t/2));
      ry = (1 & (t^rx));
      rot(s,x,y,rx,ry);
      x += s*rx;
      y += s*ry;
      t /= 4;
    }
  }

  void xy2d(const int m, int x, int y, int &d) {
    int n(0), rx(0), ry(0);
    d = 0;
    n = std::pow(2,m);
    for (int s = n/2; s > 0; s /= 2) {
      rx = (x&s) > 0;
      ry = (y&s) > 0;
      d += s*s*((3*rx)^ry);
      rot(s,x,y,rx,ry);
    }
  }
}

#endif
