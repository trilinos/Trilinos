// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef CodeGenHelper_hpp
#define CodeGenHelper_hpp

namespace percept {

#if 0
  inline double MyHeaviside(double x, double y)
  {
    if (x < y)
      return 0.0;
    else
      return 1.0;
  }
  inline double MyMax(double x, double y) { return std::max(x,y); }
  inline double MyPow(double x, double n) { return std::pow(x,n); }
  inline double MyPow2(double x) {return x*x; }
  inline double MyPow3(double x) {return x*MyPow2(x); }
  inline double MyPow4(double x) {return MyPow2(x)*MyPow2(x); }
  inline double MyInverse(double x) { return 1.0/x; }
#else
#define MyHeaviside(x, y) ((x) < (y) ? 0.0 : 1.0)
#define MyMax(x, y) ((x) > (y) ? (x) : (y))
#define MyPow(x, n) std::pow((x), (n))
#define MyPow2(x) ((x)*(x))
#define MyPow3(x) ((x)*(x)*(x))
#define MyPow4(x) ((x)*(x)*(x)*(x))
#define MyInverse(x) (1.0/(x))
#endif
}

#endif
