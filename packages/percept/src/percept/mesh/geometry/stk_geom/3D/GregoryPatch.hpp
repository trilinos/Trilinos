// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_GregoryPatch_hpp
#define percept_GregoryPatch_hpp

#include <percept/PerceptMesh.hpp>

#include <cmath>
#include <iostream>

// #define MyPow(x,n) std::pow(x,n)
// #define MyPow2(x) ((x)*(x))
// #define MyPow3(x) ((x)*(x)*(x))
// #define MyPow4(x) MyPow2(x)*MyPow2(x)
// #define MyInverse(x) (1.0/(x))


namespace percept {

  inline double MyPow(double x, double n) { return std::pow(x,n); }
  inline double MyPow2(double x) {return x*x; }
  inline double MyPow3(double x) {return x*MyPow2(x); }
  inline double MyPow4(double x) {return MyPow2(x)*MyPow2(x); }
  inline double MyInverse(double x) { return 1.0/x; }

  class GregoryPatch {
  public:
    static double
    evalQuad(double u, double v, const MDArray& Cp);

    static double
    evalTri(double u, double v, const MDArray& Cp);

    static void
    evalQuadGrad(double u, double v, const MDArray& Cp, double *grad);

    static void
    evalTriGrad(double u, double v, const MDArray& Cp, double *grad);

    static void
    evalQuadHessian(double u, double v, const MDArray& Cp, double H[2][2]);

    static void
    evalTriHessian(double u, double v, const MDArray& Cp, double H[2][2]);

    static void
    extract_degree_lower_on_boundary(const MDArray& qh, MDArray& q);

    static void
    degree_elevate(MDArray& q, MDArray& qh);

    static void
    fitCubic(MDArray& c, const MDArray& pi, const MDArray& pj, const MDArray& ni, const MDArray& nj);

    static bool
    fitRibbon(MDArray& p, MDArray& q, MDArray& r, MDArray& qh, bool pIsTri, bool rIsTri);

    static void
    fitRibbonNoNeighbor(MDArray& p, MDArray& q, MDArray& qh, bool pIsTri);

  };
}
#endif
