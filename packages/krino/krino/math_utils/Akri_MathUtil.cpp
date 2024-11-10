// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_DiagWriter.hpp>
#include <Akri_MathUtil.hpp>
#include <stk_math/StkVector.hpp>
#include <boost/math/tools/toms748_solve.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <vector>

namespace krino {

double compute_parametric_square_distance(const stk::math::Vector3d childPCoords)
{
  double dist2 = 0.0;
  if (childPCoords[0] < 0.0) dist2 += childPCoords[0]*childPCoords[0];
  if (childPCoords[1] < 0.0) dist2 += childPCoords[1]*childPCoords[1];
  if (childPCoords[2] < 0.0) dist2 += childPCoords[2]*childPCoords[2];
  const double zeta = 1.0 - childPCoords[0] - childPCoords[1] - childPCoords[2];
  if (zeta < 0.0) dist2 += zeta*zeta;
  return dist2;
}

stk::math::Vector3d
get_parametric_coordinates_of_point(const std::vector<stk::math::Vector3d> & nodeCoords, const stk::math::Vector3d & pt)
{
  if (nodeCoords.size() == 3)
  {
    const stk::math::Vector3d relativeCoords1 = nodeCoords[1] - nodeCoords[0];
    const stk::math::Vector3d relativeCoords2 = nodeCoords[2] - nodeCoords[0];
    const stk::math::Vector3d relativeCoords = pt - nodeCoords[0];

    const double a00 = relativeCoords1[0];
    const double a01 = relativeCoords2[0];
    const double a10 = relativeCoords1[1];
    const double a11 = relativeCoords2[1];
    const double b0 = relativeCoords[0];
    const double b1 = relativeCoords[1];
    const double det = a00*a11-a01*a10;
    const double x = (b0*a11-b1*a01)/det;
    const double y = (-b0*a10+b1*a00)/det;
    return stk::math::Vector3d(x,y,0.);
  }
  else
  {
    const stk::math::Vector3d relativeCoords1 = nodeCoords[1] - nodeCoords[0];
    const stk::math::Vector3d relativeCoords2 = nodeCoords[2] - nodeCoords[0];
    const stk::math::Vector3d relativeCoords3 = nodeCoords[3] - nodeCoords[0];
    const stk::math::Vector3d relativeCoords = pt - nodeCoords[0];
    const double a00 = relativeCoords1[0];
    const double a01 = relativeCoords2[0];
    const double a02 = relativeCoords3[0];
    const double a10 = relativeCoords1[1];
    const double a11 = relativeCoords2[1];
    const double a12 = relativeCoords3[1];
    const double a20 = relativeCoords1[2];
    const double a21 = relativeCoords2[2];
    const double a22 = relativeCoords3[2];
    const double b0 = relativeCoords[0];
    const double b1 = relativeCoords[1];
    const double b2 = relativeCoords[2];
    const double det =  a00*(a22*a11-a21*a12)-a10*(a22*a01-a21*a02)+a20*(a12*a01-a11*a02);
    const double x =( b0*(a22*a11-a21*a12)-b1*(a22*a01-a21*a02)+b2*(a12*a01-a11*a02))/det;
    const double y =(-b0*(a22*a10-a20*a12)+b1*(a22*a00-a20*a02)-b2*(a12*a00-a10*a02))/det;
    const double z =( b0*(a21*a10-a20*a11)-b1*(a21*a00-a20*a01)+b2*(a11*a00-a10*a01))/det;
    return stk::math::Vector3d(x,y,z);
  }
}

double compute_linear_root(const std::function<double(const double)> & f,
    const double xa,
    const double xb)
{
  // Does not check assumptions that xb >= xa and that root is bracketed (f(x)=0 for some x, xa <= x <= xb)
  const double xDiff = xb - xa;
  if (xDiff == 0.)
    return xa;

  const double fa = f(xa);
  if (fa == 0.)
    return xa;

  const double fb = f(xb);
  const double fDiff = fb - fa;
  if (fDiff == 0.)
    return xa;

  return xa - xDiff*fa/fDiff;
}

std::pair<bool, double> find_root( const std::function<double(const double)> & f,
    const double xa,
    const double xb,
    const double fa,
    const double fb,
    const unsigned maxIters,
    const double xTol)
{
  // Note:  Seems pretty strange, but the boost::toms748_solve appears to do better on the interval -1->1 than from 0->1
  boost::uintmax_t iterCount = maxIters;
  auto tol_function = [&xTol](const double a, const double b) { return std::abs(b-a) <= xTol; };
  auto result = boost::math::tools::toms748_solve(f, xa, xb, fa, fb, tol_function, iterCount);
  const bool success = iterCount < maxIters;
  return {success, compute_linear_root(f, result.first, result.second)};
}

std::pair<bool, double> find_bracketed_root_newton_raphson( const std::function<std::pair<double,double>(const double)> & f,
    double x,
    double fx,
    double dfx,
    double xa,
    double xb,
    double fa,
    double fb,
    const unsigned maxIters,
    const double fTol)
{
  unsigned iter = 0;
  while (iter++ < maxIters)
  {
    if (dfx == 0.)
    {
      x = 0.5*(xa+xb); // zero slope so use bisection
    }
    else
    {
      x -= fx / dfx;
      if (x < xa || x > xb)
        x = 0.5*(xa+xb); // Newton-Raphson step went out of bounds so use bisection
    }

    std::tie(fx, dfx) = f(x);

    if (std::abs(fx) <= fTol)
      return {true, x};

    if (fx*fb < 0.)
    {
      xa = x;
      fa = fx;
    }
    else
    {
      xb = x;
      fb = fx;
    }

    const double xTol = 4*std::numeric_limits<double>::epsilon()*(std::abs(xa)+std::abs(xb));
    if (xb-xa < xTol)
      return {false, x};
  }
  return {true, x};
}

void attempt_to_bracket_root_newton_raphson( const std::function<std::pair<double,double>(const double)> & f,
    const double guess,
    double & x,
    double & fx,
    double & dfx,
    double & xa,
    double & xb,
    double & fa,
    double & fb,
    bool & solnIsConvergedAtX,
    bool & solnIsBracketed,
    const unsigned maxIters,
    const double fTol)
{
  x = guess;
  std::tie(fx,dfx) = f(x);

  if (std::abs(fx) <= fTol)
  {
    solnIsConvergedAtX = true;
    solnIsBracketed = false;
    return;
  }

  xa = x;
  xb = x;
  fa = fx;
  fb = fx;

  unsigned iter = 0;
  while (iter++ < maxIters)
  {
    if (dfx == 0.)
    {
      // FAILED: zero slope
      solnIsConvergedAtX = false;
      solnIsBracketed = false;
      return;
    }

    xa = x;
    fa = fx;

    x -= fx / dfx;
    std::tie(fx,dfx) = f(x);

    if (std::abs(fx) <= fTol)
    {
      solnIsConvergedAtX = true;
      solnIsBracketed = false;
      return;
    }

    if (fa*fx > 0. && std::abs(fx) > std::abs(fa))
    {
      // Jacobian is poor, try line search
      const double tau = 0.5;
      const int maxLineSearchIters = 5;
      unsigned lineSearchIters = 0;
      while (lineSearchIters++ < maxLineSearchIters && fa*fx > 0. && std::abs(fx) > std::abs(fa))
      {
        x = xa + tau*(x-xa);
        std::tie(fx,dfx) = f(x);
      }
    }

    if (fa*fx < 0.)
    {
      solnIsConvergedAtX = false;
      solnIsBracketed = true;
      xb = x;
      fb = fx;
      if (xa > xb)
      {
        std::swap(xa,xb);
        std::swap(fa,fb);
      }
      return;
    }
  }

  //FAILED: did not bracket root
  solnIsConvergedAtX = false;
  solnIsBracketed = false;
}

std::pair<bool, double> find_root_newton_raphson( const std::function<std::pair<double,double>(const double)> & f,
    const double guess,
    const unsigned maxIters,
    const double fTol)
{
  double x, fx, dfx, xa, xb, fa, fb;
  bool solnIsConvergedAtX, solnIsBracketed;
  attempt_to_bracket_root_newton_raphson(f, guess, x, fx, dfx, xa, xb, fa, fb, solnIsConvergedAtX, solnIsBracketed, maxIters, fTol);

  if (solnIsConvergedAtX)
    return {true, x};

  if (!solnIsBracketed)
    return {false, guess};

  return find_bracketed_root_newton_raphson(f, x, fx, dfx, xa, xb, fa, fb, maxIters, fTol);
}

double find_quadratic_crossing( const double d0, const double d1, const double d2 )
{
  const double epsilon = std::numeric_limits<double>::epsilon()*std::sqrt(d0*d0 + d1*d1 + d2*d2);
  if ( std::fabs(d0) < epsilon ) return 0.0;
  if ( std::fabs(d1) < epsilon ) return 1.0;
  if ( std::fabs(d2) < epsilon ) return 0.5;

  STK_ThrowAssert(d0*d1 < 0.0 && (d0*d2 < 0.0 || d1*d2 < 0.0)); // Insist on one and only one crossing

  const double a = 2.0*(d0 - 2.0*d2 + d1);
  const double b = -3.0*d0 - d1 + 4.0*d2;
  const double c = d0;
  const int sign_b = ( b < 0.0 ) ? -1 : 1;
  const double q = -0.5*(b + sign_b*std::sqrt(b*b-4.0*a*c));

  const int sign_a = ( a < 0.0 ) ? -1 : 1;

  if (q*sign_a > 0.0 && q*sign_a < a*sign_a)
    {
      STK_ThrowAssert(!(c*(( q < 0.0 ) ? -1 : 1) > 0.0 && c*(( q < 0.0 ) ? -1 : 1) < q*(( q < 0.0 ) ? -1 : 1))); // Insist on only one crossing
      return (q/a);
    }
  else
    {
      STK_ThrowAssert(c*(( q < 0.0 ) ? -1 : 1) > 0.0 && c*(( q < 0.0 ) ? -1 : 1) < q*(( q < 0.0 ) ? -1 : 1));
      return (c/q);
    }
}

}


