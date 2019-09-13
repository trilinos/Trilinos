// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "GregoryPatch.hpp"

namespace percept {

double GregoryPatch::
evalQuad(double u, double v, const MDArray& Cp)
{
  /*
     b[0, 0][u, v]
b[1, 0][u, v]
b[2, 0][u, v]
b[3, 0][u, v]
b[0, 1][u, v]
b[1, 1, 1, 0][u, v]
b[1, 1, 0, 1][u, v]
b[2, 1, 1, 0][u, v]
b[2, 1, 0, 1][u, v]
b[3, 1][u, v]
b[0, 2][u, v]
b[1, 2, 1, 0][u, v]
b[1, 2, 0, 1][u, v]
b[2, 2, 1, 0][u, v]
b[2, 2, 0, 1][u, v]
b[3, 2][u, v]
b[0, 3][u, v]
b[1, 3][u, v]
b[2, 3][u, v]
b[3, 3][u, v]

     0:  b[0, 0][u, v]
     1:  b[1, 0][u, v]
     2:  b[2, 0][u, v]
     3:  b[3, 0][u, v]

     4:  b[0, 1][u, v]
     5:  b[1, 1, 1, 0][u, v]
     6:  b[1, 1, 0, 1][u, v]
     7:  b[2, 1, 1, 0][u, v]
     8:  b[2, 1, 0, 1][u, v]
     9:  b[3, 1][u, v]

     10: b[0, 2][u, v]
     11: b[1, 2, 1, 0][u, v]
     12: b[1, 2, 0, 1][u, v]
     13: b[2, 2, 1, 0][u, v]
     14: b[2, 2, 0, 1][u, v]
     15: b[3, 2][u, v]

     16: b[0, 3][u, v]
     17: b[1, 3][u, v]
     18: b[2, 3][u, v]
     19: b[3, 3][u, v]

     edges:
     v = 0
     {0,1,2,3}
     u = 1
     {3,9,15,19}
     v = 1
     {19,18,17,16}
     u = 0
     {16, 10, 4, 0}

     ribbons: (convention is interior is on the left, 1st 4 are "q" boundary points, 2nd 4 are interior "p" points)
     v = 0
     {{0,1,2,3},{4,5,7,9}}
     u = 1
     {{3,9,15,19},{2, 8, 14, 18}}
     v = 1
     {{19,18,17,16},{15,13,11,10}}
     u = 0
     {{16, 10, 4, 0}, {17, 12, 6, 1}}

   */
  double p = 0;
#define NEAR(u,val) (std::fabs((u)-(val)) < 1.e-5)
  if (NEAR(u,0.0))
    {
      
      p = v*(v*(-3*(-1 + v)*Cp(10) + v*Cp(16)) + 3*Cp(4)*MyPow2(-1 + v)) - Cp(0)*MyPow3(-1 + v) ;
    }
  else if (NEAR(u, 1))
    {
      
      p = v*(v*(-3*(-1 + v)*Cp(15) + v*Cp(19)) + 3*Cp(9)*MyPow2(-1 + v)) - Cp(3)*MyPow3(-1 + v) ;
    }
  else if (NEAR(v, 0))
    {
      
      p = 3*u*Cp(1)*MyPow2(-1 + u) - 3*(-1 + u)*Cp(2)*MyPow2(u) - Cp(0)*MyPow3(-1 + u) + Cp(3)*MyPow3(u) ;
    }
  else if (NEAR(v, 1))
    {
      
      p = 3*u*Cp(17)*MyPow2(-1 + u) - 3*(-1 + u)*Cp(18)*MyPow2(u) - Cp(16)*MyPow3(-1 + u) + Cp(19)*MyPow3(u) ;
    }
  else
    {
      
      p = 9*u*v*(u*Cp(5) + v*Cp(6))*MyInverse(u + v)*MyPow2(-1 + u)*MyPow2(-1 + v) - 9*(-1 + u)*v*((-1 + u)*Cp(7) -\
 
  v*Cp(8))*MyInverse(-1 + u - v)*MyPow2(u)*MyPow2(-1 + v) - 9*u*(-1 + v)*(u*Cp(11) - (-1 + v)*Cp(12))*MyInverse(1 + u -\
 
  v)*MyPow2(-1 + u)*MyPow2(v) + 9*(-1 + u)*(-1 + v)*((-1 + u)*Cp(13) + (-1 + v)*Cp(14))*MyInverse(-2 + u +\
 
  v)*MyPow2(u)*MyPow2(v) - 3*v*Cp(4)*MyPow2(-1 + v)*MyPow3(-1 + u) + 3*(-1 + v)*Cp(10)*MyPow2(v)*MyPow3(-1 + u) +\
 
  3*v*Cp(9)*MyPow2(-1 + v)*MyPow3(u) - 3*(-1 + v)*Cp(15)*MyPow2(v)*MyPow3(u) - 3*u*Cp(1)*MyPow2(-1 + u)*MyPow3(-1 + v)\
 
  + 3*(-1 + u)*Cp(2)*MyPow2(u)*MyPow3(-1 + v) + Cp(0)*MyPow3(-1 + u)*MyPow3(-1 + v) - Cp(3)*MyPow3(u)*MyPow3(-1 + v) +\
 
  3*u*Cp(17)*MyPow2(-1 + u)*MyPow3(v) - 3*(-1 + u)*Cp(18)*MyPow2(u)*MyPow3(v) - Cp(16)*MyPow3(-1 + u)*MyPow3(v) +\
 
  Cp(19)*MyPow3(u)*MyPow3(v) ;
    }
  return p;
}

/// @param Cp is for a single coordinate (x,y, or z)
void  GregoryPatch::
evalQuadGrad(double u, double v, const MDArray& Cp, double *grad)
{
  double gu = 0.0, gv=0.0;
  if (NEAR(u,0))
    {
      
      gu = 3*(-3*v*Cp(4)*MyPow2(-1 + v) + 3*v*Cp(6)*MyPow2(-1 + v) + 3*(-1 + v)*Cp(10)*MyPow2(v) - 3*(-1 +\
 
  v)*Cp(12)*MyPow2(v) + Cp(0)*MyPow3(-1 + v) - Cp(1)*MyPow3(-1 + v) - Cp(16)*MyPow3(v) + Cp(17)*MyPow3(v)) ;

      gv = 3*(2*(-1 + v)*v*Cp(4) - 2*(-1 + v)*v*Cp(10) - Cp(0)*MyPow2(-1 + v) + Cp(4)*MyPow2(-1 + v) - Cp(10)*MyPow2(v)\
 
  + Cp(16)*MyPow2(v)) ;
    }
  else if (NEAR(u , 1))
    {
      
      gu = 3*(-3*v*Cp(8)*MyPow2(-1 + v) + 3*v*Cp(9)*MyPow2(-1 + v) + 3*(-1 + v)*Cp(14)*MyPow2(v) - 3*(-1 +\
 
  v)*Cp(15)*MyPow2(v) + Cp(2)*MyPow3(-1 + v) - Cp(3)*MyPow3(-1 + v) - Cp(18)*MyPow3(v) + Cp(19)*MyPow3(v)) ;

      gv = 3*(2*(-1 + v)*v*Cp(9) - 2*(-1 + v)*v*Cp(15) - Cp(3)*MyPow2(-1 + v) + Cp(9)*MyPow2(-1 + v) - Cp(15)*MyPow2(v)\
 
  + Cp(19)*MyPow2(v)) ;
    }
  else if (NEAR(v, 0))
    {
      
      gu = 3*(2*(-1 + u)*u*Cp(1) - 2*(-1 + u)*u*Cp(2) - Cp(0)*MyPow2(-1 + u) + Cp(1)*MyPow2(-1 + u) - Cp(2)*MyPow2(u) +\
 
  Cp(3)*MyPow2(u)) ;

      gv = 3*(-3*u*Cp(1)*MyPow2(-1 + u) + 3*u*Cp(5)*MyPow2(-1 + u) + 3*(-1 + u)*Cp(2)*MyPow2(u) - 3*(-1 +\
 
  u)*Cp(7)*MyPow2(u) + Cp(0)*MyPow3(-1 + u) - Cp(4)*MyPow3(-1 + u) - Cp(3)*MyPow3(u) + Cp(9)*MyPow3(u)) ;
    }
  else if (NEAR(v, 1))
    {
      
      gu = 3*(2*(-1 + u)*u*Cp(17) - 2*(-1 + u)*u*Cp(18) - Cp(16)*MyPow2(-1 + u) + Cp(17)*MyPow2(-1 + u) -\
 
  Cp(18)*MyPow2(u) + Cp(19)*MyPow2(u)) ;

      gv = 3*(-3*u*Cp(11)*MyPow2(-1 + u) + 3*u*Cp(17)*MyPow2(-1 + u) + 3*(-1 + u)*Cp(13)*MyPow2(u) - 3*(-1 +\
 
  u)*Cp(18)*MyPow2(u) + Cp(10)*MyPow3(-1 + u) - Cp(16)*MyPow3(-1 + u) - Cp(15)*MyPow3(u) + Cp(19)*MyPow3(u)) ;
    }
  else
    {
      
      gu = 3*(-6*(-1 + u)*u*v*((-1 + u)*Cp(7) - v*Cp(8))*MyInverse(-1 + u - v)*MyPow2(-1 + v) + 6*(-1 + u)*u*v*(u*Cp(5)\
 
  + v*Cp(6))*MyInverse(u + v)*MyPow2(-1 + v) - 3*v*Cp(4)*MyPow2(-1 + u)*MyPow2(-1 + v) + 3*u*v*Cp(5)*MyInverse(u +\
 
  v)*MyPow2(-1 + u)*MyPow2(-1 + v) + 3*v*(u*Cp(5) + v*Cp(6))*MyInverse(u + v)*MyPow2(-1 + u)*MyPow2(-1 + v) +\
 
  3*v*Cp(9)*MyPow2(u)*MyPow2(-1 + v) - 3*(-1 + u)*v*Cp(7)*MyInverse(-1 + u - v)*MyPow2(u)*MyPow2(-1 + v) - 3*v*((-1 +\
 
  u)*Cp(7) - v*Cp(8))*MyInverse(-1 + u - v)*MyPow2(u)*MyPow2(-1 + v) - 6*(-1 + u)*u*(-1 + v)*(u*Cp(11) - (-1 +\
 
  v)*Cp(12))*MyInverse(1 + u - v)*MyPow2(v) + 6*(-1 + u)*u*(-1 + v)*((-1 + u)*Cp(13) + (-1 + v)*Cp(14))*MyInverse(-2 +\
 
  u + v)*MyPow2(v) + 3*(-1 + v)*Cp(10)*MyPow2(-1 + u)*MyPow2(v) - 3*u*(-1 + v)*Cp(11)*MyInverse(1 + u - v)*MyPow2(-1 +\
 
  u)*MyPow2(v) - 3*(-1 + v)*(u*Cp(11) - (-1 + v)*Cp(12))*MyInverse(1 + u - v)*MyPow2(-1 + u)*MyPow2(v) - 3*(-1 +\
 
  v)*Cp(15)*MyPow2(u)*MyPow2(v) + 3*(-1 + u)*(-1 + v)*Cp(13)*MyInverse(-2 + u + v)*MyPow2(u)*MyPow2(v) + 3*(-1 +\
 
  v)*((-1 + u)*Cp(13) + (-1 + v)*Cp(14))*MyInverse(-2 + u + v)*MyPow2(u)*MyPow2(v) + 3*u*(-1 + v)*(u*Cp(11) - (-1 +\
 
  v)*Cp(12))*MyPow2(-1 + u)*MyPow2(v)*MyPow2(MyInverse(1 + u - v)) + 3*(-1 + u)*v*((-1 + u)*Cp(7) -\
 
  v*Cp(8))*MyPow2(u)*MyPow2(-1 + v)*MyPow2(MyInverse(1 - u + v)) - 3*(-1 + u)*(-1 + v)*((-1 + u)*Cp(13) + (-1 +\
 
  v)*Cp(14))*MyPow2(u)*MyPow2(v)*MyPow2(MyInverse(-2 + u + v)) - 3*u*v*(u*Cp(5) + v*Cp(6))*MyPow2(-1 + u)*MyPow2(-1 +\
 
  v)*MyPow2(MyInverse(u + v)) - 2*(-1 + u)*u*Cp(1)*MyPow3(-1 + v) + 2*(-1 + u)*u*Cp(2)*MyPow3(-1 + v) + Cp(0)*MyPow2(-1\
 
  + u)*MyPow3(-1 + v) - Cp(1)*MyPow2(-1 + u)*MyPow3(-1 + v) + Cp(2)*MyPow2(u)*MyPow3(-1 + v) -\
 
  Cp(3)*MyPow2(u)*MyPow3(-1 + v) + 2*(-1 + u)*u*Cp(17)*MyPow3(v) - 2*(-1 + u)*u*Cp(18)*MyPow3(v) - Cp(16)*MyPow2(-1 +\
 
  u)*MyPow3(v) + Cp(17)*MyPow2(-1 + u)*MyPow3(v) - Cp(18)*MyPow2(u)*MyPow3(v) + Cp(19)*MyPow2(u)*MyPow3(v)) ;

      gv = 3*(-6*u*(-1 + v)*v*(u*Cp(11) - (-1 + v)*Cp(12))*MyInverse(1 + u - v)*MyPow2(-1 + u) + 6*u*(-1 +\
 
  v)*v*(u*Cp(5) + v*Cp(6))*MyInverse(u + v)*MyPow2(-1 + u) - 6*(-1 + u)*(-1 + v)*v*((-1 + u)*Cp(7) -\
 
  v*Cp(8))*MyInverse(-1 + u - v)*MyPow2(u) + 6*(-1 + u)*(-1 + v)*v*((-1 + u)*Cp(13) + (-1 + v)*Cp(14))*MyInverse(-2 + u\
 
  + v)*MyPow2(u) - 3*u*Cp(1)*MyPow2(-1 + u)*MyPow2(-1 + v) + 3*u*v*Cp(6)*MyInverse(u + v)*MyPow2(-1 + u)*MyPow2(-1 + v)\
 
  + 3*u*(u*Cp(5) + v*Cp(6))*MyInverse(u + v)*MyPow2(-1 + u)*MyPow2(-1 + v) + 3*(-1 + u)*Cp(2)*MyPow2(u)*MyPow2(-1 + v)\
 
  + 3*(-1 + u)*v*Cp(8)*MyInverse(-1 + u - v)*MyPow2(u)*MyPow2(-1 + v) - 3*(-1 + u)*((-1 + u)*Cp(7) -\
 
  v*Cp(8))*MyInverse(-1 + u - v)*MyPow2(u)*MyPow2(-1 + v) + 3*u*Cp(17)*MyPow2(-1 + u)*MyPow2(v) + 3*u*(-1 +\
 
  v)*Cp(12)*MyInverse(1 + u - v)*MyPow2(-1 + u)*MyPow2(v) - 3*u*(u*Cp(11) - (-1 + v)*Cp(12))*MyInverse(1 + u -\
 
  v)*MyPow2(-1 + u)*MyPow2(v) - 3*(-1 + u)*Cp(18)*MyPow2(u)*MyPow2(v) + 3*(-1 + u)*(-1 + v)*Cp(14)*MyInverse(-2 + u +\
 
  v)*MyPow2(u)*MyPow2(v) + 3*(-1 + u)*((-1 + u)*Cp(13) + (-1 + v)*Cp(14))*MyInverse(-2 + u + v)*MyPow2(u)*MyPow2(v) -\
 
  3*u*(-1 + v)*(u*Cp(11) - (-1 + v)*Cp(12))*MyPow2(-1 + u)*MyPow2(v)*MyPow2(MyInverse(1 + u - v)) - 3*(-1 + u)*v*((-1 +\
 
  u)*Cp(7) - v*Cp(8))*MyPow2(u)*MyPow2(-1 + v)*MyPow2(MyInverse(1 - u + v)) - 3*(-1 + u)*(-1 + v)*((-1 + u)*Cp(13) +\
 
  (-1 + v)*Cp(14))*MyPow2(u)*MyPow2(v)*MyPow2(MyInverse(-2 + u + v)) - 3*u*v*(u*Cp(5) + v*Cp(6))*MyPow2(-1 +\
 
  u)*MyPow2(-1 + v)*MyPow2(MyInverse(u + v)) - 2*(-1 + v)*v*Cp(4)*MyPow3(-1 + u) + 2*(-1 + v)*v*Cp(10)*MyPow3(-1 + u) +\
 
  Cp(0)*MyPow2(-1 + v)*MyPow3(-1 + u) - Cp(4)*MyPow2(-1 + v)*MyPow3(-1 + u) + Cp(10)*MyPow2(v)*MyPow3(-1 + u) -\
 
  Cp(16)*MyPow2(v)*MyPow3(-1 + u) + 2*(-1 + v)*v*Cp(9)*MyPow3(u) - 2*(-1 + v)*v*Cp(15)*MyPow3(u) - Cp(3)*MyPow2(-1 +\
 
  v)*MyPow3(u) + Cp(9)*MyPow2(-1 + v)*MyPow3(u) - Cp(15)*MyPow2(v)*MyPow3(u) + Cp(19)*MyPow2(v)*MyPow3(u)) ;
    }
  grad[0] = gu;
  grad[1] = gv;
}

/// @param Cp is for a single coordinate (x,y, or z)
void  GregoryPatch::
evalQuadHessian(double u, double v, const MDArray& Cp, double H[2][2])
{
  double H00, H01, H10, H11;
  if (NEAR(u,0))
    {
      
      H00 = 3*(6*v*Cp(4)*MyPow2(-1 + v) + 6*Cp(5)*MyPow2(-1 + v) - 6*Cp(6)*MyPow2(-1 + v) - 12*v*Cp(6)*MyPow2(-1 + v) +\
 
  6*v*(-Cp(7) - v*Cp(8))*MyInverse(-1 - v)*MyPow2(-1 + v) - 6*(-1 + v)*Cp(10)*MyPow2(v) + 6*Cp(11)*MyPow2(v) -\
 
  6*Cp(12)*MyPow2(v) + 12*(-1 + v)*Cp(12)*MyPow2(v) - 6*(-1 + v)*(-Cp(13) + (-1 + v)*Cp(14))*MyInverse(-2 +\
 
  v)*MyPow2(v) - 2*Cp(0)*MyPow3(-1 + v) + 4*Cp(1)*MyPow3(-1 + v) - 2*Cp(2)*MyPow3(-1 + v) + 2*Cp(16)*MyPow3(v) -\
 
  4*Cp(17)*MyPow3(v) + 2*Cp(18)*MyPow3(v)) ;

      H01 = 3*(-6*(-1 + v)*v*Cp(4) + 6*(-1 + v)*v*Cp(6) + 6*(-1 + v)*v*Cp(10) - 6*(-1 + v)*v*Cp(12) + 3*Cp(0)*MyPow2(-1\
 
  + v) - 3*Cp(1)*MyPow2(-1 + v) - 3*Cp(4)*MyPow2(-1 + v) + 3*Cp(6)*MyPow2(-1 + v) + 3*Cp(10)*MyPow2(v) -\
 
  3*Cp(12)*MyPow2(v) - 3*Cp(16)*MyPow2(v) + 3*Cp(17)*MyPow2(v)) ;

      H10 = 3*(-6*(-1 + v)*v*Cp(4) + 6*(-1 + v)*v*Cp(6) + 6*(-1 + v)*v*Cp(10) - 6*(-1 + v)*v*Cp(12) + 3*Cp(0)*MyPow2(-1\
 
  + v) - 3*Cp(1)*MyPow2(-1 + v) - 3*Cp(4)*MyPow2(-1 + v) + 3*Cp(6)*MyPow2(-1 + v) + 3*Cp(10)*MyPow2(v) -\
 
  3*Cp(12)*MyPow2(v) - 3*Cp(16)*MyPow2(v) + 3*Cp(17)*MyPow2(v)) ;

      H11 = 3*(-2*(-1 + v)*Cp(0) + 4*(-1 + v)*Cp(4) + 2*v*Cp(4) - 2*(-1 + v)*Cp(10) - 4*v*Cp(10) + 2*v*Cp(16)) ;
    }
  else if (NEAR(u , 1))
    {
      
      H00 = 3*(6*Cp(7)*MyPow2(-1 + v) - 6*Cp(8)*MyPow2(-1 + v) - 12*v*Cp(8)*MyPow2(-1 + v) + 6*v*Cp(9)*MyPow2(-1 + v) +\
 
  6*v*(Cp(5) + v*Cp(6))*MyInverse(1 + v)*MyPow2(-1 + v) + 6*Cp(13)*MyPow2(v) - 6*Cp(14)*MyPow2(v) + 12*(-1 +\
 
  v)*Cp(14)*MyPow2(v) - 6*(-1 + v)*Cp(15)*MyPow2(v) - 6*(-1 + v)*(Cp(11) - (-1 + v)*Cp(12))*MyInverse(2 - v)*MyPow2(v)\
 
  - 2*Cp(1)*MyPow3(-1 + v) + 4*Cp(2)*MyPow3(-1 + v) - 2*Cp(3)*MyPow3(-1 + v) + 2*Cp(17)*MyPow3(v) - 4*Cp(18)*MyPow3(v)\
 
  + 2*Cp(19)*MyPow3(v)) ;

      H01 = 3*(-6*(-1 + v)*v*Cp(8) + 6*(-1 + v)*v*Cp(9) + 6*(-1 + v)*v*Cp(14) - 6*(-1 + v)*v*Cp(15) + 3*Cp(2)*MyPow2(-1\
 
  + v) - 3*Cp(3)*MyPow2(-1 + v) - 3*Cp(8)*MyPow2(-1 + v) + 3*Cp(9)*MyPow2(-1 + v) + 3*Cp(14)*MyPow2(v) -\
 
  3*Cp(15)*MyPow2(v) - 3*Cp(18)*MyPow2(v) + 3*Cp(19)*MyPow2(v)) ;

      H10 = 3*(-6*(-1 + v)*v*Cp(8) + 6*(-1 + v)*v*Cp(9) + 6*(-1 + v)*v*Cp(14) - 6*(-1 + v)*v*Cp(15) + 3*Cp(2)*MyPow2(-1\
 
  + v) - 3*Cp(3)*MyPow2(-1 + v) - 3*Cp(8)*MyPow2(-1 + v) + 3*Cp(9)*MyPow2(-1 + v) + 3*Cp(14)*MyPow2(v) -\
 
  3*Cp(15)*MyPow2(v) - 3*Cp(18)*MyPow2(v) + 3*Cp(19)*MyPow2(v)) ;

      H11 = 3*(-2*(-1 + v)*Cp(3) + 4*(-1 + v)*Cp(9) + 2*v*Cp(9) - 2*(-1 + v)*Cp(15) - 4*v*Cp(15) + 2*v*Cp(19)) ;
    }
  else if (NEAR(v, 0))
    {
      
      H00 = 3*(-2*(-1 + u)*Cp(0) + 4*(-1 + u)*Cp(1) + 2*u*Cp(1) - 2*(-1 + u)*Cp(2) - 4*u*Cp(2) + 2*u*Cp(3)) ;

      H01 = 3*(-6*(-1 + u)*u*Cp(1) + 6*(-1 + u)*u*Cp(2) + 6*(-1 + u)*u*Cp(5) - 6*(-1 + u)*u*Cp(7) + 3*Cp(0)*MyPow2(-1 +\
 
  u) - 3*Cp(1)*MyPow2(-1 + u) - 3*Cp(4)*MyPow2(-1 + u) + 3*Cp(5)*MyPow2(-1 + u) + 3*Cp(2)*MyPow2(u) - 3*Cp(3)*MyPow2(u)\
 
  - 3*Cp(7)*MyPow2(u) + 3*Cp(9)*MyPow2(u)) ;

      H10 = 3*(-6*(-1 + u)*u*Cp(1) + 6*(-1 + u)*u*Cp(2) + 6*(-1 + u)*u*Cp(5) - 6*(-1 + u)*u*Cp(7) + 3*Cp(0)*MyPow2(-1 +\
 
  u) - 3*Cp(1)*MyPow2(-1 + u) - 3*Cp(4)*MyPow2(-1 + u) + 3*Cp(5)*MyPow2(-1 + u) + 3*Cp(2)*MyPow2(u) - 3*Cp(3)*MyPow2(u)\
 
  - 3*Cp(7)*MyPow2(u) + 3*Cp(9)*MyPow2(u)) ;

      H11 = 3*(6*u*Cp(1)*MyPow2(-1 + u) - 6*Cp(5)*MyPow2(-1 + u) - 12*u*Cp(5)*MyPow2(-1 + u) + 6*Cp(6)*MyPow2(-1 + u) +\
 
  6*u*(u*Cp(11) + Cp(12))*MyInverse(1 + u)*MyPow2(-1 + u) - 6*(-1 + u)*Cp(2)*MyPow2(u) - 6*Cp(7)*MyPow2(u) + 12*(-1 +\
 
  u)*Cp(7)*MyPow2(u) + 6*Cp(8)*MyPow2(u) - 6*(-1 + u)*((-1 + u)*Cp(13) - Cp(14))*MyInverse(-2 + u)*MyPow2(u) -\
 
  2*Cp(0)*MyPow3(-1 + u) + 4*Cp(4)*MyPow3(-1 + u) - 2*Cp(10)*MyPow3(-1 + u) + 2*Cp(3)*MyPow3(u) - 4*Cp(9)*MyPow3(u) +\
 
  2*Cp(15)*MyPow3(u)) ;
    }
  else if (NEAR(v, 1))
    {
      
      H00 = 3*(-2*(-1 + u)*Cp(16) + 4*(-1 + u)*Cp(17) + 2*u*Cp(17) - 2*(-1 + u)*Cp(18) - 4*u*Cp(18) + 2*u*Cp(19)) ;

      H01 = 3*(-6*(-1 + u)*u*Cp(11) + 6*(-1 + u)*u*Cp(13) + 6*(-1 + u)*u*Cp(17) - 6*(-1 + u)*u*Cp(18) +\
 
  3*Cp(10)*MyPow2(-1 + u) - 3*Cp(11)*MyPow2(-1 + u) - 3*Cp(16)*MyPow2(-1 + u) + 3*Cp(17)*MyPow2(-1 + u) +\
 
  3*Cp(13)*MyPow2(u) - 3*Cp(15)*MyPow2(u) - 3*Cp(18)*MyPow2(u) + 3*Cp(19)*MyPow2(u)) ;

      H10 = 3*(-6*(-1 + u)*u*Cp(11) + 6*(-1 + u)*u*Cp(13) + 6*(-1 + u)*u*Cp(17) - 6*(-1 + u)*u*Cp(18) +\
 
  3*Cp(10)*MyPow2(-1 + u) - 3*Cp(11)*MyPow2(-1 + u) - 3*Cp(16)*MyPow2(-1 + u) + 3*Cp(17)*MyPow2(-1 + u) +\
 
  3*Cp(13)*MyPow2(u) - 3*Cp(15)*MyPow2(u) - 3*Cp(18)*MyPow2(u) + 3*Cp(19)*MyPow2(u)) ;

      H11 = 3*(-6*Cp(11)*MyPow2(-1 + u) - 12*u*Cp(11)*MyPow2(-1 + u) + 6*Cp(12)*MyPow2(-1 + u) + 6*u*Cp(17)*MyPow2(-1 +\
 
  u) + 6*u*(u*Cp(5) + Cp(6))*MyInverse(1 + u)*MyPow2(-1 + u) - 6*Cp(13)*MyPow2(u) + 12*(-1 + u)*Cp(13)*MyPow2(u) +\
 
  6*Cp(14)*MyPow2(u) - 6*(-1 + u)*Cp(18)*MyPow2(u) - 6*(-1 + u)*((-1 + u)*Cp(7) - Cp(8))*MyInverse(-2 + u)*MyPow2(u) -\
 
  2*Cp(4)*MyPow3(-1 + u) + 4*Cp(10)*MyPow3(-1 + u) - 2*Cp(16)*MyPow3(-1 + u) + 2*Cp(9)*MyPow3(u) - 4*Cp(15)*MyPow3(u) +\
 
  2*Cp(19)*MyPow3(u)) ;
    }
  else
    {
      
      H00 = 3*(-6*(-1 + u)*v*Cp(4)*MyPow2(-1 + v) + 6*u*v*Cp(9)*MyPow2(-1 + v) - 12*(-1 + u)*u*v*Cp(7)*MyInverse(-1 + u\
 
  - v)*MyPow2(-1 + v) - 6*(-1 + u)*v*((-1 + u)*Cp(7) - v*Cp(8))*MyInverse(-1 + u - v)*MyPow2(-1 + v) - 12*u*v*((-1 +\
 
  u)*Cp(7) - v*Cp(8))*MyInverse(-1 + u - v)*MyPow2(-1 + v) + 12*(-1 + u)*u*v*Cp(5)*MyInverse(u + v)*MyPow2(-1 + v) +\
 
  12*(-1 + u)*v*(u*Cp(5) + v*Cp(6))*MyInverse(u + v)*MyPow2(-1 + v) + 6*u*v*(u*Cp(5) + v*Cp(6))*MyInverse(u +\
 
  v)*MyPow2(-1 + v) + 6*v*Cp(5)*MyInverse(u + v)*MyPow2(-1 + u)*MyPow2(-1 + v) - 6*v*Cp(7)*MyInverse(-1 + u -\
 
  v)*MyPow2(u)*MyPow2(-1 + v) + 6*(-1 + u)*(-1 + v)*Cp(10)*MyPow2(v) - 6*u*(-1 + v)*Cp(15)*MyPow2(v) - 12*(-1 +\
 
  u)*u*(-1 + v)*Cp(11)*MyInverse(1 + u - v)*MyPow2(v) - 12*(-1 + u)*(-1 + v)*(u*Cp(11) - (-1 + v)*Cp(12))*MyInverse(1 +\
 
  u - v)*MyPow2(v) - 6*u*(-1 + v)*(u*Cp(11) - (-1 + v)*Cp(12))*MyInverse(1 + u - v)*MyPow2(v) + 12*(-1 + u)*u*(-1 +\
 
  v)*Cp(13)*MyInverse(-2 + u + v)*MyPow2(v) + 6*(-1 + u)*(-1 + v)*((-1 + u)*Cp(13) + (-1 + v)*Cp(14))*MyInverse(-2 + u\
 
  + v)*MyPow2(v) + 12*u*(-1 + v)*((-1 + u)*Cp(13) + (-1 + v)*Cp(14))*MyInverse(-2 + u + v)*MyPow2(v) - 6*(-1 +\
 
  v)*Cp(11)*MyInverse(1 + u - v)*MyPow2(-1 + u)*MyPow2(v) + 6*(-1 + v)*Cp(13)*MyInverse(-2 + u + v)*MyPow2(u)*MyPow2(v)\
 
  + 6*(-1 + u)*u*v*((-1 + u)*Cp(7) - v*Cp(8))*MyPow2(-1 + v)*MyPow2(MyInverse(-1 + u - v)) + 3*(-1 +\
 
  u)*v*Cp(7)*MyPow2(u)*MyPow2(-1 + v)*MyPow2(MyInverse(-1 + u - v)) + 3*v*((-1 + u)*Cp(7) -\
 
  v*Cp(8))*MyPow2(u)*MyPow2(-1 + v)*MyPow2(MyInverse(-1 + u - v)) + 12*(-1 + u)*u*(-1 + v)*(u*Cp(11) - (-1 +\
 
  v)*Cp(12))*MyPow2(v)*MyPow2(MyInverse(1 + u - v)) + 6*u*(-1 + v)*Cp(11)*MyPow2(-1 + u)*MyPow2(v)*MyPow2(MyInverse(1 +\
 
  u - v)) + 6*(-1 + v)*(u*Cp(11) - (-1 + v)*Cp(12))*MyPow2(-1 + u)*MyPow2(v)*MyPow2(MyInverse(1 + u - v)) + 6*(-1 +\
 
  u)*u*v*((-1 + u)*Cp(7) - v*Cp(8))*MyPow2(-1 + v)*MyPow2(MyInverse(1 - u + v)) + 3*(-1 +\
 
  u)*v*Cp(7)*MyPow2(u)*MyPow2(-1 + v)*MyPow2(MyInverse(1 - u + v)) + 3*v*((-1 + u)*Cp(7) - v*Cp(8))*MyPow2(u)*MyPow2(-1\
 
  + v)*MyPow2(MyInverse(1 - u + v)) - 12*(-1 + u)*u*(-1 + v)*((-1 + u)*Cp(13) + (-1 +\
 
  v)*Cp(14))*MyPow2(v)*MyPow2(MyInverse(-2 + u + v)) - 6*(-1 + u)*(-1 +\
 
  v)*Cp(13)*MyPow2(u)*MyPow2(v)*MyPow2(MyInverse(-2 + u + v)) - 6*(-1 + v)*((-1 + u)*Cp(13) + (-1 +\
 
  v)*Cp(14))*MyPow2(u)*MyPow2(v)*MyPow2(MyInverse(-2 + u + v)) - 12*(-1 + u)*u*v*(u*Cp(5) + v*Cp(6))*MyPow2(-1 +\
 
  v)*MyPow2(MyInverse(u + v)) - 6*u*v*Cp(5)*MyPow2(-1 + u)*MyPow2(-1 + v)*MyPow2(MyInverse(u + v)) - 6*v*(u*Cp(5) +\
 
  v*Cp(6))*MyPow2(-1 + u)*MyPow2(-1 + v)*MyPow2(MyInverse(u + v)) + 2*(-1 + u)*Cp(0)*MyPow3(-1 + v) - 4*(-1 +\
 
  u)*Cp(1)*MyPow3(-1 + v) - 2*u*Cp(1)*MyPow3(-1 + v) + 2*(-1 + u)*Cp(2)*MyPow3(-1 + v) + 4*u*Cp(2)*MyPow3(-1 + v) -\
 
  2*u*Cp(3)*MyPow3(-1 + v) - 2*(-1 + u)*Cp(16)*MyPow3(v) + 4*(-1 + u)*Cp(17)*MyPow3(v) + 2*u*Cp(17)*MyPow3(v) - 2*(-1 +\
 
  u)*Cp(18)*MyPow3(v) - 4*u*Cp(18)*MyPow3(v) + 2*u*Cp(19)*MyPow3(v) - 6*u*(-1 + v)*(u*Cp(11) - (-1 +\
 
  v)*Cp(12))*MyPow2(-1 + u)*MyPow2(v)*MyPow3(MyInverse(1 + u - v)) + 6*(-1 + u)*v*((-1 + u)*Cp(7) -\
 
  v*Cp(8))*MyPow2(u)*MyPow2(-1 + v)*MyPow3(MyInverse(1 - u + v)) + 6*(-1 + u)*(-1 + v)*((-1 + u)*Cp(13) + (-1 +\
 
  v)*Cp(14))*MyPow2(u)*MyPow2(v)*MyPow3(MyInverse(-2 + u + v)) + 6*u*v*(u*Cp(5) + v*Cp(6))*MyPow2(-1 + u)*MyPow2(-1 +\
 
  v)*MyPow3(MyInverse(u + v))) ;

      H01 = 3*(-12*(-1 + u)*u*(-1 + v)*v*((-1 + u)*Cp(7) - v*Cp(8))*MyInverse(-1 + u - v) - 12*(-1 + u)*u*(-1 +\
 
  v)*v*(u*Cp(11) - (-1 + v)*Cp(12))*MyInverse(1 + u - v) + 12*(-1 + u)*u*(-1 + v)*v*((-1 + u)*Cp(13) + (-1 +\
 
  v)*Cp(14))*MyInverse(-2 + u + v) + 12*(-1 + u)*u*(-1 + v)*v*(u*Cp(5) + v*Cp(6))*MyInverse(u + v) - 6*(-1 +\
 
  v)*v*Cp(4)*MyPow2(-1 + u) + 6*(-1 + v)*v*Cp(10)*MyPow2(-1 + u) - 6*u*(-1 + v)*v*Cp(11)*MyInverse(1 + u - v)*MyPow2(-1\
 
  + u) - 6*(-1 + v)*v*(u*Cp(11) - (-1 + v)*Cp(12))*MyInverse(1 + u - v)*MyPow2(-1 + u) + 6*u*(-1 +\
 
  v)*v*Cp(5)*MyInverse(u + v)*MyPow2(-1 + u) + 6*(-1 + v)*v*(u*Cp(5) + v*Cp(6))*MyInverse(u + v)*MyPow2(-1 + u) + 6*(-1\
 
  + v)*v*Cp(9)*MyPow2(u) - 6*(-1 + v)*v*Cp(15)*MyPow2(u) - 6*(-1 + u)*(-1 + v)*v*Cp(7)*MyInverse(-1 + u - v)*MyPow2(u)\
 
  - 6*(-1 + v)*v*((-1 + u)*Cp(7) - v*Cp(8))*MyInverse(-1 + u - v)*MyPow2(u) + 6*(-1 + u)*(-1 + v)*v*Cp(13)*MyInverse(-2\
 
  + u + v)*MyPow2(u) + 6*(-1 + v)*v*((-1 + u)*Cp(13) + (-1 + v)*Cp(14))*MyInverse(-2 + u + v)*MyPow2(u) - 6*(-1 +\
 
  u)*u*Cp(1)*MyPow2(-1 + v) + 6*(-1 + u)*u*Cp(2)*MyPow2(-1 + v) + 6*(-1 + u)*u*v*Cp(8)*MyInverse(-1 + u - v)*MyPow2(-1\
 
  + v) - 6*(-1 + u)*u*((-1 + u)*Cp(7) - v*Cp(8))*MyInverse(-1 + u - v)*MyPow2(-1 + v) + 6*(-1 +\
 
  u)*u*v*Cp(6)*MyInverse(u + v)*MyPow2(-1 + v) + 6*(-1 + u)*u*(u*Cp(5) + v*Cp(6))*MyInverse(u + v)*MyPow2(-1 + v) +\
 
  3*Cp(0)*MyPow2(-1 + u)*MyPow2(-1 + v) - 3*Cp(1)*MyPow2(-1 + u)*MyPow2(-1 + v) - 3*Cp(4)*MyPow2(-1 + u)*MyPow2(-1 + v)\
 
  + 3*u*Cp(5)*MyInverse(u + v)*MyPow2(-1 + u)*MyPow2(-1 + v) + 3*v*Cp(6)*MyInverse(u + v)*MyPow2(-1 + u)*MyPow2(-1 + v)\
 
  + 3*(u*Cp(5) + v*Cp(6))*MyInverse(u + v)*MyPow2(-1 + u)*MyPow2(-1 + v) + 3*Cp(2)*MyPow2(u)*MyPow2(-1 + v) -\
 
  3*Cp(3)*MyPow2(u)*MyPow2(-1 + v) + 3*Cp(9)*MyPow2(u)*MyPow2(-1 + v) - 3*(-1 + u)*Cp(7)*MyInverse(-1 + u -\
 
  v)*MyPow2(u)*MyPow2(-1 + v) + 3*v*Cp(8)*MyInverse(-1 + u - v)*MyPow2(u)*MyPow2(-1 + v) - 3*((-1 + u)*Cp(7) -\
 
  v*Cp(8))*MyInverse(-1 + u - v)*MyPow2(u)*MyPow2(-1 + v) + 6*(-1 + u)*u*Cp(17)*MyPow2(v) - 6*(-1 +\
 
  u)*u*Cp(18)*MyPow2(v) + 6*(-1 + u)*u*(-1 + v)*Cp(12)*MyInverse(1 + u - v)*MyPow2(v) - 6*(-1 + u)*u*(u*Cp(11) - (-1 +\
 
  v)*Cp(12))*MyInverse(1 + u - v)*MyPow2(v) + 6*(-1 + u)*u*(-1 + v)*Cp(14)*MyInverse(-2 + u + v)*MyPow2(v) + 6*(-1 +\
 
  u)*u*((-1 + u)*Cp(13) + (-1 + v)*Cp(14))*MyInverse(-2 + u + v)*MyPow2(v) + 3*Cp(10)*MyPow2(-1 + u)*MyPow2(v) -\
 
  3*Cp(16)*MyPow2(-1 + u)*MyPow2(v) + 3*Cp(17)*MyPow2(-1 + u)*MyPow2(v) - 3*u*Cp(11)*MyInverse(1 + u - v)*MyPow2(-1 +\
 
  u)*MyPow2(v) + 3*(-1 + v)*Cp(12)*MyInverse(1 + u - v)*MyPow2(-1 + u)*MyPow2(v) - 3*(u*Cp(11) - (-1 +\
 
  v)*Cp(12))*MyInverse(1 + u - v)*MyPow2(-1 + u)*MyPow2(v) - 3*Cp(15)*MyPow2(u)*MyPow2(v) -\
 
  3*Cp(18)*MyPow2(u)*MyPow2(v) + 3*Cp(19)*MyPow2(u)*MyPow2(v) + 3*(-1 + u)*Cp(13)*MyInverse(-2 + u +\
 
  v)*MyPow2(u)*MyPow2(v) + 3*(-1 + v)*Cp(14)*MyInverse(-2 + u + v)*MyPow2(u)*MyPow2(v) + 3*((-1 + u)*Cp(13) + (-1 +\
 
  v)*Cp(14))*MyInverse(-2 + u + v)*MyPow2(u)*MyPow2(v) - 6*(-1 + u)*u*v*((-1 + u)*Cp(7) - v*Cp(8))*MyPow2(-1 +\
 
  v)*MyPow2(MyInverse(-1 + u - v)) - 3*(-1 + u)*v*Cp(7)*MyPow2(u)*MyPow2(-1 + v)*MyPow2(MyInverse(-1 + u - v)) -\
 
  3*v*((-1 + u)*Cp(7) - v*Cp(8))*MyPow2(u)*MyPow2(-1 + v)*MyPow2(MyInverse(-1 + u - v)) + 6*u*(-1 + v)*v*(u*Cp(11) -\
 
  (-1 + v)*Cp(12))*MyPow2(-1 + u)*MyPow2(MyInverse(1 + u - v)) - 6*(-1 + u)*u*(-1 + v)*(u*Cp(11) - (-1 +\
 
  v)*Cp(12))*MyPow2(v)*MyPow2(MyInverse(1 + u - v)) - 3*u*(-1 + v)*Cp(11)*MyPow2(-1 + u)*MyPow2(v)*MyPow2(MyInverse(1 +\
 
  u - v)) - 3*u*(-1 + v)*Cp(12)*MyPow2(-1 + u)*MyPow2(v)*MyPow2(MyInverse(1 + u - v)) + 3*u*(u*Cp(11) - (-1 +\
 
  v)*Cp(12))*MyPow2(-1 + u)*MyPow2(v)*MyPow2(MyInverse(1 + u - v)) - 3*(-1 + v)*(u*Cp(11) - (-1 + v)*Cp(12))*MyPow2(-1\
 
  + u)*MyPow2(v)*MyPow2(MyInverse(1 + u - v)) + 6*(-1 + u)*(-1 + v)*v*((-1 + u)*Cp(7) -\
 
  v*Cp(8))*MyPow2(u)*MyPow2(MyInverse(1 - u + v)) - 3*(-1 + u)*v*Cp(8)*MyPow2(u)*MyPow2(-1 + v)*MyPow2(MyInverse(1 - u\
 
  + v)) + 3*(-1 + u)*((-1 + u)*Cp(7) - v*Cp(8))*MyPow2(u)*MyPow2(-1 + v)*MyPow2(MyInverse(1 - u + v)) - 6*(-1 + u)*(-1\
 
  + v)*v*((-1 + u)*Cp(13) + (-1 + v)*Cp(14))*MyPow2(u)*MyPow2(MyInverse(-2 + u + v)) - 6*(-1 + u)*u*(-1 + v)*((-1 +\
 
  u)*Cp(13) + (-1 + v)*Cp(14))*MyPow2(v)*MyPow2(MyInverse(-2 + u + v)) - 3*(-1 + u)*(-1 +\
 
  v)*Cp(13)*MyPow2(u)*MyPow2(v)*MyPow2(MyInverse(-2 + u + v)) - 3*(-1 + u)*(-1 +\
 
  v)*Cp(14)*MyPow2(u)*MyPow2(v)*MyPow2(MyInverse(-2 + u + v)) - 3*(-1 + u)*((-1 + u)*Cp(13) + (-1 +\
 
  v)*Cp(14))*MyPow2(u)*MyPow2(v)*MyPow2(MyInverse(-2 + u + v)) - 3*(-1 + v)*((-1 + u)*Cp(13) + (-1 +\
 
  v)*Cp(14))*MyPow2(u)*MyPow2(v)*MyPow2(MyInverse(-2 + u + v)) - 6*u*(-1 + v)*v*(u*Cp(5) + v*Cp(6))*MyPow2(-1 +\
 
  u)*MyPow2(MyInverse(u + v)) - 6*(-1 + u)*u*v*(u*Cp(5) + v*Cp(6))*MyPow2(-1 + v)*MyPow2(MyInverse(u + v)) -\
 
  3*u*v*Cp(5)*MyPow2(-1 + u)*MyPow2(-1 + v)*MyPow2(MyInverse(u + v)) - 3*u*v*Cp(6)*MyPow2(-1 + u)*MyPow2(-1 +\
 
  v)*MyPow2(MyInverse(u + v)) - 3*u*(u*Cp(5) + v*Cp(6))*MyPow2(-1 + u)*MyPow2(-1 + v)*MyPow2(MyInverse(u + v)) -\
 
  3*v*(u*Cp(5) + v*Cp(6))*MyPow2(-1 + u)*MyPow2(-1 + v)*MyPow2(MyInverse(u + v)) + 6*u*(-1 + v)*(u*Cp(11) - (-1 +\
 
  v)*Cp(12))*MyPow2(-1 + u)*MyPow2(v)*MyPow3(MyInverse(1 + u - v)) - 6*(-1 + u)*v*((-1 + u)*Cp(7) -\
 
  v*Cp(8))*MyPow2(u)*MyPow2(-1 + v)*MyPow3(MyInverse(1 - u + v)) + 6*(-1 + u)*(-1 + v)*((-1 + u)*Cp(13) + (-1 +\
 
  v)*Cp(14))*MyPow2(u)*MyPow2(v)*MyPow3(MyInverse(-2 + u + v)) + 6*u*v*(u*Cp(5) + v*Cp(6))*MyPow2(-1 + u)*MyPow2(-1 +\
 
  v)*MyPow3(MyInverse(u + v))) ;

      H10 = 3*(-12*(-1 + u)*u*(-1 + v)*v*((-1 + u)*Cp(7) - v*Cp(8))*MyInverse(-1 + u - v) - 12*(-1 + u)*u*(-1 +\
 
  v)*v*(u*Cp(11) - (-1 + v)*Cp(12))*MyInverse(1 + u - v) + 12*(-1 + u)*u*(-1 + v)*v*((-1 + u)*Cp(13) + (-1 +\
 
  v)*Cp(14))*MyInverse(-2 + u + v) + 12*(-1 + u)*u*(-1 + v)*v*(u*Cp(5) + v*Cp(6))*MyInverse(u + v) - 6*(-1 +\
 
  v)*v*Cp(4)*MyPow2(-1 + u) + 6*(-1 + v)*v*Cp(10)*MyPow2(-1 + u) - 6*u*(-1 + v)*v*Cp(11)*MyInverse(1 + u - v)*MyPow2(-1\
 
  + u) - 6*(-1 + v)*v*(u*Cp(11) - (-1 + v)*Cp(12))*MyInverse(1 + u - v)*MyPow2(-1 + u) + 6*u*(-1 +\
 
  v)*v*Cp(5)*MyInverse(u + v)*MyPow2(-1 + u) + 6*(-1 + v)*v*(u*Cp(5) + v*Cp(6))*MyInverse(u + v)*MyPow2(-1 + u) + 6*(-1\
 
  + v)*v*Cp(9)*MyPow2(u) - 6*(-1 + v)*v*Cp(15)*MyPow2(u) - 6*(-1 + u)*(-1 + v)*v*Cp(7)*MyInverse(-1 + u - v)*MyPow2(u)\
 
  - 6*(-1 + v)*v*((-1 + u)*Cp(7) - v*Cp(8))*MyInverse(-1 + u - v)*MyPow2(u) + 6*(-1 + u)*(-1 + v)*v*Cp(13)*MyInverse(-2\
 
  + u + v)*MyPow2(u) + 6*(-1 + v)*v*((-1 + u)*Cp(13) + (-1 + v)*Cp(14))*MyInverse(-2 + u + v)*MyPow2(u) - 6*(-1 +\
 
  u)*u*Cp(1)*MyPow2(-1 + v) + 6*(-1 + u)*u*Cp(2)*MyPow2(-1 + v) + 6*(-1 + u)*u*v*Cp(8)*MyInverse(-1 + u - v)*MyPow2(-1\
 
  + v) - 6*(-1 + u)*u*((-1 + u)*Cp(7) - v*Cp(8))*MyInverse(-1 + u - v)*MyPow2(-1 + v) + 6*(-1 +\
 
  u)*u*v*Cp(6)*MyInverse(u + v)*MyPow2(-1 + v) + 6*(-1 + u)*u*(u*Cp(5) + v*Cp(6))*MyInverse(u + v)*MyPow2(-1 + v) +\
 
  3*Cp(0)*MyPow2(-1 + u)*MyPow2(-1 + v) - 3*Cp(1)*MyPow2(-1 + u)*MyPow2(-1 + v) - 3*Cp(4)*MyPow2(-1 + u)*MyPow2(-1 + v)\
 
  + 3*u*Cp(5)*MyInverse(u + v)*MyPow2(-1 + u)*MyPow2(-1 + v) + 3*v*Cp(6)*MyInverse(u + v)*MyPow2(-1 + u)*MyPow2(-1 + v)\
 
  + 3*(u*Cp(5) + v*Cp(6))*MyInverse(u + v)*MyPow2(-1 + u)*MyPow2(-1 + v) + 3*Cp(2)*MyPow2(u)*MyPow2(-1 + v) -\
 
  3*Cp(3)*MyPow2(u)*MyPow2(-1 + v) + 3*Cp(9)*MyPow2(u)*MyPow2(-1 + v) - 3*(-1 + u)*Cp(7)*MyInverse(-1 + u -\
 
  v)*MyPow2(u)*MyPow2(-1 + v) + 3*v*Cp(8)*MyInverse(-1 + u - v)*MyPow2(u)*MyPow2(-1 + v) - 3*((-1 + u)*Cp(7) -\
 
  v*Cp(8))*MyInverse(-1 + u - v)*MyPow2(u)*MyPow2(-1 + v) + 6*(-1 + u)*u*Cp(17)*MyPow2(v) - 6*(-1 +\
 
  u)*u*Cp(18)*MyPow2(v) + 6*(-1 + u)*u*(-1 + v)*Cp(12)*MyInverse(1 + u - v)*MyPow2(v) - 6*(-1 + u)*u*(u*Cp(11) - (-1 +\
 
  v)*Cp(12))*MyInverse(1 + u - v)*MyPow2(v) + 6*(-1 + u)*u*(-1 + v)*Cp(14)*MyInverse(-2 + u + v)*MyPow2(v) + 6*(-1 +\
 
  u)*u*((-1 + u)*Cp(13) + (-1 + v)*Cp(14))*MyInverse(-2 + u + v)*MyPow2(v) + 3*Cp(10)*MyPow2(-1 + u)*MyPow2(v) -\
 
  3*Cp(16)*MyPow2(-1 + u)*MyPow2(v) + 3*Cp(17)*MyPow2(-1 + u)*MyPow2(v) - 3*u*Cp(11)*MyInverse(1 + u - v)*MyPow2(-1 +\
 
  u)*MyPow2(v) + 3*(-1 + v)*Cp(12)*MyInverse(1 + u - v)*MyPow2(-1 + u)*MyPow2(v) - 3*(u*Cp(11) - (-1 +\
 
  v)*Cp(12))*MyInverse(1 + u - v)*MyPow2(-1 + u)*MyPow2(v) - 3*Cp(15)*MyPow2(u)*MyPow2(v) -\
 
  3*Cp(18)*MyPow2(u)*MyPow2(v) + 3*Cp(19)*MyPow2(u)*MyPow2(v) + 3*(-1 + u)*Cp(13)*MyInverse(-2 + u +\
 
  v)*MyPow2(u)*MyPow2(v) + 3*(-1 + v)*Cp(14)*MyInverse(-2 + u + v)*MyPow2(u)*MyPow2(v) + 3*((-1 + u)*Cp(13) + (-1 +\
 
  v)*Cp(14))*MyInverse(-2 + u + v)*MyPow2(u)*MyPow2(v) + 6*(-1 + u)*(-1 + v)*v*((-1 + u)*Cp(7) -\
 
  v*Cp(8))*MyPow2(u)*MyPow2(MyInverse(-1 + u - v)) - 3*(-1 + u)*v*Cp(8)*MyPow2(u)*MyPow2(-1 + v)*MyPow2(MyInverse(-1 +\
 
  u - v)) + 3*(-1 + u)*((-1 + u)*Cp(7) - v*Cp(8))*MyPow2(u)*MyPow2(-1 + v)*MyPow2(MyInverse(-1 + u - v)) + 6*u*(-1 +\
 
  v)*v*(u*Cp(11) - (-1 + v)*Cp(12))*MyPow2(-1 + u)*MyPow2(MyInverse(1 + u - v)) - 6*(-1 + u)*u*(-1 + v)*(u*Cp(11) - (-1\
 
  + v)*Cp(12))*MyPow2(v)*MyPow2(MyInverse(1 + u - v)) - 3*u*(-1 + v)*Cp(11)*MyPow2(-1 + u)*MyPow2(v)*MyPow2(MyInverse(1\
 
  + u - v)) - 3*u*(-1 + v)*Cp(12)*MyPow2(-1 + u)*MyPow2(v)*MyPow2(MyInverse(1 + u - v)) + 3*u*(u*Cp(11) - (-1 +\
 
  v)*Cp(12))*MyPow2(-1 + u)*MyPow2(v)*MyPow2(MyInverse(1 + u - v)) - 3*(-1 + v)*(u*Cp(11) - (-1 + v)*Cp(12))*MyPow2(-1\
 
  + u)*MyPow2(v)*MyPow2(MyInverse(1 + u - v)) - 6*(-1 + u)*u*v*((-1 + u)*Cp(7) - v*Cp(8))*MyPow2(-1 +\
 
  v)*MyPow2(MyInverse(1 - u + v)) - 3*(-1 + u)*v*Cp(7)*MyPow2(u)*MyPow2(-1 + v)*MyPow2(MyInverse(1 - u + v)) - 3*v*((-1\
 
  + u)*Cp(7) - v*Cp(8))*MyPow2(u)*MyPow2(-1 + v)*MyPow2(MyInverse(1 - u + v)) - 6*(-1 + u)*(-1 + v)*v*((-1 + u)*Cp(13)\
 
  + (-1 + v)*Cp(14))*MyPow2(u)*MyPow2(MyInverse(-2 + u + v)) - 6*(-1 + u)*u*(-1 + v)*((-1 + u)*Cp(13) + (-1 +\
 
  v)*Cp(14))*MyPow2(v)*MyPow2(MyInverse(-2 + u + v)) - 3*(-1 + u)*(-1 +\
 
  v)*Cp(13)*MyPow2(u)*MyPow2(v)*MyPow2(MyInverse(-2 + u + v)) - 3*(-1 + u)*(-1 +\
 
  v)*Cp(14)*MyPow2(u)*MyPow2(v)*MyPow2(MyInverse(-2 + u + v)) - 3*(-1 + u)*((-1 + u)*Cp(13) + (-1 +\
 
  v)*Cp(14))*MyPow2(u)*MyPow2(v)*MyPow2(MyInverse(-2 + u + v)) - 3*(-1 + v)*((-1 + u)*Cp(13) + (-1 +\
 
  v)*Cp(14))*MyPow2(u)*MyPow2(v)*MyPow2(MyInverse(-2 + u + v)) - 6*u*(-1 + v)*v*(u*Cp(5) + v*Cp(6))*MyPow2(-1 +\
 
  u)*MyPow2(MyInverse(u + v)) - 6*(-1 + u)*u*v*(u*Cp(5) + v*Cp(6))*MyPow2(-1 + v)*MyPow2(MyInverse(u + v)) -\
 
  3*u*v*Cp(5)*MyPow2(-1 + u)*MyPow2(-1 + v)*MyPow2(MyInverse(u + v)) - 3*u*v*Cp(6)*MyPow2(-1 + u)*MyPow2(-1 +\
 
  v)*MyPow2(MyInverse(u + v)) - 3*u*(u*Cp(5) + v*Cp(6))*MyPow2(-1 + u)*MyPow2(-1 + v)*MyPow2(MyInverse(u + v)) -\
 
  3*v*(u*Cp(5) + v*Cp(6))*MyPow2(-1 + u)*MyPow2(-1 + v)*MyPow2(MyInverse(u + v)) + 6*u*(-1 + v)*(u*Cp(11) - (-1 +\
 
  v)*Cp(12))*MyPow2(-1 + u)*MyPow2(v)*MyPow3(MyInverse(1 + u - v)) - 6*(-1 + u)*v*((-1 + u)*Cp(7) -\
 
  v*Cp(8))*MyPow2(u)*MyPow2(-1 + v)*MyPow3(MyInverse(1 - u + v)) + 6*(-1 + u)*(-1 + v)*((-1 + u)*Cp(13) + (-1 +\
 
  v)*Cp(14))*MyPow2(u)*MyPow2(v)*MyPow3(MyInverse(-2 + u + v)) + 6*u*v*(u*Cp(5) + v*Cp(6))*MyPow2(-1 + u)*MyPow2(-1 +\
 
  v)*MyPow3(MyInverse(u + v))) ;

      H11 = 3*(-6*u*(-1 + v)*Cp(1)*MyPow2(-1 + u) + 6*u*v*Cp(17)*MyPow2(-1 + u) + 12*u*(-1 + v)*v*Cp(12)*MyInverse(1 +\
 
  u - v)*MyPow2(-1 + u) - 6*u*(-1 + v)*(u*Cp(11) - (-1 + v)*Cp(12))*MyInverse(1 + u - v)*MyPow2(-1 + u) -\
 
  12*u*v*(u*Cp(11) - (-1 + v)*Cp(12))*MyInverse(1 + u - v)*MyPow2(-1 + u) + 12*u*(-1 + v)*v*Cp(6)*MyInverse(u +\
 
  v)*MyPow2(-1 + u) + 12*u*(-1 + v)*(u*Cp(5) + v*Cp(6))*MyInverse(u + v)*MyPow2(-1 + u) + 6*u*v*(u*Cp(5) +\
 
  v*Cp(6))*MyInverse(u + v)*MyPow2(-1 + u) + 6*(-1 + u)*(-1 + v)*Cp(2)*MyPow2(u) - 6*(-1 + u)*v*Cp(18)*MyPow2(u) +\
 
  12*(-1 + u)*(-1 + v)*v*Cp(8)*MyInverse(-1 + u - v)*MyPow2(u) - 12*(-1 + u)*(-1 + v)*((-1 + u)*Cp(7) -\
 
  v*Cp(8))*MyInverse(-1 + u - v)*MyPow2(u) - 6*(-1 + u)*v*((-1 + u)*Cp(7) - v*Cp(8))*MyInverse(-1 + u - v)*MyPow2(u) +\
 
  12*(-1 + u)*(-1 + v)*v*Cp(14)*MyInverse(-2 + u + v)*MyPow2(u) + 6*(-1 + u)*(-1 + v)*((-1 + u)*Cp(13) + (-1 +\
 
  v)*Cp(14))*MyInverse(-2 + u + v)*MyPow2(u) + 12*(-1 + u)*v*((-1 + u)*Cp(13) + (-1 + v)*Cp(14))*MyInverse(-2 + u +\
 
  v)*MyPow2(u) + 6*u*Cp(6)*MyInverse(u + v)*MyPow2(-1 + u)*MyPow2(-1 + v) + 6*(-1 + u)*Cp(8)*MyInverse(-1 + u -\
 
  v)*MyPow2(u)*MyPow2(-1 + v) + 6*u*Cp(12)*MyInverse(1 + u - v)*MyPow2(-1 + u)*MyPow2(v) + 6*(-1 +\
 
  u)*Cp(14)*MyInverse(-2 + u + v)*MyPow2(u)*MyPow2(v) - 6*(-1 + u)*(-1 + v)*v*((-1 + u)*Cp(7) -\
 
  v*Cp(8))*MyPow2(u)*MyPow2(MyInverse(-1 + u - v)) + 3*(-1 + u)*v*Cp(8)*MyPow2(u)*MyPow2(-1 + v)*MyPow2(MyInverse(-1 +\
 
  u - v)) - 3*(-1 + u)*((-1 + u)*Cp(7) - v*Cp(8))*MyPow2(u)*MyPow2(-1 + v)*MyPow2(MyInverse(-1 + u - v)) - 12*u*(-1 +\
 
  v)*v*(u*Cp(11) - (-1 + v)*Cp(12))*MyPow2(-1 + u)*MyPow2(MyInverse(1 + u - v)) + 6*u*(-1 + v)*Cp(12)*MyPow2(-1 +\
 
  u)*MyPow2(v)*MyPow2(MyInverse(1 + u - v)) - 6*u*(u*Cp(11) - (-1 + v)*Cp(12))*MyPow2(-1 +\
 
  u)*MyPow2(v)*MyPow2(MyInverse(1 + u - v)) - 6*(-1 + u)*(-1 + v)*v*((-1 + u)*Cp(7) -\
 
  v*Cp(8))*MyPow2(u)*MyPow2(MyInverse(1 - u + v)) + 3*(-1 + u)*v*Cp(8)*MyPow2(u)*MyPow2(-1 + v)*MyPow2(MyInverse(1 - u\
 
  + v)) - 3*(-1 + u)*((-1 + u)*Cp(7) - v*Cp(8))*MyPow2(u)*MyPow2(-1 + v)*MyPow2(MyInverse(1 - u + v)) - 12*(-1 + u)*(-1\
 
  + v)*v*((-1 + u)*Cp(13) + (-1 + v)*Cp(14))*MyPow2(u)*MyPow2(MyInverse(-2 + u + v)) - 6*(-1 + u)*(-1 +\
 
  v)*Cp(14)*MyPow2(u)*MyPow2(v)*MyPow2(MyInverse(-2 + u + v)) - 6*(-1 + u)*((-1 + u)*Cp(13) + (-1 +\
 
  v)*Cp(14))*MyPow2(u)*MyPow2(v)*MyPow2(MyInverse(-2 + u + v)) - 12*u*(-1 + v)*v*(u*Cp(5) + v*Cp(6))*MyPow2(-1 +\
 
  u)*MyPow2(MyInverse(u + v)) - 6*u*v*Cp(6)*MyPow2(-1 + u)*MyPow2(-1 + v)*MyPow2(MyInverse(u + v)) - 6*u*(u*Cp(5) +\
 
  v*Cp(6))*MyPow2(-1 + u)*MyPow2(-1 + v)*MyPow2(MyInverse(u + v)) + 2*(-1 + v)*Cp(0)*MyPow3(-1 + u) - 4*(-1 +\
 
  v)*Cp(4)*MyPow3(-1 + u) - 2*v*Cp(4)*MyPow3(-1 + u) + 2*(-1 + v)*Cp(10)*MyPow3(-1 + u) + 4*v*Cp(10)*MyPow3(-1 + u) -\
 
  2*v*Cp(16)*MyPow3(-1 + u) - 2*(-1 + v)*Cp(3)*MyPow3(u) + 4*(-1 + v)*Cp(9)*MyPow3(u) + 2*v*Cp(9)*MyPow3(u) - 2*(-1 +\
 
  v)*Cp(15)*MyPow3(u) - 4*v*Cp(15)*MyPow3(u) + 2*v*Cp(19)*MyPow3(u) - 6*u*(-1 + v)*(u*Cp(11) - (-1 +\
 
  v)*Cp(12))*MyPow2(-1 + u)*MyPow2(v)*MyPow3(MyInverse(1 + u - v)) + 6*(-1 + u)*v*((-1 + u)*Cp(7) -\
 
  v*Cp(8))*MyPow2(u)*MyPow2(-1 + v)*MyPow3(MyInverse(1 - u + v)) + 6*(-1 + u)*(-1 + v)*((-1 + u)*Cp(13) + (-1 +\
 
  v)*Cp(14))*MyPow2(u)*MyPow2(v)*MyPow3(MyInverse(-2 + u + v)) + 6*u*v*(u*Cp(5) + v*Cp(6))*MyPow2(-1 + u)*MyPow2(-1 +\
 
  v)*MyPow3(MyInverse(u + v))) ;
    }
  H[0][0] = H00;
  H[0][1] = H01;
  H[1][0] = H10;
  H[1][1] = H11;
}

#undef NEAR

}
