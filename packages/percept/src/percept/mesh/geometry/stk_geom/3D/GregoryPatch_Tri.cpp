// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "GregoryPatch.hpp"

namespace percept {

// our tri patches are degree-elevated cubics on the boundary, quartics in the interior
//   (and with extra Gregory points in interior) - since we store the quartics,
//   this function extracts the cubics on the boundary
void  GregoryPatch::
extract_degree_lower_on_boundary(const MDArray& qh, MDArray& q)
{
  
      q(0,0) = qh(0,0) ;

      q(0,1) = qh(0,1) ;

      q(0,2) = qh(0,2) ;

      q(1,0) = (-qh(0,0) + 4*qh(1,0))/3. ;

      q(1,1) = (-qh(0,1) + 4*qh(1,1))/3. ;

      q(1,2) = (-qh(0,2) + 4*qh(1,2))/3. ;

      q(2,0) = (qh(0,0) - 4*qh(1,0) + 6*qh(2,0))/3. ;

      q(2,1) = (qh(0,1) - 4*qh(1,1) + 6*qh(2,1))/3. ;

      q(2,2) = (qh(0,2) - 4*qh(1,2) + 6*qh(2,2))/3. ;

      q(3,0) = -qh(0,0) + 4*qh(1,0) - 6*qh(2,0) + 4*qh(3,0) ;

      q(3,1) = -qh(0,1) + 4*qh(1,1) - 6*qh(2,1) + 4*qh(3,1) ;

      q(3,2) = -qh(0,2) + 4*qh(1,2) - 6*qh(2,2) + 4*qh(3,2) ;
}

void  GregoryPatch::
degree_elevate(MDArray& q, MDArray& qh)
{
  
      qh(0,0) = q(0,0) ;

      qh(0,1) = q(0,1) ;

      qh(0,2) = q(0,2) ;

      qh(1,0) = q(0,0)/4. + (3*q(1,0))/4. ;

      qh(1,1) = q(0,1)/4. + (3*q(1,1))/4. ;

      qh(1,2) = q(0,2)/4. + (3*q(1,2))/4. ;

      qh(2,0) = q(1,0)/2. + q(2,0)/2. ;

      qh(2,1) = q(1,1)/2. + q(2,1)/2. ;

      qh(2,2) = q(1,2)/2. + q(2,2)/2. ;

      qh(3,0) = (3*q(2,0))/4. + q(3,0)/4. ;

      qh(3,1) = (3*q(2,1))/4. + q(3,1)/4. ;

      qh(3,2) = (3*q(2,2))/4. + q(3,2)/4. ;

      qh(4,0) = q(3,0) ;

      qh(4,1) = q(3,1) ;

      qh(4,2) = q(3,2) ;
}
/// @param Cp is for a single coordinate (x,y, or z)
double  GregoryPatch::
evalTri(double u, double v, const MDArray& Cp)
{
  /*
     b[0, 0][u, v]
b[1, 0][u, v]
b[2, 0][u, v]
b[3, 0][u, v]
b[4, 0][u, v]
b[0, 1][u, v]
b[1, 1, 1, 0][u, v]
b[1, 1, 0, 1][u, v]
b[2, 1, 1, 0][u, v]
b[2, 1, 1, 1][u, v]
b[3, 1][u, v]
b[0, 2][u, v]
b[1, 2, 1, 1][u, v]
b[1, 2, 0, 1][u, v]
b[2, 2][u, v]
b[0, 3][u, v]
b[1, 3][u, v]
b[0, 4][u, v]

     0  : b[0, 0][u, v]
     1  : b[1, 0][u, v]
     2  : b[2, 0][u, v]
     3  : b[3, 0][u, v]
     4  : b[4, 0][u, v]
     5  : b[0, 1][u, v]
     6  : b[1, 1, 1, 0][u, v]
     7  : b[1, 1, 0, 1][u, v]
     8  : b[2, 1, 1, 0][u, v]
     9  : b[2, 1, 1, 1][u, v]
     10 : b[3, 1][u, v]
     11 : b[0, 2][u, v]
     12 : b[1, 2, 1, 1][u, v]
     13 : b[1, 2, 0, 1][u, v]
     14 : b[2, 2][u, v]
     15 : b[0, 3][u, v]
     16 : b[1, 3][u, v]
     17 : b[0, 4][u, v]

     edges
     edge 0: {0,1,2,3,4}
     edge 1: {4,10,14,16,17}
     edge 2: {17,15,11,5,0}

     ribbons: (convention is interior is on the left, 1st 4 are "q" boundary points, 2nd 4 are interior "p" points)

     v = 0 (edge 0)
     { {0, DL(0,1), DL(4,3), 4}, {DL(0,5), 6, 8, DL(4,10)} }
     w = 0 (edge 1)
     { {4, DL(4,10), DL(17,16), 17}, {DL(4,3), 9, 12, DL(17,15) } }
     u = 0 (edge 2)
     { {17, DL(17,15), DL(0,5), 0}, {DL(17,16), 13, 7, DL(0,1)}}

     Note: DL(q[0],q[1]) = ql[1] = degree lower of quartic given by q[i],i=0,4 = 1/3 (4 q[1] - q[0]) (Note: ql[0] = q[0])
     Note: DH(ql[0],ql[1]) = qh[1] = degree higher of cubic ql[i],i=0,3 = (ql[0] + 3 ql[1]) / 4
   */
  double p = 0;
#define NEAR(u,val) (std::fabs((u)-(val)) < 1.e-5)
  if (NEAR(u,0))
    {
      
      p = 6*Cp(11)*MyPow2(-1 + v)*MyPow2(v) - 4*v*Cp(5)*MyPow3(-1 + v) - 4*(-1 + v)*Cp(15)*MyPow3(v) + Cp(0)*MyPow4(-1\
 
  + v) + Cp(17)*MyPow4(v) ;
    }
  else if (NEAR(u , 1-v))
    {
      
      p = 6*Cp(14)*MyPow2(1 - u)*MyPow2(u) + 4*u*Cp(16)*MyPow3(1 - u) + 4*(1 - u)*Cp(10)*MyPow3(u) + Cp(17)*MyPow4(1 -\
 
  u) + Cp(4)*MyPow4(u) ;
    }
  else if (NEAR(v, 0))
    {
      
      p = 6*Cp(2)*MyPow2(-1 + u)*MyPow2(u) - 4*u*Cp(1)*MyPow3(-1 + u) - 4*(-1 + u)*Cp(3)*MyPow3(u) + Cp(0)*MyPow4(-1 +\
 
  u) + Cp(4)*MyPow4(u) ;
    }
  else
    {
      
      p = -12*v*(-1 + u + v)*((-1 + u + v)*Cp(8) - v*Cp(9))*MyInverse(-1 + u)*MyPow2(u) - 12*u*(-1 + u +\
 
  v)*(-(u*Cp(12)) + (-1 + u + v)*Cp(13))*MyInverse(-1 + v)*MyPow2(v) + 6*Cp(14)*MyPow2(u)*MyPow2(v) + 12*u*v*(u*Cp(6) +\
 
  v*Cp(7))*MyInverse(u + v)*MyPow2(-1 + u + v) + 6*Cp(2)*MyPow2(u)*MyPow2(-1 + u + v) + 6*Cp(11)*MyPow2(v)*MyPow2(-1 +\
 
  u + v) - 4*(-1 + u + v)*Cp(3)*MyPow3(u) + 4*v*Cp(10)*MyPow3(u) - 4*(-1 + u + v)*Cp(15)*MyPow3(v) +\
 
  4*u*Cp(16)*MyPow3(v) - 4*u*Cp(1)*MyPow3(-1 + u + v) - 4*v*Cp(5)*MyPow3(-1 + u + v) + Cp(4)*MyPow4(u) +\
 
  Cp(17)*MyPow4(v) + Cp(0)*MyPow4(-1 + u + v) ;
    }
  return p;
}

/// @param Cp is for a single coordinate (x,y, or z)
void  GregoryPatch::
evalTriGrad(double u, double v, const MDArray& Cp, double *grad)
{
  double gu = 0.0, gv=0.0;
  if (NEAR(u,0))
    {
      
      gu = 4*(-3*v*Cp(5)*MyPow2(-1 + v) + 3*v*Cp(7)*MyPow2(-1 + v) + 3*(-1 + v)*Cp(11)*MyPow2(v) - 3*(-1 +\
 
  v)*Cp(13)*MyPow2(v) + Cp(0)*MyPow3(-1 + v) - Cp(1)*MyPow3(-1 + v) - Cp(15)*MyPow3(v) + Cp(16)*MyPow3(v)) ;

      gv = 4*(-3*v*Cp(5)*MyPow2(-1 + v) + 3*v*Cp(11)*MyPow2(-1 + v) + 3*(-1 + v)*Cp(11)*MyPow2(v) - 3*(-1 +\
 
  v)*Cp(15)*MyPow2(v) + Cp(0)*MyPow3(-1 + v) - Cp(5)*MyPow3(-1 + v) - Cp(15)*MyPow3(v) + Cp(17)*MyPow3(v)) ;
    }
  else if (NEAR(u , 1-v))
    {
      
      gu = 4*(-3*v*Cp(9)*MyPow2(-1 + v) + 3*v*Cp(10)*MyPow2(-1 + v) + 3*(-1 + v)*Cp(12)*MyPow2(v) - 3*(-1 +\
 
  v)*Cp(14)*MyPow2(v) + Cp(3)*MyPow3(-1 + v) - Cp(4)*MyPow3(-1 + v) - Cp(15)*MyPow3(v) + Cp(16)*MyPow3(v)) ;

      gv = 4*(-3*v*Cp(9)*MyPow2(-1 + v) + 3*v*Cp(14)*MyPow2(-1 + v) + 3*(-1 + v)*Cp(12)*MyPow2(v) - 3*(-1 +\
 
  v)*Cp(16)*MyPow2(v) + Cp(3)*MyPow3(-1 + v) - Cp(10)*MyPow3(-1 + v) - Cp(15)*MyPow3(v) + Cp(17)*MyPow3(v)) ;
    }
  else if (NEAR(v, 0))
    {
      
      gu = 4*(-3*u*Cp(1)*MyPow2(-1 + u) + 3*u*Cp(2)*MyPow2(-1 + u) + 3*(-1 + u)*Cp(2)*MyPow2(u) - 3*(-1 +\
 
  u)*Cp(3)*MyPow2(u) + Cp(0)*MyPow3(-1 + u) - Cp(1)*MyPow3(-1 + u) - Cp(3)*MyPow3(u) + Cp(4)*MyPow3(u)) ;

      gv = 4*(-3*u*Cp(1)*MyPow2(-1 + u) + 3*u*Cp(6)*MyPow2(-1 + u) + 3*(-1 + u)*Cp(2)*MyPow2(u) - 3*(-1 +\
 
  u)*Cp(8)*MyPow2(u) + Cp(0)*MyPow3(-1 + u) - Cp(5)*MyPow3(-1 + u) - Cp(3)*MyPow3(u) + Cp(10)*MyPow3(u)) ;
    }
  else
    {
      
      gu = 4*(-6*u*v*(-1 + u + v)*((-1 + u + v)*Cp(8) - v*Cp(9))*MyInverse(-1 + u) + 6*u*v*(-1 + u + v)*(u*Cp(6) +\
 
  v*Cp(7))*MyInverse(u + v) + 3*(-1 + u + v)*Cp(2)*MyPow2(u) - 3*(-1 + u + v)*Cp(3)*MyPow2(u) + 3*v*Cp(10)*MyPow2(u) -\
 
  3*v*(-1 + u + v)*Cp(8)*MyInverse(-1 + u)*MyPow2(u) - 3*v*((-1 + u + v)*Cp(8) - v*Cp(9))*MyInverse(-1 + u)*MyPow2(u) +\
 
  3*(-1 + u + v)*Cp(11)*MyPow2(v) + 3*u*Cp(14)*MyPow2(v) + 3*u*(-1 + u + v)*(Cp(12) - Cp(13))*MyInverse(-1 +\
 
  v)*MyPow2(v) - 3*u*(-(u*Cp(12)) + (-1 + u + v)*Cp(13))*MyInverse(-1 + v)*MyPow2(v) - 3*(-1 + u + v)*(-(u*Cp(12)) +\
 
  (-1 + u + v)*Cp(13))*MyInverse(-1 + v)*MyPow2(v) - 3*u*Cp(1)*MyPow2(-1 + u + v) + 3*u*Cp(2)*MyPow2(-1 + u + v) -\
 
  3*v*Cp(5)*MyPow2(-1 + u + v) + 3*u*v*Cp(6)*MyInverse(u + v)*MyPow2(-1 + u + v) + 3*v*(u*Cp(6) + v*Cp(7))*MyInverse(u\
 
  + v)*MyPow2(-1 + u + v) + 3*v*(-1 + u + v)*((-1 + u + v)*Cp(8) - v*Cp(9))*MyPow2(u)*MyPow2(MyInverse(-1 + u)) -\
 
  3*u*v*(u*Cp(6) + v*Cp(7))*MyPow2(-1 + u + v)*MyPow2(MyInverse(u + v)) - Cp(3)*MyPow3(u) + Cp(4)*MyPow3(u) -\
 
  Cp(15)*MyPow3(v) + Cp(16)*MyPow3(v) + Cp(0)*MyPow3(-1 + u + v) - Cp(1)*MyPow3(-1 + u + v)) ;

      gv = 4*(-6*u*v*(-1 + u + v)*(-(u*Cp(12)) + (-1 + u + v)*Cp(13))*MyInverse(-1 + v) + 6*u*v*(-1 + u + v)*(u*Cp(6) +\
 
  v*Cp(7))*MyInverse(u + v) + 3*(-1 + u + v)*Cp(2)*MyPow2(u) + 3*v*Cp(14)*MyPow2(u) - 3*v*(-1 + u + v)*(Cp(8) -\
 
  Cp(9))*MyInverse(-1 + u)*MyPow2(u) - 3*v*((-1 + u + v)*Cp(8) - v*Cp(9))*MyInverse(-1 + u)*MyPow2(u) - 3*(-1 + u +\
 
  v)*((-1 + u + v)*Cp(8) - v*Cp(9))*MyInverse(-1 + u)*MyPow2(u) + 3*(-1 + u + v)*Cp(11)*MyPow2(v) - 3*(-1 + u +\
 
  v)*Cp(15)*MyPow2(v) + 3*u*Cp(16)*MyPow2(v) - 3*u*(-1 + u + v)*Cp(13)*MyInverse(-1 + v)*MyPow2(v) - 3*u*(-(u*Cp(12)) +\
 
  (-1 + u + v)*Cp(13))*MyInverse(-1 + v)*MyPow2(v) - 3*u*Cp(1)*MyPow2(-1 + u + v) - 3*v*Cp(5)*MyPow2(-1 + u + v) +\
 
  3*v*Cp(11)*MyPow2(-1 + u + v) + 3*u*v*Cp(7)*MyInverse(u + v)*MyPow2(-1 + u + v) + 3*u*(u*Cp(6) + v*Cp(7))*MyInverse(u\
 
  + v)*MyPow2(-1 + u + v) + 3*u*(-1 + u + v)*(-(u*Cp(12)) + (-1 + u + v)*Cp(13))*MyPow2(v)*MyPow2(MyInverse(-1 + v)) -\
 
  3*u*v*(u*Cp(6) + v*Cp(7))*MyPow2(-1 + u + v)*MyPow2(MyInverse(u + v)) - Cp(3)*MyPow3(u) + Cp(10)*MyPow3(u) -\
 
  Cp(15)*MyPow3(v) + Cp(17)*MyPow3(v) + Cp(0)*MyPow3(-1 + u + v) - Cp(5)*MyPow3(-1 + u + v)) ;
    }
  grad[0] = gu;
  grad[1] = gv;
}

void  GregoryPatch::
evalTriHessian(double u, double v, const MDArray& Cp, double H[2][2])
{
  double H00, H01, H10, H11;
  if (NEAR(u,0))
    {
      
      H00 = 4*(-6*(-1 + v)*v*Cp(5) + 12*(-1 + v)*v*Cp(7) + 6*(-1 + v)*v*((-1 + v)*Cp(8) - v*Cp(9)) + 3*Cp(0)*MyPow2(-1\
 
  + v) - 6*Cp(1)*MyPow2(-1 + v) + 3*Cp(2)*MyPow2(-1 + v) + 6*Cp(6)*MyPow2(-1 + v) - 6*Cp(7)*MyPow2(-1 + v) +\
 
  3*Cp(11)*MyPow2(v) + 3*(Cp(12) - Cp(13))*MyPow2(v) - 6*Cp(13)*MyPow2(v) - 3*(-Cp(12) + Cp(13))*MyPow2(v) +\
 
  3*Cp(14)*MyPow2(v)) ;

      H01 = 4*(-6*(-1 + v)*v*Cp(5) + 6*(-1 + v)*v*Cp(7) + 6*(-1 + v)*v*Cp(11) - 6*(-1 + v)*v*Cp(13) + 3*Cp(0)*MyPow2(-1\
 
  + v) - 3*Cp(1)*MyPow2(-1 + v) - 3*Cp(5)*MyPow2(-1 + v) + 3*Cp(7)*MyPow2(-1 + v) + 3*Cp(11)*MyPow2(v) -\
 
  3*Cp(13)*MyPow2(v) - 3*Cp(15)*MyPow2(v) + 3*Cp(16)*MyPow2(v)) ;

      H10 = 4*(-6*(-1 + v)*v*Cp(5) + 6*(-1 + v)*v*Cp(7) + 6*(-1 + v)*v*Cp(11) - 6*(-1 + v)*v*Cp(13) + 3*Cp(0)*MyPow2(-1\
 
  + v) - 3*Cp(1)*MyPow2(-1 + v) - 3*Cp(5)*MyPow2(-1 + v) + 3*Cp(7)*MyPow2(-1 + v) + 3*Cp(11)*MyPow2(v) -\
 
  3*Cp(13)*MyPow2(v) - 3*Cp(15)*MyPow2(v) + 3*Cp(16)*MyPow2(v)) ;

      H11 = 4*(-6*(-1 + v)*v*Cp(5) + 12*(-1 + v)*v*Cp(11) - 6*(-1 + v)*v*Cp(15) + 3*Cp(0)*MyPow2(-1 + v) -\
 
  6*Cp(5)*MyPow2(-1 + v) + 3*Cp(11)*MyPow2(-1 + v) + 3*Cp(11)*MyPow2(v) - 6*Cp(15)*MyPow2(v) + 3*Cp(17)*MyPow2(v)) ;
    }
  else if (NEAR(u , 1-v))
    {
      
      H00 = 4*(-6*(-1 + v)*v*(-((-1 + v)*Cp(6)) + v*Cp(7)) + 12*(-1 + v)*v*Cp(9) - 6*(-1 + v)*v*Cp(10) +\
 
  3*Cp(2)*MyPow2(-1 + v) - 6*Cp(3)*MyPow2(-1 + v) + 3*Cp(4)*MyPow2(-1 + v) + 6*Cp(8)*MyPow2(-1 + v) - 6*Cp(9)*MyPow2(-1\
 
  + v) + 3*Cp(11)*MyPow2(v) - 6*Cp(12)*MyPow2(v) - 3*(Cp(12) - Cp(13))*MyPow2(v) + 3*(-Cp(12) + Cp(13))*MyPow2(v) +\
 
  3*Cp(14)*MyPow2(v)) ;

      H01 = 4*(-6*(-1 + v)*v*(-((-1 + v)*Cp(6)) + v*Cp(7)) + 6*(-1 + v)*v*Cp(9) + 6*(-1 + v)*v*Cp(12) - 6*(-1 +\
 
  v)*v*Cp(14) + 3*Cp(2)*MyPow2(-1 + v) - 3*Cp(3)*MyPow2(-1 + v) + 3*Cp(8)*MyPow2(-1 + v) + 3*(Cp(8) - Cp(9))*MyPow2(-1\
 
  + v) - 6*Cp(9)*MyPow2(-1 + v) + 3*Cp(10)*MyPow2(-1 + v) + 3*Cp(11)*MyPow2(v) - 6*Cp(12)*MyPow2(v) - 3*(Cp(12) -\
 
  Cp(13))*MyPow2(v) + 3*Cp(13)*MyPow2(v) - 3*Cp(15)*MyPow2(v) + 3*Cp(16)*MyPow2(v)) ;

      H10 = 4*(-6*(-1 + v)*v*(-((-1 + v)*Cp(6)) + v*Cp(7)) + 6*(-1 + v)*v*Cp(9) + 6*(-1 + v)*v*Cp(12) - 6*(-1 +\
 
  v)*v*Cp(14) + 3*Cp(2)*MyPow2(-1 + v) - 3*Cp(3)*MyPow2(-1 + v) + 3*Cp(8)*MyPow2(-1 + v) + 3*(Cp(8) - Cp(9))*MyPow2(-1\
 
  + v) - 6*Cp(9)*MyPow2(-1 + v) + 3*Cp(10)*MyPow2(-1 + v) + 3*Cp(11)*MyPow2(v) - 6*Cp(12)*MyPow2(v) +\
 
  3*Cp(13)*MyPow2(v) + 3*(-Cp(12) + Cp(13))*MyPow2(v) - 3*Cp(15)*MyPow2(v) + 3*Cp(16)*MyPow2(v)) ;

      H11 = 4*(-6*(-1 + v)*v*(-((-1 + v)*Cp(6)) + v*Cp(7)) + 12*(-1 + v)*v*Cp(12) - 6*(-1 + v)*v*Cp(16) +\
 
  3*Cp(2)*MyPow2(-1 + v) + 6*(Cp(8) - Cp(9))*MyPow2(-1 + v) - 6*Cp(9)*MyPow2(-1 + v) + 3*Cp(14)*MyPow2(-1 + v) +\
 
  3*Cp(11)*MyPow2(v) - 6*Cp(12)*MyPow2(v) + 6*Cp(13)*MyPow2(v) - 6*Cp(15)*MyPow2(v) + 3*Cp(17)*MyPow2(v)) ;
    }
  else if (NEAR(v, 0))
    {
      
      H00 = 4*(-6*(-1 + u)*u*Cp(1) + 12*(-1 + u)*u*Cp(2) - 6*(-1 + u)*u*Cp(3) + 3*Cp(0)*MyPow2(-1 + u) -\
 
  6*Cp(1)*MyPow2(-1 + u) + 3*Cp(2)*MyPow2(-1 + u) + 3*Cp(2)*MyPow2(u) - 6*Cp(3)*MyPow2(u) + 3*Cp(4)*MyPow2(u)) ;

      H01 = 4*(-6*(-1 + u)*u*Cp(1) + 6*(-1 + u)*u*Cp(2) + 6*(-1 + u)*u*Cp(6) - 6*(-1 + u)*u*Cp(8) + 3*Cp(0)*MyPow2(-1 +\
 
  u) - 3*Cp(1)*MyPow2(-1 + u) - 3*Cp(5)*MyPow2(-1 + u) + 3*Cp(6)*MyPow2(-1 + u) + 3*Cp(2)*MyPow2(u) - 3*Cp(3)*MyPow2(u)\
 
  - 3*Cp(8)*MyPow2(u) + 3*Cp(10)*MyPow2(u)) ;

      H10 = 4*(-6*(-1 + u)*u*Cp(1) + 6*(-1 + u)*u*Cp(2) + 6*(-1 + u)*u*Cp(6) - 6*(-1 + u)*u*Cp(8) + 3*Cp(0)*MyPow2(-1 +\
 
  u) - 3*Cp(1)*MyPow2(-1 + u) - 3*Cp(5)*MyPow2(-1 + u) + 3*Cp(6)*MyPow2(-1 + u) + 3*Cp(2)*MyPow2(u) - 3*Cp(3)*MyPow2(u)\
 
  - 3*Cp(8)*MyPow2(u) + 3*Cp(10)*MyPow2(u)) ;

      H11 = 4*(-6*(-1 + u)*u*Cp(1) + 12*(-1 + u)*u*Cp(6) + 6*(-1 + u)*u*(-(u*Cp(12)) + (-1 + u)*Cp(13)) +\
 
  3*Cp(0)*MyPow2(-1 + u) - 6*Cp(5)*MyPow2(-1 + u) - 6*Cp(6)*MyPow2(-1 + u) + 6*Cp(7)*MyPow2(-1 + u) +\
 
  3*Cp(11)*MyPow2(-1 + u) + 3*Cp(2)*MyPow2(u) - 6*Cp(8)*MyPow2(u) - 6*(Cp(8) - Cp(9))*MyPow2(u) + 3*Cp(14)*MyPow2(u)) ;
    }
  else
    {
      
      H00 = 4*(-6*u*(-1 + u + v)*Cp(1) + 12*u*(-1 + u + v)*Cp(2) - 6*u*(-1 + u + v)*Cp(3) - 6*v*(-1 + u + v)*Cp(5) +\
 
  6*u*v*Cp(10) - 12*u*v*(-1 + u + v)*Cp(8)*MyInverse(-1 + u) - 12*u*v*((-1 + u + v)*Cp(8) - v*Cp(9))*MyInverse(-1 + u)\
 
  - 6*v*(-1 + u + v)*((-1 + u + v)*Cp(8) - v*Cp(9))*MyInverse(-1 + u) + 12*u*v*(-1 + u + v)*Cp(6)*MyInverse(u + v) +\
 
  6*u*v*(u*Cp(6) + v*Cp(7))*MyInverse(u + v) + 12*v*(-1 + u + v)*(u*Cp(6) + v*Cp(7))*MyInverse(u + v) +\
 
  3*Cp(2)*MyPow2(u) - 6*Cp(3)*MyPow2(u) + 3*Cp(4)*MyPow2(u) - 6*v*Cp(8)*MyInverse(-1 + u)*MyPow2(u) +\
 
  3*Cp(11)*MyPow2(v) + 3*Cp(14)*MyPow2(v) + 3*u*(Cp(12) - Cp(13))*MyInverse(-1 + v)*MyPow2(v) + 3*(-1 + u + v)*(Cp(12)\
 
  - Cp(13))*MyInverse(-1 + v)*MyPow2(v) - 3*u*(-Cp(12) + Cp(13))*MyInverse(-1 + v)*MyPow2(v) - 3*(-1 + u + v)*(-Cp(12)\
 
  + Cp(13))*MyInverse(-1 + v)*MyPow2(v) - 6*(-(u*Cp(12)) + (-1 + u + v)*Cp(13))*MyInverse(-1 + v)*MyPow2(v) +\
 
  3*Cp(0)*MyPow2(-1 + u + v) - 6*Cp(1)*MyPow2(-1 + u + v) + 3*Cp(2)*MyPow2(-1 + u + v) + 6*v*Cp(6)*MyInverse(u +\
 
  v)*MyPow2(-1 + u + v) + 12*u*v*(-1 + u + v)*((-1 + u + v)*Cp(8) - v*Cp(9))*MyPow2(MyInverse(-1 + u)) + 6*v*(-1 + u +\
 
  v)*Cp(8)*MyPow2(u)*MyPow2(MyInverse(-1 + u)) + 6*v*((-1 + u + v)*Cp(8) - v*Cp(9))*MyPow2(u)*MyPow2(MyInverse(-1 + u))\
 
  - 12*u*v*(-1 + u + v)*(u*Cp(6) + v*Cp(7))*MyPow2(MyInverse(u + v)) - 6*u*v*Cp(6)*MyPow2(-1 + u +\
 
  v)*MyPow2(MyInverse(u + v)) - 6*v*(u*Cp(6) + v*Cp(7))*MyPow2(-1 + u + v)*MyPow2(MyInverse(u + v)) - 6*v*(-1 + u +\
 
  v)*((-1 + u + v)*Cp(8) - v*Cp(9))*MyPow2(u)*MyPow3(MyInverse(-1 + u)) + 6*u*v*(u*Cp(6) + v*Cp(7))*MyPow2(-1 + u +\
 
  v)*MyPow3(MyInverse(u + v))) ;

      H01 = 4*(-6*u*(-1 + u + v)*Cp(1) + 6*u*(-1 + u + v)*Cp(2) - 6*v*(-1 + u + v)*Cp(5) + 6*v*(-1 + u + v)*Cp(11) +\
 
  6*u*v*Cp(14) - 6*u*v*(-1 + u + v)*(Cp(8) - Cp(9))*MyInverse(-1 + u) - 6*u*v*((-1 + u + v)*Cp(8) -\
 
  v*Cp(9))*MyInverse(-1 + u) - 6*u*(-1 + u + v)*((-1 + u + v)*Cp(8) - v*Cp(9))*MyInverse(-1 + u) + 6*u*v*(-1 + u +\
 
  v)*(Cp(12) - Cp(13))*MyInverse(-1 + v) - 6*u*v*(-(u*Cp(12)) + (-1 + u + v)*Cp(13))*MyInverse(-1 + v) - 6*v*(-1 + u +\
 
  v)*(-(u*Cp(12)) + (-1 + u + v)*Cp(13))*MyInverse(-1 + v) + 6*u*v*(-1 + u + v)*Cp(6)*MyInverse(u + v) + 6*u*v*(-1 + u\
 
  + v)*Cp(7)*MyInverse(u + v) + 6*u*v*(u*Cp(6) + v*Cp(7))*MyInverse(u + v) + 6*u*(-1 + u + v)*(u*Cp(6) +\
 
  v*Cp(7))*MyInverse(u + v) + 6*v*(-1 + u + v)*(u*Cp(6) + v*Cp(7))*MyInverse(u + v) + 3*Cp(2)*MyPow2(u) -\
 
  3*Cp(3)*MyPow2(u) + 3*Cp(10)*MyPow2(u) - 3*v*Cp(8)*MyInverse(-1 + u)*MyPow2(u) - 3*(-1 + u + v)*Cp(8)*MyInverse(-1 +\
 
  u)*MyPow2(u) - 3*v*(Cp(8) - Cp(9))*MyInverse(-1 + u)*MyPow2(u) - 3*((-1 + u + v)*Cp(8) - v*Cp(9))*MyInverse(-1 +\
 
  u)*MyPow2(u) + 3*Cp(11)*MyPow2(v) - 3*Cp(15)*MyPow2(v) + 3*Cp(16)*MyPow2(v) + 3*u*(Cp(12) - Cp(13))*MyInverse(-1 +\
 
  v)*MyPow2(v) - 3*u*Cp(13)*MyInverse(-1 + v)*MyPow2(v) - 3*(-1 + u + v)*Cp(13)*MyInverse(-1 + v)*MyPow2(v) -\
 
  3*(-(u*Cp(12)) + (-1 + u + v)*Cp(13))*MyInverse(-1 + v)*MyPow2(v) + 3*Cp(0)*MyPow2(-1 + u + v) - 3*Cp(1)*MyPow2(-1 +\
 
  u + v) - 3*Cp(5)*MyPow2(-1 + u + v) + 3*u*Cp(6)*MyInverse(u + v)*MyPow2(-1 + u + v) + 3*v*Cp(7)*MyInverse(u +\
 
  v)*MyPow2(-1 + u + v) + 3*(u*Cp(6) + v*Cp(7))*MyInverse(u + v)*MyPow2(-1 + u + v) + 3*v*(-1 + u + v)*(Cp(8) -\
 
  Cp(9))*MyPow2(u)*MyPow2(MyInverse(-1 + u)) + 3*v*((-1 + u + v)*Cp(8) - v*Cp(9))*MyPow2(u)*MyPow2(MyInverse(-1 + u)) +\
 
  3*(-1 + u + v)*((-1 + u + v)*Cp(8) - v*Cp(9))*MyPow2(u)*MyPow2(MyInverse(-1 + u)) - 3*u*(-1 + u + v)*(Cp(12) -\
 
  Cp(13))*MyPow2(v)*MyPow2(MyInverse(-1 + v)) + 3*u*(-(u*Cp(12)) + (-1 + u + v)*Cp(13))*MyPow2(v)*MyPow2(MyInverse(-1 +\
 
  v)) + 3*(-1 + u + v)*(-(u*Cp(12)) + (-1 + u + v)*Cp(13))*MyPow2(v)*MyPow2(MyInverse(-1 + v)) - 12*u*v*(-1 + u +\
 
  v)*(u*Cp(6) + v*Cp(7))*MyPow2(MyInverse(u + v)) - 3*u*v*Cp(6)*MyPow2(-1 + u + v)*MyPow2(MyInverse(u + v)) -\
 
  3*u*v*Cp(7)*MyPow2(-1 + u + v)*MyPow2(MyInverse(u + v)) - 3*u*(u*Cp(6) + v*Cp(7))*MyPow2(-1 + u +\
 
  v)*MyPow2(MyInverse(u + v)) - 3*v*(u*Cp(6) + v*Cp(7))*MyPow2(-1 + u + v)*MyPow2(MyInverse(u + v)) + 6*u*v*(u*Cp(6) +\
 
  v*Cp(7))*MyPow2(-1 + u + v)*MyPow3(MyInverse(u + v))) ;

      H10 = 4*(-6*u*(-1 + u + v)*Cp(1) + 6*u*(-1 + u + v)*Cp(2) - 6*v*(-1 + u + v)*Cp(5) + 6*v*(-1 + u + v)*Cp(11) +\
 
  6*u*v*Cp(14) - 6*u*v*(-1 + u + v)*(Cp(8) - Cp(9))*MyInverse(-1 + u) - 6*u*v*((-1 + u + v)*Cp(8) -\
 
  v*Cp(9))*MyInverse(-1 + u) - 6*u*(-1 + u + v)*((-1 + u + v)*Cp(8) - v*Cp(9))*MyInverse(-1 + u) - 6*u*v*(-1 + u +\
 
  v)*(-Cp(12) + Cp(13))*MyInverse(-1 + v) - 6*u*v*(-(u*Cp(12)) + (-1 + u + v)*Cp(13))*MyInverse(-1 + v) - 6*v*(-1 + u +\
 
  v)*(-(u*Cp(12)) + (-1 + u + v)*Cp(13))*MyInverse(-1 + v) + 6*u*v*(-1 + u + v)*Cp(6)*MyInverse(u + v) + 6*u*v*(-1 + u\
 
  + v)*Cp(7)*MyInverse(u + v) + 6*u*v*(u*Cp(6) + v*Cp(7))*MyInverse(u + v) + 6*u*(-1 + u + v)*(u*Cp(6) +\
 
  v*Cp(7))*MyInverse(u + v) + 6*v*(-1 + u + v)*(u*Cp(6) + v*Cp(7))*MyInverse(u + v) + 3*Cp(2)*MyPow2(u) -\
 
  3*Cp(3)*MyPow2(u) + 3*Cp(10)*MyPow2(u) - 3*v*Cp(8)*MyInverse(-1 + u)*MyPow2(u) - 3*(-1 + u + v)*Cp(8)*MyInverse(-1 +\
 
  u)*MyPow2(u) - 3*v*(Cp(8) - Cp(9))*MyInverse(-1 + u)*MyPow2(u) - 3*((-1 + u + v)*Cp(8) - v*Cp(9))*MyInverse(-1 +\
 
  u)*MyPow2(u) + 3*Cp(11)*MyPow2(v) - 3*Cp(15)*MyPow2(v) + 3*Cp(16)*MyPow2(v) - 3*u*Cp(13)*MyInverse(-1 + v)*MyPow2(v)\
 
  - 3*(-1 + u + v)*Cp(13)*MyInverse(-1 + v)*MyPow2(v) - 3*u*(-Cp(12) + Cp(13))*MyInverse(-1 + v)*MyPow2(v) -\
 
  3*(-(u*Cp(12)) + (-1 + u + v)*Cp(13))*MyInverse(-1 + v)*MyPow2(v) + 3*Cp(0)*MyPow2(-1 + u + v) - 3*Cp(1)*MyPow2(-1 +\
 
  u + v) - 3*Cp(5)*MyPow2(-1 + u + v) + 3*u*Cp(6)*MyInverse(u + v)*MyPow2(-1 + u + v) + 3*v*Cp(7)*MyInverse(u +\
 
  v)*MyPow2(-1 + u + v) + 3*(u*Cp(6) + v*Cp(7))*MyInverse(u + v)*MyPow2(-1 + u + v) + 3*v*(-1 + u + v)*(Cp(8) -\
 
  Cp(9))*MyPow2(u)*MyPow2(MyInverse(-1 + u)) + 3*v*((-1 + u + v)*Cp(8) - v*Cp(9))*MyPow2(u)*MyPow2(MyInverse(-1 + u)) +\
 
  3*(-1 + u + v)*((-1 + u + v)*Cp(8) - v*Cp(9))*MyPow2(u)*MyPow2(MyInverse(-1 + u)) + 3*u*(-1 + u + v)*(-Cp(12) +\
 
  Cp(13))*MyPow2(v)*MyPow2(MyInverse(-1 + v)) + 3*u*(-(u*Cp(12)) + (-1 + u + v)*Cp(13))*MyPow2(v)*MyPow2(MyInverse(-1 +\
 
  v)) + 3*(-1 + u + v)*(-(u*Cp(12)) + (-1 + u + v)*Cp(13))*MyPow2(v)*MyPow2(MyInverse(-1 + v)) - 12*u*v*(-1 + u +\
 
  v)*(u*Cp(6) + v*Cp(7))*MyPow2(MyInverse(u + v)) - 3*u*v*Cp(6)*MyPow2(-1 + u + v)*MyPow2(MyInverse(u + v)) -\
 
  3*u*v*Cp(7)*MyPow2(-1 + u + v)*MyPow2(MyInverse(u + v)) - 3*u*(u*Cp(6) + v*Cp(7))*MyPow2(-1 + u +\
 
  v)*MyPow2(MyInverse(u + v)) - 3*v*(u*Cp(6) + v*Cp(7))*MyPow2(-1 + u + v)*MyPow2(MyInverse(u + v)) + 6*u*v*(u*Cp(6) +\
 
  v*Cp(7))*MyPow2(-1 + u + v)*MyPow3(MyInverse(u + v))) ;

      H11 = 4*(-6*u*(-1 + u + v)*Cp(1) - 6*v*(-1 + u + v)*Cp(5) + 12*v*(-1 + u + v)*Cp(11) - 6*v*(-1 + u + v)*Cp(15) +\
 
  6*u*v*Cp(16) - 12*u*v*(-1 + u + v)*Cp(13)*MyInverse(-1 + v) - 12*u*v*(-(u*Cp(12)) + (-1 + u + v)*Cp(13))*MyInverse(-1\
 
  + v) - 6*u*(-1 + u + v)*(-(u*Cp(12)) + (-1 + u + v)*Cp(13))*MyInverse(-1 + v) + 12*u*v*(-1 + u + v)*Cp(7)*MyInverse(u\
 
  + v) + 6*u*v*(u*Cp(6) + v*Cp(7))*MyInverse(u + v) + 12*u*(-1 + u + v)*(u*Cp(6) + v*Cp(7))*MyInverse(u + v) +\
 
  3*Cp(2)*MyPow2(u) + 3*Cp(14)*MyPow2(u) - 6*v*(Cp(8) - Cp(9))*MyInverse(-1 + u)*MyPow2(u) - 6*(-1 + u + v)*(Cp(8) -\
 
  Cp(9))*MyInverse(-1 + u)*MyPow2(u) - 6*((-1 + u + v)*Cp(8) - v*Cp(9))*MyInverse(-1 + u)*MyPow2(u) +\
 
  3*Cp(11)*MyPow2(v) - 6*Cp(15)*MyPow2(v) + 3*Cp(17)*MyPow2(v) - 6*u*Cp(13)*MyInverse(-1 + v)*MyPow2(v) +\
 
  3*Cp(0)*MyPow2(-1 + u + v) - 6*Cp(5)*MyPow2(-1 + u + v) + 3*Cp(11)*MyPow2(-1 + u + v) + 6*u*Cp(7)*MyInverse(u +\
 
  v)*MyPow2(-1 + u + v) + 12*u*v*(-1 + u + v)*(-(u*Cp(12)) + (-1 + u + v)*Cp(13))*MyPow2(MyInverse(-1 + v)) + 6*u*(-1 +\
 
  u + v)*Cp(13)*MyPow2(v)*MyPow2(MyInverse(-1 + v)) + 6*u*(-(u*Cp(12)) + (-1 + u +\
 
  v)*Cp(13))*MyPow2(v)*MyPow2(MyInverse(-1 + v)) - 12*u*v*(-1 + u + v)*(u*Cp(6) + v*Cp(7))*MyPow2(MyInverse(u + v)) -\
 
  6*u*v*Cp(7)*MyPow2(-1 + u + v)*MyPow2(MyInverse(u + v)) - 6*u*(u*Cp(6) + v*Cp(7))*MyPow2(-1 + u +\
 
  v)*MyPow2(MyInverse(u + v)) - 6*u*(-1 + u + v)*(-(u*Cp(12)) + (-1 + u + v)*Cp(13))*MyPow2(v)*MyPow3(MyInverse(-1 +\
 
  v)) + 6*u*v*(u*Cp(6) + v*Cp(7))*MyPow2(-1 + u + v)*MyPow3(MyInverse(u + v))) ;
    }
  H[0][0] = H00;
  H[0][1] = H01;
  H[1][0] = H10;
  H[1][1] = H11;
}

#undef NEAR
}
