// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "GregoryPatch.hpp"

#define DEBUG_PRINT 0

namespace percept {
/** given two points on an edge of a face, @param {pi,pj} and their associated
 *    @param normals  {ni,nj}, fit a cubic with normal to the normals and
 *     passing through the points, returning the control points in @param c
 */

void  GregoryPatch::
fitCubic(MDArray& c, const MDArray& pi, const MDArray& pj, const MDArray& ni, const MDArray& nj)
{
  MDArray ci("ci",3), cj("cj",3);

  
      ci(0) = (2*pi(0) + pj(0))/3. ;

      ci(1) = (2*pi(1) + pj(1))/3. ;

      ci(2) = (2*pi(2) + pj(2))/3. ;

      cj(0) = (pi(0) + 2*pj(0))/3. ;

      cj(1) = (pi(1) + 2*pj(1))/3. ;

      cj(2) = (pi(2) + 2*pj(2))/3. ;

      c(0,0) = pi(0) ;

      c(0,1) = pi(1) ;

      c(0,2) = pi(2) ;

      c(1,0) = ci(0) - ni(0)*(ni(0)*(ci(0) - pi(0)) + ni(1)*(ci(1) - pi(1)) + ni(2)*(ci(2) - pi(2))) ;

      c(1,1) = ci(1) - ni(1)*(ni(0)*(ci(0) - pi(0)) + ni(1)*(ci(1) - pi(1)) + ni(2)*(ci(2) - pi(2))) ;

      c(1,2) = ci(2) - ni(2)*(ni(0)*(ci(0) - pi(0)) + ni(1)*(ci(1) - pi(1)) + ni(2)*(ci(2) - pi(2))) ;

      c(2,0) = cj(0) - nj(0)*(nj(0)*(cj(0) - pj(0)) + nj(1)*(cj(1) - pj(1)) + nj(2)*(cj(2) - pj(2))) ;

      c(2,1) = cj(1) - nj(1)*(nj(0)*(cj(0) - pj(0)) + nj(1)*(cj(1) - pj(1)) + nj(2)*(cj(2) - pj(2))) ;

      c(2,2) = cj(2) - nj(2)*(nj(0)*(cj(0) - pj(0)) + nj(1)*(cj(1) - pj(1)) + nj(2)*(cj(2) - pj(2))) ;

      c(3,0) = pj(0) ;

      c(3,1) = pj(1) ;

      c(3,2) = pj(2) ;
}

/** left side of ribbon = p, center (shared edge) = q, right side is r
 * Note: qh will be degree-elevated for quads from input of q array; however,
 * tri elements will pass in qh and q will be calculated from qh FIXME
 */

bool  GregoryPatch::
fitRibbon(MDArray& pin, MDArray& q, MDArray& rin, MDArray& qh, bool pIsTri, bool rIsTri)
{
  MDArray p("ploc",pin.layout()), r("rloc",rin.layout());
  Kokkos::deep_copy(p,pin);
  Kokkos::deep_copy(r,rin);
  double Lam0, Lam1, Mu0, Mu1, DetAAtLhs;
  MDArray pe("pe",4,3), re("re",4,3);
  MDArray m0("m0",3,3), m1("m1",3,3), R0("R0",3), R1("R1",3);
  MDArray b("b",3,3), a("a",3,4), aat("aat",3,3), aatI("aatI",3,3), d("d",3,3), d1("d1",3,3);
  double mpe = 1.e-6, Anorm=0;

  //
  if (pIsTri)
    {
      
      p(0,0) = (3*p(0,0) + q(0,0))/4. ;

      p(0,1) = (3*p(0,1) + q(0,1))/4. ;

      p(0,2) = (3*p(0,2) + q(0,2))/4. ;

      p(3,0) = (3*p(3,0) + q(3,0))/4. ;

      p(3,1) = (3*p(3,1) + q(3,1))/4. ;

      p(3,2) = (3*p(3,2) + q(3,2))/4. ;
    }
  if (rIsTri)
    {
      
      r(0,0) = (q(0,0) + 3*r(0,0))/4. ;

      r(0,1) = (q(0,1) + 3*r(0,1))/4. ;

      r(0,2) = (q(0,2) + 3*r(0,2))/4. ;

      r(3,0) = (q(3,0) + 3*r(3,0))/4. ;

      r(3,1) = (q(3,1) + 3*r(3,1))/4. ;

      r(3,2) = (q(3,2) + 3*r(3,2))/4. ;
    }

  // get the estimate for the ribbon
  
      pe(0,0) = p(0,0) ;

      pe(0,1) = p(0,1) ;

      pe(0,2) = p(0,2) ;

      pe(3,0) = p(3,0) ;

      pe(3,1) = p(3,1) ;

      pe(3,2) = p(3,2) ;

      pe(1,0) = (2*(p(0,0) - q(0,0)))/3. + q(1,0) + (p(3,0) - q(3,0))/3. ;

      pe(1,1) = (2*(p(0,1) - q(0,1)))/3. + q(1,1) + (p(3,1) - q(3,1))/3. ;

      pe(1,2) = (2*(p(0,2) - q(0,2)))/3. + q(1,2) + (p(3,2) - q(3,2))/3. ;

      pe(2,0) = (p(0,0) - q(0,0))/3. + q(2,0) + (2*(p(3,0) - q(3,0)))/3. ;

      pe(2,1) = (p(0,1) - q(0,1))/3. + q(2,1) + (2*(p(3,1) - q(3,1)))/3. ;

      pe(2,2) = (p(0,2) - q(0,2))/3. + q(2,2) + (2*(p(3,2) - q(3,2)))/3. ;

      re(0,0) = r(0,0) ;

      re(0,1) = r(0,1) ;

      re(0,2) = r(0,2) ;

      re(3,0) = r(3,0) ;

      re(3,1) = r(3,1) ;

      re(3,2) = r(3,2) ;

      re(1,0) = q(1,0) + (2*(-q(0,0) + r(0,0)))/3. + (-q(3,0) + r(3,0))/3. ;

      re(1,1) = q(1,1) + (2*(-q(0,1) + r(0,1)))/3. + (-q(3,1) + r(3,1))/3. ;

      re(1,2) = q(1,2) + (2*(-q(0,2) + r(0,2)))/3. + (-q(3,2) + r(3,2))/3. ;

      re(2,0) = q(2,0) + (-q(0,0) + r(0,0))/3. + (2*(-q(3,0) + r(3,0)))/3. ;

      re(2,1) = q(2,1) + (-q(0,1) + r(0,1))/3. + (2*(-q(3,1) + r(3,1)))/3. ;

      re(2,2) = q(2,2) + (-q(0,2) + r(0,2))/3. + (2*(-q(3,2) + r(3,2)))/3. ;

  // solve

  
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

      m0(0,0) = -p(0,0) + r(0,0) ;

      m0(0,1) = qh(0,0) - qh(1,0) ;

      m0(1,0) = -p(0,1) + r(0,1) ;

      m0(1,1) = qh(0,1) - qh(1,1) ;

      m0(2,0) = -p(0,2) + r(0,2) ;

      m0(2,1) = qh(0,2) - qh(1,2) ;

      m1(0,0) = -p(3,0) + r(3,0) ;

      m1(0,1) = qh(3,0) - qh(4,0) ;

      m1(1,0) = -p(3,1) + r(3,1) ;

      m1(1,1) = qh(3,1) - qh(4,1) ;

      m1(2,0) = -p(3,2) + r(3,2) ;

      m1(2,1) = qh(3,2) - qh(4,2) ;

      R0(0) = -p(0,0) + qh(0,0) ;

      R0(1) = -p(0,1) + qh(0,1) ;

      R0(2) = -p(0,2) + qh(0,2) ;

      R1(0) = -p(3,0) + qh(3,0) ;

      R1(1) = -p(3,1) + qh(3,1) ;

      R1(2) = -p(3,2) + qh(3,2) ;

     double dd1, dd2, diffnc, det0, det1, localScaleFactor;
  double bdd1, bdd2, bdiffnc;
  MDArray nc1("nc1",3), nc2("nc2",3);
  MDArray bnc1("bnc1",3), bnc2("bnc2",3);
  
      localScaleFactor = 0.16666666666666666*(MyPow(MyPow2(p(0,0) - qh(0,0)) + MyPow2(p(0,1) - qh(0,1)) + MyPow2(p(0,2)\
 
  - qh(0,2)),0.5) + MyPow(MyPow2(-qh(0,0) + qh(1,0)) + MyPow2(-qh(0,1) + qh(1,1)) + MyPow2(-qh(0,2) + qh(1,2)),0.5) +\
 
  MyPow(MyPow2(p(3,0) - qh(4,0)) + MyPow2(p(3,1) - qh(4,1)) + MyPow2(p(3,2) - qh(4,2)),0.5) + MyPow(MyPow2(qh(3,0) -\
 
  qh(4,0)) + MyPow2(qh(3,1) - qh(4,1)) + MyPow2(qh(3,2) - qh(4,2)),0.5) + MyPow(MyPow2(-qh(0,0) + r(0,0)) +\
 
  MyPow2(-qh(0,1) + r(0,1)) + MyPow2(-qh(0,2) + r(0,2)),0.5) + MyPow(MyPow2(-qh(4,0) + r(3,0)) + MyPow2(-qh(4,1) +\
 
  r(3,1)) + MyPow2(-qh(4,2) + r(3,2)),0.5)) ;

      nc1(0) = p(0,2)*qh(0,1) - p(0,1)*qh(0,2) - p(0,2)*qh(1,1) + qh(0,2)*qh(1,1) + p(0,1)*qh(1,2) - qh(0,1)*qh(1,2) ;

      nc1(1) = -(p(0,2)*qh(0,0)) + p(0,0)*qh(0,2) + p(0,2)*qh(1,0) - qh(0,2)*qh(1,0) - p(0,0)*qh(1,2) + qh(0,0)*qh(1,2)\
 
  ;

      nc1(2) = p(0,1)*qh(0,0) - p(0,0)*qh(0,1) - p(0,1)*qh(1,0) + qh(0,1)*qh(1,0) + p(0,0)*qh(1,1) - qh(0,0)*qh(1,1) ;

      dd1 = MyPow(MyPow2(nc1(0)) + MyPow2(nc1(1)) + MyPow2(nc1(2)),0.5) ;

      nc1(0) = MyInverse(dd1)*nc1(0) ;

      nc1(1) = MyInverse(dd1)*nc1(1) ;

      nc1(2) = MyInverse(dd1)*nc1(2) ;

      nc2(0) = -(qh(0,2)*qh(1,1)) + qh(0,1)*qh(1,2) + qh(0,2)*r(0,1) - qh(1,2)*r(0,1) - qh(0,1)*r(0,2) + qh(1,1)*r(0,2)\
 
  ;

      nc2(1) = qh(0,2)*qh(1,0) - qh(0,0)*qh(1,2) - qh(0,2)*r(0,0) + qh(1,2)*r(0,0) + qh(0,0)*r(0,2) - qh(1,0)*r(0,2) ;

      nc2(2) = -(qh(0,1)*qh(1,0)) + qh(0,0)*qh(1,1) + qh(0,1)*r(0,0) - qh(1,1)*r(0,0) - qh(0,0)*r(0,1) + qh(1,0)*r(0,1)\
 
  ;

      dd2 = MyPow(MyPow2(nc2(0)) + MyPow2(nc2(1)) + MyPow2(nc2(2)),0.5) ;

      nc2(0) = MyInverse(dd2)*nc2(0) ;

      nc2(1) = MyInverse(dd2)*nc2(1) ;

      nc2(2) = MyInverse(dd2)*nc2(2) ;

      diffnc = MyPow2(-nc1(0) + nc2(0)) + MyPow2(-nc1(1) + nc2(1)) + MyPow2(-nc1(2) + nc2(2)) ;

      bnc1(0) = p(3,2)*qh(3,1) - p(3,1)*qh(3,2) - p(3,2)*qh(4,1) + qh(3,2)*qh(4,1) + p(3,1)*qh(4,2) - qh(3,1)*qh(4,2) ;

      bnc1(1) = -(p(3,2)*qh(3,0)) + p(3,0)*qh(3,2) + p(3,2)*qh(4,0) - qh(3,2)*qh(4,0) - p(3,0)*qh(4,2) +\
 
  qh(3,0)*qh(4,2) ;

      bnc1(2) = p(3,1)*qh(3,0) - p(3,0)*qh(3,1) - p(3,1)*qh(4,0) + qh(3,1)*qh(4,0) + p(3,0)*qh(4,1) - qh(3,0)*qh(4,1) ;

      bdd1 = MyPow(MyPow2(bnc1(0)) + MyPow2(bnc1(1)) + MyPow2(bnc1(2)),0.5) ;

      bnc1(0) = bnc1(0)*MyInverse(bdd1) ;

      bnc1(1) = bnc1(1)*MyInverse(bdd1) ;

      bnc1(2) = bnc1(2)*MyInverse(bdd1) ;

      bnc2(0) = -(qh(3,2)*qh(4,1)) + qh(3,1)*qh(4,2) + qh(3,2)*r(3,1) - qh(4,2)*r(3,1) - qh(3,1)*r(3,2) +\
 
  qh(4,1)*r(3,2) ;

      bnc2(1) = qh(3,2)*qh(4,0) - qh(3,0)*qh(4,2) - qh(3,2)*r(3,0) + qh(4,2)*r(3,0) + qh(3,0)*r(3,2) - qh(4,0)*r(3,2) ;

      bnc2(2) = -(qh(3,1)*qh(4,0)) + qh(3,0)*qh(4,1) + qh(3,1)*r(3,0) - qh(4,1)*r(3,0) - qh(3,0)*r(3,1) +\
 
  qh(4,0)*r(3,1) ;

      bdd2 = MyPow(MyPow2(bnc2(0)) + MyPow2(bnc2(1)) + MyPow2(bnc2(2)),0.5) ;

      bnc2(0) = bnc2(0)*MyInverse(bdd2) ;

      bnc2(1) = bnc2(1)*MyInverse(bdd2) ;

      bnc2(2) = bnc2(2)*MyInverse(bdd2) ;

      bdiffnc = MyPow2(-bnc1(0) + bnc2(0)) + MyPow2(-bnc1(1) + bnc2(1)) + MyPow2(-bnc1(2) + bnc2(2)) ;

      det0 = -2*m0(0,0)*m0(0,1)*m0(1,0)*m0(1,1) - 2*m0(0,0)*m0(0,1)*m0(2,0)*m0(2,1) - 2*m0(1,0)*m0(1,1)*m0(2,0)*m0(2,1)\
 
  + MyPow2(m0(0,1))*MyPow2(m0(1,0)) + MyPow2(m0(0,0))*MyPow2(m0(1,1)) + MyPow2(m0(0,1))*MyPow2(m0(2,0)) +\
 
  MyPow2(m0(1,1))*MyPow2(m0(2,0)) + MyPow2(m0(0,0))*MyPow2(m0(2,1)) + MyPow2(m0(1,0))*MyPow2(m0(2,1)) ;

      det1 = -2*m1(0,0)*m1(0,1)*m1(1,0)*m1(1,1) - 2*m1(0,0)*m1(0,1)*m1(2,0)*m1(2,1) - 2*m1(1,0)*m1(1,1)*m1(2,0)*m1(2,1)\
 
  + MyPow2(m1(0,1))*MyPow2(m1(1,0)) + MyPow2(m1(0,0))*MyPow2(m1(1,1)) + MyPow2(m1(0,1))*MyPow2(m1(2,0)) +\
 
  MyPow2(m1(1,1))*MyPow2(m1(2,0)) + MyPow2(m1(0,0))*MyPow2(m1(2,1)) + MyPow2(m1(1,0))*MyPow2(m1(2,1)) ;
     if (DEBUG_PRINT) std::cout << "\n\n ======================================================================= \n"
                                << "\npIsTri= " << pIsTri << " rIsTri= " << rIsTri << std::endl;
    if (DEBUG_PRINT) std::cout << "m0=\n" << printContainer(m0) << " m1=\n" << printContainer(m1) << std::endl;

    if (DEBUG_PRINT) std::cout << "det0 = " << det0 << " det1= " << det1 << std::endl;

    double scf1 = 0.0; //1.e-20*localScaleFactor;
    double scf2 = 1.e-4*localScaleFactor;
    VERIFY_OP_ON(localScaleFactor, >, 0.0, "localScaleFactor");
    VERIFY_OP_ON(det0, >, scf1, "det0");
    VERIFY_OP_ON(det1, >, scf1, "det1");
    VERIFY_OP_ON(dd1, >, scf1, "dd1");
    VERIFY_OP_ON(dd2, >, scf1, "dd2");
    if (0) VERIFY_OP_ON(std::fabs(diffnc), <, scf2, "diffnc");
    VERIFY_OP_ON(bdd1, >, scf1, "bdd1");
    VERIFY_OP_ON(bdd2, >, scf1, "bdd2");
    if (0) VERIFY_OP_ON(std::fabs(bdiffnc), <, scf2, "diffnc");
 if (DEBUG_PRINT) std::cout << "dd1= " << dd1 << " dd2= " << dd2 << " diffnc= " << diffnc << std::endl;
 if (DEBUG_PRINT) std::cout << "bdd1= " << bdd1 << " bdd2= " << bdd2 << " bdiffnc= " << bdiffnc << std::endl;

 
      Lam0 = MyInverse(-2*m0(0,0)*m0(0,1)*(m0(1,0)*m0(1,1) + m0(2,0)*m0(2,1)) + MyPow2(m0(0,1))*(MyPow2(m0(1,0)) +\
 
  MyPow2(m0(2,0))) + MyPow2(m0(0,0))*(MyPow2(m0(1,1)) + MyPow2(m0(2,1))) + MyPow2(m0(1,1)*m0(2,0) -\
 
  m0(1,0)*m0(2,1)))*(m0(0,0)*(MyPow2(m0(1,1)) + MyPow2(m0(2,1)))*R0(0) + (m0(1,1)*m0(2,0) -\
 
  m0(1,0)*m0(2,1))*(-(m0(2,1)*R0(1)) + m0(1,1)*R0(2)) + MyPow2(m0(0,1))*(m0(1,0)*R0(1) + m0(2,0)*R0(2)) -\
 
  m0(0,1)*(m0(1,0)*m0(1,1)*R0(0) + m0(2,0)*m0(2,1)*R0(0) + m0(0,0)*(m0(1,1)*R0(1) + m0(2,1)*R0(2)))) ;

      Mu0 = MyInverse(-2*m0(0,0)*m0(0,1)*(m0(1,0)*m0(1,1) + m0(2,0)*m0(2,1)) + MyPow2(m0(0,1))*(MyPow2(m0(1,0)) +\
 
  MyPow2(m0(2,0))) + MyPow2(m0(0,0))*(MyPow2(m0(1,1)) + MyPow2(m0(2,1))) + MyPow2(m0(1,1)*m0(2,0) -\
 
  m0(1,0)*m0(2,1)))*(-(m0(0,0)*(m0(1,0)*m0(1,1) + m0(2,0)*m0(2,1))*R0(0)) + (m0(1,1)*m0(2,0) -\
 
  m0(1,0)*m0(2,1))*(m0(2,0)*R0(1) - m0(1,0)*R0(2)) + MyPow2(m0(0,0))*(m0(1,1)*R0(1) + m0(2,1)*R0(2)) +\
 
  m0(0,1)*(MyPow2(m0(1,0))*R0(0) - m0(0,0)*m0(1,0)*R0(1) + m0(2,0)*(m0(2,0)*R0(0) - m0(0,0)*R0(2)))) ;

      Lam1 = MyInverse(-2*m1(0,0)*m1(0,1)*(m1(1,0)*m1(1,1) + m1(2,0)*m1(2,1)) + MyPow2(m1(0,1))*(MyPow2(m1(1,0)) +\
 
  MyPow2(m1(2,0))) + MyPow2(m1(0,0))*(MyPow2(m1(1,1)) + MyPow2(m1(2,1))) + MyPow2(m1(1,1)*m1(2,0) -\
 
  m1(1,0)*m1(2,1)))*(m1(0,0)*(MyPow2(m1(1,1)) + MyPow2(m1(2,1)))*R1(0) + (m1(1,1)*m1(2,0) -\
 
  m1(1,0)*m1(2,1))*(-(m1(2,1)*R1(1)) + m1(1,1)*R1(2)) + MyPow2(m1(0,1))*(m1(1,0)*R1(1) + m1(2,0)*R1(2)) -\
 
  m1(0,1)*(m1(1,0)*m1(1,1)*R1(0) + m1(2,0)*m1(2,1)*R1(0) + m1(0,0)*(m1(1,1)*R1(1) + m1(2,1)*R1(2)))) ;

      Mu1 = MyInverse(-2*m1(0,0)*m1(0,1)*(m1(1,0)*m1(1,1) + m1(2,0)*m1(2,1)) + MyPow2(m1(0,1))*(MyPow2(m1(1,0)) +\
 
  MyPow2(m1(2,0))) + MyPow2(m1(0,0))*(MyPow2(m1(1,1)) + MyPow2(m1(2,1))) + MyPow2(m1(1,1)*m1(2,0) -\
 
  m1(1,0)*m1(2,1)))*(-(m1(0,0)*(m1(1,0)*m1(1,1) + m1(2,0)*m1(2,1))*R1(0)) + (m1(1,1)*m1(2,0) -\
 
  m1(1,0)*m1(2,1))*(m1(2,0)*R1(1) - m1(1,0)*R1(2)) + MyPow2(m1(0,0))*(m1(1,1)*R1(1) + m1(2,1)*R1(2)) +\
 
  m1(0,1)*(MyPow2(m1(1,0))*R1(0) - m1(0,0)*m1(1,0)*R1(1) + m1(2,0)*(m1(2,0)*R1(0) - m1(0,0)*R1(2)))) ;

      b(0,0) = (1 - Lam1)*p(0,0) - (1 - Mu1)*qh(0,0) - Mu1*qh(1,0) - 3*((1 - Mu0)*qh(1,0) + Mu0*qh(2,0)) + Lam1*r(0,0)\
 
  ;

      b(0,1) = (1 - Lam1)*p(0,1) - (1 - Mu1)*qh(0,1) - Mu1*qh(1,1) - 3*((1 - Mu0)*qh(1,1) + Mu0*qh(2,1)) + Lam1*r(0,1)\
 
  ;

      b(0,2) = (1 - Lam1)*p(0,2) - (1 - Mu1)*qh(0,2) - Mu1*qh(1,2) - 3*((1 - Mu0)*qh(1,2) + Mu0*qh(2,2)) + Lam1*r(0,2)\
 
  ;

      b(1,0) = (1 - Mu1)*qh(1,0) + (1 - Mu0)*qh(2,0) + Mu1*qh(2,0) + Mu0*qh(3,0) ;

      b(1,1) = (1 - Mu1)*qh(1,1) + (1 - Mu0)*qh(2,1) + Mu1*qh(2,1) + Mu0*qh(3,1) ;

      b(1,2) = (1 - Mu1)*qh(1,2) + (1 - Mu0)*qh(2,2) + Mu1*qh(2,2) + Mu0*qh(3,2) ;

      b(2,0) = (1 - Lam0)*p(3,0) - (1 - Mu0)*qh(3,0) - 3*((1 - Mu1)*qh(2,0) + Mu1*qh(3,0)) - Mu0*qh(4,0) + Lam0*r(3,0)\
 
  ;

      b(2,1) = (1 - Lam0)*p(3,1) - (1 - Mu0)*qh(3,1) - 3*((1 - Mu1)*qh(2,1) + Mu1*qh(3,1)) - Mu0*qh(4,1) + Lam0*r(3,1)\
 
  ;

      b(2,2) = (1 - Lam0)*p(3,2) - (1 - Mu0)*qh(3,2) - 3*((1 - Mu1)*qh(2,2) + Mu1*qh(3,2)) - Mu0*qh(4,2) + Lam0*r(3,2)\
 
  ;

      a(0,0) = -3*(1 - Lam0) ;

      a(0,1) = -3*Lam0 ;

      a(0,2) = 0 ;

      a(0,3) = 0 ;

      a(1,0) = 1 - Lam1 ;

      a(1,1) = Lam1 ;

      a(1,2) = 1 - Lam0 ;

      a(1,3) = Lam0 ;

      a(2,0) = 0 ;

      a(2,1) = 0 ;

      a(2,2) = -3*(1 - Lam1) ;

      a(2,3) = -3*Lam1 ;

      Anorm = MyPow(MyPow2(a(0,0)) + MyPow2(a(0,1)) + MyPow2(a(0,2)) + MyPow2(a(0,3)) + MyPow2(a(1,0)) + MyPow2(a(1,1))\
 
  + MyPow2(a(1,2)) + MyPow2(a(1,3)) + MyPow2(a(2,0)) + MyPow2(a(2,1)) + MyPow2(a(2,2)) + MyPow2(a(2,3)),0.5) ;

      aat(0,0) = Anorm*mpe + MyPow2(a(0,0)) + MyPow2(a(0,1)) + MyPow2(a(0,2)) + MyPow2(a(0,3)) ;

      aat(0,1) = a(0,0)*a(1,0) + a(0,1)*a(1,1) + a(0,2)*a(1,2) + a(0,3)*a(1,3) ;

      aat(0,2) = a(0,0)*a(2,0) + a(0,1)*a(2,1) + a(0,2)*a(2,2) + a(0,3)*a(2,3) ;

      aat(1,0) = a(0,0)*a(1,0) + a(0,1)*a(1,1) + a(0,2)*a(1,2) + a(0,3)*a(1,3) ;

      aat(1,1) = Anorm*mpe + MyPow2(a(1,0)) + MyPow2(a(1,1)) + MyPow2(a(1,2)) + MyPow2(a(1,3)) ;

      aat(1,2) = a(1,0)*a(2,0) + a(1,1)*a(2,1) + a(1,2)*a(2,2) + a(1,3)*a(2,3) ;

      aat(2,0) = a(0,0)*a(2,0) + a(0,1)*a(2,1) + a(0,2)*a(2,2) + a(0,3)*a(2,3) ;

      aat(2,1) = a(1,0)*a(2,0) + a(1,1)*a(2,1) + a(1,2)*a(2,2) + a(1,3)*a(2,3) ;

      aat(2,2) = Anorm*mpe + MyPow2(a(2,0)) + MyPow2(a(2,1)) + MyPow2(a(2,2)) + MyPow2(a(2,3)) ;

      DetAAtLhs = aat(0,2)*(-(aat(1,1)*aat(2,0)) + aat(1,0)*aat(2,1)) + aat(0,1)*(aat(1,2)*aat(2,0) -\
 
  aat(1,0)*aat(2,2)) + aat(0,0)*(-(aat(1,2)*aat(2,1)) + aat(1,1)*aat(2,2)) ;

      aatI(0,0) = (-(aat(1,2)*aat(2,1)) + aat(1,1)*aat(2,2))*MyInverse(DetAAtLhs) ;

      aatI(0,1) = (aat(0,2)*aat(2,1) - aat(0,1)*aat(2,2))*MyInverse(DetAAtLhs) ;

      aatI(0,2) = (-(aat(0,2)*aat(1,1)) + aat(0,1)*aat(1,2))*MyInverse(DetAAtLhs) ;

      aatI(1,0) = (aat(1,2)*aat(2,0) - aat(1,0)*aat(2,2))*MyInverse(DetAAtLhs) ;

      aatI(1,1) = (-(aat(0,2)*aat(2,0)) + aat(0,0)*aat(2,2))*MyInverse(DetAAtLhs) ;

      aatI(1,2) = (aat(0,2)*aat(1,0) - aat(0,0)*aat(1,2))*MyInverse(DetAAtLhs) ;

      aatI(2,0) = (-(aat(1,1)*aat(2,0)) + aat(1,0)*aat(2,1))*MyInverse(DetAAtLhs) ;

      aatI(2,1) = (aat(0,1)*aat(2,0) - aat(0,0)*aat(2,1))*MyInverse(DetAAtLhs) ;

      aatI(2,2) = (-(aat(0,1)*aat(1,0)) + aat(0,0)*aat(1,1))*MyInverse(DetAAtLhs) ;

      d1(0,0) = -((-1 + Lam1)*p(0,0)) - a(0,0)*pe(1,0) - a(0,2)*pe(2,0) - qh(0,0) + Mu1*qh(0,0) - 3*qh(1,0) +\
 
  3*Mu0*qh(1,0) - Mu1*qh(1,0) - 3*Mu0*qh(2,0) + Lam1*r(0,0) - a(0,1)*re(1,0) - a(0,3)*re(2,0) ;

      d1(0,1) = -((-1 + Lam1)*p(0,1)) - a(0,0)*pe(1,1) - a(0,2)*pe(2,1) - qh(0,1) + Mu1*qh(0,1) - 3*qh(1,1) +\
 
  3*Mu0*qh(1,1) - Mu1*qh(1,1) - 3*Mu0*qh(2,1) + Lam1*r(0,1) - a(0,1)*re(1,1) - a(0,3)*re(2,1) ;

      d1(0,2) = -((-1 + Lam1)*p(0,2)) - a(0,0)*pe(1,2) - a(0,2)*pe(2,2) - qh(0,2) + Mu1*qh(0,2) - 3*qh(1,2) +\
 
  3*Mu0*qh(1,2) - Mu1*qh(1,2) - 3*Mu0*qh(2,2) + Lam1*r(0,2) - a(0,1)*re(1,2) - a(0,3)*re(2,2) ;

      d1(1,0) = -(a(1,0)*pe(1,0)) - a(1,2)*pe(2,0) + qh(1,0) - Mu1*qh(1,0) + qh(2,0) - Mu0*qh(2,0) + Mu1*qh(2,0) +\
 
  Mu0*qh(3,0) - a(1,1)*re(1,0) - a(1,3)*re(2,0) ;

      d1(1,1) = -(a(1,0)*pe(1,1)) - a(1,2)*pe(2,1) + qh(1,1) - Mu1*qh(1,1) + qh(2,1) - Mu0*qh(2,1) + Mu1*qh(2,1) +\
 
  Mu0*qh(3,1) - a(1,1)*re(1,1) - a(1,3)*re(2,1) ;

      d1(1,2) = -(a(1,0)*pe(1,2)) - a(1,2)*pe(2,2) + qh(1,2) - Mu1*qh(1,2) + qh(2,2) - Mu0*qh(2,2) + Mu1*qh(2,2) +\
 
  Mu0*qh(3,2) - a(1,1)*re(1,2) - a(1,3)*re(2,2) ;

      d1(2,0) = -((-1 + Lam0)*p(3,0)) - a(2,0)*pe(1,0) - a(2,2)*pe(2,0) - 3*qh(2,0) + 3*Mu1*qh(2,0) - qh(3,0) +\
 
  Mu0*qh(3,0) - 3*Mu1*qh(3,0) - Mu0*qh(4,0) + Lam0*r(3,0) - a(2,1)*re(1,0) - a(2,3)*re(2,0) ;

      d1(2,1) = -((-1 + Lam0)*p(3,1)) - a(2,0)*pe(1,1) - a(2,2)*pe(2,1) - 3*qh(2,1) + 3*Mu1*qh(2,1) - qh(3,1) +\
 
  Mu0*qh(3,1) - 3*Mu1*qh(3,1) - Mu0*qh(4,1) + Lam0*r(3,1) - a(2,1)*re(1,1) - a(2,3)*re(2,1) ;

      d1(2,2) = -((-1 + Lam0)*p(3,2)) - a(2,0)*pe(1,2) - a(2,2)*pe(2,2) - 3*qh(2,2) + 3*Mu1*qh(2,2) - qh(3,2) +\
 
  Mu0*qh(3,2) - 3*Mu1*qh(3,2) - Mu0*qh(4,2) + Lam0*r(3,2) - a(2,1)*re(1,2) - a(2,3)*re(2,2) ;

      d(0,0) = aatI(0,0)*d1(0,0) + aatI(0,1)*d1(1,0) + aatI(0,2)*d1(2,0) ;

      d(0,1) = aatI(0,0)*d1(0,1) + aatI(0,1)*d1(1,1) + aatI(0,2)*d1(2,1) ;

      d(0,2) = aatI(0,0)*d1(0,2) + aatI(0,1)*d1(1,2) + aatI(0,2)*d1(2,2) ;

      d(1,0) = aatI(1,0)*d1(0,0) + aatI(1,1)*d1(1,0) + aatI(1,2)*d1(2,0) ;

      d(1,1) = aatI(1,0)*d1(0,1) + aatI(1,1)*d1(1,1) + aatI(1,2)*d1(2,1) ;

      d(1,2) = aatI(1,0)*d1(0,2) + aatI(1,1)*d1(1,2) + aatI(1,2)*d1(2,2) ;

      d(2,0) = aatI(2,0)*d1(0,0) + aatI(2,1)*d1(1,0) + aatI(2,2)*d1(2,0) ;

      d(2,1) = aatI(2,0)*d1(0,1) + aatI(2,1)*d1(1,1) + aatI(2,2)*d1(2,1) ;

      d(2,2) = aatI(2,0)*d1(0,2) + aatI(2,1)*d1(1,2) + aatI(2,2)*d1(2,2) ;

      p(1,0) = a(0,0)*d(0,0) + a(1,0)*d(1,0) + a(2,0)*d(2,0) + pe(1,0) ;

      p(1,1) = a(0,0)*d(0,1) + a(1,0)*d(1,1) + a(2,0)*d(2,1) + pe(1,1) ;

      p(1,2) = a(0,0)*d(0,2) + a(1,0)*d(1,2) + a(2,0)*d(2,2) + pe(1,2) ;

      r(1,0) = a(0,1)*d(0,0) + a(1,1)*d(1,0) + a(2,1)*d(2,0) + re(1,0) ;

      r(1,1) = a(0,1)*d(0,1) + a(1,1)*d(1,1) + a(2,1)*d(2,1) + re(1,1) ;

      r(1,2) = a(0,1)*d(0,2) + a(1,1)*d(1,2) + a(2,1)*d(2,2) + re(1,2) ;

      p(2,0) = a(0,2)*d(0,0) + a(1,2)*d(1,0) + a(2,2)*d(2,0) + pe(2,0) ;

      p(2,1) = a(0,2)*d(0,1) + a(1,2)*d(1,1) + a(2,2)*d(2,1) + pe(2,1) ;

      p(2,2) = a(0,2)*d(0,2) + a(1,2)*d(1,2) + a(2,2)*d(2,2) + pe(2,2) ;

      r(2,0) = a(0,3)*d(0,0) + a(1,3)*d(1,0) + a(2,3)*d(2,0) + re(2,0) ;

      r(2,1) = a(0,3)*d(0,1) + a(1,3)*d(1,1) + a(2,3)*d(2,1) + re(2,1) ;

      r(2,2) = a(0,3)*d(0,2) + a(1,3)*d(1,2) + a(2,3)*d(2,2) + re(2,2) ;

 if (DEBUG_PRINT) std::cout << "Lam0= " << Lam0 << " Lam1= " << Lam1 << " Mu0= " << Mu0 << " Mu1= " << Mu1 << std::endl;
 if (DEBUG_PRINT) std::cout << "\nMatrix= " << printForMathematica(a) << "\nrhs= " << printForMathematica(b) << std::endl;
 if (DEBUG_PRINT) std::cout << "\nDetAAtLhs= " << DetAAtLhs << std::endl;
 if (DEBUG_PRINT) std::cout << "\npIsTri= " << pIsTri << " rIsTri= " << rIsTri
                            << "\n\npin = " << printForMathematica(pin, false)
                            << "\nq = " << printForMathematica(q, false)
                            << "\nrin = " << printForMathematica(rin, false)
                            // << "\npt = " << printForMathematica(pt, false)
                            // << "\nrt = " << printForMathematica(rt, false)
                            << "\n\npe = " << printForMathematica(pe, false)
                            << "\nre = " << printForMathematica(re, false)
                            << "\n\nqh = " << printForMathematica(qh, false)
                            << "\n\nm0 = " << printForMathematica(m0, false)
                            << "\nR0 = " << printForMathematica(R0, false)
                            << "\n\nm1 = " << printForMathematica(m1, false)
                            << "\nR1 = " << printForMathematica(R1, false)
                            << "\n\npout = " << printForMathematica(p, false)
                            << "\nq = " << printForMathematica(q, false)
                            << "\nrout = " << printForMathematica(r, false)
                            << std::endl;

   if (DEBUG_PRINT) std::cout << "\n\n ======================================================================= \n" << std::endl;

  pin = p;
  rin = r;

  return false;
}

void  GregoryPatch::
fitRibbonNoNeighbor(MDArray& pin, MDArray& q, MDArray& qh, bool pIsTri)
{
  MDArray p("ploc",pin.layout());
  Kokkos::deep_copy(p,pin);
  MDArray pe("pe",4,3);

  //
  if (pIsTri)
    {
      
      p(0,0) = (3*p(0,0) + q(0,0))/4. ;

      p(0,1) = (3*p(0,1) + q(0,1))/4. ;

      p(0,2) = (3*p(0,2) + q(0,2))/4. ;

      p(3,0) = (3*p(3,0) + q(3,0))/4. ;

      p(3,1) = (3*p(3,1) + q(3,1))/4. ;

      p(3,2) = (3*p(3,2) + q(3,2))/4. ;
    }

  // get the estimate for the ribbon
  
      pe(0,0) = p(0,0) ;

      pe(0,1) = p(0,1) ;

      pe(0,2) = p(0,2) ;

      pe(3,0) = p(3,0) ;

      pe(3,1) = p(3,1) ;

      pe(3,2) = p(3,2) ;

      pe(1,0) = (2*(p(0,0) - q(0,0)))/3. + q(1,0) + (p(3,0) - q(3,0))/3. ;

      pe(1,1) = (2*(p(0,1) - q(0,1)))/3. + q(1,1) + (p(3,1) - q(3,1))/3. ;

      pe(1,2) = (2*(p(0,2) - q(0,2)))/3. + q(1,2) + (p(3,2) - q(3,2))/3. ;

      pe(2,0) = (p(0,0) - q(0,0))/3. + q(2,0) + (2*(p(3,0) - q(3,0)))/3. ;

      pe(2,1) = (p(0,1) - q(0,1))/3. + q(2,1) + (2*(p(3,1) - q(3,1)))/3. ;

      pe(2,2) = (p(0,2) - q(0,2))/3. + q(2,2) + (2*(p(3,2) - q(3,2)))/3. ;

  for (int ip = 1; ip < 3; ++ip)
    {
      for (int jc = 0; jc < 3; ++jc)
        {
          p(ip,jc) = pe(ip,jc);
        }
    }
  pin = p;
}

}
