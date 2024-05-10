/*
// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

///////////////////////////////////////////////////////////////////////////////
//
// File: Intrepid_PolylibDef.hpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description:
// This file is redistributed with the Intrepid package. It should be used
// in accordance with the above MIT license, at the request of the authors.
// This file is NOT covered by the usual Intrepid/Trilinos LGPL license.
//
// Origin: Nektar++ library, http://www.nektar.info, downloaded on
//         March 10, 2009.
//
///////////////////////////////////////////////////////////////////////////////


/** \file   Intrepid_PolylibDef.hpp
    \brief  Definition file for a set of functions providing orthogonal polynomial
            polynomial calculus and interpolation.
    \author Created by Spencer Sherwin, Aeronautics, Imperial College London,
            modified and redistributed by D. Ridzal.
*/

namespace Intrepid {

/// Maximum number of iterations in polynomial defalation routine Jacobz
#define INTREPID_POLYLIB_STOP 50

/// Define whether to use polynomial deflation (1) or tridiagonal solver (0).
#define INTREPID_POLYLIB_POLYNOMIAL_DEFLATION 0

#ifdef INTREPID_POLYLIB_POLYNOMIAL_DEFLATION
    /// zero determination using Newton iteration with polynomial deflation
#define jacobz(n,z,alpha,beta) Jacobz(n,z,alpha,beta)
#else
    /// zero determination using eigenvalues of tridiagonal matrix
#define jacobz(n,z,alpha,beta) JacZeros(n,z,alpha,beta)
#endif


template <class Scalar>
void IntrepidPolylib::zwgj (Scalar *z, Scalar *w, const int np, const Scalar alpha, const Scalar beta){
  int i;
  Scalar fac, one = 1.0, two = 2.0, apb = alpha + beta;

  IntrepidPolylib::jacobz (np,z,alpha,beta);
  IntrepidPolylib::jacobd (np,z,w,np,alpha,beta);

  fac  = std::pow(two,apb + one)*IntrepidPolylib::gammaF(alpha + np + one)*IntrepidPolylib::gammaF(beta + np + one);
  fac /= IntrepidPolylib::gammaF((Scalar)(np + one))*IntrepidPolylib::gammaF(apb + np + one);

  for(i = 0; i < np; ++i) w[i] = fac/(w[i]*w[i]*(one-z[i]*z[i]));

  return;
}


template <class Scalar>
void IntrepidPolylib::zwgrjm(Scalar *z, Scalar *w, const int np, const Scalar alpha, const Scalar beta){

  if(np == 1){
    z[0] = 0.0;
    w[0] = 2.0;
  }
  else{
    int i;
    Scalar fac, one = 1.0, two = 2.0, apb = alpha + beta;

    z[0] = -one;
    IntrepidPolylib::jacobz  (np-1,z+1,alpha,beta+1);
    IntrepidPolylib::jacobfd (np,z,w,(Scalar*)0,np-1,alpha,beta);

    fac  = std::pow(two,apb)*IntrepidPolylib::gammaF(alpha + np)*IntrepidPolylib::gammaF(beta + np);
    fac /= IntrepidPolylib::gammaF((Scalar)np)*(beta + np)*IntrepidPolylib::gammaF(apb + np + 1);

    for(i = 0; i < np; ++i) w[i] = fac*(1-z[i])/(w[i]*w[i]);
    w[0] *= (beta + one);
  }

  return;
}


template <class Scalar>
void IntrepidPolylib::zwgrjp(Scalar *z, Scalar *w, const int np, const Scalar alpha, const Scalar beta){

  if(np == 1){
    z[0] = 0.0;
    w[0] = 2.0;
  }
  else{
    int i;
    Scalar fac, one = 1.0, two = 2.0, apb = alpha + beta;

    IntrepidPolylib::jacobz  (np-1,z,alpha+1,beta);
    z[np-1] = one;
    IntrepidPolylib::jacobfd (np,z,w,(Scalar*)0,np-1,alpha,beta);

    fac  = std::pow(two,apb)*IntrepidPolylib::gammaF(alpha + np)*IntrepidPolylib::gammaF(beta + np);
    fac /= IntrepidPolylib::gammaF((Scalar)np)*(alpha + np)*IntrepidPolylib::gammaF(apb + np + 1);

    for(i = 0; i < np; ++i) w[i] = fac*(1+z[i])/(w[i]*w[i]);
    w[np-1] *= (alpha + one);
  }

  return;
}


template <class Scalar>
void IntrepidPolylib::zwglj(Scalar *z, Scalar *w, const int np, const Scalar alpha, const Scalar beta){

  if( np == 1 ){
    z[0] = 0.0;
    w[0] = 2.0;
  }
  else{
    int i;
    Scalar   fac, one = 1.0, apb = alpha + beta, two = 2.0;

    z[0]    = -one;
    z[np-1] =  one;
    IntrepidPolylib::jacobz  (np-2,z + 1,alpha + one,beta + one);
    IntrepidPolylib::jacobfd (np,z,w,(Scalar*)0,np-1,alpha,beta);

    fac  = std::pow(two,apb + 1)*IntrepidPolylib::gammaF(alpha + np)*IntrepidPolylib::gammaF(beta + np);
    fac /= (np-1)*IntrepidPolylib::gammaF((Scalar)np)*IntrepidPolylib::gammaF(alpha + beta + np + one);

    for(i = 0; i < np; ++i) w[i] = fac/(w[i]*w[i]);
    w[0]    *= (beta  + one);
    w[np-1] *= (alpha + one);
  }

  return;
}


template <class Scalar>
void IntrepidPolylib::Dgj(Scalar *D,  const Scalar *z, const int np, const Scalar alpha, const Scalar beta)
{

    Scalar one = 1.0, two = 2.0;

    if (np <= 0){
        D[0] = 0.0;
    }
    else{
        int i,j; 
        Scalar *pd;

        pd = (Scalar *)malloc(np*sizeof(Scalar));
        IntrepidPolylib::jacobd(np,z,pd,np,alpha,beta);

        for (i = 0; i < np; i++){
            for (j = 0; j < np; j++){

                if (i != j)
                    //D[i*np+j] = pd[j]/(pd[i]*(z[j]-z[i])); <--- This is either a bug, or the derivative matrix is not defined consistently.
                    D[i*np+j] = pd[i]/(pd[j]*(z[i]-z[j]));
                else
                    D[i*np+j] = (alpha - beta + (alpha + beta + two)*z[j])/
                    (two*(one - z[j]*z[j]));
            }
        }
        free(pd);
    }
    return;
}


template <class Scalar>
void IntrepidPolylib::Dgrjm(Scalar *D, const Scalar *z, const int np, const Scalar alpha, const Scalar beta)
{

    if (np <= 0){
        D[0] = 0.0;
    }
    else{
        int i, j; 
        Scalar   one = 1.0, two = 2.0;
        Scalar   *pd;

        pd  = (Scalar *)malloc(np*sizeof(Scalar));

        pd[0] = std::pow(-one,np-1)*IntrepidPolylib::gammaF((Scalar)np+beta+one);
        pd[0] /= IntrepidPolylib::gammaF((Scalar)np)*IntrepidPolylib::gammaF(beta+two);
        IntrepidPolylib::jacobd(np-1,z+1,pd+1,np-1,alpha,beta+1);
        for(i = 1; i < np; ++i) pd[i] *= (1+z[i]);

        for (i = 0; i < np; i++) {
            for (j = 0; j < np; j++){
                if (i != j)
                    //D[i*np+j] = pd[j]/(pd[i]*(z[j]-z[i])); <--- This is either a bug, or the derivative matrix is not defined consistently.
                    D[i*np+j] = pd[i]/(pd[j]*(z[i]-z[j]));
                else {
                    if(j == 0)
                        D[i*np+j] = -(np + alpha + beta + one)*(np - one)/
                        (two*(beta + two));
                    else
                        D[i*np+j] = (alpha - beta + one + (alpha + beta + one)*z[j])/
                        (two*(one - z[j]*z[j]));
                }
            }
        }
        free(pd);
    }
    return;
}


template <class Scalar>
void IntrepidPolylib::Dgrjp(Scalar *D, const Scalar *z, const int np, const Scalar alpha, const Scalar beta)
{

    if (np <= 0){
        D[0] = 0.0;
    }
    else{
        int i, j; 
        Scalar   one = 1.0, two = 2.0;
        Scalar   *pd;

        pd  = (Scalar *)malloc(np*sizeof(Scalar));


        IntrepidPolylib::jacobd(np-1,z,pd,np-1,alpha+1,beta);
        for(i = 0; i < np-1; ++i) pd[i] *= (1-z[i]);
        pd[np-1] = -IntrepidPolylib::gammaF((Scalar)np+alpha+one);
        pd[np-1] /= IntrepidPolylib::gammaF((Scalar)np)*IntrepidPolylib::gammaF(alpha+two);

        for (i = 0; i < np; i++) {
            for (j = 0; j < np; j++){
                if (i != j)
                    //D[i*np+j] = pd[j]/(pd[i]*(z[j]-z[i])); <--- This is either a bug, or the derivative matrix is not defined consistently.
                    D[i*np+j] = pd[i]/(pd[j]*(z[i]-z[j]));
                else {
                    if(j == np-1)
                        D[i*np+j] = (np + alpha + beta + one)*(np - one)/
                        (two*(alpha + two));
                    else
                        D[i*np+j] = (alpha - beta - one + (alpha + beta + one)*z[j])/
                        (two*(one - z[j]*z[j]));
                }
            }
        }
        free(pd);
    }
    return;
}


template <class Scalar>
void IntrepidPolylib::Dglj(Scalar *D, const Scalar *z, const int np, const Scalar alpha, const Scalar beta)
{

    if (np <= 0){
        D[0] = 0.0;
    }
    else{
        int i, j; 
        Scalar   one = 1.0, two = 2.0;
        Scalar   *pd;

        pd  = (Scalar *)malloc(np*sizeof(Scalar));

        pd[0]  = two*std::pow(-one,np)*IntrepidPolylib::gammaF((Scalar)np + beta);
        pd[0] /= IntrepidPolylib::gammaF((Scalar)np - one)*IntrepidPolylib::gammaF(beta + two);
        IntrepidPolylib::jacobd(np-2,z+1,pd+1,np-2,alpha+1,beta+1);
        for(i = 1; i < np-1; ++i) pd[i] *= (one-z[i]*z[i]);
        pd[np-1]  = -two*IntrepidPolylib::gammaF((Scalar)np + alpha);
        pd[np-1] /= IntrepidPolylib::gammaF((Scalar)np - one)*IntrepidPolylib::gammaF(alpha + two);

        for (i = 0; i < np; i++) {
            for (j = 0; j < np; j++){
                if (i != j)
                    //D[i*np+j] = pd[j]/(pd[i]*(z[j]-z[i])); <--- This is either a bug, or the derivative matrix is not defined consistently.
                    D[i*np+j] = pd[i]/(pd[j]*(z[i]-z[j]));
                else {
                    if (j == 0)
                        D[i*np+j] = (alpha - (np-1)*(np + alpha + beta))/(two*(beta+ two));
                    else if (j == np-1)
                        D[i*np+j] =-(beta - (np-1)*(np + alpha + beta))/(two*(alpha+ two));
                    else
                        D[i*np+j] = (alpha - beta + (alpha + beta)*z[j])/
                        (two*(one - z[j]*z[j]));
                }
            }
        }
        free(pd);
    }
    return;
}


template <class Scalar>
Scalar IntrepidPolylib::hgj (const int i, const Scalar z, const Scalar *zgj,
                             const int np, const Scalar alpha, const Scalar beta)
{

    Scalar zi, dz, p, pd, h;

    zi  = *(zgj+i);
    dz  = z - zi;
    if (std::abs(dz) < INTREPID_TOL) return 1.0;

    IntrepidPolylib::jacobd (1, &zi, &pd , np, alpha, beta);
    IntrepidPolylib::jacobfd(1, &z , &p, (Scalar*)0 , np, alpha, beta);
    h = p/(pd*dz);

    return h;
}


template <class Scalar>
Scalar IntrepidPolylib::hgrjm (const int i, const Scalar z, const Scalar *zgrj,
                               const int np, const Scalar alpha, const Scalar beta)
{

    Scalar zi, dz, p, pd, h;

    zi  = *(zgrj+i);
    dz  = z - zi;
    if (std::abs(dz) < INTREPID_TOL) return 1.0;

    IntrepidPolylib::jacobfd (1, &zi, &p , (Scalar*)0, np-1, alpha, beta + 1);
    // need to use this routine in case zi = -1 or 1
    IntrepidPolylib::jacobd  (1, &zi, &pd, np-1, alpha, beta + 1);
    h = (1.0 + zi)*pd + p;
    IntrepidPolylib::jacobfd (1, &z, &p, (Scalar*)0,  np-1, alpha, beta + 1);
    h = (1.0 + z )*p/(h*dz);

    return h;
}


template <class Scalar>
Scalar IntrepidPolylib::hgrjp (const int i, const Scalar z, const Scalar *zgrj,
                               const int np, const Scalar alpha, const Scalar beta)
{

    Scalar zi, dz, p, pd, h;

    zi  = *(zgrj+i);
    dz  = z - zi;
    if (std::abs(dz) < INTREPID_TOL) return 1.0;

    IntrepidPolylib::jacobfd (1, &zi, &p , (Scalar*)0, np-1, alpha+1, beta );
    // need to use this routine in case z = -1 or 1
    IntrepidPolylib::jacobd  (1, &zi, &pd, np-1, alpha+1, beta );
    h = (1.0 - zi)*pd - p;
    IntrepidPolylib::jacobfd (1, &z, &p, (Scalar*)0,  np-1, alpha+1, beta);
    h = (1.0 - z )*p/(h*dz);

    return h;
}


template <class Scalar>
Scalar IntrepidPolylib::hglj (const int i, const Scalar z, const Scalar *zglj,
                              const int np, const Scalar alpha, const Scalar beta)
{
    Scalar one = 1., two = 2.;
    Scalar zi, dz, p, pd, h;

    zi  = *(zglj+i);
    dz  = z - zi;
    if (std::abs(dz) < INTREPID_TOL) return 1.0;

    IntrepidPolylib::jacobfd(1, &zi, &p , (Scalar*)0, np-2, alpha + one, beta + one);
    // need to use this routine in case z = -1 or 1
    IntrepidPolylib::jacobd (1, &zi, &pd, np-2, alpha + one, beta + one);
    h = (one - zi*zi)*pd - two*zi*p;
    IntrepidPolylib::jacobfd(1, &z, &p, (Scalar*)0, np-2, alpha + one, beta + one);
    h = (one - z*z)*p/(h*dz);

    return h;
}


template <class Scalar>
void IntrepidPolylib::Imgj(Scalar *im, const Scalar *zgj, const Scalar *zm, const int nz,
                           const int mz, const Scalar alpha, const Scalar beta){
        Scalar zp;
        int i, j;

        /* old Polylib code */
        for (i = 0; i < mz; ++i) {
            zp = zm[i];
            for (j = 0; j < nz; ++j) {
                im[i*nz+j] = IntrepidPolylib::hgj(j, zp, zgj, nz, alpha, beta);
            }
        }

        /* original Nektar++ code: is this correct???
        for (i = 0; i < nz; ++i) {
            for (j = 0; j < mz; ++j)
            {
                zp = zm[j];
                im [i*mz+j] = IntrepidPolylib::hgj(i, zp, zgj, nz, alpha, beta);
            }
        }
        */

        return;
}


template <class Scalar>
void IntrepidPolylib::Imgrjm(Scalar *im, const Scalar *zgrj, const Scalar *zm, const int nz,
                             const int mz, const Scalar alpha, const Scalar beta)
{
    Scalar zp;
    int i, j;

    for (i = 0; i < mz; i++) {
      zp = zm[i];
      for (j = 0; j < nz; j++) {
        im[i*nz+j] = IntrepidPolylib::hgrjm(j, zp, zgrj, nz, alpha, beta);
      }
    }

    /* original Nektar++ code: is this correct???
    for (i = 0; i < nz; i++) {
        for (j = 0; j < mz; j++)
        {
            zp = zm[j];
            im [i*mz+j] = IntrepidPolylib::hgrjm(i, zp, zgrj, nz, alpha, beta);
        }
    }
    */

    return;
}


template <class Scalar>
void IntrepidPolylib::Imgrjp(Scalar *im, const Scalar *zgrj, const Scalar *zm, const int nz,
                             const int mz, const Scalar alpha, const Scalar beta)
{
        Scalar zp;
        int i, j;

        for (i = 0; i < mz; i++) {
          zp = zm[i];
          for (j = 0; j < nz; j++) {
            im [i*nz+j] = IntrepidPolylib::hgrjp(j, zp, zgrj, nz, alpha, beta);
          }
        }

        /* original Nektar++ code: is this correct?
        for (i = 0; i < nz; i++) {
            for (j = 0; j < mz; j++)
            {
                zp = zm[j];
                im [i*mz+j] = IntrepidPolylib::hgrjp(i, zp, zgrj, nz, alpha, beta);
            }
        }
        */

        return;
}


template <class Scalar>
void IntrepidPolylib::Imglj(Scalar *im, const Scalar *zglj, const Scalar *zm, const int nz,
                            const int mz, const Scalar alpha, const Scalar beta)
{
    Scalar zp;
    int i, j;

    for (i = 0; i < mz; i++) {
      zp = zm[i];
      for (j = 0; j < nz; j++) {
        im[i*nz+j] = IntrepidPolylib::hglj(j, zp, zglj, nz, alpha, beta);
      }
    }

    /* original Nektar++ code: is this correct?
    for (i = 0; i < nz; i++) {
        for (j = 0; j < mz; j++)
        {
            zp = zm[j];
            im[i*mz+j] = IntrepidPolylib::hglj(i, zp, zglj, nz, alpha, beta);
        }
    }
    */

    return;
}


template <class Scalar>
void
IntrepidPolylib::
jacobfd (const int np, const Scalar *z, Scalar *poly_in, Scalar *polyd,
         const int n, const Scalar alpha, const Scalar beta)
{
  const Scalar zero = 0.0, one = 1.0, two = 2.0;

  if (! np) {
    return;
  }

  if (n == 0) {
    if (poly_in) {
      for (int i = 0; i < np; ++i) {
        poly_in[i] = one;
      }
    }
    if (polyd) {
      for (int i = 0; i < np; ++i) {
        polyd[i] = zero;
      }
    }
  }
  else if (n == 1) {
    if (poly_in) {
      for (int i = 0; i < np; ++i) {
        poly_in[i] = 0.5*(alpha - beta + (alpha + beta + two)*z[i]);
      }
    }
    if (polyd) {
      for (int i = 0; i < np; ++i) {
        polyd[i] = 0.5*(alpha + beta + two);
      }
    }
  }
  else {
    Scalar   a1,a2,a3,a4;
    Scalar apb = alpha + beta;
    Scalar   *poly, *polyn1,*polyn2;

    if (poly_in) { // switch for case of no poynomial function return
      polyn1 = (Scalar *)malloc(2*np*sizeof(Scalar));
      polyn2 = polyn1+np;
      poly   = poly_in;
    }
    else{
      polyn1 = (Scalar *)malloc(3*np*sizeof(Scalar));
      polyn2 = polyn1+np;
      poly   = polyn2+np;
    }

    for (int i = 0; i < np; ++i) {
      polyn2[i] = one;
      polyn1[i] = 0.5*(alpha - beta + (alpha + beta + two)*z[i]);
    }

    for (int k = 2; k <= n; ++k) {
      a1 =  two*k*(k + apb)*(two*k + apb - two);
      a2 = (two*k + apb - one)*(alpha*alpha - beta*beta);
      a3 = (two*k + apb - two)*(two*k + apb - one)*(two*k + apb);
      a4 =  two*(k + alpha - one)*(k + beta - one)*(two*k + apb);

      a2 /= a1;
      a3 /= a1;
      a4 /= a1;

      for (int i = 0; i < np; ++i) {
        poly  [i] = (a2 + a3*z[i])*polyn1[i] - a4*polyn2[i];
        polyn2[i] = polyn1[i];
        polyn1[i] = poly  [i];
      }
    }

    if (polyd) {
      a1 = n*(alpha - beta);
      a2 = n*(two*n + alpha + beta);
      a3 = two*(n + alpha)*(n + beta);
      a4 = (two*n + alpha + beta);
      a1 /= a4;  a2 /= a4;   a3 /= a4;

      // note polyn2 points to polyn1 at end of poly iterations
      for (int i = 0; i < np; ++i) {
        polyd[i]  = (a1- a2*z[i])*poly[i] + a3*polyn2[i];
        polyd[i] /= (one - z[i]*z[i]);
      }
    }

    free(polyn1);
  }
}


template <class Scalar>
void IntrepidPolylib::jacobd(const int np, const Scalar *z, Scalar *polyd, const int n,
                             const Scalar alpha, const Scalar beta)
{
  int i;
  Scalar one = 1.0;
  if(n == 0)
    for(i = 0; i < np; ++i) polyd[i] = 0.0;
  else{
    //jacobf(np,z,polyd,n-1,alpha+one,beta+one);
    IntrepidPolylib::jacobfd(np,z,polyd,(Scalar*)0,n-1,alpha+one,beta+one);
    for(i = 0; i < np; ++i) polyd[i] *= 0.5*(alpha + beta + (Scalar)n + one);
  }
  return;
}


template <class Scalar>
void IntrepidPolylib::Jacobz(const int n, Scalar *z, const Scalar alpha, const Scalar beta){
    int i,j,k;
    Scalar   dth = M_PI/(2.0*(Scalar)n);
    Scalar   poly,pder,rlast=0.0;
    Scalar   sum,delr,r;
    Scalar one = 1.0, two = 2.0;

    if(!n)
        return;

    for(k = 0; k < n; ++k){
        r = -std::cos((two*(Scalar)k + one) * dth);
        if(k) r = 0.5*(r + rlast);

        for(j = 1; j < INTREPID_POLYLIB_STOP; ++j){
            IntrepidPolylib::jacobfd(1,&r,&poly, &pder, n, alpha, beta);

            for(i = 0, sum = 0.0; i < k; ++i) sum += one/(r - z[i]);

            delr = -poly / (pder - sum * poly);
            r   += delr;
            if( std::abs(delr) < INTREPID_TOL ) break;
        }
        z[k]  = r;
        rlast = r;
    }
    return;
}


template <class Scalar>
void IntrepidPolylib::JacZeros(const int n, Scalar *a, const Scalar alpha, const Scalar beta){
  int i;
  Scalar apb, apbi,a2b2;
  Scalar *b;

  if(!n)
    return;

  b = (Scalar *) malloc(n*sizeof(Scalar));

  // generate normalised terms
  apb  = alpha + beta;
  apbi = 2.0 + apb;

  b[n-1] = std::pow(2.0,apb+1.0)*IntrepidPolylib::gammaF(alpha+1.0)*IntrepidPolylib::gammaF(beta+1.0)/gammaF(apbi);
  a[0]   = (beta-alpha)/apbi;
  b[0]   = std::sqrt(4.0*(1.0+alpha)*(1.0+beta)/((apbi+1.0)*apbi*apbi));

  a2b2 = beta*beta-alpha*alpha;
  for(i = 1; i < n-1; ++i){
    apbi = 2.0*(i+1) + apb;
    a[i] = a2b2/((apbi-2.0)*apbi);
    b[i] = std::sqrt(4.0*(i+1)*(i+1+alpha)*(i+1+beta)*(i+1+apb)/
                     ((apbi*apbi-1)*apbi*apbi));
  }

  apbi   = 2.0*n + apb;
  //a[n-1] = a2b2/((apbi-2.0)*apbi); // THIS IS A BUG!!!
  if (n>1) a[n-1] = a2b2/((apbi-2.0)*apbi);

  // find eigenvalues
  IntrepidPolylib::TriQL(n, a, b);

  free(b);
  return;
}


template <class Scalar>
void IntrepidPolylib::TriQL(const int n, Scalar *d,Scalar *e) {
  int m,l,iter,i,k;
  Scalar s,r,p,g,f,dd,c,b;

  for (l=0;l<n;l++) {
    iter=0;
    do {
      for (m=l;m<n-1;m++) {
        dd=std::abs(d[m])+std::abs(d[m+1]);
        if (std::abs(e[m])+dd == dd) break;
      }
      if (m != l) {
        if (iter++ == 50){
          TEUCHOS_TEST_FOR_EXCEPTION((1),
                             std::runtime_error,
                             ">>> ERROR (IntrepidPolylib): Too many iterations in TQLI.");
        }
        g=(d[l+1]-d[l])/(2.0*e[l]);
        r=std::sqrt((g*g)+1.0);
        //g=d[m]-d[l]+e[l]/(g+sign(r,g));
        g=d[m]-d[l]+e[l]/(g+((g)<0 ? -std::abs(r) : std::abs(r)));
        s=c=1.0;
        p=0.0;
        for (i=m-1;i>=l;i--) {
          f=s*e[i];
          b=c*e[i];
          if (std::abs(f) >= std::abs(g)) {
            c=g/f;
            r=std::sqrt((c*c)+1.0);
            e[i+1]=f*r;
            c *= (s=1.0/r);
          } else {
            s=f/g;
            r=std::sqrt((s*s)+1.0);
            e[i+1]=g*r;
            s *= (c=1.0/r);
          }
          g=d[i+1]-p;
          r=(d[i]-g)*s+2.0*c*b;
          p=s*r;
          d[i+1]=g+p;
          g=c*r-b;
        }
        d[l]=d[l]-p;
        e[l]=g;
        e[m]=0.0;
      }
    } while (m != l);
  }

  // order eigenvalues
  for(i = 0; i < n-1; ++i){
    k = i;
    p = d[i];
    for(l = i+1; l < n; ++l)
      if (d[l] < p) {
        k = l;
        p = d[l];
      }
    d[k] = d[i];
    d[i] = p;
  }
}


template <class Scalar>
Scalar IntrepidPolylib::gammaF(const Scalar x){
  Scalar gamma = 1.0;

  if     (x == -0.5) gamma = -2.0*std::sqrt(M_PI);
  else if (!x) return gamma;
  else if ((x-(int)x) == 0.5){
    int n = (int) x;
    Scalar tmp = x;

    gamma = std::sqrt(M_PI);
    while(n--){
      tmp   -= 1.0;
      gamma *= tmp;
    }
  }
  else if ((x-(int)x) == 0.0){
    int n = (int) x;
    Scalar tmp = x;

    while(--n){
      tmp   -= 1.0;
      gamma *= tmp;
    }
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION((1),
                       std::invalid_argument,
                       ">>> ERROR (IntrepidPolylib): Argument is not of integer or half order.");
  return gamma;
}

} // end of namespace Intrepid


#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

