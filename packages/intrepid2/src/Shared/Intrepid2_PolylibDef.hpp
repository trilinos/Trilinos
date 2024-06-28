/*
// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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


/** \file   Intrepid2_PolylibDef.hpp
    \brief  Definition file for a set of functions providing orthogonal polynomial
    calculus and interpolation.
    \author Created by Spencer Sherwin, Aeronautics, Imperial College London,
            modified and redistributed by D. Ridzal.
            modified and Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_POLYLIB_DEF_HPP__
#define __INTREPID2_POLYLIB_DEF_HPP__

#if defined(_MSC_VER) || defined(_WIN32) && defined(__ICL)
// M_PI, M_SQRT2, etc. are hidden in MSVC by #ifdef _USE_MATH_DEFINES
  #ifndef _USE_MATH_DEFINES
  #define _USE_MATH_DEFINES
  #endif
  #include <math.h>
#endif

namespace Intrepid2 {

  // -----------------------------------------------------------------------
  // Points and Weights
  // -----------------------------------------------------------------------

  template<>
  template<typename zViewType,
           typename wViewType>
  KOKKOS_INLINE_FUNCTION
  void
  Polylib::Serial::Cubature<POLYTYPE_GAUSS>::
  getValues(      zViewType z,
                  wViewType w,
            const ordinal_type np,
            const double alpha,
            const double beta) {
    const double one = 1.0, two = 2.0, apb = alpha + beta;
    double fac;

    JacobiZeros(z, np, alpha, beta);
    JacobiPolynomialDerivative(np, z, w, np, alpha, beta);
    
    fac  = pow(two, apb + one)*GammaFunction(alpha + np + one)*GammaFunction(beta + np + one);
    fac /= GammaFunction((np + one))*GammaFunction(apb + np + one);

    for (ordinal_type i = 0; i < np; ++i)
      w(i) = fac/(w(i)*w(i)*(one-z(i)*z(i)));
  }

  template<>
  template<typename zViewType,
           typename wViewType>
  KOKKOS_INLINE_FUNCTION
  void
  Polylib::Serial::Cubature<POLYTYPE_GAUSS_RADAU_LEFT>::
  getValues(      zViewType z,
                  wViewType w,
            const ordinal_type np,
            const double alpha,
            const double beta) {
    if (np == 1) {
      z(0) = 0.0;
      w(0) = 2.0;
    } else {
      const double one = 1.0, two = 2.0, apb = alpha + beta;
      double fac;

      z(0) = -one;

      auto z_plus_1 = Kokkos::subview(z, Kokkos::pair<ordinal_type,ordinal_type>(1, z.extent(0)));      
      JacobiZeros(z_plus_1, np-1, alpha, beta+1);

      Kokkos::View<typename zViewType::value_type*,Kokkos::AnonymousSpace,Kokkos::MemoryUnmanaged> null;
      JacobiPolynomial(np, z, w, null, np-1, alpha, beta);

      fac  = pow(two, apb)*GammaFunction(alpha + np)*GammaFunction(beta + np);
      fac /= GammaFunction(np)*(beta + np)*GammaFunction(apb + np + 1);

      for (ordinal_type i = 0; i < np; ++i)
        w(i) = fac*(1-z(i))/(w(i)*w(i));
      w(0) *= (beta + one);
    }
  }

  template<>
  template<typename zViewType,
           typename wViewType>
  KOKKOS_INLINE_FUNCTION
  void
  Polylib::Serial::Cubature<POLYTYPE_GAUSS_RADAU_RIGHT>::
  getValues(      zViewType z,
                  wViewType w,
            const ordinal_type np,
            const double alpha,
            const double beta) {
    if (np == 1) {
      z(0) = 0.0;
      w(0) = 2.0;
    } else {
      const double one = 1.0, two = 2.0, apb = alpha + beta;
      double fac;

      JacobiZeros(z, np-1, alpha+1, beta);
      z(np-1) = one;

      Kokkos::View<typename zViewType::value_type*,Kokkos::AnonymousSpace,Kokkos::MemoryUnmanaged> null;
      JacobiPolynomial(np, z, w, null, np-1, alpha, beta);

      fac  = pow(two,apb)*GammaFunction(alpha + np)*GammaFunction(beta + np);
      fac /= GammaFunction(np)*(alpha + np)*GammaFunction(apb + np + 1);

      for (ordinal_type i = 0; i < np; ++i)
        w(i) = fac*(1+z(i))/(w(i)*w(i));
      w(np-1) *= (alpha + one);
    }
  }

  template<>
  template<typename zViewType,
           typename wViewType>
  KOKKOS_INLINE_FUNCTION
  void
  Polylib::Serial::Cubature<POLYTYPE_GAUSS_LOBATTO>::
  getValues(      zViewType z,
                  wViewType w,
            const ordinal_type np,
            const double alpha,
            const double beta) {
    if ( np == 1 ) {
      z(0) = 0.0;
      w(0) = 2.0;
    } else {
      const double one = 1.0, apb = alpha + beta, two = 2.0;
      double fac;

      z(0)    = -one;
      z(np-1) =  one;

      auto z_plus_1 = Kokkos::subview(z, Kokkos::pair<ordinal_type,ordinal_type>(1, z.extent(0)));      
      JacobiZeros(z_plus_1, np-2, alpha+one, beta+one);

      Kokkos::View<typename zViewType::value_type*,Kokkos::AnonymousSpace,Kokkos::MemoryUnmanaged> null;
      JacobiPolynomial(np, z, w, null, np-1, alpha, beta);

      fac  = pow(two, apb + 1)*GammaFunction(alpha + np)*GammaFunction(beta + np);
      fac /= (np-1)*GammaFunction(np)*GammaFunction(alpha + beta + np + one);

      for (ordinal_type i = 0; i < np; ++i)
        w(i) = fac/(w(i)*w(i));

      w(0)    *= (beta  + one);
      w(np-1) *= (alpha + one);
    }
  }

  // -----------------------------------------------------------------------
  // Derivatives
  // -----------------------------------------------------------------------

  template<>
  template<typename DViewType,
           typename zViewType>
  KOKKOS_INLINE_FUNCTION
  void
  Polylib::Serial::Derivative<POLYTYPE_GAUSS>::
  getValues(      DViewType D,
            const zViewType z,
            const ordinal_type np,
            const double alpha,
            const double beta) {
    if (np <= 0) {
      D(0,0) = 0.0;
    } else {
      const double one = 1.0, two = 2.0;

      typename zViewType::value_type pd_buf[MaxPolylibPoint];
      Kokkos::View<typename zViewType::value_type*,
        Kokkos::AnonymousSpace,Kokkos::MemoryUnmanaged> 
        pd((typename zViewType::pointer_type)&pd_buf[0], MaxPolylibPoint);

      JacobiPolynomialDerivative(np, z, pd, np, alpha, beta);

      for (ordinal_type i = 0; i < np; ++i)
        for (ordinal_type j = 0; j < np; ++j)
          if (i != j)
            //D(i*np+j) = pd(j)/(pd(i)*(z(j)-z(i))); <--- This is either a bug, or the derivative matrix is not defined consistently.
            D(i,j) = pd(i)/(pd(j)*(z(i)-z(j)));
          else
            D(i,j) = (alpha - beta + (alpha + beta + two)*z(j))/
              (two*(one - z(j)*z(j)));
    }
  }

  template<>
  template<typename DViewType,
           typename zViewType>
  KOKKOS_INLINE_FUNCTION
  void
  Polylib::Serial::Derivative<POLYTYPE_GAUSS_RADAU_LEFT>::
  getValues(      DViewType D,
            const zViewType z,
            const ordinal_type np,
            const double alpha,
            const double beta) {
    if (np <= 0) {
      D(0,0) = 0.0;
    } else {
      const double one = 1.0, two = 2.0;

      typename zViewType::value_type pd_buf[MaxPolylibPoint];
      Kokkos::View<typename zViewType::value_type*,
        Kokkos::AnonymousSpace,Kokkos::MemoryUnmanaged>
        pd((typename zViewType::pointer_type)&pd_buf[0], MaxPolylibPoint);
      
      pd(0) = pow(-one,np-1)*GammaFunction(np+beta+one);
      pd(0) /= GammaFunction(np)*GammaFunction(beta+two);

      auto pd_plus_1 = Kokkos::subview(pd, Kokkos::pair<ordinal_type,ordinal_type>(1, pd.extent(0)));
      auto  z_plus_1 = Kokkos::subview( z, Kokkos::pair<ordinal_type,ordinal_type>(1,  z.extent(0)));

      JacobiPolynomialDerivative(np-1, z_plus_1, pd_plus_1, np-1, alpha, beta+1);
      for(ordinal_type i = 1; i < np; ++i)
        pd(i) *= (1+z(i));

      for (ordinal_type i = 0; i < np; ++i) 
        for (ordinal_type j = 0; j < np; ++j) 
          if (i != j)
            D(i,j) = pd(i)/(pd(j)*(z(i)-z(j)));
          else 
            if (j == 0)
              D(i,j) = -(np + alpha + beta + one)*(np - one)/
                (two*(beta + two));
            else
              D(i,j) = (alpha - beta + one + (alpha + beta + one)*z(j))/
                (two*(one - z(j)*z(j)));
    }
  }

  template<>
  template<typename DViewType,
           typename zViewType>
  KOKKOS_INLINE_FUNCTION
  void
  Polylib::Serial::Derivative<POLYTYPE_GAUSS_RADAU_RIGHT>::
  getValues(      DViewType D,
            const zViewType z,
            const ordinal_type np,
            const double alpha,
            const double beta) {
    if (np <= 0) {
      D(0,0) = 0.0;
    } else {
      const double one = 1.0, two = 2.0;

      typename zViewType::value_type pd_buf[MaxPolylibPoint];
      Kokkos::View<typename zViewType::value_type*,
        Kokkos::AnonymousSpace,Kokkos::MemoryUnmanaged>
        pd((typename zViewType::pointer_type)&pd_buf[0], MaxPolylibPoint);

      JacobiPolynomialDerivative(np-1, z, pd, np-1, alpha+1, beta);
      for (ordinal_type i = 0; i < np-1; ++i)
        pd(i) *= (1-z(i));

      pd(np-1) = -GammaFunction(np+alpha+one);
      pd(np-1) /= GammaFunction(np)*GammaFunction(alpha+two);

      for (ordinal_type i = 0; i < np; ++i) 
        for (ordinal_type j = 0; j < np; ++j) 
          if (i != j)
            D(i,j) = pd(i)/(pd(j)*(z(i)-z(j)));
          else 
            if (j == np-1)
              D(i,j) = (np + alpha + beta + one)*(np - one)/
                (two*(alpha + two));
            else
              D(i,j) = (alpha - beta - one + (alpha + beta + one)*z(j))/
                (two*(one - z(j)*z(j)));
    }
  }

  template<>
  template<typename DViewType,
           typename zViewType>
  KOKKOS_INLINE_FUNCTION
  void
  Polylib::Serial::Derivative<POLYTYPE_GAUSS_LOBATTO>::
  getValues(      DViewType D,
            const zViewType z,
            const ordinal_type np,
            const double alpha,
            const double beta) {
    if (np <= 0) {
      D(0,0) = 0.0;
    } else {
      const double one = 1.0, two = 2.0;

      typename zViewType::value_type pd_buf[MaxPolylibPoint];
      Kokkos::View<typename zViewType::value_type*,
        Kokkos::AnonymousSpace,Kokkos::MemoryUnmanaged>
        pd((typename zViewType::pointer_type)&pd_buf[0], MaxPolylibPoint);

      pd(0)  = two*pow(-one,np)*GammaFunction(np + beta);
      pd(0) /= GammaFunction(np - one)*GammaFunction(beta + two);

      auto pd_plus_1 = Kokkos::subview(pd, Kokkos::pair<ordinal_type,ordinal_type>(1, pd.extent(0)));
      auto  z_plus_1 = Kokkos::subview( z, Kokkos::pair<ordinal_type,ordinal_type>(1,  z.extent(0)));

      JacobiPolynomialDerivative(np-2, z_plus_1, pd_plus_1, np-2, alpha+1, beta+1);
      for (ordinal_type i = 1; i < np-1; ++i) 
        pd(i) *= (one-z(i)*z(i));

      pd(np-1)  = -two*GammaFunction(np + alpha);
      pd(np-1) /= GammaFunction(np - one)*GammaFunction(alpha + two);

      for (ordinal_type i = 0; i < np; ++i) 
        for (ordinal_type j = 0; j < np; ++j) 
          if (i != j)
            D(i,j) = pd(i)/(pd(j)*(z(i)-z(j)));
          else 
            if (j == 0)
              D(i,j) = (alpha - (np-1)*(np + alpha + beta))/(two*(beta+ two));
            else if (j == np-1)
              D(i,j) =-(beta - (np-1)*(np + alpha + beta))/(two*(alpha+ two));
            else
              D(i,j) = (alpha - beta + (alpha + beta)*z(j))/
                (two*(one - z(j)*z(j)));
    }
  }
  
  // -----------------------------------------------------------------------
  // Lagrangian Interpolants
  // -----------------------------------------------------------------------

  template<>
  template<typename zViewType>
  KOKKOS_INLINE_FUNCTION
  typename zViewType::value_type
  Polylib::Serial::LagrangianInterpolant<POLYTYPE_GAUSS>::
  getValue(const ordinal_type i,
           const typename zViewType::value_type z,
           const zViewType zgj,
           const ordinal_type np,
           const double alpha,
           const double beta) {
    const double tol = tolerence();

    typedef typename zViewType::value_type value_type;
    typedef typename zViewType::pointer_type pointer_type;

    value_type h, p_buf, pd_buf, zi_buf = zgj(i);
    Kokkos::View<value_type*,Kokkos::AnonymousSpace,Kokkos::MemoryUnmanaged> 
      p((pointer_type)&p_buf, 1), pd((pointer_type)&pd_buf, 1), zi((pointer_type)&zi_buf, 1), 
      zv(const_cast<pointer_type>(&z), 1), null;
    
    const auto dz  = z - zi(0);
    if (Util<value_type>::abs(dz) < tol) 
      return 1.0;

    JacobiPolynomialDerivative(1, zi, pd, np, alpha, beta);
    JacobiPolynomial(1, zv, p, null , np, alpha, beta);

    h = p(0)/(pd(0)*dz);

    return h;
  }

  template<>
  template<typename zViewType>
  KOKKOS_INLINE_FUNCTION
  typename zViewType::value_type
  Polylib::Serial::LagrangianInterpolant<POLYTYPE_GAUSS_RADAU_LEFT>::
  getValue(const ordinal_type i,
           const typename zViewType::value_type z,
           const zViewType zgrj,
           const ordinal_type np,
           const double alpha,
           const double beta) {
    const double tol = tolerence();

    typedef typename zViewType::value_type value_type;
    typedef typename zViewType::pointer_type pointer_type;

    value_type h, p_buf, pd_buf, zi_buf = zgrj(i);
    Kokkos::View<value_type*,Kokkos::AnonymousSpace,Kokkos::MemoryUnmanaged> 
      p((pointer_type)&p_buf, 1), pd((pointer_type)&pd_buf, 1), zi((pointer_type)&zi_buf, 1), 
      zv(const_cast<pointer_type>(&z), 1), null;
    
    const auto dz  = z - zi(0);
    if (Util<value_type>::abs(dz) < tol) 
      return 1.0;
    
    JacobiPolynomial(1, zi, p , null, np-1, alpha, beta + 1);

    // need to use this routine in case zi = -1 or 1
    JacobiPolynomialDerivative(1, zi, pd, np-1, alpha, beta + 1);
    h = (1.0 + zi(0))*pd(0) + p(0);

    JacobiPolynomial(1, zv, p, null,  np-1, alpha, beta + 1);
    h = (1.0 + z)*p(0)/(h*dz);

    return h;
  }


  template<>
  template<typename zViewType>
  KOKKOS_INLINE_FUNCTION
  typename zViewType::value_type
  Polylib::Serial::LagrangianInterpolant<POLYTYPE_GAUSS_RADAU_RIGHT>::
  getValue(const ordinal_type i,
           const typename zViewType::value_type z,
           const zViewType zgrj,
           const ordinal_type np,
           const double alpha,
           const double beta) {
    const double tol = tolerence();

    typedef typename zViewType::value_type value_type;
    typedef typename zViewType::pointer_type pointer_type;

    value_type h, p_buf, pd_buf, zi_buf = zgrj(i);
    Kokkos::View<value_type*,Kokkos::AnonymousSpace,Kokkos::MemoryUnmanaged> 
      p((pointer_type)&p_buf, 1), pd((pointer_type)&pd_buf, 1), zi((pointer_type)&zi_buf, 1), 
      zv(const_cast<pointer_type>(&z), 1), null;

    const auto dz  = z - zi(0);
    if (Util<value_type>::abs(dz) < tol) 
      return 1.0;

    JacobiPolynomial(1, zi, p , null, np-1, alpha+1, beta);

    // need to use this routine in case z = -1 or 1
    JacobiPolynomialDerivative(1, zi, pd, np-1, alpha+1, beta);
    h = (1.0 - zi(0))*pd(0) - p(0);

    JacobiPolynomial (1, zv, p, null,  np-1, alpha+1, beta);
    h = (1.0 - z)*p(0)/(h*dz);

    return h;
  }


  template<>
  template<typename zViewType>
  KOKKOS_INLINE_FUNCTION
  typename zViewType::value_type
  Polylib::Serial::LagrangianInterpolant<POLYTYPE_GAUSS_LOBATTO>::
  getValue(const ordinal_type i,
           const typename zViewType::value_type z,
           const zViewType zglj,
           const ordinal_type np,
           const double alpha,
           const double beta) {
    const double one = 1.0, two = 2.0, tol = tolerence();

    typedef typename zViewType::value_type value_type;
    typedef typename zViewType::pointer_type pointer_type;
    
    value_type h, p_buf, pd_buf, zi_buf = zglj(i);
    Kokkos::View<value_type*,Kokkos::AnonymousSpace,Kokkos::MemoryUnmanaged> 
      p((pointer_type)&p_buf, 1), pd((pointer_type)&pd_buf, 1), zi((pointer_type)&zi_buf, 1), 
      zv(const_cast<pointer_type>(&z), 1), null;

    const auto dz = z - zi(0);
    if (Util<value_type>::abs(dz) < tol) 
      return 1.0;

    JacobiPolynomial(1, zi, p , null, np-2, alpha + one, beta + one);

    // need to use this routine in case z = -1 or 1
    JacobiPolynomialDerivative(1, zi, pd, np-2, alpha + one, beta + one);
    h = (one - zi(0)*zi(0))*pd(0) - two*zi(0)*p(0);

    JacobiPolynomial(1, zv, p, null, np-2, alpha + one, beta + one);
    h = (one - z*z)*p(0)/(h*dz);

    return h;
  }

  // -----------------------------------------------------------------------
  // Interpolation Operator
  // -----------------------------------------------------------------------

  template<EPolyType polyType>
  template<typename imViewType,
           typename zgrjViewType,
           typename zmViewType>
  KOKKOS_INLINE_FUNCTION
  void
  Polylib::Serial::InterpolationOperator<polyType>::
  getMatrix(      imViewType im,
            const zgrjViewType zgrj,
            const zmViewType zm,
            const ordinal_type nz,
            const ordinal_type mz,
            const double alpha,
            const double beta) {
    for (ordinal_type i = 0; i < mz; ++i) {
      const auto zp = zm(i);
      for (ordinal_type j = 0; j < nz; ++j)
        im(i, j) = LagrangianInterpolant<polyType>::getValue(j, zp, zgrj, nz, alpha, beta);
    }
  }

  // -----------------------------------------------------------------------

  template<typename zViewType,
           typename polyiViewType,
           typename polydViewType>
  KOKKOS_INLINE_FUNCTION
  void
  Polylib::Serial::
  JacobiPolynomial(const ordinal_type np,
                   const zViewType z,
                         polyiViewType polyi,
                         polydViewType polyd,
                   const ordinal_type n,
                   const double alpha,
                   const double beta) {
    const double zero = 0.0, one = 1.0, two = 2.0;

    if (! np) {
      return;
    }

    if (n == 0) {
      if (polyi.data()) 
        for (ordinal_type i = 0; i < np; ++i) 
          polyi(i) = one;
      if (polyd.data()) 
        for (ordinal_type i = 0; i < np; ++i) 
          polyd(i) = zero;
    } else if (n == 1) {
      if (polyi.data()) 
        for (ordinal_type i = 0; i < np; ++i) 
          polyi(i) = 0.5*(alpha - beta + (alpha + beta + two)*z(i));
      if (polyd.data()) 
        for (ordinal_type i = 0; i < np; ++i) 
          polyd(i) = 0.5*(alpha + beta + two);
    } else {
      double a1, a2, a3, a4;
      const double apb = alpha + beta;

      typename polyiViewType::value_type
        poly[MaxPolylibPoint]={}, polyn1[MaxPolylibPoint]={}, polyn2[MaxPolylibPoint]={};

      if (polyi.data()) 
        for (ordinal_type i=0;i<np;++i)
          poly[i] = polyi(i);

      for (ordinal_type i = 0; i < np; ++i) {
        polyn2[i] = one;
        polyn1[i] = 0.5*(alpha - beta + (alpha + beta + two)*z(i));
      }

      for (auto k = 2; k <= n; ++k) {
        a1 =  two*k*(k + apb)*(two*k + apb - two);
        a2 = (two*k + apb - one)*(alpha*alpha - beta*beta);
        a3 = (two*k + apb - two)*(two*k + apb - one)*(two*k + apb);
        a4 =  two*(k + alpha - one)*(k + beta - one)*(two*k + apb);

        a2 /= a1;
        a3 /= a1;
        a4 /= a1;

        for (ordinal_type i = 0; i < np; ++i) {
          poly  [i] = (a2 + a3*z(i))*polyn1[i] - a4*polyn2[i];
          polyn2[i] = polyn1[i];
          polyn1[i] = poly  [i];
        }
      }

      if (polyd.data()) {
        a1 = n*(alpha - beta);
        a2 = n*(two*n + alpha + beta);
        a3 = two*(n + alpha)*(n + beta);
        a4 = (two*n + alpha + beta);
        a1 /= a4;
        a2 /= a4;
        a3 /= a4;

        // note polyn2 points to polyn1 at end of poly iterations
        for (ordinal_type i = 0; i < np; ++i) {
          polyd(i)  = (a1- a2*z(i))*poly[i] + a3*polyn2[i];
          polyd(i) /= (one - z(i)*z(i));
        }
      }

      if (polyi.data()) 
        for (ordinal_type i=0;i<np;++i)
          polyi(i) = poly[i];
    }
  }

  template<typename zViewType,
           typename polydViewType>
  KOKKOS_INLINE_FUNCTION
  void
  Polylib::Serial::
  JacobiPolynomialDerivative(const ordinal_type np,
                             const zViewType z,
                                   polydViewType polyd,
                             const ordinal_type n,
                             const double alpha,
                             const double beta) {
    const double one = 1.0;
    if (n == 0)
      for(ordinal_type i = 0; i < np; ++i)
        polyd(i) = 0.0;
    else {
      Kokkos::View<typename polydViewType::value_type*,Kokkos::AnonymousSpace,Kokkos::MemoryUnmanaged> null;
      JacobiPolynomial(np, z, polyd, null, n-1, alpha+one, beta+one);
      for(ordinal_type i = 0; i < np; ++i)
        polyd(i) *= 0.5*(alpha + beta + n + one);
    }
  }

  // -----------------------------------------------------------------------

  template<typename zViewType,
           bool DeflationEnabled>
  KOKKOS_INLINE_FUNCTION
  void
  Polylib::Serial::
  JacobiZeros(         zViewType z,
              const ordinal_type n,
              const double alpha,
              const double beta) {
    if (DeflationEnabled)
      JacobiZerosPolyDeflation(z, n, alpha, beta);
    else
      JacobiZerosTriDiagonal(z, n, alpha, beta);
  }

  template<typename zViewType>
  KOKKOS_INLINE_FUNCTION
  void
  Polylib::Serial::
  JacobiZerosPolyDeflation(      zViewType z,
                           const ordinal_type n,
                           const double alpha,
                           const double beta) {
    if(!n)
      return;

    const double dth = M_PI/(2.0*n);
    const double one = 1.0, two = 2.0;
    const double tol = tolerence();

    typedef typename zViewType::value_type value_type;
    typedef typename zViewType::pointer_type pointer_type;

    value_type r_buf, poly_buf, pder_buf;
    Kokkos::View<value_type*,Kokkos::AnonymousSpace,Kokkos::MemoryUnmanaged> 
      poly((pointer_type)&poly_buf, 1), pder((pointer_type)&pder_buf, 1), r((pointer_type)&r_buf, 1); 

    value_type rlast = 0.0;
    for (auto k = 0; k < n; ++k) {
      r(0) = -cos((two*(double)k + one) * dth);
      if (k) 
        r(0) = 0.5*(r(0) + rlast);

      for (ordinal_type j = 1; j < MaxPolylibIteration; ++j) {
        JacobiPolynomial(1, r, poly, pder, n, alpha, beta);

        value_type sum = 0;
        for (ordinal_type i = 0; i < k; ++i) 
          sum += one/(r(0) - z(i));

        const value_type delr = -poly(0) / (pder(0) - sum * poly(0));
        r(0) += delr;

        if( Util<value_type>::abs(delr) < tol ) 
          break;
      }
      z(k)  = r(0);
      rlast = r(0);
    }
  }
  
  template<typename aViewType>
  KOKKOS_INLINE_FUNCTION
  void
  Polylib::Serial::
  JacobiZerosTriDiagonal(      aViewType a,
                         const ordinal_type n,
                         const double alpha,
                         const double beta) {
    if(!n)
      return;

    typedef typename aViewType::value_type value_type;
    typedef typename aViewType::pointer_type pointer_type;

    value_type b_buf[MaxPolylibPoint];
    Kokkos::View<value_type*,Kokkos::AnonymousSpace,Kokkos::MemoryUnmanaged> 
      b((pointer_type)&b_buf[0], MaxPolylibPoint);

    // generate normalised terms
    const double apb  = alpha + beta;
    double apbi = 2.0 + apb;

    b(n-1) = pow(2.0,apb+1.0)*GammaFunction(alpha+1.0)*GammaFunction(beta+1.0)/GammaFunction(apbi);
    a(0)   = (beta-alpha)/apbi;
    b(0)   = sqrt(4.0*(1.0+alpha)*(1.0+beta)/((apbi+1.0)*apbi*apbi));

    const double a2b2 = beta*beta-alpha*alpha;
    for (ordinal_type i = 1; i < n-1; ++i) {
      apbi = 2.0*(i+1) + apb;
      a(i) = a2b2/((apbi-2.0)*apbi);
      b(i) = sqrt(4.0*(i+1)*(i+1+alpha)*(i+1+beta)*(i+1+apb)/
                  ((apbi*apbi-1)*apbi*apbi));
    }

    apbi   = 2.0*n + apb;
    //a(n-1) = a2b2/((apbi-2.0)*apbi); // THIS IS A BUG!!!
    if (n>1) a(n-1) = a2b2/((apbi-2.0)*apbi);

    // find eigenvalues
    TriQL(a, b, n);
  }


  template<typename dViewType,
           typename eViewType>
  KOKKOS_INLINE_FUNCTION
  void
  Polylib::Serial::
  TriQL(      dViewType d,
              eViewType e,
        const ordinal_type n) {
    ordinal_type m,l,iter,i,k;

    typedef typename dViewType::value_type value_type;
    value_type s,r,p,g,f,dd,c,b;

    for (l=0; l<n; ++l) {
      iter=0;
      do {
        for (m=l; m<n-1; ++m) {
          dd=Util<value_type>::abs(d(m))+Util<value_type>::abs(d(m+1));
          if (Util<value_type>::abs(e(m))+dd == dd) break;
        }
        if (m != l) {
          if (iter++ == MaxPolylibIteration) {
            INTREPID2_TEST_FOR_ABORT(true,
                                     ">>> ERROR (Polylib::Serial): Too many iterations in TQLI.");
          }
          g=(d(l+1)-d(l))/(2.0*e(l));
          r=sqrt((g*g)+1.0);
          //g=d(m)-d(l)+e(l)/(g+sign(r,g));
          g=d(m)-d(l)+e(l)/(g+((g)<0 ? value_type(-Util<value_type>::abs(r)) : Util<value_type>::abs(r)));
          s=c=1.0;
          p=0.0;
          for (i=m-1; i>=l; i--) {
            f=s*e(i);
            b=c*e(i);
            if (Util<value_type>::abs(f) >= Util<value_type>::abs(g)) {
              c=g/f;
              r=sqrt((c*c)+1.0);
              e(i+1)=f*r;
              c *= (s=1.0/r);
            } else {
              s=f/g;
              r=sqrt((s*s)+1.0);
              e(i+1)=g*r;
              s *= (c=1.0/r);
            }
            g=d(i+1)-p;
            r=(d(i)-g)*s+2.0*c*b;
            p=s*r;
            d(i+1)=g+p;
            g=c*r-b;
          }
          d(l)=d(l)-p;
          e(l)=g;
          e(m)=0.0;
        }
      } while (m != l);
    }

    // order eigenvalues
    for (i = 0; i < n-1; ++i) {
      k = i;
      p = d(i);
      for (l = i+1; l < n; ++l)
        if (d(l) < p) {
          k = l;
          p = d(l);
        }
      d(k) = d(i);
      d(i) = p;
    }
  }

  KOKKOS_INLINE_FUNCTION
  double
  Polylib::Serial::
  GammaFunction(const double x) {
    double gamma(0);

    if      (x == -0.5) gamma = -2.0*sqrt(M_PI);
    else if (x ==  0.0) gamma = 1.0;
    else if ((x-(ordinal_type)x) == 0.5) {
      ordinal_type n = (ordinal_type) x;
      auto tmp = x;

      gamma = sqrt(M_PI);
      while (n--) {
        tmp   -= 1.0;
        gamma *= tmp;
      }
    } else if ((x-(ordinal_type)x) == 0.0) {
      ordinal_type n = (ordinal_type) x;
      auto tmp = x;

      gamma = 1.0;
      while (--n) {
        tmp   -= 1.0;
        gamma *= tmp;
      }
    } else {
      INTREPID2_TEST_FOR_ABORT(true,
                               ">>> ERROR (Polylib::Serial): Argument is not of integer or half order.");
    }
    return gamma;
  }

} // end of namespace Intrepid2

#endif
