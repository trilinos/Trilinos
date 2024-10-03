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
// File: Intrepid_Polylib.hpp
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


/** \file   Intrepid2_Polylib.hpp
    \brief  Header file for Intrepid2::Polylib class providing orthogonal polynomial
    calculus and interpolation.
    \author Created by Spencer Sherwin, Aeronautics, Imperial College London,
    modified and redistributed by D. Ridzal.
    Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_POLYLIB_HPP__
#define __INTREPID2_POLYLIB_HPP__

#include "Intrepid2_ConfigDefs.hpp"

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Kokkos_Core.hpp"

namespace Intrepid2 {

  /**
     \page pagePolylib The Polylib library
     \section sectionPolyLib Routines For Orthogonal Polynomial Calculus and Interpolation

     Spencer Sherwin, 
     Aeronautics, Imperial College London

     Based on codes by Einar Ronquist and Ron Henderson

     Abbreviations
     - z    -   Set of collocation/quadrature points
     - w    -   Set of quadrature weights
     - D    -   Derivative matrix
     - h    -   Lagrange Interpolant
     - I    -   Interpolation matrix
     - g    -   Gauss
     - gr   -   Gauss-Radau
     - gl   -   Gauss-Lobatto
     - j    -   Jacobi
     - m    -   point at minus 1 in Radau rules
     - p    -   point at plus  1 in Radau rules

     -----------------------------------------------------------------------\n
     MAIN     ROUTINES\n
     -----------------------------------------------------------------------\n

     Points and Weights:

     - zwgj        Compute Gauss-Jacobi         points and weights
     - zwgrjm      Compute Gauss-Radau-Jacobi   points and weights (z=-1)
     - zwgrjp      Compute Gauss-Radau-Jacobi   points and weights (z= 1)
     - zwglj       Compute Gauss-Lobatto-Jacobi points and weights

     Derivative Matrices:

     - Dgj         Compute Gauss-Jacobi         derivative matrix
     - Dgrjm       Compute Gauss-Radau-Jacobi   derivative matrix (z=-1)
     - Dgrjp       Compute Gauss-Radau-Jacobi   derivative matrix (z= 1)
     - Dglj        Compute Gauss-Lobatto-Jacobi derivative matrix

     Lagrange Interpolants:

     - hgj         Compute Gauss-Jacobi         Lagrange interpolants
     - hgrjm       Compute Gauss-Radau-Jacobi   Lagrange interpolants (z=-1)
     - hgrjp       Compute Gauss-Radau-Jacobi   Lagrange interpolants (z= 1)
     - hglj        Compute Gauss-Lobatto-Jacobi Lagrange interpolants

     Interpolation Operators:

     - Imgj        Compute interpolation operator gj->m
     - Imgrjm      Compute interpolation operator grj->m (z=-1)
     - Imgrjp      Compute interpolation operator grj->m (z= 1)
     - Imglj       Compute interpolation operator glj->m

     Polynomial Evaluation:

     - jacobfd     Returns value and derivative of Jacobi poly. at point z
     - jacobd      Returns derivative of Jacobi poly. at point z (valid at z=-1,1)

     -----------------------------------------------------------------------\n
     LOCAL      ROUTINES\n
     -----------------------------------------------------------------------\n

     - jacobz      Returns Jacobi polynomial zeros
     - gammaf      Gamma function for integer values and halves



     ------------------------------------------------------------------------\n

     Useful references:

     - [1] Gabor Szego: Orthogonal Polynomials, American Mathematical Society,
     Providence, Rhode Island, 1939.
     - [2] Abramowitz \& Stegun: Handbook of Mathematical Functions,
     Dover, New York, 1972.
     - [3] Canuto, Hussaini, Quarteroni \& Zang: Spectral Methods in Fluid
     Dynamics, Springer-Verlag, 1988.
     - [4] Ghizzetti \& Ossicini: Quadrature Formulae, Academic Press, 1970.
     - [5] Karniadakis \& Sherwin: Spectral/hp element methods for CFD, 1999


     NOTES

     -# Legendre  polynomial \f$ \alpha = \beta = 0 \f$ 
     -# Chebychev polynomial \f$ \alpha = \beta = -0.5 \f$
     -# All array subscripts start from zero, i.e. vector[0..N-1] 
  */


  /** \class Intrepid2::Polylib
      \brief Providing orthogonal polynomial calculus and interpolation,
      created by Spencer Sherwin, Aeronautics, Imperial College London,
      modified and redistributed by D. Ridzal.

      See \ref pagePolylib "original Polylib documentation".
  */
  class Polylib {
  public:

    static constexpr ordinal_type  MaxPolylibIteration = 50;
    static constexpr ordinal_type  MaxPolylibOrder =
        (Parameters::MaxOrder > Parameters::MaxCubatureDegreeEdge) ? Parameters::MaxOrder :
                                                                     Parameters::MaxCubatureDegreeEdge;
    // NVR: HVOL bases on tri/tet use Polylib with order + spaceDim + 2 points; in 3D this can be up to Parameters::MaxOrder + 5.
    static constexpr ordinal_type  MaxPolylibPoint = (MaxPolylibOrder/2+2 > Parameters::MaxOrder + 5) ? MaxPolylibOrder/2+2 : Parameters::MaxOrder + 5;

    struct Serial {

      // -----------------------------------------------------------------------
      // Points and Weights
      // -----------------------------------------------------------------------
      
      /** \brief  Gauss-Jacobi/Gauss-Radau-Jacobi/Gauss-Lobatto zeros and weights.
          
          \li Generate \a np Gauss Jacobi zeros, \a z, and weights,\a w,
          associated with the Jacobi polynomial \f$ P^{\alpha,\beta}_{np}(z)\f$,
          
          \li Exact for polynomials of order \a 2np-1 or less

          \li Generate \a np Gauss-Radau-Jacobi zeros, \a z, and weights,\a w,
          associated with the polynomial \f$(1+z) P^{\alpha,\beta+1}_{np-1}(z)\f$.

          \li  Exact for polynomials of order \a 2np-2 or less

          \li Generate \a np Gauss-Lobatto-Jacobi points, \a z, and weights, \a w,
          associated with polynomial \f$ (1-z)(1+z) P^{\alpha+1,\beta+1}_{np-2}(z) \f$

          \li Exact for polynomials of order \a 2np-3 or less
      */
      template<EPolyType polyType>
      struct Cubature {
        template<typename zViewType,
                 typename wViewType>
        KOKKOS_INLINE_FUNCTION
        static void
        getValues(      zViewType z, 
                        wViewType w, 
                  const ordinal_type np, 
                  const double alpha, 
                  const double beta);
      };

      template<typename zViewType,
               typename wViewType>
      KOKKOS_INLINE_FUNCTION
      static void
      getCubature(      zViewType z, 
                        wViewType w, 
                  const ordinal_type np, 
                  const double alpha, 
                  const double beta,
                  const EPolyType poly) {
        switch (poly) {
        case POLYTYPE_GAUSS:             Polylib::Serial::Cubature<POLYTYPE_GAUSS>            ::getValues(z, w, np, alpha, beta); break;
        case POLYTYPE_GAUSS_RADAU_LEFT:  Polylib::Serial::Cubature<POLYTYPE_GAUSS_RADAU_LEFT> ::getValues(z, w, np, alpha, beta); break;
        case POLYTYPE_GAUSS_RADAU_RIGHT: Polylib::Serial::Cubature<POLYTYPE_GAUSS_RADAU_RIGHT>::getValues(z, w, np, alpha, beta); break;
        case POLYTYPE_GAUSS_LOBATTO:     Polylib::Serial::Cubature<POLYTYPE_GAUSS_LOBATTO>    ::getValues(z, w, np, alpha, beta); break;
        default:
          INTREPID2_TEST_FOR_ABORT(true,
                                   ">>> ERROR (Polylib::Serial::getCubature): Not supported poly type.");
          break;
        }
      }
     
      // -----------------------------------------------------------------------
      // Derivative Matrix
      // -----------------------------------------------------------------------

      /** \brief Compute the Derivative Matrix and its transpose associated
          with the Gauss-Jacobi/Gauss-Radau-Jacobi/Gauss-Lobatto-Jacobi zeros.

          \li Compute the derivative matrix \a D associated with the n_th order Lagrangian
          interpolants through the \a np Gauss-Jacobi points \a z such that \n
          \f$  \frac{du}{dz}(z[i]) =  \sum_{j=0}^{np-1} D[i*np+j] u(z[j]) \f$

      */
      template<EPolyType polyType>
      struct Derivative {
        template<typename DViewType,
                 typename zViewType>
        KOKKOS_INLINE_FUNCTION
        static void
        getValues(      DViewType D,  
                  const zViewType z, 
                  const ordinal_type np, 
                  const double alpha, 
                  const double beta);
      };

      template<typename DViewType,
               typename zViewType>
      KOKKOS_INLINE_FUNCTION
      static void
      getDerivative(      DViewType D,  
                    const zViewType z, 
                    const ordinal_type np, 
                    const double alpha, 
                    const double beta,
                    const EPolyType poly) {
        switch (poly) {
        case POLYTYPE_GAUSS:             Polylib::Serial::Derivative<POLYTYPE_GAUSS>            ::getValues(D, z, np, alpha, beta); break;
        case POLYTYPE_GAUSS_RADAU_LEFT:  Polylib::Serial::Derivative<POLYTYPE_GAUSS_RADAU_LEFT> ::getValues(D, z, np, alpha, beta); break;
        case POLYTYPE_GAUSS_RADAU_RIGHT: Polylib::Serial::Derivative<POLYTYPE_GAUSS_RADAU_RIGHT>::getValues(D, z, np, alpha, beta); break;
        case POLYTYPE_GAUSS_LOBATTO:     Polylib::Serial::Derivative<POLYTYPE_GAUSS_LOBATTO>    ::getValues(D, z, np, alpha, beta); break;
        default:
          INTREPID2_TEST_FOR_ABORT(true,
                                   ">>> ERROR (Polylib::Serial::getDerivative): Not supported poly type.");
          break;
        }
      }

      // -----------------------------------------------------------------------
      // Lagrangian Interpolants 
      // -----------------------------------------------------------------------

      /** \brief Compute the value of the \a i th Lagrangian interpolant through
          the \a np Gauss-Jacobi/Gauss-Radau-Jacobi/Gauss-Lobatto points \a zgj 
          at the arbitrary location \a z.

          \li \f$ -1 \leq z \leq 1 \f$

          \li POLYTYPE_GAUSS Uses the defintion of the Lagrangian interpolant:\n
          \f$
          \begin{array}{rcl}
          h_j(z) =  \left\{ \begin{array}{ll}
          \displaystyle \frac{P_{np}^{\alpha,\beta}(z)}
          {[P_{np}^{\alpha,\beta}(z_j)]^\prime
          (z-z_j)} & \mbox{if $z \ne z_j$}\\
          & \\
          1 & \mbox{if $z=z_j$}
          \end{array}
          \right.
          \end{array}
          \f$

          \li Gauss-Radau-Jacobi (Left) Uses the defintion of the Lagrangian interpolant:\n
          %
          \f$ \begin{array}{rcl}
          h_j(z) = \left\{ \begin{array}{ll}
          \displaystyle \frac{(1+z) P_{np-1}^{\alpha,\beta+1}(z)}
          {((1+z_j) [P_{np-1}^{\alpha,\beta+1}(z_j)]^\prime +
          P_{np-1}^{\alpha,\beta+1}(z_j) ) (z-z_j)} & \mbox{if $z \ne z_j$}\\
          & \\
          1 & \mbox{if $z=z_j$}
          \end{array}
          \right.
          \end{array}   \f$

          \li Gauss-Radau-Jacobi (Right) Uses the defintion of the Lagrangian interpolant:\n
          %
          \f$ \begin{array}{rcl}
          h_j(z) = \left\{ \begin{array}{ll}
          \displaystyle \frac{(1-z) P_{np-1}^{\alpha+1,\beta}(z)}
          {((1-z_j) [P_{np-1}^{\alpha+1,\beta}(z_j)]^\prime -
          P_{np-1}^{\alpha+1,\beta}(z_j) ) (z-z_j)} & \mbox{if $z \ne z_j$}\\
          & \\
          1 & \mbox{if $z=z_j$}
          \end{array}
          \right.
          \end{array}   \f$

          \li Gauss-Lobatto Uses the defintion of the Lagrangian interpolant:\n
          %
          \f$ \begin{array}{rcl}
          h_j(z) = \left\{ \begin{array}{ll}
          \displaystyle \frac{(1-z^2) P_{np-2}^{\alpha+1,\beta+1}(z)}
          {((1-z^2_j) [P_{np-2}^{\alpha+1,\beta+1}(z_j)]^\prime -
          2 z_j P_{np-2}^{\alpha+1,\beta+1}(z_j) ) (z-z_j)}&\mbox{if $z \ne z_j$}\\
          & \\
          1 & \mbox{if $z=z_j$}
          \end{array}
          \right.
          \end{array}   \f$
      */
      template<EPolyType polyType>
      struct LagrangianInterpolant {
        template<typename zViewType>
        KOKKOS_INLINE_FUNCTION
        static typename zViewType::value_type
        getValue(const ordinal_type i, 
                 const typename zViewType::value_type z, 
                 const zViewType zgj,
                 const ordinal_type np,
                 const double alpha,
                 const double beta);
      };

      template<typename zViewType>
      KOKKOS_INLINE_FUNCTION
      static typename zViewType::value_type
      getLagrangianInterpolant(const ordinal_type i, 
                               const typename zViewType::value_type z, 
                               const zViewType zgj,
                               const ordinal_type np,
                               const double alpha,
                               const double beta,
                               const EPolyType poly) {
        typename zViewType::value_type r_val = 0;
        switch (poly) {
        case POLYTYPE_GAUSS:             r_val = Polylib::Serial::LagrangianInterpolant<POLYTYPE_GAUSS>            ::getValue(i, z, zgj, np, alpha, beta); break;
        case POLYTYPE_GAUSS_RADAU_LEFT:  r_val = Polylib::Serial::LagrangianInterpolant<POLYTYPE_GAUSS_RADAU_LEFT> ::getValue(i, z, zgj, np, alpha, beta); break;
        case POLYTYPE_GAUSS_RADAU_RIGHT: r_val = Polylib::Serial::LagrangianInterpolant<POLYTYPE_GAUSS_RADAU_RIGHT>::getValue(i, z, zgj, np, alpha, beta); break;
        case POLYTYPE_GAUSS_LOBATTO:     r_val = Polylib::Serial::LagrangianInterpolant<POLYTYPE_GAUSS_LOBATTO>    ::getValue(i, z, zgj, np, alpha, beta); break;
        default:
          INTREPID2_TEST_FOR_ABORT(true,
                                   ">>> ERROR (Polylib::Serial::getLagrangianInterpolant): Not supported poly type.");
          break;
        }
        return r_val;
      }

      // -----------------------------------------------------------------------
      // Interpolation operators 
      // -----------------------------------------------------------------------

      /** \brief Interpolation Operator from Gauss-Jacobi points to an
          arbitrary distribution at points \a zm

          \li Computes the one-dimensional interpolation matrix, \a im, to
          interpolate a function from at Gauss-Jacobi distribution of \a nz
          zeros \a zgj to an arbitrary distribution of \a mz points \a zm, i.e.\n
          \f$
          u(zm[i]) = \sum_{j=0}^{nz-1} im[i*nz+j] \ u(zgj[j])
          \f$
      */
      template<EPolyType polyType>
      struct InterpolationOperator {
        template<typename imViewType,
                 typename zgrjViewType,
                 typename zmViewType>
        KOKKOS_INLINE_FUNCTION
        static void 
        getMatrix(      imViewType im,
                  const zgrjViewType zgrj,
                  const zmViewType zm,
                  const ordinal_type nz,
                  const ordinal_type mz,
                  const double alpha,
                  const double beta);
      };

      template<typename imViewType,
               typename zgrjViewType,
               typename zmViewType>
      KOKKOS_INLINE_FUNCTION
      static void
      getInterpolationOperator(      imViewType im,
                               const zgrjViewType zgrj,
                               const zmViewType zm,
                               const ordinal_type nz,
                               const ordinal_type mz,
                               const double alpha,
                               const double beta,
                               const EPolyType poly) {
        switch (poly) {
        case POLYTYPE_GAUSS:             Polylib::Serial::InterpolationOperator<POLYTYPE_GAUSS>            ::getMatrix(im, zgrj, zm, nz, mz, alpha, beta); break;
        case POLYTYPE_GAUSS_RADAU_LEFT:  Polylib::Serial::InterpolationOperator<POLYTYPE_GAUSS_RADAU_LEFT> ::getMatrix(im, zgrj, zm, nz, mz, alpha, beta); break;
        case POLYTYPE_GAUSS_RADAU_RIGHT: Polylib::Serial::InterpolationOperator<POLYTYPE_GAUSS_RADAU_RIGHT>::getMatrix(im, zgrj, zm, nz, mz, alpha, beta); break;
        case POLYTYPE_GAUSS_LOBATTO:     Polylib::Serial::InterpolationOperator<POLYTYPE_GAUSS_LOBATTO>    ::getMatrix(im, zgrj, zm, nz, mz, alpha, beta); break;
        default:
          INTREPID2_TEST_FOR_ABORT(true,
                                   ">>> ERROR (Polylib::Serial::getInterpolationOperator): Not supported poly type.");
          break;
        }
      }

      
      // -----------------------------------------------------------------------
      // Polynomial functions 
      // -----------------------------------------------------------------------

      /** \brief Routine to calculate Jacobi polynomials, \f$
          P^{\alpha,\beta}_n(z) \f$, and their first derivative, \f$
          \frac{d}{dz} P^{\alpha,\beta}_n(z) \f$.

          \li This function returns the vectors \a poly_in and \a poly_d
          containing the value of the \a n-th order Jacobi polynomial
          \f$ P^{\alpha,\beta}_n(z) \alpha > -1, \beta > -1 \f$ and its
          derivative at the \a np points in \a z[i]

          - If \a poly_in = NULL then only calculate derivative

          - If \a polyd   = NULL then only calculate polynomial

          - To calculate the polynomial this routine uses the recursion
          relationship (see appendix A ref [4]) :
          \f$ \begin{array}{rcl}
          P^{\alpha,\beta}_0(z) &=& 1 \\
          P^{\alpha,\beta}_1(z) &=& \frac{1}{2} [ \alpha-\beta+(\alpha+\beta+2)z] \\
          a^1_n P^{\alpha,\beta}_{n+1}(z) &=& (a^2_n + a^3_n z)
          P^{\alpha,\beta}_n(z) - a^4_n P^{\alpha,\beta}_{n-1}(z) \\
          a^1_n &=& 2(n+1)(n+\alpha + \beta + 1)(2n + \alpha + \beta) \\
          a^2_n &=& (2n + \alpha + \beta + 1)(\alpha^2 - \beta^2)  \\
          a^3_n &=& (2n + \alpha + \beta)(2n + \alpha + \beta + 1)
          (2n + \alpha + \beta + 2)  \\
          a^4_n &=& 2(n+\alpha)(n+\beta)(2n + \alpha + \beta + 2)
          \end{array} \f$

          - To calculate the derivative of the polynomial this routine uses
          the relationship (see appendix A ref [4]) :
          \f$ \begin{array}{rcl}
          b^1_n(z)\frac{d}{dz} P^{\alpha,\beta}_n(z)&=&b^2_n(z)P^{\alpha,\beta}_n(z)
          + b^3_n(z) P^{\alpha,\beta}_{n-1}(z) \hspace{2.2cm} \\
          b^1_n(z) &=& (2n+\alpha + \beta)(1-z^2) \\
          b^2_n(z) &=& n[\alpha - \beta - (2n+\alpha + \beta)z]\\
          b^3_n(z) &=& 2(n+\alpha)(n+\beta)
          \end{array} \f$

          - Note the derivative from this routine is only valid for -1 < \a z < 1.
      */
      template<typename zViewType,
               typename polyiViewType,
               typename polydViewType>
      KOKKOS_INLINE_FUNCTION
      static void 
      JacobiPolynomial(const ordinal_type np, 
                       const zViewType z, 
                             polyiViewType poly_in, 
                             polydViewType polyd,
                       const ordinal_type n, 
                       const double alpha, 
                       const double beta);


      /** \brief Calculate the  derivative of Jacobi polynomials

          \li Generates a vector \a poly of values of the derivative of the
          \a n-th order Jacobi polynomial \f$ P^(\alpha,\beta)_n(z)\f$ at the
          \a np points \a z.

          \li To do this we have used the relation
          \n
          \f$ \frac{d}{dz} P^{\alpha,\beta}_n(z)
          = \frac{1}{2} (\alpha + \beta + n + 1)  P^{\alpha,\beta}_n(z) \f$

          \li This formulation is valid for \f$ -1 \leq z \leq 1 \f$
      */
      template<typename zViewType,
               typename polydViewType>
      KOKKOS_INLINE_FUNCTION
      static void 
      JacobiPolynomialDerivative(const ordinal_type np, 
                                 const zViewType z, 
                                       polydViewType polyd, 
                                 const ordinal_type n,
                                 const double alpha, 
                                 const double beta);

      // -----------------------------------------------------------------------
      // Helper functions. 
      // -----------------------------------------------------------------------
      
      /** \brief  Calculate the \a n zeros, \a z, of the Jacobi polynomial, i.e.
          \f$ P_n^{\alpha,\beta}(z) = 0 \f$

          This routine is only valid for \f$( \alpha > -1, \beta > -1)\f$
          and uses polynomial deflation in a Newton iteration
      */
      template<typename zViewType,
               bool DeflationEnabled = false>
      KOKKOS_INLINE_FUNCTION
      static void   
      JacobiZeros (      zViewType z,
                   const ordinal_type n, 
                   const double alpha, 
                   const double beta);

      template<typename zViewType>
      KOKKOS_INLINE_FUNCTION
      static void
      JacobiZerosPolyDeflation(      zViewType z,
                               const ordinal_type n,
                               const double alpha,
                               const double beta);
      
      template<typename aViewType>
      KOKKOS_INLINE_FUNCTION
      static void
      JacobiZerosTriDiagonal(      aViewType a,
                             const ordinal_type n,
                             const double alpha,
                             const double beta); 

      /** \brief Zero determination through the eigenvalues of a tridiagonal
          matrix from the three term recursion relationship.

          Set up a symmetric tridiagonal matrix

          \f$ \left [  \begin{array}{ccccc}
          a[0] & b[0]   &        &        & \\
          b[0] & a[1]   & b[1]   &        & \\
          0   & \ddots & \ddots & \ddots &  \\
          &        & \ddots & \ddots & b[n-2] \\
          &        &        & b[n-2] & a[n-1] \end{array} \right ] \f$

          Where the coefficients a[n], b[n] come from the  recurrence relation

          \f$  b_j p_j(z) = (z - a_j ) p_{j-1}(z) - b_{j-1}   p_{j-2}(z) \f$

          where \f$ j=n+1\f$ and \f$p_j(z)\f$ are the Jacobi (normalized)
          orthogonal polynomials \f$ \alpha,\beta > -1\f$( integer values and
          halves). Since the polynomials are orthonormalized, the tridiagonal
          matrix is guaranteed to be symmetric. The eigenvalues of this
          matrix are the zeros of the Jacobi polynomial.
      */
      template<typename aViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static void   
      JacobiZeros(aViewType a,
                  bViewType b,
                  const ordinal_type n, 
                  const double alpha, 
                  const double beta);


      /** \brief QL algorithm for symmetric tridiagonal matrix

          This subroutine is a translation of an algol procedure,
          num. math. \b 12, 377-383(1968) by martin and wilkinson, as modified
          in num. math. \b 15, 450(1970) by dubrulle.  Handbook for
          auto. comp., vol.ii-linear algebra, 241-248(1971).  This is a
          modified version from numerical recipes.

          This subroutine finds the eigenvalues and first components of the
          eigenvectors of a symmetric tridiagonal matrix by the implicit QL
          method.

          on input:
          - n is the order of the matrix;
          - d contains the diagonal elements of the input matrix;
          - e contains the subdiagonal elements of the input matrix
          in its first n-1 positions. e(n) is arbitrary;

          on output:

          - d contains the eigenvalues in ascending order.
          - e has been destroyed;
      */
      template<typename dViewType,
               typename eViewType>
      KOKKOS_INLINE_FUNCTION      
      static void   
      TriQL(      dViewType d, 
                  eViewType e,
            const ordinal_type n);


      /** \brief Calculate the Gamma function , \f$ \Gamma(x)\f$, for integer
          values \a x and halves.

          Determine the value of \f$\Gamma(x)\f$ using:

          \f$ \Gamma(x) = (x-1)!  \mbox{ or  }  \Gamma(x+1/2) = (x-1/2)\Gamma(x-1/2)\f$

          where \f$ \Gamma(1/2) = \sqrt{\pi}\f$
      */
      KOKKOS_INLINE_FUNCTION      
      static double
      GammaFunction(const double x);

    };

    // -----------------------------------------------------------------------
  };
  
} // end of Intrepid namespace

  // include templated definitions
#include <Intrepid2_PolylibDef.hpp>

#endif
