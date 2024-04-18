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


/** \file   Intrepid_Polylib.hpp
    \brief  Header file for a set of functions providing orthogonal polynomial
            polynomial calculus and interpolation.
    \author Created by Spencer Sherwin, Aeronautics, Imperial College London,
            modified and redistributed by D. Ridzal.
*/

#ifndef INTREPID_POLYLIB_HPP
#define INTREPID_POLYLIB_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Types.hpp"
#include "Teuchos_Assert.hpp"

namespace Intrepid {

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


  /** \enum  Intrepid::EIntrepidPLPoly
      \brief Enumeration of coordinate frames (reference/ambient) for geometrical entities (cells, points).
  */
  enum EIntrepidPLPoly {
    PL_GAUSS=0,
    PL_GAUSS_RADAU_LEFT,
    PL_GAUSS_RADAU_RIGHT,
    PL_GAUSS_LOBATTO,
    PL_MAX
  };

  inline EIntrepidPLPoly & operator++(EIntrepidPLPoly &type) {
    return type = static_cast<EIntrepidPLPoly>(type+1);
  }

  inline EIntrepidPLPoly operator++(EIntrepidPLPoly &type, int) {
    EIntrepidPLPoly oldval = type;
    ++type;
    return oldval;
  }


  /** \class Intrepid::IntrepidPolylib
      \brief Providing orthogonal polynomial calculus and interpolation,
             created by Spencer Sherwin, Aeronautics, Imperial College London,
             modified and redistributed by D. Ridzal.

             See \ref pagePolylib "original Polylib documentation".
  */
  class IntrepidPolylib {

    public:

    /* Points and weights */

    /** \brief  Gauss-Jacobi zeros and weights.

    \li Generate \a np Gauss Jacobi zeros, \a z, and weights,\a w,
    associated with the Jacobi polynomial \f$ P^{\alpha,\beta}_{np}(z)\f$,

    \li Exact for polynomials of order \a 2np-1 or less
    */
    template<class Scalar>
    static void   zwgj   (Scalar *z, Scalar *w, const int np, const Scalar alpha, const Scalar beta);


    /** \brief  Gauss-Radau-Jacobi zeros and weights with end point at \a z=-1.

    \li Generate \a np Gauss-Radau-Jacobi zeros, \a z, and weights,\a w,
    associated with the polynomial \f$(1+z) P^{\alpha,\beta+1}_{np-1}(z)\f$.

    \li  Exact for polynomials of order \a 2np-2 or less
    */
    template<class Scalar>
    static void   zwgrjm (Scalar *z, Scalar *w, const int np, const Scalar alpha, const Scalar beta);


    /** \brief  Gauss-Radau-Jacobi zeros and weights with end point at \a z=1

    \li Generate \a np Gauss-Radau-Jacobi zeros, \a z, and weights,\a w,
    associated with the  polynomial \f$(1-z) P^{\alpha+1,\beta}_{np-1}(z)\f$.

    \li Exact for polynomials of order \a 2np-2 or less
    */
    template<class Scalar>
    static void   zwgrjp (Scalar *z, Scalar *w, const int np, const Scalar alpha, const Scalar beta);


    /** \brief  Gauss-Lobatto-Jacobi zeros and weights with end point at \a z=-1,\a 1

    \li Generate \a np Gauss-Lobatto-Jacobi points, \a z, and weights, \a w,
    associated with polynomial \f$ (1-z)(1+z) P^{\alpha+1,\beta+1}_{np-2}(z) \f$

    \li Exact for polynomials of order \a 2np-3 or less
    */
    template<class Scalar>
    static void   zwglj  (Scalar *z, Scalar *w, const int np, const Scalar alpha, const Scalar beta);



    /* Derivative operators */

    /** \brief Compute the Derivative Matrix and its transpose associated
               with the Gauss-Jacobi zeros.

    \li Compute the derivative matrix \a D associated with the n_th order Lagrangian
        interpolants through the \a np Gauss-Jacobi points \a z such that \n
        \f$  \frac{du}{dz}(z[i]) =  \sum_{j=0}^{np-1} D[i*np+j] u(z[j]) \f$

    */
    template<class Scalar>
    static void   Dgj    (Scalar *D,  const Scalar *z, const int np, const Scalar alpha, const Scalar beta);


    /** \brief Compute the Derivative Matrix and its transpose associated
               with the Gauss-Radau-Jacobi zeros with a zero at \a z=-1.

    \li Compute the derivative matrix \a D associated with the n_th
        order Lagrangian interpolants through the \a np Gauss-Radau-Jacobi
        points \a z such that \n \f$ \frac{du}{dz}(z[i]) =
        \sum_{j=0}^{np-1} D[i*np+j] u(z[j]) \f$
    */
    template<class Scalar>
    static void   Dgrjm  (Scalar *D, const Scalar *z, const int np, const Scalar alpha, const Scalar beta);


    /** \brief Compute the Derivative Matrix  associated with the
               Gauss-Radau-Jacobi zeros with a zero at \a z=1.

    \li Compute the derivative matrix \a D associated with the n_th
        order Lagrangian interpolants through the \a np Gauss-Radau-Jacobi
        points \a z such that \n \f$ \frac{du}{dz}(z[i]) =
        \sum_{j=0}^{np-1} D[i*np+j] u(z[j]) \f$
    */
    template<class Scalar>
    static void   Dgrjp  (Scalar *D, const Scalar *z, const int np, const Scalar alpha, const Scalar beta);


    /** \brief Compute the Derivative Matrix associated with the
               Gauss-Lobatto-Jacobi zeros.

    \li Compute the derivative matrix \a D associated with the n_th
        order Lagrange interpolants through the \a np
        Gauss-Lobatto-Jacobi points \a z such that \n \f$
        \frac{du}{dz}(z[i]) = \sum_{j=0}^{np-1} D[i*np+j] u(z[j]) \f$
    */
    template<class Scalar>
    static void   Dglj   (Scalar *D, const Scalar *z, const int np, const Scalar alpha, const Scalar beta);



    /* Lagrangian interpolants */

    /** \brief Compute the value of the \a i th Lagrangian interpolant through
               the \a np Gauss-Jacobi points \a zgj at the arbitrary location \a z.

    \li \f$ -1 \leq z \leq 1 \f$

    \li Uses the defintion of the Lagrangian interpolant:\n
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
    */
    template<class Scalar>
    static Scalar hgj     (const int i, const Scalar z, const Scalar *zgj,
                           const int np, const Scalar alpha, const Scalar beta);


    /** \brief Compute the value of the \a i th Lagrangian interpolant through the
               \a np Gauss-Radau-Jacobi points \a zgrj at the arbitrary location
               \a z. This routine assumes \a zgrj includes the point \a -1.

    \li \f$ -1 \leq z \leq 1 \f$

    \li Uses the defintion of the Lagrangian interpolant:\n
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
    */
    template<class Scalar>
    static Scalar hgrjm   (const int i, const Scalar z, const Scalar *zgrj,
                           const int np, const Scalar alpha, const Scalar beta);


    /** \brief Compute the value of the \a i th Lagrangian interpolant through the
               \a np Gauss-Radau-Jacobi points \a zgrj at the arbitrary location
               \a z. This routine assumes \a zgrj includes the point \a +1.

    \li \f$ -1 \leq z \leq 1 \f$

    \li Uses the defintion of the Lagrangian interpolant:\n
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
    */
    template<class Scalar>
    static Scalar hgrjp   (const int i, const Scalar z, const Scalar *zgrj,
                           const int np, const Scalar alpha, const Scalar beta);


    /** \brief Compute the value of the \a i th Lagrangian interpolant through the
               \a np Gauss-Lobatto-Jacobi points \a zglj at the arbitrary location
               \a z.

    \li \f$ -1 \leq z \leq 1 \f$

    \li Uses the defintion of the Lagrangian interpolant:\n
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
    template<class Scalar>
    static Scalar hglj    (const int i, const Scalar z, const Scalar *zglj,
                           const int np, const Scalar alpha, const Scalar beta);



    /* Interpolation operators */

    /** \brief Interpolation Operator from Gauss-Jacobi points to an
        arbitrary distribution at points \a zm

    \li Computes the one-dimensional interpolation matrix, \a im, to
    interpolate a function from at Gauss-Jacobi distribution of \a nz
    zeros \a zgj to an arbitrary distribution of \a mz points \a zm, i.e.\n
    \f$
    u(zm[i]) = \sum_{j=0}^{nz-1} im[i*nz+j] \ u(zgj[j])
    \f$
    */
    template<class Scalar>
    static void  Imgj  (Scalar *im, const Scalar *zgj, const Scalar *zm, const int nz,
                        const int mz, const Scalar alpha, const Scalar beta);


    /** \brief Interpolation Operator from Gauss-Radau-Jacobi points
               (including \a z=-1) to an arbitrary distrubtion at points \a zm

    \li Computes the one-dimensional interpolation matrix, \a im, to
    interpolate a function from at Gauss-Radau-Jacobi distribution of
    \a nz zeros \a zgrj (where \a zgrj[0]=-1) to an arbitrary
    distribution of \a mz points \a zm, i.e.
    \n
    \f$ u(zm[i]) =    \sum_{j=0}^{nz-1} im[i*nz+j] \ u(zgrj[j]) \f$
    */
    template<class Scalar>
    static void  Imgrjm(Scalar *im, const Scalar *zgrj, const Scalar *zm, const int nz,
                        const int mz, const Scalar alpha, const Scalar beta);


    /** \brief Interpolation Operator from Gauss-Radau-Jacobi points
               (including \a z=1) to an arbitrary distrubtion at points \a zm

    \li Computes the one-dimensional interpolation matrix, \a im, to
    interpolate a function from at Gauss-Radau-Jacobi distribution of
    \a nz zeros \a zgrj (where \a zgrj[nz-1]=1) to an arbitrary
    distribution of \a mz points \a zm, i.e.
    \n
    \f$ u(zm[i]) =    \sum_{j=0}^{nz-1} im[i*nz+j] \ u(zgrj[j]) \f$
    */
    template<class Scalar>
    static void  Imgrjp(Scalar *im, const Scalar *zgrj, const Scalar *zm, const int nz,
                        const int mz, const Scalar alpha, const Scalar beta);


    /** \brief Interpolation Operator from Gauss-Lobatto-Jacobi points
               to an arbitrary distrubtion at points \a zm

    \li Computes the one-dimensional interpolation matrix, \a im, to
    interpolate a function from at Gauss-Lobatto-Jacobi distribution of
    \a nz zeros \a zglj (where \a zglj[0]=-1 , \a  zglj[nz-1]=1) to an arbitrary
    distribution of \a mz points \a zm, i.e.
    \n
    \f$ u(zm[i]) =    \sum_{j=0}^{nz-1} im[i*nz+j] \ u(zglj[j]) \f$
    */
    template<class Scalar>
    static void  Imglj (Scalar *im, const Scalar *zglj, const Scalar *zm, const int nz,
                        const int mz, const Scalar alpha, const Scalar beta);


    /* Polynomial functions */

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
    template<class Scalar>
    static void jacobfd (const int np, const Scalar *z, Scalar *poly_in, Scalar *polyd,
                         const int n, const Scalar alpha, const Scalar beta);


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
    template<class Scalar>
    static void jacobd  (const int np, const Scalar *z, Scalar *polyd, const int n,
                         const Scalar alpha, const Scalar beta);



    /* Helper functions. */

    /** \brief  Calculate the \a n zeros, \a z, of the Jacobi polynomial, i.e.
                \f$ P_n^{\alpha,\beta}(z) = 0 \f$

    This routine is only valid for \f$( \alpha > -1, \beta > -1)\f$
    and uses polynomial deflation in a Newton iteration
    */
    template<class Scalar>
    static void   Jacobz (const int n, Scalar *z, const Scalar alpha, const Scalar beta);


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
    template<class Scalar>
    static void   JacZeros (const int n, Scalar *a, const Scalar alpha, const Scalar beta);


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
    template<class Scalar>
    static void   TriQL    (const int n, Scalar *d, Scalar *e);


    /** \brief Calculate the Gamma function , \f$ \Gamma(x)\f$, for integer
               values \a x and halves.

    Determine the value of \f$\Gamma(x)\f$ using:

    \f$ \Gamma(x) = (x-1)!  \mbox{ or  }  \Gamma(x+1/2) = (x-1/2)\Gamma(x-1/2)\f$

    where \f$ \Gamma(1/2) = \sqrt{\pi}\f$
    */
    template<class Scalar>
    static Scalar gammaF (const Scalar x);


  }; // class IntrepidPolylib

} // end of Intrepid namespace

// include templated definitions
#include <Intrepid_PolylibDef.hpp>

#endif

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

