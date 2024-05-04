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

/** \file   Intrepid_HGRAD_LINE_Cn_FEM_JACOBI.hpp
    \brief  Header file for the Intrepid::HGRAD_LINE_Cn_FEM_JACOBI class.
    \author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_HGRAD_LINE_Cn_FEM_JACOBI_HPP
#define INTREPID_HGRAD_LINE_Cn_FEM_JACOBI_HPP

#include "Intrepid_Basis.hpp"
#include "Intrepid_Polylib.hpp"

namespace Intrepid {
  
/** \class  Intrepid::Basis_HGRAD_LINE_Cn_FEM_JACOBI
    \brief  Implementation of the locally H(grad)-compatible FEM basis of variable order
            on the [-1,1] reference line cell, using Jacobi polynomials. 
  
            Implements Jacobi basis of variable order \f$n\f$ on
            the reference [-1,1] line cell. Jacobi polynomials depend on three parameters
            \f$ \alpha \f$, \f$ \beta \f$, and \f$ n \f$ and are defined via the so-called
            Gamma function by
            \f[
              P_n^{(\alpha,\beta)} (z) = 
                \frac{\Gamma (\alpha+n+1)}{n!\Gamma (\alpha+\beta+n+1)}
                \sum_{m=0}^n {n\choose m}
                \frac{\Gamma (\alpha + \beta + n + m + 1)}{\Gamma (\alpha + m + 1)} \left(\frac{z-1}{2}\right)^m
            \f]
            The basis has cardinality \f$n+1\f$ and spans a COMPLETE linear polynomial space.
            Basis functions are dual to a unisolvent set of degrees of freedom (DoF) enumerated as follows:

  <table>
    <tr>
    <th rowspan="2"> Basis order <th colspan="4"> DoF tag table <th rowspan="2"> DoF definition
    </tr>
    <tr>
    <th> subc dim <th> subc ordinal <th> subc DoF tag <th> subc num DoFs
    </tr>

    <tr align="center"> <td> 0 <td> 1 <td> 0 <td> 0   <td> 1
      <td align="left"> \f$ P_0^{(\alpha,\beta)} \f$ </tr>
    <tr align="center"> <td> 1 <td> 1 <td> 0 <td> 0-1 <td> 2
      <td align="left"> \f$ P_0^{(\alpha,\beta)}, P_1^{(\alpha,\beta)} \f$ </tr>
    <tr align="center"> <td> 2 <td> 1 <td> 0 <td> 0-2 <td> 3
      <td align="left"> \f$ P_0^{(\alpha,\beta)}, P_1^{(\alpha,\beta)}, P_2^{(\alpha,\beta)} \f$ </tr>
    <tr align="center"> <td> 3 <td> 1 <td> 0 <td> 0-3 <td> 4
      <td align="left"> \f$ P_0^{(\alpha,\beta)}, P_1^{(\alpha,\beta)}, ..., P_3^{(\alpha,\beta)} \f$ </tr>
    <tr align="center"> <td> ... <td> 1 <td> 0 <td> ... <td> ...
      <td align="left"> ... </tr>
    <tr align="center"> <td> n <td> 1 <td> 0 <td> 0-n <td> n+1
      <td align="left"> \f$ P_0^{(\alpha,\beta)}, P_1^{(\alpha,\beta)}, ..., P_n^{(\alpha,\beta)} \f$ </tr>
  </table>
  
            For example, for Legendre polynomials (\f$\alpha=\beta=0\f$), the first 11 bases are given by

  <table>
    <tr>
    <th rowspan="2"> Basis order <th colspan="4"> DoF tag table <th rowspan="2"> DoF definition
    </tr>
    <tr>
    <th> subc dim <th> subc ordinal <th> subc DoF tag <th> subc num DoFs
    </tr>

    <tr align="center"> <td> 0 <td> 1 <td> 0 <td> 0   <td> 1
      <td align="left"> \f$ 1 \f$ </tr>
    <tr align="center"> <td> 1 <td> 1 <td> 0 <td> 0-1 <td> 2
      <td align="char" char=":"> and: \f$ x \f$ </tr>
    <tr align="center"> <td> 2 <td> 1 <td> 0 <td> 0-2 <td> 3
      <td align="char" char=":"> and: \f$ \frac{1}{2} (3x^2-1) \f$ </tr>
    <tr align="center"> <td> 3 <td> 1 <td> 0 <td> 0-3 <td> 4
      <td align="char" char=":"> and: \f$ \frac{1}{2} (5x^3-3x) \f$ </tr>
    <tr align="center"> <td> 4 <td> 1 <td> 0 <td> 0-4 <td> 5
      <td align="char" char=":"> and: \f$ \frac{1}{8} (35x^4-30x^2+3) \f$ </tr>
    <tr align="center"> <td> 5 <td> 1 <td> 0 <td> 0-5 <td> 6
      <td align="char" char=":"> and: \f$ \frac{1}{8} (63x^5-70x^3+15x) \f$ </tr>
    <tr align="center"> <td> 6 <td> 1 <td> 0 <td> 0-6 <td> 7
      <td align="char" char=":"> and: \f$ \frac{1}{16} (231x^6-315x^4+105x^2-5) \f$ </tr>
    <tr align="center"> <td> 7 <td> 1 <td> 0 <td> 0-7 <td> 8
      <td align="char" char=":"> and: \f$ \frac{1}{16} (429x^7-693x^5+315x^3-35x) \f$ </tr>
    <tr align="center"> <td> 8 <td> 1 <td> 0 <td> 0-8 <td> 9
      <td align="char" char=":"> and: \f$ \frac{1}{128} (6435x^8-12012x^6+6930x^4-1260x^2+35) \f$ </tr>
    <tr align="center"> <td> 9 <td> 1 <td> 0 <td> 0-9 <td> 10
      <td align="char" char=":"> and: \f$ \frac{1}{128} (12155x^9-25740x^7+18018x^5-4620x^3+315x) \f$ </tr>
    <tr align="center"> <td>10 <td> 1 <td> 0 <td> 0-10<td> 11
      <td align="char" char=":"> and: \f$ \frac{1}{128} (46189x^{10}-109395x^8+90090x^6-30030x^4+3465x^2-63) \f$ </tr>
  </table>
*/
template<class Scalar, class ArrayScalar> 
class Basis_HGRAD_LINE_Cn_FEM_JACOBI: public Basis<Scalar, ArrayScalar> {
private:

  Scalar jacobiAlpha_;

  Scalar jacobiBeta_;
  
  /** \brief  Initializes <var>tagToOrdinal_</var> and <var>ordinalToTag_</var> lookup arrays.
  */
  void initializeTags();
  
public:
  
  /** \brief  Constructor.
  */
  Basis_HGRAD_LINE_Cn_FEM_JACOBI(int order, Scalar alpha = 0, Scalar beta = 0);  
  
  
  /** \brief  Evaluation of a FEM basis on a <strong>reference Line</strong> cell. 
  
              Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
              points in the <strong>reference Line</strong> cell. For rank and dimensions of
              I/O array arguments see Section \ref basis_md_array_sec .

      \param  outputValues      [out] - variable rank array with the basis values
      \param  inputPoints       [in]  - rank-2 array (P,D) with the evaluation points
      \param  operatorType      [in]  - the operator acting on the basis functions    
  */
  void getValues(ArrayScalar &          outputValues,
                 const ArrayScalar &    inputPoints,
                 const EOperator        operatorType) const;
  
  
  /**  \brief  FVD basis evaluation: invocation of this method throws an exception.
  */
  void getValues(ArrayScalar &          outputValues,
                 const ArrayScalar &    inputPoints,
                 const ArrayScalar &    cellVertices,
                 const EOperator        operatorType = OPERATOR_VALUE) const;

  /** \brief Sets private data \a <b>basisDegree_</b>, \a <b>basisCardinality_</b>,
             \a <b>jacobiAlpha_</b>, and \a <b>jacobiBeta_</b>, to
             \a <b>n</b>, \a <b>n+1</b>, \a <b>alpha</b>, and \a <b>beta</b>, respectively.
  */
  void setBasisParameters(int n, Scalar alpha = 0, Scalar beta = 0);
};

}// namespace Intrepid

#include "Intrepid_HGRAD_LINE_Cn_FEM_JACOBIDef.hpp"

#endif

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

