// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_HGRAD_WEDGE_C2_FEMDef.hpp
    \brief  Definition file for bi-quadratic FEM basis functions for H(grad) functions on WEDGE cells.
    \author Created by P. Bochev and D. Ridzal.
    Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_HGRAD_WEDGE_C2_FEM_DEF_HPP__
#define __INTREPID2_HGRAD_WEDGE_C2_FEM_DEF_HPP__

namespace Intrepid2 {

  // -------------------------------------------------------------------------------------
  template<typename SpT, typename OT, typename PT>
  template<EOperator opType>
  template<typename outputValueValueType, class ...outputValueProperties,
           typename inputPointValueType,  class ...inputPointProperties>
  KOKKOS_INLINE_FUNCTION
  void
  Basis_HGRAD_WEDGE_C2_FEM<SpT,OT,PT>::Serial<opType>::
  getValues( /**/  Kokkos::DynRankView<outputValueValueType,outputValueProperties...> output,
             const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  input ) {
    switch (opType) {
    case OPERATOR_VALUE: {
      const auto x = input(0);
      const auto y = input(1);
      const auto z = input(2);

      // output is a rank-2 array with dimensions (basisCardinality_, dim0)
      output(0) =  ((-1. + x + y)*(-1. + 2.*x + 2.*y)*(-1. + z)*z)/2.;
      output(1) =  (x*(-1. + 2.*x)*(-1. + z)*z)/2.;
      output(2) =  (y*(-1. + 2.*y)*(-1. + z)*z)/2.;
      output(3) =  ((-1. + x + y)*(-1. + 2.*x + 2.*y)*z*(1. + z))/2.;
      output(4) =  (x*(-1. + 2.*x)*z*(1. + z))/2.;
      output(5) =  (y*(-1. + 2.*y)*z*(1. + z))/2.;

      output(6) = -2.*x*(-1. + x + y)*(-1. + z)*z;
      output(7) =  2.*x*y*(-1. + z)*z;
      output(8) = -2.*y*(-1. + x + y)*(-1. + z)*z;
      output(9) = -((-1. + x + y)*(-1. + 2.*x + 2.*y)*(-1. + z)*(1. + z));
      output(10) = -(x*(-1. + 2.*x)*(-1. + z)*(1. + z));
      output(11) = -(y*(-1. + 2.*y)*(-1. + z)*(1. + z));
      output(12) = -2.*x*(-1. + x + y)*z*(1. + z);
      output(13) =  2.*x*y*z*(1. + z);
      output(14) = -2.*y*(-1. + x + y)*z*(1. + z);
      output(15) =  4.*x*(-1. + x + y)*(-1. + z)*(1. + z);
      output(16) = -4.*x*y*(-1. + z)*(1. + z);
      output(17) =  4.*y*(-1. + x + y)*(-1. + z)*(1. + z);
      break;
    }
    case OPERATOR_GRAD: {
      const auto x = input(0);
      const auto y = input(1);
      const auto z = input(2);

      // output is a rank-3 array with dimensions (basisCardinality_, dim0, spaceDim)
      output(0, 0) = ((-3 + 4*x + 4*y)*(-1 + z)*z)/2.;
      output(0, 1) = ((-3 + 4*x + 4*y)*(-1 + z)*z)/2.;
      output(0, 2) = ((-1 + x + y)*(-1 + 2*x + 2*y)*(-1 + 2*z))/2.;

      output(1, 0) = ((-1 + 4*x)*(-1 + z)*z)/2.;
      output(1, 1) = 0.;
      output(1, 2) = (x*(-1 + 2*x)*(-1 + 2*z))/2.;

      output(2, 0) = 0.;
      output(2, 1) = ((-1 + 4*y)*(-1 + z)*z)/2.;
      output(2, 2) = (y*(-1 + 2*y)*(-1 + 2*z))/2.;

      output(3, 0) = ((-3 + 4*x + 4*y)*z*(1 + z))/2.;
      output(3, 1) = ((-3 + 4*x + 4*y)*z*(1 + z))/2.;
      output(3, 2) = ((-1 + x + y)*(-1 + 2*x + 2*y)*(1 + 2*z))/2.;

      output(4, 0) = ((-1 + 4*x)*z*(1 + z))/2.;
      output(4, 1) = 0.;
      output(4, 2) = (x*(-1 + 2*x)*(1 + 2*z))/2.;

      output(5, 0) = 0.;
      output(5, 1) = ((-1 + 4*y)*z*(1 + z))/2.;
      output(5, 2) = (y*(-1 + 2*y)*(1 + 2*z))/2.;

      output(6, 0) = -2*(-1 + 2*x + y)*(-1 + z)*z;
      output(6, 1) = -2*x*(-1 + z)*z;
      output(6, 2) = 2*x*(-1 + x + y)*(1 - 2*z);

      output(7, 0) = 2*y*(-1 + z)*z;
      output(7, 1) = 2*x*(-1 + z)*z;
      output(7, 2) = 2*x*y*(-1 + 2*z);

      output(8, 0) = -2*y*(-1 + z)*z;
      output(8, 1) = -2*(-1 + x + 2*y)*(-1 + z)*z;
      output(8, 2) = 2*y*(-1 + x + y)*(1 - 2*z);

      output(9, 0) = -(-3 + 4*x + 4*y)*(-1 + z*z);
      output(9, 1) = -(-3 + 4*x + 4*y)*(-1 + z*z);
      output(9, 2) = -2*(1 + 2*x*x - 3*y + 2*y*y + x*(-3 + 4*y))*z;

      output(10, 0) = -(-1 + 4*x)*(-1 + z*z);
      output(10, 1) =  0;
      output(10, 2) =  2*(1 - 2*x)*x*z;

      output(11, 0) =  0;
      output(11, 1) =  -(-1 + 4*y)*(-1 + z*z);
      output(11, 2) =  2*(1 - 2*y)*y*z;

      output(12, 0) = -2*(-1 + 2*x + y)*z*(1 + z);
      output(12, 1) = -2*x*z*(1 + z);
      output(12, 2) = -2*x*(-1 + x + y)*(1 + 2*z);

      output(13, 0) =  2*y*z*(1 + z);
      output(13, 1) =  2*x*z*(1 + z);
      output(13, 2) =  2*x*y*(1 + 2*z);

      output(14, 0) = -2*y*z*(1 + z);
      output(14, 1) = -2*(-1 + x + 2*y)*z*(1 + z);
      output(14, 2) = -2*y*(-1 + x + y)*(1 + 2*z);

      output(15, 0) =  4*(-1 + 2*x + y)*(-1 + z*z);
      output(15, 1) =  4*x*(-1 + z)*(1 + z);
      output(15, 2) =  8*x*(-1 + x + y)*z;

      output(16, 0) = -4*y*(-1 + z)*(1 + z);
      output(16, 1) = -4*x*(-1 + z)*(1 + z);
      output(16, 2) = -8*x*y*z;

      output(17, 0) =  4*y*(-1 + z)*(1 + z);
      output(17, 1) =  4*(-1 + x + 2*y)*(-1 + z*z);
      output(17, 2) =  8*y*(-1 + x + y)*z;
      break;
    }
    case OPERATOR_D2: {
      const auto x = input(0);
      const auto y = input(1);
      const auto z = input(2);

      output(0, 0) =  2.*(-1. + z)*z;
      output(0, 1) =  2.*(-1. + z)*z;
      output(0, 2) =  ((-3. + 4.*x + 4.*y)*(-1. + 2.*z))/2.;
      output(0, 3) =  2.*(-1. + z)*z;
      output(0, 4) =  ((-3. + 4.*x + 4.*y)*(-1. + 2.*z))/2.;
      output(0, 5) =  (-1. + x + y)*(-1. + 2.*x + 2.*y);

      output(1, 0) =  2.*(-1. + z)*z;
      output(1, 1) =  0.;
      output(1, 2) =  ((-1. + 4.*x)*(-1. + 2.*z))/2.;
      output(1, 3) =  0.;
      output(1, 4) =  0.;
      output(1, 5) =  x*(-1. + 2.*x);

      output(2, 0) =  0.;
      output(2, 1) =  0.;
      output(2, 2) =  0.;
      output(2, 3) =  2.*(-1. + z)*z;
      output(2, 4) =  ((-1. + 4.*y)*(-1. + 2.*z))/2.;
      output(2, 5) =  y*(-1. + 2.*y);

      output(3, 0) =  2.*z*(1. + z);
      output(3, 1) =  2.*z*(1. + z);
      output(3, 2) =  ((-3. + 4.*x + 4.*y)*(1. + 2.*z))/2.;
      output(3, 3) =  2.*z*(1. + z);
      output(3, 4) =  ((-3. + 4.*x + 4.*y)*(1. + 2.*z))/2.;
      output(3, 5) =  (-1. + x + y)*(-1. + 2.*x + 2.*y);

      output(4, 0) =  2.*z*(1. + z);
      output(4, 1) =  0.;
      output(4, 2) =  ((-1. + 4.*x)*(1. + 2.*z))/2.;;
      output(4, 3) =  0.;
      output(4, 4) =  0.;
      output(4, 5) =  x*(-1. + 2.*x);

      output(5, 0) =  0.;
      output(5, 1) =  0.;
      output(5, 2) =  0.;
      output(5, 3) =  2.*z*(1. + z);
      output(5, 4) =  ((-1. + 4.*y)*(1. + 2.*z))/2.;
      output(5, 5) =  y*(-1. + 2.*y);

      output(6, 0) = -4.*(-1. + z)*z;
      output(6, 1) = -2.*(-1. + z)*z;
      output(6, 2) = -2.*(-1. + 2.*x + y)*(-1. + 2.*z);
      output(6, 3) =  0.;
      output(6, 4) =  x*(2. - 4.*z);
      output(6, 5) = -4.*x*(-1. + x + y);

      output(7, 0) =  0.;
      output(7, 1) =  2.*(-1. + z)*z;
      output(7, 2) =  2.*y*(-1. + 2.*z);
      output(7, 3) =  0.;
      output(7, 4) =  2.*x*(-1. + 2.*z);
      output(7, 5) =  4.*x*y;

      output(8, 0) =  0.;
      output(8, 1) = -2.*(-1. + z)*z;
      output(8, 2) =  y*(2. - 4.*z);
      output(8, 3) = -4.*(-1. + z)*z;
      output(8, 4) = -2.*(-1. + x + 2.*y)*(-1. + 2.*z);
      output(8, 5) = -4.*y*(-1. + x + y);

      output(9, 0) =  4. - 4.*z*z;
      output(9, 1) =  4. - 4.*z*z;
      output(9, 2) = -2.*(-3. + 4.*x + 4.*y)*z;
      output(9, 3) =  4. - 4.*z*z;
      output(9, 4) = -2.*(-3. + 4.*x + 4.*y)*z;
      output(9, 5) = -2.*(-1. + x + y)*(-1. + 2.*x + 2.*y);

      output(10, 0) =  4. - 4.*z*z;
      output(10, 1) =  0.;
      output(10, 2) =  (2. - 8.*x)*z;
      output(10, 3) =  0.;
      output(10, 4) =  0.;
      output(10, 5) = -2.*x*(-1. + 2.*x);

      output(11, 0) =  0.;
      output(11, 1) =  0.;
      output(11, 2) =  0.;
      output(11, 3) =  4. - 4.*z*z;
      output(11, 4) =  (2. - 8.*y)*z;
      output(11, 5) = -2.*y*(-1. + 2.*y);

      output(12, 0) = -4.*z*(1. + z);
      output(12, 1) = -2.*z*(1. + z);
      output(12, 2) = -2.*(-1. + 2.*x + y)*(1. + 2.*z);
      output(12, 3) =  0.;
      output(12, 4) = -2.*(x + 2.*x*z);
      output(12, 5) = -4.*x*(-1. + x + y);

      output(13, 0) =  0.;
      output(13, 1) =  2.*z*(1. + z);
      output(13, 2) =  2.*(y + 2.*y*z);
      output(13, 3) =  0.;
      output(13, 4) =  2.*(x + 2.*x*z);
      output(13, 5) =  4.*x*y;

      output(14, 0) =  0.;
      output(14, 1) = -2.*z*(1. + z);
      output(14, 2) = -2.*(y + 2.*y*z);
      output(14, 3) = -4.*z*(1. + z);
      output(14, 4) = -2.*(-1. + x + 2.*y)*(1. + 2.*z);
      output(14, 5) = -4.*y*(-1. + x + y);

      output(15, 0) =  8.*(-1. + z*z);
      output(15, 1) =  4.*(-1. + z*z);
      output(15, 2) =  8.*(-1. + 2.*x + y)*z;
      output(15, 3) =  0.;
      output(15, 4) =  8.*x*z;
      output(15, 5) =  8.*x*(-1. + x + y);

      output(16, 0) =  0.;
      output(16, 1) =  4. - 4.*z*z;
      output(16, 2) = -8.*y*z;
      output(16, 3) =  0.;
      output(16, 4) = -8.*x*z;
      output(16, 5) = -8.*x*y;


      output(17, 0) =  0.;
      output(17, 1) =  4.*(-1. + z*z);
      output(17, 2) =  8.*y*z;
      output(17, 3) =  8.*(-1. + z*z);
      output(17, 4) =  8.*(-1. + x + 2.*y)*z;
      output(17, 5) =  8.*y*(-1. + x + y);
      break;
    }
    case OPERATOR_D3: {
      const auto x = input(0);
      const auto y = input(1);
      const auto z = input(2);

      output(0, 0) =  0.;
      output(0, 1) =  0.;
      output(0, 2) = -2. + 4.*z;
      output(0, 3) =  0.;
      output(0, 4) = -2. + 4.*z;
      output(0, 5) = -3. + 4.*x + 4.*y;
      output(0, 6) =  0.;
      output(0, 7) = -2. + 4.*z;
      output(0, 8) = -3. + 4.*x + 4.*y;
      output(0, 9) =  0.;

      output(1, 0) =  0.;
      output(1, 1) =  0.;
      output(1, 2) = -2. + 4.*z;
      output(1, 3) =  0.;
      output(1, 4) =  0.;
      output(1, 5) = -1 + 4.*x;
      output(1, 6) =  0.;
      output(1, 7) =  0.;
      output(1, 8) =  0.;
      output(1, 9) =  0.;

      output(2, 0) =  0.;
      output(2, 1) =  0.;
      output(2, 2) =  0.;
      output(2, 3) =  0.;
      output(2, 4) =  0.;
      output(2, 5) =  0.;
      output(2, 6) =  0.;
      output(2, 7) = -2. + 4.*z;
      output(2, 8) = -1 + 4.*y;
      output(2, 9) =  0.;

      output(3, 0) =  0.;
      output(3, 1) =  0.;
      output(3, 2) =  2. + 4.*z;
      output(3, 3) =  0.;
      output(3, 4) =  2. + 4.*z;
      output(3, 5) = -3. + 4.*x + 4.*y;
      output(3, 6) =  0.;
      output(3, 7) =  2. + 4.*z;
      output(3, 8) = -3. + 4.*x + 4.*y;
      output(3, 9) =  0.;

      output(4, 0) =  0.;
      output(4, 1) =  0.;
      output(4, 2) =  2. + 4.*z;
      output(4, 3) =  0.;
      output(4, 4) =  0.;
      output(4, 5) = -1 + 4.*x;
      output(4, 6) =  0.;
      output(4, 7) =  0.;
      output(4, 8) =  0.;
      output(4, 9) =  0.;

      output(5, 0) =  0.;
      output(5, 1) =  0.;
      output(5, 2) =  0.;
      output(5, 3) =  0.;
      output(5, 4) =  0.;
      output(5, 5) =  0.;
      output(5, 6) =  0.;
      output(5, 7) =  2. + 4.*z;
      output(5, 8) = -1 + 4.*y;
      output(5, 9) =  0.;

      output(6, 0) =  0.;
      output(6, 1) =  0.;
      output(6, 2) =  4. - 8.*z;
      output(6, 3) =  0.;
      output(6, 4) =  2. - 4.*z;
      output(6, 5) = -4.*(-1 + 2*x + y);
      output(6, 6) =  0.;
      output(6, 7) =  0.;
      output(6, 8) = -4.*x;
      output(6, 9) =  0.;

      output(7, 0) =  0.;
      output(7, 1) =  0.;
      output(7, 2) =  0.;
      output(7, 3) =  0.;
      output(7, 4) = -2. + 4.*z;
      output(7, 5) =  4.*y;
      output(7, 6) =  0.;
      output(7, 7) =  0.;
      output(7, 8) =  4.*x;
      output(7, 9) =  0.;

      output(8, 0) =  0.;
      output(8, 1) =  0.;
      output(8, 2) =  0.;
      output(8, 3) =  0.;
      output(8, 4) =  2. - 4.*z;
      output(8, 5) = -4.*y;
      output(8, 6) =  0.;
      output(8, 7) =  4. - 8.*z;
      output(8, 8) = -4.*(-1 + x + 2*y);
      output(8, 9) =  0.;

      output(9, 0) =  0.;
      output(9, 1) =  0.;
      output(9, 2) = -8.*z;
      output(9, 3) =  0.;
      output(9, 4) = -8.*z;
      output(9, 5) =  6. - 8.*x - 8.*y;
      output(9, 6) =  0.;
      output(9, 7) = -8.*z;
      output(9, 8) =  6. - 8.*x - 8.*y;
      output(9, 9) =  0.;

      output(10, 0) =  0.;
      output(10, 1) =  0.;
      output(10, 2) = -8.*z;
      output(10, 3) =  0.;
      output(10, 4) =  0.;
      output(10, 5) =  2. - 8.*x;
      output(10, 6) =  0.;
      output(10, 7) =  0.;
      output(10, 8) =  0.;
      output(10, 9) =  0.;

      output(11, 0) =  0.;
      output(11, 1) =  0.;
      output(11, 2) =  0.;
      output(11, 3) =  0.;
      output(11, 4) =  0.;
      output(11, 5) =  0.;
      output(11, 6) =  0.;
      output(11, 7) = -8.*z;
      output(11, 8) =  2. - 8.*y;
      output(11, 9) =  0.;

      output(12, 0) =  0.;
      output(12, 1) =  0.;
      output(12, 2) = -4. - 8.*z;
      output(12, 3) =  0.;
      output(12, 4) = -2. - 4.*z;
      output(12, 5) = -4.*(-1 + 2*x + y);
      output(12, 6) =  0.;
      output(12, 7) =  0.;
      output(12, 8) = -4.*x;
      output(12, 9) =  0.;

      output(13, 0) =  0.;
      output(13, 1) =  0.;
      output(13, 2) =  0.;
      output(13, 3) =  0.;
      output(13, 4) =  2. + 4.*z;
      output(13, 5) =  4.*y;
      output(13, 6) =  0.;
      output(13, 7) =  0.;
      output(13, 8) =  4.*x;
      output(13, 9) =  0.;

      output(14, 0) =  0.;
      output(14, 1) =  0.;
      output(14, 2) =  0.;
      output(14, 3) =  0.;
      output(14, 4) = -2. - 4.*z;
      output(14, 5) = -4.*y;
      output(14, 6) =  0.;
      output(14, 7) = -4. - 8.*z;
      output(14, 8) = -4.*(-1 + x + 2*y);
      output(14, 9) =  0.;

      output(15, 0) =  0.;
      output(15, 1) =  0.;
      output(15, 2) =  16.*z;
      output(15, 3) =  0.;
      output(15, 4) =  8.*z;
      output(15, 5) =  8.*(-1 + 2*x + y);
      output(15, 6) =  0.;
      output(15, 7) =  0.;
      output(15, 8) =  8.*x;
      output(15, 9) =  0.;

      output(16, 0) =  0.;
      output(16, 1) =  0.;
      output(16, 2) =  0.;
      output(16, 3) =  0.;
      output(16, 4) = -8.*z;
      output(16, 5) = -8.*y;
      output(16, 6) =  0.;
      output(16, 7) =  0.;
      output(16, 8) = -8.*x;
      output(16, 9) =  0.;

      output(17, 0) =  0.;
      output(17, 1) =  0.;
      output(17, 2) =  0.;
      output(17, 3) =  0.;
      output(17, 4) =  8.*z;
      output(17, 5) =  8.*y;
      output(17, 6) =  0.;
      output(17, 7) =  16.*z;
      output(17, 8) =  8.*(-1 + x + 2*y);
      output(17, 9) =  0.;
      break;
    }
    case OPERATOR_D4: {
      const auto jend = output.dimension(1);
      const auto iend = output.dimension(0);

      for (size_type j=0;j<jend;++j)
        for (size_type i=0;i<iend;++i)
          output(i, j) = 0.0;

      output(0, 5) = 4.;
      output(0, 8) = 4.;
      output(0,12) = 4.;

      output(1, 5) = 4.;

      output(2,12) = 4.;

      output(3, 5) = 4.;
      output(3, 8) = 4.;
      output(3,12) = 4.;

      output(4, 5) = 4.0;

      output(5,12) = 4.0;

      output(6, 5) =-8.;
      output(6, 8) =-4.;

      output(7, 8) = 4.;

      output(8, 8) =-4.;
      output(8,12) =-8.;

      output(9, 5) =-8.;
      output(9, 8) =-8.;
      output(9,12) =-8.;

      output(10, 5) =-8.;

      output(11,12) =-8.;

      output(12, 5) =-8.;
      output(12, 8) =-4.;

      output(13, 8) = 4.;

      output(14, 8) =-4;
      output(14,12) =-8.;

      output(15, 5) =16.;
      output(15, 8) = 8.;

      output(16, 8) =-8.;


      output(17, 8) = 8.;
      output(17,12) =16.;
      break;
    }
    case OPERATOR_MAX : {
      const auto jend = output.dimension(1);
      const auto iend = output.dimension(0);

      for (size_type j=0;j<jend;++j)
        for (size_type i=0;i<iend;++i)
          output(i, j) = 0.0;
      break;
    }

    default: {
      INTREPID2_TEST_FOR_ABORT( opType != OPERATOR_VALUE &&
                                opType != OPERATOR_GRAD &&
                                opType != OPERATOR_D2 &&
                                opType != OPERATOR_D3 &&
                                opType != OPERATOR_D4 &&
                                opType != OPERATOR_MAX,
                                ">>> ERROR: (Intrepid2::Basis_HGRAD_WEDGE_C2_FEM::Serial::getValues) operator is not supported");
    }
    }
  }

  // -------------------------------------------------------------------------------------

  template<typename SpT, typename OT, typename PT>
  Basis_HGRAD_WEDGE_C2_FEM<SpT,OT,PT>::
  Basis_HGRAD_WEDGE_C2_FEM()
    : impl_(this) {
    this->basisCardinality_  = 18;
    this->basisDegree_       = 2;
    this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Wedge<6> >() );
    this->basisType_         = BASIS_FEM_DEFAULT;
    this->basisCoordinates_  = COORDINATES_CARTESIAN;

    // initialize tags
    {
      // Basis-dependent intializations
      const ordinal_type tagSize  = 4;        // size of DoF tag
      const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
      const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
      const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

      // An array with local DoF tags assigned to basis functions, in the order of their local enumeration
      ordinal_type tags[72]  = { 0, 0, 0, 1,
                                 0, 1, 0, 1,
                                 0, 2, 0, 1,
                                 0, 3, 0, 1,
                                 0, 4, 0, 1,
                                 0, 5, 0, 1,
                                 1, 0, 0, 1,
                                 1, 1, 0, 1,
                                 1, 2, 0, 1,
                                 1, 6, 0, 1,
                                 1, 7, 0, 1,
                                 1, 8, 0, 1,
                                 1, 3, 0, 1,
                                 1, 4, 0, 1,
                                 1, 5, 0, 1,
                                 2, 0, 0, 1,
                                 2, 1, 0, 1,
                                 2, 2, 0, 1
      };

      // host tags
      ordinal_type_array_1d_host tagView(&tags[0], 72);

      // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
      //ordinal_type_array_2d_host ordinalToTag;
      //ordinal_type_array_3d_host tagToOrdinal;
      this->setOrdinalTagData(this->tagToOrdinal_,
                              this->ordinalToTag_,
                              tagView,
                              this->basisCardinality_,
                              tagSize,
                              posScDim,
                              posScOrd,
                              posDfOrd);
    }

    // dofCoords on host and create its mirror view to device
    Kokkos::DynRankView<PT,typename SpT::array_layout,Kokkos::HostSpace>
      dofCoords("dofCoordsHost", this->basisCardinality_,this->basisCellTopology_.getDimension());

    dofCoords(0,0) =  0.0;  dofCoords(0,1) =  0.0;  dofCoords(0,2) = -1.0;
    dofCoords(1,0) =  1.0;  dofCoords(1,1) =  0.0;  dofCoords(1,2) = -1.0;
    dofCoords(2,0) =  0.0;  dofCoords(2,1) =  1.0;  dofCoords(2,2) = -1.0;
    dofCoords(3,0) =  0.0;  dofCoords(3,1) =  0.0;  dofCoords(3,2) =  1.0;
    dofCoords(4,0) =  1.0;  dofCoords(4,1) =  0.0;  dofCoords(4,2) =  1.0;
    dofCoords(5,0) =  0.0;  dofCoords(5,1) =  1.0;  dofCoords(5,2) =  1.0;

    dofCoords(6,0) =  0.5;  dofCoords(6,1) =  0.0;  dofCoords(6,2) = -1.0;
    dofCoords(7,0) =  0.5;  dofCoords(7,1) =  0.5;  dofCoords(7,2) = -1.0;
    dofCoords(8,0) =  0.0;  dofCoords(8,1) =  0.5;  dofCoords(8,2) = -1.0;
    dofCoords(9,0) =  0.0;  dofCoords(9,1) =  0.0;  dofCoords(9,2) =  0.0;
    dofCoords(10,0)=  1.0;  dofCoords(10,1)=  0.0;  dofCoords(10,2)=  0.0;
    dofCoords(11,0)=  0.0;  dofCoords(11,1)=  1.0;  dofCoords(11,2)=  0.0;

    dofCoords(12,0)=  0.5;  dofCoords(12,1)=  0.0;  dofCoords(12,2)=  1.0;
    dofCoords(13,0)=  0.5;  dofCoords(13,1)=  0.5;  dofCoords(13,2)=  1.0;
    dofCoords(14,0)=  0.0;  dofCoords(14,1)=  0.5;  dofCoords(14,2)=  1.0;
    dofCoords(15,0)=  0.5;  dofCoords(15,1)=  0.0;  dofCoords(15,2)=  0.0;
    dofCoords(16,0)=  0.5;  dofCoords(16,1)=  0.5;  dofCoords(16,2)=  0.0;
    dofCoords(17,0)=  0.0;  dofCoords(17,1)=  0.5;  dofCoords(17,2)=  0.0;

    this->dofCoords_ = Kokkos::create_mirror_view(typename SpT::memory_space(), dofCoords);
    Kokkos::deep_copy(this->dofCoords_, dofCoords);
  }

  template<typename SpT, typename OT, typename PT>
  template<typename outputValueValueType, class ...outputValueProperties,
           typename inputPointValueType,  class ...inputPointProperties>
  void
  Basis_HGRAD_WEDGE_C2_FEM<SpT,OT,PT>::Internal::
  getValues( /**/  Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
             const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
             const EOperator operatorType ) const {
#ifdef HAVE_INTREPID2_DEBUG
    Intrepid2::getValues_HGRAD_Args(outputValues,
                                    inputPoints,
                                    operatorType,
                                    obj_->getBaseCellTopology(),
                                    obj_->getCardinality() );
#endif

    typedef          Kokkos::DynRankView<outputValueValueType,outputValueProperties...>         outputValueViewType;
    typedef          Kokkos::DynRankView<inputPointValueType, inputPointProperties...>          inputPointViewType;
    typedef typename ExecSpace<typename inputPointViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;

    // Number of evaluation points = dim 0 of inputPoints
    const auto loopSize = inputPoints.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    switch (operatorType) {

    case OPERATOR_VALUE: {
      typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_VALUE> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
      break;
    }
    case OPERATOR_GRAD:
    case OPERATOR_D1: {
      typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_GRAD> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
      break;
    }
    case OPERATOR_CURL: {
      INTREPID2_TEST_FOR_EXCEPTION( operatorType == OPERATOR_CURL, std::invalid_argument,
                                    ">>> ERROR (Basis_HGRAD_WEDGE_C2_FEM): CURL is invalid operator for rank-0 (scalar) functions in 3D");
      break;
    }
    case OPERATOR_DIV: {
      INTREPID2_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_DIV), std::invalid_argument,
                                    ">>> ERROR (Basis_HGRAD_WEDGE_C2_FEM): DIV is invalid operator for rank-0 (scalar) functions in 3D");
      break;
    }
    case OPERATOR_D2: {
      typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_D2> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
      break;
    }
    case OPERATOR_D3: {
      typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_D3> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
      break;
    }
    case OPERATOR_D4: {
      typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_D3> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
      break;
    }
    case OPERATOR_D5:
    case OPERATOR_D6:
    case OPERATOR_D7:
    case OPERATOR_D8:
    case OPERATOR_D9:
    case OPERATOR_D10: {
      typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_MAX> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
      break;
    }
    default: {
      INTREPID2_TEST_FOR_EXCEPTION( !( Intrepid2::isValidOperator(operatorType) ), std::invalid_argument,
                                    ">>> ERROR (Basis_HGRAD_WEDGE_C2_FEM): Invalid operator type");
    }
    }
  }

  template<typename SpT, typename OT, typename PT>
  template<typename dofCoordValueType, class ...dofCoordProperties>
  void
  Basis_HGRAD_WEDGE_C2_FEM<SpT,OT,PT>::Internal::
  getDofCoords( Kokkos::DynRankView<dofCoordValueType,dofCoordProperties...> dofCoords ) const {
#ifdef HAVE_INTREPID2_DEBUG
    // Verify rank of output array.
    INTREPID2_TEST_FOR_EXCEPTION( dofCoords.rank() != 2, std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::Basis_HGRAD_WEDGE_C2_FEM::getDofCoords) rank = 2 required for dofCoords array");
    // Verify 0th dimension of output array.
    INTREPID2_TEST_FOR_EXCEPTION( dofCoords.dimension(0) != obj_->basisCardinality_, std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::Basis_HGRAD_WEDGE_C2_FEM::getDofCoords) mismatch in number of dof and 0th dimension of dofCoords array");
    // Verify 1st dimension of output array.
    INTREPID2_TEST_FOR_EXCEPTION( dofCoords.dimension(1) != obj_->basisCellTopology_.getDimension(), std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::Basis_HGRAD_WEDGE_C2_FEM::getDofCoords) incorrect reference cell (1st) dimension in dofCoords array");
#endif
    Kokkos::deep_copy(dofCoords, obj_->dofCoords_);
  }

}// namespace Intrepid2
#endif
