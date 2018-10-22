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

/** \file   Intrepid2_HGRAD_HEX_C2_FEMDef.hpp
    \brief  Definition file for FEM basis functions of degree 2 for H(grad) functions on HEX cells.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
 */

#ifndef __INTREPID2_HGRAD_HEX_C2_FEM_DEF_HPP__
#define __INTREPID2_HGRAD_HEX_C2_FEM_DEF_HPP__

namespace Intrepid2 {

  // -------------------------------------------------------------------------------------
  namespace Impl {
    
    template<EOperator opType>
    template<typename outputViewType,
             typename inputViewType>
    KOKKOS_INLINE_FUNCTION
    void
    Basis_HGRAD_HEX_C2_FEM::Serial<opType>::
    getValues(       outputViewType output,
               const inputViewType input ) {
      switch (opType) {
      case OPERATOR_VALUE : {
        const auto x = input(0);
        const auto y = input(1);
        const auto z = input(2);
        
        // output is a rank-2 array with dimensions (basisCardinality_, dim0)
        output.access( 0) = 0.125*(-1. + x)*x*(-1. + y)*y*(-1. + z)*z;
        output.access( 1) = 0.125*x*(1.+ x)*(-1. + y)*y*(-1. + z)*z;
        output.access( 2) = 0.125*x*(1.+ x)*y*(1.+ y)*(-1. + z)*z;
        output.access( 3) = 0.125*(-1. + x)*x*y*(1.+ y)*(-1. + z)*z;
        output.access( 4) = 0.125*(-1. + x)*x*(-1. + y)*y*z*(1.+ z);
        output.access( 5) = 0.125*x*(1.+ x)*(-1. + y)*y*z*(1.+ z);
        output.access( 6) = 0.125*x*(1.+ x)*y*(1.+ y)*z*(1.+ z);
        output.access( 7) = 0.125*(-1. + x)*x*y*(1.+ y)*z*(1.+ z);
        output.access( 8) = 0.25*(1. - x)*(1. + x)*(-1. + y)*y*(-1. + z)*z;
        output.access( 9) = 0.25*x*(1.+ x)*(1. - y)*(1. + y)*(-1. + z)*z;
        output.access(10) = 0.25*(1. - x)*(1. + x)*y*(1.+ y)*(-1. + z)*z;
        output.access(11) = 0.25*(-1. + x)*x*(1. - y)*(1. + y)*(-1. + z)*z;
        output.access(12) = 0.25*(-1. + x)*x*(-1. + y)*y*(1. - z)*(1. + z);
        output.access(13) = 0.25*x*(1.+ x)*(-1. + y)*y*(1. - z)*(1. + z);
        output.access(14) = 0.25*x*(1.+ x)*y*(1.+ y)*(1. - z)*(1. + z);
        output.access(15) = 0.25*(-1. + x)*x*y*(1.+ y)*(1. - z)*(1. + z);
        output.access(16) = 0.25*(1. - x)*(1. + x)*(-1. + y)*y*z*(1.+ z);
        output.access(17) = 0.25*x*(1.+ x)*(1. - y)*(1. + y)*z*(1.+ z);
        output.access(18) = 0.25*(1. - x)*(1. + x)*y*(1.+ y)*z*(1.+ z);
        output.access(19) = 0.25*(-1. + x)*x*(1. - y)*(1. + y)*z*(1.+ z);
        output.access(20) = (1. - x)*(1. + x)*(1. - y)*(1. + y)*(1. - z)*(1. + z);
        output.access(21) = 0.5*(1. - x)*(1. + x)*(1. - y)*(1. + y)*(-1. + z)*z;
        output.access(22) = 0.5*(1. - x)*(1. + x)*(1. - y)*(1. + y)*z*(1.+ z);
        output.access(23) = 0.5*(-1. + x)*x*(1. - y)*(1. + y)*(1. - z)*(1. + z);
        output.access(24) = 0.5*x*(1.+ x)*(1. - y)*(1. + y)*(1. - z)*(1. + z);
        output.access(25) = 0.5*(1. - x)*(1. + x)*(-1. + y)*y*(1. - z)*(1. + z);
        output.access(26) = 0.5*(1. - x)*(1. + x)*y*(1.+ y)*(1. - z)*(1. + z);
        break;
      }
      case OPERATOR_GRAD :
      case OPERATOR_D1 : {
        const auto x = input(0);
        const auto y = input(1);
        const auto z = input(2);

        // output.access is a rank-3 array with dimensions (basisCardinality_, dim0, spaceDim)
        output.access(0,  0) = (-0.125 + 0.25*x)*(-1. + y)*y*(-1. + z)*z;
        output.access(0,  1) = (-1. + x)*x*(-0.125 + 0.25*y)*(-1. + z)*z;
        output.access(0,  2) = (-1. + x)*x*(-1. + y)*y*(-0.125 + 0.25*z);

        output.access(1,  0) = (0.125 + 0.25*x)*(-1. + y)*y*(-1. + z)*z;
        output.access(1,  1) = x*(1. + x)*(-0.125 + 0.25*y)*(-1. + z)*z;
        output.access(1,  2) = x*(1. + x)*(-1. + y)*y*(-0.125 + 0.25*z);

        output.access(2,  0) = (0.125 + 0.25*x)*y*(1. + y)*(-1. + z)*z;
        output.access(2,  1) = x*(1. + x)*(0.125 + 0.25*y)*(-1. + z)*z;
        output.access(2,  2) = x*(1. + x)*y*(1. + y)*(-0.125 + 0.25*z);

        output.access(3,  0) = (-0.125 + 0.25*x)*y*(1. + y)*(-1. + z)*z;
        output.access(3,  1) = (-1. + x)*x*(0.125 + 0.25*y)*(-1. + z)*z;
        output.access(3,  2) = (-1. + x)*x*y*(1. + y)*(-0.125 + 0.25*z);

        output.access(4,  0) = (-0.125 + 0.25*x)*(-1. + y)*y*z*(1. + z);
        output.access(4,  1) = (-1. + x)*x*(-0.125 + 0.25*y)*z*(1. + z);
        output.access(4,  2) = (-1. + x)*x*(-1. + y)*y*(0.125 + 0.25*z);

        output.access(5,  0) = (0.125 + 0.25*x)*(-1. + y)*y*z*(1. + z);
        output.access(5,  1) = x*(1. + x)*(-0.125 + 0.25*y)*z*(1. + z);
        output.access(5,  2) = x*(1. + x)*(-1. + y)*y*(0.125 + 0.25*z);

        output.access(6,  0) = (0.125 + 0.25*x)*y*(1. + y)*z*(1. + z);
        output.access(6,  1) = x*(1. + x)*(0.125 + 0.25*y)*z*(1. + z);
        output.access(6,  2) = x*(1. + x)*y*(1. + y)*(0.125 + 0.25*z);

        output.access(7,  0) = (-0.125 + 0.25*x)*y*(1. + y)*z*(1. + z);
        output.access(7,  1) = (-1. + x)*x*(0.125 + 0.25*y)*z*(1. + z);
        output.access(7,  2) = (-1. + x)*x*y*(1. + y)*(0.125 + 0.25*z);

        output.access(8,  0) = -0.5*x*(-1. + y)*y*(-1. + z)*z;
        output.access(8,  1) = (1. - x)*(1. + x)*(-0.25 + 0.5*y)*(-1. + z)*z;
        output.access(8,  2) = (1. - x)*(1. + x)*(-1. + y)*y*(-0.25 + 0.5*z);

        output.access(9,  0) = (0.25 + 0.5*x)*(1. - y)*(1. + y)*(-1. + z)*z;
        output.access(9,  1) = x*(1. + x)*(-0.5*y)*(-1. + z)*z;
        output.access(9,  2) = x*(1. + x)*(1. - y)*(1. + y)*(-0.25 + 0.5*z);

        output.access(10, 0) = -0.5*x*y*(1. + y)*(-1. + z)*z;
        output.access(10, 1) = (1. - x)*(1. + x)*(0.25 + 0.5*y)*(-1. + z)*z;
        output.access(10, 2) = (1. - x)*(1. + x)*y*(1. + y)*(-0.25 + 0.5*z);

        output.access(11, 0) = (-0.25 + 0.5*x)*(1. - y)*(1. + y)*(-1. + z)*z;
        output.access(11, 1) = (-1. + x)*x*(-0.5*y)*(-1. + z)*z;
        output.access(11, 2) = (-1. + x)*x*(1. - y)*(1. + y)*(-0.25 + 0.5*z);

        output.access(12, 0) = (-0.25 + 0.5*x)*(-1. + y)*y*(1. - z)*(1. + z);
        output.access(12, 1) = (-1. + x)*x*(-0.25 + 0.5*y)*(1. - z)*(1. + z);
        output.access(12, 2) = (-1. + x)*x*(-1. + y)*y*(-0.5*z);

        output.access(13, 0) = (0.25 + 0.5*x)*(-1. + y)*y*(1. - z)*(1. + z);
        output.access(13, 1) = x*(1. + x)*(-0.25 + 0.5*y)*(1. - z)*(1. + z);
        output.access(13, 2) = x*(1. + x)*(-1. + y)*y*(-0.5*z);

        output.access(14, 0) = (0.25 + 0.5*x)*y*(1. + y)*(1. - z)*(1. + z);
        output.access(14, 1) = x*(1. + x)*(0.25 + 0.5*y)*(1. - z)*(1. + z);
        output.access(14, 2) = x*(1. + x)*y*(1. + y)*(-0.5*z);

        output.access(15, 0) = (-0.25 + 0.5*x)*y*(1. + y)*(1. - z)*(1. + z);
        output.access(15, 1) = (-1. + x)*x*(0.25 + 0.5*y)*(1. - z)*(1. + z);
        output.access(15, 2) = (-1. + x)*x*y*(1. + y)*(-0.5*z);

        output.access(16, 0) = -0.5*x*(-1. + y)*y*z*(1. + z);
        output.access(16, 1) = (1. - x)*(1. + x)*(-0.25 + 0.5*y)*z*(1. + z);
        output.access(16, 2) = (1. - x)*(1. + x)*(-1. + y)*y*(0.25 + 0.5*z);

        output.access(17, 0) = (0.25 + 0.5*x)*(1. - y)*(1. + y)*z*(1. + z);
        output.access(17, 1) = x*(1. + x)*(-0.5*y)*z*(1. + z);
        output.access(17, 2) = x*(1. + x)*(1. - y)*(1. + y)*(0.25 + 0.5*z);

        output.access(18, 0) = -0.5*x*y*(1. + y)*z*(1. + z);
        output.access(18, 1) = (1. - x)*(1. + x)*(0.25 + 0.5*y)*z*(1. + z);
        output.access(18, 2) = (1. - x)*(1. + x)*y*(1. + y)*(0.25 + 0.5*z);

        output.access(19, 0) = (-0.25 + 0.5*x)*(1. - y)*(1. + y)*z*(1. + z);
        output.access(19, 1) = (-1. + x)*x*(-0.5*y)*z*(1. + z);
        output.access(19, 2) = (-1. + x)*x*(1. - y)*(1. + y)*(0.25 + 0.5*z);

        output.access(20, 0) = -2.*x*(1. - y)*(1. + y)*(1. - z)*(1. + z);
        output.access(20, 1) = (1. - x)*(1. + x)*(-2.*y)*(1. - z)*(1. + z);
        output.access(20, 2) = (1. - x)*(1. + x)*(1. - y)*(1. + y)*(-2.*z);

        output.access(21, 0) = -x*(1. - y)*(1. + y)*(-1. + z)*z;
        output.access(21, 1) = (1. - x)*(1. + x)*(-y)*(-1. + z)*z;
        output.access(21, 2) = (1. - x)*(1. + x)*(1. - y)*(1. + y)*(-0.5 + z);

        output.access(22, 0) = -x*(1. - y)*(1. + y)*z*(1. + z);
        output.access(22, 1) = (1. - x)*(1. + x)*(-y)*z*(1. + z);
        output.access(22, 2) = (1. - x)*(1. + x)*(1. - y)*(1. + y)*(0.5 + z);

        output.access(23, 0) = (-0.5 + x)*(1. - y)*(1. + y)*(1. - z)*(1. + z);
        output.access(23, 1) = (-1. + x)*x*(-y)*(1. - z)*(1. + z);
        output.access(23, 2) = (-1. + x)*x*(1. - y)*(1. + y)*(-z);

        output.access(24, 0) = (0.5 + x)*(1. - y)*(1. + y)*(1. - z)*(1. + z);
        output.access(24, 1) = x*(1. + x)*(-y)*(1. - z)*(1. + z);
        output.access(24, 2) = x*(1. + x)*(1. - y)*(1. + y)*(-z);

        output.access(25, 0) = -x*(-1. + y)*y*(1. - z)*(1. + z);
        output.access(25, 1) = (1. - x)*(1. + x)*(-0.5 + y)*(1. - z)*(1. + z);
        output.access(25, 2) = (1. - x)*(1. + x)*(-1. + y)*y*(-z);

        output.access(26, 0) = -x*y*(1. + y)*(1. - z)*(1. + z);
        output.access(26, 1) = (1. - x)*(1. + x)*(0.5 + y)*(1. - z)*(1. + z);
        output.access(26, 2) = (1. - x)*(1. + x)*y*(1. + y)*(-z);
        break;
      }
      case OPERATOR_D2 : {
        const auto x = input(0);
        const auto y = input(1);
        const auto z = input(2);

        // output.access is a rank-3 array with dimensions (basisCardinality_, dim0, D2Cardinality = 6)
        output.access(0,  0) = 0.25*(-1. + y)*y*(-1. + z)*z;
        output.access(0,  1) = (-0.125 + y*(0.25 - 0.25*z) + x*(0.25 + y*(-0.5 + 0.5*z) - 0.25*z) + 0.125*z)*z;
        output.access(0,  2) = y*(-0.125 + x*(0.25 + y*(-0.25 + 0.5*z) - 0.5*z) + y*(0.125 - 0.25*z) + 0.25*z);
        output.access(0,  3) = 0.25*(-1. + x)*x*(-1. + z)*z;
        output.access(0,  4) = x*(-0.125 + y*(0.25 - 0.5*z) + x*(0.125 + y*(-0.25 + 0.5*z) - 0.25*z) + 0.25*z);
        output.access(0,  5) = 0.25*(-1. + x)*x*(-1. + y)*y;

        output.access(1,  0) = 0.25*(-1. + y)*y*(-1. + z)*z;
        output.access(1,  1) = (0.125 + x*(0.25 + y*(-0.5 + 0.5*z) - 0.25*z) + y*(-0.25 + 0.25*z) - 0.125*z)*z;
        output.access(1,  2) = y*(0.125 + x*(0.25 + y*(-0.25 + 0.5*z) - 0.5*z) + y*(-0.125 + 0.25*z) - 0.25*z);
        output.access(1,  3) = 0.25*x*(1 + x)*(-1. + z)*z;
        output.access(1,  4) = x*(1. + x)*(0.125 + y*(-0.25 + 0.5*z) - 0.25*z);
        output.access(1,  5) = 0.25*x*(1 + x)*(-1. + y)*y;

        output.access(2,  0) = 0.25*y*(1 + y)*(-1. + z)*z;
        output.access(2,  1) = (0.125 + x*(0.25 + 0.5*y) + 0.25*y)*(-1. + z)*z;
        output.access(2,  2) = y*(1. + y)*(-0.125 + x*(-0.25 + 0.5*z) + 0.25*z);
        output.access(2,  3) = 0.25*x*(1 + x)*(-1. + z)*z;
        output.access(2,  4) = x*(1. + x)*(-0.125 + y*(-0.25 + 0.5*z) + 0.25*z);
        output.access(2,  5) = 0.25*x*(1 + x)*y*(1 + y);

        output.access(3,  0) = 0.25*y*(1 + y)*(-1. + z)*z;
        output.access(3,  1) = (0.125 + y*(0.25 - 0.25*z) + x*(-0.25 + y*(-0.5 + 0.5*z) + 0.25*z) - 0.125*z)*z;
        output.access(3,  2) = y*(1. + y)*(0.125 + x*(-0.25 + 0.5*z) - 0.25*z);
        output.access(3,  3) = 0.25*(-1. + x)*x*(-1. + z)*z;
        output.access(3,  4) = x*(0.125 + y*(0.25 - 0.5*z) + x*(-0.125 + y*(-0.25 + 0.5*z) + 0.25*z) - 0.25*z);
        output.access(3,  5) = 0.25*(-1. + x)*x*y*(1 + y);

        output.access(4,  0) = 0.25*(-1. + y)*y*z*(1 + z);
        output.access(4,  1) = (0.125 + x*(-0.25 + 0.5*y) - 0.25*y)*z*(1. + z);
        output.access(4,  2) = y*(0.125 + x*(-0.25 + y*(0.25 + 0.5*z) - 0.5*z) + y*(-0.125 - 0.25*z) + 0.25*z);
        output.access(4,  3) = 0.25*(-1. + x)*x*z*(1 + z);
        output.access(4,  4) = x*(0.125 + y*(-0.25 - 0.5*z) + x*(-0.125 + y*(0.25 + 0.5*z) - 0.25*z) + 0.25*z);
        output.access(4,  5) = 0.25*(-1. + x)*x*(-1. + y)*y;

        output.access(5,  0) = 0.25*(-1. + y)*y*z*(1 + z);
        output.access(5,  1) = (-0.125 + x*(-0.25 + 0.5*y) + 0.25*y)*z*(1. + z);
        output.access(5,  2) = (-1. + y)*y*(0.125 + x*(0.25 + 0.5*z) + 0.25*z);
        output.access(5,  3) = 0.25*x*(1 + x)*z*(1 + z);
        output.access(5,  4) = x*(1. + x)*(-0.125 + y*(0.25 + 0.5*z) - 0.25*z);
        output.access(5,  5) = 0.25*x*(1 + x)*(-1. + y)*y;

        output.access(6,  0) = 0.25*y*(1 + y)*z*(1 + z);
        output.access(6,  1) = (0.125 + x*(0.25 + 0.5*y) + 0.25*y)*z*(1. + z);
        output.access(6,  2) = y*(1. + y)*(0.125 + x*(0.25 + 0.5*z) + 0.25*z);
        output.access(6,  3) = 0.25*x*(1 + x)*z*(1 + z);
        output.access(6,  4) = x*(1. + x)*(0.125 + y*(0.25 + 0.5*z) + 0.25*z);
        output.access(6,  5) = 0.25*x*(1 + x)*y*(1 + y);

        output.access(7,  0) = 0.25*y*(1 + y)*z*(1 + z);
        output.access(7,  1) = (-0.125 + x*(0.25 + 0.5*y) - 0.25*y)*z*(1. + z);
        output.access(7,  2) = y*(1. + y)*(-0.125 + x*(0.25 + 0.5*z) - 0.25*z);
        output.access(7,  3) = 0.25*(-1. + x)*x*z*(1 + z);
        output.access(7,  4) = (-1. + x)*x*(0.125 + y*(0.25 + 0.5*z) + 0.25*z);
        output.access(7,  5) = 0.25*(-1. + x)*x*y*(1 + y);

        output.access(8,  0) = -0.5*(-1. + y)*y*(-1. + z)*z;
        output.access(8,  1) = (0.  + x*(-0.5 + y))*z + (x*(0.5 - y) )*(z*z);
        output.access(8,  2) = (y*y)*(x*(0.5 - z) ) + y*(x*(-0.5 + z));
        output.access(8,  3) = 0.5*(1. - x)*(1. + x)*(-1. + z)*z;
        output.access(8,  4) = 0.25 + (x*x)*(-0.25 + y*(0.5 - z) + 0.5*z) - 0.5*z + y*(-0.5 + z);
        output.access(8,  5) = 0.5*(1. - x)*(1. + x)*(-1. + y)*y;

        output.access(9,  0) = 0.5*(1. - y)*(1. + y)*(-1. + z)*z;
        output.access(9,  1) = (0.5*y + x*(y))*z + (x*(-y) - 0.5*y)*(z*z);
        output.access(9,  2) = -0.25 + (y*y)*(0.25 - 0.5*z) + 0.5*z + x*(-0.5 + (y*y)*(0.5 - z) + z);
        output.access(9,  3) = -0.5*x*(1 + x)*(-1. + z)*z;
        output.access(9,  4) = x*(y*(0.5 - z) ) + (x*x)*(y*(0.5 - z) );
        output.access(9,  5) = 0.5*x*(1 + x)*(1. - y)*(1. + y);

        output.access(10, 0) = -0.5*y*(1 + y)*(-1. + z)*z;
        output.access(10, 1) = (0.  + x*(0.5 + y))*z + (x*(-0.5 - y) )*(z*z);
        output.access(10, 2) = y*(x*(0.5 - z) ) + (y*y)*(x*(0.5 - z) );
        output.access(10, 3) = 0.5*(1. - x)*(1. + x)*(-1. + z)*z;
        output.access(10, 4) = -0.25 + (x*x)*(0.25 + y*(0.5 - z) - 0.5*z) + 0.5*z + y*(-0.5 + z);
        output.access(10, 5) = 0.5*(1. - x)*(1. + x)*y*(1 + y);

        output.access(11, 0) = 0.5*(1. - y)*(1. + y)*(-1. + z)*z;
        output.access(11, 1) = (-0.5*y + x*(y))*z + (x*(-y) + 0.5*y)*(z*z);
        output.access(11, 2) = 0.25 + (y*y)*(-0.25 + 0.5*z) - 0.5*z + x*(-0.5 + (y*y)*(0.5 - z) + z);
        output.access(11, 3) = -0.5*(-1. + x)*x*(-1. + z)*z;
        output.access(11, 4) = (x*x)*(y*(0.5 - z) ) + x*(y*(-0.5 + z));
        output.access(11, 5) = 0.5*(-1. + x)*x*(1. - y)*(1. + y);

        output.access(12, 0) = 0.5*(-1. + y)*y*(1. - z)*(1. + z);
        output.access(12, 1) = 0.25  - 0.25*(z*z) + y*(-0.5  + 0.5*(z*z)) + x*(-0.5  + 0.5*(z*z) + y*(1.  - (z*z)));
        output.access(12, 2) = (y*y)*(x*(-z) + 0.5*z) + y*(-0.5*z + x*(z));
        output.access(12, 3) = 0.5*(-1. + x)*x*(1. - z)*(1. + z);
        output.access(12, 4) = (x*x)*(y*(-z) + 0.5*z) + x*(-0.5*z + y*(z));
        output.access(12, 5) = -0.5*(-1. + x)*x*(-1. + y)*y;

        output.access(13, 0) = 0.5*(-1. + y)*y*(1. - z)*(1. + z);
        output.access(13, 1) = -0.25  + 0.25*(z*z) + y*(0.5  - 0.5*(z*z)) + x*(-0.5  + 0.5*(z*z) + y*(1.  - (z*z)));
        output.access(13, 2) = (y*y)*(x*(-z) - 0.5*z) + y*(0.5*z + x*(z));
        output.access(13, 3) = 0.5*x*(1 + x)*(1. - z)*(1. + z);
        output.access(13, 4) = x*(y*(-z) + 0.5*z) + (x*x)*(y*(-z) + 0.5*z);
        output.access(13, 5) = -0.5*x*(1 + x)*(-1. + y)*y;

        output.access(14, 0) = 0.5*y*(1 + y)*(1. - z)*(1. + z);
        output.access(14, 1) = 0.25  - 0.25*(z*z) + y*(0.5  - 0.5*(z*z)) + x*(0.5  - 0.5*(z*z) + y*(1.  - (z*z)));
        output.access(14, 2) = y*(x*(-z) - 0.5*z) + (y*y)*(x*(-z) - 0.5*z);
        output.access(14, 3) = 0.5*x*(1 + x)*(1. - z)*(1. + z);
        output.access(14, 4) = x*(y*(-z) - 0.5*z) + (x*x)*(y*(-z) - 0.5*z);
        output.access(14, 5) = -0.5*x*(1 + x)*y*(1 + y);

        output.access(15, 0) = 0.5*y*(1 + y)*(1. - z)*(1. + z);
        output.access(15, 1) = -0.25  + 0.25*(z*z) + y*(-0.5  + 0.5*(z*z)) + x*(0.5  - 0.5*(z*z) + y*(1.  - (z*z)));
        output.access(15, 2) = y*(x*(-z) + 0.5*z) + (y*y)*(x*(-z) + 0.5*z);
        output.access(15, 3) = 0.5*(-1. + x)*x*(1. - z)*(1. + z);
        output.access(15, 4) = (x*x)*(y*(-z) - 0.5*z) + x*(0.5*z + y*(z));
        output.access(15, 5) = -0.5*(-1. + x)*x*y*(1 + y);

        output.access(16, 0) = -0.5*(-1. + y)*y*z*(1 + z);
        output.access(16, 1) = (x*(0.5 - y) )*z + (x*(0.5 - y) )*(z*z);
        output.access(16, 2) = (y*y)*(x*(-0.5 - z) ) + y*(x*(0.5 + z));
        output.access(16, 3) = 0.5*(1. - x)*(1. + x)*z*(1 + z);
        output.access(16, 4) = -0.25 + (x*x)*(0.25 + y*(-0.5 - z) + 0.5*z) - 0.5*z + y*(0.5 + z);
        output.access(16, 5) = 0.5*(1. - x)*(1. + x)*(-1. + y)*y;

        output.access(17, 0) = 0.5*(1. - y)*(1. + y)*z*(1 + z);
        output.access(17, 1) = (x*(-y) - 0.5*y)*z + (x*(-y) - 0.5*y)*(z*z);
        output.access(17, 2) = 0.25 + (y*y)*(-0.25 - 0.5*z) + 0.5*z + x*(0.5 + (y*y)*(-0.5 - z) + z);
        output.access(17, 3) = -0.5*x*(1 + x)*z*(1 + z);
        output.access(17, 4) = x*(y*(-0.5 - z) ) + (x*x)*(y*(-0.5 - z) );
        output.access(17, 5) = 0.5*x*(1 + x)*(1. - y)*(1. + y);

        output.access(18, 0) = -0.5*y*(1 + y)*z*(1 + z);
        output.access(18, 1) = (x*(-0.5 - y) )*z + (x*(-0.5 - y) )*(z*z);
        output.access(18, 2) = y*(x*(-0.5 - z) ) + (y*y)*(x*(-0.5 - z) );
        output.access(18, 3) = 0.5*(1. - x)*(1. + x)*z*(1 + z);
        output.access(18, 4) = 0.25 + (x*x)*(-0.25 + y*(-0.5 - z) - 0.5*z) + 0.5*z + y*(0.5 + z);
        output.access(18, 5) = 0.5*(1. - x)*(1. + x)*y*(1 + y);

        output.access(19, 0) = 0.5*(1. - y)*(1. + y)*z*(1 + z);
        output.access(19, 1) = (x*(-y) + 0.5*y)*z + (x*(-y) + 0.5*y)*(z*z);
        output.access(19, 2) = -0.25 + (y*y)*(0.25 + 0.5*z) - 0.5*z + x*(0.5 + (y*y)*(-0.5 - z) + z);
        output.access(19, 3) = -0.5*(-1. + x)*x*z*(1 + z);
        output.access(19, 4) = (x*x)*(y*(-0.5 - z) ) + x*(y*(0.5 + z));
        output.access(19, 5) = 0.5*(-1. + x)*x*(1. - y)*(1. + y);

        output.access(20, 0) = -2.*(1. - y)*(1. + y)*(1. - z)*(1. + z);
        output.access(20, 1) = -4.*x*y*(-1. + z*z);
        output.access(20, 2) = x*((y*y)*(-4.*z) + 4.*z);
        output.access(20, 3) = -2.*(1. - x)*(1. + x)*(1. - z)*(1. + z);
        output.access(20, 4) = (x*x)*(y*(-4.*z) ) + y*(4.*z);
        output.access(20, 5) = -2.*(1. - x)*(1. + x)*(1. - y)*(1. + y);

        output.access(21, 0) = -(1. - y)*(1. + y)*(-1. + z)*z;
        output.access(21, 1) = (x*(-2.*y) )*z + (0.  + x*(2.*y))*(z*z);
        output.access(21, 2) =  x*(1. - 2.*z + (y*y)*(-1. + 2.*z));
        output.access(21, 3) = -(1. - x)*(1. + x)*(-1. + z)*z;
        output.access(21, 4) = y*(1. - 2.*z)  + (x*x)*(y*(-1. + 2.*z));
        output.access(21, 5) = (1. - x)*(1. + x)*(1. - y)*(1. + y);

        output.access(22, 0) = -(1. - y)*(1. + y)*z*(1 + z);
        output.access(22, 1) = (0.  + x*(2.*y))*z + (0.  + x*(2.*y))*(z*z);
        output.access(22, 2) = x*(-1. - 2.*z + (y*y)*(1. + 2.*z));
        output.access(22, 3) = -(1. - x)*(1. + x)*z*(1 + z);
        output.access(22, 4) = y*(-1. - 2.*z) + (x*x)*(y*(1. + 2.*z));
        output.access(22, 5) = (1. - x)*(1. + x)*(1. - y)*(1. + y);

        output.access(23, 0) = (1. - y)*(1. + y)*(1. - z)*(1. + z);
        output.access(23, 1) = (-1. + 2.*x)*y*(-1. + z*z);
        output.access(23, 2) = (-1. + 2.*x)*(-1. + y*y)*z;
        output.access(23, 3) =-(-1. + x)*x*(1. - z)*(1. + z);
        output.access(23, 4) =  2.*(-1. + x)*x*y*z;
        output.access(23, 5) =-(-1. + x)*x*(1. - y)*(1. + y);

        output.access(24, 0) = (1. - y)*(1. + y)*(1. - z)*(1. + z);
        output.access(24, 1) = (1. + 2.*x)*y*(-1. + z*z);
        output.access(24, 2) = (1. + 2.*x)*(-1. + y*y)*z;
        output.access(24, 3) = x*(1. + x)*(-1. + z)*(1. + z);
        output.access(24, 4) = 2.*x*(1. + x)*y*z;
        output.access(24, 5) = x*(1. + x)*(-1. + y)*(1. + y);

        output.access(25, 0) = -(-1. + y)*y*(1. - z)*(1. + z);
        output.access(25, 1) = x*(-1. + 2.*y)*(-1. + z*z);
        output.access(25, 2) = 2.*x*(-1. + y)*y*z;
        output.access(25, 3) = (1. - x)*(1. + x)*(1. - z)*(1. + z);
        output.access(25, 4) = (-1. + x*x)*(-1. + 2.*y)*z;
        output.access(25, 5) =-(1. - x)*(1. + x)*(-1. + y)*y;

        output.access(26, 0) =  y*(1. + y)*(-1. + z)*(1. + z);
        output.access(26, 1) =  x*(1. + 2.*y)*(-1. + z*z);
        output.access(26, 2) =  2.*x*y*(1. + y)*z;
        output.access(26, 3) =  (-1. + x)*(1. + x)*(-1. + z)*(1. + z);
        output.access(26, 4) =  (-1. + x*x)*(1. + 2.*y)*z;
        output.access(26, 5) =  (-1. + x)*(1. + x)*y*(1. + y);
        break;
      }
      case OPERATOR_D3 : {
        const auto x = input(0);
        const auto y = input(1);
        const auto z = input(2);

        output.access(0, 0) = 0.;
        output.access(0, 1) = ((-1.+ 2.*y)*(-1.+ z)*z)/4.;
        output.access(0, 2) = ((-1.+ y)*y*(-1.+ 2.*z))/4.;
        output.access(0, 3) = ((-1.+ 2.*x)*(-1.+ z)*z)/4.;
        output.access(0, 4) = ((-1.+ 2.*x)*(-1.+ 2.*y)*(-1.+ 2.*z))/8.;
        output.access(0, 5) = ((-1.+ 2.*x)*(-1.+ y)*y)/4.;
        output.access(0, 6) = 0.;
        output.access(0, 7) = ((-1.+ x)*x*(-1.+ 2.*z))/4.;
        output.access(0, 8) = ((-1.+ x)*x*(-1.+ 2.*y))/4.;
        output.access(0, 9) = 0.;

        output.access(1,  0) = 0.;
        output.access(1,  1) = ((-1.+ 2.*y)*(-1.+ z)*z)/4.;
        output.access(1,  2) = ((-1.+ y)*y*(-1.+ 2.*z))/4.;
        output.access(1,  3) = ((1.+ 2.*x)*(-1.+ z)*z)/4.;
        output.access(1,  4) = ((1.+ 2.*x)*(-1.+ 2.*y)*(-1.+ 2.*z))/8.;
        output.access(1,  5) = ((1.+ 2.*x)*(-1.+ y)*y)/4.;
        output.access(1,  6) = 0.;
        output.access(1,  7) = (x*(1.+ x)*(-1.+ 2.*z))/4.;
        output.access(1,  8) = (x*(1.+ x)*(-1.+ 2.*y))/4.;
        output.access(1,  9) = 0.;

        output.access(2,  0) = 0.;
        output.access(2,  1) = ((1.+ 2.*y)*(-1.+ z)*z)/4.;
        output.access(2,  2) = (y*(1.+ y)*(-1.+ 2.*z))/4.;
        output.access(2,  3) = ((1.+ 2.*x)*(-1.+ z)*z)/4.;
        output.access(2,  4) = ((1.+ 2.*x)*(1.+ 2.*y)*(-1.+ 2.*z))/8.;
        output.access(2,  5) = ((1.+ 2.*x)*y*(1.+ y))/4.;
        output.access(2,  6) = 0.;
        output.access(2,  7) = (x*(1.+ x)*(-1.+ 2.*z))/4.;
        output.access(2,  8) = (x*(1.+ x)*(1.+ 2.*y))/4.;
        output.access(2,  9) = 0.;

        output.access(3,  0) = 0.;
        output.access(3,  1) = ((1.+ 2.*y)*(-1.+ z)*z)/4.;
        output.access(3,  2) = (y*(1.+ y)*(-1.+ 2.*z))/4.;
        output.access(3,  3) = ((-1.+ 2.*x)*(-1.+ z)*z)/4.;
        output.access(3,  4) = ((-1.+ 2.*x)*(1.+ 2.*y)*(-1.+ 2.*z))/8.;
        output.access(3,  5) = ((-1.+ 2.*x)*y*(1.+ y))/4.;
        output.access(3,  6) = 0.;
        output.access(3,  7) = ((-1.+ x)*x*(-1.+ 2.*z))/4.;
        output.access(3,  8) = ((-1.+ x)*x*(1.+ 2.*y))/4.;
        output.access(3,  9) = 0.;

        output.access(4,  0) = 0.;
        output.access(4,  1) = ((-1.+ 2.*y)*z*(1.+ z))/4.;
        output.access(4,  2) = ((-1.+ y)*y*(1.+ 2.*z))/4.;
        output.access(4,  3) = ((-1.+ 2.*x)*z*(1.+ z))/4.;
        output.access(4,  4) = ((-1.+ 2.*x)*(-1.+ 2.*y)*(1.+ 2.*z))/8.;
        output.access(4,  5) = ((-1.+ 2.*x)*(-1.+ y)*y)/4.;
        output.access(4,  6) = 0.;
        output.access(4,  7) = ((-1.+ x)*x*(1.+ 2.*z))/4.;
        output.access(4,  8) = ((-1.+ x)*x*(-1.+ 2.*y))/4.;
        output.access(4,  9) = 0.;

        output.access(5,  0) = 0.;
        output.access(5,  1) = ((-1.+ 2.*y)*z*(1.+ z))/4.;
        output.access(5,  2) = ((-1.+ y)*y*(1.+ 2.*z))/4.;
        output.access(5,  3) = ((1.+ 2.*x)*z*(1.+ z))/4.;
        output.access(5,  4) = ((1.+ 2.*x)*(-1.+ 2.*y)*(1.+ 2.*z))/8.;
        output.access(5,  5) = ((1.+ 2.*x)*(-1.+ y)*y)/4.;
        output.access(5,  6) = 0.;
        output.access(5,  7) = (x*(1.+ x)*(1.+ 2.*z))/4.;
        output.access(5,  8) = (x*(1.+ x)*(-1.+ 2.*y))/4.;
        output.access(5,  9) = 0.;

        output.access(6,  0) = 0.;
        output.access(6,  1) = ((1.+ 2.*y)*z*(1.+ z))/4.;
        output.access(6,  2) = (y*(1.+ y)*(1.+ 2.*z))/4.;
        output.access(6,  3) = ((1.+ 2.*x)*z*(1.+ z))/4.;
        output.access(6,  4) = ((1.+ 2.*x)*(1.+ 2.*y)*(1.+ 2.*z))/8.;
        output.access(6,  5) = ((1.+ 2.*x)*y*(1.+ y))/4.;
        output.access(6,  6) = 0.;
        output.access(6,  7) = (x*(1.+ x)*(1.+ 2.*z))/4.;
        output.access(6,  8) = (x*(1.+ x)*(1.+ 2.*y))/4.;
        output.access(6,  9) = 0.;

        output.access(7,  0) = 0.;
        output.access(7,  1) = ((1.+ 2.*y)*z*(1.+ z))/4.;
        output.access(7,  2) = (y*(1.+ y)*(1.+ 2.*z))/4.;
        output.access(7,  3) = ((-1.+ 2.*x)*z*(1.+ z))/4.;
        output.access(7,  4) = ((-1.+ 2.*x)*(1.+ 2.*y)*(1.+ 2.*z))/8.;
        output.access(7,  5) = ((-1.+ 2.*x)*y*(1.+ y))/4.;
        output.access(7,  6) = 0.;
        output.access(7,  7) = ((-1.+ x)*x*(1.+ 2.*z))/4.;
        output.access(7,  8) = ((-1.+ x)*x*(1.+ 2.*y))/4.;
        output.access(7,  9) = 0.;

        output.access(8,  0) = 0.;
        output.access(8,  1) = -((-1.+ 2.*y)*(-1.+ z)*z)/2.;
        output.access(8,  2) = -((-1.+ y)*y*(-1.+ 2.*z))/2.;
        output.access(8,  3) = -(x*(-1.+ z)*z);
        output.access(8,  4) = -(x*(-1.+ 2.*y)*(-1.+ 2.*z))/2.;
        output.access(8,  5) = -(x*(-1.+ y)*y);
        output.access(8,  6) = 0.;
        output.access(8,  7) = -((-1.+ (x*x))*(-1.+ 2.*z))/2.;
        output.access(8,  8) = -((-1.+ (x*x))*(-1.+ 2.*y))/2.;
        output.access(8,  9) = 0.;

        output.access(9,  0) = 0.;
        output.access(9,  1) = -(y*(-1.+ z)*z);
        output.access(9,  2) = -((-1.+ (y*y))*(-1.+ 2.*z))/2.;
        output.access(9,  3) = -((1.+ 2.*x)*(-1.+ z)*z)/2.;
        output.access(9,  4) = -((1.+ 2.*x)*y*(-1.+ 2.*z))/2.;
        output.access(9,  5) = -((1.+ 2.*x)*(-1.+ (y*y)))/2.;
        output.access(9,  6) = 0.;
        output.access(9,  7) = -(x*(1.+ x)*(-1.+ 2.*z))/2.;
        output.access(9,  8) = -(x*(1.+ x)*y);
        output.access(9,  9) = 0.;

        output.access(10, 0) = 0.;
        output.access(10, 1) = -((1.+ 2.*y)*(-1.+ z)*z)/2.;
        output.access(10, 2) = -(y*(1.+ y)*(-1.+ 2.*z))/2.;
        output.access(10, 3) = -(x*(-1.+ z)*z);
        output.access(10, 4) = -(x*(1.+ 2.*y)*(-1.+ 2.*z))/2.;
        output.access(10, 5) = -(x*y*(1.+ y));
        output.access(10, 6) = 0.;
        output.access(10, 7) =  -((-1.+ (x*x))*(-1.+ 2.*z))/2.;
        output.access(10, 8) = -((-1.+ (x*x))*(1.+ 2.*y))/2.;
        output.access(10, 9) = 0.;

        output.access(11, 0) = 0.;
        output.access(11, 1) = -(y*(-1.+ z)*z);
        output.access(11, 2) = -((-1.+ (y*y))*(-1.+ 2.*z))/2.;
        output.access(11, 3) = -((-1.+ 2.*x)*(-1.+ z)*z)/2.;
        output.access(11, 4) = -((-1.+ 2.*x)*y*(-1.+ 2.*z))/2.;
        output.access(11, 5) = -((-1.+ 2.*x)*(-1.+ (y*y)))/2.;
        output.access(11, 6) = 0.;
        output.access(11, 7) = -((-1.+ x)*x*(-1.+ 2.*z))/2.;
        output.access(11, 8) = -((-1.+ x)*x*y);
        output.access(11, 9) = 0.;

        output.access(12, 0) = 0.;
        output.access(12, 1) = -((-1.+ 2.*y)*(-1.+ (z*z)))/2.;
        output.access(12, 2) = -((-1.+ y)*y*z);
        output.access(12, 3) = -((-1.+ 2.*x)*(-1.+ (z*z)))/2.;
        output.access(12, 4) = -((-1.+ 2.*x)*(-1.+ 2.*y)*z)/2.;
        output.access(12, 5) = -((-1.+ 2.*x)*(-1.+ y)*y)/2.;
        output.access(12, 6) = 0.;
        output.access(12, 7) = -((-1.+ x)*x*z);
        output.access(12, 8) = -((-1.+ x)*x*(-1.+ 2.*y))/2.;
        output.access(12, 9) = 0.;

        output.access(13, 0) = 0.;
        output.access(13, 1) = -((-1.+ 2.*y)*(-1.+ (z*z)))/2.;
        output.access(13, 2) = -((-1.+ y)*y*z);
        output.access(13, 3) = -((1.+ 2.*x)*(-1.+ (z*z)))/2.;
        output.access(13, 4) = -((1.+ 2.*x)*(-1.+ 2.*y)*z)/2.;
        output.access(13, 5) = -((1.+ 2.*x)*(-1.+ y)*y)/2.;
        output.access(13, 6) = 0.;
        output.access(13, 7) = -(x*(1.+ x)*z);
        output.access(13, 8) = -(x*(1.+ x)*(-1.+ 2.*y))/2.;
        output.access(13, 9) = 0.;

        output.access(14, 0) = 0.;
        output.access(14, 1) = -((1.+ 2.*y)*(-1.+ (z*z)))/2.;
        output.access(14, 2) = -(y*(1.+ y)*z);
        output.access(14, 3) = -((1.+ 2.*x)*(-1.+ (z*z)))/2.;
        output.access(14, 4) = -((1.+ 2.*x)*(1.+ 2.*y)*z)/2.;
        output.access(14, 5) = -((1.+ 2.*x)*y*(1.+ y))/2.;
        output.access(14, 6) = 0.;
        output.access(14, 7) = -(x*(1.+ x)*z);
        output.access(14, 8) = -(x*(1.+ x)*(1.+ 2.*y))/2.;
        output.access(14, 9) = 0.;

        output.access(15, 0) = 0.;
        output.access(15, 1) = -((1.+ 2.*y)*(-1.+ (z*z)))/2.;
        output.access(15, 2) = -(y*(1.+ y)*z);
        output.access(15, 3) = -((-1.+ 2.*x)*(-1.+ (z*z)))/2.;
        output.access(15, 4) = -((-1.+ 2.*x)*(1.+ 2.*y)*z)/2.;
        output.access(15, 5) = -((-1.+ 2.*x)*y*(1.+ y))/2.;
        output.access(15, 6) = 0.;
        output.access(15, 7) = -((-1.+ x)*x*z);
        output.access(15, 8) = -((-1.+ x)*x*(1.+ 2.*y))/2.;
        output.access(15, 9) = 0.;

        output.access(16, 0) = 0.;
        output.access(16, 1) = -((-1.+ 2.*y)*z*(1.+ z))/2.;
        output.access(16, 2) = -((-1.+ y)*y*(1.+ 2.*z))/2.;
        output.access(16, 3) = -(x*z*(1.+ z));
        output.access(16, 4) = -(x*(-1.+ 2.*y)*(1.+ 2.*z))/2.;
        output.access(16, 5) = -(x*(-1.+ y)*y);
        output.access(16, 6) = 0.;
        output.access(16, 7) = -((-1.+ (x*x))*(1.+ 2.*z))/2.;
        output.access(16, 8) = -((-1.+ (x*x))*(-1.+ 2.*y))/2.;
        output.access(16, 9) = 0.;

        output.access(17, 0) = 0.;
        output.access(17, 1) = -(y*z*(1.+ z));
        output.access(17, 2) = -((-1.+ (y*y))*(1.+ 2.*z))/2.;
        output.access(17, 3) = -((1.+ 2.*x)*z*(1.+ z))/2.;
        output.access(17, 4) = -((1.+ 2.*x)*y*(1.+ 2.*z))/2.;
        output.access(17, 5) = -((1.+ 2.*x)*(-1.+ (y*y)))/2.;
        output.access(17, 6) = 0.;
        output.access(17, 7) = -(x*(1.+ x)*(1.+ 2.*z))/2.;
        output.access(17, 8) = -(x*(1.+ x)*y);
        output.access(17, 9) = 0.;

        output.access(18, 0) = 0.;
        output.access(18, 1) = -((1.+ 2.*y)*z*(1.+ z))/2.;
        output.access(18, 2) = -(y*(1.+ y)*(1.+ 2.*z))/2.;
        output.access(18, 3) = -(x*z*(1.+ z));
        output.access(18, 4) = -(x*(1.+ 2.*y)*(1.+ 2.*z))/2.;
        output.access(18, 5) = -(x*y*(1.+ y));
        output.access(18, 6) = 0.;
        output.access(18, 7) = -((-1.+ (x*x))*(1.+ 2.*z))/2.;
        output.access(18, 8) = -((-1.+ (x*x))*(1.+ 2.*y))/2.;
        output.access(18, 9) = 0.;

        output.access(19, 0) = 0.;
        output.access(19, 1) = -(y*z*(1.+ z));
        output.access(19, 2) = -((-1.+ (y*y))*(1.+ 2.*z))/2.;
        output.access(19, 3) = -((-1.+ 2.*x)*z*(1.+ z))/2.;
        output.access(19, 4) = -((-1.+ 2.*x)*y*(1.+ 2.*z))/2.;
        output.access(19, 5) = -((-1.+ 2.*x)*(-1.+ (y*y)))/2.;
        output.access(19, 6) = 0.;
        output.access(19, 7) = -((-1.+ x)*x*(1.+ 2.*z))/2.;
        output.access(19, 8) = -((-1.+ x)*x*y);
        output.access(19, 9) = 0.;

        output.access(20, 0) = 0.;
        output.access(20, 1) = -4*y*(-1.+ (z*z));
        output.access(20, 2) = -4*(-1.+ (y*y))*z;
        output.access(20, 3) = -4*x*(-1.+ (z*z));
        output.access(20, 4) = -8*x*y*z;
        output.access(20, 5) = -4*x*(-1.+ (y*y));
        output.access(20, 6) = 0.;
        output.access(20, 7) = -4*(-1.+ (x*x))*z;
        output.access(20, 8) = -4*(-1.+ (x*x))*y;
        output.access(20, 9) = 0.;

        output.access(21, 0) = 0.;
        output.access(21, 1) = 2.*y*(-1.+ z)*z;
        output.access(21, 2) = (-1.+ (y*y))*(-1.+ 2.*z);
        output.access(21, 3) = 2.*x*(-1.+ z)*z;
        output.access(21, 4) = 2.*x*y*(-1.+ 2.*z);
        output.access(21, 5) = 2.*x*(-1.+ (y*y));
        output.access(21, 6) = 0.;
        output.access(21, 7) = (-1.+ (x*x))*(-1.+ 2.*z);
        output.access(21, 8) = 2.*(-1.+ (x*x))*y;
        output.access(21, 9) = 0.;

        output.access(22, 0) = 0.;
        output.access(22, 1) = 2.*y*z*(1.+ z);
        output.access(22, 2) = (-1.+ (y*y))*(1.+ 2.*z);
        output.access(22, 3) = 2.*x*z*(1.+ z);
        output.access(22, 4) = 2.*x*y*(1.+ 2.*z);
        output.access(22, 5) = 2.*x*(-1.+ (y*y));
        output.access(22, 6) = 0.;
        output.access(22, 7) = (-1.+ (x*x))*(1.+ 2.*z);
        output.access(22, 8) = 2.*(-1.+ (x*x))*y;
        output.access(22, 9) = 0.;

        output.access(23, 0) = 0.;
        output.access(23, 1) = 2.*y*(-1.+ (z*z));
        output.access(23, 2) = 2.*(-1.+ (y*y))*z;
        output.access(23, 3) = (-1.+ 2.*x)*(-1.+ (z*z));
        output.access(23, 4) = 2.*(-1.+ 2.*x)*y*z;
        output.access(23, 5) = (-1.+ 2.*x)*(-1.+ (y*y));
        output.access(23, 6) = 0.;
        output.access(23, 7) = 2.*(-1.+ x)*x*z;
        output.access(23, 8) = 2.*(-1.+ x)*x*y;
        output.access(23, 9) = 0.;

        output.access(24, 0) = 0.;
        output.access(24, 1) = 2.*y*(-1.+ (z*z));
        output.access(24, 2) = 2.*(-1.+ (y*y))*z;
        output.access(24, 3) = (1.+ 2.*x)*(-1.+ (z*z));
        output.access(24, 4) = 2.*(1.+ 2.*x)*y*z;
        output.access(24, 5) = (1.+ 2.*x)*(-1.+ (y*y));
        output.access(24, 6) = 0.;
        output.access(24, 7) = 2.*x*(1.+ x)*z;
        output.access(24, 8) = 2.*x*(1.+ x)*y;
        output.access(24, 9) = 0.;

        output.access(25, 0) = 0.;
        output.access(25, 1) = (-1.+ 2.*y)*(-1.+ (z*z));
        output.access(25, 2) = 2.*(-1.+ y)*y*z;
        output.access(25, 3) = 2.*x*(-1.+ (z*z));
        output.access(25, 4) = 2.*x*(-1.+ 2.*y)*z;
        output.access(25, 5) = 2.*x*(-1.+ y)*y;
        output.access(25, 6) = 0.;
        output.access(25, 7) = 2.*(-1.+ (x*x))*z;
        output.access(25, 8) = (-1.+ (x*x))*(-1.+ 2.*y);
        output.access(25, 9) = 0.;

        output.access(26, 0) = 0.;
        output.access(26, 1) = (1.+ 2.*y)*(-1.+ (z*z));
        output.access(26, 2) = 2.*y*(1.+ y)*z;
        output.access(26, 3) = 2.*x*(-1.+ (z*z));
        output.access(26, 4) = 2.*x*(1.+ 2.*y)*z;
        output.access(26, 5) = 2.*x*y*(1.+ y);
        output.access(26, 6) = 0.;
        output.access(26, 7) = 2.*(-1.+ (x*x))*z;
        output.access(26, 8) = (-1.+ (x*x))*(1.+ 2.*y);
        output.access(26, 9) = 0.;
        break;
      }
      case OPERATOR_D4 : {
        // Non-zero entries have Dk (derivative cardinality) indices {3,4,5,7,8,12}, all other entries are 0.
        // Intitialize array by zero and then fill only non-zero entries.
        const ordinal_type jend = output.extent(1);
        const ordinal_type iend = output.extent(0);

        for (ordinal_type j=0;j<jend;++j)
          for (ordinal_type i=0;i<iend;++i)
            output.access(i, j) = 0.0;

        const auto x = input(0);
        const auto y = input(1);
        const auto z = input(2);
        output.access(0,  3) = ((-1.+ z)*z)/2.;
        output.access(0,  4) = ((-1.+ 2.*y)*(-1.+ 2.*z))/4.;
        output.access(0,  5) = ((-1.+ y)*y)/2.;
        output.access(0,  7) = ((-1.+ 2.*x)*(-1.+ 2.*z))/4.;
        output.access(0,  8) = ((-1.+ 2.*x)*(-1.+ 2.*y))/4.;
        output.access(0,  12)= ((-1.+ x)*x)/2.;

        output.access(1,  3) = ((-1.+ z)*z)/2.;
        output.access(1,  4) = ((-1.+ 2.*y)*(-1.+ 2.*z))/4.;
        output.access(1,  5) = ((-1.+ y)*y)/2.;
        output.access(1,  7) = ((1. + 2.*x)*(-1.+ 2.*z))/4.;
        output.access(1,  8) = ((1. + 2.*x)*(-1.+ 2.*y))/4.;
        output.access(1,  12)= (x*(1. + x))/2.;

        output.access(2,  3) = ((-1.+ z)*z)/2.;
        output.access(2,  4) = ((1. + 2.*y)*(-1.+ 2.*z))/4.;
        output.access(2,  5) = (y*(1. + y))/2.;
        output.access(2,  7) = ((1. + 2.*x)*(-1.+ 2.*z))/4.;
        output.access(2,  8) =  ((1. + 2.*x)*(1. + 2.*y))/4.;
        output.access(2,  12)= (x*(1. + x))/2.;

        output.access(3,  3) = ((-1.+ z)*z)/2.;
        output.access(3,  4) = ((1. + 2.*y)*(-1.+ 2.*z))/4.;
        output.access(3,  5) = (y*(1. + y))/2.;
        output.access(3,  7) = ((-1.+ 2.*x)*(-1.+ 2.*z))/4.;
        output.access(3,  8) = ((-1.+ 2.*x)*(1. + 2.*y))/4.;
        output.access(3,  12)= ((-1.+ x)*x)/2.;

        output.access(4,  3) = (z*(1. + z))/2.;
        output.access(4,  4) = ((-1.+ 2.*y)*(1. + 2.*z))/4.;
        output.access(4,  5) = ((-1.+ y)*y)/2.;
        output.access(4,  7) = ((-1.+ 2.*x)*(1. + 2.*z))/4.;
        output.access(4,  8) = ((-1.+ 2.*x)*(-1.+ 2.*y))/4.;
        output.access(4,  12)= ((-1.+ x)*x)/2.;

        output.access(5,  3) = (z*(1. + z))/2.;
        output.access(5,  4) = ((-1.+ 2.*y)*(1. + 2.*z))/4.;
        output.access(5,  5) = ((-1.+ y)*y)/2.;
        output.access(5,  7) = ((1. + 2.*x)*(1. + 2.*z))/4.;
        output.access(5,  8) = ((1. + 2.*x)*(-1.+ 2.*y))/4.;
        output.access(5,  12)= (x*(1. + x))/2.;

        output.access(6,  3) = (z*(1. + z))/2.;
        output.access(6,  4) = ((1. + 2.*y)*(1. + 2.*z))/4.;
        output.access(6,  5) = (y*(1. + y))/2.;
        output.access(6,  7) = ((1. + 2.*x)*(1. + 2.*z))/4.;
        output.access(6,  8) = ((1. + 2.*x)*(1. + 2.*y))/4.;
        output.access(6,  12)=  (x*(1. + x))/2.;

        output.access(7,  3) = (z*(1. + z))/2.;
        output.access(7,  4) = ((1. + 2.*y)*(1. + 2.*z))/4.;
        output.access(7,  5) = (y*(1. + y))/2.;
        output.access(7,  7) = ((-1.+ 2.*x)*(1. + 2.*z))/4.;
        output.access(7,  8) = ((-1.+ 2.*x)*(1. + 2.*y))/4.;
        output.access(7,  12)= ((-1.+ x)*x)/2.;

        output.access(8,  3) = -((-1.+ z)*z);
        output.access(8,  4) = -0.5 + y + z - 2.*y*z;
        output.access(8,  5) = -((-1.+ y)*y);
        output.access(8,  7) = x - 2.*x*z;
        output.access(8,  8) = x - 2.*x*y;
        output.access(8,  12)= 1. - x*x;

        output.access(9,  3) = -((-1.+ z)*z);
        output.access(9,  4) = y - 2.*y*z;
        output.access(9,  5) = 1 - y*y;
        output.access(9,  7) = 0.5 + x - z - 2.*x*z;
        output.access(9,  8) = -((1. + 2.*x)*y);
        output.access(9,  12)= -(x*(1. + x));

        output.access(10, 3) = -((-1.+ z)*z);
        output.access(10, 4) = 0.5 + y - z - 2.*y*z;
        output.access(10, 5) = -(y*(1. + y));
        output.access(10, 7) = x - 2.*x*z;
        output.access(10, 8) = -(x*(1. + 2.*y));
        output.access(10, 12)=  1. - x*x;

        output.access(11, 3) = -((-1.+ z)*z);
        output.access(11, 4) =  y - 2.*y*z;
        output.access(11, 5) =  1. - y*y;
        output.access(11, 7) = -0.5 + x + z - 2.*x*z;
        output.access(11, 8) =  y - 2.*x*y;
        output.access(11, 12)= -((-1.+ x)*x);

        output.access(12, 3) = 1. - z*z;
        output.access(12, 4) = z - 2.*y*z;
        output.access(12, 5) = -((-1.+ y)*y);
        output.access(12, 7) =  z - 2.*x*z;
        output.access(12, 8) = -0.5 + x + y - 2.*x*y;
        output.access(12, 12)= -((-1.+ x)*x);

        output.access(13, 3) =  1. - z*z;
        output.access(13, 4) = z - 2.*y*z;
        output.access(13, 5) = -((-1.+ y)*y);
        output.access(13, 7) =  -((1. + 2.*x)*z);
        output.access(13, 8) = 0.5 + x - y - 2.*x*y;
        output.access(13, 12)= -(x*(1. + x));

        output.access(14, 3) = 1. - z*z;
        output.access(14, 4) = -((1. + 2.*y)*z);
        output.access(14, 5) = -(y*(1. + y));
        output.access(14, 7) = -((1. + 2.*x)*z);
        output.access(14, 8) = -((1. + 2.*x)*(1. + 2.*y))/2.;
        output.access(14, 12)= -(x*(1. + x));

        output.access(15, 3) =  1. - z*z;
        output.access(15, 4) = -((1. + 2.*y)*z);
        output.access(15, 5) = -(y*(1. + y));
        output.access(15, 7) = z - 2.*x*z;
        output.access(15, 8) = 0.5 + y - x*(1. + 2.*y);
        output.access(15, 12)= -((-1.+ x)*x);

        output.access(16, 3) = -(z*(1. + z));
        output.access(16, 4) = 0.5 + z - y*(1. + 2.*z);
        output.access(16, 5) = -((-1.+ y)*y);
        output.access(16, 7) = -(x*(1. + 2.*z));
        output.access(16, 8) = x - 2.*x*y;
        output.access(16, 12)= 1. - x*x;

        output.access(17, 3) = -(z*(1. + z));
        output.access(17, 4) = -(y*(1. + 2.*z));
        output.access(17, 5) = 1. - y*y;
        output.access(17, 7) = -((1. + 2.*x)*(1. + 2.*z))/2.;
        output.access(17, 8) = -((1. + 2.*x)*y);
        output.access(17, 12)= -(x*(1. + x));

        output.access(18, 3) = -(z*(1. + z));
        output.access(18, 4) = -((1. + 2.*y)*(1. + 2.*z))/2.;
        output.access(18, 5) = -(y*(1. + y));
        output.access(18, 7) =  -(x*(1. + 2.*z));
        output.access(18, 8) =  -(x*(1. + 2.*y));
        output.access(18, 12)= 1. - x*x;

        output.access(19, 3) = -(z*(1. + z));
        output.access(19, 4) = -(y*(1. + 2.*z));
        output.access(19, 5) = 1. - y*y;
        output.access(19, 7) = 0.5 + z - x*(1. + 2.*z);
        output.access(19, 8) = y - 2.*x*y;
        output.access(19, 12)= -((-1.+ x)*x);

        output.access(20, 3) = 4. - 4.*z*z;
        output.access(20, 4) = -8.*y*z;
        output.access(20, 5) = 4. - 4.*y*y;
        output.access(20, 7) = -8.*x*z;
        output.access(20, 8) = -8.*x*y;
        output.access(20, 12)= 4. - 4.*x*x;

        output.access(21, 3) = 2.*(-1.+ z)*z;
        output.access(21, 4) = 2.*y*(-1.+ 2.*z);
        output.access(21, 5) = 2.*(-1.+ y*y);
        output.access(21, 7) = 2.*x*(-1.+ 2.*z);
        output.access(21, 8) = 4.*x*y;
        output.access(21, 12)= 2.*(-1.+ x*x);

        output.access(22, 3) = 2.*z*(1. + z);
        output.access(22, 4) = 2.*(y + 2.*y*z);
        output.access(22, 5) = 2.*(-1.+ y*y);
        output.access(22, 7) = 2.*(x + 2.*x*z);
        output.access(22, 8) = 4.*x*y;
        output.access(22, 12)= 2.*(-1.+ x*x);

        output.access(23, 3) = 2.*(-1.+ z*z);
        output.access(23, 4) = 4.*y*z;
        output.access(23, 5) = 2.*(-1.+ y*y);
        output.access(23, 7) = 2.*(-1.+ 2.*x)*z;
        output.access(23, 8) = 2.*(-1.+ 2.*x)*y;
        output.access(23, 12)= 2.*(-1.+ x)*x;

        output.access(24, 3) = 2.*(-1.+ z*z);
        output.access(24, 4) = 4.*y*z;
        output.access(24, 5) = 2.*(-1.+ y*y);
        output.access(24, 7) = 2.*(z + 2.*x*z);
        output.access(24, 8) = 2.*(y + 2.*x*y);
        output.access(24, 12)= 2.*x*(1. + x);

        output.access(25, 3) =  2.*(-1.+ z*z);
        output.access(25, 4) = 2.*(-1.+ 2.*y)*z;
        output.access(25, 5) = 2.*(-1.+ y)*y;
        output.access(25, 7) = 4.*x*z;
        output.access(25, 8) = 2.*x*(-1.+ 2.*y);
        output.access(25, 12)= 2.*(-1.+ x*x);

        output.access(26, 3) = 2.*(-1.+ z*z);
        output.access(26, 4) = 2.*(z + 2.*y*z);
        output.access(26, 5) = 2.*y*(1. + y);
        output.access(26, 7) =  4.*x*z;
        output.access(26, 8) = 2.*(x + 2.*x*y);
        output.access(26, 12)= 2.*(-1.+ x*x);

        break;
      }
      case OPERATOR_MAX : {
        const ordinal_type jend = output.extent(1);
        const ordinal_type iend = output.extent(0);

        for (ordinal_type j=0;j<jend;++j)
          for (ordinal_type i=0;i<iend;++i)
            output.access(i, j) = 0.0;
        break;
      }
      default: {
        INTREPID2_TEST_FOR_ABORT( opType != OPERATOR_VALUE &&
                                  opType != OPERATOR_GRAD &&
                                  opType != OPERATOR_CURL &&
                                  opType != OPERATOR_D1 &&
                                  opType != OPERATOR_D2 &&
                                  opType != OPERATOR_D3 &&
                                  opType != OPERATOR_D4 &&
                                  opType != OPERATOR_MAX,
                                  ">>> ERROR: (Intrepid2::Basis_HGRAD_HEX_C2_FEM::Serial::getValues) operator is not supported");

      }
      }
    }

    template<typename SpT, 
             typename outputValueValueType, class ...outputValueProperties,
             typename inputPointValueType,  class ...inputPointProperties>
    void
    Basis_HGRAD_HEX_C2_FEM::
    getValues(       Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
               const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
               const EOperator operatorType ) {
      typedef          Kokkos::DynRankView<outputValueValueType,outputValueProperties...>         outputValueViewType;
      typedef          Kokkos::DynRankView<inputPointValueType, inputPointProperties...>          inputPointViewType;
      typedef typename ExecSpace<typename inputPointViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;

      // Number of evaluation points = dim 0 of inputPoints
      const auto loopSize = inputPoints.extent(0);
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
      case OPERATOR_CURL:
        INTREPID2_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_CURL), std::invalid_argument,
                                      ">>> ERROR (Basis_HGRAD_HEX_C2_FEM): CURL is invalid operator for rank-0 (scalar) functions in 3D");
        break;

      case OPERATOR_DIV:
        INTREPID2_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_DIV), std::invalid_argument,
                                      ">>> ERROR (Basis_HGRAD_HEX_C2_FEM): DIV is invalid operator for rank-0 (scalar) functions in 3D");
        break;

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
        typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_D4> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
        break;
      }
      case OPERATOR_D5:
      case OPERATOR_D6: {
        INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                                      ">>> ERROR (Basis_HGRAD_HEX_C2_FEM): operator not supported");
        // break; commented out becuase this always throws
      }
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
                                      ">>> ERROR (Basis_HGRAD_HEX_C2_FEM): Invalid operator type");
      }
      }
    }

  }

  // -------------------------------------------------------------------------------------

  template<typename SpT, typename OT, typename PT>
  Basis_HGRAD_HEX_C2_FEM<SpT,OT,PT>::
  Basis_HGRAD_HEX_C2_FEM() {
    this -> basisCardinality_  = 27;
    this -> basisDegree_       = 2;
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >() );
    this -> basisType_         = BASIS_FEM_DEFAULT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;

    // initialize tags
    {
      // Basis-dependent intializations
      const ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
      const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
      const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
      const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

      // An array with local DoF tags assigned to basis functions, in the order of their local enumeration
      ordinal_type tags[108]  = { 0, 0, 0, 1,     // Nodes 0 to 7 follow vertex order of the topology
                                  0, 1, 0, 1,
                                  0, 2, 0, 1,
                                  0, 3, 0, 1,
                                  0, 4, 0, 1,
                                  0, 5, 0, 1,
                                  0, 6, 0, 1,
                                  0, 7, 0, 1,
                                  1, 0, 0, 1,      // Node 8  -> edge 0
                                  1, 1, 0, 1,      // Node 9  -> edge 1
                                  1, 2, 0, 1,      // Node 10 -> edge 2
                                  1, 3, 0, 1,      // Node 11 -> edge 3
                                  1, 8, 0, 1,      // Node 12 -> edge 8
                                  1, 9, 0, 1,      // Node 13 -> edge 9
                                  1,10, 0, 1,      // Node 14 -> edge 10
                                  1,11, 0, 1,      // Node 15 -> edge 11
                                  1, 4, 0, 1,      // Node 16 -> edge 4
                                  1, 5, 0, 1,      // Node 17 -> edge 5
                                  1, 6, 0, 1,      // Node 18 -> edge 6
                                  1, 7, 0, 1,      // Node 19 -> edge 7
                                  3, 0, 0, 1,      // Node 20 -> Hexahedron
                                  2, 4, 0, 1,      // Node 21 -> face 4
                                  2, 5, 0, 1,      // Node 22 -> face 5
                                  2, 3, 0, 1,      // Node 23 -> face 3
                                  2, 1, 0, 1,      // Node 24 -> face 1
                                  2, 0, 0, 1,      // Node 25 -> face 0
                                  2, 2, 0, 1,      // Node 26 -> face 2
      };

      // host tags
      ordinal_type_array_1d_host tagView(&tags[0], 108);

      // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
      this->setOrdinalTagData(this -> tagToOrdinal_,
                              this -> ordinalToTag_,
                              tagView,
                              this -> basisCardinality_,
                              tagSize,
                              posScDim,
                              posScOrd,
                              posDfOrd);
    }
    // dofCoords on host and create its mirror view to device
    Kokkos::DynRankView<typename scalarViewType::value_type,typename SpT::array_layout,Kokkos::HostSpace>
      dofCoords("dofCoordsHost", this->basisCardinality_,this->basisCellTopology_.getDimension());

    dofCoords(0,0) = -1.0;   dofCoords(0,1) = -1.0; dofCoords(0,2) = -1.0;
    dofCoords(1,0) =  1.0;   dofCoords(1,1) = -1.0; dofCoords(1,2) = -1.0;
    dofCoords(2,0) =  1.0;   dofCoords(2,1) =  1.0; dofCoords(2,2) = -1.0;
    dofCoords(3,0) = -1.0;   dofCoords(3,1) =  1.0; dofCoords(3,2) = -1.0;
    dofCoords(4,0) = -1.0;   dofCoords(4,1) = -1.0; dofCoords(4,2) =  1.0;
    dofCoords(5,0) =  1.0;   dofCoords(5,1) = -1.0; dofCoords(5,2) =  1.0;
    dofCoords(6,0) =  1.0;   dofCoords(6,1) =  1.0; dofCoords(6,2) =  1.0;
    dofCoords(7,0) = -1.0;   dofCoords(7,1) =  1.0; dofCoords(7,2) =  1.0;

    dofCoords(8,0) =   0.0;   dofCoords(8,1) =  -1.0; dofCoords(8,2) =  -1.0;
    dofCoords(9,0) =   1.0;   dofCoords(9,1) =   0.0; dofCoords(9,2) =  -1.0;
    dofCoords(10,0) =  0.0;   dofCoords(10,1) =  1.0; dofCoords(10,2) = -1.0;
    dofCoords(11,0) = -1.0;   dofCoords(11,1) =  0.0; dofCoords(11,2) = -1.0;
    dofCoords(12,0) = -1.0;   dofCoords(12,1) = -1.0; dofCoords(12,2) =  0.0;
    dofCoords(13,0) =  1.0;   dofCoords(13,1) = -1.0; dofCoords(13,2) =  0.0;
    dofCoords(14,0) =  1.0;   dofCoords(14,1) =  1.0; dofCoords(14,2) =  0.0;
    dofCoords(15,0) = -1.0;   dofCoords(15,1) =  1.0; dofCoords(15,2) =  0.0;
    dofCoords(16,0) =  0.0;   dofCoords(16,1) = -1.0; dofCoords(16,2) =  1.0;
    dofCoords(17,0) =  1.0;   dofCoords(17,1) =  0.0; dofCoords(17,2) =  1.0;
    dofCoords(18,0) =  0.0;   dofCoords(18,1) =  1.0; dofCoords(18,2) =  1.0;
    dofCoords(19,0) = -1.0;   dofCoords(19,1) =  0.0; dofCoords(19,2) =  1.0;

    dofCoords(20,0) =  0.0;   dofCoords(20,1) =  0.0; dofCoords(20,2) =  0.0;

    dofCoords(21,0) =  0.0;   dofCoords(21,1) =  0.0; dofCoords(21,2) = -1.0;
    dofCoords(22,0) =  0.0;   dofCoords(22,1) =  0.0; dofCoords(22,2) =  1.0;
    dofCoords(23,0) = -1.0;   dofCoords(23,1) =  0.0; dofCoords(23,2) =  0.0;
    dofCoords(24,0) =  1.0;   dofCoords(24,1) =  0.0; dofCoords(24,2) =  0.0;
    dofCoords(25,0) =  0.0;   dofCoords(25,1) = -1.0; dofCoords(25,2) =  0.0;
    dofCoords(26,0) =  0.0;   dofCoords(26,1) =  1.0; dofCoords(26,2) =  0.0;

    this->dofCoords_ = Kokkos::create_mirror_view(typename SpT::memory_space(), dofCoords);
    Kokkos::deep_copy(this->dofCoords_, dofCoords);

  }

}// namespace Intrepid2
#endif
