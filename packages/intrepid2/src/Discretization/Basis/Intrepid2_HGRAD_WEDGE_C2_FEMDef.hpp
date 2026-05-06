// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_HGRAD_WEDGE_C2_FEMDef.hpp
    \brief  Definition file for FEM basis functions of degree 2 for H(grad) functions on WEDGE cells.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_HGRAD_WEDGE_C2_FEM_DEF_HPP__
#define __INTREPID2_HGRAD_WEDGE_C2_FEM_DEF_HPP__

namespace Intrepid2 {

  // -------------------------------------------------------------------------------------
  namespace Impl {
    
    template<bool serendipity>
    template<EOperator OpType>
    template<typename OutputViewType,
             typename inputViewType>
    KOKKOS_INLINE_FUNCTION
    void
    Basis_HGRAD_WEDGE_DEG2_FEM<serendipity>::Serial<OpType>::
    getValues(       OutputViewType output,
               const inputViewType input ) {
      typedef typename inputViewType::value_type value_type;

      if constexpr (OpType == OPERATOR_VALUE) {
        const value_type x = input(0);
        const value_type y = input(1);
        const value_type z = input(2);
        const value_type w = 1.0 - x - y;

        // output is a rank-1 array with dimensions (basisCardinality_)
        if constexpr (!serendipity) {
          output(0) =  -(w*(-1. + 2.*x + 2.*y)*(-1. + z)*z)/2.;
          output(1) =  (x*(-1. + 2.*x)*(-1. + z)*z)/2.;
          output(2) =  (y*(-1. + 2.*y)*(-1. + z)*z)/2.;
          output(3) =  -(w*(-1. + 2.*x + 2.*y)*z*(1. + z))/2.;
          output(4) =  (x*(-1. + 2.*x)*z*(1. + z))/2.;
          output(5) =  (y*(-1. + 2.*y)*z*(1. + z))/2.;

          output(6) = 2.*x*w*(-1. + z)*z;
          output(7) =  2.*x*y*(-1. + z)*z;
          output(8) = 2.*y*w*(-1. + z)*z;
          output(9) = (w*(-1. + 2.*x + 2.*y)*(-1. + z)*(1. + z));
          output(10) = -(x*(-1. + 2.*x)*(-1. + z)*(1. + z));
          output(11) = -(y*(-1. + 2.*y)*(-1. + z)*(1. + z));
          output(12) = 2.*x*w*z*(1. + z);
          output(13) =  2.*x*y*z*(1. + z);
          output(14) = 2.*y*w*z*(1. + z);
        
          output(15) = -4.*x*w*(-1. + z)*(1. + z);
          output(16) = -4.*x*y*(-1. + z)*(1. + z);
          output(17) = -4.*y*w*(-1. + z)*(1. + z);
        } else {
          output(0) =  0.5 * w *(2.0*w - 2.0 - z)*(1.0 - z);
          output(1) =  0.5 * x *(2.0*x - 2.0 - z)*(1.0 - z);
          output(2) =  0.5 * y *(2.0*y - 2.0 - z)*(1.0 - z);
          output(3) =  0.5 * w *(2.0*w - 2.0 + z)*(1.0 + z);
          output(4) =  0.5 * x *(2.0*x - 2.0 + z)*(1.0 + z);
          output(5) =  0.5 * y *(2.0*y - 2.0 + z)*(1.0 + z);
          
          output(6) = 2.0 * w * x * (1.0 - z); 
          output(7) = 2.0 * x * y * (1.0 - z); 
          output(8) = 2.0 * y * w * (1.0 - z);
          output(9) = w * (1.0 - z*z);
          output(10) = x * (1.0 - z*z);
          output(11) = y * (1.0 - z*z);
          output(12) = 2.0 * w * x * (1.0 + z);
          output(13) = 2.0 * x * y * (1.0 + z);
          output(14) = 2.0 * y * w * (1.0 + z);
        }
      }
      else if constexpr (OpType == OPERATOR_GRAD) {
        const value_type x = input(0);
        const value_type y = input(1);
        const value_type z = input(2);

        if constexpr (!serendipity) {
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
        } else {
          const value_type w = 1.0 - x - y;
      
          output(0, 0) = -(2.0*w - 1.0 - 0.5*z)*(1.0 - z);
          output(0, 1) = -(2.0*w - 1.0 - 0.5*z)*(1.0 - z);
          output(0, 2) = -0.5*w*(2.0*w - 1.0 - 2.0*z);
          
          output(1, 0) = (2.0*x - 1.0 - 0.5*z)*(1.0 - z);
          output(1, 1) = 0.0;
          output(1, 2) = -0.5*x*(2.0*x - 1.0 - 2.0*z);
          
          output(2, 0) = 0.0;
          output(2, 1) = (2.0*y - 1.0 - 0.5*z)*(1.0 - z);
          output(2, 2) = -0.5*y*(2.0*y - 1.0 - 2.0*z);
          
          output(3, 0) = -(2.0*w - 1.0 + 0.5*z)*(1.0 + z);
          output(3, 1) = -(2.0*w - 1.0 + 0.5*z)*(1.0 + z);
          output(3, 2) = 0.5*w*(2.0*w - 1.0 + 2.0*z);
          
          output(4, 0) = (2.0*x - 1.0 + 0.5*z)*(1.0 + z);
          output(4, 1) = 0.0;
          output(4, 2) = 0.5*x*(2.0*x - 1.0 + 2.0*z);
          
          output(5, 0) = 0.0;
          output(5, 1) = (2.0*y - 1.0 + 0.5*z)*(1.0 + z);
          output(5, 2) = 0.5*y*(2.0*y - 1.0 + 2.0*z);
          
          output(6, 0) =  2.0*(w - x)*(1.0 - z);
          output(6, 1) = -2.0*x*(1.0 - z);
          output(6, 2) = -2.0*x*w;
          
          output(7, 0) =  2.0*y*(1.0 - z);
          output(7, 1) =  2.0*x*(1.0 - z);
          output(7, 2) = -2.0*x*y;
          
          output(8, 0) = -2.0*y*(1.0 - z);
          output(8, 1) =  2.0*(w - y)*(1.0 - z);
          output(8, 2) = -2.0*w*y;
          
          output(9, 0) = -(1.0 - z*z);
          output(9, 1) = -(1.0 - z*z);
          output(9, 2) = -2.0*z*w;
          
          output(10,0) = (1.0 - z*z);
          output(10,1) =  0.0;
          output(10,2) = -2.0*z*x;
          
          output(11,0) =  0.0;
          output(11,1) = (1.0 - z*z);
          output(11,2) = -2.0*z*y;
          
          output(12,0) =  2.0*(w - x)*(1.0 + z);
          output(12,1) = -2.0*x*(1.0 + z);
          output(12,2) =  2.0*x*w;
          
          output(13,0) =  2.0*y*(1.0 + z); 
          output(13,1) =  2.0*x*(1.0 + z); 
          output(13,2) =  2.0*x*y;
          
          output(14,0) = -2.0*y*(1.0 + z);
          output(14,1) =  2.0*(w - y)*(1.0 + z);
          output(14,2) =  2.0*w*y;
        }
      }
      else if constexpr (OpType == OPERATOR_D2) {
        const value_type x = input(0);
        const value_type y = input(1);
        const value_type z = input(2);

        if constexpr (!serendipity) {
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
          output(4, 2) =  ((-1. + 4.*x)*(1. + 2.*z))/2.;
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

        } else { //serendipity element
        
          const value_type w = 1.0 - x - y;
          output(0, 0) =  2.0*(1.0 - z);     
          output(0, 1) =  2.0*(1.0 - z);     
          output(0, 2) =  2.0*w - 0.5 - z;
          output(0, 3) =  2.0*(1.0-z);
          output(0, 4) =  2.0*w -0.5 - z;
          output(0, 5) =  w;
          
          output(1, 0) =  2.0*(1.0 - z);
          output(1, 1) =  0.;     
          output(1, 2) = -2.0*x + 0.5 + z;     
          output(1, 3) =  0.;     
          output(1, 4) =  0.;     
          output(1, 5) =  x;     
          
          output(2, 0) =  0.;     
          output(2, 1) =  0.;     
          output(2, 2) =  0.;     
          output(2, 3) =  2.0*(1.0 - z);      
          output(2, 4) = -2.0*y + 0.5 + z;
          output(2, 5) =  y;     
          
          output(3, 0) =  2.0*(1.0 + z); 
          output(3, 1) =  2.0*(1.0 + z);
          output(3, 2) = -2.0*w + 0.5 - z;
          output(3, 3) =  2.0*(1.0 + z);
          output(3, 4) = -2.0*w + 0.5 - z;
          output(3, 5) =  w;
          
          output(4, 0) =  2.0*(1.0 + z);
          output(4, 1) =  0.; 
          output(4, 2) =  2.0*x - 0.5 + z;
          output(4, 3) =  0.;
          output(4, 4) =  0.;
          output(4, 5) =  x;
          
          output(5, 0) =  0.;
          output(5, 1) =  0.;
          output(5, 2) =  0.;
          output(5, 3) =  2.0*(1.0 + z);
          output(5, 4) =  2.0*y - 0.5 + z;
          output(5, 5) =  y;         
          
          output(6, 0) = -4.0*(1.0 - z);
          output(6, 1) = -2.0*(1.0 - z);
          output(6, 2) = -2.0*(w - x);
          output(6, 3) =  0.;
          output(6, 4) =  2.0*x;     
          output(6, 5) =  0.;
          
          output(7, 0) =  0.;
          output(7, 1) =  2.0*(1.0 - z);
          output(7, 2) = -2.0*y;
          output(7, 3) =  0.;
          output(7, 4) = -2.0*x;
          output(7, 5) =  0.;
          
          output(8, 0) =  0.;
          output(8, 1) = -2.0*(1.0 - z);
          output(8, 2) =  2.0*y;
          output(8, 3) = -4.0*(1.0 - z);
          output(8, 4) = -2.0*(w - y);
          output(8, 5) =  0.;

          output(9, 0) =  0.;
          output(9, 1) =  0.;
          output(9, 2) =  2.0*z;     
          output(9, 3) =  0.;
          output(9, 4) =  2.0*z;
          output(9, 5) = -2.0*w;
          
          output(10,0) =  0.;
          output(10,1) =  0.;
          output(10,2) = -2.0*z;
          output(10,3) =  0.;
          output(10,4) =  0.;
          output(10,5) = -2.0*x;     
          
          output(11,0) =  0.;
          output(11,1) =  0.;
          output(11,2) =  0.;
          output(11,3) =  0.;
          output(11,4) = -2.0*z;
          output(11,5) = -2.0*y;     
          
          output(12,0) = -4.0*(1.0 + z);
          output(12,1) = -2.0*(1.0 + z);
          output(12,2) =  2.0*(w - x);
          output(12,3) =  0.;
          output(12,4) = -2.0*x;     
          output(12,5) =  0.;
          
          output(13,0) =  0.;
          output(13,1) =  2.0*(1.0 + z);
          output(13,2) =  2.0*y;
          output(13,3) =  0.;
          output(13,4) =  2.0*x;
          output(13,5) =  0.;     
          
          output(14,0) =  0.;
          output(14,1) = -2.0*(1.0 + z);
          output(14,2) = -2.0*y;
          output(14,3) = -4.0*(1.0 + z);
          output(14,4) =  2.0*(w - y);
          output(14,5) =  0.;
        }
      }
      else if constexpr (OpType == OPERATOR_D3) {
        if constexpr (!serendipity) {
          const value_type x = input(0);
          const value_type y = input(1);
          const value_type z = input(2);

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
        
        } else {  //serendipity 
        
          output(0, 0) =  0.0;
          output(0, 1) =  0.0;
          output(0, 2) = -2.0;
          output(0, 3) =  0.0;
          output(0, 4) = -2.0;
          output(0, 5) = -1.0;
          output(0, 6) =  0.0;
          output(0, 7) = -2.0;
          output(0, 8) = -1.0;
          output(0, 9) =  0.0;  
            
          output(1, 0) =  0.0;
          output(1, 1) =  0.0;
          output(1, 2) = -2.0;
          output(1, 3) =  0.0;
          output(1, 4) =  0.0;
          output(1, 5) = -1.0;
          output(1, 6) =  0.0;
          output(1, 7) =  0.0;
          output(1, 8) =  0.0;
          output(1, 9) =  0.0;  
          
          output(2, 0) =  0.0;
          output(2, 1) =  0.0;
          output(2, 2) =  0.0;
          output(2, 3) =  0.0;
          output(2, 4) =  0.0;
          output(2, 5) =  0.0;
          output(2, 6) =  0.0;
          output(2, 7) = -2.0;
          output(2, 8) = -1.0;
          output(2, 9) =  0.0;  
          
          output(3, 0) =  0.0;
          output(3, 1) =  0.0;
          output(3, 2) =  2.0;
          output(3, 3) =  0.0;
          output(3, 4) =  2.0;
          output(3, 5) = -1.0;
          output(3, 6) =  0.0;
          output(3, 7) =  2.0;
          output(3, 8) = -1.0;
          output(3, 9) =  0.0;  
          
          output(4, 0) =  0.0;
          output(4, 1) =  0.0;
          output(4, 2) =  2.0;
          output(4, 3) =  0.0;
          output(4, 4) =  0.0;
          output(4, 5) =  1.0;
          output(4, 6) =  0.0;
          output(4, 7) =  0.0;
          output(4, 8) =  0.0;
          output(4, 9) =  0.0;  
          
          output(5, 0) =  0.0;
          output(5, 1) =  0.0;
          output(5, 2) =  0.0;
          output(5, 3) =  0.0;
          output(5, 4) =  0.0;
          output(5, 5) =  0.0;
          output(5, 6) =  0.0;
          output(5, 7) =  2.0;
          output(5, 8) =  1.0;
          output(5, 9) =  0.0;  
          
          output(6, 0) =  0.0;
          output(6, 1) =  0.0;
          output(6, 2) =  4.0;
          output(6, 3) =  0.0;
          output(6, 4) =  2.0;
          output(6, 5) =  0.0;
          output(6, 6) =  0.0;
          output(6, 7) =  0.0;
          output(6, 8) =  0.0;
          output(6, 9) =  0.0;  
          
          output(7, 0) =  0.0;
          output(7, 1) =  0.0;
          output(7, 2) =  0.0;
          output(7, 3) =  0.0;
          output(7, 4) = -2.0;
          output(7, 5) =  0.0;
          output(7, 6) =  0.0;
          output(7, 7) =  0.0;
          output(7, 8) =  0.0;
          output(7, 9) =  0.0;  
          
          output(8, 0) =  0.0;
          output(8, 1) =  0.0;
          output(8, 2) =  0.0;
          output(8, 3) =  0.0;
          output(8, 4) =  2.0;
          output(8, 5) =  0.0;
          output(8, 6) =  0.0;
          output(8, 7) =  4.0;
          output(8, 8) =  0.0;
          output(8, 9) =  0.0;  
          
          output(9, 0) =  0.0;
          output(9, 1) =  0.0;
          output(9, 2) =  0.0;
          output(9, 3) =  0.0;
          output(9, 4) =  0.0;
          output(9, 5) =  2.0;
          output(9, 6) =  0.0;
          output(9, 7) =  0.0;
          output(9, 8) =  2.0;
          output(9, 9) =  0.0;  
          
          output(10,0) =  0.0;
          output(10,1) =  0.0;
          output(10,2) =  0.0;
          output(10,3) =  0.0;
          output(10,4) =  0.0;
          output(10,5) = -2.0;
          output(10,6) =  0.0;
          output(10,7) =  0.0;
          output(10,8) =  0.0;
          output(10,9) =  0.0;  
          
          output(11,0) =  0.0;
          output(11,1) =  0.0;
          output(11,2) =  0.0;
          output(11,3) =  0.0;
          output(11,4) =  0.0;
          output(11,5) =  0.0;
          output(11,6) =  0.0;
          output(11,7) =  0.0;
          output(11,8) = -2.0;
          output(11,9) =  0.0;  
          
          output(12,0) =  0.0;
          output(12,1) =  0.0;
          output(12,2) = -4.0;
          output(12,3) =  0.0;
          output(12,4) = -2.0;
          output(12,5) =  0.0;
          output(12,6) =  0.0;
          output(12,7) =  0.0;
          output(12,8) =  0.0;
          output(12,9) =  0.0;  
          
          output(13,0) =  0.0;
          output(13,1) =  0.0;
          output(13,2) =  0.0;
          output(13,3) =  0.0;
          output(13,4) =  2.0;
          output(13,5) =  0.0;
          output(13,6) =  0.0;
          output(13,7) =  0.0;
          output(13,8) =  0.0;
          output(13,9) =  0.0;  
          
          output(14,0) =  0.0;
          output(14,1) =  0.0;
          output(14,2) =  0.0;
          output(14,3) =  0.0;
          output(14,4) = -2.0;
          output(14,5) =  0.0;
          output(14,6) =  0.0;
          output(14,7) = -4.0;
          output(14,8) =  0.0;
          output(14,9) =  0.0; 
        }
      }
      else if constexpr (OpType == OPERATOR_D4) {
        const ordinal_type jend = output.extent(1);
        const ordinal_type iend = output.extent(0);

        for (ordinal_type j=0;j<jend;++j)
          for (ordinal_type i=0;i<iend;++i)
            output(i, j) = 0.0;

        if constexpr (!serendipity) {

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
        }
      }
      else if constexpr (OpType == OPERATOR_MAX) {
        const ordinal_type jend = output.extent(1);
        const ordinal_type iend = output.extent(0);

        for (ordinal_type j=0;j<jend;++j)
          for (ordinal_type i=0;i<iend;++i)
            output(i, j) = 0.0;
      }

      else {
        INTREPID2_TEST_FOR_ABORT( OpType != OPERATOR_VALUE &&
                                  OpType != OPERATOR_GRAD &&
                                  OpType != OPERATOR_D2 &&
                                  OpType != OPERATOR_D3 &&
                                  OpType != OPERATOR_D4 &&
                                  OpType != OPERATOR_MAX,
                                  ">>> ERROR: (Intrepid2::Basis_HGRAD_WEDGE_DEG2_FEM::Serial::getValues) operator is not supported");
      }
    }

    template<bool serendipity>
    template<typename DT, 
             typename outputValueValueType, class ...outputValueProperties,
             typename inputPointValueType,  class ...inputPointProperties>
    void
    Basis_HGRAD_WEDGE_DEG2_FEM<serendipity>::
    getValues(       Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
               const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
               const EOperator operatorType ) {
      typedef          Kokkos::DynRankView<outputValueValueType,outputValueProperties...>         outputValueViewType;
      typedef          Kokkos::DynRankView<inputPointValueType, inputPointProperties...>          inputPointViewType;
      typedef typename ExecSpace<typename inputPointViewType::execution_space,typename DT::execution_space>::ExecSpaceType ExecSpaceType;

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
      case OPERATOR_CURL: {
        INTREPID2_TEST_FOR_EXCEPTION( operatorType == OPERATOR_CURL, std::invalid_argument,
                                      ">>> ERROR (Basis_HGRAD_WEDGE_DEG2_FEM): CURL is invalid operator for rank-0 (scalar) functions in 3D");
        break;
      }
      case OPERATOR_DIV: {
        INTREPID2_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_DIV), std::invalid_argument,
                                      ">>> ERROR (Basis_HGRAD_WEDGE_DEG2_FEM): DIV is invalid operator for rank-0 (scalar) functions in 3D");
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
        typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_D4> FunctorType;
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
                                      ">>> ERROR (Basis_HGRAD_WEDGE_DEG2_FEM): Invalid operator type");
      }
      }
    }

  }

  // -------------------------------------------------------------------------------------

  template<bool serendipity, typename DT, typename OT, typename PT>
  Basis_HGRAD_WEDGE_DEG2_FEM<serendipity,DT,OT,PT>::
  Basis_HGRAD_WEDGE_DEG2_FEM() {
    const ordinal_type spaceDim = 3;
    this->basisCardinality_     = serendipity ? 15 : 18;
    this->basisDegree_          = 2;
    this->basisCellTopologyKey_ = shards::Wedge<6>::key;
    this->basisType_            = BASIS_FEM_DEFAULT;
    this->basisCoordinates_     = COORDINATES_CARTESIAN;
    this->functionSpace_        = FUNCTION_SPACE_HGRAD;

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
                                 // following entries not used for serendipity elements
                                 2, 0, 0, 1,
                                 2, 1, 0, 1,
                                 2, 2, 0, 1
      };

      // host tags
      OrdinalTypeArray1DHost tagView(&tags[0], serendipity ? 60 : 72);

      // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
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
    Kokkos::DynRankView<typename ScalarViewType::value_type,typename DT::execution_space::array_layout,Kokkos::HostSpace>
      dofCoords("dofCoordsHost", this->basisCardinality_,spaceDim);

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

    if constexpr (!serendipity) {
      dofCoords(15,0)=  0.5;  dofCoords(15,1)=  0.0;  dofCoords(15,2)=  0.0;
      dofCoords(16,0)=  0.5;  dofCoords(16,1)=  0.5;  dofCoords(16,2)=  0.0;
      dofCoords(17,0)=  0.0;  dofCoords(17,1)=  0.5;  dofCoords(17,2)=  0.0;
    }

    this->dofCoords_ = Kokkos::create_mirror_view(typename DT::memory_space(), dofCoords);
    Kokkos::deep_copy(this->dofCoords_, dofCoords);
  }

  template<bool serendipity, typename DT, typename OT, typename PT>
  void 
  Basis_HGRAD_WEDGE_DEG2_FEM<serendipity,DT,OT,PT>::getScratchSpaceSize(       
                                    ordinal_type& perThreadSpaceSize,
                              const PointViewType inputPoints,
                              const EOperator operatorType) const {
    perThreadSpaceSize = 0;
  }

  template<bool serendipity, typename DT, typename OT, typename PT>
  KOKKOS_INLINE_FUNCTION
  void 
  Basis_HGRAD_WEDGE_DEG2_FEM<serendipity,DT,OT,PT>::getValues(       
          OutputViewType outputValues,
      const PointViewType  inputPoints,
      const EOperator operatorType,
      const typename Kokkos::TeamPolicy<typename DT::execution_space>::member_type& team_member,
      const int threadScratchLevel, 
      const ordinal_type subcellDim,
      const ordinal_type subcellOrdinal) const {

      INTREPID2_TEST_FOR_ABORT( !((subcellDim <= 0) && (subcellOrdinal == -1)),
        ">>> ERROR: (Intrepid2::Basis_HGRAD_WEDGE_DEG2_FEM::getValues), The capability of selecting subsets of basis functions has not been implemented yet.");

      (void) threadScratchLevel; //avoid unused variable warning

      const int numPoints = inputPoints.extent(0);

      switch(operatorType) {
        case OPERATOR_VALUE:
          Kokkos::parallel_for (Kokkos::TeamThreadRange (team_member, numPoints), [=] (ordinal_type& pt) {
            auto       output = Kokkos::subview( outputValues, Kokkos::ALL(), pt, Kokkos::ALL() );
            const auto input  = Kokkos::subview( inputPoints,                 pt, Kokkos::ALL() );
            using SerialValue = typename Impl::Basis_HGRAD_WEDGE_DEG2_FEM<serendipity>::template Serial<OPERATOR_VALUE>;
            SerialValue::getValues( output, input);
          });
          break;
        case OPERATOR_GRAD:
          Kokkos::parallel_for (Kokkos::TeamThreadRange (team_member, numPoints), [=] (ordinal_type& pt) {
            auto       output = Kokkos::subview( outputValues, Kokkos::ALL(), pt, Kokkos::ALL() );
            const auto input  = Kokkos::subview( inputPoints,                 pt, Kokkos::ALL() );
            using SerialGrad = typename Impl::Basis_HGRAD_WEDGE_DEG2_FEM<serendipity>::template Serial<OPERATOR_GRAD>;
            SerialGrad::getValues( output, input);
          });
          break;
        default: {}
    }
  }

}// namespace Intrepid2
#endif
