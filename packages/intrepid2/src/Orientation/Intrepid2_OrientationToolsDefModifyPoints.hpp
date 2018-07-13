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


/** \file   Intrepid2_OrientationToolsDefModifyPoints.hpp
    \brief  Definition file for functions that modify points due to orientation in the Intrepid2::Impl::OrientationTools class.
    \author Created by Kyungjoo Kim
*/
#ifndef __INTREPID2_ORIENTATIONTOOLS_DEF_MODIFY_POINTS_HPP__
#define __INTREPID2_ORIENTATIONTOOLS_DEF_MODIFY_POINTS_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {

  namespace Impl {

    // ------------------------------------------------------------------------------------
    // Modified points according to orientations
    //
    //
    template<typename VT>
    KOKKOS_INLINE_FUNCTION
    void
    OrientationTools::
    getModifiedLinePoint(VT &ot,
                         const VT pt,
                         const ordinal_type ort) {
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_ABORT( !( -1.0 <= pt && pt <= 1.0 ), 
                                ">>> ERROR (Intrepid::OrientationTools::getModifiedLinePoint): "
                                "Input point is out of range [-1, 1].");
#endif
      
      switch (ort) {
      case 0: ot =  pt; break;
      case 1: ot = -pt; break;
      default:
        INTREPID2_TEST_FOR_ABORT( true, 
                                  ">>> ERROR (Intrepid2::OrientationTools::getModifiedLinePoint): "
                                  "Orientation is invalid (0--1)." );
      }
    }
    
    template<typename JacobianViewType>
    KOKKOS_INLINE_FUNCTION
    void
    OrientationTools::
    getLineJacobian(JacobianViewType jacobian, const ordinal_type ort) {
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_ABORT( ort <0 || ort >1,
                                ">>> ERROR (Intrepid2::OrientationTools::getLineJacobian): " \
                                "Orientation is invalid (0--1)." );
#endif

      ordinal_type jac[2] = {  1,  -1 };

      jacobian(0,0) = jac[ort];
    }

    template<typename VT>
    KOKKOS_INLINE_FUNCTION
    void
    OrientationTools::
    getModifiedTrianglePoint(VT &ot0,
                             VT &ot1,
                             const VT pt0,
                             const VT pt1,
                             const ordinal_type ort) {
      const VT lambda[3] = { 1.0 - pt0 - pt1,
                             pt0,
                             pt1 };
      
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_ABORT( !( 0.0 <= lambda[0] && lambda[0] <= 1.0 ), 
                                ">>> ERROR (Intrepid::OrientationTools::getModifiedTrianglePoint): " \
                                "Computed bicentric coordinate (lamba[0]) is out of range [0, 1].");
      
      INTREPID2_TEST_FOR_ABORT( !( 0.0 <= lambda[1] && lambda[1] <= 1.0 ), 
                                ">>> ERROR (Intrepid::OrientationTools::getModifiedTrianglePoint): " \
                                "Computed bicentric coordinate (lamba[1]) is out of range [0, 1].");
      
      INTREPID2_TEST_FOR_ABORT( !( 0.0 <= lambda[2] && lambda[2] <= 1.0 ), 
                                ">>> ERROR (Intrepid::OrientationTools::getModifiedTrianglePoint): "
                                "Computed bicentric coordinate (lamba[2]) is out of range [0, 1].");
#endif
      
      switch (ort) {
      case 0: ot0 = lambda[1]; ot1 = lambda[2]; break;
      case 1: ot0 = lambda[0]; ot1 = lambda[1]; break;
      case 2: ot0 = lambda[2]; ot1 = lambda[0]; break;
        
      case 3: ot0 = lambda[2]; ot1 = lambda[1]; break;
      case 4: ot0 = lambda[0]; ot1 = lambda[2]; break;
      case 5: ot0 = lambda[1]; ot1 = lambda[0]; break;
      default:
        INTREPID2_TEST_FOR_ABORT( true, 
                                  ">>> ERROR (Intrepid2::OrientationTools::getModifiedTrianglePoint): " \
                                  "Orientation is invalid (0--5)." );
      }
    }

    template<typename JacobianViewType>
    KOKKOS_INLINE_FUNCTION
    void
    OrientationTools::
    getTriangleJacobian(JacobianViewType jacobian, const ordinal_type ort) {
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_ABORT( ort <0 || ort >5,
                                ">>> ERROR (Intrepid2::OrientationTools::getTriangleJacobian): " \
                                "Orientation is invalid (0--5)." );
#endif

      ordinal_type jac[6][2][2] = { { {  1,  0 },
                                      {  0,  1 } }, // 0
                                    { { -1, -1 },
                                      {  1,  0 } }, // 1
                                    { {  0,  1 },
                                      { -1, -1 } }, // 2
                                    { {  0,  1 },
                                      {  1,  0 } }, // 3
                                    { { -1, -1 },
                                      {  0,  1 } }, // 4
                                    { {  1,  0 },
                                      { -1, -1 } } }; // 5

      jacobian(0,0) = jac[ort][0][0];
      jacobian(0,1) = jac[ort][0][1];
      jacobian(1,0) = jac[ort][1][0];
      jacobian(1,1) = jac[ort][1][1];
    }

    template<typename VT>
    KOKKOS_INLINE_FUNCTION
    void
    OrientationTools::
    getModifiedQuadrilateralPoint(VT &ot0,
                                  VT &ot1,
                                  const VT pt0,
                                  const VT pt1,
                                  const ordinal_type ort) {
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_ABORT( !( -1.0 <= pt0 && pt0 <= 1.0 ), 
                                ">>> ERROR (Intrepid::OrientationTools::getModifiedQuadrilateralPoint): " \
                                "Input point(0) is out of range [-1, 1].");
      
      INTREPID2_TEST_FOR_ABORT( !( -1.0 <= pt1 && pt1 <= 1.0 ), 
                                ">>> ERROR (Intrepid::OrientationTools::getModifiedQuadrilateralPoint): " \
                                "Input point(1) is out of range [-1, 1].");
#endif
      
      const VT lambda[2][2] = { { pt0, -pt0 },
                                { pt1, -pt1 } };
      
      switch (ort) {
      case 0: ot0 = lambda[0][0]; ot1 = lambda[1][0]; break;
      case 1: ot0 = lambda[1][1]; ot1 = lambda[0][0]; break;
      case 2: ot0 = lambda[0][1]; ot1 = lambda[1][1]; break;
      case 3: ot0 = lambda[1][0]; ot1 = lambda[0][1]; break;
        // flip
      case 4: ot0 = lambda[1][0]; ot1 = lambda[0][0]; break;
      case 5: ot0 = lambda[0][1]; ot1 = lambda[1][0]; break;
      case 6: ot0 = lambda[1][1]; ot1 = lambda[0][1]; break;
      case 7: ot0 = lambda[0][0]; ot1 = lambda[1][1]; break;
      default:
        INTREPID2_TEST_FOR_ABORT( true, 
                                  ">>> ERROR (Intrepid2::OrientationTools::getModifiedQuadrilateralPoint): " \
                                  "Orientation is invalid (0--7)." );
      }
    }

    template<typename JacobianViewType>
    KOKKOS_INLINE_FUNCTION
    void
    OrientationTools::
    getQuadrilateralJacobian(JacobianViewType jacobian, const ordinal_type ort) {
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_ABORT( ort <0 || ort >7,
                                ">>> ERROR (Intrepid2::OrientationTools::getQuadrilateralJacobian): " \
                                "Orientation is invalid (0--7)." );
#endif

      ordinal_type jac[8][2][2] = { { {  1,  0 },
                                      {  0,  1 } }, // 0
                                    { {  0, -1 },
                                      {  1,  0 } }, // 1
                                    { { -1,  0 },
                                      {  0, -1 } }, // 2
                                    { {  0,  1 },
                                      { -1,  0 } }, // 3
                                    { {  0,  1 },
                                      {  1,  0 } }, // 4
                                    { { -1,  0 },
                                      {  0,  1 } }, // 5
                                    { {  0, -1 },
                                      { -1,  0 } }, // 6
                                    { {  1,  0 },
                                      {  0, -1 } } }; // 7

      jacobian(0,0) = jac[ort][0][0];
      jacobian(0,1) = jac[ort][0][1];
      jacobian(1,0) = jac[ort][1][0];
      jacobian(1,1) = jac[ort][1][1];
    }

    template<typename outPointViewType,
             typename refPointViewType>
    inline
    void 
    OrientationTools::
    mapToModifiedReference(outPointViewType outPoints,
                           const refPointViewType refPoints,
                           const shards::CellTopology cellTopo,
                           const ordinal_type cellOrt) {
#ifdef HAVE_INTREPID2_DEBUG
      {
        const auto cellDim = cellTopo.getDimension();
        INTREPID2_TEST_FOR_EXCEPTION( !( (1 <= cellDim) && (cellDim <= 2 ) ), std::invalid_argument,
                                      ">>> ERROR (Intrepid::OrientationTools::mapToModifiedReference): " \
                                      "Method defined only for 1 and 2-dimensional subcells.");
        
        INTREPID2_TEST_FOR_EXCEPTION( !( outPoints.extent(0) == refPoints.extent(0) ), std::invalid_argument,
                                      ">>> ERROR (Intrepid::OrientationTools::mapToModifiedReference): " \
                                      "Size of input and output point arrays does not match each other.");
      }
#endif
      
      // Apply the parametrization map to every point in parameter domain
      const ordinal_type numPts = outPoints.extent(0);
      const auto key = cellTopo.getBaseCellTopologyData()->key;
      switch (key) {
      case shards::Line<>::key : {
        for (ordinal_type pt=0;pt<numPts;++pt)
          getModifiedLinePoint(outPoints(pt, 0),
                               refPoints(pt, 0),
                               cellOrt);
        break;
      }
      case shards::Triangle<>::key : {
        for (ordinal_type pt=0;pt<numPts;++pt)
          getModifiedTrianglePoint(outPoints(pt, 0), outPoints(pt, 1),
                                   refPoints(pt, 0), refPoints(pt, 1),
                                   cellOrt);
        break;
      }
      case shards::Quadrilateral<>::key : {
        for (ordinal_type pt=0;pt<numPts;++pt)
          getModifiedQuadrilateralPoint(outPoints(pt, 0), outPoints(pt, 1),
                                        refPoints(pt, 0), refPoints(pt, 1),
                                        cellOrt);
        break;
      }
      default: {
        INTREPID2_TEST_FOR_WARNING( true, 
                                    ">>> ERROR (Intrepid2::OrientationTools::mapToModifiedReference): " \
                                    "Invalid cell topology." );
        break;
      }
      }
    }


    template<typename outPointViewType>
    inline
    void
    OrientationTools::
    getJacobianOfOrientationMap(outPointViewType jacobian,
                           const shards::CellTopology cellTopo,
                           const ordinal_type cellOrt) {
#ifdef HAVE_INTREPID2_DEBUG
      {
        const auto cellDim = cellTopo.getDimension();
        INTREPID2_TEST_FOR_EXCEPTION( !( (1 <= cellDim) && (cellDim <= 2 ) ), std::invalid_argument,
                                      ">>> ERROR (Intrepid::OrientationTools::getJacobianOfOrientationMap): " \
                                      "Method defined only for 1 and 2-dimensional subcells.");

        INTREPID2_TEST_FOR_ABORT( jacobian.rank() != 2,
                                  ">>> ERROR (Intrepid2::OrientationTools::getJacobianOfOrientationMap): " \
                                  "Jacobian should have rank 2" );

        INTREPID2_TEST_FOR_EXCEPTION( ((jacobian.extent(0) != cellDim) || (jacobian.extent(1) != cellDim)), std::invalid_argument,
                                      ">>> ERROR (Intrepid::OrientationTools::getJacobianOfOrientationMap): " \
                                      "Size of jacobian is not compatible with cell topology.");
      }
#endif

      const auto key = cellTopo.getBaseCellTopologyData()->key;
      switch (key) {
      case shards::Line<>::key :
        getLineJacobian(jacobian, cellOrt);
        break;
      case shards::Triangle<>::key :
          getTriangleJacobian(jacobian, cellOrt);
        break;
      case shards::Quadrilateral<>::key :
          getQuadrilateralJacobian(jacobian, cellOrt);
        break;
      default: {
        INTREPID2_TEST_FOR_WARNING( true,
                                    ">>> ERROR (Intrepid2::OrientationTools::mapToModifiedReference): " \
                                    "Invalid cell topology." );
        break;
      }
      }
    }



  }
}

#endif
