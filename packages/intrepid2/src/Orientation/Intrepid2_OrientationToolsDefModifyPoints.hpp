// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
      // hard-coded epsilon to avoid CUDA issue with calling host code.  Would be nice to have a portable way to define this, but probably not worth the effort at the moment.
      VT eps = 1e-14; // = 10.0*std::numeric_limits<VT>::epsilon();
      INTREPID2_TEST_FOR_ABORT( !( -eps <= lambda[0] && lambda[0] <= 1.0+eps ),
                                ">>> ERROR (Intrepid::OrientationTools::getModifiedTrianglePoint): " \
                                "Computed bicentric coordinate (lamba[0]) is out of range [0, 1].");
      
      INTREPID2_TEST_FOR_ABORT( !( -eps <= lambda[1] && lambda[1] <= 1.0+eps ),
                                ">>> ERROR (Intrepid::OrientationTools::getModifiedTrianglePoint): " \
                                "Computed bicentric coordinate (lamba[1]) is out of range [0, 1].");
      
      INTREPID2_TEST_FOR_ABORT( !( -eps <= lambda[2] && lambda[2] <= 1.0+eps ),
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
      
      return mapToModifiedReference(outPoints, refPoints, cellTopo.getKey(), cellOrt);
    }


    template<typename outPointViewType,
             typename refPointViewType>
    KOKKOS_INLINE_FUNCTION
    void
    OrientationTools::
    mapToModifiedReference(outPointViewType outPoints,
                           const refPointViewType refPoints,
                           const unsigned cellTopoKey,
                           const ordinal_type cellOrt) {
      // Apply the parametrization map to every point in parameter domain
      const ordinal_type numPts = outPoints.extent(0);
      switch (cellTopoKey) {
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
        INTREPID2_TEST_FOR_ABORT( true,
                                    ">>> ERROR (Intrepid2::OrientationTools::mapToModifiedReference): " \
                                    "Invalid cell topology key." );
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

      return getJacobianOfOrientationMap(jacobian, cellTopo.getKey(), cellOrt);
    }

    template<typename outPointViewType>
    KOKKOS_INLINE_FUNCTION
    void
    OrientationTools::
    getJacobianOfOrientationMap(outPointViewType jacobian,
                           const unsigned cellTopoKey,
                           const ordinal_type cellOrt) {
      switch (cellTopoKey) {
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
        INTREPID2_TEST_FOR_ABORT( true,
                                    ">>> ERROR (Intrepid2::OrientationTools::mapToModifiedReference): " \
                                    "Invalid cell topology key." );
        break;
      }
      }
    }


    template<typename TanViewType, typename ParamViewType>
    KOKKOS_INLINE_FUNCTION
    void
    OrientationTools::getRefSubcellTangents(TanViewType tangents,
                                        const ParamViewType subcellParametrization,
                                        const unsigned subcellTopoKey,
                                        const ordinal_type subcellOrd,
                                        const ordinal_type ort){
      typename ParamViewType::non_const_value_type data[4];
      typename ParamViewType::non_const_type jac(data, 2, 2);

      ordinal_type cellDim = subcellParametrization.extent(1);
      ordinal_type numTans = subcellParametrization.extent(2)-1;

      getJacobianOfOrientationMap(jac,subcellTopoKey,ort);
      for(ordinal_type d=0; d<cellDim; ++d)
        for(ordinal_type j=0; j<numTans; ++j) {
          tangents(j,d) = 0;
          for(ordinal_type k=0; k<numTans; ++k)
            tangents(j,d) += subcellParametrization(subcellOrd,d,k+1)*jac(k,j);
        }
    }


    template<typename TanNormViewType, typename ParamViewType>
    KOKKOS_INLINE_FUNCTION
    void
    OrientationTools::getRefSideTangentsAndNormal(TanNormViewType tangentsAndNormal,
                                        const ParamViewType subcellParametrization,
                                        const unsigned subcellTopoKey,
                                        const ordinal_type subcellOrd,
                                        const ordinal_type ort){
      ordinal_type cellDim = subcellParametrization.extent(1);
      auto range = Kokkos::pair<ordinal_type, ordinal_type>(0,cellDim-1);
      auto tangents = Kokkos::subview(tangentsAndNormal,range,Kokkos::ALL());
      getRefSubcellTangents(tangents,subcellParametrization,subcellTopoKey,subcellOrd,ort);

      //compute normal
      if(cellDim==2){
        tangentsAndNormal(1,0) = tangents(0,1);
        tangentsAndNormal(1,1) = -tangents(0,0);
      } else { // cellDim==3
          //cross-product
        tangentsAndNormal(2,0) = tangents(0,1)*tangents(1,2) - tangents(0,2)*tangents(1,1);
        tangentsAndNormal(2,1) = tangents(0,2)*tangents(1,0) - tangents(0,0)*tangents(1,2);
        tangentsAndNormal(2,2) = tangents(0,0)*tangents(1,1) - tangents(0,1)*tangents(1,0);
      }
    }

    template<typename coordsViewType, typename subcellCoordsViewType, typename ParamViewType>
    KOKKOS_INLINE_FUNCTION
    void
    OrientationTools::mapSubcellCoordsToRefCell(coordsViewType cellCoords,
                                        const subcellCoordsViewType subcellCoords,
                                        const ParamViewType subcellParametrization,
                                        const unsigned subcellTopoKey,
                                        const ordinal_type subcellOrd,
                                        const ordinal_type ort){

      ordinal_type cellDim = subcellParametrization.extent(1);
      ordinal_type numCoords = subcellCoords.extent(0);
      ordinal_type subcellDim = subcellCoords.extent(1);
      auto range = Kokkos::pair<ordinal_type, ordinal_type>(0,subcellDim);
      auto modSubcellCoords = Kokkos::subview(cellCoords, Kokkos::ALL(),range);
      mapToModifiedReference(modSubcellCoords,subcellCoords,subcellTopoKey,ort);
      typename coordsViewType::value_type subCoord[2];

      for(ordinal_type i=0; i<numCoords; ++i) {
        for(ordinal_type d=0; d<subcellDim; ++d)
          subCoord[d] = modSubcellCoords(i,d);

        for(ordinal_type d=0; d<cellDim; ++d) {
          cellCoords(i,d) = subcellParametrization(subcellOrd, d, 0);
          for(ordinal_type k=0; k<subcellDim; ++k)
            cellCoords(i,d) += subcellParametrization(subcellOrd, d, k+1)*subCoord[k];
        }
      }
    }
  }
}

#endif
