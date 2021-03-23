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


/** \file   Intrepid2_CellToolsDefInclusion.hpp
    \brief  Definition file for point inclusion functions of the Intrepid2::CellTools class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/
#ifndef __INTREPID2_CELLTOOLS_DEF_INCLUSION_HPP__
#define __INTREPID2_CELLTOOLS_DEF_INCLUSION_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {

  //============================================================================================//
  //                                                                                            //
  //                                        Inclusion tests                                     //
  //                                                                                            //
  //============================================================================================//

  
  template<typename DeviceType>
  template<typename pointValueType, class ...pointProperties>
  bool 
  CellTools<DeviceType>::
  checkPointInclusion( const Kokkos::DynRankView<pointValueType,pointProperties...> point,
                       const shards::CellTopology cellTopo,
                       const double               threshold) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( point.rank() != 1, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::checkPointInclusion): Point must have rank 1. ");
    INTREPID2_TEST_FOR_EXCEPTION( point.extent(0) != cellTopo.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::checkPointInclusion): Point and cell dimensions do not match. ");
#endif
    bool testResult = true;
  
    const double t = threshold; //(threshold < 0 ? threshold() : threshold); 

    // Using these values in the tests effectievly inflates the reference element to a larger one
    const double minus_one =  -1.0 - t;
    const double plus_one  =   1.0 + t;
    const double minus_zero = -t;

    // A cell with extended topology has the same reference cell as a cell with base topology. 
    // => testing for inclusion in a reference Triangle<> and a reference Triangle<6> relies on 
    // on the same set of inequalities. To eliminate unnecessary cases we switch on the base topology
    const auto key = cellTopo.getBaseKey();
    switch (key) {
    
    case shards::Line<>::key :
      if( !(minus_one <= point(0) && point(0) <= plus_one))  testResult = false;
      break;
      
    case shards::Triangle<>::key : {
      const auto distance = Util<pointValueType>::max( std::max( -point(0), -point(1) ), point(0) + point(1) - 1.0 );
      if( distance > threshold ) testResult = false;
      break;
    }
      
    case shards::Quadrilateral<>::key :
      if(!( (minus_one <= point(0) && point(0) <= plus_one) &&          
            (minus_one <= point(1) && point(1) <= plus_one) ) ) testResult = false;   
      break;
      
    case shards::Tetrahedron<>::key : {
      const auto distance = Util<pointValueType>::max(  Util<pointValueType>::max(-point(0),-point(1)), 
                                        Util<pointValueType>::max(-point(2), point(0) + point(1) + point(2) - 1)  );
      if( distance > threshold ) testResult = false;
      break;
    }
      
    case shards::Hexahedron<>::key :
      if(!((minus_one <= point(0) && point(0) <= plus_one) && 
           (minus_one <= point(1) && point(1) <= plus_one) && 
           (minus_one <= point(2) && point(2) <= plus_one)))  
        testResult = false;
      break;
      
      // The base of the reference prism is the same as the reference triangle => apply triangle test
      // to X and Y coordinates and test whether Z is in [-1,1]
    case shards::Wedge<>::key : {
      const auto distance = Util<pointValueType>::max( Util<pointValueType>::max( -point(0), -point(1) ), point(0) + point(1) - 1 );
      if( distance > threshold ||                     
          point(2) < minus_one || point(2) > plus_one) 
        testResult = false;
      break;
    }

      // The base of the reference pyramid is the same as the reference quad cell => a horizontal plane
      // through a point P(x,y,z) intersects the pyramid at a quadrilateral that equals the base quad 
      // scaled by (1-z) => P(x,y,z) is inside the pyramid <=> (x,y) is in [-1+z,1-z]^2 && 0 <= Z <= 1 
    case shards::Pyramid<>::key : {
      const auto left  = minus_one + point(2);
      const auto right = plus_one  - point(2);
      if(!( (left       <= point(0) && point(0) <= right) && \
            (left       <= point(1) && point(1) <= right) && 
            (minus_zero <= point(2) && point(2) <= plus_one) ) )  \
        testResult = false;  
      break;
    }
      
    default:
      INTREPID2_TEST_FOR_EXCEPTION( !( (key == shards::Line<>::key ) ||
                                       (key == shards::Triangle<>::key)  ||
                                       (key == shards::Quadrilateral<>::key) ||
                                       (key == shards::Tetrahedron<>::key)  ||
                                       (key == shards::Hexahedron<>::key)  ||
                                       (key == shards::Wedge<>::key)  ||
                                       (key == shards::Pyramid<>::key) ),
                                    std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::checkPointInclusion): Invalid cell topology. ");
    }
    return testResult;
  }



//   template<class Scalar>
//   template<class ArrayPoint>
//   ordinal_type CellTools<Scalar>::checkPointsetInclusion(const ArrayPoint&             points,
//                                                 const shards::CellTopology &  cellTopo, 
//                                                 const double &                threshold) {
  
//     ordinal_type rank = points.rank();  
  
// #ifdef HAVE_INTREPID2_DEBUG
//     INTREPID2_TEST_FOR_EXCEPTION( !( (1 <=getrank(points) ) && (getrank(points) <= 3) ), std::invalid_argument,
//                                   ">>> ERROR (Intrepid2::CellTools::checkPointsetInclusion): rank-1, 2 or 3 required for input points array. ");

//     // The last dimension of points array at (rank - 1) is the spatial dimension. Must equal the cell dimension.
//     INTREPID2_TEST_FOR_EXCEPTION( !((index_type) points.extent(rank - 1) == (index_type)cellTopo.getDimension() ), std::invalid_argument,
//                                   ">>> ERROR (Intrepid2::CellTools::checkPointsetInclusion): Point and cell dimensions do not match. ");
// #endif
  
//     // create temp output array depending on the rank of the input array 
//     FieldContainer<ordinal_type> inRefCell;
//     index_type dim0(0), dim1(0);
//     switch(rank) {
//     case 1: 
//       inRefCell.resize(1); 
//       break;
//     case 2: 
//       dim0 = static_cast<index_type>(points.extent(0)); 
//       inRefCell.resize(dim0); 
//       break;
//     case 3: 
//       dim0 = static_cast<index_type>(points.extent(0)); 
//       dim1 = static_cast<index_type>(points.extent(1)); 
//       inRefCell.resize(dim0, dim1); 
//       break;
//     }

//     // Call the inclusion method which returns inclusion results for all points
//     checkPointwiseInclusion(inRefCell, points, cellTopo, threshold);
  
//     // Check if any points were outside, return 0 after finding one
  
//     switch(rank) {
//     case 1:  
//       if (inRefCell(0) == 0) 
//         return 0;
//       break;
//     case 2:
//       for(index_type i = 0; i < dim0; i++ )
//         if (inRefCell(i) == 0) 
//           return 0;
//       break;
    
//     case 3: 
//       for(index_type i = 0; i < dim0; i++ )
//         for(index_type j = 0; j < dim1; j++ )
//           if (inRefCell(i,j) == 0)
//             return 0;
//       break;
//     }
  
//     return 1; //all points are inside
//   }


  template<typename DeviceType>
  template<typename inCellValueType, class ...inCellProperties,
           typename pointValueType, class ...pointProperties>
  void
  CellTools<DeviceType>::
  checkPointwiseInclusion(       Kokkos::DynRankView<inCellValueType,inCellProperties...> inCell,
                           const Kokkos::DynRankView<pointValueType,pointProperties...> points,
                           const shards::CellTopology cellTopo,
                           const double threshold ) {
#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( inCell.rank() != (points.rank()-1), std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): rank difference between inCell and points is 1.");  
      const ordinal_type iend = inCell.rank();
      for (ordinal_type i=0;i<iend;++i) {
        INTREPID2_TEST_FOR_EXCEPTION( inCell.extent(i) != points.extent(i), std::invalid_argument, 
                                      ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): dimension mismatch between inCell and points.");  
      }
    }
#endif

    // do we really need to support 3 ranks ? 
    switch (points.rank()) {
    case 2: {
      const ordinal_type iend = points.extent(0);
      for (ordinal_type i=0;i<iend;++i) {
        const auto point = Kokkos::subview(points, i, Kokkos::ALL());
        inCell(i) = checkPointInclusion(point, cellTopo, threshold);
      }
      break;
    }
    case 3: {
      const ordinal_type 
        iend = points.extent(0), 
        jend = points.extent(1); 
      for (ordinal_type i=0;i<iend;++i) 
        for (ordinal_type j=0;j<jend;++j) {
          const auto point = Kokkos::subview(points, i, j, Kokkos::ALL());
          inCell(i, j) = checkPointInclusion(point, cellTopo, threshold);
        }
      break;
    }
    }  
  }

  template<typename DeviceType>
  template<typename inCellValueType, class ...inCellProperties,
           typename pointValueType, class ...pointProperties,
           typename cellWorksetValueType, class ...cellWorksetProperties>
  void
  CellTools<DeviceType>::
  checkPointwiseInclusion(       Kokkos::DynRankView<inCellValueType,inCellProperties...> inCell,
                           const Kokkos::DynRankView<pointValueType,pointProperties...> points,
                           const Kokkos::DynRankView<cellWorksetValueType,cellWorksetProperties...> cellWorkset,
                           const shards::CellTopology cellTopo,
                           const double threshold ) {
#ifdef HAVE_INTREPID2_DEBUG
    {
      const auto key = cellTopo.getBaseKey();
      INTREPID2_TEST_FOR_EXCEPTION( key != shards::Line<>::key &&
                                    key != shards::Triangle<>::key &&
                                    key != shards::Quadrilateral<>::key &&
                                    key != shards::Tetrahedron<>::key &&                                                                                                                                               
                                    key != shards::Hexahedron<>::key &&                                                                                                                                                
                                    key != shards::Wedge<>::key &&                                                                                                                                                     
                                    key != shards::Pyramid<>::key, 
                                    std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): cell topology not supported");

      INTREPID2_TEST_FOR_EXCEPTION( points.rank() != 3, std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): Points must have rank 3. ");
      INTREPID2_TEST_FOR_EXCEPTION( cellWorkset.rank() != 3, std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): cellWorkset must have rank 3. ");
      INTREPID2_TEST_FOR_EXCEPTION( points.extent(2) != cellTopo.getDimension(), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::checkPointInclusion): Points and cell dimensions do not match. ");
      INTREPID2_TEST_FOR_EXCEPTION( cellWorkset.extent(2) != cellTopo.getDimension(), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::checkPointInclusion): cellWorkset and cell dimensions do not match. ");
      INTREPID2_TEST_FOR_EXCEPTION( points.extent(0) != cellWorkset.extent(0) , std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::checkPointInclusion): cellWorkset and points dimension(0) does not match. ");
    }
#endif    
    const ordinal_type 
      numCells = cellWorkset.extent(0),
      numPoints = points.extent(1), 
      spaceDim = cellTopo.getDimension();

    using result_layout = typename DeduceLayout< decltype(points) >::result_layout;
    using device_type = typename decltype(points)::device_type;
    auto vcprop = Kokkos::common_view_alloc_prop(points);
    using common_value_type = typename decltype(vcprop)::value_type;
    Kokkos::DynRankView< common_value_type, result_layout, device_type > refPoints ( Kokkos::view_alloc("CellTools::checkPointwiseInclusion::refPoints", vcprop), numCells, numPoints, spaceDim);
    
    // expect refPoints(CPD), points(CPD), cellWorkset(CND) 
    mapToReferenceFrame(refPoints, points, cellWorkset, cellTopo);
    checkPointwiseInclusion(inCell, refPoints, cellTopo, threshold);
  }

}

#endif
