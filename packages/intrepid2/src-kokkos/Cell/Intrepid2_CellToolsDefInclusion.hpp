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


/** \file   Intrepid_CellToolsDef.hpp
    \brief  Definition file for the Intrepid2::CellTools class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/
#ifndef __INTREPID2_CELLTOOLS_DEF_HPP__
#define __INTREPID2_CELLTOOLS_DEF_HPP__

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


  template<class Scalar>
  int CellTools<Scalar>::checkPointInclusion(const Scalar*                 point,
                                             const int                     pointDim,
                                             const shards::CellTopology &  cellTopo,
                                             const double &                threshold) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( !(pointDim == (int)cellTopo.getDimension() ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::checkPointInclusion): Point and cell dimensions do not match. ");
#endif
    int testResult = 1;
  
    // Using these values in the tests effectievly inflates the reference element to a larger one
    double minus_one = -1.0 - threshold;
    double plus_one  =  1.0 + threshold;
    double minus_zero = - threshold;
  
    // A cell with extended topology has the same reference cell as a cell with base topology. 
    // => testing for inclusion in a reference Triangle<> and a reference Triangle<6> relies on 
    // on the same set of inequalities. To eliminate unnecessary cases we switch on the base topology
    unsigned key = cellTopo.getBaseCellTopologyData() -> key ;
    switch( key ) {
    
    case shards::Line<>::key :
      if( !(minus_one <= point[0] && point[0] <= plus_one))  testResult = 0;
      break;
      
    case shards::Triangle<>::key : {
      Scalar distance = std::max( std::max( -point[0], -point[1] ), point[0] + point[1] - 1.0 );
      if( distance > threshold ) testResult = 0;
      break;
    }
      
    case shards::Quadrilateral<>::key :
      if(!( (minus_one <= point[0] && point[0] <= plus_one) && \
            (minus_one <= point[1] && point[1] <= plus_one) ) ) testResult = 0;   
      break;
      
    case shards::Tetrahedron<>::key : {
      Scalar distance = std::max(  std::max(-point[0],-point[1]), \
                                   std::max(-point[2], point[0] + point[1] + point[2] - 1)  );
      if( distance > threshold ) testResult = 0;
      break;
    }
      
    case shards::Hexahedron<>::key :
      if(!((minus_one <= point[0] && point[0] <= plus_one) && \
           (minus_one <= point[1] && point[1] <= plus_one) && \
           (minus_one <= point[2] && point[2] <= plus_one)))  \
        testResult = 0;
      break;
      
      // The base of the reference prism is the same as the reference triangle => apply triangle test
      // to X and Y coordinates and test whether Z is in [-1,1]
    case shards::Wedge<>::key : {
      Scalar distance = std::max( std::max( -point[0], -point[1] ), point[0] + point[1] - 1 );
      if( distance > threshold  || \
          point[2] < minus_one || point[2] > plus_one) \
        testResult = 0;
      break;
    }

      // The base of the reference pyramid is the same as the reference quad cell => a horizontal plane
      // through a point P(x,y,z) intersects the pyramid at a quadrilateral that equals the base quad 
      // scaled by (1-z) => P(x,y,z) is inside the pyramid <=> (x,y) is in [-1+z,1-z]^2 && 0 <= Z <= 1 
    case shards::Pyramid<>::key : {
      Scalar left  = minus_one + point[2];
      Scalar right = plus_one  - point[2];
      if(!( (left       <= point[0] && point[0] <= right) && \
            (left       <= point[1] && point[1] <= right) && 
            (minus_zero <= point[2] && point[2] <= plus_one) ) )  \
        testResult = 0;  
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



  template<class Scalar>
  template<class ArrayPoint>
  int CellTools<Scalar>::checkPointsetInclusion(const ArrayPoint&             points,
                                                const shards::CellTopology &  cellTopo, 
                                                const double &                threshold) {
  
    int rank = points.rank();  
  
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( !( (1 <=getrank(points) ) && (getrank(points) <= 3) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::checkPointsetInclusion): rank-1, 2 or 3 required for input points array. ");

    // The last dimension of points array at (rank - 1) is the spatial dimension. Must equal the cell dimension.
    INTREPID2_TEST_FOR_EXCEPTION( !((index_type) points.dimension(rank - 1) == (index_type)cellTopo.getDimension() ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::checkPointsetInclusion): Point and cell dimensions do not match. ");
#endif
  
    // create temp output array depending on the rank of the input array 
    FieldContainer<int> inRefCell;
    index_type dim0(0), dim1(0);
    switch(rank) {
    case 1: 
      inRefCell.resize(1); 
      break;
    case 2: 
      dim0 = static_cast<index_type>(points.dimension(0)); 
      inRefCell.resize(dim0); 
      break;
    case 3: 
      dim0 = static_cast<index_type>(points.dimension(0)); 
      dim1 = static_cast<index_type>(points.dimension(1)); 
      inRefCell.resize(dim0, dim1); 
      break;
    }

    // Call the inclusion method which returns inclusion results for all points
    checkPointwiseInclusion(inRefCell, points, cellTopo, threshold);
  
    // Check if any points were outside, return 0 after finding one
  
    switch(rank) {
    case 1:  
      if (inRefCell(0) == 0) 
        return 0;
      break;
    case 2:
      for(index_type i = 0; i < dim0; i++ )
        if (inRefCell(i) == 0) 
          return 0;
      break;
    
    case 3: 
      for(index_type i = 0; i < dim0; i++ )
        for(index_type j = 0; j < dim1; j++ )
          if (inRefCell(i,j) == 0)
            return 0;
      break;
    }
  
    return 1; //all points are inside
  }



  template<class Scalar>
  template<class ArrayIncl, class ArrayPoint>
  void CellTools<Scalar>::checkPointwiseInclusion(ArrayIncl &                   inRefCell,
                                                  const ArrayPoint &            points,
                                                  const shards::CellTopology &  cellTopo, 
                                                  const double &                threshold) {
    int apRank   = points.rank();
  
#ifdef HAVE_INTREPID2_DEBUG
  
    // Verify that points and inRefCell have correct ranks and dimensions
    std::string errmsg = ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion):";
    if(getrank(points) == 1) {
      INTREPID2_TEST_FOR_EXCEPTION( !(getrank(inRefCell) == 1 ), std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): rank-1 input array requires rank-1 output array.");  
      INTREPID2_TEST_FOR_EXCEPTION( !(static_cast<index_type>(inRefCell.dimension(0)) == 1), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): rank-1 input array requires dim0 = 1 for output array.");  
    }
    else if(getrank(points) == 2){
      INTREPID2_TEST_FOR_EXCEPTION( !(getrank(inRefCell) == 1 ), std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): rank-2 input array requires rank-1 output array.");  
      // dimension 0 of the arrays must match
      INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionMatch( errmsg, inRefCell, 0,  points, 0), std::invalid_argument, errmsg);
    }
    else if (getrank(points) == 3) {
      INTREPID2_TEST_FOR_EXCEPTION( !(getrank(inRefCell) == 2 ), std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): rank-3 input array requires rank-2 output array.");  
      // dimensions 0 and 1 of the arrays must match
      INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionMatch( errmsg, inRefCell, 0,1,  points, 0,1), std::invalid_argument, errmsg);
    }
    else{
      INTREPID2_TEST_FOR_EXCEPTION( !( (getrank(points) == 1) || (getrank(points) == 2) || (getrank(points) == 3) ), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): rank-1, 2 or 3 required for input points array. ");      
    }    
  
    // The last dimension of points array at (rank - 1) is the spatial dimension. Must equal the cell dimension.
    INTREPID2_TEST_FOR_EXCEPTION( !((index_type)points.dimension(apRank - 1) == (index_type)cellTopo.getDimension() ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): Point and cell dimensions do not match. ");
  
#endif
  
    // Initializations
    index_type dim0     = 1;
    index_type dim1     = 1;
    index_type pointDim = 0;
    switch(apRank) {
    case 1:
      pointDim = static_cast<index_type>(points.dimension(0));
      break;
    case 2:
      dim1     = static_cast<index_type>(points.dimension(0));
      pointDim = static_cast<index_type>(points.dimension(1));
      break;
    case 3:
      dim0     = static_cast<index_type>(points.dimension(0));
      dim1     = static_cast<index_type>(points.dimension(1));
      pointDim = static_cast<index_type>(points.dimension(2));
      break;
    default:
      INTREPID2_TEST_FOR_EXCEPTION( !( (1 <= getrank(points) ) && (getrank(points) <= 3) ), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): rank-1, 2 or 3 required for input points array. ");      
    }// switch
  
  
    // This method can handle up to rank-3 input arrays. The spatial dim must be the last dimension. 
    // The method uses [] accessor because array rank is determined at runtime and the appropriate
    // (i,j,..,k) accessor is not known. Use of [] requires the following offsets:
    //    for input array  = i0*dim1*pointDim + i1*dim1  (computed in 2 pieces: inPtr0 and inPtr1, resp)
    //    for output array = i0*dim1                     (computed in one piece: outPtr0)
    Scalar point[3] = {0.0, 0.0, 0.0};

    INTREPID2_TEST_FOR_EXCEPTION( !( (1 <= pointDim) && (pointDim <= 3)), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): Input array specifies invalid point dimension ");

    switch(apRank) {
    case 1:
      for(index_type i2 = 0; i2 < pointDim; i2++)
        point[i2] = points(i2);
      inRefCell(0) = checkPointInclusion(point, pointDim, cellTopo, threshold);
      break;
    case 2:
      for(index_type i1 = 0; i1 < dim1; i1++) {
        for(index_type i2 = 0; i2 < pointDim; i2++)
          point[i2] = points(i1,i2);
        inRefCell(i1) = checkPointInclusion(point, pointDim, cellTopo, threshold);
      }
      break;
    case 3:
      for(index_type i0 = 0; i0 < dim0; i0++){
        for(index_type i1 = 0; i1 < dim1; i1++) {
          for(index_type i2 = 0; i2 < pointDim; i2++)
            point[i2] = points(i0,i1,i2);
          inRefCell(i0,i1) = checkPointInclusion(point, pointDim, cellTopo, threshold);
        }
      }
      break;
    }
  }  


  template<class Scalar>
  template<class ArrayIncl, class ArrayPoint, class ArrayCell>
  void CellTools<Scalar>::checkPointwiseInclusion(ArrayIncl &                   inCell,
                                                  const ArrayPoint &            points,
                                                  const ArrayCell &             cellWorkset,
                                                  const shards::CellTopology &  cell,
                                                  const int &                   whichCell, 
                                                  const double &                threshold)
  {
    INTREPID2_VALIDATE( validateArguments_checkPointwiseInclusion(inCell, points, cellWorkset, whichCell, cell) );
  
    // For cell topologies with reference cells this test maps the points back to the reference cell
    // and uses the method for reference cells
    unsigned baseKey = cell.getBaseCellTopologyData() -> key;
  
    switch(baseKey){
    
    case shards::Line<>::key :
    case shards::Triangle<>::key:
    case shards::Quadrilateral<>::key :
    case shards::Tetrahedron<>::key :
    case shards::Hexahedron<>::key :
    case shards::Wedge<>::key :
    case shards::Pyramid<>::key :
      {
        FieldContainer<Scalar> refPoints;
        
        if(getrank(points) == 2){
          refPoints.resize(static_cast<index_type>(points.dimension(0)), static_cast<index_type>(points.dimension(1)) );
          mapToReferenceFrame(refPoints, points, cellWorkset, cell, whichCell);
          checkPointwiseInclusion(inCell, refPoints, cell, threshold );
        }
        else if(getrank(points) == 3){
          refPoints.resize(static_cast<index_type>(points.dimension(0)), static_cast<index_type>(points.dimension(1)), static_cast<index_type>(points.dimension(2)) );
          mapToReferenceFrame(refPoints, points, cellWorkset, cell, whichCell);
          checkPointwiseInclusion(inCell, refPoints, cell, threshold );          
        }
        break;
      }
    default: 
      INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): cell topology not supported");
    }// switch
  
  }

}

#endif
