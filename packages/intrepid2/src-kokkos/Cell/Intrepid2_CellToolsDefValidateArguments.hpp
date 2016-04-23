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
#ifndef __INTREPID2_CELLTOOLS_DEF_VALIDATE_ARGUMENTS_HPP__
#define __INTREPID2_CELLTOOLS_DEF_VALIDATE_ARGUMENTS_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {

  //============================================================================================//
  //                                                                                            //
  //                  Validation of input/output arguments for CellTools methods                //
  //                                                                                            //
  //============================================================================================//

  template<typename jacobianViewType, 
           typename pointViewType,
           typename worksetCellViewType>
  void 
  CellTools_setJacobianArgs( const jacobianViewType     jacobian,
                             const pointViewType        points,
                             const worksetCellViewType  worksetCell,
                             const shards::CellTopology cellTopo ) {
    // Validate worksetCell array
    INTREPID2_TEST_FOR_EXCEPTION( worksetCell.rank() != 3, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): rank = 3 required for worksetCell array." );
  
    INTREPID2_TEST_FOR_EXCEPTION( worksetCell.dimension(1) != cellTopo.getSubcellCount(0), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 1 (number of cell nodes) of worksetCell array does not match cell topology." );
  
    INTREPID2_TEST_FOR_EXCEPTION( worksetCell.dimension(2) != cellTopo.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 2 (spatial dimension) of worksetCell array  does not match cell dimension." );
    
    // Validate points array: can be rank-2 (P,D) or rank-3 (C,P,D)
    // If rank-2: admissible jacobians: rank-3 (P,D,D) or rank-4 (C,P,D,D); admissible whichCell: -1 (default) or cell ordinal.
    const auto pointRank = points.rank();
    INTREPID2_TEST_FOR_EXCEPTION( pointRank != 2 &&
                                  pointRank != 3, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): points must have rank 2 or 3." );

    switch (pointRank) {
    case 2: {
      INTREPID2_TEST_FOR_EXCEPTION( points.dimension(1) != cellTopo.getDimension(), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 1 (spatial dimension) of points array does not match cell dimension." );
      
      INTREPID2_TEST_FOR_EXCEPTION( jacobian.rank() != 4, std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): rank = 4 required for jacobian array." );
      
      INTREPID2_TEST_FOR_EXCEPTION( jacobian.dimension(0) != worksetCell.dimension(0), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 0 (number of cells) of jacobian array must equal dim 0 of worksetCell array." );
      
      INTREPID2_TEST_FOR_EXCEPTION( jacobian.dimension(1) != points.dimension(0), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 1 (number of points) of jacobian array must equal dim 0 of points array." );
      
      INTREPID2_TEST_FOR_EXCEPTION( jacobian.dimension(2) != points.dimension(1), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 2 (spatial dimension) of jacobian array must equal dim 1 of points array." );
      
      INTREPID2_TEST_FOR_EXCEPTION( jacobian.dimension(2) != jacobian.dimension(3), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 2 = dim 3 (same spatial dimensions) required for jacobian array." );
      
      INTREPID2_TEST_FOR_EXCEPTION( jacobian.dimension(3) < 1 || jacobian.dimension(3) > 3, std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 2 and dim 3 (spatial dimensions) must be between 1 and 3." );
      break;
    }
    case 3: {
      INTREPID2_TEST_FOR_EXCEPTION( points.dimension(0) != cellWorkset.dimension(0), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 0 (number of cells) of points array must equal dim 0 of cellWorkset array");

      INTREPID2_TEST_FOR_EXCEPTION( points.dimension(2) != cellTopo.getDimension(), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 2 (spatial dimension) of points array does not match cell dimension");
      
      // rank-4 (C,P,D,D) jacobian required for rank-3 (C,P,D) input points
      INTREPID2_TEST_FOR_EXCEPTION( jacobian.rank() != 4, std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): rank = 4 required for jacobian array." );                                    

      INTREPID2_TEST_FOR_EXCEPTION( jacobian.dimension(0) != points.dimension(0), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 0 (number of cells) of jacobian array must equal dim 0 of points array");

      INTREPID2_TEST_FOR_EXCEPTION( jacobian.dimension(1) != points.dimension(1), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 1 (number of points) of jacobian array must equal dim 1 of points array");
      
      INTREPID2_TEST_FOR_EXCEPTION( jacobian.dimension(2) != points.dimension(2), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 2 (spatial dimension) of jacobian array must equal dim 2 of points array");
      
      INTREPID2_TEST_FOR_EXCEPTION( jacobian.dimension(2) != jacobian.dimension(3), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 2 = dim 3 (same spatial dimensions) required for jacobian array. ");

      INTREPID2_TEST_FOR_EXCEPTION( jacobian.dimension(3) < 1 || jacobian.dimension(3) > 3, std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 2 and dim 3 (spatial dimensions) must be between 1 and 3." );
      break;
    }
    }
  }

  template<typename jacobianInvViewType, 
           typename jacobianViewType>
  void 
  CellTools_setJacobianInvArgs( const jacobianInvViewType jacobianInv,
                                const jacobianViewType    jacobian ) {
    // Validate input jacobian array: admissible ranks & dimensions are: 
    // - rank-4 with dimensions (C,P,D,D), or rank-3 with dimensions (P,D,D).
    const auto jacobRank = jacobian.rank();
    INTREPID2_TEST_FOR_EXCEPTION( jacobRank != 4 && 
                                  jacobRank != 3, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobianInv): rank = 4 or 3 required for jacobian array." );
  
    // Verify correctness of spatial dimensions - they are the last two dimensions of the array: rank-2 and rank-1
    INTREPID2_TEST_FOR_EXCEPTION( jacobian.dimension(jacobRank - 1) != jacobian.dimension(jacobRank - 2), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobianInv): dim(rank-2) = dim(rank-2) (same spatial dimensions) required for jacobian array." );
    
    INTREPID2_TEST_FOR_EXCEPTION( jacobian.dimension(jacobRank - 1) < 1 || 
                                  jacobian.dimension(jacobRank - 1) > 3, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobianInv): dim(rank-1) and dim(rank-2) (spatial dimensions) must be between 1 and 3." );
    
    // Validate output jacobianInv array: must have the same rank and dimensions as the input array.
    const auto jacoInvRank = jacobianInv.rank();
    INTREPID2_TEST_FOR_EXCEPTION( jacoInvRank != jacoRank, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobianInv): jacobian rank does not match to jacobianInv." );
  
    for (auto i=0;i<jacoRank;++i) {
      INTREPID2_TEST_FOR_EXCEPTION( jacobianInv.dimension(i) != jacobian.dimension(i), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobianInv): jacobian dimension (i) does not match to jacobianInv dimension (i)." );
    }
  }
  
  
  template<typename jacobianDetViewType, 
           typename jacobianViewType>
  void 
  CellTools_setJacobianDetArgs( const jacobianDetViewType jacobianDet,
                                const jacobianViewType    jacobian ) {
    // Validate input jacobian array: admissible ranks & dimensions are: 
    // - rank-4 with dimensions (C,P,D,D), or rank-3 with dimensions (P,D,D).
    const auto jacobRank = jacobian.rank();
    INTREPID2_TEST_FOR_EXCEPTION( jacobRank != 4 &&
                                  jacobRank != 3, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobianInv): rank = 4 or 3 required for jacobian array." );
  
    // Verify correctness of spatial dimensions - they are the last two dimensions of the array: rank-2 and rank-1
    INTREPID2_TEST_FOR_EXCEPTION( jacobian.dimension(jacobRank - 1) != jacobian.dimension(jacobRank - 2), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobianInv): dim(rank-2) = dim(rank-2) (same spatial dimensions) required for jacobian array." );
  
    INTREPID2_TEST_FOR_EXCEPTION( jacobian.dimension(jacobRank - 1) < 1 || 
                                  jacobian.dimension(jacobRank - 1) > 3, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobianInv): dim(rank-1) and dim(rank-2) (spatial dimensions) must be between 1 and 3." );
    
    // Validate output jacobianDet array
    const auto jacoDetRank = jacobianDet.rank();
    //  must be rank-2 with dimensions (C,P) if jacobian was rank-4
    // must be rank-1 with dimension (P) if jacobian was rank-3
    INTREPID2_TEST_FOR_EXCEPTION( jacoDetRank != (jacoRank-2), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobianDetArgs): rank = 2 required for jacobianDet if jacobian is rank-4." );
    
    for (auto i=0;i<jacoDetRank;++i) {
      INTREPID2_TEST_FOR_EXCEPTION( jacobianDet.dimension(i) != jacobian.dimension(i), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobianDetArgs): jacobianDet dimension (i) does not match to jacobian dimension (i)." );
    }
  }



  template<typename physPointViewType,
           typename refPointViewType,
           typename worksetCellViewType>
  void 
  CellTools_mapToPhysicalFrameArgs( const physPointViewType    physPoints,
                                    const refPointViewType     refPoints,
                                    const worksetCellViewType  worksetCell,
                                    const shards::CellTopology cellTopo ) {
    // Validate worksetCell array
    INTREPID2_TEST_FOR_EXCEPTION( worksetCell.rank() != 3, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): rank = 3 required for worksetCell array." );
  
    INTREPID2_TEST_FOR_EXCEPTION( worksetCell.dimension(1) != cellTopo.getSubcellCount(0), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): dim 1 (number of cell nodes) of worksetCell array does not match cell topology." );
  
    INTREPID2_TEST_FOR_EXCEPTION( worksetCell.dimension(2) != cellTopo.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): dim 2 (spatial dimension) of worksetCell array  does not match cell dimension." );


    // Validate refPoints array: can be rank-2 (P,D) or rank-3 (C,P,D) array
    const auto refPointRank = refPoints.rank();
    const auto physPointRank = physPoints.rank();

    INTREPID2_TEST_FOR_EXCEPTION( refPointRank != 2 &&
                                  refPointRank != 3, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): refPoints requires rank 2 or 3." );
    
    switch (refPointRank) {
    case 2: {
      // If rank-2: admissible output array is (P,D) or (C,P,D)
      INTREPID2_TEST_FOR_EXCEPTION( refPoints.dimension(1) != cellTopo.getDimension(), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): dim 1 (spatial dimension) of refPoints array does not match cell dimension." );
      
      INTREPID2_TEST_FOR_EXCEPTION( physPoints.rank() != 3, std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): rank = 3 required for physPoints array for the default whichCell value." );
      
      INTREPID2_TEST_FOR_EXCEPTION( physPoints.dimension(0) != worksetCell.dimension(0), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): dim 0 (number of cells) of physPoints array must equal dim 0 of worksetCell array." );
      
      INTREPID2_TEST_FOR_EXCEPTION( physPoints.dimension(1) != refPoints.dimension(0), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): dim 1 (number of points) of physPoints array must equal dim 0 of refPoints array." ); 
      
      INTREPID2_TEST_FOR_EXCEPTION( physPoints.dimension(2) != cellTopo.getDimension(), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): dim 2 (spatial dimension) does not match cell dimension." );  
      break;
    }
    case 3: {
      // refPoints is (C,P,D): requires physPoints to be (C,P,D) and whichCell=-1  (because all cell mappings are applied)
      // validate refPoints dimensions and rank
      INTREPID2_TEST_FOR_EXCEPTION( refPoints.dimension(0) != worksetCell.dimension(0), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): dim 0 (number of cells) of refPoints and worksetCell arraya are required to match." );
      
      INTREPID2_TEST_FOR_EXCEPTION( refPoints.dimension(2) != cellTopo.getDimension(), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): dim 2 (spatial dimension) of refPoints array does not match cell dimension." );
    
      // physPoints must match rank and dimensions of refPoints
      INTREPID2_TEST_FOR_EXCEPTION( refPointRank != physPointRank, std::invalid_argument, 
                                    " >>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): refPoints rank does not match to physPoints rank." );
      
      for (auto i=0;i<refPointRank;++i) {
        INTREPID2_TEST_FOR_EXCEPTION( refPoints.dimension(i) != physPoints.dimension(i), std::invalid_argument, 
                                      " >>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): refPoints dimension(i) does not match to physPoints dimension(i)." );
      }
      break;
    }
    }
  }
}

#endif




// template<typename refPointViewType, 
//          typename physPointViewType,
//          typename worksetCellViewType>
// void 
// CellTools_mapToReferenceFrameArgs( const refPointViewType     refPoints,
//                                    const physPointViewType    physPoints,
//                                    const worksetCellViewType  worksetCell,
//                                    const shards::CellTopology cellTopo ) {
//   // Validate worksetCell array
//   const auto worksetCellRank = worksetCell.rank();
//   INTREPID2_TEST_FOR_EXCEPTION( worksetCellRank != 3, std::invalid_argument,
//                                 ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToReferenceFrame): rank = 3 required for worksetCell array" );
    
//   INTREPID2_TEST_FOR_EXCEPTION( worksetCell.dimension(1) != cellTopo.getSubcellCount(0), std::invalid_argument,
//                                 ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToReferenceFrame): dim 1 (number of cell nodes) of worksetCell array does not match cell topology" );
  
//   INTREPID2_TEST_FOR_EXCEPTION( worksetCell.dimension(2) != cellTopo.getDimension(), std::invalid_argument,
//                                 ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToReferenceFrame): dim 2 (spatial dimension) of worksetCell array  does not match cell dimension" );
    
//   // Admissible ranks and dimensions of refPoints and physPoints depend on whichCell value:
//   // default is to map multiple sets of points to multiple sets of points. (C,P,D) arrays required

//   const auto physPointRank = physPoints.rank();
//   const auto refPointRank = refPoints.rank();

//   INTREPID2_TEST_FOR_EXCEPTION( refPointRank != 2 &&
//                                 refPointRank != 3, std::invalid_argument, 
//                                 ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToReferenceFrame): refPoint must have rank 2 or 3." );
    
//   INTREPID2_TEST_FOR_EXCEPTION( physPointRank != refPointRank, std::invalid_argument, 
//                                 ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToReferenceFrame): physPoints rank does not match refPoints rank." );
//   for (auto i=0;i<refPointRank;++i) {
//     INTREPID2_TEST_FOR_EXCEPTION( refPoints.dimension(i) != physPoints.dimension(i), std::invalid_argument, 
//                                   ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToReferenceFrame): physPoints dimension (i) does not match refPoints dimension (i)." );
//   }
// }



// template<class ArrayRefPoint, class ArrayInitGuess, class ArrayPhysPoint, class ArrayCell>
// void CellTools_mapToReferenceFrameArgs(const ArrayRefPoint  &        refPoints,
//                                        const ArrayInitGuess &        initGuess,
//                                        const ArrayPhysPoint &        physPoints,
//                                        const ArrayCell      &        worksetCell,
//                                        const shards::CellTopology &  cellTopo,
//                                        const int&                    whichCell)
// {
//   // Call the method that validates arguments with the default initial guess selection
//   validateArguments_mapToReferenceFrame(refPoints, physPoints, worksetCell, cellTopo, whichCell);
  
//   // Then check initGuess: its rank and dimensions must match those of physPoints.
//   std::string errmsg = ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToReferenceFrame):";
//   INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, initGuess, physPoints), std::invalid_argument, errmsg);  
// }


// template<class Scalar>
// template<class ArrayIncl, class ArrayPoint, class ArrayCell>
// void CellTools<Scalar>::validateArguments_checkPointwiseInclusion(ArrayIncl &                   inCell,
//                                                                   const ArrayPoint &            physPoints,
//                                                                   const ArrayCell &             worksetCell,
//                                                                   const int &                   whichCell,
//                                                                   const shards::CellTopology &  cell)
// {
//   // Validate worksetCell array
//   INTREPID2_TEST_FOR_EXCEPTION( (getrank(worksetCell) != 3), std::invalid_argument,
//                                 ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): rank = 3 required for worksetCell array" );
  
//   INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(worksetCell.dimension(1)) != (index_type)cell.getSubcellCount(0) ), std::invalid_argument,
//                                 ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): dim 1 (number of cell nodes) of worksetCell array does not match cell topology" );
  
//   INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(worksetCell.dimension(2)) != (index_type)cell.getDimension() ), std::invalid_argument,
//                                 ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): dim 2 (spatial dimension) of worksetCell array  does not match cell dimension" );
  
  
//   // Validate whichCell It can be either -1 (default value) or a valid cell ordinal.
//   INTREPID2_TEST_FOR_EXCEPTION( !( ( (0 <= whichCell ) && (whichCell < worksetCell.dimension(0) ) ) || (whichCell == -1) ), std::invalid_argument,
//                                 ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): whichCell = -1 or a valid cell ordinal is required." );  
  
//   // Validate points array: can be rank-2 (P,D) or rank-3 (C,P,D)
//   // If rank-2: admissible inCell is rank-1 (P); admissible whichCell is valid cell ordinal but not -1.
//   if(getrank(physPoints) == 2) {
    
//     INTREPID2_TEST_FOR_EXCEPTION( (whichCell == -1), std::invalid_argument,
//                                   ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): whichCell = a valid cell ordinal is required with rank-2 input array." );

//     INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(physPoints.dimension(1)) != (index_type)cell.getDimension() ), std::invalid_argument,
//                                   ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): dim 1 (spatial dimension) of physPoints array does not match cell dimension" );
    
//     // Validate inCell
//     INTREPID2_TEST_FOR_EXCEPTION( (getrank(inCell) != 1), std::invalid_argument, 
//                                   ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): rank = 1 required for inCell array" );
    
//     INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(inCell.dimension(0)) != static_cast<index_type>(physPoints.dimension(0))), std::invalid_argument,
//                                   ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): dim 0 (number of points) of inCell array must equal dim 0 of physPoints array" );
//   }
//   // If rank-3: admissible inCell is rank-2 (C,P); admissible whichCell = -1.
//   else if (getrank(physPoints) == 3){
    
//     INTREPID2_TEST_FOR_EXCEPTION( !(whichCell == -1), std::invalid_argument,
//                                   ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): whichCell = -1 is required with rank-3 input array." );
    
//     INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(physPoints.dimension(0)) != static_cast<index_type>(worksetCell.dimension(0)) ), std::invalid_argument,
//                                   ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): dim 0 (number of cells)  of physPoints array must equal dim 0 of worksetCell array " );

//     INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(physPoints.dimension(2)) != (index_type)cell.getDimension() ), std::invalid_argument,
//                                   ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): dim 2 (spatial dimension) of physPoints array does not match cell dimension" );
    
//     // Validate inCell
//     INTREPID2_TEST_FOR_EXCEPTION( (getrank(inCell) != 2), std::invalid_argument, 
//                                   ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): rank = 2 required for inCell array" );
    
//     INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(inCell.dimension(0)) != static_cast<index_type>(physPoints.dimension(0))), std::invalid_argument,
//                                   ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): dim 0 (number of cells) of inCell array must equal dim 0 of physPoints array" );    

//     INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(inCell.dimension(1)) != static_cast<index_type>(physPoints.dimension(1))), std::invalid_argument,
//                                   ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): dim 1 (number of points) of inCell array must equal dim 1 of physPoints array" );    
//   }
//   else {
//     INTREPID2_TEST_FOR_EXCEPTION( !( (getrank(physPoints) == 2) && (getrank(physPoints) ==3) ), std::invalid_argument,
//                                   ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): rank = 2 or 3 required for points array" );
//   }
// }
