// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file   Intrepid2_CellToolsDefValidateArguments.hpp
    \brief  Definition file for the validate arguments functions of the Intrepid2::CellTools class.
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
           typename PointViewType,
           typename worksetCellViewType>
  void 
  CellTools_setJacobianArgs( const jacobianViewType     jacobian,
                             const PointViewType        points,
                             const worksetCellViewType  worksetCell,
                             const shards::CellTopology cellTopo,
                             const int startCell, const int endCell) {
    // Validate worksetCell array
    INTREPID2_TEST_FOR_EXCEPTION( worksetCell.rank() != 3, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::setJacobian): rank = 3 required for worksetCell array." );
    //TODO: check this. not working for composite tet
    //INTREPID2_TEST_FOR_EXCEPTION( worksetCell.extent(1) != cellTopo.getSubcellCount(0), std::invalid_argument,
    //                              ">>> ERROR (Intrepid2::CellTools::setJacobian): dim 1 (number of cell nodes) of worksetCell array does not match cell topology." );
  
    INTREPID2_TEST_FOR_EXCEPTION( worksetCell.extent(2) != cellTopo.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::setJacobian): dim 2 (spatial dimension) of worksetCell array  does not match cell dimension." );
    
    // Validate points array: can be rank-2 (P,D) or rank-3 (C,P,D)
    // If rank-2: admissible jacobians: rank-3 (P,D,D) or rank-4 (C,P,D,D); admissible whichCell: -1 (default) or cell ordinal.
    const auto pointRank = points.rank();
    INTREPID2_TEST_FOR_EXCEPTION( pointRank != 2 &&
                                  pointRank != 3, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::setJacobian): points must have rank 2 or 3." );

    const int endCellResolved = (endCell == -1) ? worksetCell.extent_int(0) : endCell;
    const int numCells = endCellResolved - startCell;
    
    INTREPID2_TEST_FOR_EXCEPTION(startCell < 0, std::invalid_argument, "Invalid startCell");
    INTREPID2_TEST_FOR_EXCEPTION(startCell >= worksetCell.extent_int(0), std::invalid_argument, "startCell is out of bounds in workset.");
    INTREPID2_TEST_FOR_EXCEPTION(endCellResolved > worksetCell.extent_int(0), std::invalid_argument, "resolved endCell is out of bounds in workset.");
    
    switch (pointRank) {
    case 2: {
      INTREPID2_TEST_FOR_EXCEPTION( points.extent(1) != cellTopo.getDimension(), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::setJacobian): dim 1 (spatial dimension) of points array does not match cell dimension." );
      
      INTREPID2_TEST_FOR_EXCEPTION( jacobian.rank() != 4, std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::setJacobian): rank = 4 required for jacobian array." );
      
      INTREPID2_TEST_FOR_EXCEPTION( jacobian.extent_int(0) != numCells, std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::setJacobian): dim 0 (number of cells) of jacobian array must equal number of cells requested from in the workset." );
      
      INTREPID2_TEST_FOR_EXCEPTION( jacobian.extent(1) != points.extent(0), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::setJacobian): dim 1 (number of points) of jacobian array must equal dim 0 of points array." );
      
      INTREPID2_TEST_FOR_EXCEPTION( jacobian.extent(2) != points.extent(1), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::setJacobian): dim 2 (spatial dimension) of jacobian array must equal dim 1 of points array." );
      
      INTREPID2_TEST_FOR_EXCEPTION( jacobian.extent(2) != jacobian.extent(3), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::setJacobian): dim 2 = dim 3 (same spatial dimensions) required for jacobian array." );
      
      INTREPID2_TEST_FOR_EXCEPTION( jacobian.extent(3) < 1 || jacobian.extent(3) > 3, std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::setJacobian): dim 2 and dim 3 (spatial dimensions) must be between 1 and 3." );
      break;
    }
    case 3: {
      INTREPID2_TEST_FOR_EXCEPTION( points.extent_int(0) != numCells, std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::setJacobian): dim 0 (number of cells) of points array must equal number of cells requested from in the workset.");

      INTREPID2_TEST_FOR_EXCEPTION( points.extent(2) != cellTopo.getDimension(), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::setJacobian): dim 2 (spatial dimension) of points array does not match cell dimension");
      
      // rank-4 (C,P,D,D) jacobian required for rank-3 (C,P,D) input points
      INTREPID2_TEST_FOR_EXCEPTION( jacobian.rank() != 4, std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::setJacobian): rank = 4 required for jacobian array." );                                    

      INTREPID2_TEST_FOR_EXCEPTION( jacobian.extent(0) != points.extent(0), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::setJacobian): dim 0 (number of cells) of jacobian array must equal dim 0 of points array");

      INTREPID2_TEST_FOR_EXCEPTION( jacobian.extent(1) != points.extent(1), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::setJacobian): dim 1 (number of points) of jacobian array must equal dim 1 of points array");
      
      INTREPID2_TEST_FOR_EXCEPTION( jacobian.extent(2) != points.extent(2), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::setJacobian): dim 2 (spatial dimension) of jacobian array must equal dim 2 of points array");
      
      INTREPID2_TEST_FOR_EXCEPTION( jacobian.extent(2) != jacobian.extent(3), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::setJacobian): dim 2 = dim 3 (same spatial dimensions) required for jacobian array. ");

      INTREPID2_TEST_FOR_EXCEPTION( jacobian.extent(3) < 1 || jacobian.extent(3) > 3, std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::setJacobian): dim 2 and dim 3 (spatial dimensions) must be between 1 and 3." );
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
    const ordinal_type jacoRank = jacobian.rank();
    INTREPID2_TEST_FOR_EXCEPTION( jacoRank != 4 && 
                                  jacoRank != 3, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::setJacobianInv): rank = 4 or 3 required for jacobian array." );
  
    // Verify correctness of spatial dimensions - they are the last two dimensions of the array: rank-2 and rank-1
    INTREPID2_TEST_FOR_EXCEPTION( jacobian.extent(jacoRank - 1) != jacobian.extent(jacoRank - 2), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::setJacobianInv): dim(rank-2) = dim(rank-2) (same spatial dimensions) required for jacobian array." );
    
    INTREPID2_TEST_FOR_EXCEPTION( jacobian.extent(jacoRank - 1) < 1 || 
                                  jacobian.extent(jacoRank - 1) > 3, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::setJacobianInv): dim(rank-1) and dim(rank-2) (spatial dimensions) must be between 1 and 3." );
    
    // Validate output jacobianInv array: must have the same rank and dimensions as the input array.
    const ordinal_type jacoInvRank = jacobianInv.rank();
    INTREPID2_TEST_FOR_EXCEPTION( jacoInvRank != jacoRank, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::setJacobianInv): jacobian rank does not match to jacobianInv." );
  
    for (ordinal_type i=0;i<jacoRank;++i) {
      INTREPID2_TEST_FOR_EXCEPTION( jacobianInv.extent(i) != jacobian.extent(i), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::setJacobianInv): jacobian dimension (i) does not match to jacobianInv dimension (i)." );
    }
  }
  
  
  template<typename jacobianDetViewType, 
           typename jacobianViewType>
  void 
  CellTools_setJacobianDetArgs( const jacobianDetViewType jacobianDet,
                                const jacobianViewType    jacobian ) {
    // Validate input jacobian array: admissible ranks & dimensions are: 
    // - rank-4 with dimensions (C,P,D,D), or rank-3 with dimensions (P,D,D).
    const ordinal_type jacoRank = jacobian.rank();
    INTREPID2_TEST_FOR_EXCEPTION( jacoRank != 4 &&
                                  jacoRank != 3, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::setJacobianInv): rank = 4 or 3 required for jacobian array." );
  
    // Verify correctness of spatial dimensions - they are the last two dimensions of the array: rank-2 and rank-1
    INTREPID2_TEST_FOR_EXCEPTION( jacobian.extent(jacoRank - 1) != jacobian.extent(jacoRank - 2), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::setJacobianInv): dim(rank-2) = dim(rank-2) (same spatial dimensions) required for jacobian array." );
  
    INTREPID2_TEST_FOR_EXCEPTION( jacobian.extent(jacoRank - 1) < 1 || 
                                  jacobian.extent(jacoRank - 1) > 3, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::setJacobianInv): dim(rank-1) and dim(rank-2) (spatial dimensions) must be between 1 and 3." );
    
    // Validate output jacobianDet array
    const ordinal_type jacoDetRank = jacobianDet.rank();
    //  must be rank-2 with dimensions (C,P) if jacobian was rank-4
    // must be rank-1 with dimension (P) if jacobian was rank-3
    INTREPID2_TEST_FOR_EXCEPTION( jacoDetRank != (jacoRank-2), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::setJacobianDetArgs): rank = 2 required for jacobianDet if jacobian is rank-4." );
    
    for (ordinal_type i=0;i<jacoDetRank;++i) {
      INTREPID2_TEST_FOR_EXCEPTION( jacobianDet.extent(i) != jacobian.extent(i), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::setJacobianDetArgs): jacobianDet dimension (i) does not match to jacobian dimension (i)." );
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
    INTREPID2_TEST_FOR_EXCEPTION( (worksetCell.rank() != 3) || (physPoints.rank() != 3), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::mapToPhysicalFrame): rank = 3 required for worksetCell and physPoints arrays." );
 
    //TODO: check this, not working for tria6 
    //INTREPID2_TEST_FOR_EXCEPTION( worksetCell.extent(1) != cellTopo.getSubcellCount(0), std::invalid_argument,
    //                              ">>> ERROR (Intrepid2::CellTools::mapToPhysicalFrame): dim 1 (number of cell nodes) of worksetCell array does not match cell topology." );
  
    //we allow cells immersed in a higher-dimensional space  (e.g. 2d cell in a 3d space)
    INTREPID2_TEST_FOR_EXCEPTION( worksetCell.extent(2) < cellTopo.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::mapToPhysicalFrame): dim 2 (spatial dimension) of worksetCell array is smaller than the cell dimension." );

    INTREPID2_TEST_FOR_EXCEPTION( physPoints.extent(2) !=  worksetCell.extent(2), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::mapToPhysicalFrame): physPoints and worksetCell should have the same spatial dimension." );


    // Validate refPoints array: can be rank-2 (P,D) or rank-3 (C,P,D) array
    const ordinal_type refPointRank = refPoints.rank();
    const ordinal_type physPointRank = physPoints.rank();

    INTREPID2_TEST_FOR_EXCEPTION( refPointRank != 2 &&
                                  refPointRank != 3, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::mapToPhysicalFrame): refPoints requires rank 2 or 3." );
    
    switch (refPointRank) {
    case 2: {
      // If rank-2: admissible output array is (P,D) or (C,P,D)
      INTREPID2_TEST_FOR_EXCEPTION( refPoints.extent(1) != cellTopo.getDimension(), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::mapToPhysicalFrame): dim 1 (spatial dimension) of refPoints array does not match cell dimension." );
      
      INTREPID2_TEST_FOR_EXCEPTION( physPoints.rank() != 3, std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::mapToPhysicalFrame): rank = 3 required for physPoints array for the default whichCell value." );
      
      INTREPID2_TEST_FOR_EXCEPTION( physPoints.extent(0) != worksetCell.extent(0), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::mapToPhysicalFrame): dim 0 (number of cells) of physPoints array must equal dim 0 of worksetCell array." );
      
      INTREPID2_TEST_FOR_EXCEPTION( physPoints.extent(1) != refPoints.extent(0), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::mapToPhysicalFrame): dim 1 (number of points) of physPoints array must equal dim 0 of refPoints array." ); 
      
  
      break;
    }
    case 3: {
      // refPoints is (C,P,D): requires physPoints to be (C,P,D) and whichCell=-1  (because all cell mappings are applied)
      // validate refPoints dimensions and rank
      INTREPID2_TEST_FOR_EXCEPTION( refPoints.extent(0) != worksetCell.extent(0), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::mapToPhysicalFrame): dim 0 (number of cells) of refPoints and worksetCell arraya are required to match." );
      
      INTREPID2_TEST_FOR_EXCEPTION( refPoints.extent(2) != cellTopo.getDimension(), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::mapToPhysicalFrame): dim 2 (spatial dimension) of refPoints array does not match cell dimension." );
    
      // physPoints must match rank and dimensions of refPoints
      INTREPID2_TEST_FOR_EXCEPTION( refPointRank != physPointRank, std::invalid_argument, 
                                    " >>> ERROR (Intrepid2::CellTools::mapToPhysicalFrame): refPoints rank does not match to physPoints rank." );
      
      for (ordinal_type i=0;i<refPointRank-1;++i) {
        INTREPID2_TEST_FOR_EXCEPTION( refPoints.extent(i) != physPoints.extent(i), std::invalid_argument, 
                                      " >>> ERROR (Intrepid2::CellTools::mapToPhysicalFrame): refPoints dimension(i) does not match to physPoints dimension(i)." );
      }
      break;
    }
    }
  }

  template<typename refPointViewType, 
           typename physPointViewType,
           typename worksetCellViewType>
  void 
  CellTools_mapToReferenceFrameArgs( const refPointViewType     refPoints,
                                     const physPointViewType    physPoints,
                                     const worksetCellViewType  worksetCell,
                                     const shards::CellTopology cellTopo ) {
    // Validate worksetCell array
    const ordinal_type worksetCellRank = worksetCell.rank();
    INTREPID2_TEST_FOR_EXCEPTION( worksetCellRank != 3, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::mapToReferenceFrame): rank = 3 required for worksetCell array" );
    // TODO: check this. 
    // INTREPID2_TEST_FOR_EXCEPTION( worksetCell.extent(1) != cellTopo.getSubcellCount(0), std::invalid_argument,
    //                              ">>> ERROR (Intrepid2::CellTools::mapToReferenceFrame): dim 1 (number of cell nodes) of worksetCell array does not match cell topology" );
  
    INTREPID2_TEST_FOR_EXCEPTION( worksetCell.extent(2) != cellTopo.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::mapToReferenceFrame): dim 2 (spatial dimension) of worksetCell array  does not match cell dimension" );
    
    // Admissible ranks and dimensions of refPoints and physPoints depend on whichCell value:
    // default is to map multiple sets of points to multiple sets of points. (C,P,D) arrays required
    
    const ordinal_type physPointRank = physPoints.rank();
    const ordinal_type refPointRank = refPoints.rank();
    
    INTREPID2_TEST_FOR_EXCEPTION( refPointRank != 2 &&
                                  refPointRank != 3, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::mapToReferenceFrame): refPoint must have rank 2 or 3." );
    
    INTREPID2_TEST_FOR_EXCEPTION( physPointRank != refPointRank, std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::mapToReferenceFrame): physPoints rank does not match refPoints rank." );
    for (ordinal_type i=0;i<refPointRank;++i) {
      INTREPID2_TEST_FOR_EXCEPTION( refPoints.extent(i) != physPoints.extent(i), std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::mapToReferenceFrame): physPoints dimension (i) does not match refPoints dimension (i)." );
    }
  }

  template<typename refPointViewType, 
           typename initGuessViewType,
           typename physPointViewType,
           typename worksetCellViewType>
  void CellTools_mapToReferenceFrameInitGuessArgs( const refPointViewType     refPoints,
                                                   const initGuessViewType    initGuess,
                                                   const physPointViewType    physPoints,
                                                   const worksetCellViewType  worksetCell,
                                                   const shards::CellTopology cellTopo ) {
    // Call the method that validates arguments with the default initial guess selection
    CellTools_mapToReferenceFrameArgs(refPoints, physPoints, worksetCell, cellTopo);
  
    // Then check initGuess: its rank and dimensions must match those of physPoints.
    INTREPID2_TEST_FOR_EXCEPTION( initGuess.rank() != physPoints.rank(), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::mapToReferenceFrame): InitGuess must have the same rank as physPoints");         

    const ordinal_type r = initGuess.rank();
    for (ordinal_type i=0;i<r;++i) {
      INTREPID2_TEST_FOR_EXCEPTION( initGuess.extent(i) != physPoints.extent(i), std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::mapToReferenceFrame): InitGuess dimension (i) does not match ot physPoints dimension(i).");
    }         
  }


} // end of intrepid2

#endif





// template<class Scalar>
// template<class ArrayIncl, class ArrayPoint, class ArrayCell>
// void CellTools<Scalar>::checkPointwiseInclusion(ArrayIncl &                   inCell,
//                                                                   const ArrayPoint &            physPoints,
//                                                                   const ArrayCell &             worksetCell,
//                                                                   const ordinal_type &                   whichCell,
//                                                                   const shards::CellTopology &  cell)
// {
//   // Validate worksetCell array
//   INTREPID2_TEST_FOR_EXCEPTION( (getrank(worksetCell) != 3), std::invalid_argument,
//                                 ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): rank = 3 required for worksetCell array" );
  
//   INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(worksetCell.extent(1)) != (index_type)cell.getSubcellCount(0) ), std::invalid_argument,
//                                 ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): dim 1 (number of cell nodes) of worksetCell array does not match cell topology" );
  
//   INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(worksetCell.extent(2)) != (index_type)cell.getDimension() ), std::invalid_argument,
//                                 ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): dim 2 (spatial dimension) of worksetCell array  does not match cell dimension" );
  
  
//   // Validate whichCell It can be either -1 (default value) or a valid cell ordinal.
//   INTREPID2_TEST_FOR_EXCEPTION( !( ( (0 <= whichCell ) && (whichCell < worksetCell.extent(0) ) ) || (whichCell == -1) ), std::invalid_argument,
//                                 ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): whichCell = -1 or a valid cell ordinal is required." );  
  
//   // Validate points array: can be rank-2 (P,D) or rank-3 (C,P,D)
//   // If rank-2: admissible inCell is rank-1 (P); admissible whichCell is valid cell ordinal but not -1.
//   if(getrank(physPoints) == 2) {
    
//     INTREPID2_TEST_FOR_EXCEPTION( (whichCell == -1), std::invalid_argument,
//                                   ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): whichCell = a valid cell ordinal is required with rank-2 input array." );

//     INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(physPoints.extent(1)) != (index_type)cell.getDimension() ), std::invalid_argument,
//                                   ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): dim 1 (spatial dimension) of physPoints array does not match cell dimension" );
    
//     // Validate inCell
//     INTREPID2_TEST_FOR_EXCEPTION( (getrank(inCell) != 1), std::invalid_argument, 
//                                   ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): rank = 1 required for inCell array" );
    
//     INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(inCell.extent(0)) != static_cast<index_type>(physPoints.extent(0))), std::invalid_argument,
//                                   ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): dim 0 (number of points) of inCell array must equal dim 0 of physPoints array" );
//   }
//   // If rank-3: admissible inCell is rank-2 (C,P); admissible whichCell = -1.
//   else if (getrank(physPoints) == 3){
    
//     INTREPID2_TEST_FOR_EXCEPTION( !(whichCell == -1), std::invalid_argument,
//                                   ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): whichCell = -1 is required with rank-3 input array." );
    
//     INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(physPoints.extent(0)) != static_cast<index_type>(worksetCell.extent(0)) ), std::invalid_argument,
//                                   ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): dim 0 (number of cells)  of physPoints array must equal dim 0 of worksetCell array " );

//     INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(physPoints.extent(2)) != (index_type)cell.getDimension() ), std::invalid_argument,
//                                   ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): dim 2 (spatial dimension) of physPoints array does not match cell dimension" );
    
//     // Validate inCell
//     INTREPID2_TEST_FOR_EXCEPTION( (getrank(inCell) != 2), std::invalid_argument, 
//                                   ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): rank = 2 required for inCell array" );
    
//     INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(inCell.extent(0)) != static_cast<index_type>(physPoints.extent(0))), std::invalid_argument,
//                                   ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): dim 0 (number of cells) of inCell array must equal dim 0 of physPoints array" );    

//     INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(inCell.extent(1)) != static_cast<index_type>(physPoints.extent(1))), std::invalid_argument,
//                                   ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): dim 1 (number of points) of inCell array must equal dim 1 of physPoints array" );    
//   }
//   else {
//     INTREPID2_TEST_FOR_EXCEPTION( !( (getrank(physPoints) == 2) && (getrank(physPoints) ==3) ), std::invalid_argument,
//                                   ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): rank = 2 or 3 required for points array" );
//   }
// }
