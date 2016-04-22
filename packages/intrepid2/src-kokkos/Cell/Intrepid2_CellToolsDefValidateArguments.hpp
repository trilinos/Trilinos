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

  template<typename SpT>
  template<typename jacobianViewType, 
           typename pointViewType,
           typename worksetCellViewType>
  void CellTools<Scalar>::validateArguments_setJacobian( const jacobianViewType      jacobian,
                                                         const pointViewType         points,
                                                         const worksetCellViewType   worksetCell,
                                                         const shards::CellTopology  cellTopo ){
 
    // Validate worksetCell array
    INTREPID2_TEST_FOR_EXCEPTION( worksetCell.rank() != 3, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): rank = 3 required for worksetCell array");
  
    INTREPID2_TEST_FOR_EXCEPTION( worksetCell.dimension(0) <= 0, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 0 (number of cells) >= 1 required for worksetCell array");
  
    INTREPID2_TEST_FOR_EXCEPTION( worksetCell.dimension(1) != cellTopo.getSubcellCount(0), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 1 (number of cell nodes) of worksetCell array does not match cell topology");
  
    INTREPID2_TEST_FOR_EXCEPTION( worksetCell.dimension(2) != cellTopo.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 2 (spatial dimension) of worksetCell array  does not match cell dimension");
    
    // validate whichCell. It can be either -1 (default value) or a valid cell ordinal.
    //INTREPID2_TEST_FOR_EXCEPTION( !( ( (0 <= whichCell ) && (static_cast<index_type>(whichCell) < static_cast<index_type>(worksetCell.dimension(0)) ) ) || (whichCell == -1) ), std::invalid_argument,
    //                              ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): whichCell = -1 or a valid cell ordinal is required.");
  
  
    // Validate points array: can be rank-2 (P,D) or rank-3 (C,P,D)
    // If rank-2: admissible jacobians: rank-3 (P,D,D) or rank-4 (C,P,D,D); admissible whichCell: -1 (default) or cell ordinal.
    if (points.rank() == 2) {
      INTREPID2_TEST_FOR_EXCEPTION( points.dimension(0) <= 0, std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 0 (number of points) >= 1 required for points array ");
    
      INTREPID2_TEST_FOR_EXCEPTION( points.dimension(1) != cellTopo.getDimension(), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 1 (spatial dimension) of points array does not match cell dimension");
    
      // Validate the output array for the Jacobian: if whichCell == -1 all Jacobians are computed, rank-4 (C,P,D,D) required
      //if (whichCell == -1) 
      {
        INTREPID2_TEST_FOR_EXCEPTION( jacobian.rank() != 4, std::invalid_argument, 
                                      ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): rank = 4 required for jacobian array");
        
        INTREPID2_TEST_FOR_EXCEPTION( jacobian.dimension(0) != worksetCell.dimension(0), std::invalid_argument,
                                      ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 0 (number of cells) of jacobian array must equal dim 0 of worksetCell array");
      
        INTREPID2_TEST_FOR_EXCEPTION( jacobian.dimension(1) != points.dimension(0), std::invalid_argument,
                                      ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 1 (number of points) of jacobian array must equal dim 0 of points array");

        INTREPID2_TEST_FOR_EXCEPTION( jacobian.dimension(2) != points.dimension(1), std::invalid_argument,
                                      ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 2 (spatial dimension) of jacobian array must equal dim 1 of points array");
      
        INTREPID2_TEST_FOR_EXCEPTION( jacobian.dimension(2) != jacobian.dimension(3), std::invalid_argument,
                                      ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 2 = dim 3 (same spatial dimensions) required for jacobian array. ");
      
        INTREPID2_TEST_FOR_EXCEPTION( jacobian.dimension(3) < 1 || jacobian.dimension(3) > 3, std::invalid_argument,
                                      ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 2 and dim 3 (spatial dimensions) must be between 1 and 3. ");
      }     
      // // A single cell Jacobian is computed when whichCell != -1 (whichCell has been already validated), rank-3 (P,D,D) required
      // else {
      //   INTREPID2_TEST_FOR_EXCEPTION( (getrank(jacobian) != 3), std::invalid_argument, 
      //                                 ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): rank = 3 required for jacobian array");
      
      //   INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(jacobian.dimension(0)) != static_cast<index_type>(points.dimension(0))), std::invalid_argument,
      //                                 ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 0 (number of points) of jacobian array must equal dim 0 of points array");

      //   INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(jacobian.dimension(1)) != static_cast<index_type>(points.dimension(1))), std::invalid_argument,
      //                                 ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 1 (spatial dimension) of jacobian array must equal dim 1 of points array");
      
      //   INTREPID2_TEST_FOR_EXCEPTION( !(static_cast<index_type>(jacobian.dimension(1)) == static_cast<index_type>(jacobian.dimension(2)) ), std::invalid_argument,
      //                                 ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 1 = dim 2 (same spatial dimensions) required for jacobian array. ");
      
      //   INTREPID2_TEST_FOR_EXCEPTION( !( (0 < static_cast<index_type>(jacobian.dimension(1)) ) && (static_cast<index_type>(jacobian.dimension(1)) < 4) ), std::invalid_argument,
      //                                 ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 1 and dim 2 (spatial dimensions) must be between 1 and 3. ");
      // }
    }
    // Point array is rank-3 (C,P,D): requires whichCell = -1 and rank-4 (C,P,D,D) jacobians
    else if(getrank(points) ==3){
      std::string errmsg  = ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian):";
      INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(points.dimension(0)) != static_cast<index_type>(worksetCell.dimension(0)) ), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 0 (number of cells) of points array must equal dim 0 of worksetCell array");

      INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(points.dimension(1)) <= 0), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 1 (number of points) >= 1 required for points array ");
    
      INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(points.dimension(2)) != (index_type)cellTopo.getDimension() ), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 2 (spatial dimension) of points array does not match cell dimension");
    
      INTREPID2_TEST_FOR_EXCEPTION( (whichCell != -1), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): default value whichCell=-1 required for rank-3 input points");
    
      // rank-4 (C,P,D,D) jacobian required for rank-3 (C,P,D) input points
      INTREPID2_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, jacobian,  4, 4), std::invalid_argument,errmsg);
    
      INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(jacobian.dimension(0)) != static_cast<index_type>(points.dimension(0))), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 0 (number of cells) of jacobian array must equal dim 0 of points array");
    
      INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(jacobian.dimension(1)) != static_cast<index_type>(points.dimension(1))), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 1 (number of points) of jacobian array must equal dim 1 of points array");
  
      INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(jacobian.dimension(2)) != static_cast<index_type>(points.dimension(2))), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 2 (spatial dimension) of jacobian array must equal dim 2 of points array");
    
      INTREPID2_TEST_FOR_EXCEPTION( !(static_cast<index_type>(jacobian.dimension(2)) == static_cast<index_type>(jacobian.dimension(3)) ), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 2 = dim 3 (same spatial dimensions) required for jacobian array. ");
    
      INTREPID2_TEST_FOR_EXCEPTION( !( (0 < static_cast<index_type>(jacobian.dimension(3)) ) && (static_cast<index_type>(jacobian.dimension(3)) < 4) ), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): dim 2 and dim 3 (spatial dimensions) must be between 1 and 3. ");
    }
    else {
      INTREPID2_TEST_FOR_EXCEPTION( !( (getrank(points) == 2) && (getrank(points) ==3) ), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobian): rank = 2 or 3 required for points array");
    }  
  }



  template<class Scalar>
  template<class ArrayJacInv, class ArrayJac>
  void CellTools<Scalar>::validateArguments_setJacobianInv(const ArrayJacInv & jacobianInv,
                                                           const ArrayJac &    jacobian)
  {
    // Validate input jacobian array: admissible ranks & dimensions are: 
    // - rank-4 with dimensions (C,P,D,D), or rank-3 with dimensions (P,D,D).
    if(!CheckType<ArrayJac>::value || !CheckType<ArrayJacInv>::value ){
      std::cout <<std::endl<<"WARNING:: A Nonsupported Container is Being used with Intrepid2::CellTools<Scalar>::setJacobianInv"<<std::endl;
    }
    int jacobRank = getrank(jacobian);
    INTREPID2_TEST_FOR_EXCEPTION( !( (jacobRank == 4) || (jacobRank == 3) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobianInv): rank = 4 or 3 required for jacobian array. ");
  
    // Verify correctness of spatial dimensions - they are the last two dimensions of the array: rank-2 and rank-1
    INTREPID2_TEST_FOR_EXCEPTION( !(jacobian.dimension(jacobRank - 1) == jacobian.dimension(jacobRank - 2) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobianInv): dim(rank-2) = dim(rank-2) (same spatial dimensions) required for jacobian array. ");
  
    INTREPID2_TEST_FOR_EXCEPTION( !( (0 < jacobian.dimension(jacobRank - 1) ) && (jacobian.dimension(jacobRank - 1) < 4) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobianInv): dim(rank-1) and dim(rank-2) (spatial dimensions) must be between 1 and 3. ");
  
    // Validate output jacobianInv array: must have the same rank and dimensions as the input array.
    std::string errmsg = ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobianInv):";

    INTREPID2_TEST_FOR_EXCEPTION( !(requireRankMatch(errmsg, jacobianInv, jacobian) ), std::invalid_argument, errmsg);
  
    INTREPID2_TEST_FOR_EXCEPTION( !(requireDimensionMatch(errmsg, jacobianInv, jacobian) ), std::invalid_argument, errmsg);
  }



  template<class Scalar>
  template<class ArrayJacDet, class ArrayJac>
  void CellTools<Scalar>::validateArguments_setJacobianDetArgs(const ArrayJacDet &  jacobianDet,
                                                               const ArrayJac    &  jacobian)
  {
    // Validate input jacobian array: admissible ranks & dimensions are: 
    // - rank-4 with dimensions (C,P,D,D), or rank-3 with dimensions (P,D,D).
    if(!CheckType<ArrayJac>::value || !CheckType<ArrayJacDet>::value){
      std::cout <<std::endl<<"WARNING:: A Nonsupported Container is Being used with Intrepid2::CellTools<Scalar>::setJacobianDet"<<std::endl;
    }
    int jacobRank = getrank(jacobian);
    INTREPID2_TEST_FOR_EXCEPTION( !( (jacobRank == 4) || (jacobRank == 3) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobianInv): rank = 4 or 3 required for jacobian array. ");
  
    // Verify correctness of spatial dimensions - they are the last two dimensions of the array: rank-2 and rank-1
    INTREPID2_TEST_FOR_EXCEPTION( !(jacobian.dimension(jacobRank - 1) == jacobian.dimension(jacobRank - 2) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobianInv): dim(rank-2) = dim(rank-2) (same spatial dimensions) required for jacobian array. ");
  
    INTREPID2_TEST_FOR_EXCEPTION( !( (0 < jacobian.dimension(jacobRank - 1) ) && (jacobian.dimension(jacobRank - 1) < 4) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobianInv): dim(rank-1) and dim(rank-2) (spatial dimensions) must be between 1 and 3. ");

  
    // Validate output jacobianDet array: must be rank-2 with dimensions (C,P) if jacobian was rank-4:
    if(jacobRank == 4){
      INTREPID2_TEST_FOR_EXCEPTION( !(getrank(jacobianDet) == 2), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobianDetArgs): rank = 2 required for jacobianDet if jacobian is rank-4. ");
    
      INTREPID2_TEST_FOR_EXCEPTION( !(static_cast<index_type>(jacobianDet.dimension(0)) == static_cast<index_type>(jacobian.dimension(0)) ), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobianDetArgs): dim 0 (number of cells) of jacobianDet array must equal dim 0 of jacobian array. ");
    
      INTREPID2_TEST_FOR_EXCEPTION( !(static_cast<index_type>(jacobianDet.dimension(1)) == static_cast<index_type>(jacobian.dimension(1)) ), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobianDetArgs): dim 1 (number of points) of jacobianDet array must equal dim 1 of jacobian array.");  
    }
  
    // must be rank-1 with dimension (P) if jacobian was rank-3
    else {
      INTREPID2_TEST_FOR_EXCEPTION( !(getrank(jacobianDet) == 1), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobianDetArgs): rank = 1 required for jacobianDet if jacobian is rank-3. ");
    
      INTREPID2_TEST_FOR_EXCEPTION( !(static_cast<index_type>(jacobianDet.dimension(0)) == static_cast<index_type>(jacobian.dimension(0)) ), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_setJacobianDetArgs): dim 0 (number of points) of jacobianDet array must equal dim 0 of jacobian array.");  
    }
  }



  template<class Scalar>
  template<class ArrayPhysPoint, class ArrayRefPoint, class ArrayCell>
  void CellTools<Scalar>::validateArguments_mapToPhysicalFrame(const ArrayPhysPoint &        physPoints,
                                                               const ArrayRefPoint  &        refPoints,
                                                               const ArrayCell      &        worksetCell,
                                                               const shards::CellTopology &  cellTopo,
                                                               const int&                    whichCell)
  {
    if(!CheckType<ArrayPhysPoint>::value || !CheckType<ArrayRefPoint>::value || !CheckType<ArrayCell>::value){
      std::cout <<std::endl<<"WARNING:: A Nonsupported Container is Being used with Intrepid2::CellTools<Scalar>::mapToPhysicalFrame"<<std::endl;
    }
    std::string errmsg = ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame):";
  
    // Validate worksetCell array
    INTREPID2_TEST_FOR_EXCEPTION( (getrank(worksetCell) != 3), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): rank = 3 required for worksetCell array");
  
    INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(worksetCell.dimension(0)) <= 0), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): dim 0 (number of cells) >= 1 required for worksetCell array");
  
    INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(worksetCell.dimension(1)) != (index_type)cellTopo.getSubcellCount(0) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): dim 1 (number of cell nodes) of worksetCell array does not match cell topology");
  
    INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(worksetCell.dimension(2)) != (index_type)cellTopo.getDimension() ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): dim 2 (spatial dimension) of worksetCell array  does not match cell dimension");
  
    


    INTREPID2_TEST_FOR_EXCEPTION( !( ( (0 <= whichCell ) && ((index_type)whichCell < (index_type)worksetCell.dimension(0) ) ) || (whichCell == -1) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): whichCell = -1 or a valid cell ordinal is required.");
  
    // Validate refPoints array: can be rank-2 (P,D) or rank-3 (C,P,D) array
    // If rank-2: admissible output array is (P,D) or (C,P,D); admissible whichCell: -1 (default) or cell ordinal
    if(getrank(refPoints) == 2) {
      INTREPID2_TEST_FOR_EXCEPTION( (refPoints.dimension(0) <= 0), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): dim 0 (number of points) >= 1 required for refPoints array ");
    
      INTREPID2_TEST_FOR_EXCEPTION( ((index_type)refPoints.dimension(1) != (index_type)cellTopo.getDimension() ), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): dim 1 (spatial dimension) of refPoints array does not match cell dimension");

      // Validate output array: whichCell = -1 requires rank-3 array with dimensions (C,P,D)  
      if(whichCell == -1) {
        INTREPID2_TEST_FOR_EXCEPTION( ( (getrank(physPoints) != 3) && (whichCell == -1) ), std::invalid_argument,
                                      ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): rank = 3 required for physPoints array for the default whichCell value");
      
        INTREPID2_TEST_FOR_EXCEPTION( ((index_type)physPoints.dimension(0) != (index_type)worksetCell.dimension(0)), std::invalid_argument,
                                      ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): dim 0 (number of cells) of physPoints array must equal dim 0 of worksetCell array");
      
        INTREPID2_TEST_FOR_EXCEPTION( ((index_type)physPoints.dimension(1) != (index_type)refPoints.dimension(0)), std::invalid_argument,
                                      ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): dim 1 (number of points) of physPoints array must equal dim 0 of refPoints array"); 
      
        INTREPID2_TEST_FOR_EXCEPTION( ((index_type)physPoints.dimension(2) != (index_type)cellTopo.getDimension()), std::invalid_argument,
                                      ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): dim 2 (spatial dimension) does not match cell dimension ");  
      }
      // 0 <= whichCell < num cells requires rank-2 (P,D) arrays for both refPoints and physPoints
      else{
        INTREPID2_TEST_FOR_EXCEPTION( (getrank(physPoints) != 2), std::invalid_argument,
                                      ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): rank = 2 required for physPoints array");
      
        INTREPID2_TEST_FOR_EXCEPTION( ((index_type)physPoints.dimension(0) != (index_type)refPoints.dimension(0)), std::invalid_argument,
                                      ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): dim 0 (number of points) of physPoints array must equal dim 0 of refPoints array"); 
      
        INTREPID2_TEST_FOR_EXCEPTION( ((index_type)physPoints.dimension(1) != (index_type)cellTopo.getDimension()), std::invalid_argument,
                                      ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): dim 1 (spatial dimension) does not match cell dimension ");      
      }
    }
    // refPoints is (C,P,D): requires physPoints to be (C,P,D) and whichCell=-1  (because all cell mappings are applied)
    else if(getrank(refPoints) == 3) {
    
      // 1. validate refPoints dimensions and rank
      INTREPID2_TEST_FOR_EXCEPTION( ((index_type)refPoints.dimension(0) !=(index_type) worksetCell.dimension(0) ), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): dim 0 (number of cells) of refPoints and worksetCell arraya are required to match ");

      INTREPID2_TEST_FOR_EXCEPTION( (refPoints.dimension(1) <= 0), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): dim 1 (number of points) >= 1 required for refPoints array ");
    
      INTREPID2_TEST_FOR_EXCEPTION( ((index_type)refPoints.dimension(2) != (index_type)cellTopo.getDimension() ), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): dim 2 (spatial dimension) of refPoints array does not match cell dimension");
    
      // 2. whichCell  must be -1
      INTREPID2_TEST_FOR_EXCEPTION( (whichCell != -1), std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): default value is required for rank-3 refPoints array");

      // 3.  physPoints must match rank and dimensions of refPoints
      INTREPID2_TEST_FOR_EXCEPTION( !requireRankMatch(errmsg, refPoints, physPoints), std::invalid_argument, errmsg );
      INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, refPoints, physPoints), std::invalid_argument, errmsg);
    }
    // if rank is not 2 or 3 throw exception
    else {
      INTREPID2_TEST_FOR_EXCEPTION( !( (getrank(refPoints) == 2) || (getrank(refPoints) == 3) ), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToPhysicalFrame): rank = 2 or 3 required for refPoints array");
    }
  }
  template<class Scalar>
  template<class ArrayRefPoint, class ArrayPhysPoint, class ArrayCell>
  void CellTools<Scalar>::validateArguments_mapToReferenceFrame(const ArrayRefPoint  &        refPoints,
                                                                const ArrayPhysPoint &        physPoints,
                                                                const ArrayCell      &        worksetCell,
                                                                const shards::CellTopology &  cellTopo,
                                                                const int&                    whichCell)
  {
    std::string errmsg  = ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToReferenceFrame):";
    std::string errmsg1 = ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToReferenceFrame):";
  
    // Validate worksetCell array
    INTREPID2_TEST_FOR_EXCEPTION( (getrank(worksetCell) != 3), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToReferenceFrame): rank = 3 required for worksetCell array");
  
    INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(worksetCell.dimension(0)) <= 0), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToReferenceFrame): dim 0 (number of cells) >= 1 required for worksetCell array");
  
    INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(worksetCell.dimension(1)) != (index_type)cellTopo.getSubcellCount(0) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToReferenceFrame): dim 1 (number of cell nodes) of worksetCell array does not match cell topology");
  
    INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(worksetCell.dimension(2)) != (index_type)cellTopo.getDimension() ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToReferenceFrame): dim 2 (spatial dimension) of worksetCell array  does not match cell dimension");
    
    // Validate whichCell. It can be either -1 (default value) or a valid celli ordinal.
    INTREPID2_TEST_FOR_EXCEPTION( !( ( (0 <= whichCell ) && ((index_type)whichCell <(index_type) worksetCell.dimension(0) ) ) || (whichCell == -1) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToReferenceFrame): whichCell = -1 or a valid cell ordinal is required.");  
    // Admissible ranks and dimensions of refPoints and physPoints depend on whichCell value:
    // default is to map multiple sets of points to multiple sets of points. (C,P,D) arrays required
    int validRank;
    if(whichCell == -1) {
      validRank = 3;
      errmsg1 += " default value of whichCell requires rank-3 arrays:";
    }
    // whichCell is valid cell ordinal => we map single set of pts to a single set of pts. (P,D) arrays required
    else{
      errmsg1 += " rank-2 arrays required when whichCell is valid cell ordinal";
      validRank = 2;
    }
    INTREPID2_TEST_FOR_EXCEPTION( !requireRankRange(errmsg1, refPoints,  validRank,validRank), std::invalid_argument, errmsg1);
    INTREPID2_TEST_FOR_EXCEPTION( !requireRankMatch(errmsg1, physPoints, refPoints),           std::invalid_argument, errmsg1);
    INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg1, refPoints, physPoints),      std::invalid_argument, errmsg1);
  }



  template<class Scalar>
  template<class ArrayRefPoint, class ArrayInitGuess, class ArrayPhysPoint, class ArrayCell>
  void CellTools<Scalar>::validateArguments_mapToReferenceFrame(const ArrayRefPoint  &        refPoints,
                                                                const ArrayInitGuess &        initGuess,
                                                                const ArrayPhysPoint &        physPoints,
                                                                const ArrayCell      &        worksetCell,
                                                                const shards::CellTopology &  cellTopo,
                                                                const int&                    whichCell)
  {
    // Call the method that validates arguments with the default initial guess selection
    validateArguments_mapToReferenceFrame(refPoints, physPoints, worksetCell, cellTopo, whichCell);
  
    // Then check initGuess: its rank and dimensions must match those of physPoints.
    std::string errmsg = ">>> ERROR (Intrepid2::CellTools::validateArguments_mapToReferenceFrame):";
    INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, initGuess, physPoints), std::invalid_argument, errmsg);  
  }


  template<class Scalar>
  template<class ArrayIncl, class ArrayPoint, class ArrayCell>
  void CellTools<Scalar>::validateArguments_checkPointwiseInclusion(ArrayIncl &                   inCell,
                                                                    const ArrayPoint &            physPoints,
                                                                    const ArrayCell &             worksetCell,
                                                                    const int &                   whichCell,
                                                                    const shards::CellTopology &  cell)
  {
    // Validate worksetCell array
    INTREPID2_TEST_FOR_EXCEPTION( (getrank(worksetCell) != 3), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): rank = 3 required for worksetCell array");
  
    INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(worksetCell.dimension(0)) <= 0), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): dim 0 (number of cells) >= 1 required for worksetCell array");
  
    INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(worksetCell.dimension(1)) != (index_type)cell.getSubcellCount(0) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): dim 1 (number of cell nodes) of worksetCell array does not match cell topology");
  
    INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(worksetCell.dimension(2)) != (index_type)cell.getDimension() ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): dim 2 (spatial dimension) of worksetCell array  does not match cell dimension");
  
  
    // Validate whichCell It can be either -1 (default value) or a valid cell ordinal.
    INTREPID2_TEST_FOR_EXCEPTION( !( ( (0 <= whichCell ) && (whichCell < worksetCell.dimension(0) ) ) || (whichCell == -1) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): whichCell = -1 or a valid cell ordinal is required.");  
  
    // Validate points array: can be rank-2 (P,D) or rank-3 (C,P,D)
    // If rank-2: admissible inCell is rank-1 (P); admissible whichCell is valid cell ordinal but not -1.
    if(getrank(physPoints) == 2) {
    
      INTREPID2_TEST_FOR_EXCEPTION( (whichCell == -1), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): whichCell = a valid cell ordinal is required with rank-2 input array.");

      INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(physPoints.dimension(0)) <= 0), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): dim 0 (number of points) >= 1 required for physPoints array ");
    
      INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(physPoints.dimension(1)) != (index_type)cell.getDimension() ), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): dim 1 (spatial dimension) of physPoints array does not match cell dimension");
    
      // Validate inCell
      INTREPID2_TEST_FOR_EXCEPTION( (getrank(inCell) != 1), std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): rank = 1 required for inCell array");
    
      INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(inCell.dimension(0)) != static_cast<index_type>(physPoints.dimension(0))), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): dim 0 (number of points) of inCell array must equal dim 0 of physPoints array");
    }
    // If rank-3: admissible inCell is rank-2 (C,P); admissible whichCell = -1.
    else if (getrank(physPoints) == 3){
    
      INTREPID2_TEST_FOR_EXCEPTION( !(whichCell == -1), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): whichCell = -1 is required with rank-3 input array.");
    
      INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(physPoints.dimension(0)) != static_cast<index_type>(worksetCell.dimension(0)) ), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): dim 0 (number of cells)  of physPoints array must equal dim 0 of worksetCell array ");

      INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(physPoints.dimension(1)) <= 0), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): dim 1 (number of points) >= 1 required for physPoints array ");
    
      INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(physPoints.dimension(2)) != (index_type)cell.getDimension() ), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): dim 2 (spatial dimension) of physPoints array does not match cell dimension");
    
      // Validate inCell
      INTREPID2_TEST_FOR_EXCEPTION( (getrank(inCell) != 2), std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): rank = 2 required for inCell array");
    
      INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(inCell.dimension(0)) != static_cast<index_type>(physPoints.dimension(0))), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): dim 0 (number of cells) of inCell array must equal dim 0 of physPoints array");    

      INTREPID2_TEST_FOR_EXCEPTION( (static_cast<index_type>(inCell.dimension(1)) != static_cast<index_type>(physPoints.dimension(1))), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): dim 1 (number of points) of inCell array must equal dim 1 of physPoints array");    
    }
    else {
      INTREPID2_TEST_FOR_EXCEPTION( !( (getrank(physPoints) == 2) && (getrank(physPoints) ==3) ), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::validateArguments_checkPointwiseInclusion): rank = 2 or 3 required for points array");
    }
  }
}

#endif
