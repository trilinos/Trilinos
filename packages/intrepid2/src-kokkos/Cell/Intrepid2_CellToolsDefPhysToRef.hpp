1;2c// @HEADER
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
#ifndef __INTREPID2_CELLTOOLS_DEF_PHYS_TO_REF_HPP__
#define __INTREPID2_CELLTOOLS_DEF_PHYS_TO_REF_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {

  
  //============================================================================================//          
  //                                                                                            //          
  //                      Reference-to-physical frame mapping and its inverse                   //          
  //                                                                                            //          
  //============================================================================================//   
  

  template<class Scalar>
  template<class ArrayRefPoint, class ArrayPhysPoint, class ArrayCell>
  void CellTools<Scalar>::mapToReferenceFrame(ArrayRefPoint        &        refPoints,
                                              const ArrayPhysPoint &        physPoints,
                                              const ArrayCell      &        cellWorkset,
                                              const shards::CellTopology &  cellTopo,
                                              const int &                   whichCell)
  {
    INTREPID2_VALIDATE( validateArguments_mapToReferenceFrame(refPoints, physPoints, cellWorkset, cellTopo, whichCell) );
  
    index_type spaceDim  = (index_type)cellTopo.getDimension();
    index_type numPoints;
    index_type numCells;

    // Define initial guesses to be  the Cell centers of the reference cell topology
    FieldContainer<Scalar> cellCenter(spaceDim);
    switch( cellTopo.getKey() ){
      // Standard Base topologies (number of cellWorkset = number of vertices)
    case shards::Line<2>::key:
      cellCenter(0) = 0.0;    break;

    case shards::Triangle<3>::key:
    case shards::Triangle<6>::key:    
      cellCenter(0) = 1./3.;    cellCenter(1) = 1./3.;  break;
      
    case shards::Quadrilateral<4>::key:
    case shards::Quadrilateral<9>::key:
      cellCenter(0) = 0.0;      cellCenter(1) = 0.0;    break;
      
    case shards::Tetrahedron<4>::key:
    case shards::Tetrahedron<10>::key:
    case shards::Tetrahedron<11>::key:
      cellCenter(0) = 1./6.;    cellCenter(1) =  1./6.;    cellCenter(2) =  1./6.;  break;
      
    case shards::Hexahedron<8>::key:
    case shards::Hexahedron<20>::key:
    case shards::Hexahedron<27>::key:
      cellCenter(0) = 0.0;      cellCenter(1) =  0.0;       cellCenter(2) =  0.0;   break;

    case shards::Wedge<6>::key:
    case shards::Wedge<15>::key:
    case shards::Wedge<18>::key:
      cellCenter(0) = 1./3.;    cellCenter(1) =  1./3.;     cellCenter(2) = 0.0;    break;

    case shards::Pyramid<5>::key:
    case shards::Pyramid<13>::key:
      cellCenter(0) = 0.;       cellCenter(1) = 0.;         cellCenter(2) = 0.25;    break;

      // These extended topologies are not used for mapping purposes
    case shards::Quadrilateral<8>::key:
      INTREPID2_TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::mapToReferenceFrame): Cell topology not supported. ");
      break;

      // Base and Extended Line, Beam and Shell topologies  
    case shards::Line<3>::key:
    case shards::Beam<2>::key:
    case shards::Beam<3>::key:
    case shards::ShellLine<2>::key:
    case shards::ShellLine<3>::key:
    case shards::ShellTriangle<3>::key:
    case shards::ShellTriangle<6>::key:
    case shards::ShellQuadrilateral<4>::key:
    case shards::ShellQuadrilateral<8>::key:
    case shards::ShellQuadrilateral<9>::key:
      INTREPID2_TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::mapToReferenceFrame): Cell topology not supported. ");
      break;
    default:
      INTREPID2_TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::mapToReferenceFrame): Cell topology not supported.");        
    }// switch key 
  
    // Resize initial guess depending on the rank of the physical points array
    FieldContainer<Scalar> initGuess;
  
    // Default: map (C,P,D) array of physical pt. sets to (C,P,D) array. Requires (C,P,D) initial guess.
    if(whichCell == -1){
      numPoints = static_cast<index_type>(physPoints.dimension(1));
      numCells = static_cast<index_type>(cellWorkset.dimension(0));
      initGuess.resize(numCells, numPoints, spaceDim);
      // Set initial guess:
      for(index_type c = 0; c < numCells; c++){
        for(index_type p = 0; p < numPoints; p++){
          for(index_type d = 0; d < spaceDim; d++){
            initGuess(c, p, d) = cellCenter(d);
          }// d
        }// p
      }// c
    }
    // Custom: map (P,D) array of physical pts. to (P,D) array. Requires (P,D) initial guess.
    else {
      numPoints = static_cast<index_type>(physPoints.dimension(0));
      initGuess.resize(numPoints, spaceDim);
      // Set initial guess:
      for(index_type p = 0; p < numPoints; p++){
        for(index_type d = 0; d < spaceDim; d++){
          initGuess(p, d) = cellCenter(d);
        }// d
      }// p
    }
    // Call method with initial guess
    mapToReferenceFrameInitGuess(refPoints, initGuess, physPoints, cellWorkset, cellTopo, whichCell);  

  }
  
  
  template<class Scalar>
  template<class ArrayRefPoint, class ArrayInitGuess, class ArrayPhysPoint, class ArrayCell>
  void CellTools<Scalar>::mapToReferenceFrameInitGuess(ArrayRefPoint        &        refPoints,
                                                       const ArrayInitGuess &        initGuess,
                                                       const ArrayPhysPoint &        physPoints,
                                                       const ArrayCell      &        cellWorkset,
                                                       const Teuchos::RCP<Basis<Scalar, FieldContainer<Scalar> > > HGRAD_Basis,
                                                       const int &                   whichCell)
  {
    ArrayWrapper<Scalar,ArrayInitGuess, Rank<ArrayInitGuess >::value, true>initGuessWrap(initGuess);
    ArrayWrapper<Scalar,ArrayRefPoint, Rank<ArrayRefPoint >::value, false>refPointsWrap(refPoints);
    // INTREPID2_VALIDATE( validateArguments_mapToReferenceFrame(refPoints, initGuess, physPoints, cellWorkset, cellTopo, whichCell) );
    index_type spaceDim  = (index_type)HGRAD_Basis->cellTopology().getDimension();
    index_type numPoints;
    index_type numCells=0;
  
    // Temp arrays for Newton iterates and Jacobians. Resize according to rank of ref. point array
    FieldContainer<Scalar> xOld;
    FieldContainer<Scalar> xTem;  
    FieldContainer<Scalar> jacobian;
    FieldContainer<Scalar> jacobInv;
    FieldContainer<Scalar> error; 
    FieldContainer<Scalar> cellCenter(spaceDim);
  
    // Default: map (C,P,D) array of physical pt. sets to (C,P,D) array. Requires (C,P,D) temp arrays and (C,P,D,D) Jacobians.
    if(whichCell == -1){
      numPoints = static_cast<index_type>(physPoints.dimension(1));
      numCells = static_cast<index_type>(cellWorkset.dimension(0));
      xOld.resize(numCells, numPoints, spaceDim);
      xTem.resize(numCells, numPoints, spaceDim);  
      jacobian.resize(numCells,numPoints, spaceDim, spaceDim);
      jacobInv.resize(numCells,numPoints, spaceDim, spaceDim);
      error.resize(numCells,numPoints); 
      // Set initial guess to xOld
      for(index_type c = 0; c < numCells; c++){
        for(index_type p = 0; p < numPoints; p++){
          for(index_type d = 0; d < spaceDim; d++){
            xOld(c, p, d) = initGuessWrap(c, p, d);
          }// d
        }// p
      }// c
    }
    // Custom: map (P,D) array of physical pts. to (P,D) array. Requires (P,D) temp arrays and (P,D,D) Jacobians.
    else {
      numPoints = static_cast<index_type>(physPoints.dimension(0));
      xOld.resize(numPoints, spaceDim);
      xTem.resize(numPoints, spaceDim);  
      jacobian.resize(numPoints, spaceDim, spaceDim);
      jacobInv.resize(numPoints, spaceDim, spaceDim);
      error.resize(numPoints); 
      // Set initial guess to xOld
      for(index_type p = 0; p < numPoints; p++){
        for(index_type d = 0; d < spaceDim; d++){
          xOld(p, d) = initGuessWrap(p, d);
        }// d
      }// p
    }
  
    // Newton method to solve the equation F(refPoints) - physPoints = 0:
    // refPoints = xOld - DF^{-1}(xOld)*(F(xOld) - physPoints) = xOld + DF^{-1}(xOld)*(physPoints - F(xOld))
    for(int iter = 0; iter < INTREPID2_MAX_NEWTON; ++iter) {
    
      // Jacobians at the old iterates and their inverses. 
      setJacobian(jacobian, xOld, cellWorkset, HGRAD_Basis, whichCell);
      setJacobianInv(jacobInv, jacobian);
      // The Newton step.
      mapToPhysicalFrame( xTem, xOld, cellWorkset, HGRAD_Basis->cellTopology(), whichCell );      // xTem <- F(xOld)
      RealSpaceTools<Scalar>::subtract( xTem, physPoints, xTem );        // xTem <- physPoints - F(xOld)
      RealSpaceTools<Scalar>::matvec( refPoints, jacobInv, xTem);        // refPoints <- DF^{-1}( physPoints - F(xOld) )
      RealSpaceTools<Scalar>::add( refPoints, xOld );                    // refPoints <- DF^{-1}( physPoints - F(xOld) ) + xOld

      // l2 error (Euclidean distance) between old and new iterates: |xOld - xNew|
      RealSpaceTools<Scalar>::subtract( xTem, xOld, refPoints );
      RealSpaceTools<Scalar>::vectorNorm( error, xTem, NORM_TWO );

      // Average L2 error for a multiple sets of physical points: error is rank-2 (C,P) array 
      Scalar totalError;
      if(whichCell == -1) {
        FieldContainer<Scalar> cellWiseError(numCells);
        // error(C,P) -> cellWiseError(P)

        RealSpaceTools<Scalar>::vectorNorm( cellWiseError, error, NORM_ONE );
        totalError = RealSpaceTools<Scalar>::vectorNorm( cellWiseError, NORM_ONE );
      }
      //Average L2 error for a single set of physical points: error is rank-1 (P) array
      else{

        totalError = RealSpaceTools<Scalar>::vectorNorm( error, NORM_ONE ); 
        totalError = totalError;
      }
    
      // Stopping criterion:
      if (totalError < INTREPID_TOL) {
        break;
      } 
      else if ( iter > INTREPID2_MAX_NEWTON) {
        INTREPID2_VALIDATE(std::cout << " Intrepid2::CellTools::mapToReferenceFrameInitGuess failed to converge to desired tolerance within " 
                           << INTREPID2_MAX_NEWTON  << " iterations\n" );
        break;
      }

      // initialize next Newton step
      //    xOld = refPoints;
      int refPointsRank=getrank(refPoints);
      if (refPointsRank==3){
        for(index_type i=0;i<static_cast<index_type>(refPoints.dimension(0));i++){
          for(index_type j=0;j<static_cast<index_type>(refPoints.dimension(1));j++){
            for(index_type k=0;k<static_cast<index_type>(refPoints.dimension(2));k++){
              xOld(i,j,k) = refPointsWrap(i,j,k);
            }
          }
        }
      }else if(refPointsRank==2){
        for(index_type i=0;i<static_cast<index_type>(refPoints.dimension(0));i++){
          for(index_type j=0;j<static_cast<index_type>(refPoints.dimension(1));j++){
            xOld(i,j) = refPointsWrap(i,j);
          }
        }

      }



    } // for(iter)
  }


  template<class Scalar>
  template<class ArrayRefPoint, class ArrayInitGuess, class ArrayPhysPoint, class ArrayCell>
  void CellTools<Scalar>::mapToReferenceFrameInitGuess(ArrayRefPoint        &        refPoints,
                                                       const ArrayInitGuess &        initGuess,
                                                       const ArrayPhysPoint &        physPoints,
                                                       const ArrayCell      &        cellWorkset,
                                                       const shards::CellTopology &  cellTopo,
                                                       const int &                   whichCell)
  {
    ArrayWrapper<Scalar,ArrayInitGuess, Rank<ArrayInitGuess >::value, true>initGuessWrap(initGuess);
    ArrayWrapper<Scalar,ArrayRefPoint, Rank<ArrayRefPoint >::value, false>refPointsWrap(refPoints);
    INTREPID2_VALIDATE( validateArguments_mapToReferenceFrame(refPoints, initGuess, physPoints, cellWorkset, cellTopo, whichCell) );
    index_type spaceDim  = (index_type)cellTopo.getDimension();
    index_type numPoints;
    index_type numCells=0;
  
    // Temp arrays for Newton iterates and Jacobians. Resize according to rank of ref. point array
    FieldContainer<Scalar> xOld;
    FieldContainer<Scalar> xTem;  
    FieldContainer<Scalar> jacobian;
    FieldContainer<Scalar> jacobInv;
    FieldContainer<Scalar> error; 
    FieldContainer<Scalar> cellCenter(spaceDim);
  
    // Default: map (C,P,D) array of physical pt. sets to (C,P,D) array. Requires (C,P,D) temp arrays and (C,P,D,D) Jacobians.
    if(whichCell == -1){
      numPoints = static_cast<index_type>(physPoints.dimension(1));
      numCells = static_cast<index_type>(cellWorkset.dimension(0));
      xOld.resize(numCells, numPoints, spaceDim);
      xTem.resize(numCells, numPoints, spaceDim);  
      jacobian.resize(numCells,numPoints, spaceDim, spaceDim);
      jacobInv.resize(numCells,numPoints, spaceDim, spaceDim);
      error.resize(numCells,numPoints); 
      // Set initial guess to xOld
      for(index_type c = 0; c < numCells; c++){
        for(index_type p = 0; p < numPoints; p++){
          for(index_type d = 0; d < spaceDim; d++){
            xOld(c, p, d) = initGuessWrap(c, p, d);
          }// d
        }// p
      }// c
    }
    // Custom: map (P,D) array of physical pts. to (P,D) array. Requires (P,D) temp arrays and (P,D,D) Jacobians.
    else {
      numPoints = static_cast<index_type>(physPoints.dimension(0));
      xOld.resize(numPoints, spaceDim);
      xTem.resize(numPoints, spaceDim);  
      jacobian.resize(numPoints, spaceDim, spaceDim);
      jacobInv.resize(numPoints, spaceDim, spaceDim);
      error.resize(numPoints); 
      // Set initial guess to xOld
      for(index_type p = 0; p < numPoints; p++){
        for(index_type d = 0; d < spaceDim; d++){
          xOld(p, d) = initGuessWrap(p, d);
        }// d
      }// p
    }
  
    // Newton method to solve the equation F(refPoints) - physPoints = 0:
    // refPoints = xOld - DF^{-1}(xOld)*(F(xOld) - physPoints) = xOld + DF^{-1}(xOld)*(physPoints - F(xOld))
    for(int iter = 0; iter < INTREPID2_MAX_NEWTON; ++iter) {
    
      // Jacobians at the old iterates and their inverses. 
      setJacobian(jacobian, xOld, cellWorkset, cellTopo, whichCell);
      setJacobianInv(jacobInv, jacobian);
      // The Newton step.
      mapToPhysicalFrame( xTem, xOld, cellWorkset, cellTopo, whichCell );      // xTem <- F(xOld)
      RealSpaceTools<Scalar>::subtract( xTem, physPoints, xTem );        // xTem <- physPoints - F(xOld)
      RealSpaceTools<Scalar>::matvec( refPoints, jacobInv, xTem);        // refPoints <- DF^{-1}( physPoints - F(xOld) )
      RealSpaceTools<Scalar>::add( refPoints, xOld );                    // refPoints <- DF^{-1}( physPoints - F(xOld) ) + xOld

      // l2 error (Euclidean distance) between old and new iterates: |xOld - xNew|
      RealSpaceTools<Scalar>::subtract( xTem, xOld, refPoints );
      RealSpaceTools<Scalar>::vectorNorm( error, xTem, NORM_TWO );

      // Average L2 error for a multiple sets of physical points: error is rank-2 (C,P) array 
      Scalar totalError;
      if(whichCell == -1) {
        FieldContainer<Scalar> cellWiseError(numCells);
        // error(C,P) -> cellWiseError(P)

        RealSpaceTools<Scalar>::vectorNorm( cellWiseError, error, NORM_ONE );
        totalError = RealSpaceTools<Scalar>::vectorNorm( cellWiseError, NORM_ONE );
      }
      //Average L2 error for a single set of physical points: error is rank-1 (P) array
      else{

        totalError = RealSpaceTools<Scalar>::vectorNorm( error, NORM_ONE ); 
        totalError = totalError;
      }
    
      // Stopping criterion:
      if (totalError < INTREPID_TOL) {
        break;
      } 
      else if ( iter > INTREPID2_MAX_NEWTON) {
        INTREPID2_VALIDATE(std::cout << " Intrepid2::CellTools::mapToReferenceFrameInitGuess failed to converge to desired tolerance within " 
                           << INTREPID2_MAX_NEWTON  << " iterations\n" );
        break;
      }

      // initialize next Newton step
      //    xOld = refPoints;
      int refPointsRank=getrank(refPoints);
      if (refPointsRank==3){
        for(index_type i=0;i<static_cast<index_type>(refPoints.dimension(0));i++){
          for(index_type j=0;j<static_cast<index_type>(refPoints.dimension(1));j++){
            for(index_type k=0;k<static_cast<index_type>(refPoints.dimension(2));k++){
              xOld(i,j,k) = refPointsWrap(i,j,k);
            }
          }
        }
      }else if(refPointsRank==2){
        for(index_type i=0;i<static_cast<index_type>(refPoints.dimension(0));i++){
          for(index_type j=0;j<static_cast<index_type>(refPoints.dimension(1));j++){
            xOld(i,j) = refPointsWrap(i,j);
          }
        }

      }



    } // for(iter)
  }


}

#endif
