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

  template<typename SpT>
  template<typename refPointValueType,    class ...refPointProperties,
           typename physPointValueType,   class ...physPointProperties,
           typename worksetCellValueType, class ...worksetCellProperties>
  void
  CellTools<SpT>::
  mapToReferenceFrame( /**/  Kokkos::DynRankView<refPointValueType,refPointProperties...>       refPoints,
                       const Kokkos::DynRankView<physPointValueType,physPointProperties...>     physPoints,
                       const Kokkos::DynRankView<worksetCellValueType,worksetCellProperties...> worksetCell,
                       const shards::CellTopology cellTopo ) {
#ifdef HAVE_INTREPID2_DEBUG
    CellTools_mapToReferenceFrameArgs(refPoints, physPoints, worksetCell, cellTopo);
#endif  
    typedef RealSpaceTools<SpT> rst;
    typedef Kokkos::DynRankView<refPointValueType,refPointProperties...> refPointViewType;
    typedef Kokkos::DynRankView<refPointValueType,SpT> refPointViewSpType;

    const auto spaceDim  = cellTopo.getDimension();
    refPointViewSpType 
      cellCenter("CellTools::mapToReferenceFrame::cellCenter", spaceDim), 
      cellVertex("CellTools::mapToReferenceFrame::cellCenter", spaceDim);
    getReferenceCellCenter(cellCenter, cellVertex, cellTopo);

    // Default: map (C,P,D) array of physical pt. sets to (C,P,D) array. Requires (C,P,D) initial guess.
    const auto numCells = worksetCell.dimension(0);
    const auto numPoints = physPoints.dimension(1);
    
    refPointViewSpType initGuess("CellTools::mapToReferenceFrame::initGuess", numCells, numPoints, spaceDim);
    rst::clone(initGuess, cellCenter);
    
    mapToReferenceFrameInitGuess(refPoints, initGuess, physPoints, worksetCell, cellTopo);  
  }
  

  template<typename SpT>
  template<typename refPointValueType,    class ...refPointProperties,
           typename initGuessValueType,   class ...initGuessProperties,
           typename physPointValueType,   class ...physPointProperties,
           typename worksetCellValueType, class ...worksetCellProperties,
           typename HGradBasisPtrType>
  void
  CellTools<SpT>::
  mapToReferenceFrameInitGuess( /**/  Kokkos::DynRankView<refPointValueType,refPointProperties...>       refPoints,
                                const Kokkos::DynRankView<initGuessValueType,initGuessProperties...>     initGuess,
                                const Kokkos::DynRankView<physPointValueType,physPointProperties...>     physPoints,
                                const Kokkos::DynRankView<worksetCellValueType,worksetCellProperties...> worksetCell,
                                const HGradBasisPtrType basis ) {
#ifdef HAVE_INTREPID2_DEBUG
    CellTools_mapToReferenceFrameInitGuessArgs(refPoints, initGuess, physPoints, worksetCell, 
                                               basis->getBaseCellTopology());
#endif
    const auto cellTopo = basis->getBaseCellTopology();
    const auto spaceDim = cellTopo.getDimension();

    // Default: map (C,P,D) array of physical pt. sets to (C,P,D) array. 
    // Requires (C,P,D) temp arrays and (C,P,D,D) Jacobians.
    const auto numCells = worksetCell.dimension(0);
    const auto numPoints = physPoints.dimension(1);

    typedef RealSpaceTools<SpT> rst;
    const auto tol = Parameters::Tolerence;

    // Temp arrays for Newton iterates and Jacobians. Resize according to rank of ref. point array
    typedef Kokkos::DynRankView<initGuessValueType,SpT> xViewType;
    xViewType 
      xOld("CellTools::mapToReferenceFrameInitGuess::xOld", numCells, numPoints, spaceDim),
      xTmp("CellTools::mapToReferenceFrameInitGuess::xTmp", numCells, numPoints, spaceDim);

    // Set initial guess to xOld; we can do soft copy with initGuess but leave it not touched 
    // user may want to keep the same init guess
    Kokkos::deep_copy(xOld, initGuess);
    
    typedef Kokkos::DynRankView<physPointValueType,SpT> jacobianViewType;
    jacobianViewType
      jacobian   ("CellTools::mapToReferenceFrameInitGuess::jacobian", numCells, numPoints, spaceDim, spaceDim),
      jacobianInv("CellTools::mapToReferenceFrameInitGuess::jacobian", numCells, numPoints, spaceDim, spaceDim);
    
    typedef Kokkos::DynRankView<initGuessValueType,SpT> errorViewType;
    errorViewType 
      errorPointwise("CellTools::mapToReferenceFrameInitGuess::errorPointwise", numCells, numPoints),
      errorCellwise ("CellTools::mapToReferenceFrameInitGuess::errorCellwise",  numCells);

    // Newton method to solve the equation F(refPoints) - physPoints = 0:
    // refPoints = xOld - DF^{-1}(xOld)*(F(xOld) - physPoints) = xOld + DF^{-1}(xOld)*(physPoints - F(xOld))
    for (auto iter=0;iter<Parameters::MaxNewton;++iter) {
    
      // Jacobians at the old iterates and their inverses. 
      setJacobian(jacobian, xOld, worksetCell, basis);
      setJacobianInv(jacobianInv, jacobian);

      // The Newton step.
      mapToPhysicalFrame(xTmp, xOld, worksetCell, basis); // xTmp <- F(xOld)
      rst::subtract(xTmp, physPoints, xTmp);              // xTmp <- physPoints - F(xOld)
      rst::matvec(refPoints, jacobianInv, xTmp);          // refPoints <- DF^{-1}( physPoints - F(xOld) )
      rst::add(refPoints, xOld);                          // refPoints <- DF^{-1}( physPoints - F(xOld) ) + xOld

      // l2 error (Euclidean distance) between old and new iterates: |xOld - xNew|
      rst::subtract(xTmp, xOld, refPoints);
      rst::vectorNorm(errorPointwise, xTmp, NORM_TWO);

      // Average L2 error for a multiple sets of physical points: error is rank-2 (C,P) array 
      rst::vectorNorm(errorCellwise, errorPointwise, NORM_ONE);
      const auto errorTotal = rst::Serial::vectorNorm(errorCellwise, NORM_ONE);
    
      // Stopping criterion:
      if (errorTotal < tol) 
        break;

      // initialize next Newton step
      //    xOld = refPoints;
      Kokkos::deep_copy(xOld, refPoints);
    }
  }

}

#endif
