// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file   Intrepid2_CellToolsDefPhysToRef.hpp
    \brief  Definition file for the physical to reference mappings in the Intrepid2::CellTools class.
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

  template<typename DeviceType>
  template<typename refPointValueType,    class ...refPointProperties,
           typename physPointValueType,   class ...physPointProperties,
           typename worksetCellValueType, class ...worksetCellProperties>
  void
  CellTools<DeviceType>::
  mapToReferenceFrame(       Kokkos::DynRankView<refPointValueType,refPointProperties...>       refPoints,
                       const Kokkos::DynRankView<physPointValueType,physPointProperties...>     physPoints,
                       const Kokkos::DynRankView<worksetCellValueType,worksetCellProperties...> worksetCell,
                       const shards::CellTopology cellTopo ) {
    constexpr bool are_accessible =
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(refPoints)::memory_space>::accessible &&
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(physPoints)::memory_space>::accessible &&
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(worksetCell)::memory_space>::accessible;

    static_assert(are_accessible, "CellTools<DeviceType>::mapToReferenceFrame(..): input/output views' memory spaces are not compatible with DeviceType");

#ifdef HAVE_INTREPID2_DEBUG
    CellTools_mapToReferenceFrameArgs(refPoints, physPoints, worksetCell, cellTopo);
#endif  
    using deviceType = typename decltype(refPoints)::device_type;

    typedef RealSpaceTools<deviceType> rst;
    typedef Kokkos::DynRankView<typename ScalarTraits<refPointValueType>::scalar_type,deviceType> refPointViewSpType;

    const auto spaceDim  = cellTopo.getDimension();
    refPointViewSpType 
      cellCenter("CellTools::mapToReferenceFrame::cellCenter", spaceDim);
    getReferenceCellCenter(cellCenter, cellTopo);

    // Default: map (C,P,D) array of physical pt. sets to (C,P,D) array. Requires (C,P,D) initial guess.
    const auto numCells = worksetCell.extent(0);
    const auto numPoints = physPoints.extent(1);
    
    // init guess is created locally and non fad whatever refpoints type is 
    using result_layout = typename DeduceLayout< decltype(refPoints) >::result_layout;
    auto vcprop = Kokkos::common_view_alloc_prop(refPoints);
    using common_value_type = typename decltype(vcprop)::value_type;
    Kokkos::DynRankView< common_value_type, result_layout, deviceType > initGuess ( Kokkos::view_alloc("CellTools::mapToReferenceFrame::initGuess", vcprop), numCells, numPoints, spaceDim );
    //refPointViewSpType initGuess("CellTools::mapToReferenceFrame::initGuess", numCells, numPoints, spaceDim);
    rst::clone(initGuess, cellCenter);
    
    mapToReferenceFrameInitGuess(refPoints, initGuess, physPoints, worksetCell, cellTopo);  
  }
  

  template<typename DeviceType>
  template<typename refPointValueType,    class ...refPointProperties,
           typename initGuessValueType,   class ...initGuessProperties,
           typename physPointValueType,   class ...physPointProperties,
           typename worksetCellValueType, class ...worksetCellProperties,
           typename HGradBasisPtrType>
  void
  CellTools<DeviceType>::
  mapToReferenceFrameInitGuess(       Kokkos::DynRankView<refPointValueType,refPointProperties...>       refPoints,
                                const Kokkos::DynRankView<initGuessValueType,initGuessProperties...>     initGuess,
                                const Kokkos::DynRankView<physPointValueType,physPointProperties...>     physPoints,
                                const Kokkos::DynRankView<worksetCellValueType,worksetCellProperties...> worksetCell,
                                const HGradBasisPtrType basis ) {
#ifdef HAVE_INTREPID2_DEBUG
    CellTools_mapToReferenceFrameInitGuessArgs(refPoints, initGuess, physPoints, worksetCell, 
                                               basis->getBaseCellTopology());

#endif

    constexpr bool are_accessible =
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(refPoints)::memory_space>::accessible &&
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(initGuess)::memory_space>::accessible &&
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(physPoints)::memory_space>::accessible &&
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(worksetCell)::memory_space>::accessible;

    static_assert(are_accessible, "CellTools<DeviceType>::mapToReferenceFrameInitGuess(..): input/output views' memory spaces are not compatible with DeviceType");


    const auto cellTopo = basis->getBaseCellTopology();
    const auto spaceDim = cellTopo.getDimension();

    // Default: map (C,P,D) array of physical pt. sets to (C,P,D) array. 
    // Requires (C,P,D) temp arrays and (C,P,D,D) Jacobians.
    const auto numCells = worksetCell.extent(0);
    const auto numPoints = physPoints.extent(1);

    using rst = RealSpaceTools<DeviceType>;
    const auto tol = tolerence();

    using result_layout = typename DeduceLayout< decltype(refPoints) >::result_layout;
    auto vcprop = Kokkos::common_view_alloc_prop(refPoints);
    using viewType = Kokkos::DynRankView<typename decltype(vcprop)::value_type, result_layout, DeviceType >;

    // Temp arrays for Newton iterates and Jacobians. Resize according to rank of ref. point array
    viewType xOld(Kokkos::view_alloc("CellTools::mapToReferenceFrameInitGuess::xOld", vcprop), numCells, numPoints, spaceDim);
    viewType xTmp(Kokkos::view_alloc("CellTools::mapToReferenceFrameInitGuess::xTmp", vcprop), numCells, numPoints, spaceDim);

    // deep copy may not work with FAD but this is right thing to do as it can move data between devices
    Kokkos::deep_copy(xOld, initGuess);

    // jacobian should select fad dimension between xOld and worksetCell as they are input; no front interface yet
    auto vcpropJ = Kokkos::common_view_alloc_prop(refPoints, worksetCell);
    using viewTypeJ = Kokkos::DynRankView<typename decltype(vcpropJ)::value_type, result_layout, DeviceType >;
    viewTypeJ jacobian(Kokkos::view_alloc("CellTools::mapToReferenceFrameInitGuess::jacobian", vcpropJ), numCells, numPoints, spaceDim, spaceDim);
    viewTypeJ jacobianInv(Kokkos::view_alloc("CellTools::mapToReferenceFrameInitGuess::jacobianInv", vcpropJ), numCells, numPoints, spaceDim, spaceDim);
    
    using errorViewType = Kokkos::DynRankView<typename ScalarTraits<refPointValueType>::scalar_type, DeviceType>;
    errorViewType
      xScalarTmp    ("CellTools::mapToReferenceFrameInitGuess::xScalarTmp",     numCells, numPoints, spaceDim),
      errorPointwise("CellTools::mapToReferenceFrameInitGuess::errorPointwise", numCells, numPoints),
      errorCellwise ("CellTools::mapToReferenceFrameInitGuess::errorCellwise",  numCells);

    // Newton method to solve the equation F(refPoints) - physPoints = 0:
    // refPoints = xOld - DF^{-1}(xOld)*(F(xOld) - physPoints) = xOld + DF^{-1}(xOld)*(physPoints - F(xOld))
    for (ordinal_type iter=0;iter<Parameters::MaxNewton;++iter) {
      
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

      // extract values
      rst::extractScalarValues(xScalarTmp, xTmp);
      rst::vectorNorm(errorPointwise, xScalarTmp, NORM_TWO);

      // Average L2 error for a multiple sets of physical points: error is rank-2 (C,P) array 
      rst::vectorNorm(errorCellwise, errorPointwise, NORM_ONE);

      auto errorCellwise_h = Kokkos::create_mirror_view(errorCellwise);
      Kokkos::deep_copy(errorCellwise_h, errorCellwise);
      const auto errorTotal = rst::Serial::vectorNorm(errorCellwise_h, NORM_ONE);
    
      // Stopping criterion:
      if (errorTotal < tol) 
        break;

      // initialize next Newton step ( this is not device friendly )
      Kokkos::deep_copy(xOld, refPoints);
    }
  }

}

#endif
