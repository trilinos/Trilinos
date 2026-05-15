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

  namespace FunctorCellTools {

    /**
     \brief Functor for computing the scaled normal component of the residual, 3d version
     
    */

    template<typename OutViewType, typename TanViewType, typename InViewType>
    struct F_scaledResidualNormalComponent3d {
            OutViewType result_;
            TanViewType workview_;
            const TanViewType tangents_;
            const InViewType residual_;
            const typename OutViewType::value_type scaling_;

      KOKKOS_INLINE_FUNCTION
      F_scaledResidualNormalComponent3d( OutViewType result,
                           TanViewType workview,
                     const TanViewType tangents,
                     const InViewType residual,
                     const typename OutViewType::value_type scaling)
        : result_(result), workview_(workview), tangents_(tangents), residual_(residual), scaling_(scaling) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cell,
                      const ordinal_type point) const {
        
        auto t1 = Kokkos::subview( tangents_,  cell, point, Kokkos::ALL(), 0);
        auto t2 = Kokkos::subview( tangents_,  cell, point, Kokkos::ALL(), 1);
        auto n = Kokkos::subview( workview_,  cell, point, Kokkos::ALL());
        Kernels::Serial::vector_product_d3(n, t1, t2);
        result_(cell,point) = (n(0)*residual_(cell,point,0) + n(1)*residual_(cell,point,1) + n(2)*residual_(cell,point,2))*scaling_/Kernels::Serial::norm(n, NORM_TWO);
      }
    };

    /**
     \brief Functor for computing the scaled normal component of the residual, 2d version
    */
    template<typename OutViewType, typename TanViewType, typename InViewType>
    struct F_scaledResidualNormalComponent2d {
            OutViewType result_;
            const TanViewType tangents_;
            const InViewType residual_;
            const typename InViewType::value_type scaling_;

      KOKKOS_INLINE_FUNCTION
      F_scaledResidualNormalComponent2d( OutViewType result,
                     const TanViewType tangents,
                     const InViewType residual,
                     const typename InViewType::value_type scaling)
        : result_(result), tangents_(tangents), residual_(residual), scaling_(scaling) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cell,
                      const ordinal_type point) const {        
        auto t = Kokkos::subview( tangents_,  cell, point, Kokkos::ALL(), 0); // n = [-t_1, t_0]
        result_(cell,point) = (-t(1)*residual_(cell,point,0) + t(0)*residual_(cell,point,1))*scaling_/Kernels::Serial::norm(t, NORM_TWO);
      }
    };

    /**
     \brief Functor for computing the max norm of a view
    */
    template<typename ViewType>
    struct F_maxNorm {
      const ViewType view_;
      KOKKOS_INLINE_FUNCTION
      F_maxNorm( const ViewType view) 
        : view_(view) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cl,
                      const ordinal_type pt,
                      typename ViewType::value_type& lmax) const {
        if (view_(cl,pt) > lmax) lmax = view_(cl,pt);
      }
    };

  }

  template<typename DeviceType>
  template<typename OutViewType, typename TanViewType, typename InViewType>
  void
  CellTools<DeviceType>::
  scaledResidualNormalComponent(
               OutViewType result,
               const TanViewType tangents,
               const InViewType residual,
               const typename InViewType::value_type scaling)
  {
    constexpr bool is_accessible =
        Kokkos::SpaceAccessibility<MemSpaceType,
        typename decltype(tangents)::memory_space>::accessible;
    static_assert(is_accessible, "CellTools<DeviceType>::setJacobian(..): output view's memory space is not compatible with DeviceType");

    using range_policy_type = Kokkos::MDRangePolicy
      < ExecSpaceType, Kokkos::Rank<2>, Kokkos::IndexType<ordinal_type> >;
    range_policy_type policy( { 0, 0 },
                              { tangents.extent(0), tangents.extent(1) } );
    if (tangents.extent_int(2) == 2) {
      using FunctorType      = FunctorCellTools::F_scaledResidualNormalComponent2d<OutViewType, TanViewType, InViewType> ;
      Kokkos::parallel_for( policy, FunctorType(result, tangents, residual, scaling) );
    } else {
      using FunctorType      = FunctorCellTools::F_scaledResidualNormalComponent3d<OutViewType, TanViewType, InViewType> ;
      TanViewType work = Impl::createMatchingDynRankView(tangents, "work_view", tangents.extent_int(0), tangents.extent_int(1), tangents.extent_int(2));
      Kokkos::parallel_for( policy, FunctorType(result, work, tangents, residual, scaling) );
    }
  }



  template<typename DeviceType>
  template<typename refPointValueType,    class ...refPointProperties,
           typename physPointValueType,   class ...physPointProperties,
           typename worksetCellValueType, class ...worksetCellProperties>
  bool
  CellTools<DeviceType>::
  mapToReferenceFrame(       Kokkos::DynRankView<refPointValueType,refPointProperties...>       refPoints,
                       const Kokkos::DynRankView<physPointValueType,physPointProperties...>     physPoints,
                       const Kokkos::DynRankView<worksetCellValueType,worksetCellProperties...> worksetCell,
                       const shards::CellTopology cellTopo, 
                       const physPointValueType shellThickness ) {
    constexpr bool are_accessible =
        Kokkos::SpaceAccessibility<MemSpaceType,
        typename decltype(refPoints)::memory_space>::accessible &&
        Kokkos::SpaceAccessibility<MemSpaceType,
        typename decltype(physPoints)::memory_space>::accessible &&
        Kokkos::SpaceAccessibility<MemSpaceType,
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
    auto initGuess = Impl::createMatchingDynRankView(refPoints, "CellTools::mapToReferenceFrame::initGuess", numCells, numPoints, spaceDim );
    rst::clone(refPoints, cellCenter);
    
    auto basis = createHGradBasis<refPointValueType,refPointValueType>(cellTopo);
    return mapToReferenceFrame(refPoints, physPoints, worksetCell, basis, shellThickness);  
  }
  

  template<typename DeviceType>
  template<typename refPointValueType,    class ...refPointProperties,
           typename physPointValueType,   class ...physPointProperties,
           typename worksetCellValueType, class ...worksetCellProperties>
  bool
  CellTools<DeviceType>::
  mapToReferenceFrame(       Kokkos::DynRankView<refPointValueType,refPointProperties...>       refPoints,
                       const Kokkos::DynRankView<physPointValueType,physPointProperties...>     physPoints,
                       const Kokkos::DynRankView<worksetCellValueType,worksetCellProperties...> worksetCell,
                       const BasisPtr<DeviceType, refPointValueType, refPointValueType>         basis,
                       const physPointValueType shellThickness) {
    constexpr bool are_accessible =
        Kokkos::SpaceAccessibility<MemSpaceType,
        typename decltype(refPoints)::memory_space>::accessible &&
        Kokkos::SpaceAccessibility<MemSpaceType,
        typename decltype(physPoints)::memory_space>::accessible &&
        Kokkos::SpaceAccessibility<MemSpaceType,
        typename decltype(worksetCell)::memory_space>::accessible;

    static_assert(are_accessible, "CellTools<DeviceType>::mapToReferenceFrameInitGuess(..): input/output views' memory spaces are not compatible with DeviceType");

    const bool isShell = (shellThickness>0);

    //note, for Beam and Shell elements, refPoints have refDim+1 dimension 
    const ordinal_type refDim = basis->getDomainDimension();
    const ordinal_type physDim = physPoints.extent_int(2);

    // Default: map (C,P,D) array of physical pt. sets to (C,P,D) array. 
    // Requires (C,P,D) temp arrays and (C,P,D,D) Jacobians.
    const auto numCells = worksetCell.extent(0);
    const auto numPoints = physPoints.extent(1);

    using rst = RealSpaceTools<DeviceType>;
    auto refDimRange = Kokkos::pair<ordinal_type,ordinal_type>(0,refDim);

    using result_layout = typename DeduceLayout< decltype(refPoints) >::result_layout;
    // Temp arrays for Newton iterates and Jacobians. Resize according to rank of ref. point array
    auto xOld = Impl::createMatchingDynRankView(refPoints, "CellTools::mapToReferenceFrameInitGuess::xOld", numCells, numPoints, refDim);
    auto xTmp = Impl::createMatchingDynRankView(refPoints, "CellTools::mapToReferenceFrameInitGuess::xTmp", numCells, numPoints, refDim);
    auto physTmp = Impl::createMatchingDynRankView(refPoints, "CellTools::mapToReferenceFrameInitGuess::physTmp", numCells, numPoints, physDim);

    // deep copy may not work with FAD but this is right thing to do as it can move data between devices
    Kokkos::deep_copy(xOld, Kokkos::subview(refPoints, Kokkos::ALL(), Kokkos::ALL(), refDimRange));

    // jacobian should select fad dimension between xOld and worksetCell as they are input; no front interface yet
    using valueTypeJ = std::common_type_t<typename decltype(refPoints)::value_type, typename decltype(worksetCell)::value_type>;    
    using viewTypeJ = Kokkos::DynRankView<valueTypeJ, result_layout, DeviceType >;
    using view_factory = Impl::CreateViewFactory<decltype(refPoints), decltype(worksetCell)>;

    viewTypeJ metric;
    if(refDim < physDim)
      metric = view_factory::template create_view<viewTypeJ>(refPoints, worksetCell, "CellTools::mapToReferenceFrameInitGuess::metric", numCells, numPoints, refDim, refDim);
    viewTypeJ jacobian = view_factory::template create_view<viewTypeJ>(refPoints, worksetCell, "CellTools::mapToReferenceFrameInitGuess::jacobian", numCells, numPoints, physDim, refDim);
    
    //  when refDim < physDim, this will contain the inverse of the metric
    viewTypeJ jacobianInv = view_factory::template create_view<viewTypeJ>(refPoints, worksetCell, "CellTools::mapToReferenceFrameInitGuess::jacobianInv", numCells, numPoints, refDim, refDim);

    
    using errorViewType = Kokkos::DynRankView<typename ScalarTraits<refPointValueType>::scalar_type, DeviceType>;
    using errorType = typename errorViewType::value_type;
    errorViewType
      xScalarTmp    ("CellTools::mapToReferenceFrameInitGuess::xScalarTmp",     numCells, numPoints, refDim),
      errorPointwise("CellTools::mapToReferenceFrameInitGuess::errorPointwise", numCells, numPoints);

    
    const auto tol = tolerance<errorType>();

    errorType error(0), physInfNorm(0);

    rst::extractScalarValues(xScalarTmp, physPoints);
    rst::vectorNorm(errorPointwise, Kokkos::subview(xScalarTmp, Kokkos::ALL(), Kokkos::ALL(), refDimRange), NORM_TWO);
    using range_policy_type = Kokkos::MDRangePolicy<typename DeviceType::execution_space, Kokkos::Rank<2>, Kokkos::IndexType<ordinal_type> >;
    range_policy_type policy( { 0, 0 }, { numCells, numPoints } );
    using FunctorType = FunctorCellTools::F_maxNorm<errorViewType>;
    Kokkos::parallel_reduce("MaxReduction", policy, FunctorType(errorPointwise), Kokkos::Max<errorType>(physInfNorm));
    
    auto refPts = Kokkos::subview(refPoints,Kokkos::ALL(), Kokkos::ALL(), refDimRange);
      

    // Newton method to solve the equation F(refPoints) - physPoints = 0:
    // refPoints = xOld - DF^{-1}(xOld)*(F(xOld) - physPoints) = xOld + DF^{-1}(xOld)*(physPoints - F(xOld))
    bool converged(false);
    for (ordinal_type iter=0;iter<Parameters::MaxNewton;++iter) {
      
      // Jacobians at the old iterates and their inverses. 
      setJacobian(jacobian, xOld, worksetCell, basis);
      
      // The Newton step.
      mapToPhysicalFrame(physTmp, xOld, worksetCell, basis); // physTmp <- F(xOld)
      rst::subtract(physTmp, physPoints, physTmp);           // physTmp <- physPoints - F(xOld)

      // l2 error (Euclidean distance) between old and new iterates: |physPoints - F(xOld)|
      rst::extractScalarValues(xScalarTmp, physTmp);
      rst::vectorNorm(errorPointwise, Kokkos::subview(xScalarTmp, Kokkos::ALL(), Kokkos::ALL(), refDimRange), NORM_TWO);

      error = 0;
      Kokkos::parallel_reduce("MaxReduction", policy, FunctorType(errorPointwise), Kokkos::Max<errorType>(error));

      if (error < tol*physInfNorm) {
        converged = true;
        break;
      }


      if(refDim < physDim) {
        rst::matvec(xTmp, jacobian, physTmp, true);          // xTmp <- DF^{T}( physPoints - F(xOld) ))
        rst::AtA(metric,jacobian);
        setJacobianInv(jacobianInv, metric);             
        rst::matvec(refPts, jacobianInv, xTmp);              // refPoints <- (DF^{T} DF)^{-1} DF^{T}( physPoints - F(xOld) )
      }
      else { 
        setJacobianInv(jacobianInv, jacobian);
        rst::matvec(refPts, jacobianInv, physTmp);           // refPoints <- DF^{-1}( physPoints - F(xOld) )
      }
  
      // extract values
      rst::extractScalarValues(xScalarTmp, refPts);

      rst::add(refPts, xOld);                                 // refPoints <- refPoints + xOld

      
      error = 0;
      // l2 error (Euclidean distance) between old and new iterates: |xOld - xNew|
      rst::vectorNorm(errorPointwise, Kokkos::subview(xScalarTmp, Kokkos::ALL(), Kokkos::ALL(), refDimRange), NORM_TWO);
      Kokkos::parallel_reduce("MaxReduction", policy, FunctorType(errorPointwise), Kokkos::Max<errorType>(error));    
    
      // Stopping criterion:
      if (error < tol) {
        converged = true;
        break;
      }

      // initialize next Newton step ( this is not device friendly )
      Kokkos::deep_copy(xOld, refPts);
    }

    if(isShell) {
      //for shell elements, we can reconstruct the orthogonal reference component by approximating the 
      //signed distance between the physical point and the map to physical space of the reference line (in 2d) or triangle,quad (in 3d),
      //and dividing it by the shell half thickness. We do so by solving
      //(physPoints - F(xOld)) * n /|n|^2 / (H/2)
      //where H is the shell thickness and  n is the normal to the manifold defined by the map to physical frame (n is orthogonal to the tangents defined by the Jacobian)
      setJacobian(jacobian, refPts, worksetCell, basis);
      scaledResidualNormalComponent(
               Kokkos::subview(refPoints,Kokkos::ALL(), Kokkos::ALL(), refDim),
               jacobian,
               physTmp,
               2.0/shellThickness);
    } 

    return converged;
  }


  template<typename DeviceType>
  template<typename refPointValueType,    class ...refPointProperties,
           typename initGuessValueType,   class ...initGuessProperties,
           typename physPointValueType,   class ...physPointProperties,
           typename worksetCellValueType, class ...worksetCellProperties>
  bool
  CellTools<DeviceType>::
  mapToReferenceFrameInitGuess(       Kokkos::DynRankView<refPointValueType,refPointProperties...>       refPoints,
                                const Kokkos::DynRankView<initGuessValueType,initGuessProperties...>     initGuess,
                                const Kokkos::DynRankView<physPointValueType,physPointProperties...>     physPoints,
                                const Kokkos::DynRankView<worksetCellValueType,worksetCellProperties...> worksetCell,
                                const shards::CellTopology cellTopo,
                                const physPointValueType shellThickness) {
    
  #ifdef HAVE_INTREPID2_DEBUG
    CellTools_mapToReferenceFrameInitGuessArgs(refPoints, initGuess, physPoints, worksetCell, cellTopo);
  #endif
    
    Kokkos::deep_copy(refPoints,initGuess);
    auto basis = createHGradBasis<refPointValueType,refPointValueType>(cellTopo);
    return mapToReferenceFrame(refPoints,
                                 physPoints,
                                 worksetCell,
                                 basis,
                                 shellThickness);
  }

}

#endif
