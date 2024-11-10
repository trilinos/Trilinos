// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file   Intrepid2_CellToolsDefRefToPhys.hpp
    \brief  Definition file for the reference to physical mappings in the Intrepid2::CellTools class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/
#ifndef __INTREPID2_CELLTOOLS_DEF_REF_TO_PHYS_HPP__
#define __INTREPID2_CELLTOOLS_DEF_REF_TO_PHYS_HPP__

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
     \brief Functor for mapping reference points to physical frame see Intrepid2::CellTools for more
    */
    template<typename physPointViewType,
             typename worksetCellType,
             typename basisValType>
    struct F_mapToPhysicalFrame {
            physPointViewType _physPoints;
      const worksetCellType   _worksetCells;
      const basisValType      _basisVals;

      KOKKOS_INLINE_FUNCTION
      F_mapToPhysicalFrame( physPointViewType physPoints_,
                            worksetCellType   worksetCells_,
                            basisValType      basisVals_ )
        : _physPoints(physPoints_), _worksetCells(worksetCells_), _basisVals(basisVals_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) const {
        size_type cell, pt;
        unrollIndex( cell, pt,
                           _physPoints.extent(0),
                           _physPoints.extent(1),
                           iter );
              auto phys = Kokkos::subdynrankview( _physPoints, cell, pt, Kokkos::ALL());

        const auto valRank = rank(_basisVals);
        const auto val = ( valRank == 2 ? Kokkos::subdynrankview( _basisVals,       Kokkos::ALL(), pt) :
                                          Kokkos::subdynrankview( _basisVals, cell, Kokkos::ALL(), pt));

        const ordinal_type dim = phys.extent(0);
        const ordinal_type cardinality = val.extent(0);

        for (ordinal_type i=0;i<dim;++i) {
          phys(i) = 0;
          for (ordinal_type bf=0;bf<cardinality;++bf)
            phys(i) += _worksetCells(cell, bf, i)*val(bf);
        }
      }
    };


    template<typename refSubcellViewType,
             typename paramPointsViewType,
             typename subcellMapViewType>
    struct F_mapReferenceSubcell {
      refSubcellViewType refSubcellPoints_;
      const paramPointsViewType paramPoints_;
      const subcellMapViewType subcellMap_;
      const ordinal_type subcellOrd_;
      const ordinal_type sideDim_;

      KOKKOS_INLINE_FUNCTION
      F_mapReferenceSubcell( refSubcellViewType  refSubcellPoints,
                             const paramPointsViewType paramPoints,
                             const subcellMapViewType  subcellMap,
                             ordinal_type subcellOrd)
        : refSubcellPoints_(refSubcellPoints), paramPoints_(paramPoints), subcellMap_(subcellMap),
          subcellOrd_(subcellOrd), sideDim_(paramPoints_.extent(1)){};

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type pt, const size_type d) const {

        refSubcellPoints_(pt, d) = subcellMap_(subcellOrd_, d, 0);
        for(ordinal_type k=0; k<sideDim_; ++k)
          refSubcellPoints_(pt, d) += subcellMap_(subcellOrd_, d, k+1)*paramPoints_(pt, k);
      }
    };

    template<typename refSubcellViewType,
          typename paramPointsViewType,
          typename subcellMapViewType,
          typename ordViewType>
    struct F_mapReferenceSubcellBatch {
      refSubcellViewType refSubcellPoints_;
      const paramPointsViewType paramPoints_;
      const subcellMapViewType subcellMap_;
      const ordViewType subcellOrd_;
      const ordinal_type sideDim_;

      KOKKOS_INLINE_FUNCTION
      F_mapReferenceSubcellBatch( refSubcellViewType  refSubcellPoints,
                             const paramPointsViewType paramPoints,
                             const subcellMapViewType  subcellMap,
                             ordViewType subcellOrd)
        : refSubcellPoints_(refSubcellPoints), paramPoints_(paramPoints), subcellMap_(subcellMap),
          subcellOrd_(subcellOrd), sideDim_(paramPoints_.extent(1)){};

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type isc, const size_type pt, const size_type d) const {

        refSubcellPoints_(isc, pt, d) = subcellMap_(subcellOrd_(isc), d, 0);
        for(ordinal_type k=0; k<sideDim_; ++k)
          refSubcellPoints_(isc, pt, d) += subcellMap_(subcellOrd_(isc), d, k+1)*paramPoints_(pt, k);
      }
    };
  }

  template<typename DeviceType>
  template<typename PhysPointViewType,
           typename RefPointViewType,
           typename WorksetType,
           typename HGradBasisPtrType>
  void
  CellTools<DeviceType>::
  mapToPhysicalFrame(      PhysPointViewType physPoints,
                     const RefPointViewType  refPoints,
                     const WorksetType       worksetCell,
                     const HGradBasisPtrType basis ) {
#ifdef HAVE_INTREPID2_DEBUG
    CellTools_mapToPhysicalFrameArgs( physPoints, refPoints, worksetCell, basis->getBaseCellTopology() );
#endif
    constexpr bool are_accessible =
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(physPoints)::memory_space>::accessible &&
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(refPoints)::memory_space>::accessible;

    static_assert(are_accessible, "CellTools<DeviceType>::mapToPhysicalFrame(..): input/output views' memory spaces are not compatible with DeviceType");

    const auto cellTopo = basis->getBaseCellTopology();
    const auto numCells = worksetCell.extent(0);

    //points can be rank-2 (P,D), or rank-3 (C,P,D)
    const auto refPointRank = refPoints.rank();
    const auto numPoints = (refPointRank == 2 ? refPoints.extent(0) : refPoints.extent(1));
    const auto basisCardinality = basis->getCardinality();
    auto vcprop = Kokkos::common_view_alloc_prop(physPoints);

    using valViewType = Kokkos::DynRankView<decltype(basis->getDummyOutputValue()),DeviceType>;

    valViewType vals;

    switch (refPointRank) {
    case 2: {
      // refPoints is (P,D): single set of ref. points is mapped to one or multiple physical cells
      vals = valViewType(Kokkos::view_alloc("CellTools::mapToPhysicalFrame::vals", vcprop), basisCardinality, numPoints);
      basis->getValues(vals,
                       refPoints,
                       OPERATOR_VALUE);
      break;
    }
    case 3: {
      // refPoints is (C,P,D): multiple sets of ref. points are mapped to matching number of physical cells.
      //vals = valViewType("CellTools::mapToPhysicalFrame::vals", numCells, basisCardinality, numPoints);
      vals = valViewType(Kokkos::view_alloc("CellTools::mapToPhysicalFrame::vals", vcprop), numCells, basisCardinality, numPoints);
      for (size_type cell=0;cell<numCells;++cell)
        basis->getValues(Kokkos::subdynrankview( vals,      cell, Kokkos::ALL(), Kokkos::ALL() ),
                         Kokkos::subdynrankview( refPoints, cell, Kokkos::ALL(), Kokkos::ALL() ),
                         OPERATOR_VALUE);
      break;
    }
    }
    
    using FunctorType    = FunctorCellTools::F_mapToPhysicalFrame<PhysPointViewType,WorksetType,valViewType>;

    const auto loopSize = physPoints.extent(0)*physPoints.extent(1);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(physPoints, worksetCell, vals) );
  }

  template<typename DeviceType>
  template<typename refSubcellViewType, typename paramViewType>
  void
  CellTools<DeviceType>::
  mapToReferenceSubcell(       refSubcellViewType refSubcellPoints,
                         const paramViewType  paramPoints,
                         const ordinal_type subcellDim,
                         const ordinal_type subcellOrd,
                         const shards::CellTopology parentCell ) {
#ifdef HAVE_INTREPID2_DEBUG
    ordinal_type parentCellDim = parentCell.getDimension();
    INTREPID2_TEST_FOR_EXCEPTION( !hasReferenceCell(parentCell), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::mapToReferenceSubcell): the specified cell topology does not have a reference cell.");

    INTREPID2_TEST_FOR_EXCEPTION( subcellDim < 1 ||
                                  subcellDim > parentCellDim-1, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::mapToReferenceSubcell): method defined only for subcells with dimension greater than 0 and less than the cell dimension");

    INTREPID2_TEST_FOR_EXCEPTION( subcellOrd <  0 ||
                                  subcellOrd >= static_cast<ordinal_type>(parentCell.getSubcellCount(subcellDim)), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::mapToReferenceSubcell): subcell ordinal out of range.");

    // refSubcellPoints is rank-2 (P,D1), D1 = cell dimension
    INTREPID2_TEST_FOR_EXCEPTION( refSubcellPoints.rank() != 2, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::mapToReferenceSubcell): refSubcellPoints must have rank 2.");
    INTREPID2_TEST_FOR_EXCEPTION( refSubcellPoints.extent(1) != parentCell.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::mapToReferenceSubcell): refSubcellPoints dimension (1) does not match to parent cell dimension.");

    // paramPoints is rank-2 (P,D2) with D2 = subcell dimension
    INTREPID2_TEST_FOR_EXCEPTION( paramPoints.rank() != 2, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::mapToReferenceSubcell): paramPoints must have rank 2.");
    INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(paramPoints.extent(1)) != subcellDim, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::mapToReferenceSubcell): paramPoints dimension (1) does not match to subcell dimension.");

    // cross check: refSubcellPoints and paramPoints: dimension 0 must match
    INTREPID2_TEST_FOR_EXCEPTION( refSubcellPoints.extent(0) != paramPoints.extent(0), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::mapToReferenceSubcell): refSubcellPoints dimension (0) does not match to paramPoints dimension(0).");
#endif

    // Get the subcell map, i.e., the coefficients of the parametrization function for the subcell
    const auto subcellMap = RefSubcellParametrization<DeviceType>::get(subcellDim, parentCell.getKey());

    mapToReferenceSubcell(refSubcellPoints, paramPoints, subcellMap, subcellOrd);
  }

  template<typename DeviceType>
  template<typename refSubcellViewType, typename paramViewType>
  void
  CellTools<DeviceType>::
  mapToReferenceSubcell(       refSubcellViewType                                            refSubcellPoints,
                         const paramViewType                                                 paramPoints,
                         const typename RefSubcellParametrization<DeviceType>::ConstViewType subcellParametrization,
                         const ordinal_type subcellOrd) {

    constexpr bool are_accessible =
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(refSubcellPoints)::memory_space>::accessible &&
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(paramPoints)::memory_space>::accessible;
    static_assert(are_accessible, "CellTools<DeviceType>::mapToReferenceSubcell(..): input/output views' memory spaces are not compatible with DeviceType");


  #ifdef HAVE_INTREPID2_DEBUG
    const bool ranks_and_dims_compatible = (refSubcellPoints.rank() == 2) && (paramPoints.rank() == 2) && (subcellParametrization.rank() == 3) && 
                                           (refSubcellPoints.extent(0) == paramPoints.extent(0)) && 
                                           (refSubcellPoints.extent(1) == subcellParametrization.extent(1)) && 
                                           (paramPoints.extent(1) == subcellParametrization.extent(2)-1);
    
    INTREPID2_TEST_FOR_EXCEPTION(!ranks_and_dims_compatible, std::invalid_argument, "CellTools<DeviceType>::mapToReferenceSubcell(..): input/output views' ranks and dimensions are not compatible");
  #endif

    const ordinal_type parentCellDim = subcellParametrization.extent(1);
    const ordinal_type numPts  = paramPoints.extent(0);
    
    //Note: this function has several template parameters and the compiler gets confused if using a lambda function
    using FunctorType = FunctorCellTools::F_mapReferenceSubcell<decltype(refSubcellPoints), decltype(paramPoints), decltype(subcellParametrization)>;

    Kokkos::MDRangePolicy<ExecSpaceType,Kokkos::Rank<2>> rangePolicy({0,0},{numPts,parentCellDim});
    
    // Apply the parametrization map to every point in parameter domain
    Kokkos::parallel_for( rangePolicy, FunctorType(refSubcellPoints, paramPoints, subcellParametrization, subcellOrd) );
  }

  template<typename DeviceType>
  template<typename refSubcellViewType, typename paramViewType, typename ordViewType>
  void
  CellTools<DeviceType>::
  mapToReferenceSubcell(       refSubcellViewType                                            refSubcellPoints,
                         const paramViewType                                                 paramPoints,
                         const typename RefSubcellParametrization<DeviceType>::ConstViewType subcellParametrization,
                         const ordViewType                                                   subcellOrd) {

    constexpr bool are_accessible =
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(refSubcellPoints)::memory_space>::accessible &&
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(paramPoints)::memory_space>::accessible && 
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(subcellOrd)::memory_space>::accessible;
    static_assert(are_accessible, "CellTools<DeviceType>::mapToReferenceSubcell(..): input/output views' memory spaces are not compatible with DeviceType");
    

#ifdef HAVE_INTREPID2_DEBUG
    const bool ranks_and_dims_compatible = (refSubcellPoints.rank() == 3) && (paramPoints.rank() == 2) && (subcellParametrization.rank() == 3) && 
                                           (refSubcellPoints.extent(0) == subcellOrd.extent(0)) && 
                                           (refSubcellPoints.extent(1) == paramPoints.extent(0)) && 
                                           (refSubcellPoints.extent(2) == subcellParametrization.extent(1)) && 
                                           (paramPoints.extent(1) == subcellParametrization.extent(2)-1);

    INTREPID2_TEST_FOR_EXCEPTION(!ranks_and_dims_compatible, std::invalid_argument, "CellTools<DeviceType>::mapToReferenceSubcell(..): input/output views' ranks and dimensions are not compatible");
#endif

    const ordinal_type numSubCells = refSubcellPoints.extent(0);
    const ordinal_type parentCellDim = subcellParametrization.extent(1);
    const ordinal_type numPts  = paramPoints.extent(0);

    
    //Note: this function has several template parameters and the compiler gets confused if using a lambda function
    using FunctorType = FunctorCellTools::F_mapReferenceSubcellBatch<decltype(refSubcellPoints), decltype(paramPoints), decltype(subcellParametrization), decltype(subcellOrd)>;

    Kokkos::MDRangePolicy<ExecSpaceType,Kokkos::Rank<3>> rangePolicy({0,0,0},{numSubCells,numPts,parentCellDim});
    
    // Apply the parametrization map to every point in parameter domain
    Kokkos::parallel_for( rangePolicy, FunctorType(refSubcellPoints, paramPoints, subcellParametrization, subcellOrd) );
  }

}

#endif
