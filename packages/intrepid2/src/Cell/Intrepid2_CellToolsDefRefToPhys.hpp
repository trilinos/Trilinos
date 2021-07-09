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

        const auto valRank = _basisVals.rank();
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
    struct F_mapReferenceSubcell2 {
      refSubcellViewType refSubcellPoints_;
      const paramPointsViewType paramPoints_;
      const subcellMapViewType subcellMap_;
      ordinal_type subcellOrd_, dim_;

      KOKKOS_INLINE_FUNCTION
      F_mapReferenceSubcell2( refSubcellViewType  refSubcellPoints,
                             const paramPointsViewType paramPoints,
                             const subcellMapViewType  subcellMap,
                             ordinal_type subcellOrd,
                             ordinal_type dim)
        : refSubcellPoints_(refSubcellPoints), paramPoints_(paramPoints), subcellMap_(subcellMap),
          subcellOrd_(subcellOrd), dim_(dim){};

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type pt) const {

        const auto u = paramPoints_(pt, 0);
        const auto v = paramPoints_(pt, 1);

        // map_dim(u,v) = c_0(dim) + c_1(dim)*u + c_2(dim)*v because both Quad and Tri ref faces are affine!
        for (ordinal_type i=0;i<dim_;++i)
          refSubcellPoints_(pt, i) = subcellMap_(subcellOrd_, i, 0) + ( subcellMap_(subcellOrd_, i, 1)*u +
                                                                 subcellMap_(subcellOrd_, i, 2)*v );
      }
    };

    template<typename refSubcellViewType,
             typename paramPointsViewType,
             typename subcellMapViewType>
    struct F_mapReferenceSubcell1 {
      refSubcellViewType refSubcellPoints_;
      const paramPointsViewType paramPoints_;
      const subcellMapViewType subcellMap_;
      ordinal_type subcellOrd_, dim_;

      KOKKOS_INLINE_FUNCTION
      F_mapReferenceSubcell1( refSubcellViewType  refSubcellPoints,
                             const paramPointsViewType paramPoints,
                             const subcellMapViewType  subcellMap,
                             ordinal_type subcellOrd,
                             ordinal_type dim)
        : refSubcellPoints_(refSubcellPoints), paramPoints_(paramPoints), subcellMap_(subcellMap),
          subcellOrd_(subcellOrd), dim_(dim){};

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type pt) const {
        const auto u = paramPoints_(pt, 0);
        for (ordinal_type i=0;i<dim_;++i)
          refSubcellPoints_(pt, i) = subcellMap_(subcellOrd_, i, 0) + ( subcellMap_(subcellOrd_, i, 1)*u );
      }
    };
  }


/*
  template<typename DeviceType>
  template<typename physPointValueType,   class ...physPointProperties,
           typename refPointValueType,    class ...refPointProperties,
           typename worksetCellValueType, class ...worksetCellProperties>
  void
  CellTools<DeviceType>::
  mapToPhysicalFrame(       Kokkos::DynRankView<physPointValueType,physPointProperties...>     physPoints,
                      const Kokkos::DynRankView<refPointValueType,refPointProperties...>       refPoints,
                      const Kokkos::DynRankView<worksetCellValueType,worksetCellProperties...> worksetCell,
                      const shards::CellTopology cellTopo ) {

   auto basis = createHGradBasis<refPointValueType,refPointValueType>(cellTopo);
   mapToPhysicalFrame(physPoints,
                      refPoints,
                      worksetCell,
                      basis);
  }
*/

  template<typename DeviceType>
  template<typename physPointValueType,   class ...physPointProperties,
           typename refPointValueType,    class ...refPointProperties,
           typename WorksetType,
           typename HGradBasisPtrType>
  void
  CellTools<DeviceType>::
  mapToPhysicalFrame(      Kokkos::DynRankView<physPointValueType,physPointProperties...>     physPoints,
                     const Kokkos::DynRankView<refPointValueType,refPointProperties...>       refPoints,
                     const WorksetType worksetCell,
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

    using physPointViewType =Kokkos::DynRankView<physPointValueType,physPointProperties...>;
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
    
    using FunctorType    = FunctorCellTools::F_mapToPhysicalFrame<physPointViewType,WorksetType,valViewType>;

    const auto loopSize = physPoints.extent(0)*physPoints.extent(1);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(physPoints, worksetCell, vals) );
  }

  template<typename DeviceType>
  template<typename refSubcellPointValueType, class ...refSubcellPointProperties,
           typename paramPointValueType, class ...paramPointProperties>
  void
  CellTools<DeviceType>::
  mapToReferenceSubcell(       Kokkos::DynRankView<refSubcellPointValueType,refSubcellPointProperties...> refSubcellPoints,
                         const Kokkos::DynRankView<paramPointValueType,paramPointProperties...>           paramPoints,
                         const ordinal_type subcellDim,
                         const ordinal_type subcellOrd,
                         const shards::CellTopology parentCell ) {
    ordinal_type parentCellDim = parentCell.getDimension();
#ifdef HAVE_INTREPID2_DEBUG
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
    INTREPID2_TEST_FOR_EXCEPTION( refSubcellPoints.extent(0) < paramPoints.extent(0), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::mapToReferenceSubcell): refSubcellPoints dimension (0) does not match to paramPoints dimension(0).");
#endif

    constexpr bool are_accessible =
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(refSubcellPoints)::memory_space>::accessible &&
        Kokkos::Impl::MemorySpaceAccess<MemSpaceType,
        typename decltype(paramPoints)::memory_space>::accessible;

    static_assert(are_accessible, "CellTools<DeviceType>::mapToReferenceSubcell(..): input/output views' memory spaces are not compatible with DeviceType");

    // Get the subcell map, i.e., the coefficients of the parametrization function for the subcell
    const auto subcellMap = RefSubcellParametrization<DeviceType>::get(subcellDim, parentCell.getKey());
    
    const ordinal_type numPts  = paramPoints.extent(0);
    
    //Note: this function has several template parameters and the compiler gets confused if using a lambda function
    using FunctorType1 = FunctorCellTools::F_mapReferenceSubcell1<decltype(refSubcellPoints), decltype(paramPoints), decltype(subcellMap)>;
    using FunctorType2 = FunctorCellTools::F_mapReferenceSubcell2<decltype(refSubcellPoints), decltype(paramPoints), decltype(subcellMap)>;

    Kokkos::RangePolicy<ExecSpaceType> policy(0, numPts);
    // Apply the parametrization map to every point in parameter domain
    switch (subcellDim) {
    case 2: {
      Kokkos::parallel_for( policy, FunctorType2(refSubcellPoints, paramPoints, subcellMap, subcellOrd, parentCellDim) );
      break;
    }
    case 1: {
      Kokkos::parallel_for( policy, FunctorType1(refSubcellPoints, paramPoints, subcellMap, subcellOrd, parentCellDim) );
      break;
    }
    default: {}
    }
  }
}

#endif
