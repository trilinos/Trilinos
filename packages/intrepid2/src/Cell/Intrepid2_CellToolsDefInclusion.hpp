// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file   Intrepid2_CellToolsDefInclusion.hpp
    \brief  Definition file for point inclusion functions of the Intrepid2::CellTools class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/
#ifndef __INTREPID2_CELLTOOLS_DEF_INCLUSION_HPP__
#define __INTREPID2_CELLTOOLS_DEF_INCLUSION_HPP__

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

  
  template<typename DeviceType>
  template<typename PointViewType>
  bool 
  CellTools<DeviceType>::
  checkPointInclusion( const PointViewType          point,
                       const shards::CellTopology   cellTopo,
                       const typename ScalarTraits<typename PointViewType::value_type>::scalar_type threshold) {
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( point.rank() != 1, std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::checkPointInclusion): Point must have rank 1. ");
    INTREPID2_TEST_FOR_EXCEPTION( point.extent(0) != cellTopo.getDimension(), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::checkPointInclusion): Point and cell dimensions do not match. ");
#endif
    bool testResult = true;

    // A cell with extended topology has the same reference cell as a cell with base topology. 
    // => testing for inclusion in a reference Triangle<> and a reference Triangle<6> relies on 
    // on the same set of inequalities. To eliminate unnecessary cases we switch on the base topology
    const auto key = cellTopo.getBaseKey();
    switch (key) {
    
    case shards::Line<>::key :
      testResult = PointInclusion<shards::Line<>::key>::check(point, threshold);
      break;
    case shards::Triangle<>::key :
      testResult = PointInclusion<shards::Triangle<>::key>::check(point, threshold);
      break;      
    case shards::Quadrilateral<>::key :
      testResult = PointInclusion<shards::Quadrilateral<>::key>::check(point, threshold);
      break;      
    case shards::Tetrahedron<>::key :
      testResult = PointInclusion<shards::Tetrahedron<>::key>::check(point, threshold);
      break;
    case shards::Hexahedron<>::key :
      testResult = PointInclusion<shards::Hexahedron<>::key>::check(point, threshold);
      break;
    case shards::Wedge<>::key :
      testResult = PointInclusion<shards::Wedge<>::key>::check(point, threshold);
      break;
    case shards::Pyramid<>::key : 
      testResult = PointInclusion<shards::Pyramid<>::key>::check(point, threshold);
      break;      
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



  template<unsigned cellTopologyKey,
             typename OutputViewType,
             typename InputViewType>
  struct checkPointInclusionFunctor {
    OutputViewType output_;
    InputViewType input_;
    using ScalarType = typename ScalarTraits<typename InputViewType::value_type>::scalar_type;
    ScalarType threshold_;

    KOKKOS_INLINE_FUNCTION
    checkPointInclusionFunctor(       OutputViewType                      output,
                               const  InputViewType                       input,
                               const  ScalarType threshold)
      : output_(output), 
        input_(input),
        threshold_(threshold) {}

    KOKKOS_INLINE_FUNCTION
    void
    operator()(const ordinal_type i) const {
      const auto in = Kokkos::subview(input_,i,Kokkos::ALL());
      const auto check = PointInclusion<cellTopologyKey>::check(in, threshold_);
      output_(i) = check;        
    }
    
    KOKKOS_INLINE_FUNCTION
    void
    operator()(const ordinal_type i, const ordinal_type j) const {
      const auto in = Kokkos::subview(input_,i,j,Kokkos::ALL());
      const auto check = PointInclusion<cellTopologyKey>::check(in, threshold_);
      output_(i,j) = check;        
    }
  };


  template<typename DeviceType>
  template<unsigned cellTopologyKey,
           typename OutputViewType,
           typename InputViewType>
  void CellTools<DeviceType>::
  checkPointwiseInclusion(      OutputViewType inCell, 
                          const InputViewType  points,
                          const typename ScalarTraits<typename InputViewType::value_type>::scalar_type threshold) {     

     using FunctorType = checkPointInclusionFunctor<cellTopologyKey,decltype(inCell),decltype(points)>;
    if (points.rank() == 2) {     // inCell.rank() == 1
      Kokkos::RangePolicy<typename DeviceType::execution_space> policy(0, points.extent(0));
      Kokkos::parallel_for(policy, FunctorType(inCell, points, threshold));
    } else { //points.rank() == 3, inCell.rank() == 2
      Kokkos::MDRangePolicy<typename DeviceType::execution_space,Kokkos::Rank<2>> policy({0,0},{points.extent(0),points.extent(1)});
      Kokkos::parallel_for(policy, FunctorType(inCell, points, threshold));
    }
  }


  template<typename DeviceType>
  template<typename InCellViewType,
           typename InputViewType>
  void
  CellTools<DeviceType>::
  checkPointwiseInclusion(       InCellViewType         inCell,
                           const InputViewType          points,
                           const shards::CellTopology   cellTopo,
                           const typename ScalarTraits<typename InputViewType::value_type>::scalar_type threshold ) {
#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( (inCell.rank() != 1) && (inCell.rank() != 2), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): InCell must have rank 1 or 2. ");
      INTREPID2_TEST_FOR_EXCEPTION( points.extent(points.rank()-1) != cellTopo.getDimension(), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::checkPointInclusion): Points and cell dimensions do not match. ");
      INTREPID2_TEST_FOR_EXCEPTION( inCell.rank() != (points.rank()-1), std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): rank difference between inCell and points is 1.");  
      const ordinal_type iend = inCell.rank();
      
      for (ordinal_type i=0;i<iend;++i) {
        INTREPID2_TEST_FOR_EXCEPTION( inCell.extent(i) != points.extent(i), std::invalid_argument, 
                                      ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): dimension mismatch between inCell and points. " );  
      }
    }
#endif

   const auto key = cellTopo.getBaseKey();
   switch (key) {
    
    case shards::Line<>::key :
      checkPointwiseInclusion<shards::Line<>::key,decltype(inCell),decltype(points)>(inCell, points, threshold);
      break;

    case shards::Triangle<>::key :
      checkPointwiseInclusion<shards::Triangle<>::key,decltype(inCell),decltype(points)>(inCell, points, threshold);
      break;

    case shards::Quadrilateral<>::key :
      checkPointwiseInclusion<shards::Quadrilateral<>::key,decltype(inCell),decltype(points)>(inCell, points, threshold);
      break;

    case shards::Tetrahedron<>::key :
      checkPointwiseInclusion<shards::Tetrahedron<>::key,decltype(inCell),decltype(points)>(inCell, points, threshold);
      break;

    case shards::Hexahedron<>::key :
      checkPointwiseInclusion<shards::Hexahedron<>::key,decltype(inCell),decltype(points)>(inCell, points, threshold);
      break;

    case shards::Wedge<>::key :
      checkPointwiseInclusion<shards::Wedge<>::key,decltype(inCell),decltype(points)>(inCell, points, threshold);
      break;

    case shards::Pyramid<>::key :
      checkPointwiseInclusion<shards::Pyramid<>::key,decltype(inCell),decltype(points)>(inCell, points, threshold);
      break;
      
    default:
      INTREPID2_TEST_FOR_EXCEPTION( false,
                                    std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::checkPointInclusion): Invalid cell topology. ");
    }
  }

  template<typename DeviceType>
  template<typename inCellValueType, class ...inCellProperties,
           typename pointValueType, class ...pointProperties,
           typename cellWorksetValueType, class ...cellWorksetProperties>
  void
  CellTools<DeviceType>::
  checkPointwiseInclusion(       Kokkos::DynRankView<inCellValueType,inCellProperties...> inCell,
                           const Kokkos::DynRankView<pointValueType,pointProperties...> points,
                           const Kokkos::DynRankView<cellWorksetValueType,cellWorksetProperties...> cellWorkset,
                           const shards::CellTopology cellTopo,
                           const typename ScalarTraits<pointValueType>::scalar_type threshold ) {
#ifdef HAVE_INTREPID2_DEBUG
    {
      const auto key = cellTopo.getBaseKey();
      INTREPID2_TEST_FOR_EXCEPTION( key != shards::Line<>::key &&
                                    key != shards::Triangle<>::key &&
                                    key != shards::Quadrilateral<>::key &&
                                    key != shards::Tetrahedron<>::key &&                                                                                                                                               
                                    key != shards::Hexahedron<>::key &&                                                                                                                                                
                                    key != shards::Wedge<>::key &&                                                                                                                                                     
                                    key != shards::Pyramid<>::key, 
                                    std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): cell topology not supported");
      INTREPID2_TEST_FOR_EXCEPTION( inCell.rank() != 2, std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): InCell must have rank 2. ");
      INTREPID2_TEST_FOR_EXCEPTION( cellWorkset.rank() != 3, std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::checkPointwiseInclusion): cellWorkset must have rank 3. ");
      INTREPID2_TEST_FOR_EXCEPTION( cellWorkset.extent(2) != cellTopo.getDimension(), std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::checkPointInclusion): cellWorkset and points have incompatible dimensions. ");
      INTREPID2_TEST_FOR_EXCEPTION( (points.rank() == 3) && (points.extent(0) != cellWorkset.extent(0)) , std::invalid_argument,
                                    ">>> ERROR (Intrepid2::CellTools::checkPointInclusion): cellWorkset and points have incompatible dimensions. ");
    }
#endif    
    const ordinal_type 
      numCells = cellWorkset.extent(0),
      numPoints = points.extent(points.rank()-2), 
      spaceDim = cellTopo.getDimension();

    using result_layout = typename DeduceLayout< decltype(points) >::result_layout;
    auto vcprop = Kokkos::common_view_alloc_prop(points);
    using common_value_type = typename decltype(vcprop)::value_type;
    Kokkos::DynRankView< common_value_type, result_layout, DeviceType > refPoints ( Kokkos::view_alloc("CellTools::checkPointwiseInclusion::refPoints", vcprop), numCells, numPoints, spaceDim);
    
    // expect refPoints(CPD), points (CPD or PD), cellWorkset(CND) 
    if(points.rank() == 3)  
      mapToReferenceFrame(refPoints, points, cellWorkset, cellTopo);
    else { //points.rank() == 2
      Kokkos::DynRankView< common_value_type, result_layout, DeviceType > cellPoints ( Kokkos::view_alloc("CellTools::checkPointwiseInclusion::physCellPoints", vcprop), numCells, numPoints, spaceDim);
      RealSpaceTools<DeviceType>::clone(cellPoints,points);
      mapToReferenceFrame(refPoints, cellPoints, cellWorkset, cellTopo);
    }    
    checkPointwiseInclusion(inCell, refPoints, cellTopo, threshold);
  }

}

#endif
