// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_TensorPoints.hpp
    \brief  View-like interface to tensor points; point components are stored separately; the appropriate coordinate is determined from the composite point index and requested dimension at access time.

    \author Nathan V. Roberts
*/

#ifndef Intrepid2_TensorPoints_h
#define Intrepid2_TensorPoints_h

#include <Kokkos_Core.hpp>

namespace Intrepid2 {
/** \class Intrepid2::TensorPoints
    \brief View-like interface to tensor points; point components are stored separately; the appropriate coordinate is determined from the composite point index and requested dimension at access time.
*/
  template<class PointScalar, typename DeviceType>
  class TensorPoints {
  public:
    using value_type = PointScalar;
  protected:
    Kokkos::Array< ScalarView<PointScalar,DeviceType>, Parameters::MaxTensorComponents> pointTensorComponents_; // each component has dimensions (P,D)
    ordinal_type numTensorComponents_;
    ordinal_type totalPointCount_;
    ordinal_type totalDimension_;
    Kokkos::View<ordinal_type*, DeviceType> dimToComponent_;
    Kokkos::View<ordinal_type*, DeviceType> dimToComponentDim_;
    Kokkos::Array<ordinal_type, Parameters::MaxTensorComponents> pointModulus_;
    Kokkos::Array<ordinal_type, Parameters::MaxTensorComponents> pointDivisor_;
    
    bool isValid_;
    using reference_type = typename ScalarView<PointScalar,DeviceType>::reference_type;
    
    void TEST_VALID_POINT_COMPONENTS()
    {
#ifdef HAVE_INTREPID2_DEBUG
      if (isValid_)
      {
        for (ordinal_type r=0; r<numTensorComponents_; r++)
        {
          const auto & pointComponent = pointTensorComponents_[r];
          INTREPID2_TEST_FOR_EXCEPTION(2 != pointComponent.rank(), std::invalid_argument, "Each component must have shape (P,D)");
          INTREPID2_TEST_FOR_EXCEPTION(pointComponent.extent_int(0) <= 0, std::invalid_argument, "Each component must have at least one point");
        }
      }
#endif
    }
    
    /**
     \brief Initialize members based on constructor parameters.
    */
    void initialize()
    {
      TEST_VALID_POINT_COMPONENTS();
      
      totalPointCount_ = 1;
      totalDimension_  = 0;
      for (ordinal_type r=0; r<numTensorComponents_; r++)
      {
        totalPointCount_ *= pointTensorComponents_[r].extent_int(0);
        totalDimension_  += pointTensorComponents_[r].extent_int(1);
      }
      ordinal_type pointDivisor = 1;
      for (ordinal_type r=0; r<numTensorComponents_; r++)
      {
        pointModulus_[r] = pointTensorComponents_[r].extent_int(0);
        pointDivisor_[r] = pointDivisor;
        pointDivisor *= pointTensorComponents_[r].extent_int(0);
      }
      dimToComponent_    = Kokkos::View<ordinal_type*, DeviceType>("dimToComponent_",totalDimension_);
      dimToComponentDim_ = Kokkos::View<ordinal_type*, DeviceType>("dimToComponentDim_",totalDimension_);
      ordinal_type d=0;
      ordinal_type dimsSoFar = 0;

      auto dimToComponentHost = Kokkos::create_mirror_view(dimToComponent_);
      auto dimToComponentDimHost = Kokkos::create_mirror_view(dimToComponentDim_);
      for (ordinal_type r=0; r<numTensorComponents_; r++)
      {
        const int componentDim = pointTensorComponents_[r].extent_int(1);
        for (int i=0; i<componentDim; i++)
        {
          dimToComponentHost[d]    = r;
          dimToComponentDimHost[d] = d - dimsSoFar;
          d++;
        }
        dimsSoFar += componentDim;
      }
      Kokkos::deep_copy(dimToComponent_,dimToComponentHost);
      Kokkos::deep_copy(dimToComponentDim_,dimToComponentDimHost);
    }
  public:
    /**
     \brief Constructor with fixed-length Kokkos::Array argument.
     \param [in] pointTensorComponents - the components representing the points.
     
     TensorPoints has shape (P,D), where P is the product of the first dimensions of the component points, and D is the sum of the second dimensions.
    */
    template<size_t numTensorComponents>
    TensorPoints(Kokkos::Array< ScalarView<PointScalar,DeviceType>, numTensorComponents> pointTensorComponents)
    :
    numTensorComponents_(numTensorComponents),
    isValid_(true)
    {
      for (unsigned r=0; r<numTensorComponents; r++)
      {
        pointTensorComponents_[r] = pointTensorComponents[r];
      }
      
      initialize();
    }
    
    /**
     \brief Constructor that takes a subset of the tensorial components of another points container.
     \param [in] otherPointsContainer - the original points container
     \param [in] whichDims - the tensorial component indices to take from the other container.
     
     \note this does not copy the points.
    */
    TensorPoints(TensorPoints otherPointsContainer, std::vector<int> whichDims)
    :
    numTensorComponents_(whichDims.size()),
    isValid_(true)
    {
      int r = 0;
      for (const auto & componentOrdinal : whichDims)
      {
        pointTensorComponents_[r++] = otherPointsContainer.getTensorComponent(componentOrdinal);
      }
      
      initialize();
    }
    
    /**
     \brief Constructor with variable-length std::vector argument.
     \param [in] pointTensorComponents - the components representing the points.
     
     TensorPoints has shape (P,D), where P is the product of the first dimensions of the component points, and D is the sum of the second dimensions.
    */
    TensorPoints(std::vector< ScalarView<PointScalar,DeviceType>> pointTensorComponents)
    :
    numTensorComponents_(pointTensorComponents.size()),
    isValid_(true)
    {
      for (ordinal_type r=0; r<numTensorComponents_; r++)
      {
        pointTensorComponents_[r] = pointTensorComponents[r];
      }
      
      initialize();
    }
    
    /**
     \brief Constructor for point set with trivial tensor structure.
     \param [in] points - the points, with shape (P,D).
     
     TensorPoints has the same shape (P,D) as the input points.
    */
    TensorPoints(ScalarView<PointScalar,DeviceType> points)
    :
    numTensorComponents_(1),
    isValid_(true)
    {
      pointTensorComponents_[0] = points;
      initialize();
    }

    /**
     \brief Copy from one points container, which may be an arbitrary functor, to a DynRankView.
     \param [in] toPoints - the container to copy to.
     \param [in] fromPoints - the container to copy from.
    */
    template<class OtherPointsContainer>
    void copyPointsContainer(ScalarView<PointScalar,DeviceType> toPoints, OtherPointsContainer fromPoints)
    {
      const int numPoints = fromPoints.extent_int(0);
      const int numDims   = fromPoints.extent_int(1);
      using ExecutionSpace = typename DeviceType::execution_space;
      auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{numPoints,numDims});
      Kokkos::parallel_for("copy points", policy,
      KOKKOS_LAMBDA (const int &i0, const int &i1) {
        toPoints(i0,i1) = fromPoints(i0,i1);
      });
    }
    
    //! copy-like constructor for differing device type, but same memory space.  This does a shallow copy of the underlying view.
    template<typename OtherDeviceType, class = typename std::enable_if< std::is_same<typename DeviceType::memory_space, typename OtherDeviceType::memory_space>::value>::type,
                                       class = typename std::enable_if<!std::is_same<DeviceType,OtherDeviceType>::value>::type>
    TensorPoints(const TensorPoints<PointScalar,OtherDeviceType> &tensorPoints)
    :
    numTensorComponents_(tensorPoints.numTensorComponents()),
    isValid_(tensorPoints.isValid())
    {
      if (isValid_)
      {
        for (ordinal_type r=0; r<numTensorComponents_; r++)
        {
          pointTensorComponents_[r] = tensorPoints.getTensorComponent(r);
        }
        initialize();
      }
    }
    
    //! copy-like constructor for differing memory spaces.  This does a deep_copy of the underlying view.
    template<typename OtherDeviceType, class = typename std::enable_if<!std::is_same<typename DeviceType::memory_space, typename OtherDeviceType::memory_space>::value>::type>
    TensorPoints(const TensorPoints<PointScalar,OtherDeviceType> &tensorPoints)
    :
    numTensorComponents_(tensorPoints.numTensorComponents()),
    isValid_(tensorPoints.isValid())
    {
      if (isValid_)
      {
        for (ordinal_type r=0; r<numTensorComponents_; r++)
        {
          ScalarView<PointScalar,OtherDeviceType> otherPointComponent = tensorPoints.getTensorComponent(r);
          const int numPoints = otherPointComponent.extent_int(0);
          const int numDims   = otherPointComponent.extent_int(1);
          pointTensorComponents_[r] = ScalarView<PointScalar,DeviceType>("Intrepid2 point component", numPoints, numDims);
          
          using MemorySpace = typename DeviceType::memory_space;
          auto pointComponentMirror = Kokkos::create_mirror_view_and_copy(MemorySpace(), otherPointComponent);
          
          copyPointsContainer(pointTensorComponents_[r], pointComponentMirror);
        }
        initialize();
      }
    }
    
    //! Default constructor.  TensorPoints::isValid() will return false.
    TensorPoints() :
    isValid_(false)
    // when constructed with default arguments, TensorPoints should not be used…
    // default constructor is only provided for things like CellGeometry, which has TensorPoints as a member,
    // but only uses it in certain modes.
    {}

    //! Returns the number of points in the indicated component.
    ordinal_type componentPointCount(const ordinal_type &tensorComponentOrdinal) const
    {
      return pointTensorComponents_[tensorComponentOrdinal].extent_int(0);
    }
    
    /**
     \brief Accessor that accepts a composite point index.
     \param [in] tensorPointIndex - the composite point index.
     \param [in] dim - the coordinate dimension.
     From <var>tensorPointIndex</var>, we can determine a point index for each component.  From <var>dim</var>, we can determine which component contains the coordinate of interest, and the coordinate dimension in that component.
     
     \return The component coordinate.
    */
    template <typename iType0, typename iType1>
    KOKKOS_INLINE_FUNCTION typename std::enable_if<
        (std::is_integral<iType0>::value && std::is_integral<iType1>::value),
        reference_type>::type
    operator()(const iType0& tensorPointIndex, const iType1& dim) const {
      const ordinal_type component = dimToComponent_[dim];
      const ordinal_type d = dimToComponentDim_[dim];
      const ordinal_type componentPointOrdinal = (tensorPointIndex / pointDivisor_[component]) % pointModulus_[component];
      return pointTensorComponents_[component](componentPointOrdinal,d);
    }
    
    /**
     \brief Accessor that accepts a a fixed-length array with entries corresponding to component indices.
     \param [in] pointOrdinalComponents - the component point indices.
     \param [in] dim - the coordinate dimension.
     From <var>dim</var>, we can determine which component contains the coordinate of interest, and the coordinate dimension in that component.  We use the indicated point index for that component specified in <var>pointOrdinalComponents</var> to return the appropriate coordinate.
     
     \return The component coordinate.
    */
    template <typename iType0, typename iType1, size_t numTensorComponents>
    KOKKOS_INLINE_FUNCTION typename std::enable_if<
        (std::is_integral<iType0>::value && std::is_integral<iType1>::value),
        reference_type>::type
    operator()(const Kokkos::Array<iType0,numTensorComponents>& pointOrdinalComponents, const iType1& dim) const {
      const ordinal_type component = dimToComponent_[dim];
      const ordinal_type d = dimToComponentDim_[dim];
      const ordinal_type componentPointOrdinal = pointOrdinalComponents[component];
      return pointTensorComponents_[component](componentPointOrdinal,d);
    }
    
    /** \brief  Returns the logical extent in the requested dimension.
       \param [in] r - the dimension
       
     \return the logical extent in the requested dimension, as an int.  (This will be 1 for r >= 2, as this is a rank-2 container.)
    */
    template <typename iType>
    KOKKOS_INLINE_FUNCTION
    typename std::enable_if<std::is_integral<iType>::value, int>::type
    extent_int(const iType& r) const {
      if (r == static_cast<iType>(0))
      {
        return totalPointCount_;
      }
      else if (r == static_cast<iType>(1))
      {
        return totalDimension_;
      }
      else
      {
        return 1;
      }
    }
    
    /** \brief  Returns the logical extent in the requested dimension.
       \param [in] r - the dimension
       
     \return the logical extent in the requested dimension.  (This will be 1 for r >= 2, as this is a rank-2 container.)
    */
    template <typename iType>
    KOKKOS_INLINE_FUNCTION constexpr
    typename std::enable_if<std::is_integral<iType>::value, size_t>::type
    extent(const iType& r) const {
      // C++ prior to 14 doesn't allow constexpr if statements; the compound ternary is here to keep things constexpr
      return (r == static_cast<iType>(0)) ? totalPointCount_
                                          :
             (r == static_cast<iType>(1)) ? totalDimension_ : 1;
    }
    
    //! This method is for compatibility with existing methods that take raw point views.  Note that in general it is probably better for performance to make the existing methods accept a TensorPoints object, since that can dramatically reduce memory footprint, and avoids an allocation here.
    ScalarView<PointScalar,DeviceType> allocateAndFillExpandedRawPointView() const
    {
      const int numPoints = this->extent_int(0);
      const int spaceDim  = this->extent_int(1);
      ScalarView<PointScalar,DeviceType> expandedRawPoints("expanded raw points from TensorPoints", numPoints, spaceDim);
      TensorPoints<PointScalar,DeviceType> tensorPoints(*this); // (shallow) copy for lambda capture
      using ExecutionSpace = typename DeviceType::execution_space;
      Kokkos::parallel_for(
        Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{numPoints,spaceDim}),
        KOKKOS_LAMBDA (const int &pointOrdinal, const int &d) {
          expandedRawPoints(pointOrdinal,d) = tensorPoints(pointOrdinal,d);
      });
      return expandedRawPoints;
    }
    
    /**
     \brief Returns the requested tensor component.
     \param [in] r - the tensor ordinal of the component
     \return the requested tensor component.
    */
    KOKKOS_INLINE_FUNCTION
    ScalarView<PointScalar,DeviceType> getTensorComponent(const ordinal_type &r) const
    {
      return pointTensorComponents_[r];
    }
    
    //! Returns true for containers that have data; false for those that don't (e.g., those that have been constructed by the default constructor).
    KOKKOS_INLINE_FUNCTION
    bool isValid() const
    {
      return isValid_;
    }
    
    //! Returns the number of tensorial components.
    KOKKOS_INLINE_FUNCTION
    ordinal_type numTensorComponents() const
    {
      return numTensorComponents_;
    }
    
    //! Return the rank of the container, which is 2.
    KOKKOS_INLINE_FUNCTION
    constexpr ordinal_type rank() const
    {
      return 2; // shape is (P,D)
    }
    
    // NOTE: the extractTensorPoints() code commented out below appears not to work.  We don't have tests against it, though.
    // TODO: either delete this, or re-enable, add tests, and fix.
//    template<class ViewType>
//    static TensorPoints<PointScalar,DeviceType> extractTensorPoints( ViewType expandedPoints, const std::vector<ordinal_type> &dimensionExtents )
//    {
//      const ordinal_type numComponents = dimensionExtents.size();
//      const ordinal_type numPoints     = expandedPoints.extent_int(0);
//      Kokkos::Array<ordinal_type,Parameters::MaxTensorComponents> componentPointCounts;
//      Kokkos::Array<ordinal_type,Parameters::MaxTensorComponents> componentPointTensorStride;
//      Kokkos::Array<ordinal_type,Parameters::MaxTensorComponents> componentPointDimOffsets;
//
//      // for simplicity of implementation, we copy to host:
//      auto hostExpandedPoints = Kokkos::create_mirror_view_and_copy(typename Kokkos::HostSpace::memory_space(), expandedPoints);
//
//      ordinal_type dimOffset = 0;
//      ordinal_type tensorPointStride = 1; // increases with componentOrdinal
//
//      TensorPoints<PointScalar,DeviceType> invalidTensorData; // will be returned if extraction does not succeed.
//
//      for (ordinal_type componentOrdinal=0; componentOrdinal<numComponents; componentOrdinal++)
//      {
//        componentPointDimOffsets[componentOrdinal]   = dimOffset;
//        componentPointTensorStride[componentOrdinal] = tensorPointStride;
//        const ordinal_type numDimsForComponent = dimensionExtents[componentOrdinal];
//        std::vector<PointScalar> firstPoint(numDimsForComponent);
//        for (ordinal_type d=0; d<numDimsForComponent; d++)
//        {
//          firstPoint[d] = hostExpandedPoints(0,d+dimOffset);
//        }
//
//        // we assume that once we see the same point twice, we've found the cycle length:
//        componentPointCounts[componentOrdinal] = -1;
//        for (ordinal_type pointOrdinal=1; pointOrdinal<numPoints; pointOrdinal += tensorPointStride)
//        {
//          bool matches = true;
//          for (ordinal_type d=0; d<numDimsForComponent; d++)
//          {
//            matches = matches && (firstPoint[d] == hostExpandedPoints(pointOrdinal,d+dimOffset));
//          }
//          if (matches)
//          {
//            componentPointCounts[componentOrdinal] = pointOrdinal;
//            break;
//          }
//        }
//        if (componentPointCounts[componentOrdinal] == -1)
//        {
//          // no matches found -> no tensor decomposition available
//          return invalidTensorData;
//        }
//
//        // check that we got the cycle length correct:
//        for (ordinal_type pointOrdinal=0; pointOrdinal<componentPointCounts[componentOrdinal]; pointOrdinal += tensorPointStride)
//        {
//          std::vector<PointScalar> point(numDimsForComponent);
//          for (ordinal_type d=0; d<numDimsForComponent; d++)
//          {
//            point[d] = hostExpandedPoints(pointOrdinal,d+dimOffset);
//          }
//          // each of the following points should match:
//          for (ordinal_type secondPointOrdinal=0; secondPointOrdinal<numPoints; secondPointOrdinal += tensorPointStride*componentPointCounts[componentOrdinal])
//          {
//            bool matches = true;
//            for (ordinal_type d=0; d<numDimsForComponent; d++)
//            {
//              matches = matches && (point[d] == hostExpandedPoints(secondPointOrdinal,d+dimOffset));
//            }
//            if (!matches)
//            {
//              // fail:
//              return invalidTensorData;
//            }
//          }
//        }
//
//        dimOffset += numDimsForComponent;
//        tensorPointStride *= componentPointCounts[componentOrdinal];
//      }
//
//      std::vector< ScalarView<PointScalar,DeviceType> > componentPoints(numComponents);
//      std::vector< ScalarView<PointScalar,Kokkos::HostSpace> > componentPointsHost(numComponents);
//      for (ordinal_type componentOrdinal=0; componentOrdinal<numComponents; componentOrdinal++)
//      {
//        const ordinal_type numPointsForComponent = componentPointCounts[componentOrdinal];
//        const ordinal_type dimForComponent       = dimensionExtents[componentOrdinal];
//        componentPoints[componentOrdinal] = ScalarView<PointScalar,DeviceType>("extracted tensor components", numPointsForComponent, dimForComponent);
//        componentPointsHost[componentOrdinal] = Kokkos::create_mirror(componentPoints[componentOrdinal]);
//      }
//
//      for (ordinal_type componentOrdinal=0; componentOrdinal<numComponents; componentOrdinal++)
//      {
//        const ordinal_type numComponentPoints = componentPointCounts[componentOrdinal];
//
//        auto hostView = componentPointsHost[componentOrdinal];
//        auto deviceView = componentPoints[componentOrdinal];
//        const ordinal_type tensorPointStride   = componentPointTensorStride[componentOrdinal];
//        const ordinal_type dimOffset           = componentPointDimOffsets[componentOrdinal];
//        const ordinal_type numDimsForComponent = dimensionExtents[componentOrdinal];
//        for (ordinal_type componentPointOrdinal=0; componentPointOrdinal<numComponentPoints; componentPointOrdinal++)
//        {
//          const ordinal_type expandedPointOrdinal = componentPointOrdinal*tensorPointStride;
//          for (ordinal_type d=0; d<numDimsForComponent; d++)
//          {
//            hostView(componentPointOrdinal,d) = hostExpandedPoints(expandedPointOrdinal,d+dimOffset);
//          }
//        }
//        Kokkos::deep_copy(deviceView, hostView);
//      }
//
//      // prior to return, check all points agree in all dimensions with the input points
//      // for convenience, we do this check on host, too
//      TensorPoints<PointScalar,Kokkos::HostSpace> hostTensorPoints(componentPointsHost);
//      const ordinal_type totalDim = expandedPoints.extent_int(1);
//      bool matches = true;
//      for (ordinal_type pointOrdinal=0; pointOrdinal<numPoints; pointOrdinal++)
//      {
//        for (ordinal_type d=0; d<totalDim; d++)
//        {
//          const auto &originalCoord = hostExpandedPoints(pointOrdinal,d);
//          const auto &tensorCoord   = hostTensorPoints(pointOrdinal,d);
//          if (originalCoord != tensorCoord)
//          {
//            matches = false;
//          }
//        }
//      }
//      if (!matches)
//      {
//        return invalidTensorData;
//      }
//
//      return TensorPoints<PointScalar,DeviceType>(componentPoints);
//    }
  };

}

#endif /* Intrepid2_TensorPoints_h */
