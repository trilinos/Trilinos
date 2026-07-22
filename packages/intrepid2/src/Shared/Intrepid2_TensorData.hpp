// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_TensorData.hpp
    \brief  View-like interface to tensor data; tensor components are stored separately and multiplied together at access time.

    \author Nathan V. Roberts
*/

#ifndef Intrepid2_TensorData_h
#define Intrepid2_TensorData_h

#include "Intrepid2_Data.hpp"

#include "Intrepid2_ScalarView.hpp"
#include "Intrepid2_Types.hpp"

namespace Intrepid2
{
  /** \class Intrepid2::TensorData
      \brief View-like interface to tensor data; tensor components are stored separately and multiplied together at access time.
  */
  template<class Scalar, typename DeviceType>
  class TensorData {
  public:
    using value_type = Scalar;
    using execution_space = typename DeviceType::execution_space;
  protected:
    Kokkos::Array< Data<Scalar,DeviceType>, Parameters::MaxTensorComponents> tensorComponents_;
    Kokkos::Array<ordinal_type, 7> extents_;
    Kokkos::Array<Kokkos::Array<ordinal_type, Parameters::MaxTensorComponents>, 7> entryModulus_;
    ordinal_type rank_;
    bool separateFirstComponent_ = false; // true supported only for rank 1 components; uses a rank-2 operator() in that case, where the first argument corresponds to the index in the first component, while the second argument corresponds to the tensorially-multiplied remaining components
    ordinal_type numTensorComponents_ = 0;
    
    /**
     \brief Initialize members based on constructor parameters.
    */
    void initialize()
    {
      ordinal_type maxComponentRank = -1;
      for (ordinal_type r=0; r<numTensorComponents_; r++)
      {
        const ordinal_type componentRank = tensorComponents_[r].rank();
        maxComponentRank = (maxComponentRank > componentRank) ? maxComponentRank : componentRank;
      }
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(separateFirstComponent_ && (maxComponentRank != 1), std::invalid_argument, "separateFirstComponent = true only supported if all components have rank 1");
      ordinal_type rd_start; // where to begin the extents_ and entryModulus_ loops below
      if ((maxComponentRank == 1) && separateFirstComponent_)
      {
        rank_ = 2;
        rd_start = 1;
        extents_[0] = tensorComponents_[0].extent_int(0);
        entryModulus_[0][0] = extents_[0]; // should not be used
      }
      else
      {
        rank_ = maxComponentRank;
        rd_start = 0;
      }
      
      for (ordinal_type d=rd_start; d<7; d++)
      {
        extents_[d] = 1;
        for (ordinal_type r=rd_start; r<numTensorComponents_; r++)
        {
          extents_[d] *= tensorComponents_[r].extent_int(d-rd_start);
        }
        ordinal_type entryModulus = extents_[d];
        for (ordinal_type r=rd_start; r<numTensorComponents_; r++)
        {
          entryModulus /= tensorComponents_[r].extent_int(d-rd_start);
          entryModulus_[d][r] = entryModulus;
        }
      }
    }
  public:
    /**
     \brief Constructor with fixed-length Kokkos::Array argument.
     \param [in] tensorComponents - the data components that will be multiplied together.  May not have more than Parameters::MaxTensorComponents entries.
     \param [in] separateFirstComponent - if true, indicates that the first component will be indexed separately (this is used when the first index corresponds to a cell ordinal)
     
     When <var>separateFirstComponent</var> is false, TensorData has rank equal to the maximum rank of the components in <var>tensorComponents</var>, and the logical index in rank <var>r</var> is a function of the indices in rank <var>r</var> of its components, where the function is such that the fastest-moving component index is the one for the final component.  Components that have rank less than <var>r</var> are understood to have index 0 in that dimension.
     
    When <var>separateFirstComponent</var> is true, all components are required to have rank 1, and TensorData has rank 2, with the first argument reserved for the first component.  The second argument is indexed precisely as described above, omitting the first component.
    */
    template<size_t numTensorComponents>
    TensorData(Kokkos::Array< Data<Scalar,DeviceType>, numTensorComponents> tensorComponents, bool separateFirstComponent = false)
    :
    separateFirstComponent_(separateFirstComponent),
    numTensorComponents_(numTensorComponents)
    {
      for (size_t r=0; r< numTensorComponents; r++)
      {
        tensorComponents_[r] = tensorComponents[r];
      }
      
      initialize();
    }
    
    /**
     \brief Constructor with variable-length std::vector containing the components.
     \param [in] tensorComponents - the data components that will be multiplied together.  May not have more than Parameters::MaxTensorComponents entries.
     \param [in] separateFirstComponent - if true, indicates that the first component will be indexed separately (this is used when the first index corresponds to a cell ordinal)
     
     When <var>separateFirstComponent</var> is false, TensorData has rank equal to the maximum rank of the components in <var>tensorComponents</var>, and the logical index in rank <var>r</var> is a function of the indices in rank <var>r</var> of its components, where the function is such that the fastest-moving component index is the one for the final component.  Components that have rank less than <var>r</var> are understood to have index 0 in that dimension.
     
    When <var>separateFirstComponent</var> is true, all components are required to have rank 1, and TensorData has rank 2, with the first argument reserved for the first component.  The second argument is indexed precisely as described above, omitting the first component.
    */
    TensorData(std::vector< Data<Scalar,DeviceType> > tensorComponents, bool separateFirstComponent = false)
    :
    separateFirstComponent_(separateFirstComponent),
    numTensorComponents_(tensorComponents.size())
    {
      for (ordinal_type r=0; r<numTensorComponents_; r++)
      {
        tensorComponents_[r] = tensorComponents[r];
      }
      
      initialize();
    }
    
    /**
     \brief Constructor to combine two other TensorData objects
     \param [in] first - TensorData object with the components for the first dimension(s)
     \param [in] second - TensorData object with the components for the remaining dimension(s).
     \param [in] separateFirstComponent - if true, indicates that the first component (from the first TensorData object) will be indexed separately (this is used when the first index corresponds to a cell ordinal)
     
       When <var>separateFirstComponent</var> is true, all components are required to have rank 1, and TensorData has rank 2, with the first argument reserved for the first component.  The second argument is indexed precisely as described above, omitting the first component.
    */
    TensorData(const TensorData &first, const TensorData &second, bool separateFirstComponent = false)
    :
    separateFirstComponent_(separateFirstComponent),
    numTensorComponents_(first.numTensorComponents() + second.numTensorComponents())
    {
      ordinal_type r = 0;
      for (ordinal_type r1=0; r1<first.numTensorComponents(); r1++, r++)
      {
        tensorComponents_[r] = first.getTensorComponent(r1);
      }
      for (ordinal_type r2=0; r2<second.numTensorComponents(); r2++, r++)
      {
        tensorComponents_[r] = second.getTensorComponent(r2);
      }
      
      initialize();
    }
    
    /**
     \brief Simple constructor for the case of trivial tensor-product structure (single component)
     \param [in] tensorComponent - the data component.
     
      Simple constructor for trivial tensor-product structure.  The TensorData object will have precisely the same logical data layout as the provided <var>tensorComponent</var>.
    */
    TensorData(Data<Scalar,DeviceType> tensorComponent)
    :
    TensorData(Kokkos::Array< Data<Scalar,DeviceType>, 1>({tensorComponent}), false)
    {}
    
    /**
     \brief Default constructor.
     
      Default constructor provided to allow an indication of empty/zero data.  TensorData::isValid() will return false.
    */
    TensorData()
    :
    extents_({0,0,0,0,0,0,0}),
    rank_(0)
    {}
    
    /**
     \brief Constructor that takes a subset of the tensorial components of another TensorData container.
     \param [in] otherTensorData - the original TensorData container
     \param [in] whichComps - the tensorial component indices to take from the other container.
     
     \note this does not copy the data.
    */
    TensorData(TensorData otherTensorData, std::vector<int> whichComps)
    :
    numTensorComponents_(whichComps.size())
    {
      int r = 0;
      for (const auto & componentOrdinal : whichComps)
      {
        tensorComponents_[r++] = otherTensorData.getTensorComponent(componentOrdinal);
      }
      
      initialize();
    }
    
    //! copy-like constructor for differing device type, but same memory space.  This does a shallow copy of the underlying view.
    template<typename OtherDeviceType, class = typename std::enable_if< std::is_same<typename DeviceType::memory_space, typename OtherDeviceType::memory_space>::value>::type,
                                       class = typename std::enable_if<!std::is_same<DeviceType,OtherDeviceType>::value>::type>
    TensorData(const TensorData<Scalar,OtherDeviceType> &tensorData)
    {
      if (tensorData.isValid())
      {
        numTensorComponents_ = tensorData.numTensorComponents();
        for (ordinal_type r=0; r<numTensorComponents_; r++)
        {
          Data<Scalar,OtherDeviceType> otherTensorComponent = tensorData.getTensorComponent(r);
          tensorComponents_[r] = Data<Scalar,DeviceType>(otherTensorComponent);
        }
        initialize();
      }
      else
      {
        extents_ = Kokkos::Array<ordinal_type,7>{0,0,0,0,0,0,0};
        rank_    = 0;
      }
    }
    
    /**
     \brief Copy-like constructor for differing execution spaces.  This performs a deep copy of the underlying data.
    */
    template<typename OtherDeviceType, class = typename std::enable_if<!std::is_same<typename DeviceType::memory_space, typename OtherDeviceType::memory_space>::value>::type>
    TensorData(const TensorData<Scalar,OtherDeviceType> &tensorData)
    {
      if (tensorData.isValid())
      {
        numTensorComponents_ = tensorData.numTensorComponents();
        for (ordinal_type r=0; r<numTensorComponents_; r++)
        {
          Data<Scalar,OtherDeviceType> otherTensorComponent = tensorData.getTensorComponent(r);
          tensorComponents_[r] = Data<Scalar,DeviceType>(otherTensorComponent);
        }
        initialize();
      }
      else
      {
        extents_ = Kokkos::Array<ordinal_type,7>{0,0,0,0,0,0,0};
        rank_    = 0;
      }
    }
    
    /**
     \brief Returns the requested tensor component.
     \param [in] r - the tensor ordinal of the component
     \return the requested tensor component.
    */
    KOKKOS_INLINE_FUNCTION
    const Data<Scalar,DeviceType> & getTensorComponent(const ordinal_type &r) const
    {
      return tensorComponents_[r];
    }
    
    /**
     \brief Accessor for rank-1 objects.
     \param [in] tensorEntryIndex - the composite entry index.
     \return The product of tensorComponent values with component entry indices corresponding to <var>tensorEntryIndex</var>.
    */
    template <typename iType0>
    KOKKOS_INLINE_FUNCTION typename std::enable_if<std::is_integral<iType0>::value, Scalar>::type
    operator()(const iType0& tensorEntryIndex) const {
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(rank_ != 1, std::invalid_argument, "This method is only valid for rank 1 containers.");
#endif
      Scalar value = 1.0;
      iType0 remainingEntryOrdinal = tensorEntryIndex;
      for (ordinal_type r=0; r<numTensorComponents_; r++)
      {
        const ordinal_type componentEntryCount   = tensorComponents_[r].extent_int(0);
        const ordinal_type componentEntryOrdinal = remainingEntryOrdinal % componentEntryCount;
        remainingEntryOrdinal /= componentEntryCount;
        
        value *= tensorComponents_[r](componentEntryOrdinal);
      }
      
      return value;
    }
    
    /**
     \brief Accessor that accepts a fixed-length array with entries corresponding to component indices.
     \param [in] entryComponents - an array with one entry per tensorial component, each entry indicating the requested index into that component.
     \return The product of tensor component values with the specified component entry indices.
    */
    template <typename iType0, ordinal_type numTensorComponents>
    KOKKOS_INLINE_FUNCTION typename std::enable_if<std::is_integral<iType0>::value, Scalar>::type
    operator()(const Kokkos::Array<iType0,numTensorComponents>& entryComponents) const {
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(rank_ != 1, std::invalid_argument, "This method is only valid for rank 1 containers.");
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(numTensorComponents_ != numTensorComponents, std::invalid_argument, "Tensorial component count mismatch");
#endif
      Scalar value = 1.0;
      for (ordinal_type r=0; r<numTensorComponents; r++)
      {
        value *= tensorComponents_[r](entryComponents[r]);
      }
      return value;
    }
    
    /**
     \brief Accessor for rank-2 objects.
     \param [in] tensorEntryIndex0 - the composite tensor index in the first dimension
     \param [in] tensorEntryIndex1 - the composite tensor index in the second dimension
     If constructed with separateFirstComponent = true, <var>tensorEntryIndex0</var> corresponds to the index to be used for the first component.  Otherwise, it corresponds to an enumeration of all valid combinations of first arguments to the tensorial components.
      
     If constructed with separateFirstComponent = true, <var>tensorEntryIndex1</var> corresponds to an enumeration of all valid combinations of second arguments to components after the first.  Otherwise, it corresponds to an enumeration of all valid combinations of second arguments to all the tensorial components.
     
     \return The product of tensor component values with the specified component entry indices.
    */
    template <typename iType0, typename iType1>
    KOKKOS_INLINE_FUNCTION typename std::enable_if<
        (std::is_integral<iType0>::value && std::is_integral<iType1>::value),
        Scalar>::type
    operator()(const iType0& tensorEntryIndex0, const iType1& tensorEntryIndex1) const {
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(rank_ != 2, std::invalid_argument, "This method is only valid for rank 2 containers.");
#endif
      
      if (numTensorComponents_ == 1)
      {
        return tensorComponents_[0](tensorEntryIndex0,tensorEntryIndex1);
      }
      
      if (!separateFirstComponent_)
      {
        Scalar value = 1.0;
        iType0 remainingEntryOrdinal0 = tensorEntryIndex0;
        iType1 remainingEntryOrdinal1 = tensorEntryIndex1;
        for (ordinal_type r=0; r<numTensorComponents_; r++)
        {
          auto & component = tensorComponents_[r];
          const ordinal_type componentEntryCount0   = component.extent_int(0);
          const ordinal_type componentEntryCount1   = component.extent_int(1);
          const iType0 componentEntryOrdinal0 = remainingEntryOrdinal0 % componentEntryCount0;
          const iType1 componentEntryOrdinal1 = remainingEntryOrdinal1 % componentEntryCount1;
          remainingEntryOrdinal0 /= componentEntryCount0;
          remainingEntryOrdinal1 /= componentEntryCount1;
          
          const ordinal_type componentRank = component.rank();
          
          if (componentRank == 2)
          {
            value *= component(componentEntryOrdinal0,componentEntryOrdinal1);
          }
          else if (componentRank == 1)
          {
            value *= component(componentEntryOrdinal0);
          }
          else
          {
            INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::logic_error, "unsupported component rank encountered");
          }
        }
        
        return value;
      }
      else
      {
        Scalar value = tensorComponents_[0](tensorEntryIndex0);
        iType0 remainingEntryOrdinal = tensorEntryIndex1;
        for (ordinal_type r=1; r<numTensorComponents_; r++)
        {
          const ordinal_type componentEntryCount   = tensorComponents_[r].extent_int(0);
          const ordinal_type componentEntryOrdinal = remainingEntryOrdinal % componentEntryCount;
          remainingEntryOrdinal /= componentEntryCount;
          
          value *= tensorComponents_[r](componentEntryOrdinal);
        }
        return value;
      }
    }
    
    //! return the index into the specified tensorial component in the dimension specified corresponding to the enumerationIndex given for that dimension.
    KOKKOS_INLINE_FUNCTION
    ordinal_type getTensorComponentIndex(const ordinal_type &tensorComponent, const ordinal_type &dim, const ordinal_type &enumerationIndex) const
    {
      ordinal_type remainingEntryOrdinal = enumerationIndex;
      for (ordinal_type r=0; r<tensorComponent; r++)
      {
        const auto & component = tensorComponents_[r];
        const ordinal_type & componentEntryCount = component.extent_int(dim);
        
        remainingEntryOrdinal /= componentEntryCount;
      }
      return remainingEntryOrdinal % tensorComponents_[tensorComponent].extent_int(dim);
    }
    
    /**
     \brief Accessor for rank-3 objects.
     \param [in] tensorEntryIndex0 - the composite tensor index in the first dimension
     \param [in] tensorEntryIndex1 - the composite tensor index in the second dimension
     \param [in] tensorEntryIndex2 - the composite tensor index in the third dimension
     <var>tensorEntryIndex0</var> corresponds to an enumeration of all valid combinations of first arguments to the tensorial components.  Similarly, <var>tensorEntryIndex1</var> and <var>tensorEntryIndex2</var> correspond to an enumeration of all valid second and third arguments to the tensorial components.
     
     \return The product of tensor component values with the specified component entry indices.
    */
    template <typename iType0, typename iType1, typename iType2>
    KOKKOS_INLINE_FUNCTION typename std::enable_if<
        (std::is_integral<iType0>::value && std::is_integral<iType1>::value && std::is_integral<iType2>::value),
        Scalar>::type
    operator()(const iType0& tensorEntryIndex0, const iType1& tensorEntryIndex1, const iType2& tensorEntryIndex2) const
    {
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(rank_ != 3, std::invalid_argument, "This method is only valid for rank 3 containers.");
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(separateFirstComponent_, std::logic_error, "This method should never be called when separateFirstComponent is true");
#endif
      
      Scalar value = 1.0;
      Kokkos::Array<ordinal_type,3> remainingEntryOrdinal {tensorEntryIndex0, tensorEntryIndex1, tensorEntryIndex2};
      for (ordinal_type r=0; r<numTensorComponents_; r++)
      {
        auto & component = tensorComponents_[r];
        const ordinal_type componentEntryCount0   = component.extent_int(0);
        const ordinal_type componentEntryCount1   = component.extent_int(1);
        const ordinal_type componentEntryCount2   = component.extent_int(2);
        const ordinal_type componentEntryOrdinal0 = remainingEntryOrdinal[0] % componentEntryCount0;
        const ordinal_type componentEntryOrdinal1 = remainingEntryOrdinal[1] % componentEntryCount1;
        const ordinal_type componentEntryOrdinal2 = remainingEntryOrdinal[2] % componentEntryCount2;
        remainingEntryOrdinal[0] /= componentEntryCount0;
        remainingEntryOrdinal[1] /= componentEntryCount1;
        remainingEntryOrdinal[2] /= componentEntryCount2;
        
        const ordinal_type componentRank = component.rank();
        
        if (componentRank == 3)
        {
          value *= component(componentEntryOrdinal0,componentEntryOrdinal1,componentEntryOrdinal2);
        }
        else if (componentRank == 2)
        {
          value *= component(componentEntryOrdinal0,componentEntryOrdinal1);
        }
        else if (componentRank == 1)
        {
          value *= component(componentEntryOrdinal0);
        }
        else
        {
          INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::logic_error, "unsupported component rank encountered");
        }
      }
      return value;
    }
    
    /**
     \brief Accessor for rank-2 objects that accepts fixed-length arrays with entries corresponding to component indices.
     \param [in] entryComponents0 - an array with one entry per tensorial component, each entry indicating the requested index into that component's first dimension
     \param [in] entryComponents1 - an array with one entry per tensorial component, each entry indicating the requested index into that component's second dimension
     
     \return The product of tensor component values with the specified component entry indices.
    */
    template <typename iType0, typename iType1, ordinal_type numTensorComponents>
    KOKKOS_INLINE_FUNCTION typename std::enable_if<
        (std::is_integral<iType0>::value && std::is_integral<iType1>::value),
        Scalar>::type
    operator()(const Kokkos::Array<iType0,numTensorComponents>& entryComponents0, const Kokkos::Array<iType1,numTensorComponents>& entryComponents1) const {
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(numTensorComponents_ != numTensorComponents, std::invalid_argument, "Tensorial component count mismatch");
#endif
      Scalar value = 1.0;
      for (ordinal_type r=0; r<numTensorComponents; r++)
      {
        auto & component = tensorComponents_[r];
        const ordinal_type componentRank = component.rank();
        if (componentRank == 2)
        {
          value *= component(entryComponents0[r],entryComponents1[r]);
        }
        else if (componentRank == 1)
        {
          value *= component(entryComponents0[r]);
        }
        else
        {
          INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::logic_error, "unsupported component rank encountered");
        }
      }
      return value;
    }
    
  /** \brief  Returns the logical extent in the requested dimension
       \param [in] d - the dimension
       
     \return logical extent as an integer
    */
    template <typename iType>
    KOKKOS_INLINE_FUNCTION
    typename std::enable_if<std::is_integral<iType>::value, ordinal_type>::type
    extent_int(const iType& d) const {
      return extents_[d];
    }
    
    /** \brief  Returns the logical extent in the requested dimension
         \param [in] d - the dimension
         
       \return the logical extent in the requested dimension.
      */
    template <typename iType>
    KOKKOS_INLINE_FUNCTION constexpr
    typename std::enable_if<std::is_integral<iType>::value, size_t>::type
    extent(const iType& d) const {
      return extents_[d];
    }
    
    //! Returns true for containers that have data; false for those that don't (e.g., those that have been constructed by the default constructor).
    KOKKOS_INLINE_FUNCTION constexpr bool isValid() const
    {
      return extents_[0] > 0;
    }
    
    //! Returns the rank of the container.
    KOKKOS_INLINE_FUNCTION
    ordinal_type rank() const
    {
      return rank_;
    }
    
    //! Return the number of tensorial components.
    KOKKOS_INLINE_FUNCTION
    ordinal_type numTensorComponents() const
    {
      return numTensorComponents_;
    }
    
    //! Returns true if the first component is indexed separately; false if not.
    KOKKOS_INLINE_FUNCTION
    bool separateFirstComponent() const
    {
      return separateFirstComponent_;
    }
    
    //! Sets the extent of the first component.  Only valid when either there is only one component, or when separateFirstComponent() returns true.  The intended use case is when the 0 dimension in first component represents a cell index, and the container is resized to match a workset size that does not evenly divide the number of cells.
    void setFirstComponentExtentInDimension0(const ordinal_type &newExtent)
    {
      INTREPID2_TEST_FOR_EXCEPTION(!separateFirstComponent_ && (numTensorComponents_ != 1), std::invalid_argument, "setFirstComponentExtent() is only allowed when separateFirstComponent_ is true, or there is only one component");
      tensorComponents_[0].setExtent(0,newExtent);
      extents_[0] = newExtent;
    }
  };
}

#endif /* Intrepid2_TensorData_h */
