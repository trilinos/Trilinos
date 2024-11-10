// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_VectorData.hpp
    \brief  Reference-space field values for a basis, designed to support typical vector-valued bases.

    \author Nathan V. Roberts
*/

#ifndef Intrepid2_VectorData_h
#define Intrepid2_VectorData_h

namespace Intrepid2 {
/** \class Intrepid2::VectorData
    \brief Reference-space field values for a basis, designed to support typical vector-valued bases.
 
  VectorData is designed with typical HDIV/HCURL bases in mind; these often involve reference-space basis functions
  each of which is zero in every dimension but one.  Moreover, on tensor product topologies each nonzero scalar function can be
  expressed as a product of functions defined on the topological components.  Also supported: gradients of typical H^1 bases,
  which have a full vector, each component of which is a product of basis functions defined on lower-dimensional topologies.
 
 Typically, HDIV and HCURL bases have multiple families corresponding to the nonzero structure of the vector components.  VectorData therefore supports multiple families.  (Gradients of H^1 are expressed as a single family.)
*/
  template<class Scalar, typename DeviceType>
  class VectorData
  {
  public:
    using VectorArray = Kokkos::Array< TensorData<Scalar,DeviceType>, Parameters::MaxVectorComponents >; // for axis-aligned case, these correspond entry-wise to the axis with which the vector values align
    using FamilyVectorArray = Kokkos::Array< VectorArray, Parameters::MaxTensorComponents>;

    FamilyVectorArray vectorComponents_; // outer: family ordinal; inner: component/spatial dimension ordinal
    bool axialComponents_; // if true, each entry in vectorComponents_ is an axial component vector; for 3D: (f1,0,0); (0,f2,0); (0,0,f3).  The 0s are represented by trivial/invalid TensorData objects.  In this case, numComponents_ == numFamilies_.
     
    int totalDimension_;
    Kokkos::Array<int, Parameters::MaxVectorComponents> dimToComponent_;
    Kokkos::Array<int, Parameters::MaxVectorComponents> dimToComponentDim_;
    Kokkos::Array<int, Parameters::MaxVectorComponents> numDimsForComponent_;
    
    Kokkos::Array<int,Parameters::MaxTensorComponents> familyFieldUpperBound_; // one higher than the end of family indicated
    
    unsigned numFamilies_;   // number of valid entries in vectorComponents_
    unsigned numComponents_; // number of valid entries in each entry of vectorComponents_
    unsigned numPoints_;     // the second dimension of each (valid) TensorData entry
    
    /**
     \brief Initialize members based on constructor parameters; all constructors should call this after populating numFamilies_, numComponents_, and vectorComponents_.
    */
    void initialize()
    {
      int lastFieldUpperBound = 0;
      int numPoints = 0;
      axialComponents_ = true;  // will set to false if we find any valid entries that are not on the "diagonal" (like position for family/component)
      for (unsigned i=0; i<numFamilies_; i++)
      {
        bool validEntryFoundForFamily = false;
        int numFieldsInFamily = 0;
        for (unsigned j=0; j<numComponents_; j++)
        {
          if (vectorComponents_[i][j].isValid())
          {
            if (!validEntryFoundForFamily)
            {
              numFieldsInFamily = vectorComponents_[i][j].extent_int(0); // (F,P[,D])
              validEntryFoundForFamily = true;
            }
            else
            {
              INTREPID2_TEST_FOR_EXCEPTION(numFieldsInFamily != vectorComponents_[i][j].extent_int(0), std::invalid_argument, "Each valid TensorData entry within a family must agree with the others on the number of fields in the family");
            }
            if (numPoints == 0)
            {
              numPoints = vectorComponents_[i][j].extent_int(1); // (F,P[,D])
            }
            else
            {
              INTREPID2_TEST_FOR_EXCEPTION(numPoints != vectorComponents_[i][j].extent_int(1), std::invalid_argument, "Each valid TensorData entry must agree with the others on the number of points");
            }
            if (i != j)
            {
              // valid entry found that is not on the "diagonal": axialComponents is false
              axialComponents_ = false;
            }
          }
        }
        lastFieldUpperBound += numFieldsInFamily;
        familyFieldUpperBound_[i] = lastFieldUpperBound;
        INTREPID2_TEST_FOR_EXCEPTION(!validEntryFoundForFamily, std::invalid_argument, "Each family must have at least one valid TensorData entry");
      }

      // do a pass through components to determine total component dim (totalDimension_) and size lookups appropriately
      int currentDim = 0;
      for (unsigned j=0; j<numComponents_; j++)
      {
        bool validEntryFoundForComponent = false;
        int numDimsForComponent = 0;
        for (unsigned i=0; i<numFamilies_; i++)
        {
          if (vectorComponents_[i][j].isValid())
          {
            if (!validEntryFoundForComponent)
            {
              validEntryFoundForComponent = true;
              numDimsForComponent = vectorComponents_[i][j].extent_int(2); // (F,P,D) container or (F,P) container
            }
            else
            {
              INTREPID2_TEST_FOR_EXCEPTION(numDimsForComponent != vectorComponents_[i][j].extent_int(2), std::invalid_argument, "Components in like positions must agree across families on the number of dimensions spanned by that component position");
            }
          }
        }
        if (!validEntryFoundForComponent)
        {
          // assume that the component takes up exactly one space dim
          numDimsForComponent = 1;
        }
        
        numDimsForComponent_[j] = numDimsForComponent;
        
        currentDim += numDimsForComponent;
      }
      totalDimension_ = currentDim;
      
      currentDim = 0;
      for (unsigned j=0; j<numComponents_; j++)
      {
        int numDimsForComponent = numDimsForComponent_[j];
        for (int dim=0; dim<numDimsForComponent; dim++)
        {
          dimToComponent_[currentDim+dim]    = j;
          dimToComponentDim_[currentDim+dim] = dim;
        }
        currentDim += numDimsForComponent;
      }
      numPoints_      = numPoints;
    }
  public:
    /**
     \brief Standard constructor for the arbitrary case, accepting a fixed-length array argument.
     \param [in] vectorComponents - an array of arrays, where the outer dimension corresponds to the family, and the inner to the number of components in each vector.
     
     Outer dimension: number of families; inner dimension: number of components in each vector.  Use empty/invalid TensorData objects to indicate zeroes.  Each family, and each vector component dimension, must have at least one valid entry, and the number of points in each valid entry must agree with each other.  The field count within components of a family must also agree; across families these may differ.  Vector components are allowed to span multiple spatial dimensions, but when they do, any valid TensorData object in that vector position must agree in the dimension across families.
    */
    template<size_t numFamilies, size_t numComponents>
    VectorData(Kokkos::Array< Kokkos::Array<TensorData<Scalar,DeviceType>, numComponents>, numFamilies> vectorComponents)
    :
    numFamilies_(numFamilies),
    numComponents_(numComponents)
    {
      static_assert(numFamilies <= Parameters::MaxTensorComponents,   "numFamilies must be less than Parameters::MaxTensorComponents");
      static_assert(numComponents <= Parameters::MaxVectorComponents, "numComponents must be less than Parameters::MaxVectorComponents");
      for (unsigned i=0; i<numFamilies; i++)
      {
        for (unsigned j=0; j<numComponents; j++)
        {
          vectorComponents_[i][j] = vectorComponents[i][j];
        }
      }
      initialize();
    }
    
    /**
     \brief Standard constructor for the arbitrary case, accepting a variable-length std::vector argument.
     \param [in] vectorComponents - an array of arrays, where the outer dimension corresponds to the family, and the inner to the number of components in each vector.
     
     Outer dimension: number of families; inner dimension: number of components in each vector.  Use empty/invalid TensorData objects to indicate zeroes.  Each family, and each vector component dimension, must have at least one valid entry, and the number of points in each valid entry must agree with each other.  The field count within components of a family must also agree; across families these may differ.  Vector components are allowed to span multiple spatial dimensions, but when they do, any valid TensorData object in that vector position must agree in the dimension across families.
    */
    VectorData(const std::vector< std::vector<TensorData<Scalar,DeviceType> > > &vectorComponents)
    {
      numFamilies_ = vectorComponents.size();
      INTREPID2_TEST_FOR_EXCEPTION(numFamilies_ <= 0, std::invalid_argument, "numFamilies must be at least 1");
      numComponents_ = vectorComponents[0].size();
      for (unsigned i=1; i<numFamilies_; i++)
      {
        INTREPID2_TEST_FOR_EXCEPTION(vectorComponents[i].size() != numComponents_, std::invalid_argument, "each family must have the same number of components");
      }
      
      INTREPID2_TEST_FOR_EXCEPTION(numFamilies_ > Parameters::MaxTensorComponents,   std::invalid_argument, "numFamilies must be at most Parameters::MaxTensorComponents");
      INTREPID2_TEST_FOR_EXCEPTION(numComponents_ > Parameters::MaxVectorComponents,   std::invalid_argument, "numComponents must be at most Parameters::MaxVectorComponents");
      for (unsigned i=0; i<numFamilies_; i++)
      {
        for (unsigned j=0; j<numComponents_; j++)
        {
          vectorComponents_[i][j] = vectorComponents[i][j];
        }
      }
      initialize();
    }
    
    /**
     \brief Simplified constructor for gradients of HGRAD, and values of HDIV and HCURL vector bases.
     \param [in] vectorComponents - an array of components; see below for interpretation
     \param [in] axialComponents - use false for a "full" vector, such as the gradient of a scalar HGRAD basis; use true for all-but-one zero entries in each vector, as for the values in a typical HDIV or HCURL basis
     
     For the gradient use case, VectorData will have a single family.  For the HDIV and HCURL use cases, will have the same number of families as there are components in the vectorComponents argument; each family will consist of vectors that have one entry filled, the others zeroes; this is what we mean by "axial components."
    */
    template<size_t numComponents>
    VectorData(Kokkos::Array< TensorData<Scalar,DeviceType>, numComponents> vectorComponents, bool axialComponents)
    {
      if (axialComponents)
      {
        numFamilies_   = numComponents;
        numComponents_ = numComponents;
        for (unsigned d=0; d<numComponents_; d++)
        {
          vectorComponents_[d][d] = vectorComponents[d];
        }
      }
      else
      {
        numFamilies_   = 1;
        numComponents_ = numComponents;
        for (unsigned d=0; d<numComponents_; d++)
        {
          vectorComponents_[0][d] = vectorComponents[d];
        }
      }
      initialize();
    }
    
    /**
     \brief Simplified constructor for gradients of HGRAD, and values of HDIV and HCURL vector bases.
     \param [in] vectorComponents - an array of components; see below for interpretation
     \param [in] axialComponents - use false for a "full" vector, such as the gradient of a scalar HGRAD basis; use true for all-but-one zero entries in each vector, as for the values in a typical HDIV or HCURL basis
     
     For the gradient use case, VectorData will have a single family.  For the HDIV and HCURL use cases, will have the same number of families as there are components in the vectorComponents argument; each family will consist of vectors that have one entry filled, the others zeroes; this is what we mean by "axial components."
    */
    VectorData(std::vector< TensorData<Scalar,DeviceType> > vectorComponents, bool axialComponents)
    : numComponents_(vectorComponents.size())
    {
      if (axialComponents)
      {
        numFamilies_   = numComponents_;
        for (unsigned d=0; d<numComponents_; d++)
        {
          vectorComponents_[d][d] = vectorComponents[d];
        }
      }
      else
      {
        numFamilies_   = 1;
        for (unsigned d=0; d<numComponents_; d++)
        {
          vectorComponents_[0][d] = vectorComponents[d];
        }
      }
      initialize();
    }
    
    //! copy-like constructor for differing device type, but same memory space.  This does a shallow copy of the underlying view.
    template<typename OtherDeviceType, class = typename std::enable_if< std::is_same<typename DeviceType::memory_space, typename OtherDeviceType::memory_space>::value>::type,
                                       class = typename std::enable_if<!std::is_same<DeviceType,OtherDeviceType>::value>::type>
    VectorData(const VectorData<Scalar,OtherDeviceType> &vectorData)
    :
    numFamilies_(vectorData.numFamilies()),
    numComponents_(vectorData.numComponents())
    {
      if (vectorData.isValid())
      {
        for (unsigned i=0; i<numFamilies_; i++)
        {
          for (unsigned j=0; j<numComponents_; j++)
          {
            vectorComponents_[i][j] = vectorData.getComponent(i, j);
          }
        }
        initialize();
      }
    }
    
    //! copy-like constructor for differing execution spaces.  This does a deep copy of underlying views.
    template<typename OtherDeviceType, class = typename std::enable_if<!std::is_same<typename DeviceType::memory_space, typename OtherDeviceType::memory_space>::value>::type>
    VectorData(const VectorData<Scalar,OtherDeviceType> &vectorData)
    :
    numFamilies_(vectorData.numFamilies()),
    numComponents_(vectorData.numComponents())
    {
      if (vectorData.isValid())
      {
        for (unsigned i=0; i<numFamilies_; i++)
        {
          for (unsigned j=0; j<numComponents_; j++)
          {
            vectorComponents_[i][j] = vectorData.getComponent(i, j);
          }
        }
        initialize();
      }
    }
    
    //! Simple 1-argument constructor for the case of trivial tensor product structure.  The argument should have shape (F,P,D) where D has extent equal to the spatial dimension.
    VectorData(TensorData<Scalar,DeviceType> data)
    :
    VectorData(Kokkos::Array< TensorData<Scalar,DeviceType>, 1>(data), true)
    {}
    
    //! Simple 1-argument constructor for the case of trivial tensor product structure.  The argument should have shape (F,P,D) where D has extent equal to the spatial dimension.
    VectorData(Data<Scalar,DeviceType> data)
    :
    VectorData(Kokkos::Array< TensorData<Scalar,DeviceType>, 1>({TensorData<Scalar,DeviceType>(data)}), true)
    {}
    
    //! default constructor; results in an invalid container.
    VectorData()
    :
    numFamilies_(0), numComponents_(0)
    {}
    
    //! Returns true only if the families are so structured that the first family has nonzeros only in the x component, the second only in the y component, etc.
    KOKKOS_INLINE_FUNCTION
    bool axialComponents() const
    {
      return axialComponents_;
    }
    
    //! Returns the number of dimensions corresponding to the specified component.
    KOKKOS_INLINE_FUNCTION
    int numDimsForComponent(int &componentOrdinal) const
    {
      return numDimsForComponent_[componentOrdinal];
    }
    
    //! Returns the total number of fields; corresponds to the first dimension of this container.
    KOKKOS_INLINE_FUNCTION
    int numFields() const
    {
      return familyFieldUpperBound_[numFamilies_-1];
    }
    
    //! Returns the field ordinal offset for the specified family.
    KOKKOS_INLINE_FUNCTION
    int familyFieldOrdinalOffset(const int &familyOrdinal) const
    {
      return (familyOrdinal == 0) ? 0 : familyFieldUpperBound_[familyOrdinal-1];
    }
    
    //! Returns the number of points; corresponds to the second dimension of this container.
    KOKKOS_INLINE_FUNCTION
    int numPoints() const
    {
      return numPoints_;
    }
    
    //! Returns the spatial dimension; corresponds to the third dimension of this container.
    KOKKOS_INLINE_FUNCTION
    int spaceDim() const
    {
      return totalDimension_;
    }
    
    //! Accessor for the container, which has shape (F,P,D).
    KOKKOS_INLINE_FUNCTION
    Scalar operator()(const int &fieldOrdinal, const int &pointOrdinal, const int &dim) const
    {
      int fieldOrdinalInFamily = fieldOrdinal;
      int familyForField = 0;
      if (numFamilies_ > 1)
      {
        familyForField = -1;
        int previousFamilyEnd = -1;
        int fieldAdjustment = 0;
        // this loop is written in such a way as to avoid branching for CUDA performance
        for (unsigned family=0; family<numFamilies_; family++)
        {
          const bool fieldInRange = (fieldOrdinal > previousFamilyEnd) && (fieldOrdinal < familyFieldUpperBound_[family]);
          familyForField = fieldInRange ? family : familyForField;
          fieldAdjustment = fieldInRange ? previousFamilyEnd + 1 : fieldAdjustment;
          previousFamilyEnd = familyFieldUpperBound_[family] - 1;
        }
#ifdef HAVE_INTREPID2_DEBUG
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(familyForField == -1, std::invalid_argument, "family for field not found");
#endif
        
        fieldOrdinalInFamily = fieldOrdinal - fieldAdjustment;
      }
      
      const int componentOrdinal = dimToComponent_[dim];
      
      const auto &component = vectorComponents_[familyForField][componentOrdinal];
      if (component.isValid())
      {
        const int componentRank     = component.rank();
        if (componentRank == 2) // (F,P) container
        {
          return component(fieldOrdinalInFamily,pointOrdinal);
        }
        else if (componentRank == 3) // (F,P,D)
        {
          return component(fieldOrdinalInFamily,pointOrdinal,dimToComponentDim_[dim]);
        }
        else
        {
          INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::logic_error, "Unsupported component rank");
          return -1; // unreachable, but compilers complain otherwise...
        }
      }
      else // invalid component: placeholder means 0
      {
        return 0;
      }
    }
    
    /**
     \brief Single-argument component accessor for the axial-component or the single-family case; in this case, one argument suffices to uniquely identify the component.
     \param [in] componentOrdinal - the vector component ordinal.
    */
    KOKKOS_INLINE_FUNCTION
    const TensorData<Scalar,DeviceType> & getComponent(const int &componentOrdinal) const
    {
      if (axialComponents_)
      {
        return vectorComponents_[componentOrdinal][componentOrdinal];
      }
      else if (numFamilies_ == 1)
      {
        return vectorComponents_[0][componentOrdinal];
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "Ambiguous component request; use the two-argument getComponent()");
      }
      // nvcc warns here about a missing return.
      return vectorComponents_[6][6]; // likely this is an empty container, but anyway it's an unreachable line...
    }
    
    /**
     \brief General component accessor.
     \param [in] familyOrdinal - the family ordinal for the requested component.
     \param [in] componentOrdinal - the vector component ordinal.
    */
    KOKKOS_INLINE_FUNCTION
    const TensorData<Scalar,DeviceType> & getComponent(const int &familyOrdinal, const int &componentOrdinal) const
    {
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(familyOrdinal < 0, std::invalid_argument, "familyOrdinal must be non-negative");
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(static_cast<unsigned>(familyOrdinal) >= numFamilies_, std::invalid_argument, "familyOrdinal out of bounds");
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(componentOrdinal < 0, std::invalid_argument, "componentOrdinal must be non-negative");
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(static_cast<unsigned>(componentOrdinal) >= numComponents_, std::invalid_argument, "componentOrdinal out of bounds");
      
      return vectorComponents_[familyOrdinal][componentOrdinal];
    }
    
    //! Returns the extent in the specified dimension as an int.
    KOKKOS_INLINE_FUNCTION
    int extent_int(const int &r) const
    {
      // logically (F,P,D) container
      if      (r == 0) return numFields();
      else if (r == 1) return numPoints();
      else if (r == 2) return totalDimension_;
      else if (r  > 2) return 1;
      
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "Unsupported rank");
      return -1; // unreachable; return here included to avoid compiler warnings.
    }
    
    //! Returns the rank of this container, which is 3.
    KOKKOS_INLINE_FUNCTION
    unsigned rank() const
    {
      // logically (F,P,D) container
      return 3;
    }
    
    //! returns the number of components
    KOKKOS_INLINE_FUNCTION int numComponents() const
    {
      return numComponents_;
    }
    
    //! returns the number of families
    KOKKOS_INLINE_FUNCTION int numFamilies() const
    {
      return numFamilies_;
    }
    
    //! Returns the family ordinal corresponding to the indicated field ordinal.
    KOKKOS_INLINE_FUNCTION int familyForFieldOrdinal(const int &fieldOrdinal) const
    {
      int matchingFamily = -1;
      int fieldsSoFar = 0;
      // logic here is a little bit more complex to avoid branch divergence
      for (int i=0; i<numFamilies_; i++)
      {
        const bool fieldIsBeyondPreviousFamily = (fieldOrdinal >= fieldsSoFar);
        fieldsSoFar += numFieldsInFamily(i);
        const bool fieldIsBeforeCurrentFamily  = (fieldOrdinal < fieldsSoFar);
        const bool fieldMatchesFamily = fieldIsBeyondPreviousFamily && fieldIsBeforeCurrentFamily;
        matchingFamily = fieldMatchesFamily ? i : matchingFamily;
      }
      return matchingFamily;
    }
    
    //! returns the number of fields in the specified family
    KOKKOS_INLINE_FUNCTION int numFieldsInFamily(const unsigned &familyOrdinal) const
    {
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(familyOrdinal >= numFamilies_, std::invalid_argument, "familyOrdinal out of bounds");
      int numFields = -1;
      for (unsigned componentOrdinal=0; componentOrdinal<numComponents_; componentOrdinal++)
      {
        numFields = vectorComponents_[familyOrdinal][componentOrdinal].isValid() ? vectorComponents_[familyOrdinal][componentOrdinal].extent_int(0) : numFields;
      }
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(numFields < 0, std::logic_error, "numFields was not properly initialized");
      return numFields;
    }
    
    //! returns true for containers that have data; false for those that don't (e.g., those that have been constructed by the default constructor).
    KOKKOS_INLINE_FUNCTION constexpr bool isValid() const
    {
      return numComponents_ > 0;
    }
  };
}

#endif /* Intrepid2_VectorData_h */
