// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_TransformedBasisValues.hpp
    \brief  Structure-preserving representation of transformed basis values; reference space values and transformations are stored separately.
 
 There are effectively two modes: one for vector-valued BasisValues, and one for scalar-valued BasisValues.  In the former case the transformation is logically a matrix, with shape (C,P,D,D).  In the latter case, the transformation is logically a weight on each physical-space quadrature point, with shape (C,P).  If the transform is left unset, it is understood to be the identity.

    \author Nathan V. Roberts
*/

#ifndef Intrepid2_TransformedBasisValues_h
#define Intrepid2_TransformedBasisValues_h

#include "Intrepid2_BasisValues.hpp"
#include "Intrepid2_DataTools.hpp"
#include "Intrepid2_ScalarView.hpp"
#include "Intrepid2_Utils.hpp"

namespace Intrepid2 {
/** \class Intrepid2::TransformedBasisValues
    \brief Structure-preserving representation of transformed vector data; reference space values and transformations are stored separately.
 
 TransformedBasisValues provides a View-like interface of rank 3 or 4, with shape (C,F,P) or (C,F,P,D).  When the corresponding accessor is used, the transformed value is determined from corresponding reference space values and the transformation.
*/
  template<class Scalar, typename DeviceType>
  class TransformedBasisValues
  {
  public:
    ordinal_type numCells_;
    
    Data<Scalar,DeviceType> transform_; // vector case: (C,P,D,D) jacobian or jacobian inverse; can also be unset for identity transform.  Scalar case: (C,P), or unset for identity.  Contracted vector case: (C,P,D) transform, to be contracted with a vector field to produce a scalar result.
    
    BasisValues<Scalar, DeviceType> basisValues_;
    
    /**
     \brief Standard constructor.
     \param [in] transform - the transformation (matrix), with logical shape (C,P), (C,P,D), or (C,P,D,D)
     \param [in] basisValues - the reference-space data to be transformed, with logical shape (F,P) (for scalar values) or (F,P,D) (for vector values)
    */
    TransformedBasisValues(const Data<Scalar,DeviceType> &transform, const BasisValues<Scalar,DeviceType> &basisValues)
    :
    numCells_(transform.extent_int(0)),
    transform_(transform),
    basisValues_(basisValues)
    {
      // sanity check: when transform is diagonal, we expect there to be no pointwise variation.
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(transform_.isDiagonal() && (transform_.getVariationTypes()[1] != CONSTANT), std::invalid_argument, "When transform is diagonal, we assume in various places that there is no pointwise variation; the transform_ Data should have CONSTANT as its variation type in dimension 1.");
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE((transform_.rank() < 2) || (transform_.rank() > 4), std::invalid_argument, "Only transforms of rank 2, 3, or 4 are supported");
    }
    
    /**
     \brief Constructor for the case of an identity transform.
     \param [in] vectorData - the reference-space data, with logical shape (F,P,D)
    */
    TransformedBasisValues(const ordinal_type &numCells, const BasisValues<Scalar,DeviceType> &basisValues)
    :
    numCells_(numCells),
    basisValues_(basisValues)
    {}
    
    //! copy-like constructor for differing device types.  This may do a deep_copy of underlying views, depending on the memory spaces involved.
    template<typename OtherDeviceType, class = typename std::enable_if<!std::is_same<DeviceType, OtherDeviceType>::value>::type>
    TransformedBasisValues(const TransformedBasisValues<Scalar,OtherDeviceType> &transformedVectorData)
    :
    numCells_(transformedVectorData.numCells()),
    transform_(transformedVectorData.transform()),
    basisValues_(transformedVectorData.basisValues())
    {}
    
    /**
     \brief Default constructor; an invalid container.  Will return -1 for numCells().
     */
    TransformedBasisValues()
    :
    numCells_(-1)
    {}
    
    //! Returns true if the transformation matrix is diagonal.
    KOKKOS_INLINE_FUNCTION bool axisAligned() const
    {
      if (!transform_.isValid())
      {
        // null transform is understood as the identity
        return true;
      }
      else
      {
        return transform_.isDiagonal();
      }
    }

    BasisValues<Scalar, DeviceType> basisValues() const
    {
      return basisValues_;
    }

    //! Returns the true data extent in the cell dimension (e.g., will be 1 for transform matrices that do not vary from one cell to the next).
    KOKKOS_INLINE_FUNCTION int cellDataExtent() const
    {
      return transform_.getDataExtent(0);
    }
    
    //! Returns the variation type corresponding to the cell dimension.
    KOKKOS_INLINE_FUNCTION DataVariationType cellVariationType() const
    {
      return transform_.getVariationTypes()[0];
    }
    
    //! Replaces the internal pullback (transformation operator) with the result of the pullback multiplied by the specified (C,P) weights.  ViewType may be a rank-2 Kokkos::View, a rank-2 Kokkos::DynRankView, or a rank-2 Intrepid2::Data object.
    template<class ViewType>
    void multiplyByPointwiseWeights(const ViewType &weights)
    {
      ordinal_type weightRank = getFunctorRank(weights); // .rank() or ::rank, depending on weights type
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(weightRank != 2, std::invalid_argument, "pointwise weights must have shape (C,P).");
      
      Data<Scalar,DeviceType> weightData(weights);
      if (!transform_.isValid())
      {
        // empty transform_ is understood as identity; multiplying by weightData is thus
        // the same as transform_ = weightData
        transform_ = weightData;
        return;
      }
      else
      {
        if ((transform_.rank() == 4) || (transform_.rank() == 3))
        {
          transform_ = DataTools::multiplyByCPWeights(transform_,weightData);
        }
        else // transformRank == 2
        {
          auto result = Data<Scalar,DeviceType>::allocateInPlaceCombinationResult(weightData, transform_);
          
          result.storeInPlaceProduct(weightData,transform_);
          transform_ = result;
        }
      }
    }
    
    //! Returns the logical extent in the cell dimension, which is the 0 dimension in this container.
    KOKKOS_INLINE_FUNCTION int numCells() const
    {
      return numCells_;
    }
    
    //! Returns the logical extent in the fields dimension, which is the 1 dimension in this container.
    KOKKOS_INLINE_FUNCTION int numFields() const
    {
      return basisValues_.extent_int(0);
    }
    
    //! Returns the logical extent in the points dimension, which is the 2 dimension in this container.
    KOKKOS_INLINE_FUNCTION int numPoints() const
    {
      return basisValues_.extent_int(1);
    }
    
    //! Returns the logical extent in the space dimension, which is the 3 dimension in this container.
    KOKKOS_INLINE_FUNCTION int spaceDim() const
    {
      if ((transform_.rank() == 3) && (basisValues_.rank() == 3)) // (C,P,D) contracted in D against (F,P,D)
      {
        return 1; // spaceDim contracted away
      }
      else if ((transform_.rank() == 3) && (basisValues_.rank() == 2)) // (C,P,D) weighting (F,P)
      {
        return transform_.extent_int(2);
      }
      else if (transform_.isValid())
      {
        return transform_.extent_int(2);
      }
      else
      {
        return basisValues_.extent_int(2);
      }
    }
    
    //! Scalar accessor, with arguments (C,F,P).
    KOKKOS_INLINE_FUNCTION Scalar operator()(const int &cellOrdinal, const int &fieldOrdinal, const int &pointOrdinal) const
    {
      if (!transform_.isValid())
      {
        // null transform is understood as the identity
        return basisValues_(fieldOrdinal,pointOrdinal);
      }
      else if (transform_.rank() == 2)
      {
        return transform_(cellOrdinal,pointOrdinal) * basisValues_(fieldOrdinal,pointOrdinal);
      }
      else if (transform_.rank() == 3)
      {
        Scalar value = 0;
        for (int d=0; d<transform_.extent_int(2); d++)
        {
          value += transform_(cellOrdinal,pointOrdinal,d) * basisValues_(fieldOrdinal,pointOrdinal,d);
        }
        return value;
      }
      return 0; // should not be reachable
    }
    
    //! Vector accessor, with arguments (C,F,P,D).
    KOKKOS_INLINE_FUNCTION Scalar operator()(const int &cellOrdinal, const int &fieldOrdinal, const int &pointOrdinal, const int &dim) const
    {
      if (!transform_.isValid())
      {
        // null transform is understood as the identity
        return basisValues_(fieldOrdinal,pointOrdinal,dim);
      }
      else if (transform_.isDiagonal())
      {
        return transform_(cellOrdinal,pointOrdinal,dim,dim) * basisValues_(fieldOrdinal,pointOrdinal,dim);
      }
      else if (transform_.rank() == 4)
      {
        Scalar value = 0.0;
        for (int d2=0; d2<transform_.extent_int(2); d2++)
        {
          value += transform_(cellOrdinal,pointOrdinal,dim,d2) * basisValues_(fieldOrdinal,pointOrdinal,d2);
        }
        return value;
      }
      else if (transform_.rank() == 3)
      {
        Scalar value = transform_(cellOrdinal,pointOrdinal,dim) * basisValues_(fieldOrdinal,pointOrdinal);
        return value;
      }
      else // rank 2 transform
      {
        Scalar value = transform_(cellOrdinal,pointOrdinal) * basisValues_(fieldOrdinal,pointOrdinal,dim);
        return value;
      }
    }
    
    //! Returns the specified entry in the (scalar) transform.  (Only valid for scalar-valued BasisValues; see the four-argument transformWeight() for the vector-valued case.)
    KOKKOS_INLINE_FUNCTION Scalar transformWeight(const int &cellOrdinal, const int &pointOrdinal) const
    {
      if (!transform_.isValid())
      {
        // null transform is understood as identity
        return 1.0;
      }
      else
      {
        return transform_(cellOrdinal,pointOrdinal);
      }
    }
    
    //! Returns the specified entry in the transformation vector.
    KOKKOS_INLINE_FUNCTION Scalar transformWeight(const int &cellOrdinal, const int &pointOrdinal, const int &d) const
    {
      if (!transform_.isValid())
      {
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "three-argument transformWeight() is not supported for invalid transform_ object -- no meaningful interpretation for vector-valued identity");
      }
      else
      {
        return transform_(cellOrdinal,pointOrdinal,d);
      }
    }
    
    //! Returns the specified entry in the transform matrix.
    KOKKOS_INLINE_FUNCTION Scalar transformWeight(const int &cellOrdinal, const int &pointOrdinal, const int &dim1, const int &dim2) const
    {
      if (!transform_.isValid())
      {
        // null transform is understood as identity
        return (dim1 == dim2) ? 1.0 : 0.0;
      }
      else
      {
        return transform_(cellOrdinal,pointOrdinal,dim1,dim2);
      }
    }
    
    //! Returns the transform matrix.  An invalid/empty container indicates the identity transform.
    const Data<Scalar,DeviceType> & transform() const
    {
      return transform_;
    }
    
    //! Returns the reference-space vector data.
    const VectorData<Scalar,DeviceType> & vectorData() const
    {
      return basisValues_.vectorData();
    }
    
    //! Returns the rank of the container, which is 3 for scalar values, and 4 for vector values.
    KOKKOS_INLINE_FUNCTION
    unsigned rank() const
    {
      if ((transform_.rank() == 4) && (basisValues_.rank() == 3))
      {
        return 4; // (C,F,P,D)
      }
      else if (transform_.rank() == 2)
      {
        return basisValues_.rank() + 1; // transformation adds a cell dimension
      }
      else if (transform_.rank() == 3)
      {
        if (basisValues_.rank() == 3)
        {
          // transform contracts with basisValues in D dimension
          return 3; // (C,F,P)
        }
        else if (basisValues_.rank() == 2) // (F,P)
        {
          return 4; // (C,F,P,D)
        }
      }
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "Unhandled basisValues_/transform_ rank combination");
    }
    
    //! Returns the extent in the specified dimension as an int.
    KOKKOS_INLINE_FUNCTION
    int extent_int(const int &r) const
    {
      if      (r == 0) return numCells();
      else if (r == 1) return numFields();
      else if (r == 2) return numPoints();
      else if (r == 3) return spaceDim();
      else if (r  > 3) return 1;
      
      return -1; // unreachable return; here to avoid compiler warnings.
    }
  };
}

#endif /* Intrepid2_TransformedBasisValues_h */
