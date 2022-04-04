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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov),
//                    Mauro Perego  (mperego@sandia.gov), or
//                    Nate Roberts  (nvrober@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid2_TransformedBasisValues.hpp
    \brief  Structure-preserving representation of transformed basis values; reference space values and transformations are stored separately.
 
 There are effectively two modes: one for vector-valued BasisValues, and one for scalar-valued BasisValues.  In the former case the transformation is logically a matrix, with shape (C,P,D,D).  In the latter case, the transformation is logically a weight on each physical-space quadrature point, with shape (C,P).  If the transform is left unset, it is understood to be the identity.

    \author Nathan V. Roberts
*/

#ifndef Intrepid2_TransformedBasisValues_h
#define Intrepid2_TransformedBasisValues_h

#include "Intrepid2_BasisValues.hpp"
#include "Intrepid2_ScalarView.hpp"

namespace Intrepid2 {
/** \class Intrepid2::TransformedBasisValues
    \brief Structure-preserving representation of transformed vector data; reference space values and transformations are stored separately.
 
 TransformedBasisValues provides a View-like interface of rank 4, with shape (C,F,P,D).  When the corresponding accessor is used, the transformed value is determined from corresponding reference space values and the transformation.
*/
  template<class Scalar, typename DeviceType>
  class TransformedBasisValues
  {
  public:
    ordinal_type numCells_;
    
    Data<Scalar,DeviceType> transform_; // vector case: (C,P,D,D) jacobian or jacobian inverse; can also be unset for identity transform.  Scalar case: (C,P), or unset for identity.
    
    BasisValues<Scalar, DeviceType> basisValues_;
    
    /**
     \brief Standard constructor.
     \param [in] transform - the transformation (matrix), with logical shape (C,P) or (C,P,D,D)
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
      return basisValues_.extent_int(2);
    }
    
    //! Scalar accessor, with arguments (C,F,P).
    KOKKOS_INLINE_FUNCTION Scalar operator()(const int &cellOrdinal, const int &fieldOrdinal, const int &pointOrdinal) const
    {
      if (!transform_.isValid())
      {
        // null transform is understood as the identity
        return basisValues_(fieldOrdinal,pointOrdinal);
      }
      else
      {
        return transform_(cellOrdinal,pointOrdinal) * basisValues_(fieldOrdinal,pointOrdinal);
      }
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
      else
      {
        Scalar value = 0.0;
        for (int d2=0; d2<transform_.extent_int(2); d2++)
        {
          value += transform_(cellOrdinal,pointOrdinal,dim,d2) * basisValues_(fieldOrdinal,pointOrdinal,d2);
        }
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
      return basisValues_.rank() + 1; // transformation adds a cell dimension
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
