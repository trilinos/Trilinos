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

/** \file   Intrepid2_BasisValues.hpp
    \brief  Header file for the data-wrapper class Intrepid2::BasisValues.
    \author Nathan V. Roberts
 */

#include "Intrepid2_TensorData.hpp"
#include "Intrepid2_VectorData.hpp"

#ifndef Intrepid2_BasisValues_h
#define Intrepid2_BasisValues_h

/** \class  Intrepid2::BasisValues
    \brief  The data containers in Intrepid2 that support sum factorization and other reduced-data
            optimizations distinguish between scalar-valued data that is a simple product of elements
            in tensor components, and vector-valued data that is made up of a series of such products.
 
            BasisValues contains both a VectorData and a TensorData object; only one of them
            will be valid in any BasisValues object.  This allows us to maintain a common getValues() interface
            in the Basis class, rather than distinguishing at the interface level between operations that result
            in vector values and those that result in scalars.
*/

namespace Intrepid2
{
  template<class Scalar, typename ExecSpaceType>
  class BasisValues
  {
    using TensorDataType = TensorData<Scalar,ExecSpaceType>;
    using VectorDataType = VectorData<Scalar,ExecSpaceType>;
    
    Kokkos::Array<TensorDataType,Parameters::MaxTensorComponents> tensorDataFamilies_;
    VectorDataType vectorData_;
    
    int numTensorDataFamilies_ = -1;
  public:
    //! Constructor for scalar-valued BasisValues with a single family of values.
    BasisValues(TensorDataType tensorData)
    :
    tensorDataFamilies_({tensorData}),
    numTensorDataFamilies_(1)
    {}
    
    //! Constructor for scalar-valued BasisValues, with potentially multiple families of values.  (Used, e.g., for op = DIV and functionSpace = HDIV.)
    BasisValues(std::vector<TensorDataType> tensorDataFamilies)
    :
    numTensorDataFamilies_(tensorDataFamilies.size())
    {
      for (int family=0; family<numTensorDataFamilies_; family++)
      {
        tensorDataFamilies_[family] = tensorDataFamilies[family];
      }
    }
    
    //! Constructor for vector-valued BasisValues.
    BasisValues(VectorDataType vectorData)
    :
    vectorData_(vectorData)
    {}
    
    //! Default constructor.
    BasisValues()
    :
    numTensorDataFamilies_(0)
    {}
    
    
    //! copy-like constructor for differing execution spaces.  This does a deep copy of underlying views.
    template<typename OtherExecSpaceType, class = typename std::enable_if<!std::is_same<ExecSpaceType, OtherExecSpaceType>::value>::type>
    BasisValues(const BasisValues<Scalar,OtherExecSpaceType> &basisValues)
    :
    vectorData_(basisValues.vectorData()),
    numTensorDataFamilies_(basisValues.numTensorDataFamilies())
    {
      auto otherFamilies = basisValues.tensorDataFamilies();
      for (int family=0; family<numTensorDataFamilies_; family++)
      {
        tensorDataFamilies_[family] = TensorData<Scalar,ExecSpaceType>(otherFamilies[family]);
      }
    }
    
    //! field start and length must align with families in vectorData_ or tensorDataFamilies_ (whichever is valid).
    BasisValues<Scalar,ExecSpaceType> basisValuesForFields(const int &fieldStartOrdinal, const int &numFields)
    {
      int familyStartOrdinal = -1, familyEndOrdinal = -1;
      const int familyCount = this->numFamilies();
      int fieldsSoFar = 0;
      for (int i=0; i<familyCount; i++)
      {
        const bool startMatches = (fieldsSoFar == fieldStartOrdinal);
        familyStartOrdinal      = startMatches ? i : familyStartOrdinal;
        fieldsSoFar            += numFieldsInFamily(i);
        const bool endMatches   = (fieldsSoFar - fieldStartOrdinal == numFields);
        familyEndOrdinal        = endMatches ? i : familyEndOrdinal;
      }
      INTREPID2_TEST_FOR_EXCEPTION(familyStartOrdinal == -1, std::invalid_argument, "fieldStartOrdinal does not align with the start of a family.");
      INTREPID2_TEST_FOR_EXCEPTION(familyEndOrdinal   == -1, std::invalid_argument, "fieldStartOrdinal + numFields does not align with the end of a family.");
      
      const int numFamiliesInFieldSpan = familyEndOrdinal - familyStartOrdinal + 1;
      if (numTensorDataFamilies_ > 0)
      {
        std::vector<TensorDataType> tensorDataFamilies(numFamiliesInFieldSpan);
        for (int i=familyStartOrdinal; i<=familyEndOrdinal; i++)
        {
          tensorDataFamilies[i-familyStartOrdinal] = tensorDataFamilies_[i];
        }
        return BasisValues<Scalar,ExecSpaceType>(tensorDataFamilies);
      }
      else
      {
        const int componentCount = vectorData_.numComponents();
        std::vector< std::vector<TensorData<Scalar,ExecSpaceType> > > vectorComponents(numFamiliesInFieldSpan, std::vector<TensorData<Scalar,ExecSpaceType> >(componentCount));
        for (int i=familyStartOrdinal; i<=familyEndOrdinal; i++)
        {
          for (int j=0; j<componentCount; j++)
          {
            vectorComponents[i-familyStartOrdinal][j] = vectorData_.getComponent(i,j);
          }
        }
        return BasisValues<Scalar,ExecSpaceType>(vectorComponents);
      }
    }
    
    //! TensorData accessor for single-family scalar data
    TensorDataType & tensorData()
    {
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(numTensorDataFamilies_ != 1, std::invalid_argument, "this method is not supported when numTensorDataFamilies_ != 1");
      return tensorDataFamilies_[0];
    }
    
    //! TensorData accessor for multi-family scalar data
    TensorDataType & tensorData(const int &familyOrdinal)
    {
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(familyOrdinal >= numTensorDataFamilies_, std::invalid_argument, "familyOrdinal too large");
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(familyOrdinal < 0, std::invalid_argument, "familyOrdinal may not be less than 0");
      return tensorDataFamilies_[familyOrdinal];
    }
    
    //! For valid vectorData, returns the number of families in vectorData; otherwise, returns number of TensorData families
    KOKKOS_INLINE_FUNCTION
    int numFamilies() const
    {
      if (vectorData_.isValid())
      {
        return vectorData_.numFamilies();
      }
      else
      {
        return numTensorDataFamilies_;
      }
    }
    
    KOKKOS_INLINE_FUNCTION
    int numTensorDataFamilies() const
    {
      return numTensorDataFamilies_;
    }
    
    KOKKOS_INLINE_FUNCTION
    int numFieldsInFamily(int familyOrdinal) const
    {
      if (vectorData_.isValid())
      {
        return vectorData_.numFieldsInFamily(familyOrdinal);
      }
      else
      {
        return tensorDataFamilies_[familyOrdinal].extent_int(0); // (F,P,…)
      }
    }
    
    //! TensorDataFamilies accessor
    const Kokkos::Array<TensorDataType,Parameters::MaxTensorComponents> & tensorDataFamilies() const
    {
      return tensorDataFamilies_;
    }
    
    //! VectorData accessor
    const VectorDataType & vectorData() const
    {
      return vectorData_;
    }
    
    //! operator() for (F,P) scalar data; throws an exception if this is not a scalar-valued container
    KOKKOS_INLINE_FUNCTION
    Scalar operator()(const int &fieldOrdinal, const int &pointOrdinal) const
    {
      if (numTensorDataFamilies_ == 1)
      {
#ifdef HAVE_INTREPID2_DEBUG
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(! tensorDataFamilies_[0].isValid(), std::invalid_argument, "TensorData object not initialized!");
#endif
        return tensorDataFamilies_[0](fieldOrdinal, pointOrdinal);
      }
      else
      {
        int familyForField = -1;
        int previousFamilyEnd = -1;
        int fieldAdjustment = 0;
        // this loop is written in such a way as to avoid branching for CUDA performance
        for (int family=0; family<numTensorDataFamilies_; family++)
        {
          const int familyFieldCount = tensorDataFamilies_[family].extent_int(0);
          const bool fieldInRange    = (fieldOrdinal > previousFamilyEnd) && (fieldOrdinal <= previousFamilyEnd + familyFieldCount);
          familyForField = fieldInRange ? family : familyForField;
          fieldAdjustment = fieldInRange ? previousFamilyEnd + 1 : fieldAdjustment;
          previousFamilyEnd += familyFieldCount;
        }
#ifdef HAVE_INTREPID2_DEBUG
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE( familyForField == -1, std::invalid_argument, "fieldOrdinal appears to be out of range");
#endif
        return tensorDataFamilies_[familyForField](fieldOrdinal-fieldAdjustment,pointOrdinal);
      }
    }
    
    //! operator() for (F,P,D) vector data; throws an exception if this is not a vector-valued container
    KOKKOS_INLINE_FUNCTION
    Scalar operator()(const int &fieldOrdinal, const int &pointOrdinal, const int &dim) const
    {
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(! vectorData_.isValid(), std::invalid_argument, "VectorData object not initialized!");
#endif
      return vectorData_(fieldOrdinal, pointOrdinal, dim);
    }
    
    //! operator() for (C,F,P,D) data, which arises in CVFEM; at present unimplemented, and only declared here to allow a generic setJacobian() method in CellTools to compile.
    KOKKOS_INLINE_FUNCTION
    Scalar operator()(const int &cellOrdinal, const int &fieldOrdinal, const int &pointOrdinal, const int &dim) const
    {
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "CVFEM support not yet implemented in BasisValues");
      return 0;
    }
    
    KOKKOS_INLINE_FUNCTION
    int extent_int(const int &i) const
    {
      // shape is (F,P) or (F,P,D)
      if (i == 0) // field dimension
      {
        int numFields = 0;
        for (int familyOrdinal=0; familyOrdinal<numFamilies(); familyOrdinal++)
        {
          numFields += numFieldsInFamily(familyOrdinal);
        }
        return numFields;
      }
      else
      {
        if (vectorData_.isValid())
        {
          return vectorData_.extent_int(i);
        }
        else if (tensorDataFamilies_[0].isValid())
        {
          return tensorDataFamilies_[0].extent_int(i);
        }
        else
        {
          return 0;
        }
      }
    }
    
    
    KOKKOS_INLINE_FUNCTION
    size_t extent(const int &i) const
    {
      return static_cast<size_t>(extent_int(i));
    }
    
    KOKKOS_INLINE_FUNCTION
    size_t rank() const
    {
      if (vectorData_.isValid())
      {
        return vectorData_.rank();
      }
      else if (tensorDataFamilies_[0].isValid())
      {
        return tensorDataFamilies_[0].rank();
      }
      else
      {
        return 0;
      }
    }
  };
}

#endif /* Intrepid2_BasisValues_h */
