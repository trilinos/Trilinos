// @HEADER
// ************************************************************************
//
//                           Intrepid Package
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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_FieldContainerDef.hpp
    \brief  Definition file for utility class to provide multidimensional containers.
    \author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid {


  //--------------------------------------------------------------------------------------------//
  //                                                                                            //
  //                 Member function definitions of the class FieldContainer                    //
  //                                                                                            //
  //--------------------------------------------------------------------------------------------//


template<class Scalar, int ArrayTypeId>
FieldContainer<Scalar, ArrayTypeId>::FieldContainer(const FieldContainer<Scalar, ArrayTypeId>& right) {

  // Copy dimensions and data values from right
  dimensions_.assign(right.dimensions_.begin(),right.dimensions_.end());
  data_.assign(right.data_.begin(),right.data_.end());
  dim0_ = right.dim0_;
  dim1_ = right.dim1_;
  dim2_ = right.dim2_;
  dim3_ = right.dim3_;
  dim4_ = right.dim4_;
  data_ptr_ = data_.begin();
}

//--------------------------------------------------------------------------------------------//
//                                                                                            //
//                              Constructors of FieldContainer class                          //
//                                                                                            //
//--------------------------------------------------------------------------------------------//


template<class Scalar, int ArrayTypeId>
FieldContainer<Scalar, ArrayTypeId>::FieldContainer(const int dim0) : dim0_(dim0), dim1_(0), dim2_(0), dim3_(0), dim4_(0)
{
  using Teuchos::as;
  using Teuchos::Ordinal;
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (0 > dim0), std::invalid_argument,
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative dimension.");

#endif
  dimensions_.resize(as<Ordinal>(1));
  dimensions_[0] = dim0_;
  data_.assign(as<Ordinal>(dim0_), as<Scalar>(0));
  data_ptr_ = data_.begin();
}



template<class Scalar, int ArrayTypeId>
FieldContainer<Scalar, ArrayTypeId>::FieldContainer(const int dim0,
                                       const int dim1) : dim0_(dim0), dim1_(dim1), dim2_(0), dim3_(0), dim4_(0)
{
  using Teuchos::as;
  using Teuchos::Ordinal;
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (0 > dim0), std::invalid_argument,
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 1st dimension.");
  TEUCHOS_TEST_FOR_EXCEPTION( (0 > dim1), std::invalid_argument,
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 2nd dimension.");

#endif
  dimensions_.resize(2);
  dimensions_[0] = dim0_;
  dimensions_[1] = dim1_;
  data_.assign(as<Ordinal>(dim0_*dim1_), as<Scalar>(0));
  data_ptr_ = data_.begin();
}



template<class Scalar, int ArrayTypeId>
FieldContainer<Scalar, ArrayTypeId>::FieldContainer(const int dim0,
                                       const int dim1,
                                       const int dim2) : dim0_(dim0), dim1_(dim1), dim2_(dim2), dim3_(0), dim4_(0)
{
  using Teuchos::as;
  using Teuchos::Ordinal;
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (0 > dim0), std::invalid_argument,
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 1st dimension.");
  TEUCHOS_TEST_FOR_EXCEPTION( (0 > dim1), std::invalid_argument,
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 2nd dimension.");
  TEUCHOS_TEST_FOR_EXCEPTION( (0 > dim2), std::invalid_argument,
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 3rd dimension.");
#endif
  dimensions_.resize(3);
  dimensions_[0] = dim0_;
  dimensions_[1] = dim1_;
  dimensions_[2] = dim2_;
  data_.assign(as<Ordinal>(dim0_*dim1_*dim2_), as<Scalar>(0));
  data_ptr_ = data_.begin();
}



template<class Scalar, int ArrayTypeId>
FieldContainer<Scalar, ArrayTypeId>::FieldContainer(const int dim0,
                                       const int dim1,
                                       const int dim2,
                                       const int dim3) : dim0_(dim0), dim1_(dim1), dim2_(dim2), dim3_(dim3), dim4_(0)
{
  using Teuchos::as;
  using Teuchos::Ordinal;
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (0 > dim0), std::invalid_argument,
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 1st dimension.");
  TEUCHOS_TEST_FOR_EXCEPTION( (0 > dim1), std::invalid_argument,
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 2nd dimension.");
  TEUCHOS_TEST_FOR_EXCEPTION( (0 > dim2), std::invalid_argument,
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 3rd dimension.");
  TEUCHOS_TEST_FOR_EXCEPTION( (0 > dim3), std::invalid_argument,
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 4th dimension.");
#endif
  dimensions_.resize(4);
  dimensions_[0] = dim0_;
  dimensions_[1] = dim1_;
  dimensions_[2] = dim2_;
  dimensions_[3] = dim3_;
  data_.assign(as<Ordinal>(dim0_*dim1_*dim2_*dim3_), as<Scalar>(0));
  data_ptr_ = data_.begin();
}



template<class Scalar, int ArrayTypeId>
FieldContainer<Scalar, ArrayTypeId>::FieldContainer(const int dim0,
                                       const int dim1,
                                       const int dim2,
                                       const int dim3,
                                       const int dim4) : dim0_(dim0), dim1_(dim1), dim2_(dim2), dim3_(dim3), dim4_(dim4)
{
  using Teuchos::as;
  using Teuchos::Ordinal;
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (0 > dim0), std::invalid_argument,
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 1st dimension.");
  TEUCHOS_TEST_FOR_EXCEPTION( (0 > dim1), std::invalid_argument,
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 2nd dimension.");
  TEUCHOS_TEST_FOR_EXCEPTION( (0 > dim2), std::invalid_argument,
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 3rd dimension.");
  TEUCHOS_TEST_FOR_EXCEPTION( (0 > dim3), std::invalid_argument,
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 4th dimension.");
  TEUCHOS_TEST_FOR_EXCEPTION( (0 > dim4), std::invalid_argument,
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 5th dimension.");
#endif
  dimensions_.resize(5);
  dimensions_[0] = dim0_;
  dimensions_[1] = dim1_;
  dimensions_[2] = dim2_;
  dimensions_[3] = dim3_;
  dimensions_[4] = dim4_;
  data_.assign(as<Ordinal>(dim0_*dim1_*dim2_*dim3_*dim4_), as<Scalar>(0));
  data_ptr_ = data_.begin();
}



template<class Scalar, int ArrayTypeId>
FieldContainer<Scalar, ArrayTypeId>::FieldContainer(const Teuchos::Array<int>& dimensionsArray) {

#ifdef HAVE_INTREPID_DEBUG
// srkenno@sandia.gov 6/12/10: changed unsigned int to int - this was causing a warning on compilers that
//   signed & unsigned int's were being comparied.
  for( int dim = 0; dim < dimensionsArray.size(); dim++) {
    TEUCHOS_TEST_FOR_EXCEPTION( (0 > dimensionsArray[dim] ), std::invalid_argument,
                        ">>> ERROR (FieldContainer): One or more negative dimensions");
  }
#endif

  // Copy dimensions and resize container storage to match them
  dimensions_.assign(dimensionsArray.begin(),dimensionsArray.end());

  // Copy first 5 dimensions to optimize class for low rank containers
  unsigned int theRank = dimensions_.size();
  switch(theRank) {
    case 1:
      dim0_ = dimensions_[0];
      dim1_ = 0;
      dim2_ = 0;
      dim3_ = 0;
      dim4_ = 0;
      break;

    case 2:
      dim0_ = dimensions_[0];
      dim1_ = dimensions_[1];
      dim2_ = 0;
      dim3_ = 0;
      dim4_ = 0;
      break;

    case 3:
      dim0_ = dimensions_[0];
      dim1_ = dimensions_[1];
      dim2_ = dimensions_[2];
      dim3_ = 0;
      dim4_ = 0;
      break;

    case 4:
      dim0_ = dimensions_[0];
      dim1_ = dimensions_[1];
      dim2_ = dimensions_[2];
      dim3_ = dimensions_[3];
      dim4_ = 0;
      break;

    case 5:
    default:
      dim0_ = dimensions_[0];
      dim1_ = dimensions_[1];
      dim2_ = dimensions_[2];
      dim3_ = dimensions_[3];
      dim4_ = dimensions_[4];
  }

  // resize data array according to specified dimensions
  data_.assign( this -> size(), (Scalar)0);
  data_ptr_ = data_.begin();

}



template<class Scalar, int ArrayTypeId>
FieldContainer<Scalar, ArrayTypeId>::FieldContainer(const Teuchos::Array<int>&         dimensionsArray,
                                       const Teuchos::ArrayView<Scalar>&  data) {

  // Copy all dimensions
  dimensions_.assign(dimensionsArray.begin(),dimensionsArray.end());

  // Copy first 5 dimensions to optimize class for low rank containers
  unsigned int theRank = dimensions_.size();
  switch (theRank) {
    case 0:
      dim0_ = 0;
      dim1_ = 0;
      dim2_ = 0;
      dim3_ = 0;
      dim4_ = 0;
      break;

    case 1:
      dim0_ = dimensions_[0];
      dim1_ = 0;
      dim2_ = 0;
      dim3_ = 0;
      dim4_ = 0;
      break;

    case 2:
      dim0_ = dimensions_[0];
      dim1_ = dimensions_[1];
      dim2_ = 0;
      dim3_ = 0;
      dim4_ = 0;
      break;

    case 3:
      dim0_ = dimensions_[0];
      dim1_ = dimensions_[1];
      dim2_ = dimensions_[2];
      dim3_ = 0;
      dim4_ = 0;
      break;

    case 4:
      dim0_ = dimensions_[0];
      dim1_ = dimensions_[1];
      dim2_ = dimensions_[2];
      dim3_ = dimensions_[3];
      dim4_ = 0;
      break;

    case 5:
    default:
      dim0_ = dimensions_[0];
      dim1_ = dimensions_[1];
      dim2_ = dimensions_[2];
      dim3_ = dimensions_[3];
      dim4_ = dimensions_[4];
  }

    // Validate input: size of data array must match container size specified by its dimensions
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( (int)data.size() != this -> size() ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Size of input data does not match size of this container.");
#endif

  // Deep-copy ArrayView data.
  data_.assign(data.begin(),data.end());
  data_ptr_ = data_.begin();
}



template<class Scalar, int ArrayTypeId>
FieldContainer<Scalar, ArrayTypeId>::FieldContainer(const Teuchos::Array<int>&        dimensionsArray,
                                       const Teuchos::ArrayRCP<Scalar>&  data) {

  // Copy all dimensions
  dimensions_.assign(dimensionsArray.begin(),dimensionsArray.end());

  // Copy first 5 dimensions to optimize class for low rank containers
  unsigned int theRank = dimensions_.size();
  switch(theRank) {
    case 0:
      dim0_ = 0;
      dim1_ = 0;
      dim2_ = 0;
      dim3_ = 0;
      dim4_ = 0;
      break;

    case 1:
      dim0_ = dimensions_[0];
      dim1_ = 0;
      dim2_ = 0;
      dim3_ = 0;
      dim4_ = 0;
      break;

    case 2:
      dim0_ = dimensions_[0];
      dim1_ = dimensions_[1];
      dim2_ = 0;
      dim3_ = 0;
      dim4_ = 0;
      break;

    case 3:
      dim0_ = dimensions_[0];
      dim1_ = dimensions_[1];
      dim2_ = dimensions_[2];
      dim3_ = 0;
      dim4_ = 0;
      break;

    case 4:
      dim0_ = dimensions_[0];
      dim1_ = dimensions_[1];
      dim2_ = dimensions_[2];
      dim3_ = dimensions_[3];
      dim4_ = 0;
      break;

    case 5:
    default:
      dim0_ = dimensions_[0];
      dim1_ = dimensions_[1];
      dim2_ = dimensions_[2];
      dim3_ = dimensions_[3];
      dim4_ = dimensions_[4];
  }

    // Validate input: size of data array must match container size specified by its dimensions
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( (int)data.size() != this -> size() ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Size of input data does not match size of this container.");
#endif

  // Shallow-copy ArrayRCP data.
  data_ = data;
  data_ptr_ = data_.begin();
}



template<class Scalar, int ArrayTypeId>
FieldContainer<Scalar, ArrayTypeId>::FieldContainer(const Teuchos::Array<int>&    dimensionsArray,
                                       Scalar*                       data,
                                       const bool                    deep_copy,
                                       const bool                    owns_mem) {

  // Copy all dimensions
  dimensions_.assign(dimensionsArray.begin(),dimensionsArray.end());

  // Copy first 5 dimensions to optimize class for low rank containers
  unsigned int theRank = dimensions_.size();
  switch (theRank) {
    case 0:
      dim0_ = 0;
      dim1_ = 0;
      dim2_ = 0;
      dim3_ = 0;
      dim4_ = 0;
      break;

    case 1:
      dim0_ = dimensions_[0];
      dim1_ = 0;
      dim2_ = 0;
      dim3_ = 0;
      dim4_ = 0;
      break;

    case 2:
      dim0_ = dimensions_[0];
      dim1_ = dimensions_[1];
      dim2_ = 0;
      dim3_ = 0;
      dim4_ = 0;
      break;

    case 3:
      dim0_ = dimensions_[0];
      dim1_ = dimensions_[1];
      dim2_ = dimensions_[2];
      dim3_ = 0;
      dim4_ = 0;
      break;

    case 4:
      dim0_ = dimensions_[0];
      dim1_ = dimensions_[1];
      dim2_ = dimensions_[2];
      dim3_ = dimensions_[3];
      dim4_ = 0;
      break;

    case 5:
    default:
      dim0_ = dimensions_[0];
      dim1_ = dimensions_[1];
      dim2_ = dimensions_[2];
      dim3_ = dimensions_[3];
      dim4_ = dimensions_[4];
  }


  if (deep_copy) {
    Teuchos::ArrayRCP<Scalar> arrayrcp = Teuchos::arcp<Scalar>(data, 0, this -> size(), false);
    data_.deepCopy(arrayrcp());
    data_ptr_ = data_.begin();
  }
  else {
    data_ = Teuchos::arcp<Scalar>(data, 0, this -> size(), owns_mem);
    data_ptr_ = data_.begin();
  }
}



template<class Scalar, int ArrayTypeId>
FieldContainer<Scalar, ArrayTypeId>::FieldContainer(const shards::Array<Scalar,shards::NaturalOrder>&  data,
                                       const bool                                         deep_copy,
                                       const bool                                         owns_mem) {

  // Copy all dimensions
  dimensions_.resize(data.rank());

  // Copy first 5 dimensions to optimize class for low rank containers
  unsigned int theRank = dimensions_.size();
  switch(theRank) {
    case 1:
      dimensions_[0] = data.dimension(0);
      dim0_ = dimensions_[0];
      dim1_ = 0;
      dim2_ = 0;
      dim3_ = 0;
      dim4_ = 0;
      break;

    case 2:
      dimensions_[0] = data.dimension(0);
      dimensions_[1] = data.dimension(1);
      dim0_ = dimensions_[0];
      dim1_ = dimensions_[1];
      dim2_ = 0;
      dim3_ = 0;
      dim4_ = 0;
      break;

    case 3:
      dimensions_[0] = data.dimension(0);
      dimensions_[1] = data.dimension(1);
      dimensions_[2] = data.dimension(2);
      dim0_ = dimensions_[0];
      dim1_ = dimensions_[1];
      dim2_ = dimensions_[2];
      dim3_ = 0;
      dim4_ = 0;
      break;

    case 4:
      dimensions_[0] = data.dimension(0);
      dimensions_[1] = data.dimension(1);
      dimensions_[2] = data.dimension(2);
      dimensions_[3] = data.dimension(3);
      dim0_ = dimensions_[0];
      dim1_ = dimensions_[1];
      dim2_ = dimensions_[2];
      dim3_ = dimensions_[3];
      dim4_ = 0;
      break;

    case 5:
      dimensions_[0] = data.dimension(0);
      dimensions_[1] = data.dimension(1);
      dimensions_[2] = data.dimension(2);
      dimensions_[3] = data.dimension(3);
      dimensions_[4] = data.dimension(4);
      dim0_ = dimensions_[0];
      dim1_ = dimensions_[1];
      dim2_ = dimensions_[2];
      dim3_ = dimensions_[3];
      dim4_ = dimensions_[4];
      break;

    default:
      for (int i=0; i<data.rank(); i++) {
        dimensions_[i] = data.dimension(i);
      }
      dim0_ = dimensions_[0];
      dim1_ = dimensions_[1];
      dim2_ = dimensions_[2];
      dim3_ = dimensions_[3];
      dim4_ = dimensions_[4];
  }


  if (deep_copy) {
    Teuchos::ArrayRCP<Scalar> arrayrcp = Teuchos::arcp<Scalar>(data.contiguous_data(), 0, this -> size(), false);
    data_.deepCopy(arrayrcp());
    data_ptr_ = data_.begin();
  }
  else {
    data_ = Teuchos::arcp<Scalar>(data.contiguous_data(), 0, this -> size(), owns_mem);
    data_ptr_ = data_.begin();
  }
}



//--------------------------------------------------------------------------------------------//
//                                                                                            //
//                            Access methods of FieldContainer class                          //
//                                                                                            //
//--------------------------------------------------------------------------------------------//


template<class Scalar, int ArrayTypeId>
inline int FieldContainer<Scalar, ArrayTypeId>::rank() const {
  return dimensions_.size();
}



template<class Scalar, int ArrayTypeId>
int FieldContainer<Scalar, ArrayTypeId>::size() const {
  // Important! This method is used by constructors to find out what is the needed size of data_
  // based on the specified dimensions. Therefore, it cannot be implmented by returning data_.size
  // and must be able to compute the size of the container based only on its specified dimensions

  // Size equals product of all dimensions stored in dimensions_
  const int theRank = dimensions_.size ();

  // If container has no dimensions its size is zero
  if (theRank == 0) {
    return 0;
  }
  else {
    int theSize = dim0_;

    // Compute size directly to optimize method for low rank (<=5) containers
    switch(theRank) {
    case 5:
      theSize *= dim1_*dim2_*dim3_*dim4_;
      break;

    case 4:
      theSize *= dim1_*dim2_*dim3_;
      break;

    case 3:
      theSize *= dim1_*dim2_;
      break;

    case 2:
      theSize *= dim1_;
      break;

    case 1:
      break;

      // Compute size for containers with ranks hihger than 5
    default:
      for (int r = 1; r < theRank; ++r) {
        theSize *= dimensions_[r];
      }
    }
    return theSize;
  }
}



template<class Scalar, int ArrayTypeId>
template<class Vector>
inline void FieldContainer<Scalar, ArrayTypeId>::dimensions(Vector& dimVec) const {
  dimVec.assign(dimensions_.begin(),dimensions_.end());
}



template<class Scalar, int ArrayTypeId>
inline int FieldContainer<Scalar, ArrayTypeId>::dimension(const int whichDim) const {
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (0 > whichDim), std::invalid_argument,
                      ">>> ERROR (FieldContainer): dimension order cannot be negative");
  TEUCHOS_TEST_FOR_EXCEPTION( (whichDim >= this -> rank() ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): dimension order cannot exceed rank of the container");
#endif
  return dimensions_[whichDim];
}



template<class Scalar, int ArrayTypeId>
inline int FieldContainer<Scalar, ArrayTypeId>::getEnumeration(const int i0) const {
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( this -> rank() != 1), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i0 < 0) || ( i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): index is out of range.");
#endif
  return i0;
}



template<class Scalar, int ArrayTypeId>
inline int FieldContainer<Scalar, ArrayTypeId>::getEnumeration(const int i0,
                                                  const int i1) const {
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( this -> rank() != 2), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i0 < 0) || ( i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");
#endif
  return i0*dim1_ + i1;
}



template<class Scalar, int ArrayTypeId>
inline int FieldContainer<Scalar, ArrayTypeId>::getEnumeration(const int i0,
                                                  const int i1,
                                                  const int i2) const {
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( this -> rank() != 3), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i0 < 0) || ( i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 3rd index is out of range.");
#endif
  return (i0*dim1_ + i1)*dim2_ + i2;
}



template<class Scalar, int ArrayTypeId>
inline int FieldContainer<Scalar, ArrayTypeId>::getEnumeration(const int i0,
                                                  const int i1,
                                                  const int i2,
                                                  const int i3) const {
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( this -> rank() != 4), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i0 < 0) || ( i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 3rd index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i3 < 0) || (i3 >= dim3_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 4th index is out of range.");
#endif
  return ( (i0*dim1_ + i1 )*dim2_ + i2 )*dim3_ + i3;
}



template<class Scalar, int ArrayTypeId>
inline int FieldContainer<Scalar, ArrayTypeId>::getEnumeration(const int i0,
                                                  const int i1,
                                                  const int i2,
                                                  const int i3,
                                                  const int i4) const {
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( this -> rank() != 5), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i0 < 0) || ( i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 3rd index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i3 < 0) || (i3 >= dim3_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 4th index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i4 < 0) || (i4 >= dim4_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 5th index is out of range.");
#endif
  return ( ( (i0*dim1_ + i1 )*dim2_ + i2 )*dim3_ + i3 )*dim4_ + i4;
}




template<class Scalar, int ArrayTypeId>
int FieldContainer<Scalar, ArrayTypeId>::getEnumeration(const Teuchos::Array<int>& multiIndex) const {

#ifdef HAVE_INTREPID_DEBUG
  // Check if number of multi-indices matches rank of the FieldContainer object
  TEUCHOS_TEST_FOR_EXCEPTION( ( multiIndex.size() != dimensions_.size() ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of multi-indices does not match rank of container.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( ( multiIndex[0] < 0) || ( multiIndex[0] >= dim0_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");
#endif

  int theRank = dimensions_.size();
  int address = 0;
  switch (theRank) {

    // Optimize enumeration computation for low rank (<= 5) containers
    case 5:
#ifdef HAVE_INTREPID_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION( ( (multiIndex[4] < 0) || (multiIndex[4] >= dim4_) ),
                          std::invalid_argument,
                          ">>> ERROR (FieldContainer): 5th index is out of range.");
      TEUCHOS_TEST_FOR_EXCEPTION( ( (multiIndex[3] < 0) || (multiIndex[3] >= dim3_) ),
                          std::invalid_argument,
                          ">>> ERROR (FieldContainer): 4th index is out of range.");
      TEUCHOS_TEST_FOR_EXCEPTION( ( (multiIndex[2] < 0) || (multiIndex[2] >= dim2_) ),
                          std::invalid_argument,
                          ">>> ERROR (FieldContainer): 3rd index is out of range.");
      TEUCHOS_TEST_FOR_EXCEPTION( ( (multiIndex[1] < 0) || (multiIndex[1] >= dim1_) ),
                          std::invalid_argument,
                          ">>> ERROR (FieldContainer): 2nd index is out of range.");
#endif
      address = (((multiIndex[0]*dim1_ + multiIndex[1])*dim2_ + multiIndex[2])*dim3_ + multiIndex[3])*dim4_ + multiIndex[4];
      break;

    case 4:
#ifdef HAVE_INTREPID_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION( ( (multiIndex[3] < 0) || (multiIndex[3] >= dim3_) ),
                          std::invalid_argument,
                          ">>> ERROR (FieldContainer): 4th index is out of range.");
      TEUCHOS_TEST_FOR_EXCEPTION( ( (multiIndex[2] < 0) || (multiIndex[2] >= dim2_) ),
                          std::invalid_argument,
                          ">>> ERROR (FieldContainer): 3rd index is out of range.");
      TEUCHOS_TEST_FOR_EXCEPTION( ( (multiIndex[1] < 0) || (multiIndex[1] >= dim1_) ),
                          std::invalid_argument,
                          ">>> ERROR (FieldContainer): 2nd index is out of range.");
#endif
      address = ((multiIndex[0]*dim1_ + multiIndex[1])*dim2_ + multiIndex[2])*dim3_ + multiIndex[3];
      break;

    case 3:
#ifdef HAVE_INTREPID_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION( ( (multiIndex[2] < 0) || (multiIndex[2] >= dim2_) ),
                          std::invalid_argument,
                          ">>> ERROR (FieldContainer): 3rd index is out of range.");
      TEUCHOS_TEST_FOR_EXCEPTION( ( (multiIndex[1] < 0) || (multiIndex[1] >= dim1_) ),
                          std::invalid_argument,
                          ">>> ERROR (FieldContainer): 2nd index is out of range.");
#endif
      address = (multiIndex[0]*dim1_ + multiIndex[1])*dim2_ + multiIndex[2];
      break;

    case 2:
#ifdef HAVE_INTREPID_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION( ( (multiIndex[1] < 0) || (multiIndex[1] >= dim1_) ),
                          std::invalid_argument,
                          ">>> ERROR (FieldContainer): 2nd index is out of range.");
#endif
      address = multiIndex[0]*dim1_ + multiIndex[1];
      break;

    case 1:
      address = multiIndex[0];
      break;

    default:

      // Arbitrary rank: compute address using Horner's nested scheme: intialize address to 0th index value
      address = multiIndex[0];
      for (int r = 0; r < theRank - 1; r++){
#ifdef HAVE_INTREPID_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION( ( (multiIndex[r+1] < 0) || (multiIndex[r+1] >= dimensions_[r+1]) ),
                            std::invalid_argument,
                            ">>> ERROR (FieldContainer): Multi-index component out of range.");
#endif
        // Add increment
        address = address*dimensions_[r+1] + multiIndex[r+1];
      }
  } // end switch(theRank)

  return address;
}



template<class Scalar, int ArrayTypeId>
void FieldContainer<Scalar, ArrayTypeId>::getMultiIndex(int & i0,
                                           const int valueEnum) const
{
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( this -> rank() != 1), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (valueEnum < 0) || (valueEnum >= (int)data_.size()) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Value enumeration is out of range.");
#endif
  i0 = valueEnum;
}



template<class Scalar, int ArrayTypeId>
void FieldContainer<Scalar, ArrayTypeId>::getMultiIndex(int & i0,
                                           int & i1,
                                           const int valueEnum) const
{
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( this -> rank() != 2), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (valueEnum < 0) || (valueEnum >= (int)data_.size()) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Value enumeration is out of range.");
#endif

  i0 = valueEnum/dim1_;
  i1 = valueEnum - i0*dim1_;
}



template<class Scalar, int ArrayTypeId>
void FieldContainer<Scalar, ArrayTypeId>::getMultiIndex(int & i0,
                                           int & i1,
                                           int & i2,
                                           const int valueEnum) const
{
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( this -> rank() != 3), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (valueEnum < 0) || (valueEnum >= (int)data_.size()) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Value enumeration is out of range.");
#endif
  int tempDim = dim1_*dim2_;
  int tempEnu = valueEnum;
  i0 = tempEnu/tempDim;

  tempEnu -= i0*tempDim;
  tempDim /= dim1_;
  i1 = tempEnu/tempDim;

  tempEnu -= i1*tempDim;
  tempDim /= dim2_;
  i2 = tempEnu/tempDim;
}



template<class Scalar, int ArrayTypeId>
void FieldContainer<Scalar, ArrayTypeId>::getMultiIndex(int & i0,
                                           int & i1,
                                           int & i2,
                                           int & i3,
                                           const int valueEnum) const
{
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( this -> rank() != 4), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (valueEnum < 0) || (valueEnum >= (int)data_.size()) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Value enumeration is out of range.");
#endif
  int tempDim = dim1_*dim2_*dim3_;
  int tempEnu = valueEnum;
  i0 = tempEnu/tempDim;

  tempEnu -= i0*tempDim;
  tempDim /= dim1_;
  i1 = tempEnu/tempDim;

  tempEnu -= i1*tempDim;
  tempDim /= dim2_;
  i2 = tempEnu/tempDim;

  tempEnu -= i2*tempDim;
  tempDim /= dim3_;
  i3 = tempEnu/tempDim;
}




template<class Scalar, int ArrayTypeId>
void FieldContainer<Scalar, ArrayTypeId>::getMultiIndex(int & i0,
                                           int & i1,
                                           int & i2,
                                           int & i3,
                                           int & i4,
                                           const int valueEnum) const
{
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( this -> rank() != 5), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (valueEnum < 0) || (valueEnum >= (int)data_.size()) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Value enumeration is out of range.");
#endif
  int tempDim = dim1_*dim2_*dim3_*dim4_;
  int tempEnu = valueEnum;
  i0 = tempEnu/tempDim;

  tempEnu -= i0*tempDim;
  tempDim /= dim1_;
  i1 = tempEnu/tempDim;

  tempEnu -= i1*tempDim;
  tempDim /= dim2_;
  i2 = tempEnu/tempDim;

  tempEnu -= i2*tempDim;
  tempDim /= dim3_;
  i3 = tempEnu/tempDim;

  tempEnu -= i3*tempDim;
  tempDim /= dim4_;
  i4 = tempEnu/tempDim;
}



template<class Scalar, int ArrayTypeId>
template<class Vector>
void FieldContainer<Scalar, ArrayTypeId>::getMultiIndex(Vector &             multiIndex,
                                           const int            valueEnum) const
{

  // Verify address is in the admissible range for this FieldContainer
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( (valueEnum < 0) || (valueEnum >= (int)data_.size()) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Value enumeration is out of range.");
#endif

  // make sure multiIndex has the right size to hold all multi-indices
  const int theRank = dimensions_.size ();
  multiIndex.resize (theRank);

  // Initializations
  int temp_enum = valueEnum;
  int temp_range = 1;

  // Compute product of all but the first upper bound
  for (int r = 1; r < theRank; ++r) {
    temp_range *= dimensions_[r];
  }

  // Index 0 is computed first using integer division
  multiIndex[0] = temp_enum/temp_range;

  // Indices 1 to (theRank - 2) are computed next; will be skipped if
  // theRank <= 2
  for (int r = 1; r < theRank - 1; ++r) {
    temp_enum  -= multiIndex[r-1]*temp_range;
    temp_range /= dimensions_[r];
    multiIndex[r] = temp_enum/temp_range;
  }

  // Index (theRank - 1) is computed last, skip if theRank = 1 and
  // keep if theRank = 2
  if (theRank > 1) {
    multiIndex[theRank - 1] = temp_enum - multiIndex[theRank - 2] * temp_range;
  }
}

//--------------------------------------------------------------------------------------------//
//                                                                                            //
//                          Methods to shape (resize) a field container                       //
//                                                                                            //
//--------------------------------------------------------------------------------------------//

template<class Scalar, int ArrayTypeId>
inline void FieldContainer<Scalar, ArrayTypeId>::clear() {
  dimensions_.resize(0);

  // Reset first five dimensions:
  dim0_ = 0;
  dim1_ = 0;
  dim2_ = 0;
  dim3_ = 0;
  dim4_ = 0;

  // Clears data array and sets to zero length
  data_.clear();
  data_ptr_ = data_.begin();
}



template<class Scalar, int ArrayTypeId>
void FieldContainer<Scalar, ArrayTypeId>::resize(const Teuchos::Array<int>& newDimensions) {

  // First handle the trivial case of zero dimensions
  if( newDimensions.size() == 0) {
    dimensions_.resize(0);
    dim0_ = 0;
    dim1_ = 0;
    dim2_ = 0;
    dim3_ = 0;
    dim4_ = 0;
    data_.resize(0);
    data_ptr_ = data_.begin();
  }
  else {

    // Copy upper index bounds and resize container storage to match new upper bounds.
    dimensions_.assign(newDimensions.begin(),newDimensions.end());

    // Copy first 5 dimensions for faster access
    unsigned int theRank = dimensions_.size();
    switch (theRank) {
      case 1:
        dim0_ = dimensions_[0];
        dim1_ = 0;
        dim2_ = 0;
        dim3_ = 0;
        dim4_ = 0;
        break;

      case 2:
        dim0_ = dimensions_[0];
        dim1_ = dimensions_[1];
        dim2_ = 0;
        dim3_ = 0;
        dim4_ = 0;
        break;

      case 3:
        dim0_ = dimensions_[0];
        dim1_ = dimensions_[1];
        dim2_ = dimensions_[2];
        dim3_ = 0;
        dim4_ = 0;
        break;

      case 4:
        dim0_ = dimensions_[0];
        dim1_ = dimensions_[1];
        dim2_ = dimensions_[2];
        dim3_ = dimensions_[3];
        dim4_ = 0;
        break;

      case 5:
      default:
        dim0_ = dimensions_[0];
        dim1_ = dimensions_[1];
        dim2_ = dimensions_[2];
        dim3_ = dimensions_[3];
        dim4_ = dimensions_[4];
    }

    // Resize data array
    data_.resize(this -> size());
    data_ptr_ = data_.begin();
  }
}



template<class Scalar, int ArrayTypeId>
inline void FieldContainer<Scalar, ArrayTypeId>::resize(const int dim0) {
  dim0_ = dim0;
  dim1_ = 0;
  dim2_ = 0;
  dim3_ = 0;
  dim4_ = 0;
  dimensions_.resize(1);
  dimensions_[0] = dim0_;
  data_.resize(dim0_);
  data_ptr_ = data_.begin();
}



template<class Scalar, int ArrayTypeId>
inline void FieldContainer<Scalar, ArrayTypeId>::resize(const int dim0,
                                           const int dim1) {
  dim0_ = dim0;
  dim1_ = dim1;
  dim2_ = 0;
  dim3_ = 0;
  dim4_ = 0;
  dimensions_.resize(2);
  dimensions_[0] = dim0_;
  dimensions_[1] = dim1_;
  data_.resize(dim0_*dim1_);
  data_ptr_ = data_.begin();
}



template<class Scalar, int ArrayTypeId>
inline void FieldContainer<Scalar, ArrayTypeId>::resize(const int dim0,
                                           const int dim1,
                                           const int dim2) {
  dim0_ = dim0;
  dim1_ = dim1;
  dim2_ = dim2;
  dim3_ = 0;
  dim4_ = 0;
  dimensions_.resize(3);
  dimensions_[0] = dim0_;
  dimensions_[1] = dim1_;
  dimensions_[2] = dim2_;
  data_.resize(dim0_*dim1_*dim2_);
  data_ptr_ = data_.begin();
}



template<class Scalar, int ArrayTypeId>
inline void FieldContainer<Scalar, ArrayTypeId>::resize(const int dim0,
                                           const int dim1,
                                           const int dim2,
                                           const int dim3) {
  dim0_ = dim0;
  dim1_ = dim1;
  dim2_ = dim2;
  dim3_ = dim3;
  dim4_ = 0;
  dimensions_.resize(4);
  dimensions_[0] = dim0_;
  dimensions_[1] = dim1_;
  dimensions_[2] = dim2_;
  dimensions_[3] = dim3_;
  data_.resize(dim0_*dim1_*dim2_*dim3_);
  data_ptr_ = data_.begin();
}



template<class Scalar, int ArrayTypeId>
inline void FieldContainer<Scalar, ArrayTypeId>::resize(const int dim0,
                                           const int dim1,
                                           const int dim2,
                                           const int dim3,
                                           const int dim4) {
  dim0_ = dim0;
  dim1_ = dim1;
  dim2_ = dim2;
  dim3_ = dim3;
  dim4_ = dim4;
  dimensions_.resize(5);
  dimensions_[0] = dim0_;
  dimensions_[1] = dim1_;
  dimensions_[2] = dim2_;
  dimensions_[3] = dim3_;
  dimensions_[4] = dim4_;
  data_.resize(dim0_*dim1_*dim2_*dim3_*dim4_);
  data_ptr_ = data_.begin();
}



template<class Scalar, int ArrayTypeId>
inline void FieldContainer<Scalar, ArrayTypeId>::resize(const FieldContainer<Scalar, ArrayTypeId>& anotherContainer) {

  // Copy dimensions from the specified container
  anotherContainer.dimensions(dimensions_);
  int newRank = dimensions_.size();

  // Copy first 5 dimensions for faster access
  switch(newRank) {
    case 1:
      dim0_ = dimensions_[0];
      dim1_ = 0;
      dim2_ = 0;
      dim3_ = 0;
      dim4_ = 0;
      break;

    case 2:
      dim0_ = dimensions_[0];
      dim1_ = dimensions_[1];
      dim2_ = 0;
      dim3_ = 0;
      dim4_ = 0;
      break;

    case 3:
      dim0_ = dimensions_[0];
      dim1_ = dimensions_[1];
      dim2_ = dimensions_[2];
      dim3_ = 0;
      dim4_ = 0;
      break;

    case 4:
      dim0_ = dimensions_[0];
      dim1_ = dimensions_[1];
      dim2_ = dimensions_[2];
      dim3_ = dimensions_[3];
      dim4_ = 0;
      break;

    case 5:
    default:
      dim0_ = dimensions_[0];
      dim1_ = dimensions_[1];
      dim2_ = dimensions_[2];
      dim3_ = dimensions_[3];
      dim4_ = dimensions_[4];
  }

  // Resize data array
  data_.resize(this->size());
  data_ptr_ = data_.begin();
}


template<class Scalar, int ArrayTypeId>
void FieldContainer<Scalar, ArrayTypeId>::resize(const int             numPoints,
                                    const int             numFields,
                                    const EFunctionSpace  spaceType,
                                    const EOperator       operatorType,
                                    const int             spaceDim) {
  // Validate input
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( numPoints < 0),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of points cannot be negative!");
  TEUCHOS_TEST_FOR_EXCEPTION( ( numFields < 0),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of fields cannot be negative!");
  TEUCHOS_TEST_FOR_EXCEPTION( !( (1 <=  spaceDim ) && ( spaceDim <= 3  ) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Invalid space dimension.");
#endif

  // Find out field and operator ranks
  const int fieldRank    = getFieldRank(spaceType);
  const int operatorRank = getOperatorRank(spaceType,operatorType,spaceDim);

  // Compute rank of the container = 1(numPoints) + 1(numFields) + fieldRank + operatorRank
  const int theRank = 1 + 1 + fieldRank + operatorRank;

  // Define temp array for the dimensions
  Teuchos::Array<int> newDimensions (theRank);

  // Dimensions 0 and 1 are number of points and number of fields, resp.
  newDimensions[0] = numPoints;
  newDimensions[1] = numFields;

  // The rest of the dimensions depend on whether we had VALUE, GRAD (D1), CURL, DIV or Dk, k>1
  switch (operatorType) {

  case OPERATOR_VALUE:
  case OPERATOR_GRAD:
  case OPERATOR_D1:
  case OPERATOR_CURL:
  case OPERATOR_DIV:

    // For these operators all dimensions from 2 to 2 + fieldRank + OperatorRank are bounded by spaceDim
    for (int i = 0; i < fieldRank + operatorRank; ++i) {
      newDimensions[2 + i] = spaceDim;
    }
    break;

  case OPERATOR_D2:
  case OPERATOR_D3:
  case OPERATOR_D4:
  case OPERATOR_D5:
  case OPERATOR_D6:
  case OPERATOR_D7:
  case OPERATOR_D8:
  case OPERATOR_D9:
  case OPERATOR_D10:

    // All dimensions from 2 to 2 + fieldRank, if any, are bounded by spaceDim
    for(int i = 0; i < fieldRank; i++){
      newDimensions[2 + i] = spaceDim;
    }

    // We know that for Dk operatorRank = 1 and so there's just one more dimension left
    // given by the cardinality of the set of all derivatives of order k
    newDimensions[2 + fieldRank] = getDkCardinality(operatorType,spaceDim);
    break;

  default:
    TEUCHOS_TEST_FOR_EXCEPTION(
      !(Intrepid::isValidOperator(operatorType) ), std::invalid_argument,
      ">>> ERROR (FieldContainer): Invalid operator type");
  }

  // Resize FieldContainer using the newDimensions in the array
  this->resize (newDimensions);
}

//--------------------------------------------------------------------------------------------//
//                                                                                            //
//                     Methods to read and write values to FieldContainer                     //
//                                                                                            //
//--------------------------------------------------------------------------------------------//



template<class Scalar, int ArrayTypeId>
inline void FieldContainer<Scalar, ArrayTypeId>::initialize(const Scalar value) {
  for (int i=0; i < this->size(); i++) {
    data_[i] = value;
  }
}



template<class Scalar, int ArrayTypeId>
inline Scalar FieldContainer<Scalar, ArrayTypeId>::getValue(const Teuchos::Array<int>& multiIndex) const {
  return data_[this -> getEnumeration(multiIndex)];
}



template<class Scalar, int ArrayTypeId>
inline void FieldContainer<Scalar, ArrayTypeId>::setValue(const Scalar dataValue,
                                             const Teuchos::Array<int>& multiIndex) {
  data_[this -> getEnumeration(multiIndex)] = dataValue;
}



template<class Scalar, int ArrayTypeId>
inline void FieldContainer<Scalar, ArrayTypeId>::setValue(const Scalar dataValue,
                                             const int    order) {
  data_[order] = dataValue;
}



template<class Scalar, int ArrayTypeId>
void FieldContainer<Scalar, ArrayTypeId>::setValues(const Teuchos::ArrayView<Scalar>& dataArray) {
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (dataArray.size() != (data_.size()) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Size of argument does not match the size of container.");
#endif
  data_.assign(dataArray.begin(),dataArray.end());
  data_ptr_ = data_.begin();
}



template<class Scalar, int ArrayTypeId>
void FieldContainer<Scalar, ArrayTypeId>::setValues(const Scalar* dataPtr,
                                       const int numData)
{
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( (numData != this -> size() ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of data does not match the size of container.");

#endif
  data_.assign(dataPtr, dataPtr + numData);
  data_ptr_ = data_.begin();
}



template<class Scalar, int ArrayTypeId>
inline const Scalar& FieldContainer<Scalar, ArrayTypeId>::operator () (const int i0) const
{
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( this -> rank() != 1), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i0 < 0) || (i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): index is out of range.");
#endif
  return data_ptr_[i0];
}


template<class Scalar, int ArrayTypeId>
inline Scalar& FieldContainer<Scalar, ArrayTypeId>::operator () (const int i0)
{
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( this -> rank() != 1), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i0 < 0) || (i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): index is out of range.");
#endif
  return data_ptr_[i0];
}



template<class Scalar, int ArrayTypeId>
inline const Scalar& FieldContainer<Scalar, ArrayTypeId>::operator () (const int i0,
                                                          const int i1) const
{
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( this -> rank() != 2), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i0 < 0) || (i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");
#endif
  return data_ptr_[i0*dim1_ + i1];
}


template<class Scalar, int ArrayTypeId>
inline Scalar& FieldContainer<Scalar, ArrayTypeId>::operator () (const int i0,
                                                    const int i1)
{
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( this -> rank() != 2), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i0 < 0) || (i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");
#endif
  return data_ptr_[i0*dim1_ + i1];
}



template<class Scalar, int ArrayTypeId>
inline const Scalar& FieldContainer<Scalar, ArrayTypeId>::operator () (const int i0,
                                                          const int i1,
                                                          const int i2) const
{
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( this -> rank() != 3), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i0 < 0) || (i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 3rd index is out of range.");
#endif
  return data_ptr_[(i0*dim1_ + i1)*dim2_ + i2];
}

template<class Scalar, int ArrayTypeId>
inline Scalar& FieldContainer<Scalar, ArrayTypeId>::operator () (const int i0,
                                                    const int i1,
                                                    const int i2)
{
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( this -> rank() != 3), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i0 < 0) || (i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 3rd index is out of range.");
#endif
  return data_ptr_[(i0*dim1_ + i1)*dim2_ + i2];
}



template<class Scalar, int ArrayTypeId>
inline const Scalar& FieldContainer<Scalar, ArrayTypeId>::operator ()  (const int i0,
                                                           const int i1,
                                                           const int i2,
                                                           const int i3) const {
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( this -> rank() != 4), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i0 < 0) || (i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 3rd index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i3 < 0) || (i3 >= dim3_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 4th index is out of range.");
#endif
  return data_ptr_[( (i0*dim1_ + i1 )*dim2_ + i2 )*dim3_ + i3];
}


template<class Scalar, int ArrayTypeId>
inline Scalar& FieldContainer<Scalar, ArrayTypeId>::operator ()  (const int i0,
                                                     const int i1,
                                                     const int i2,
                                                     const int i3) {
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( this -> rank() != 4), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i0 < 0) || (i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 3rd index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i3 < 0) || (i3 >= dim3_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 4th index is out of range.");
#endif
  return data_ptr_[( (i0*dim1_ + i1 )*dim2_ + i2 )*dim3_ + i3];
}



template<class Scalar, int ArrayTypeId>
inline const Scalar& FieldContainer<Scalar, ArrayTypeId>::operator ()  (const int i0,
                                                           const int i1,
                                                           const int i2,
                                                           const int i3,
                                                           const int i4) const {
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( this -> rank() != 5), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i0 < 0) || (i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 3rd index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i3 < 0) || (i3 >= dim3_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 4th index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i4 < 0) || (i4 >= dim4_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 5th index is out of range.");
#endif
  return data_ptr_[( ( (i0*dim1_ + i1 )*dim2_ + i2 )*dim3_ + i3 )*dim4_ + i4];
}

template<class Scalar, int ArrayTypeId>
inline Scalar& FieldContainer<Scalar, ArrayTypeId>::operator ()  (const int i0,
                                                     const int i1,
                                                     const int i2,
                                                     const int i3,
                                                     const int i4) {
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( this -> rank() != 5), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i0 < 0) || (i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 3rd index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i3 < 0) || (i3 >= dim3_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 4th index is out of range.");
  TEUCHOS_TEST_FOR_EXCEPTION( ( (i4 < 0) || (i4 >= dim4_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 5th index is out of range.");
#endif
  return data_ptr_[( ( (i0*dim1_ + i1 )*dim2_ + i2 )*dim3_ + i3 )*dim4_ + i4];
}



template<class Scalar, int ArrayTypeId>
const Scalar& FieldContainer<Scalar, ArrayTypeId>::operator [] (const int address) const {
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( (address < 0) || (address >= (int)data_.size() ) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Specified address is out of range.");
#endif
  return data_ptr_[address];
}



template<class Scalar, int ArrayTypeId>
Scalar& FieldContainer<Scalar, ArrayTypeId>::operator [] (const int address) {
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( (address < 0) || (address >= (int)data_.size() ) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Specified address is out of range.");
#endif
  return data_ptr_[address];
}



template<class Scalar, int ArrayTypeId>
inline FieldContainer<Scalar, ArrayTypeId>& FieldContainer<Scalar, ArrayTypeId>::operator = (const FieldContainer<Scalar, ArrayTypeId>& right)
{
#ifdef HAVE_INTREPID_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( ( this == &right ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Invalid right-hand side to '='. Self-assignment prohibited.");
#endif
  dim0_ = right.dim0_;
  dim1_ = right.dim1_;
  dim2_ = right.dim2_;
  dim3_ = right.dim3_;
  dim4_ = right.dim4_;
  data_.deepCopy(right.data_());
  data_ptr_ = data_.begin();
  dimensions_ = right.dimensions_;
  return *this;
}


//===========================================================================//
//                                                                           //
//           END of member definitions; START friends and related            //
//                                                                           //
//===========================================================================//


template<class Scalar, int ArrayTypeId>
std::ostream& operator << (std::ostream& os, const FieldContainer<Scalar, ArrayTypeId>& container) {

  // Save the format state of the original ostream os.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(os);

  os.setf(std::ios_base::scientific, std::ios_base::floatfield);
  os.setf(std::ios_base::right);
  int myprec = os.precision();

  int size = container.size();
  int rank = container.rank();
  Teuchos::Array<int> multiIndex(rank);
  Teuchos::Array<int> dimensions;
  container.dimensions(dimensions);

  os<< "===============================================================================\n"\
    << "\t Container size = " << size << "\n"
    << "\t Container rank = " << rank << "\n" ;

  if( (rank == 0 ) && (size == 0) ) {
    os<< "====================================================================================\n"\
      << "|                        *** This is an empty container ****                       |\n";
  }
  else {
    os<< "\t Dimensions     = ";

    for(int r = 0; r < rank; r++){
      os << " (" << dimensions[r] <<") ";
    }
    os << "\n";

    os<< "====================================================================================\n"\
      << "|              Multi-index          Enumeration             Value                  |\n"\
      << "====================================================================================\n";
  }

  for(int address = 0; address < size; address++){
    container.getMultiIndex(multiIndex,address);
    std::ostringstream mistring;
    for(int r = 0; r < rank; r++){
      mistring <<  multiIndex[r] << std::dec << " ";
    }
    os.setf(std::ios::right, std::ios::adjustfield);
    os << std::setw(27) << mistring.str();
    os << std::setw(20) << address;
    os << "             ";
    os.setf(std::ios::left, std::ios::adjustfield);
    os << std::setw(myprec+8) << container[address] << "\n";
  }

  os<< "====================================================================================\n\n";

  // reset format state of os
  os.copyfmt(oldFormatState);

  return os;
}
// End member, friend, and related function definitions of class FieldContainer.

} // end namespace Intrepid

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

