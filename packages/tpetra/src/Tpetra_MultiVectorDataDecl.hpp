// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TPETRA_MULTIVECTORDATA_DECL_HPP
#define TPETRA_MULTIVECTORDATA_DECL_HPP

#include <Teuchos_ArrayRCP.hpp>

namespace Tpetra {

  template<class Scalar>
  class MultiVectorData {

  // cannot declare as friends any partial specializations; therefore, all MultiVectors and Vectors will be friends
  template <class S2, class O1, class O2>
  friend class MultiVector;
  template <class S2, class O1, class O2>
  friend class Vector;

  public:
    MultiVectorData();
    ~MultiVectorData();

  protected:

    /* only one of the following is valid:
       if constantStride_ == false then
         - contigValues_ == Teuchos::null
         - nonContigValues_.size() == ptrs_.size()
       else, if constantStride_ == true, then
         - nonContigValues_.size() == 0
         - contigValues_ == null only if stride_ == 0

       if constantStride_ == true and contigValues_ != null, then it points to the 
       first entry in the first vector and runs to the last entry of the last vector.
       if stride == myLength(), then it has length myLength() * numVectors()
    */
    Teuchos::ArrayRCP<Scalar> contigValues_;                          
    Teuchos::Array<Teuchos::ArrayRCP<Scalar> > nonContigValues_;

    /* make sure that this is an iterator of an appropriately sized view, so that the iterator is valid only for the span of the column.
       this is filled under all circumstances, although when the data is contiguous, it may be more efficient to access data through values_ */
    mutable Teuchos::Array<typename Teuchos::ArrayView<Scalar>::iterator>  ptrs_;

    bool constantStride_;
    Teuchos_Ordinal stride_;
    
  private:
    void setupPointers(Teuchos_Ordinal MyLength, Teuchos_Ordinal NumVectors); 

    // Copy constructor (declared but not defined, do not use)
    MultiVectorData(const MultiVectorData<Scalar> &source);
    // Assignment operator (declared but not defined, do not use)
    MultiVectorData<Scalar>& operator=(const MultiVectorData<Scalar> &source);

  }; // class MultiVectorData

} // namespace Tpetra

#endif // TPETRA_MULTIVECTORDATA_DECL_HPP

