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

#ifndef TPETRA_MULTIVECTORDATA_HPP
#define TPETRA_MULTIVECTORDATA_HPP

#include "Tpetra_MultiVectorDataDecl.hpp"
#include <Teuchos_Describable.hpp>
#include <Teuchos_ArrayRCP.hpp>

namespace Tpetra {

  template<class Scalar, class Node>
  MultiVectorData<Scalar>::MultiVectorData(Node &node) : node_(node) {}

  template<class Scalar, class Node>
  MultiVectorData<Scalar>::~MultiVectorData() {}

  template<class Scalar, class Node>
  void MultiVectorData<Scalar>::setupPointers(Teuchos_Ordinal MyLength, Teuchos_Ordinal NumVectors) {
  //REFACTOR//   TEST_FOR_EXCEPTION(!(NumVectors > 0) || !(MyLength >= 0), std::logic_error,
  //REFACTOR//       "MultiVectorData::setupPointers(): logic error. Please contact Tpetra team.");
  //REFACTOR//   using Teuchos::null;
  //REFACTOR//   using Teuchos::ArrayRCP;
  //REFACTOR//   // setup ptrs_ array of pointers to iterators
  //REFACTOR//   // each iterator should be an ArrayView::iterator from a properly sized ArrayView
  //REFACTOR//   // in a debug mode, these will be ArrayRCP with correct bounds
  //REFACTOR//   // in an optimized build, they are C pointers
  //REFACTOR//   ptrs_.resize(NumVectors);
  //REFACTOR//   if (constantStride_) {
  //REFACTOR//     if (contigValues_ != null) { // stride_ > 0
  //REFACTOR//       ArrayRCP<Scalar> ptr = contigValues_;
  //REFACTOR//       for (Teuchos_Ordinal j=0; j<NumVectors; ++j) {
  //REFACTOR//         ptrs_[j] = ptr(0,MyLength).begin();
  //REFACTOR//         ptr += stride_;
  //REFACTOR//       }
  //REFACTOR//     }
  //REFACTOR//   }
  //REFACTOR//   else {
  //REFACTOR//     if (MyLength > 0) {
  //REFACTOR//       for (Teuchos_Ordinal j=0; j<NumVectors; ++j) {
  //REFACTOR//         ptrs_[j] = nonContigValues_[j](0,MyLength).begin();
  //REFACTOR//       }
  //REFACTOR//     }
  //REFACTOR//   }
  }


} // namespace Tpetra

#endif // TPETRA_MULTIVECTORDATA_HPP

