// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
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

#include <Teuchos_Object.hpp>
#include <Teuchos_ArrayRCP.hpp>

namespace Tpetra {

  template<typename Ordinal, typename Scalar>
  class MultiVectorData : public Teuchos::Object {

  friend class MultiVector<Ordinal, Scalar>;
  friend class Vector<Ordinal, Scalar>;

  public:
    MultiVectorData();
    ~MultiVectorData();

  protected:
    void updateConstPointers();
    Teuchos::ArrayRCP<Scalar> values_;
    Teuchos::Array<Teuchos::ArrayView<Scalar> > ptrs_;
    Teuchos::Array<Teuchos::ArrayView<const Scalar> > cPtrs_;
    bool constantStride_;
    Ordinal stride_;
    
  private:
    //! Copy constructor (declared but not defined, do not use)
    MultiVectorData(const MultiVectorData<Ordinal,Scalar> &source);
    //! Assignment operator (declared but not defined, do not use)
    MultiVectorData<Ordinal,Scalar>& operator=(const MultiVectorData<Ordinal,Scalar> &source);

  }; // class MultiVectorData

} // namespace Tpetra

#endif // TPETRA_MULTIVECTORDATA_DECL_HPP

