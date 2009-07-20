//@HEADER
// ************************************************************************
// 
//          Kokkos: A Fast Kernel Package
//              Copyright (2004) Sandia Corporation
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
// ************************************************************************
//@HEADER

#ifndef KOKKOS_DEFAULTSPARSEMULTIPLY_H
#define KOKKOS_DEFAULTSPARSEMULTIPLY_H

#include <Teuchos_TestForException.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <stdexcept>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CrsMatrix.hpp" 
#include "Kokkos_MultiVector.hpp"

#ifndef KERNEL_PREFIX
  #define KERNEL_PREFIX
#endif

namespace Kokkos {

  template <class Scalar, class Ordinal, class Node>
  struct DefaultSparseMultiplyOp {
    typename Node::template buffer<const size_type>::buffer_t offsets;
    typename Node::template buffer<const Ordinal>::buffer_t inds;
    typename Node::template buffer<const  Scalar>::buffer_t x, vals;
    typename Node::template buffer<       Scalar>::buffer_t y;
    Scalar alpha, beta;

    inline KERNEL_PREFIX void execute(int i)
    {
      Scalar tmp = 0;
      for (size_type c=offsets[i]; c != offsets[i+1]; ++c) {
        tmp += vals[c] * x[inds[c]];
      }
      Scalar tmp2 = beta*y[i];
      y[i] = alpha*tmp + tmp2;
    }
  };


  template<class MAT, class MV>
  class DefaultSparseMultiply {
  public:
    typedef typename MAT::ScalarType  ScalarType;
    typedef typename MAT::OrdinalType OrdinalType;
    typedef typename MAT::NodeType    NodeType;

    //! @name Constructors/Destructor

    //@{

    //! DefaultSparseMultiply constuctor with variable number of indices per row.
    DefaultSparseMultiply(typename MAT::NodeType &node);

    //! DefaultSparseMultiply Destructor
    ~DefaultSparseMultiply();

    //@}

    //! @name Initialization of structure

    //@{

    //! Initialize structure of matrix
    int initializeStructure(const MAT& A, bool view);

    //! Initialize values of matrix
    int initializeValues(const MAT& A, bool view);

    //@}

    //! @name Computational methods

    //@{

    //! Applies the matrix to a MultiVector.
    int Apply(bool transpose, ScalarType alpha, const MV &X, ScalarType beta, MV &Y) const;

    //@}

  protected:
    // No data yet, because we have no implementations 
  };

  template <class MAT, class MV>
  DefaultSparseMultiply<MAT,MV>::DefaultSparseMultiply(typename MAT::NodeType &node)
  {}

  template <class MAT, class MV>
  DefaultSparseMultiply<MAT,MV>::~DefaultSparseMultiply()
  {}

  template <class MAT, class MV>
  int DefaultSparseMultiply<MAT,MV>::initializeStructure(const MAT& A, bool view)
  {
    TEST_FOR_EXCEPTION(true, std::logic_error, 
        "DefaultSparseMultiply::initializeStructure() not implemented for matrix type " << Teuchos::TypeNameTraits<MAT>::name() << std::endl;);
    return 0;
  }

  template <class MAT, class MV>
  int DefaultSparseMultiply<MAT,MV>::initializeValues(const MAT& A, bool view)
  {
    TEST_FOR_EXCEPTION(true, std::logic_error, 
        "DefaultSparseMultiply::initializeValues() not implemented for matrix type " << Teuchos::TypeNameTraits<MAT>::name() << std::endl;);
    return 0;
  }

  template <class MAT, class MV>
  int DefaultSparseMultiply<MAT,MV>::Apply(bool transpose, ScalarType alpha, const MV &X, ScalarType beta, MV &Y) const
  {
    TEST_FOR_EXCEPTION(true, std::logic_error, 
        "DefaultSparseMultiply::Apply() not implemented for matrix type " << Teuchos::TypeNameTraits<MAT>::name()
        << " and multivector type " << Teuchos::TypeNameTraits<MV>::name() << std::endl;);
    return 0;
  }


  // default implementation
  template<class Scalar, class Ordinal, class Node>
  class DefaultSparseMultiply<CrsMatrix<Scalar,Ordinal,Node>, MultiVector<Scalar,Ordinal,Node> > {
  public:
    typedef Scalar  ScalarType;
    typedef Ordinal OrdinalType;
    typedef Node    NodeType;

    //! @name Constructors/Destructor

    //@{

    //! DefaultSparseMultiply constuctor with variable number of indices per row.
    DefaultSparseMultiply(Node &node);

    //! DefaultSparseMultiply Destructor
    ~DefaultSparseMultiply();

    //@}

    //! @name Initialization of structure

    //@{

    //! Initialize structure of matrix
    int initializeStructure(const CrsMatrix<Scalar,Ordinal,Node> & A, bool view);

    //! Initialize values of matrix
    int initializeValues(const CrsMatrix<Scalar,Ordinal,Node> & A, bool view);

    //@}

    //! @name Computational methods

    //@{

    //! Applies the matrix to a MultiVector.
    inline int Apply(bool transpose, Scalar alpha, const MultiVector<Scalar,Ordinal,Node> &X, Scalar beta, MultiVector<Scalar,Ordinal,Node> &Y) const;

    //@}

  protected:
    //! Copy constructor (protected and unimplemented)
    DefaultSparseMultiply(const DefaultSparseMultiply& source);

    Node &node_;
    mutable DefaultSparseMultiplyOp<Scalar,Ordinal,Node> op_;

    inline void deleteStructure();
    inline void deleteValues();

    Ordinal numRows_;
    size_type numEntries_;

    bool storedStructure_, storedValues_;
  };

  template<class Scalar, class Ordinal, class Node>
  DefaultSparseMultiply<CrsMatrix<Scalar,Ordinal,Node>, MultiVector<Scalar,Ordinal,Node> >::DefaultSparseMultiply(Node &node)
  : node_(node), storedStructure_(false), storedValues_(false)
  {
    op_.alpha = Teuchos::ScalarTraits<Scalar>::one();
    op_.beta =  Teuchos::ScalarTraits<Scalar>::zero();
  }

  template<class Scalar, class Ordinal, class Node>
  DefaultSparseMultiply<CrsMatrix<Scalar,Ordinal,Node>, MultiVector<Scalar,Ordinal,Node> >::~DefaultSparseMultiply()
  {
    deleteStructure();
    deleteValues();
  }

  template<class Scalar, class Ordinal, class Node>
  void DefaultSparseMultiply<CrsMatrix<Scalar,Ordinal,Node>, MultiVector<Scalar,Ordinal,Node> >::deleteStructure() 
  {
    if (storedStructure_) {
      node_.template freeBuffer<const size_type>(op_.offsets); // FINISH: this can't be <const T>; what to do?
      node_.template freeBuffer<const Ordinal>(op_.inds);   // FINISH: this can't be <const T>; what to do?
    }
    storedStructure_ = false;
  }

  template<class Scalar, class Ordinal, class Node>
  void DefaultSparseMultiply<CrsMatrix<Scalar,Ordinal,Node>, MultiVector<Scalar,Ordinal,Node> >::deleteValues() 
  {
    if (storedValues_) {
      node_.template freeBuffer<const Scalar>(op_.vals);     // FINISH: this can't be <const T>; what to do?
    }
    storedValues_ = false;
  }

  template<class Scalar, class Ordinal, class Node>
  int DefaultSparseMultiply<CrsMatrix<Scalar,Ordinal,Node>, MultiVector<Scalar,Ordinal,Node> >::initializeStructure(const CrsMatrix<Scalar,Ordinal,Node> &A, bool view)
  {
    node_ = A.getNode();
    numRows_ = A.getNumRows();
    numEntries_ = A.getNumEntries();
    deleteStructure();
    if (view) {
      op_.offsets = A.const_offsets();
      op_.inds = A.const_indices();
    }
    else {
      typename NodeType::template buffer<const size_type>::buffer_t   srcoffsets;
      typename NodeType::template buffer<const Ordinal>::buffer_t srcindices;
      typename NodeType::template buffer<size_type>::buffer_t   nc_offsets;
      typename NodeType::template buffer<Ordinal>::buffer_t nc_indices;
      srcoffsets = A.const_offsets();
      srcindices = A.const_indices();
      nc_offsets = node_.template allocBuffer<size_type>(numRows_+1);
      nc_indices = node_.template allocBuffer<Ordinal>(numEntries_);
      node_.template copyBuffers<size_type>(numRows_+1,srcoffsets,0,nc_offsets,0);
      node_.template copyBuffers<Ordinal>(numEntries_,srcindices,0,nc_indices,0);
      op_.inds = nc_indices;
      op_.offsets = nc_offsets;
      storedStructure_ = true;
    }
    return 0;
  }

  template <class Scalar, class Ordinal, class Node>
  int DefaultSparseMultiply<CrsMatrix<Scalar,Ordinal,Node>, MultiVector<Scalar,Ordinal,Node> >::initializeValues(const CrsMatrix<Scalar,Ordinal,Node> &A, bool view)
  {
#ifdef HAVE_KOKKOS_DEBUG
    TEST_FOR_EXCEPTION(A.numRows() != numRows_ || A.numEntries() != numEntries_, std::runtime_error,
        Teuchos::typeName(*this) << "::initializevalues(A): structure of A does not match current structure.");
#endif
    deleteValues();
    if (view) {
      op_.vals = A.const_values();
    }
    else {
      typename NodeType::template buffer<const  Scalar>::buffer_t srcvalues;
      typename NodeType::template buffer<Scalar>::buffer_t nc_values;
      srcvalues = A.const_values();
      nc_values = node_.template allocBuffer<Scalar>(numEntries_);
      node_.template copyBuffers<Scalar>(numEntries_,srcvalues,0,nc_values,0);
      op_.vals = nc_values;
      storedValues_ = true;
    }
    return 0;
  }

  template <class Scalar, class Ordinal, class Node>
  inline int DefaultSparseMultiply<CrsMatrix<Scalar,Ordinal,Node>, MultiVector<Scalar,Ordinal,Node> >::Apply(
      bool transpose, 
      Scalar alpha, const MultiVector<Scalar,Ordinal,Node> &X, 
      Scalar beta, MultiVector<Scalar,Ordinal,Node> &Y) const 
  {
    TEST_FOR_EXCEPTION(transpose == true, std::runtime_error,
        Teuchos::typeName(*this) << "::Apply(): Operation does not currently support tranpose.");
#ifdef HAVE_KOKKOS_DEBUG
#endif
    op_.alpha = alpha;
    op_.beta  = beta;
    op_.x = X.getValues(0);
    op_.y = Y.getValues(0);
    node_.template parallel_for<DefaultSparseMultiplyOp<Scalar,Ordinal,Node> >(0,numRows_,op_);
    return 0;
  }

} // namespace Kokkos
#endif /* KOKKOS_DEFAULTSPARSEMULTIPLY_H */
