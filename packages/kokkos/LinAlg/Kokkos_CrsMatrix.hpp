//@HEADER
// ************************************************************************
// 
//                Kokkos: A Fast Kernel Package
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

#ifndef KOKKOS_CRSMATRIX_H
#define KOKKOS_CRSMATRIX_H

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"

namespace Kokkos {

//! Kokkos::CrsMatrix: Kokkos compressed index sparse matrix base class.

/*! The Kokkos::CrsMatrix implements the Kokkos::CisMatrix interface 
    using either a Harwell-Boeing matrix or generalized form of one.

*/    

  template<class Scalar, class Ordinal, class Node = DefaultNode::DefaultNodeType>
  class CrsMatrix {
  public:
    typedef Scalar  ScalarType;
    typedef Ordinal OrdinalType;
    typedef Node    NodeType;

    //! @name Constructors/Destructor

    //@{

    //! Default CrsMatrix constuctor.
    CrsMatrix(Node &node = DefaultNode::getDefaultNode());

    //! CrsMatrix Destructor
    ~CrsMatrix();

    //@}

    //! @name Harwell-Boeing Format Initialization Methods

    //@{

    //! Initialize structure of matrix
    void initializeStructure(int N, Ordinal nnzEachRow);

    //! Initialize structure of matrix
    void initializeStructure(int N, const Ordinal *nnzPerRow);

    //@}

    //! @name Matrix entry methods

    //@{

    int insertEntries(int row, int numEntries, const Ordinal *indices, const Scalar *values);

    //@}

    //! @name Matrix Attribute access methods

    //@{

    //! Compute node 
    Node & getNode() const;

    //! Number of rows
    Ordinal getNumRows() const;

    //! Number of matrix entries
    Ordinal getNumEntries() const;

    //! Compute buffers (non-const)
    void data(typename Node::template buffer<Ordinal>::buffer_t &offsets,
              typename Node::template buffer<Ordinal>::buffer_t &indices,
              typename Node::template buffer<Scalar>::buffer_t &values);

    //! Compute buffers (const)
    void data(typename Node::template buffer<const Ordinal>::buffer_t &offsets,
              typename Node::template buffer<const Ordinal>::buffer_t &indices,
              typename Node::template buffer<const Scalar>::buffer_t &values) const;

    //! Offset compute buffer (const)
    typename Node::template buffer<const Ordinal>::buffer_t const_offsets() const;

    //! Offset compute buffer (non-const)
    typename Node::template buffer<Ordinal>::buffer_t offsets();

    //! Index compute buffer (const)
    typename Node::template buffer<const Ordinal>::buffer_t const_indices() const;

    //! Index compute buffer (non-const)
    typename Node::template buffer<Ordinal>::buffer_t indices();

    //! Values compute buffer (const)
    typename Node::template buffer<const Scalar>::buffer_t const_values() const;

    //! Values compute buffer (non-const)
    typename Node::template buffer<Scalar>::buffer_t values();

    void Print(std::ostream &out);

    //@}

  protected:

    //! Copy constructor (protected and not implemented)
    CrsMatrix(const CrsMatrix& source);

    Node &node_;
    int numRows_, numEntries_;
    typename Node::template buffer<Ordinal>::buffer_t offsets_, indices_;
    typename Node::template buffer<Scalar>::buffer_t values_;
  };


  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  CrsMatrix<Scalar,Ordinal,Node>::CrsMatrix(Node &node)
  : node_(node), numRows_(0), numEntries_(0) {}

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  CrsMatrix<Scalar,Ordinal,Node>::~CrsMatrix()
  {
    if (numRows_ > 0) {
      node_.template freeBuffer<Scalar>(values_);
      node_.template freeBuffer<Ordinal>(indices_);
      node_.template freeBuffer<Ordinal>(offsets_);
    }
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  void CrsMatrix<Scalar,Ordinal,Node>::initializeStructure(int N, Ordinal nnzEachRow)
  {
    numRows_ = N;
    if (numRows_ > 0) {
      numEntries_ = 0;
      offsets_ = node_.template allocBuffer<Ordinal>(numRows_+1);
      Ordinal *h_offsets = node_.template viewBuffer<Ordinal>(true,numRows_+1,offsets_,0);
      numEntries_ = nnzEachRow*N;
      h_offsets[0] = 0;
      for (int i=1; i<numRows_; ++i) {
        h_offsets[i] = h_offsets[i-1]+nnzEachRow;
      }
      h_offsets[numRows_] = numEntries_;
      node_.template releaseView<Ordinal>(h_offsets); h_offsets = NULL;
      values_  = node_.template allocBuffer<Scalar>(numEntries_);
      indices_ = node_.template allocBuffer<Ordinal>(numEntries_);
    }
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  void CrsMatrix<Scalar,Ordinal,Node>::initializeStructure(int N, const Ordinal *nnzPerRow)
  {
    numRows_ = N;
    if (numRows_ > 0) {
      numEntries_ = 0;
      offsets_ = node_.template allocBuffer<Ordinal>(numRows_+1);
      Ordinal *h_offsets = node_.template viewBuffer<Ordinal>(true,numRows_+1,offsets_,0);
      for (int i=0; i<numRows_; ++i) {
        h_offsets[i] = numEntries_;
        numEntries_ += nnzPerRow[i];
      }
      h_offsets[numRows_] = numEntries_;
      node_.template releaseView<Ordinal>(h_offsets); h_offsets = NULL;
      values_  = node_.template allocBuffer<Scalar>(numEntries_);
      indices_ = node_.template allocBuffer<Ordinal>(numEntries_);
    }
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  int CrsMatrix<Scalar,Ordinal,Node>::insertEntries(int row, int numEntries, const Ordinal *indices, const Scalar *values)
  {
    if (row < 0 || row >= numRows_) return -1;
    const Ordinal *h_offsets = node_.template viewBufferConst<Ordinal>(2,offsets_,row);
    const int rowNNZ = h_offsets[1]-h_offsets[0];
    if (numEntries > rowNNZ) return -2;
    Ordinal *h_indices = node_.template viewBuffer<Ordinal>(true,rowNNZ,indices_,h_offsets[0]);
    Scalar  *h_values  = node_.template viewBuffer<Scalar >(true,rowNNZ,values_, h_offsets[0]);
    int e = 0;
    while (e != numEntries) {
      h_indices[e] = indices[e];
      h_values[e] = values[e];
      ++e;
    }
    while (e != rowNNZ) {
      h_indices[e] = 0;
      h_values[e] = 0;
      ++e;
    }
    node_.template releaseView<Ordinal>(h_offsets); h_offsets = NULL;
    node_.template releaseView<Ordinal>(h_indices); h_indices = NULL;
    node_.template releaseView<Scalar >(h_values);  h_values  = NULL;
    return 0;
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  Node & CrsMatrix<Scalar,Ordinal,Node>::getNode() const
  { return node_; }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  Ordinal CrsMatrix<Scalar,Ordinal,Node>::getNumRows() const
  { return numRows_; }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  Ordinal CrsMatrix<Scalar,Ordinal,Node>::getNumEntries() const
  { return numEntries_; }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  void CrsMatrix<Scalar,Ordinal,Node>::data(typename Node::template buffer<Ordinal>::buffer_t &offsets,
                       typename Node::template buffer<Ordinal>::buffer_t &indices,
                       typename Node::template buffer<Scalar>::buffer_t &values)
  {
    offsets = offsets_;
    indices = indices_;
    values  = values_;
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  void CrsMatrix<Scalar,Ordinal,Node>::data(typename Node::template buffer<const Ordinal>::buffer_t &offsets,
                       typename Node::template buffer<const Ordinal>::buffer_t &indices,
                       typename Node::template buffer<const Scalar>::buffer_t &values) const
  {
    offsets = offsets_;
    indices = indices_;
    values  = values_;
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  typename Node::template buffer<const Ordinal>::buffer_t 
  CrsMatrix<Scalar,Ordinal,Node>::const_offsets() const
  {
    return offsets_;
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  typename Node::template buffer<Ordinal>::buffer_t 
  CrsMatrix<Scalar,Ordinal,Node>::offsets()
  {
    return offsets_;
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  typename Node::template buffer<const Ordinal>::buffer_t 
  CrsMatrix<Scalar,Ordinal,Node>::const_indices() const
  {
    return indices_;
  }


  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  typename Node::template buffer<Ordinal>::buffer_t 
  CrsMatrix<Scalar,Ordinal,Node>::indices()
  {
    return indices_;
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  typename Node::template buffer<const Scalar>::buffer_t 
  CrsMatrix<Scalar,Ordinal,Node>::const_values() const
  {
    return values_;
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  typename Node::template buffer<Scalar>::buffer_t 
  CrsMatrix<Scalar,Ordinal,Node>::values()
  {
    return values_;
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  void CrsMatrix<Scalar,Ordinal,Node>::Print(std::ostream &out)
  {
    using std::endl;
    const Ordinal *h_offsets = node_.template viewBufferConst<Ordinal>(numRows_+1,offsets_,0);
    const Ordinal *h_indices = node_.template viewBufferConst<Ordinal>(numEntries_,indices_,0);
    const Scalar  *h_values  = node_.template viewBufferConst<Scalar >(numEntries_,values_ ,0);
    out << "Matrix data: " << endl;
    for (int i=0; i<numRows_; ++i) {
      for (int j=h_offsets[i]; j!=h_offsets[i+1]; ++j) {
        out << "(" << i << "," << h_indices[j] << ") : " << h_values[j] << endl;
      }
    }
    node_.template releaseView<Ordinal>(h_offsets); h_offsets = NULL;
    node_.template releaseView<Ordinal>(h_indices); h_indices = NULL;
    node_.template releaseView<Scalar >(h_values ); h_values  = NULL;
  }

} // namespace Kokkos


#endif /* KOKKOS_CRSMATRIX_H */
