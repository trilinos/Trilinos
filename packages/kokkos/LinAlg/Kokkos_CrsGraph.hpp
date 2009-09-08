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

#ifndef KOKKOS_CRSGRAPH_H
#define KOKKOS_CRSGRAPH_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"

namespace Kokkos {

  //! Kokkos::CrsGraph: Kokkos compressed index sparse class.

  template <class Ordinal, class Node = DefaultNode::DefaultNodeType>  
  class CrsGraph {
  public:

  typedef Ordinal OrdinalType;
  typedef Node    NodeType;

    //! @name Constructors/Destructor
    //@{

    //! Default CrsGraph constuctor.
    CrsGraph(size_t numRows, const Teuchos::RCP<Node> &node = DefaultNode::getDefaultNode());

    //! CrsGraph Destructor
    ~CrsGraph();

    //@}

    //! @name Accessor routines.
    //@{ 
    
    //! Node accessor.
    Teuchos::RCP<Node> getNode() const;

    //@}

    //! @name Data entry and accessor methods.
    //@{

    //! Submit the indices and offset for 1D storage.
    void setPackedStructure(const Teuchos::ArrayRCP<const size_t> &offsets,
                            const Teuchos::ArrayRCP<const Ordinal> &allinds);

    //! Submit the indices and offset for 1D storage.
    //! Submit the indices for one row of 2D storage.
    void set2DStructure(size_t row, const Teuchos::ArrayRCP<const Ordinal> &rowinds);

    //! Retrieve the offsets for 1D storage.
    Teuchos::ArrayRCP<const size_t> get1DOffsets() const;

    //! Retrieve the indices for 1D storage.
    Teuchos::ArrayRCP<const Ordinal> get1DIndices() const;

    //! Retrieve the indices for one row of 2D storage.
    Teuchos::ArrayRCP<const Ordinal> get2DIndices(size_t row) const;

    //! Indicates whether or not the graph entries are packed.
    bool isPacked() const;
  
    //! Indicates that the graph is initialized, but empty.
    bool isEmpty() const;

    //! Release data associated with this graph.
    void clear();

    //! Return the number of rows in the graph.
    size_t getNumRows() const;

    //@}

  private:
    //! Copy constructor (protected and not implemented) 
    CrsGraph(const CrsGraph& sources);

    Teuchos::RCP<Node> node_;
    size_t numRows_;
    bool isInitialized_, isPacked_, isEmpty_;

    Teuchos::ArrayRCP<const Ordinal>                      pbuf_indices1D_;
    Teuchos::ArrayRCP<const size_t>                       pbuf_offsets_;
    Teuchos::ArrayRCP< Teuchos::ArrayRCP<const Ordinal> > pbuf_indices2D_;
  };


  //==============================================================================
  template <class Ordinal, class Node>
  CrsGraph<Ordinal,Node>::CrsGraph(size_t numRows, const Teuchos::RCP<Node> &node) 
  : node_(node)
  , numRows_(numRows)
  , isInitialized_(false)
  , isPacked_(false)
  , isEmpty_(true) {
  }

  //==============================================================================
  template <class Ordinal, class Node>
  CrsGraph<Ordinal,Node>::~CrsGraph() {
  }

  //==============================================================================
  template <class Ordinal, class Node>
  Teuchos::RCP<Node> CrsGraph<Ordinal,Node>::getNode() const {
    return node_;
  }

  //==============================================================================
  template <class Ordinal, class Node>
  void CrsGraph<Ordinal,Node>::clear() { 
    pbuf_indices1D_ = Teuchos::null;
    pbuf_indices2D_ = Teuchos::null;
    pbuf_offsets_   = Teuchos::null;
    isInitialized_ = false;
    isEmpty_       = true;
  }

  //==============================================================================
  template <class Ordinal, class Node>
  void CrsGraph<Ordinal,Node>::setPackedStructure(
                        const Teuchos::ArrayRCP<const size_t> &offsets,
                        const Teuchos::ArrayRCP<const Ordinal> &allinds) {
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(isInitialized_ == true, std::runtime_error,
        Teuchos::typeName(*this) << "::setPackedStructure(): graph is already initialized. Call clear() before reinitializing.");
#endif
    isEmpty_ = (allinds == Teuchos::null);
    pbuf_indices1D_ = allinds;
    pbuf_offsets_ = offsets;
    isInitialized_ = true;
    isPacked_ = true;
  }

  //==============================================================================
  template <class Ordinal, class Node>
  void CrsGraph<Ordinal,Node>::set2DStructure(
                              size_t row, 
                              const Teuchos::ArrayRCP<const Ordinal> &rowinds) {
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(isPacked_ == true, std::runtime_error,
        Teuchos::typeName(*this) << "::set2DStructure(): graph is already initialized with 1D structure. Call clear() before reinitializing.");
#endif
    if (isInitialized_ == false) {
      pbuf_indices2D_ = Teuchos::arcp<Teuchos::ArrayRCP<const Ordinal> >(numRows_);
      isInitialized_ = true;
    }
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION((row < 1 && row != 0) || row > numRows_, std::runtime_error,
        Teuchos::typeName(*this) << ":;set2DStructure(): specified row is invalid.");
#endif
    isEmpty_ = isEmpty_ && (rowinds == Teuchos::null);
    pbuf_indices2D_[row] = rowinds;
  }

  //==============================================================================
  template <class Ordinal, class Node>
  Teuchos::ArrayRCP<const Ordinal> 
  CrsGraph<Ordinal,Node>::get1DIndices() const {
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(isPacked_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::get1DIndices(): graph is uninitialized or not packed.");
#endif
    return pbuf_indices1D_;
  }

  //==============================================================================
  template <class Ordinal, class Node>
  Teuchos::ArrayRCP<const Ordinal> 
  CrsGraph<Ordinal,Node>::get2DIndices(size_t row) const {
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(isInitialized_ == false || isPacked_ == true, std::runtime_error,
        Teuchos::typeName(*this) << "::get2DIndices(): graph is uninitialized or initialized packed.");
    TEST_FOR_EXCEPTION((row < 1 && row != 0) || row > numRows_, std::runtime_error,
        Teuchos::typeName(*this) << "::get2DIndices(): row number is invalid.");
#endif
    return pbuf_indices2D_[row];
  }

  //==============================================================================
  template <class Ordinal, class Node>
  bool CrsGraph<Ordinal,Node>::isPacked() const {
    return isPacked_;
  }

  //==============================================================================
  template <class Ordinal, class Node>
  bool CrsGraph<Ordinal,Node>::isEmpty() const {
    return isEmpty_;
  }

  //==============================================================================
  template <class Ordinal, class Node>
  Teuchos::ArrayRCP<const size_t> 
  CrsGraph<Ordinal,Node>::get1DOffsets() const {
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(isPacked_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::get1DOffsets(): graph is uninitialized or not packed.");
#endif
    return pbuf_offsets_;
  }

  //==============================================================================
  template <class Ordinal, class Node>
  size_t CrsGraph<Ordinal,Node>::getNumRows() const {
    return numRows_;
  }

} // namespace Kokkos

#endif /* KOKKOS_CRSGRAPH_H */
