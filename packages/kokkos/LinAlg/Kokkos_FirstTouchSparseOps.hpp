//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2009) Sandia Corporation
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

#ifndef KOKKOS_FIRSTTOUCHSPARSEOPS_HPP
#define KOKKOS_FIRSTTOUCHSPARSEOPS_HPP

#include "Kokkos_DefaultSparseOps.hpp"

namespace Kokkos {

  //=========================================================================================================================
  // 
  // A first-touch sparse ops
  // 
  //=========================================================================================================================

  /** \brief Sparse matrix-vector multiplication and solve routines with first-touch allocation of matrix data.
      \ingroup kokkos_crs_ops

      Effectively, this class is identical to DefaultSparseOps, except for the allocRowPtrs() and allocStorage() methods, which are 
      modified to use a first-touch allocation. The relationship is done using private inheritence, though this class should never be 
      accessed via a reference to its parent class. (Tpetra is unaware this relationship, and so will never do that.)
   */
  template <class Scalar, class Ordinal, class Node>
  class FirstTouchSparseOps : private DefaultHostSparseOps<Scalar,Ordinal,Node> {
  public:
    //@{ 
    //! @name Typedefs and structs

    //! The type of the individual entries of the sparse matrix.
    typedef Scalar  scalar_type;
    //! The type of the (local) indices describing the structure of the sparse matrix.
    typedef Ordinal ordinal_type;
    //! The Kokkos Node type.
    typedef Node    node_type;
    //! The type of this object, the sparse operator object
    typedef FirstTouchSparseOps<Scalar,Ordinal,Node> sparse_ops_type;
    
    using DefaultHostSparseOps<Scalar,Ordinal,Node>::graph;
    using DefaultHostSparseOps<Scalar,Ordinal,Node>::matrix;

    /** \brief Rebind struct, for specifying type information for a different scalar.
          
        This specifies a DefaultHostSparseOps object, regardless of scalar type.
      */
    template <class S2>
    struct bind_scalar {
      typedef FirstTouchSparseOps<S2,Ordinal,Node> other_type;
    };

    //@}
    //! @name Constructors/Destructor
    //@{

    //! \brief Constructor accepting and retaining a node object.
    FirstTouchSparseOps(const RCP<Node> &node);

    //! Destructor
    virtual ~FirstTouchSparseOps();

    //@}
    //! @name Overloaded from DefaultHostSparseOps for first-touch allocation.
    //@{

    //! \brief Allocate and initialize the storage for the matrix values.
    static ArrayRCP<size_t> allocRowPtrs(const RCP<Node> &node, const ArrayView<const size_t> &numEntriesPerRow);

    //! \brief Allocate and initialize the storage for a sparse graph.
    template <class T>
    static ArrayRCP<T> allocStorage(const RCP<Node> &node, const ArrayView<const size_t> &rowPtrs);

    using DefaultHostSparseOps<Scalar,Ordinal,Node>::finalizeGraph;
    using DefaultHostSparseOps<Scalar,Ordinal,Node>::finalizeMatrix;
    using DefaultHostSparseOps<Scalar,Ordinal,Node>::finalizeGraphAndMatrix;
    using DefaultHostSparseOps<Scalar,Ordinal,Node>::setGraphAndMatrix;
    using DefaultHostSparseOps<Scalar,Ordinal,Node>::multiply;
    using DefaultHostSparseOps<Scalar,Ordinal,Node>::solve;
    using DefaultHostSparseOps<Scalar,Ordinal,Node>::getNode;

    //@}

  protected:
    //! Copy constructor (protected and unimplemented)
    FirstTouchSparseOps(const DefaultHostSparseOps<Scalar,Ordinal,Node>& source);
  };

  namespace FirstTouchSparseOpsDetails {

    struct rowPtrsInitKernel {
      size_t       *rowPtrs;
      inline KERNEL_PREFIX void execute(int i) const {
        rowPtrs[i] = 0;
      }
    };

    template <class T>
    struct valsInitKernel {
      const size_t *rowPtrs;
      T * entries;
      inline KERNEL_PREFIX void execute(int i) const {
        size_t *beg = entries+rowPtrs[i],
               *end = entries+rowPtrs[i+1];
        while (beg != end) {
          *beg = Teuchos::ScalarTraits<T>::zero();
        }
      }
    };
  }

  // ======= pointer allocation ===========
  template <class Scalar, class Ordinal, class Node>
  ArrayRCP<size_t>
  FirstTouchSparseOps<Scalar,Ordinal,Node>::allocRowPtrs(const RCP<Node> &node, const ArrayView<const size_t> &numEntriesPerRow)
  {
    const size_t numrows = numEntriesPerRow.size();
    FirstTouchSparseOpsDetails::rowPtrsInitKernel kern;
    // allocate
    kern.rowPtrs = new size_t[numrows+1];
    // parallel first touch
    node->parallel_for(0,numrows+1,kern);
    // encapsulate
    ArrayRCP<size_t> ptrs = arcp<size_t>(kern.rowPtrs,0,numrows+1,true);
    // compute in serial. parallelize later, perhaps; it's only O(N)
    ptrs[0] = 0;
    std::partial_sum( numEntriesPerRow.getRawPtr(), numEntriesPerRow.getRawPtr()+numEntriesPerRow.size(), ptrs.begin()+1 );
    return ptrs;
  }

  // ======= other allocation ===========
  template <class Scalar, class Ordinal, class Node>
  template <class T>
  ArrayRCP<T>
  FirstTouchSparseOps<Scalar,Ordinal,Node>::allocStorage(const RCP<Node> &node, const ArrayView<const size_t> &rowPtrs)
  {
    const size_t totalNumEntries = *(rowPtrs.end()-1);
    const size_t numRows = rowPtrs.size() - 1;
    FirstTouchSparseOpsDetails::valsInitKernel<T> kern;
    ArrayRCP<T> vals;
    if (totalNumEntries > 0) {
      // allocate
      kern.entries = new T[totalNumEntries];
      kern.rowPtrs = rowPtrs.getRawPtr();
      // first touch
      node->parallel_for(0,numRows,kern);
      // encapsulate
      vals = arcp<T>(kern.entries,0,totalNumEntries,true);
    }
    return vals;
  }

  template <class Scalar, class Ordinal, class Node>
  FirstTouchSparseOps<Scalar,Ordinal,Node>::FirstTouchSparseOps(const RCP<Node> &node) 
  : DefaultHostSparseOps<Scalar,Ordinal,Node>(node) 
  {}

  template <class Scalar, class Ordinal, class Node>
  FirstTouchSparseOps<Scalar,Ordinal,Node>::~FirstTouchSparseOps()
  {}

} // end namespace Kokkos

#endif // KOKKOS_FIRSTTOUCH_SPARSEOPS_HPP
