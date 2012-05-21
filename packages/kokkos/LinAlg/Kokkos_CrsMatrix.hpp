//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef KOKKOS_CRSMATRIX_HPP
#define KOKKOS_CRSMATRIX_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_Assert.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_CrsGraph.hpp"

namespace Kokkos {

  /*! @class CrsMatrix
      @brief A compressed-row sparse matrix that lives in host memory.

      Matrix coefficients are stored in one of the following ways:
      - The matrix data is packed into a single contiguous array. The array is of length num-total-nonzeros.  Two auxiliary arrays
       give starting/ending indices for each row and column indices, respectively.  This is essentially the standard compressed sparse row (CSR) format.
      This is sometimes referred to as "1D" storage elsewhere in Kokkos.
      - Data is packed into an array of arrays.  The outer array is length number-of-rows, and each row is an array of length num-nonzeros-per-row.
      This is sometimes referred to as "2D" storage.
      @ingroup kokkos_crs_ops
  */
  template <class Scalar,
            class Ordinal,
            class Node>
  class CrsMatrix {
  public:

    typedef Scalar   ScalarType;
    typedef Ordinal  OrdinalType;
    typedef Node     NodeType;

    //! @name Constructors/Destructor
    //@{

    //! Default constructor with no graph. (Must be set later.)
    CrsMatrix();

    //! Constructor with a matrix-owned non-const graph
    CrsMatrix(CrsGraphHostCompute<Ordinal,Node> &graph);

    //! Constructor with a non-owned const graph.
    CrsMatrix(const CrsGraphHostCompute<Ordinal,Node> &graph);

    //! CrsMatrix Destructor
    virtual ~CrsMatrix();

    //@}

    //! @name Graph set routines.
    //@{

    //! Set matrix-owned graph.
    void setOwnedGraph(CrsGraphHostCompute<Ordinal,Node> &graph);

    //! Set static graph.
    void setStaticGraph(const CrsGraphHostCompute<Ordinal,Node> &graph);
    
    //@}

    //! @name Data entry and accessor methods.
    //@{

    //! Node accessor.
    RCP<Node> getNode() const;

    //! Return the number of rows in the matrix.
    size_t getNumRows() const;

    //! Return the number of entries in the matrix.
    size_t getNumEntries() const;

    //! Indicates that the graph is filled, but empty.
    bool isEmpty() const;

    //! Whether the graph has been finalized.
    bool isFinalized() const;

    /// \brief Whether the structure is 1D.
    ///
    /// It will never be the case that both \c is1DStructure() and
    /// \c is2DStructure() return true.
    bool is1DStructure() const;

    /// \brief Whether the structure is 2D.
    ///
    /// It will never be the case that both \c is1DStructure() and \c
    /// is2DStructure() return true.
    bool is2DStructure() const;

    //! Whether the sparse matrix stucture is optimized.
    bool isOptimized() const;

    //! Submit the values for 1D storage.
    /**
          \post is1DStructure() == true
     */
    void set1DValues(ArrayRCP<Scalar> vals);

    //! Submit the values for 2D storage.
    /**
          \post is2DStructure() == true
     */
    void set2DValues(ArrayRCP<ArrayRCP<Scalar> > vals);

    //! Retrieve the values for 1D storage.
    /**
          If is1DStructure() == false, then
          \post vals == null
     */
    void get1DValues(ArrayRCP<Scalar> &vals);

    //! Retrieve the structure for 2D storage.
    /**
          If is2DStructure() == false, then
          \post vals == null
     */
    void get2DValues(ArrayRCP<ArrayRCP<Scalar> > &inds);

    //! Instruct the matrix to perform any necessary manipulation, including (optionally) optimizing the storage of the matrix data.
    /**
          If the matrix is associated with a matrix-owned graph, then it will be finalized as well. A static graph will not be modified.

          @param[in] OptimizeStorage   Permit the graph to reallocate storage on the host in order to provide optimal storage and/or performance.
          \post if OptimizeStorage == true, then is2DStructure() == true
     */
    void finalize(bool OptimizeStorage);

    //! Release data associated with this graph.
    virtual void clear();

    //@}

  protected:
    //! Copy constructor (protected and not implemented)
    CrsMatrix(const CrsMatrix& sources);

    CrsGraphHostCompute<Ordinal,Node> *myGraph_;
    const CrsGraphHostCompute<Ordinal,Node> *staticGraph_;
    bool isFinalized_;

    // 2D storage
    ArrayRCP<ArrayRCP<Scalar> >  values2D_;
    // 1D storage
    ArrayRCP<Scalar>             values1D_;
  };


  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  CrsMatrix<Scalar,Ordinal,Node>::CrsMatrix()
  {
    myGraph_ = NULL;
    staticGraph_ = NULL;
    CrsMatrix<Scalar,Ordinal,Node>::clear();
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  CrsMatrix<Scalar,Ordinal,Node>::CrsMatrix(const CrsGraphHostCompute<Ordinal,Node> &graph)
  {
    setStaticGraph(graph);
    CrsMatrix<Scalar,Ordinal,Node>::clear();
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  CrsMatrix<Scalar,Ordinal,Node>::CrsMatrix(CrsGraphHostCompute<Ordinal,Node> &graph)
  {
    setOwnedGraph(graph);
    CrsMatrix<Scalar,Ordinal,Node>::clear();
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  void CrsMatrix<Scalar,Ordinal,Node>::setStaticGraph(const CrsGraphHostCompute<Ordinal,Node> &hgraph)
  {
    myGraph_ = NULL;
    staticGraph_ = &hgraph;
    TEUCHOS_TEST_FOR_EXCEPTION( ! staticGraph_->isFinalized (),
      std::runtime_error, Teuchos::typeName(*this) << ": construction with "
      "static graph requires that the graph is already finalized.");
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  void CrsMatrix<Scalar,Ordinal,Node>::setOwnedGraph(CrsGraphHostCompute<Ordinal,Node> &hgraph)
  {
    myGraph_ = &hgraph;
    staticGraph_ = &hgraph;
    TEUCHOS_TEST_FOR_EXCEPTION( myGraph_->isFinalized() == true, std::runtime_error,
        Teuchos::typeName(*this) << ": construction with matrix-owned graph requires that the graph is not yet finalized.");
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  CrsMatrix<Scalar,Ordinal,Node>::~CrsMatrix() {
  }

  // ======= clear ===========
  template <class Scalar, class Ordinal, class Node>
  void CrsMatrix<Scalar,Ordinal,Node>::clear() {
    isFinalized_   = false;
    values1D_      = null;
    values2D_      = null;
  }

  // ======= node ===========
  template <class Scalar, class Ordinal, class Node>
  RCP<Node> CrsMatrix<Scalar,Ordinal,Node>::getNode() const {
    return staticGraph_->getNode();
  }

  // ======= numrows ===========
  template <class Scalar, class Ordinal, class Node>
  size_t CrsMatrix<Scalar,Ordinal,Node>::getNumRows() const {
    return staticGraph_->getNumRows();
  }

  // ======= numentries ===========
  template <class Scalar, class Ordinal, class Node>
  size_t CrsMatrix<Scalar,Ordinal,Node>::getNumEntries() const {
    return staticGraph_->getNumEntries();
  }

  // ======= isempty ===========
  template <class Scalar, class Ordinal, class Node>
  bool CrsMatrix<Scalar,Ordinal,Node>::isEmpty() const {
    return staticGraph_->isEmpty();
  }

  // ======= isfinalized ===========
  template <class Scalar, class Ordinal, class Node>
  bool CrsMatrix<Scalar,Ordinal,Node>::isFinalized() const {
    return isFinalized_;
  }

  // ======= is1d ===========
  template <class Scalar, class Ordinal, class Node>
  bool CrsMatrix<Scalar,Ordinal,Node>::is1DStructure() const {
    return staticGraph_->is1DStructure();
  }

  // ======= is2d ===========
  template <class Scalar, class Ordinal, class Node>
  bool CrsMatrix<Scalar,Ordinal,Node>::is2DStructure() const {
    return staticGraph_->is2DStructure();
  }

  // ======= isopt ===========
  template <class Scalar, class Ordinal, class Node>
  bool CrsMatrix<Scalar,Ordinal,Node>::isOptimized() const {
    return staticGraph_->isOptimized();
  }

  // ======= get 1d ===========
  template <class Scalar, class Ordinal, class Node>
  void CrsMatrix<Scalar,Ordinal,Node>::get1DValues(ArrayRCP<Scalar> &vals)
  {
    vals = values1D_;
  }

  // ======= get 2d ===========
  template <class Scalar, class Ordinal, class Node>
  void CrsMatrix<Scalar,Ordinal,Node>::get2DValues(ArrayRCP<ArrayRCP<Scalar> > &vals)
  {
    vals = values2D_;
  }

  // ======= set 1d ===========
  template <class Scalar, class Ordinal, class Node>
  void CrsMatrix<Scalar,Ordinal,Node>::set1DValues(ArrayRCP<Scalar> vals)
  {
    this->clear();
    TEUCHOS_TEST_FOR_EXCEPTION(! is1DStructure (), std::runtime_error,
      Teuchos::typeName(*this) << "::set1DValues(vals): graph must have 1D "
      "structure and must be set before matrix.");
    values1D_ = vals;
    TEUCHOS_TEST_FOR_EXCEPTION((size_t)values1D_.size() < getNumEntries(),
      std::runtime_error, Teuchos::typeName(*this) << "::set1DValues(vals): "
      "inds must have as many entries as the number of entries in the "
      "associated graph, and the graph entries must be set first.");
  }

  // ======= set 2d ===========
  template <class Scalar, class Ordinal, class Node>
  void CrsMatrix<Scalar,Ordinal,Node>::set2DValues(ArrayRCP<ArrayRCP<Scalar> > vals)
  {
    TEUCHOS_TEST_FOR_EXCEPTION( is2DStructure() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::set2DValues(vals): graph must have 2D structure and must be set before matrix.");
    TEUCHOS_TEST_FOR_EXCEPTION( (size_t)vals.size() != getNumRows(), std::runtime_error,
        Teuchos::typeName(*this) << "::set2DValues(inds): vals must have as many entries as the number of rows in the associated graph.");
    this->clear();
    values2D_  = vals;
  }


  // ======= finalize ===========
  template <class Scalar, class Ordinal, class Node>
  void CrsMatrix<Scalar,Ordinal,Node>::finalize(bool OptimizeStorage)
  {
    if (isFinalized_ && ! (OptimizeStorage && ! staticGraph_->isOptimized ())) {
      // If we've already finalized, and if we don't need to optimize
      // the graph's storage, then we don't have to do anything.
      return;
    }
    // we only have something to do if we have a matrix-owned graph
    if (myGraph_ != NULL) {
      myGraph_->template finalize<Scalar> (OptimizeStorage, values2D_, values1D_);
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        OptimizeStorage && ! staticGraph_->isOptimized (),
        std::runtime_error, Teuchos::typeName(*this) << "::finalize(Optimize"
        "Storage == true): underlying static graph is not optimized and cannot "
        "be optimized.");
    }
    isFinalized_ = true;
  }

} // namespace Kokkos

#endif /* KOKKOS_CRSMATRIX_HPP */
