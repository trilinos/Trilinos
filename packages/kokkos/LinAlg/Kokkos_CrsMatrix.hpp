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

#ifndef KOKKOS_CRSMATRIX_HPP
#define KOKKOS_CRSMATRIX_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_CrsGraph.hpp"

namespace Kokkos {

  //=========================================================================================================================
  // 
  // A host-resident CrsMatrix
  // 
  //=========================================================================================================================

  /** \brief A default host-compute compressed-row sparse matrix.
      \ingroup kokkos_crs_ops
   */
  template <class Scalar, 
            class Ordinal, 
            class Node,
            class LocalMatOps>
  class CrsMatrixHostCompute {
  public:

    typedef Scalar                ScalarType;
    typedef Ordinal               OrdinalType;
    typedef Node                  NodeType;
    typedef LocalMatOps           LocalMatOpsType;

    //! @name Constructors/Destructor
    //@{

    //! Default constructor with no graph. (Must be set later.) 
    CrsMatrixHostCompute();

    //! Constructor with a matrix-owned non-const graph
    CrsMatrixHostCompute(CrsGraphHostCompute<Ordinal,Node,LocalMatOps> &graph);

    //! Constructor with a non-owned const graph.
    CrsMatrixHostCompute(const CrsGraphHostCompute<Ordinal,Node,LocalMatOps> &graph);

    //! CrsMatrixHostCompute Destructor
    virtual ~CrsMatrixHostCompute();

    //@}

    //! @name Graph set routines.
    //@{

    //! Set matrix-owned graph.
    void setOwnedGraph(CrsGraphHostCompute<Ordinal,Node,LocalMatOps> &graph);

    //! Set static graph.
    void setStaticGraph(const CrsGraphHostCompute<Ordinal,Node,LocalMatOps> &graph);
    
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

    //! Indicatest that the graph has been finalized.
    bool isFinalized() const;

    //! \brief Indicate that the structure is 1D.
    //! It will never be the case that both is1DStructure() and is2DStructure() return true.
    bool is1DStructure() const;

    //! \brief Indicate that the structure is 2D.
    //! It will never be the case that both is1DStructure() and is2DStructure() return true.
    bool is2DStructure() const;

    //! \brief Indicate that the stucture is optimized.
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
    CrsMatrixHostCompute(const CrsMatrixHostCompute& sources);

    CrsGraphHostCompute<Ordinal,Node,LocalMatOps> *myGraph_;
    const CrsGraphHostCompute<Ordinal,Node,LocalMatOps> *staticGraph_;
    bool isFinalized_;

    // 2D storage
    ArrayRCP<ArrayRCP<Scalar> >  values2D_;
    // 1D storage
    ArrayRCP<Scalar>             values1D_;
  };


  //==============================================================================
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::CrsMatrixHostCompute()
  {
    myGraph_ = NULL;
    staticGraph_ = NULL;
    CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::clear();
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::CrsMatrixHostCompute(const CrsGraphHostCompute<Ordinal,Node,LocalMatOps> &graph)
  {
    setStaticGraph(graph);
    CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::clear();
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::CrsMatrixHostCompute(CrsGraphHostCompute<Ordinal,Node,LocalMatOps> &graph)
  {
    setOwnedGraph(graph);
    CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::clear();
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  void CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::setStaticGraph(const CrsGraphHostCompute<Ordinal,Node,LocalMatOps> &hgraph)
  {
    myGraph_ = NULL;
    staticGraph_ = &hgraph;
    TEST_FOR_EXCEPTION( staticGraph_->isFinalized() == false, std::runtime_error, 
        Teuchos::typeName(*this) << ": construction with static graph requires that the graph is already finalized.");
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  void CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::setOwnedGraph(CrsGraphHostCompute<Ordinal,Node,LocalMatOps> &hgraph)
  {
    myGraph_ = &hgraph;
    staticGraph_ = &hgraph;
    TEST_FOR_EXCEPTION( myGraph_->isFinalized() == true, std::runtime_error,
        Teuchos::typeName(*this) << ": construction with matrix-owned graph requires that the graph is not yet finalized.");
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::~CrsMatrixHostCompute() {
  }

  // ======= clear ===========
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  void CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::clear() {
    isFinalized_   = false;
    values1D_      = null;
    values2D_      = null;
  }

  // ======= node ===========
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  RCP<Node> CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::getNode() const {
    return staticGraph_->getNode();
  }

  // ======= numrows ===========
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  size_t CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::getNumRows() const {
    return staticGraph_->getNumRows();
  }

  // ======= numentries ===========
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  size_t CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::getNumEntries() const {
    return staticGraph_->getNumEntries();
  }

  // ======= isempty ===========
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  bool CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::isEmpty() const {
    return staticGraph_->isEmpty();
  }

  // ======= isfinalized ===========
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  bool CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::isFinalized() const {
    return isFinalized_;
  }

  // ======= is1d ===========
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  bool CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::is1DStructure() const {
    return staticGraph_->is1DStructure();
  }

  // ======= is2d ===========
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  bool CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::is2DStructure() const {
    return staticGraph_->is2DStructure();
  }

  // ======= isopt ===========
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  bool CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::isOptimized() const {
    return staticGraph_->isOptimized();
  }

  // ======= get 1d ===========
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  void CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::get1DValues(ArrayRCP<Scalar> &vals)
  {
    vals = values1D_;
  }

  // ======= get 2d ===========
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  void CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::get2DValues(ArrayRCP<ArrayRCP<Scalar> > &vals)
  {
    vals = values2D_;
  }

  // ======= set 1d ===========
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  void CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::set1DValues(ArrayRCP<Scalar> vals)
  {
    this->clear();
    TEST_FOR_EXCEPTION( is1DStructure() == false, std::runtime_error, 
        Teuchos::typeName(*this) << "::set1DValues(vals): graph must have 1D structure and must be set before matrix.");
    values1D_ = vals;
    TEST_FOR_EXCEPTION( (size_t)values1D_.size() < getNumEntries(), std::runtime_error,
        Teuchos::typeName(*this) << "::set1DValues(vals): inds must have as many entries as the number of entries in the associated graph, and the graph entries must be set first.");
  }

  // ======= set 2d ===========
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  void CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::set2DValues(ArrayRCP<ArrayRCP<Scalar> > vals)
  {
    TEST_FOR_EXCEPTION( is2DStructure() == false, std::runtime_error, 
        Teuchos::typeName(*this) << "::set2DValues(vals): graph must have 2D structure and must be set before matrix.");
    TEST_FOR_EXCEPTION( (size_t)vals.size() != getNumRows(), std::runtime_error,
        Teuchos::typeName(*this) << "::set2DValues(inds): vals must have as many entries as the number of rows in the associated graph.");
    this->clear();
    values2D_  = vals;
  }


  // ======= finalize ===========
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  void CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::finalize(bool OptimizeStorage)
  {
    if (isFinalized_ && !(OptimizeStorage == true && staticGraph_->isOptimized() == false)) return;
    // we only have something to do if we have a matrix-owned graph
    if (myGraph_ != NULL) {
      myGraph_->template finalize<Scalar>(OptimizeStorage, values2D_, values1D_);
    }
    else {
      TEST_FOR_EXCEPTION(OptimizeStorage == true && staticGraph_->isOptimized() == false, std::runtime_error, 
          Teuchos::typeName(*this) << "::finalize(OptimizeStorage == true): underlying static graph is not optimized and cannot be optimized.");
    }
    isFinalized_ = true;
  }


  //=========================================================================================================================
  // 
  // A device-resident CrsMatrix
  // 
  //=========================================================================================================================


  /** \brief A default device-compute compressed-row sparse matrix.
      \ingroup kokkos_crs_ops

      This is externally identical to the host-based matrix; in fact, it
      derives from CrsMatrixHostCompute. The difference is that that it
      contains additional storage and logic for device-bound compute buffers.
   */
  template <class Scalar,
            class Ordinal, 
            class Node,
            class LocalMatOps>
  class CrsMatrixDeviceCompute : public CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps> {
  public:

    //! @name Constructors/Destructor
    //@{

    //! Default constructor with a no graph (must be set later).
    CrsMatrixDeviceCompute();

    //! Constructor with a matrix-owned non-const graph
    CrsMatrixDeviceCompute(CrsGraphDeviceCompute<Ordinal,Node,LocalMatOps> &graph);

    //! Constructor with a non-owned const graph.
    CrsMatrixDeviceCompute(const CrsGraphDeviceCompute<Ordinal,Node,LocalMatOps> &graph);

    //! CrsMatrixDeviceCompute Destructor
    ~CrsMatrixDeviceCompute();

    //@}

    //! @name Methods over-riding CrsGraphDeviceCompute.
    //@{

    //! Set matrix-owned graph.
    void setOwnedGraph(CrsGraphDeviceCompute<Ordinal,Node,LocalMatOps> &graph);

    //! Set static graph.
    void setStaticGraph(const CrsGraphDeviceCompute<Ordinal,Node,LocalMatOps> &graph);
    
    //! Instruct the matrix to perform any necessary manipulation, including (optionally) optimizing the storage of the matrix data.
    /** 
          @param[in] OptimizeStorage   Permit the matrix to reallocate storage on the host in order to provide optimal storage and/or performance.
          \post if OptimizeStorage == true, then is2DStructure() == true
     */
    void finalize(bool OptimizeStorage);

    //! Return the device-bound buffers.
    void getDeviceBuffer(ArrayRCP<Scalar> &d_vals) const;

    //! Release data associated with this graph.
    virtual void clear();

    //@}

  protected:
    //! Copy constructor (protected and not implemented) 
    CrsMatrixDeviceCompute(const CrsMatrixDeviceCompute& sources);

    // pointers to device-based graph
    CrsGraphDeviceCompute<Ordinal,Node,LocalMatOps> *myDeviceGraph_;
    const CrsGraphDeviceCompute<Ordinal,Node,LocalMatOps> *staticDeviceGraph_;
    // device storage (always 1D packed)
    ArrayRCP<Scalar> pbuf_values_;
  };

  //==============================================================================
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  CrsMatrixDeviceCompute<Scalar,Ordinal,Node,LocalMatOps>::CrsMatrixDeviceCompute()
  : CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>()
  , myDeviceGraph_(NULL)
  , staticDeviceGraph_(NULL)
  {}

  //==============================================================================
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  CrsMatrixDeviceCompute<Scalar,Ordinal,Node,LocalMatOps>::CrsMatrixDeviceCompute(const CrsGraphDeviceCompute<Ordinal,Node,LocalMatOps> &graph)
  : CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>(graph)
  , myDeviceGraph_(NULL)
  , staticDeviceGraph_(&graph)
  {}

  //==============================================================================
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  CrsMatrixDeviceCompute<Scalar,Ordinal,Node,LocalMatOps>::CrsMatrixDeviceCompute(CrsGraphDeviceCompute<Ordinal,Node,LocalMatOps> &graph)
  : CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>(graph)
  , myDeviceGraph_(&graph)
  , staticDeviceGraph_(&graph)
  {}

  //===== destructor =====
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  CrsMatrixDeviceCompute<Scalar,Ordinal,Node,LocalMatOps>::~CrsMatrixDeviceCompute() 
  {}

  //===== clear =====
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  void CrsMatrixDeviceCompute<Scalar,Ordinal,Node,LocalMatOps>::clear() { 
    CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::clear();
    pbuf_values_ = null;
  }

  // ======= get device ===========
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  void CrsMatrixDeviceCompute<Scalar,Ordinal,Node,LocalMatOps>::getDeviceBuffer(ArrayRCP<Scalar> &vals) const
  {
    vals = pbuf_values_;
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  void CrsMatrixDeviceCompute<Scalar,Ordinal,Node,LocalMatOps>::setStaticGraph(const CrsGraphDeviceCompute<Ordinal,Node,LocalMatOps> &dgraph)
  {
    CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::setStaticGraph(dgraph);
    myDeviceGraph_     = NULL;
    staticDeviceGraph_ = &dgraph;
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  void CrsMatrixDeviceCompute<Scalar,Ordinal,Node,LocalMatOps>::setOwnedGraph(CrsGraphDeviceCompute<Ordinal,Node,LocalMatOps> &dgraph)
  {
    CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::setOwnedGraph(dgraph);
    myDeviceGraph_     = &dgraph;
    staticDeviceGraph_ = &dgraph;
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  void CrsMatrixDeviceCompute<Scalar,Ordinal,Node,LocalMatOps>::finalize(bool OptimizeStorage)
  {
    if (this->isFinalized() && !(OptimizeStorage == true && staticDeviceGraph_->isOptimized() == false)) return;
    if (myDeviceGraph_ != NULL) {
      // this will finalize both classes and copy data to the device
      myDeviceGraph_->template finalize<Scalar>(OptimizeStorage, this->values2D_, this->values1D_, pbuf_values_);
    }
    else {
      TEST_FOR_EXCEPTION(OptimizeStorage == true && this->staticGraph_->isOptimized() == false, std::runtime_error, 
          Teuchos::typeName(*this) << "::finalize(OptimizeStorage == true): underlying static graph is not optimized and cannot be optimized.");
      // static graph. no changes, but need to copy values to device.
      if (this->isEmpty()) {
        pbuf_values_  = null;
      }
      else {
        // allocate space on the device and copy data there, in a packed format
        pbuf_values_ = this->getNode()->template allocBuffer<Scalar>(this->getNumEntries());
        ArrayRCP<size_t> rowBegs, rowEnds, numEntriesPerRow;
        if (this->is1DStructure()) {
          ArrayRCP<Ordinal> inds;
          const_cast<CrsGraphDeviceCompute<Ordinal,Node,LocalMatOps> *>(staticDeviceGraph_)->get1DStructure(inds, rowBegs, rowEnds);
        }
        else {
          ArrayRCP<ArrayRCP<Ordinal> > inds;
          const_cast<CrsGraphDeviceCompute<Ordinal,Node,LocalMatOps> *>(staticDeviceGraph_)->get2DStructure(inds, numEntriesPerRow);
        }
        if (this->isOptimized()) {
          // should be packed now; single copy should do, and offsets are rowBegs_
          this->getNode()->template copyToBuffer<Scalar >(this->getNumEntries(), this->values1D_(0,this->getNumEntries()), pbuf_values_);
        }
        else {
          // need to get graph structure from graph
          ArrayRCP<Scalar> view_values = this->getNode()->template viewBufferNonConst<Scalar >(WriteOnly, pbuf_values_.size(), pbuf_values_);
          typename ArrayRCP<Scalar >::iterator oldvals, newvals;
          newvals = view_values.begin();
          size_t curnuminds;
          for (size_t i=0; i < this->getNumRows(); ++i) {
            if (this->is1DStructure()) {
              curnuminds = rowEnds[i] - rowBegs[i];
              oldvals = this->values1D_.begin() + rowBegs[i];
            }
            else {
              curnuminds = numEntriesPerRow[i];
              oldvals = this->values2D_[i].begin();
            }
            std::copy(oldvals, oldvals+curnuminds, newvals);
            newvals += curnuminds;
          }
          view_values  = null;
        }
      }
    }
    this->isFinalized_ = true;
  }


  //=========================================================================================================================
  // 
  // A first-touch allocation host-resident CrsMatrix
  // 
  //=========================================================================================================================

  /** \brief A host-compute compressed-row sparse matrix with first-touch allocation.
      \ingroup kokkos_crs_ops
   */
  template <class Scalar, 
            class Ordinal, 
            class Node,
            class LocalMatOps>
  class FirstTouchHostCrsMatrix : public CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps> {
  public:

    typedef typename CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::ScalarType        ScalarType;
    typedef typename CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::OrdinalType       OrdinalType;
    typedef typename CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::NodeType          NodeType;
    typedef typename CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::LocalMatOpsType   LocalMatOpsType;

    //! @name Constructors/Destructor
    //@{

    //! Default constructor with no graph. (Must be set later.) 
    FirstTouchHostCrsMatrix();

    //! Constructor with a matrix-owned non-const graph
    FirstTouchHostCrsMatrix(FirstTouchHostCrsGraph<Ordinal,Node,LocalMatOps> &graph);

    //! Constructor with a non-owned const graph.
    FirstTouchHostCrsMatrix(const FirstTouchHostCrsGraph<Ordinal,Node,LocalMatOps> &graph);

    //! FirstTouchHostCrsMatrix Destructor
    virtual ~FirstTouchHostCrsMatrix();

    //@}

    //! @name Graph set routines.
    //@{

    //! Set matrix-owned graph.
    void setOwnedGraph(FirstTouchHostCrsGraph<Ordinal,Node,LocalMatOps> &graph);

    //! Set static graph.
    void setStaticGraph(const FirstTouchHostCrsGraph<Ordinal,Node,LocalMatOps> &graph);
    
    //@}

    //! @name Data entry and accessor methods.
    //@{

    //! Instruct the matrix to perform any necessary manipulation, including (optionally) optimizing the storage of the matrix data.
    /** 
          If the matrix is associated with a matrix-owned graph, then it will be finalized as well. A static graph will not be modified.

          @param[in] OptimizeStorage   Permit the graph to reallocate storage on the host in order to provide optimal storage and/or performance.
          \post if OptimizeStorage == true, then is2DStructure() == true
     */
    void finalize(bool OptimizeStorage);

    //@}

  protected:
    //! Copy constructor (protected and not implemented) 
    FirstTouchHostCrsMatrix(const FirstTouchHostCrsMatrix& sources);
    // pointers to first-touch graphs
    FirstTouchHostCrsGraph<Ordinal,Node,LocalMatOps>           *myFTGraph_;
    const FirstTouchHostCrsGraph<Ordinal,Node,LocalMatOps> *staticFTGraph_;
  };


  //==============================================================================
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  FirstTouchHostCrsMatrix<Scalar,Ordinal,Node,LocalMatOps>::FirstTouchHostCrsMatrix()
  : CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>() 
  , myFTGraph_(NULL)
  , staticFTGraph_(NULL)
  {}

  //==============================================================================
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  FirstTouchHostCrsMatrix<Scalar,Ordinal,Node,LocalMatOps>::FirstTouchHostCrsMatrix(const FirstTouchHostCrsGraph<Ordinal,Node,LocalMatOps> &graph)
  : CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>(graph) 
  , myFTGraph_(NULL)
  , staticFTGraph_(&graph)
  {}

  //==============================================================================
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  FirstTouchHostCrsMatrix<Scalar,Ordinal,Node,LocalMatOps>::FirstTouchHostCrsMatrix(FirstTouchHostCrsGraph<Ordinal,Node,LocalMatOps> &graph)
  : CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>(graph) 
  , myFTGraph_(&graph)
  , staticFTGraph_(&graph)
  {}

  //==============================================================================
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  FirstTouchHostCrsMatrix<Scalar,Ordinal,Node,LocalMatOps>::~FirstTouchHostCrsMatrix() {
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  void FirstTouchHostCrsMatrix<Scalar,Ordinal,Node,LocalMatOps>::setStaticGraph(const FirstTouchHostCrsGraph<Ordinal,Node,LocalMatOps> &dgraph)
  {
    CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::setStaticGraph(dgraph);
    myFTGraph_     = NULL;
    staticFTGraph_ = &dgraph;
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  void FirstTouchHostCrsMatrix<Scalar,Ordinal,Node,LocalMatOps>::setOwnedGraph(FirstTouchHostCrsGraph<Ordinal,Node,LocalMatOps> &dgraph)
  {
    CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::setOwnedGraph(dgraph);
    myFTGraph_     = &dgraph;
    staticFTGraph_ = &dgraph;
  }

  // ======= finalize ===========
  template <class Scalar, class Ordinal, class Node, class LocalMatOps>
  void FirstTouchHostCrsMatrix<Scalar,Ordinal,Node,LocalMatOps>::finalize(bool OptimizeStorage)
  {
    if (this->isFinalized_ == false || 
       (OptimizeStorage == true && staticFTGraph_->isOptimized() == false))
    {
      if (this->myFTGraph_ != NULL) {
        // call graph finalize to finalize/first-touch graph and matrix at the same time
        myFTGraph_->template finalize<Scalar>(OptimizeStorage, this->values2D_, this->values1D_);
      }
      else {
        TEST_FOR_EXCEPTION(OptimizeStorage == true && staticFTGraph_->isOptimized() == false, std::runtime_error, 
            Teuchos::typeName(*this) << "::finalize(OptimizeStorage == true): underlying static graph is not already optimized and cannot be optimized.");
      }
      this->isFinalized_ = true;
    }
#ifdef HAVE_KOKKOS_DEBUG
    TEST_FOR_EXCEPTION(this->isFinalized_ == false || staticFTGraph_->isFinalized() == false, std::logic_error,
        Teuchos::typeName(*this) << "::finalize(): logic error. Post-condition violated. Please contact Tpetra team.");
    TEST_FOR_EXCEPTION(OptimizeStorage == true && (staticFTGraph_->isOptimized() == false && staticFTGraph_->isEmpty() == false), std::logic_error,
        Teuchos::typeName(*this) << "::finalize(): logic error. Post-condition violated. Please contact Tpetra team.");
#endif
  }


  //=========================================================================================================================
  // 
  // Specializations
  // 
  //=========================================================================================================================
  // C++0x: it would be nice if these were template aliases

  /** \brief Kokkos compressed-row sparse matrix class.
      \ingroup kokkos_crs_ops

      Default specialization is a host-bound CrsMatrixHostCompute object.
    */
  template <class Scalar,
            class Ordinal, 
            class Node,
            class LocalMatOps>
  class CrsMatrix : public CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps> {
  public:
    CrsMatrix() : CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>() {}
    CrsMatrix(CrsGraph<Ordinal,Node,LocalMatOps> &graph) : CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>(graph) {}
    CrsMatrix(const CrsGraph<Ordinal,Node,LocalMatOps> &graph) : CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>(graph) {}
    void setStaticGraph(const CrsGraph<Ordinal,Node,LocalMatOps> &graph) {
      const CrsGraphHostCompute<Ordinal,Node,LocalMatOps> &hgraph = dynamic_cast<const CrsGraphHostCompute<Ordinal,Node,LocalMatOps> &>(graph);
      CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::setStaticGraph(hgraph);
    }
    void setOwnedGraph(CrsGraph<Ordinal,Node,LocalMatOps> &graph) {
      CrsGraphHostCompute<Ordinal,Node,LocalMatOps> &hgraph = dynamic_cast<CrsGraphHostCompute<Ordinal,Node,LocalMatOps> &>(graph);
      CrsMatrixHostCompute<Scalar,Ordinal,Node,LocalMatOps>::setOwnedGraph(hgraph);
    }
  private:
    CrsMatrix(const CrsMatrix<Scalar,Ordinal,Node,LocalMatOps> &mat); // not implemented
  };

#ifndef HAVE_KOKKOS_NO_FIRST_TOUCH_MATVEC_ALLOCATION
#ifdef HAVE_KOKKOS_TBB
  /** \brief Kokkos compressed-row sparse matrix class.
      \ingroup kokkos_crs_ops

      Specialization is a first-touch host-bound FirstTouchHostCrsMatrix object.
    */
  class TBBNode;
  template <class Scalar,
            class Ordinal, 
            class LocalMatOps>
  class CrsMatrix<Scalar,Ordinal,TBBNode,LocalMatOps> : public FirstTouchHostCrsMatrix<Scalar,Ordinal,TBBNode,LocalMatOps> {
  public:
    CrsMatrix() : FirstTouchHostCrsMatrix<Scalar,Ordinal,TBBNode,LocalMatOps>() {}
    CrsMatrix(CrsGraph<Ordinal,TBBNode,LocalMatOps> &graph) : FirstTouchHostCrsMatrix<Scalar,Ordinal,TBBNode,LocalMatOps>(graph) {}
    CrsMatrix(const CrsGraph<Ordinal,TBBNode,LocalMatOps> &graph) : FirstTouchHostCrsMatrix<Scalar,Ordinal,TBBNode,LocalMatOps>(graph) {}
    void setStaticGraph(const CrsGraph<Ordinal,TBBNode,LocalMatOps> &graph) {
      const FirstTouchHostCrsGraph<Ordinal,TBBNode,LocalMatOps> &hgraph = dynamic_cast<const FirstTouchHostCrsGraph<Ordinal,TBBNode,LocalMatOps> &>(graph);
      FirstTouchHostCrsMatrix<Scalar,Ordinal,TBBNode,LocalMatOps>::setStaticGraph(hgraph);
    }
    void setOwnedGraph(CrsGraph<Ordinal,TBBNode,LocalMatOps> &graph) {
      FirstTouchHostCrsGraph<Ordinal,TBBNode,LocalMatOps> &hgraph = dynamic_cast<FirstTouchHostCrsGraph<Ordinal,TBBNode,LocalMatOps> &>(graph);
      FirstTouchHostCrsMatrix<Scalar,Ordinal,TBBNode,LocalMatOps>::setOwnedGraph(hgraph);
    }
  private:
    CrsMatrix(const CrsMatrix<Scalar,Ordinal,TBBNode,LocalMatOps> &mat); // not implemented
  };
#endif
#ifdef HAVE_KOKKOS_THREADPOOL
  class TPINode;
  /** \brief Kokkos compressed-row sparse matrix class.
      \ingroup kokkos_crs_ops

      Specialization is a first-touch host-bound FirstTouchHostCrsMatrix object.
    */
  template <class Scalar,
            class Ordinal, 
            class LocalMatOps>
  class CrsMatrix<Scalar,Ordinal,TPINode,LocalMatOps> : public FirstTouchHostCrsMatrix<Scalar,Ordinal,TPINode,LocalMatOps> {
  public:
    CrsMatrix() : FirstTouchHostCrsMatrix<Scalar,Ordinal,TPINode,LocalMatOps>() {}
    CrsMatrix(CrsGraph<Ordinal,TPINode,LocalMatOps> &graph) : FirstTouchHostCrsMatrix<Scalar,Ordinal,TPINode,LocalMatOps>(graph) {}
    CrsMatrix(const CrsGraph<Ordinal,TPINode,LocalMatOps> &graph) : FirstTouchHostCrsMatrix<Scalar,Ordinal,TPINode,LocalMatOps>(graph) {}
    void setStaticGraph(const CrsGraph<Ordinal,TPINode,LocalMatOps> &graph) {
      const FirstTouchHostCrsGraph<Ordinal,TPINode,LocalMatOps> &hgraph = dynamic_cast<const FirstTouchHostCrsGraph<Ordinal,TPINode,LocalMatOps> &>(graph);
      FirstTouchHostCrsMatrix<Scalar,Ordinal,TPINode,LocalMatOps>::setStaticGraph(hgraph);
    }
    void setOwnedGraph(CrsGraph<Ordinal,TPINode,LocalMatOps> &graph) {
      FirstTouchHostCrsGraph<Ordinal,TPINode,LocalMatOps> &hgraph = dynamic_cast<FirstTouchHostCrsGraph<Ordinal,TPINode,LocalMatOps> &>(graph);
      FirstTouchHostCrsMatrix<Scalar,Ordinal,TPINode,LocalMatOps>::setOwnedGraph(hgraph);
    }
  private:
    CrsMatrix(const CrsMatrix<Scalar,Ordinal,TPINode,LocalMatOps> &mat); // not implemented
  };
#endif
#endif

  /** \brief Kokkos compressed-row sparse matrix class.
      \ingroup kokkos_crs_ops

      Specialization for device-based matrix operation is a CrsMatrixDeviceCompute.
    */
  template <class S, class O, class N> class DefaultDeviceSparseOps;
  template <class Scalar,
            class Ordinal,
            class Node>
  class CrsMatrix<Scalar,Ordinal,Node,DefaultDeviceSparseOps<void,Ordinal,Node> > : public CrsMatrixDeviceCompute<Scalar,Ordinal,Node,DefaultDeviceSparseOps<void,Ordinal,Node> > {
  public:
    typedef DefaultDeviceSparseOps<void,Ordinal,Node> LocalMatOps;
    CrsMatrix() : CrsMatrixDeviceCompute<Scalar,Ordinal,Node,LocalMatOps >() {}
    CrsMatrix(CrsGraph<Ordinal,Node,LocalMatOps > &graph) : CrsMatrixDeviceCompute<Scalar,Ordinal,Node,LocalMatOps >(graph) {}
    CrsMatrix(const CrsGraph<Ordinal,Node,LocalMatOps > &graph) : CrsMatrixDeviceCompute<Scalar,Ordinal,Node,LocalMatOps >(graph) {}
    void setStaticGraph(const CrsGraph<Ordinal,Node,LocalMatOps > &graph) {
      const CrsGraphDeviceCompute<Ordinal,Node,LocalMatOps > &dgraph = dynamic_cast<const CrsGraphDeviceCompute<Ordinal,Node,LocalMatOps > &>(graph);
      CrsMatrixDeviceCompute<Scalar,Ordinal,Node,LocalMatOps >::setStaticGraph(dgraph);
    }
    void setOwnedGraph(CrsGraph<Ordinal,Node,LocalMatOps > &graph) {
      CrsGraphDeviceCompute<Ordinal,Node,LocalMatOps > &dgraph = dynamic_cast<CrsGraphDeviceCompute<Ordinal,Node,LocalMatOps > &>(graph);
      CrsMatrixDeviceCompute<Scalar,Ordinal,Node,LocalMatOps >::setOwnedGraph(dgraph);
    }
  private:
    CrsMatrix(const CrsMatrix<Scalar,Ordinal,Node,LocalMatOps > &mat); // not implemented
  };

} // namespace Kokkos

#endif /* KOKKOS_CRSMATRIX_HPP */
