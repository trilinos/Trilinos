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

#include "Kokkos_CrsGraph.hpp"
#include "Kokkos_CrsMatrix.hpp"

namespace Kokkos {

  //=========================================================================================================================
  // 
  // A first-touch sparse ops
  // 
  //=========================================================================================================================

  /** \brief Sparse matrix-vector multiplication and solve routines with first-touch semantics. 
      \ingroup kokkos_crs_ops

      Effectively, this class is trigger so that downstream software will use the FirstTouchHostCrsGraph and FirstTOuchHostCrsMatrix
      classes. The sparse mat-vec and solve operations are inherited from DefautlHostSparseOps.
   */
  template <class Scalar, class Ordinal, class Node>
  class FirstTouchSparseOps {
  public:
    //@{ 
    //! @name Typedefs and structs

    //!
    typedef Scalar  ScalarType;
    //!
    typedef Ordinal OrdinalType;
    //!
    typedef Node    NodeType;


    /** \brief Typedef for local graph class */
    template <class O, class N>
    struct graph {
      typedef FirstTouchHostCrsGraph<O,N> type;
    };

    /** \brief Typedef for local matrix class */
    template <class S, class O, class N>
    struct matrix {
      typedef FirstTouchHostCrsMatrix<S,O,N> type;
    };

    /** \brief Rebind struct, for specifying type information for a different scalar.
          
        This specifies a DefaultHostSparseOps object, regardless of scalar type.
      */
    template <class S2>
    struct rebind {
      typedef DefaultHostSparseOps<S2,Ordinal,Node> type;
    };

    //@}
    //! @name Constructors/Destructor
    //@{

    //! \brief Constructor accepting and retaining a node object.
    FirstTouchHostSparseOps(const RCP<Node> &node) : DefaultHostSparseOps(node) {}

    //! DefaultHostSparseOps Destructor
    ~FirstTouchHostSparseOps();

    //@}
    //! @name Overloaded from DefaultHostSparseOps for first-touch CRS objects.
    //@{

    void initializeStructure(const FirstTouchHostCrsGraph<Ordinal,Node> &graph) {
      DefaultHostSparseOps<Scalar,Ordinal,Node>::initializeStructure(graph);
    }

    void initializeValues(const FirstTouchHostCrsMatrix<Scalar,Ordinal,Node> &matrix) {
      DefaultHostSparseOps<Scalar,Ordinal,Node>::initializeValues(matrix);
    }

    //@}

  protected:
    //! Copy constructor (protected and unimplemented)
    FirstTouchHostSparseOps(const DefaultHostSparseOps& source);
  };


  //=========================================================================================================================
  // 
  // A first-touch allocation host-resident CrsGraph
  // 
  //=========================================================================================================================

  /** \brief A host-compute compressed-row sparse graph with first-touch allocation.
      \ingroup kokkos_crs_ops
   */
  template <class Ordinal, 
            class Node>
  class FirstTouchHostCrsGraph : public CrsGraphHostCompute<Ordinal,Node> {
  public:

    typedef typename CrsGraphHostCompute<Ordinal,Node>::OrdinalType      OrdinalType;
    typedef typename CrsGraphHostCompute<Ordinal,Node>::NodeType         NodeType;

    //! @name Constructors/Destructor
    //@{

    //! Default FirstTouchHostCrsGraph constuctor.
    FirstTouchHostCrsGraph(size_t numRows, const RCP<Node> &node);

    //! FirstTouchHostCrsGraph Destructor
    virtual ~FirstTouchHostCrsGraph();

    //@}

    //! @name Data entry and accessor methods.
    //@{

    //! Instruct the graph to perform any necessary manipulation, including (optionally) optimizing the storage of the graph data.
    /** 
          @param[in] OptimizeStorage   Permit the graph to reallocate storage on the host in order to provide optimal storage and/or performance.
          \post if OptimizeStorage == true, then is2DStructure() == true
     */
    void finalize(bool OptimizeStorage);

    /** 
          @param[in] OptimizeStorage   Permit the graph to reallocate storage on the host in order to provide optimal storage and/or performance.
          @param[in/out] values2D      2D-structured matrix values. Required to be non-null if is2DStructure() is true. Set to null if OptimizeStorage is true.
          @param[in/out] values1D      1D-structured matrix values. Required to be non-null if is1DStructure() is true. Allocated if OptimizeStorage is true.
          \post if OptimizeStorage == true or already is2DStructure(), then is2DStructure() == true.
     */
    template <class Scalar>
    void finalize(bool OptimizeStorage, ArrayRCP<ArrayRCP<Scalar> > &values2D, ArrayRCP<Scalar> &values1D);

    //@}

  protected:
    //! Copy constructor (protected and not implemented) 
    FirstTouchHostCrsGraph(const FirstTouchHostCrsGraph& sources);
  };

#ifndef KERNEL_PREFIX
#define KERNEL_PREFIX
#endif
  /// \class FirstTouchCopyIndicesKernel
  ///
  /// Kokkos kernel for copying array indices using a first-touch
  /// initialization strategy for CPU memory.
  ///
  /// \note We have to store member data as raw pointers, rather than
  /// ArrayRCPs, because ArrayRCP are not thread safe, and the arrays
  /// get accessed inside a Kokkos kernel.
  template <class T>
  struct FirstTouchCopyIndicesKernel {
    const size_t * numEntriesPerRow;
    const size_t * offsets1D;
    T * data1D;
    const ArrayRCP<T> * data2D;

    inline KERNEL_PREFIX void execute(size_t row) {
      const size_t rowNumInds = numEntriesPerRow[row];
      const T* const oldinds = data2D[row].getRawPtr();
      T* const newinds = data1D + offsets1D[row];
      std::copy(oldinds, oldinds+rowNumInds, newinds);
    }
  };

  //==============================================================================
  template <class Ordinal, class Node>
  FirstTouchHostCrsGraph<Ordinal,Node>::FirstTouchHostCrsGraph(size_t numRows, const RCP<Node> &node) 
  : CrsGraphHostCompute<Ordinal,Node>(numRows,node) 
  {
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;
  }

  //==============================================================================
  template <class Ordinal, class Node>
  FirstTouchHostCrsGraph<Ordinal,Node>::~FirstTouchHostCrsGraph() {}

  // ======= finalize ===========
  template <class Ordinal, class Node>
  void FirstTouchHostCrsGraph<Ordinal,Node>::finalize(bool OptimizeStorage)
  {
    // allocations not done using the Node. no current need for host-based nodes, and 
    // this leads to incorrect behavior when we try to reuse this code from child CrsGraphDeviceCompute
    if (this->isFinalized() && !(OptimizeStorage == true && this->isOptimized() == false)) return;
    if ((this->indices1D_ == null && this->indices2D_ == null) || (this->getNumEntries() == 0)) {
      this->isEmpty_ = true;
    }
    else {
      this->isEmpty_ = false;
      if (OptimizeStorage) {
        // move into packed 1D storage
        ArrayRCP<size_t> offsets = arcp<size_t>(this->numRows_+1);
        if (this->is2DStructure() == true) {
          // 2D to 1D packed: first-touch allocation
          // allocate 1D storage
          this->indices1D_ = arcp<Ordinal>(this->getNumEntries());
          // compute offset array on host thread
          {
            size_t curoffset = 0;
            for (size_t i=0; i < this->numRows_; ++i) {
              offsets[i] = curoffset;
              curoffset += this->numEntriesPerRow_[i];
            }
            offsets[this->numRows_] = curoffset;
            TEST_FOR_EXCEPTION( curoffset != this->getNumEntries(), std::logic_error, 
                Teuchos::typeName(*this) << "::finalize(): Internal logic error. Please contact Kokkos team.");
          }
          FirstTouchCopyIndicesKernel<Ordinal> kern;
          kern.offsets1D = offsets.getRawPtr();
          kern.numEntriesPerRow = this->numEntriesPerRow_.getRawPtr();
          kern.data1D = this->indices1D_.getRawPtr();
          kern.data2D = this->indices2D_.getRawPtr();
          this->getNode()->template parallel_for<FirstTouchCopyIndicesKernel<Ordinal> >(0,this->numRows_,kern);
        }
        else {
          // 1D to 1D packed: no first-touch
          // copy/pack data
          size_t curoffset = 0;
          size_t curnuminds;
          typename ArrayRCP<Ordinal>::iterator oldinds, newinds;
          newinds = this->indices1D_.begin();
          for (size_t i=0; i < this->numRows_; ++i) {
            offsets[i] = curoffset;
            curnuminds = this->rowEnds_[i] - this->rowBegs_[i];
            oldinds = this->indices1D_.begin() + this->rowBegs_[i];
            std::copy(oldinds, oldinds+curnuminds, newinds);
            newinds += curnuminds;
            curoffset += curnuminds;
          }
          offsets[this->numRows_] = curoffset;
          TEST_FOR_EXCEPTION( curoffset != this->getNumEntries(), std::logic_error, 
              Teuchos::typeName(*this) << "::finalize(): Internal logic error. Please contact Kokkos team.");
        }
        // done with the original row beg/end offsets, can point to the new overlapping one
        this->rowBegs_   = offsets;
        this->rowEnds_   = offsets.persistingView(1,this->numRows_);
        this->isOpt_     = true;
        this->is1D_      = true;
        // delete 2D storage (if any)
        this->is2D_      = false;
        this->numEntriesPerRow_ = null;
        this->indices2D_        = null;
      }
    }
    this->isFinalized_ = true;
  }


  // ======= finalize ===========
  // finalize() storage for the graph with associated matrix values
  // this is called from a CrsMatrix, and we're doing the finalize the for the graph and matrix at the same time, so the matrix doesn't have to.
  template <class Ordinal, class Node>
  template <class Scalar>
  void FirstTouchHostCrsGraph<Ordinal,Node>::finalize(bool OptimizeStorage, ArrayRCP<ArrayRCP<Scalar> > &values2D, ArrayRCP<Scalar> &values1D)
  {
    if (this->isFinalized() && !(OptimizeStorage == true && this->isOptimized() == false)) return;
    if ((this->indices1D_ == null && this->indices2D_ == null) || (this->getNumEntries() == 0)) {
      this->isEmpty_ = true;
    }
    else {
      this->isEmpty_ = false;
      // move into packed 1D storage
      if (OptimizeStorage) {
        ArrayRCP<size_t> offsets = arcp<size_t>(this->numRows_+1);
        if (this->is2DStructure() == true) {
          // 2D to 1D packed: first-touch allocation
          // allocate 1D storage
          this->indices1D_ = arcp<Ordinal>(this->getNumEntries());
          values1D         = arcp<Scalar >(this->getNumEntries());
          // compute offset array on host thread
          {
            size_t curoffset = 0;
            for (size_t i=0; i < this->numRows_; ++i) {
              offsets[i] = curoffset;
              curoffset += this->numEntriesPerRow_[i];
            }
            offsets[this->numRows_] = curoffset;
            TEST_FOR_EXCEPTION( curoffset != this->getNumEntries(), std::logic_error, 
                Teuchos::typeName(*this) << "::finalize(): Internal logic error. Please contact Kokkos team.");
          }
          {
            FirstTouchCopyIndicesKernel<Ordinal> indskern;
            indskern.offsets1D = offsets.getRawPtr();
            indskern.numEntriesPerRow = this->numEntriesPerRow_.getRawPtr();
            indskern.data1D = this->indices1D_.getRawPtr();
            indskern.data2D = this->indices2D_.getRawPtr();
            this->getNode()->template parallel_for<FirstTouchCopyIndicesKernel<Ordinal> >(0,this->numRows_,indskern);
          }
          {
            FirstTouchCopyIndicesKernel<Scalar> valskern;
            valskern.offsets1D = offsets.getRawPtr();
            valskern.numEntriesPerRow = this->numEntriesPerRow_.getRawPtr();
            valskern.data1D = values1D.getRawPtr();
            valskern.data2D = values2D.getRawPtr();
            this->getNode()->template parallel_for<FirstTouchCopyIndicesKernel<Scalar> >(0,this->numRows_,valskern);
          }
        }
        else {
          // copy/pack data
          // 1D to 1D packed: no first-touch
          size_t curoffset = 0;
          size_t curnuminds;
          typename ArrayRCP<Ordinal>::iterator oldinds, newinds;
          typename ArrayRCP<Scalar >::iterator oldvals, newvals;
          newinds = this->indices1D_.begin();
          newvals = values1D.begin();
          for (size_t i=0; i < this->numRows_; ++i) {
            offsets[i] = curoffset;
            curnuminds = this->rowEnds_[i] - this->rowBegs_[i];
            oldinds = this->indices1D_.begin() + this->rowBegs_[i];
            oldvals = values1D.begin() + this->rowBegs_[i];
            std::copy(oldinds, oldinds+curnuminds, newinds);
            std::copy(oldvals, oldvals+curnuminds, newvals);
            newinds += curnuminds;
            newvals += curnuminds;
            curoffset += curnuminds;
          }
          offsets[this->numRows_] = curoffset;
          TEST_FOR_EXCEPTION( curoffset != this->getNumEntries(), std::logic_error, 
			      Teuchos::typeName(*this) << "::finalize(): "
			      "Internal logic error: curoffset (= " 
			      << curoffset << ") != this->getNumEntries() (= "
			      << this->getNumEntries() 
			      << ").  Please contact Kokkos team.");
        }
        // done with the original row beg/end offsets, can point to the new overlapping one
        this->rowBegs_   = offsets;
        this->rowEnds_   = offsets.persistingView(1,this->numRows_);
        this->is1D_      = true;
        this->isOpt_     = true;
        // delete 2D storage (if there was any)
        this->is2D_      = false;
        this->numEntriesPerRow_ = null;
        this->indices2D_        = null;
        values2D          = null;
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
            class Node>
  class FirstTouchHostCrsMatrix : public CrsMatrixHostCompute<Scalar,Ordinal,Node> {
  public:

    typedef typename CrsMatrixHostCompute<Scalar,Ordinal,Node>::ScalarType   ScalarType;
    typedef typename CrsMatrixHostCompute<Scalar,Ordinal,Node>::OrdinalType  OrdinalType;
    typedef typename CrsMatrixHostCompute<Scalar,Ordinal,Node>::NodeType     NodeType;

    //! @name Constructors/Destructor
    //@{

    //! Default constructor with no graph. (Must be set later.) 
    FirstTouchHostCrsMatrix();

    //! Constructor with a matrix-owned non-const graph
    FirstTouchHostCrsMatrix(FirstTouchHostCrsGraph<Ordinal,Node> &graph);

    //! Constructor with a non-owned const graph.
    FirstTouchHostCrsMatrix(const FirstTouchHostCrsGraph<Ordinal,Node> &graph);

    //! FirstTouchHostCrsMatrix Destructor
    virtual ~FirstTouchHostCrsMatrix();

    //@}

    //! @name Graph set routines.
    //@{

    //! Set matrix-owned graph.
    void setOwnedGraph(FirstTouchHostCrsGraph<Ordinal,Node> &graph);

    //! Set static graph.
    void setStaticGraph(const FirstTouchHostCrsGraph<Ordinal,Node> &graph);
    
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
    FirstTouchHostCrsGraph<Ordinal,Node>           *myFTGraph_;
    const FirstTouchHostCrsGraph<Ordinal,Node> *staticFTGraph_;
  };


  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  FirstTouchHostCrsMatrix<Scalar,Ordinal,Node>::FirstTouchHostCrsMatrix()
  : CrsMatrixHostCompute<Scalar,Ordinal,Node>() 
  , myFTGraph_(NULL)
  , staticFTGraph_(NULL)
  {}

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  FirstTouchHostCrsMatrix<Scalar,Ordinal,Node>::FirstTouchHostCrsMatrix(const FirstTouchHostCrsGraph<Ordinal,Node> &graph)
  : CrsMatrixHostCompute<Scalar,Ordinal,Node>(graph) 
  , myFTGraph_(NULL)
  , staticFTGraph_(&graph)
  {}

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  FirstTouchHostCrsMatrix<Scalar,Ordinal,Node>::FirstTouchHostCrsMatrix(FirstTouchHostCrsGraph<Ordinal,Node> &graph)
  : CrsMatrixHostCompute<Scalar,Ordinal,Node>(graph) 
  , myFTGraph_(&graph)
  , staticFTGraph_(&graph)
  {}

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  FirstTouchHostCrsMatrix<Scalar,Ordinal,Node>::~FirstTouchHostCrsMatrix() {
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  void FirstTouchHostCrsMatrix<Scalar,Ordinal,Node>::setStaticGraph(const FirstTouchHostCrsGraph<Ordinal,Node> &dgraph)
  {
    CrsMatrixHostCompute<Scalar,Ordinal,Node>::setStaticGraph(dgraph);
    myFTGraph_     = NULL;
    staticFTGraph_ = &dgraph;
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  void FirstTouchHostCrsMatrix<Scalar,Ordinal,Node>::setOwnedGraph(FirstTouchHostCrsGraph<Ordinal,Node> &dgraph)
  {
    CrsMatrixHostCompute<Scalar,Ordinal,Node>::setOwnedGraph(dgraph);
    myFTGraph_     = &dgraph;
    staticFTGraph_ = &dgraph;
  }

  // ======= finalize ===========
  template <class Scalar, class Ordinal, class Node>
  void FirstTouchHostCrsMatrix<Scalar,Ordinal,Node>::finalize(bool OptimizeStorage)
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

} // end namespace Kokkos

#endif // KOKKOS_FIRSTTOUCH_SPARSEOPS_HPP
