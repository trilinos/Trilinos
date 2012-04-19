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

#ifndef KOKKOS_DEFAULTRELAXATION_HPP
#define KOKKOS_DEFAULTRELAXATION_HPP

#include <stdio.h>
#include <stdexcept>

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_Assert.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CrsMatrix.hpp" 
#include "Kokkos_CrsGraph.hpp" 
#include "Kokkos_MultiVector.hpp"
#include "Kokkos_NodeHelpers.hpp"
#include "Kokkos_DefaultArithmetic.hpp"
#include "Kokkos_DefaultRelaxationKernelOps.hpp"

namespace Kokkos {

  /*!
    \class DefaultRelaxation
    \brief Various relaxation methods.

    Methods include Jacobi, Gauss-Seidel, and Chebyshev polynomial relaxation.
  */

  template <class Scalar, class Ordinal, class Node = DefaultNode::DefaultNodeType>
  class DefaultRelaxation {
  public:
    typedef Scalar  ScalarType;
    typedef Ordinal OrdinalType;
    typedef Node    NodeType;

    //! @name Constructors/Destructor

    //@{

    //! DefaultRelaxation constructor 
    DefaultRelaxation(const RCP<Node> &node = DefaultNode::getDefaultNode());

    //! DefaultRelaxation Destructor
    ~DefaultRelaxation();

    //@}

    //! @name Accessor routines.
    //@{ 
    
    //! Node accessor.
    RCP<Node> getNode() const;

    //@}

    //! @name Initialization of structure

    //@{

    /*! Initialize structure of matrix.

       @todo Not implemented for general sparse graphs.
    */
    template <class GRAPH>
    void initializeStructure(GRAPH &graph, Teuchos::DataAccess cv);

    /*! Initialize values of matrix.
    
       @todo Not implemented for general sparse matrices
    */
    template <class MATRIX>
    void initializeValues(MATRIX &matrix, Teuchos::DataAccess cv);

    //! Initialize structure of the matrix, using CrsGraph
    template <class SparseOps>
    void initializeStructure(CrsGraph<Ordinal,Node,SparseOps > &graph, Teuchos::DataAccess cv);

    //! Initialize values of the matrix, using CrsMatrix
    template <class SparseOps>
    void initializeValues(CrsMatrix<Scalar,Ordinal,Node,SparseOps > &matrix, Teuchos::DataAccess cv);

    /*! Sets the diagonal inverted for relaxation using a MultiVector

      @todo Not implemented yet.
    */
    void setDiagonal(MultiVector<Scalar,Node> & diag);    

    //! Clear all matrix structures and values.
    void clear();

    //! 

    //@}

    //! @name Computational methods

    //@{

    //! Applies a sweep of Jacobi
    void sweep_jacobi(Scalar dampingFactor_, MultiVector<Scalar,Node> &X, const MultiVector<Scalar,Node> &B) const;

#ifdef ENABLE_ALL_OTHER_RELAXATION
    //! Applies a sweep of fine-grain Hybrid Gauss-Seidel
    void sweep_fine_hybrid(Scalar dampingFactor_, MultiVector<Scalar,Node> &X, const MultiVector<Scalar,Node> &B) const;

    //! Applies a sweep of coarse-grain Hybrid Gauss-Seidel
    void sweep_coarse_hybrid(Scalar dampingFactor_, size_t num_chunks, MultiVector<Scalar,Node> &X, const MultiVector<Scalar,Node> &B) const;

#endif //ifdef ENABLE_ALL_OTHER_RELAXATION
    //! Does setup for Chebyshev
    void setup_chebyshev(const Scalar lambda_max, const Scalar lambda_min);

    //! Applies a sweep of Chebyshev iteration
    void sweep_chebyshev(MultiVector<Scalar,Node> &X, const MultiVector<Scalar,Node> &B) const;

    //@}

  protected:
    //! Copy constructor (protected and unimplemented)
    DefaultRelaxation(const DefaultRelaxation& source);

    //! Extract the diagonal from the matrix, if the user hasn't set it already.
    // NTS: This should eventually disappear into some other Kokkos class...
    void ExtractDiagonal();    

    //! Update the Jacobi temporary vector size.
    bool UpdateJacobiTemp(size_t num_vectors, size_t vec_leng) const;

    //! Update the Chebyshev temporary vector size.
    bool UpdateChebyTemp(size_t num_vectors, size_t vec_leng) const;

    //! My node
    RCP<Node> node_;

    // we do this one of two ways: 
    // 1D/packed: array of offsets, pointer for ordinals, pointer for values. obviously the smallest footprint.
    ArrayRCP<const Ordinal> pbuf_inds1D_;
    ArrayRCP<const size_t>  begs1D_, ends1D_;
    ArrayRCP<Scalar>  pbuf_vals1D_;
    // 2D: array of pointers
    ArrayRCP<const ArrayRCP<Ordinal> > pbuf_inds2D_;
    ArrayRCP<const ArrayRCP<Scalar> > pbuf_vals2D_;
    ArrayRCP<const size_t>          pbuf_numEntries_;
    ArrayRCP<const Ordinal *> indPtrs_;
    ArrayRCP<Scalar  *> valPtrs_;
    
    //! Array containing matrix diagonal for easy access
    ArrayRCP<Scalar> diagonal_;

    //! Temporary work storage for Jacobi
    mutable ArrayRCP<Scalar> tempJacobiVector_;
    mutable size_t lastNumJacobiVectors_;

    // Arrays containing temp vectors for Chebyshev
    mutable ArrayRCP<Scalar> tempChebyVectorX_;
    mutable ArrayRCP<Scalar> tempChebyVectorW_;
    mutable size_t lastNumChebyVectors_;
    mutable bool cheby_setup_done_;
    mutable bool first_cheby_iteration_;

    // Constants for Chebyshev
    mutable Scalar lmin_,lmax_,delta_,s1_,oneOverTheta_,rho_,rho_new_,dtemp1_,dtemp2_;

    size_t numRows_;
    bool indsInit_, valsInit_, isEmpty_;
  }; //class DefaultRelaxation declaration


  /**********************************************************************/
  template<class Scalar, class Ordinal, class Node>
  DefaultRelaxation<Scalar,Ordinal,Node>::DefaultRelaxation(const RCP<Node> &node)
  : node_(node)
  , lastNumJacobiVectors_(0)
  , lastNumChebyVectors_(0)
  , cheby_setup_done_(false)
  , first_cheby_iteration_(false)  
  , indsInit_(false)
  , valsInit_(false)
  , isEmpty_(false)
  {
    lmin_=lmax_=delta_=s1_=oneOverTheta_=rho_=rho_new_=dtemp1_=dtemp2_=Teuchos::ScalarTraits<Scalar>::zero();
  }

  /**********************************************************************/
  template<class Scalar, class Ordinal, class Node>
  DefaultRelaxation<Scalar,Ordinal,Node>::~DefaultRelaxation() {
    clear();
  }

  /**********************************************************************/
  template <class Scalar, class Ordinal, class Node>
  template <class GRAPH>
  void DefaultRelaxation<Scalar,Ordinal,Node>::initializeStructure(GRAPH &graph, Teuchos::DataAccess cv) {
    // not implemented for general sparse graphs
    TEUCHOS_TEST_FOR_EXCEPT(true);
  }

  /**********************************************************************/
  template <class Scalar, class Ordinal, class Node>
  template <class MATRIX>
  void DefaultRelaxation<Scalar,Ordinal,Node>::initializeValues(MATRIX &graph, Teuchos::DataAccess cv) {
    // not implemented for general sparse matrices
    TEUCHOS_TEST_FOR_EXCEPT(true);
  }

  /**********************************************************************/
  template <class Scalar, class Ordinal, class Node>
  template <class SparseOps>
  void DefaultRelaxation<Scalar,Ordinal,Node>::initializeStructure(CrsGraph<Ordinal,Node,SparseOps > &graph, Teuchos::DataAccess cv) {
    using Teuchos::arcp;
    TEUCHOS_TEST_FOR_EXCEPTION(indsInit_ == true || valsInit_ == true, std::runtime_error, Teuchos::typeName(*this) << "::initializeStructure(): structure already initialized.");
    numRows_ = graph.getNumRows();
    if (graph.isEmpty() || numRows_ == 0) {
      isEmpty_ = true;
    }
    else if (graph.is1DStructure()) {
      isEmpty_ = false;
      ArrayRCP<Ordinal> inds;
      ArrayRCP<size_t> begs, ends;
      const_cast<CrsGraph<Ordinal,Node,SparseOps > &>(graph).get1DStructure( inds, begs, ends );
      pbuf_inds1D_ = inds;
      begs1D_ = begs;
      ends1D_ = ends;
    }
    else {
      isEmpty_  = false;
      {
        ArrayRCP<ArrayRCP<Ordinal> > inds;
        ArrayRCP<size_t>  sizes;
        const_cast<CrsGraph<Ordinal,Node,SparseOps > &>(graph).get2DStructure(inds,sizes);
        pbuf_inds2D_     = inds;
        pbuf_numEntries_ = sizes;
      }
      indPtrs_    = arcp<const Ordinal *>(numRows_);
      for (size_t r=0; r < numRows_; ++r) {
        indPtrs_[r] = pbuf_inds2D_[r].getRawPtr();
      }
    }
    indsInit_ = true;
  } //initializeStructure()

  /**********************************************************************/
  template <class Scalar, class Ordinal, class Node>
  template <class SparseOps>
  void DefaultRelaxation<Scalar,Ordinal,Node>::initializeValues(CrsMatrix<Scalar,Ordinal,Node,SparseOps > &matrix, Teuchos::DataAccess cv) {
    using Teuchos::arcp;
    TEUCHOS_TEST_FOR_EXCEPTION(indsInit_ == false, std::runtime_error, Teuchos::typeName(*this) << "::initializeValues(): must initialize values after graph.");
    TEUCHOS_TEST_FOR_EXCEPTION(numRows_ != matrix.getNumRows() || isEmpty_ != matrix.isEmpty() || (pbuf_inds2D_ != null && matrix.is1DStructure()) || (pbuf_inds1D_ != null && matrix.is2DStructure()), std::runtime_error, Teuchos::typeName(*this) << "::initializeValues(): matrix not compatible with previously supplied graph.");
    if (!isEmpty_) {        
      if (matrix.is1DStructure()) {
        ArrayRCP<Scalar> vals;
        const_cast<CrsMatrix<Scalar,Ordinal,Node,SparseOps > &>(matrix).get1DValues( vals );
        pbuf_vals1D_ = vals;
      }      
      else {
        {  
          ArrayRCP<ArrayRCP<Scalar> > vals;
          const_cast<CrsMatrix<Scalar,Ordinal,Node,SparseOps > &>(matrix).get2DValues(vals);
          pbuf_vals2D_ = vals;
        }      
        valPtrs_ = arcp<Scalar *>(numRows_);
        for (size_t r=0; r < numRows_; ++r) {
          valPtrs_[r] = pbuf_vals2D_[r].getRawPtr();
        }      
      }      
    }      
    valsInit_ = true;
    ExtractDiagonal();
  } //initializeValues()

  /**********************************************************************/
  template <class Scalar, class Ordinal, class Node>
  RCP<Node> DefaultRelaxation<Scalar,Ordinal,Node>::getNode() const { 
    return node_; 
  }

  /**********************************************************************/
  template <class Scalar, class Ordinal, class Node>
  void DefaultRelaxation<Scalar,Ordinal,Node>::clear() {
    pbuf_inds1D_      = null;
    pbuf_vals1D_      = null;
    pbuf_inds2D_      = null;
    pbuf_vals2D_      = null;
    pbuf_numEntries_  = null;
    diagonal_         = null;
    tempJacobiVector_ = null;
    tempChebyVectorW_ = null;
    tempChebyVectorX_ = null;
    indsInit_ = false;
    valsInit_ = false;
    isEmpty_  = false;
    lastNumJacobiVectors_=0;
    lastNumChebyVectors_=0;
  }

  /**********************************************************************/
  // Sets the diagonal inverted for relaxation using a MultiVector
  template <class Scalar, class Ordinal, class Node>
  void DefaultRelaxation<Scalar,Ordinal,Node>::setDiagonal(MultiVector<Scalar,Node> & diag){
    // Allocate space for diagonal
    diagonal_ = node_->template allocBuffer<Scalar>(numRows_);    

    // NTS: Copy diag over

    // Make it fail for now...
    TEUCHOS_TEST_FOR_EXCEPT(true);
  }
    
  /**********************************************************************/
  template <class Scalar, class Ordinal, class Node>
  void DefaultRelaxation<Scalar,Ordinal,Node>::ExtractDiagonal(){    
    TEUCHOS_TEST_FOR_EXCEPTION(valsInit_ == false, std::runtime_error, 
        Teuchos::typeName(*this) << "::ExtractDiagonal(): initializeValues() hasn't been called.");

    // Allocate space for diagonal
    diagonal_ = node_->template allocBuffer<Scalar>(numRows_);    
      
    if (pbuf_vals1D_ != null){
      // Extract the diagonal for Type 1 storage
      typedef ExtractDiagonalOp1<Scalar,Ordinal>  Op1D;
      ReadyBufferHelper<Node> rbh(node_);
      Op1D wdp;
      rbh.begin();
      wdp.numRows = numRows_;
      wdp.begs    = rbh.template addConstBuffer<size_t>(begs1D_);
      wdp.ends    = rbh.template addConstBuffer<size_t>(ends1D_);
      wdp.inds    = rbh.template addConstBuffer<Ordinal>(pbuf_inds1D_);
      wdp.vals    = rbh.template addConstBuffer<Scalar>(pbuf_vals1D_);
      wdp.diag    = rbh.template addNonConstBuffer<Scalar>(diagonal_);
      rbh.end();
      node_->template parallel_for<Op1D>(0,numRows_,wdp);
    }
    else {
      // Extract the diagonal for Type 2 storage
      typedef ExtractDiagonalOp2<Scalar,Ordinal>  Op2D;
      ReadyBufferHelper<Node> rbh(node_);
      Op2D wdp;
      rbh.begin();
      wdp.numRows = numRows_;
      wdp.numEntries = rbh.template addConstBuffer<size_t>(pbuf_numEntries_);
      wdp.inds_beg   = rbh.template addConstBuffer<const Ordinal *>(indPtrs_);
      wdp.vals_beg   = rbh.template addConstBuffer<Scalar *>(valPtrs_);
      wdp.diag    = rbh.template addNonConstBuffer<Scalar>(diagonal_);
      rbh.end();
      rbh.end();
      node_->template parallel_for<Op2D>(0,numRows_,wdp);
      
    }
  } //ExtractDiagonal()

  /**********************************************************************/
  template <class Scalar, class Ordinal, class Node>
  bool DefaultRelaxation<Scalar,Ordinal,Node>::UpdateJacobiTemp(size_t num_vectors, size_t vec_leng) const{
    // Re-allocate memory if needed
    if(num_vectors > lastNumJacobiVectors_){
      tempJacobiVector_= null;
      lastNumJacobiVectors_=num_vectors;
      tempJacobiVector_ = node_->template allocBuffer<Scalar>(vec_leng*lastNumJacobiVectors_); 
      return true;
    }
    return false;
  } //UpdateJacobiTemp()

  /**********************************************************************/
  template <class Scalar, class Ordinal, class Node>
  bool DefaultRelaxation<Scalar,Ordinal,Node>::UpdateChebyTemp(size_t num_vectors, size_t vec_leng) const{
    // Re-allocate memory if needed
    if(num_vectors > lastNumChebyVectors_){
      tempChebyVectorW_= null;
      tempChebyVectorX_= null;
      lastNumChebyVectors_=num_vectors;
      tempChebyVectorW_ = node_->template allocBuffer<Scalar>(numRows_*lastNumChebyVectors_);
      tempChebyVectorX_ = node_->template allocBuffer<Scalar>(vec_leng*lastNumChebyVectors_); 
      return true;
    }
    return false;
  } //UpdateChebyTemp()


#ifdef ENABLE_ALL_OTHER_RELAXATION
  /**********************************************************************/
  template <class Scalar, class Ordinal, class Node>
  void DefaultRelaxation<Scalar,Ordinal,Node>::sweep_fine_hybrid(Scalar dampingFactor_,
                         MultiVector<Scalar,Node> &X, const MultiVector<Scalar,Node> &B) const{
    typedef DefaultFineGrainHybridGaussSeidelOp1<Scalar,Ordinal>  Op1D;
    typedef DefaultFineGrainHybridGaussSeidelOp2<Scalar,Ordinal>  Op2D;


    TEUCHOS_TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::sweep_fine_hybrid(): operation not fully initialized.");
    TEUCHOS_TEST_FOR_EXCEPT(X.getNumCols() != B.getNumCols());
    ReadyBufferHelper<Node> rbh(node_);
    if (isEmpty_ == true) {
      // This makes no sense to try to call ...
      TEUCHOS_TEST_FOR_EXCEPT(true);
    }
    else if (begs1D_ != null) {
      Op1D wdp;
      rbh.begin();
      wdp.numRows = numRows_;
      wdp.inds    = rbh.template addConstBuffer<Ordinal>(pbuf_inds1D_);
      wdp.vals    = rbh.template addConstBuffer<Scalar>(pbuf_vals1D_);
      wdp.x       = rbh.template addNonConstBuffer<Scalar>(X.getValuesNonConst());
      wdp.b       = rbh.template addConstBuffer<Scalar>(B.getValues());
      wdp.diag    = rbh.template addConstBuffer<Scalar>(diagonal_);
      wdp.damping_factor = dampingFactor_;
      wdp.xstride = X.getStride();
      wdp.bstride = B.getStride();
      rbh.end();
      const size_t numRHS = X.getNumCols();
      node_->template parallel_for<Op1D>(0,numRows_*numRHS,wdp);
    }
    else {
      Op2D wdp;
      rbh.begin();
      wdp.numRows = numRows_;
      wdp.numEntries = rbh.template addConstBuffer<size_t>(pbuf_numEntries_);
      wdp.inds_beg   = rbh.template addConstBuffer<const Ordinal *>(pbuf_inds2D_);
      wdp.vals_beg   = rbh.template addConstBuffer<const Scalar *>(pbuf_vals2D_);
      wdp.x       = rbh.template addNonConstBuffer<Scalar>(X.getValuesNonConst());
      wdp.b       = rbh.template addConstBuffer<Scalar>(B.getValues());
      wdp.diag    = rbh.template addConstBuffer<Scalar>(diagonal_);
      wdp.damping_factor = dampingFactor_;
      wdp.xstride = X.getStride();
      wdp.bstride = B.getStride();
      rbh.end();
      const size_t numRHS = X.getNumCols();
      node_->template parallel_for<Op2D>(0,numRows_*numRHS,wdp);
    }
    return;
  } //sweep_fine_hybrid()

  /**********************************************************************/
  template <class Scalar, class Ordinal, class Node>
  void DefaultRelaxation<Scalar,Ordinal,Node>::sweep_coarse_hybrid(Scalar dampingFactor_,size_t num_chunks,
                         MultiVector<Scalar,Node> &X, const MultiVector<Scalar,Node> &B) const{
    typedef DefaultCoarseGrainHybridGaussSeidelOp1<Scalar,Ordinal>  Op1D;
    typedef DefaultCoarseGrainHybridGaussSeidelOp2<Scalar,Ordinal>  Op2D;


    TEUCHOS_TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::sweep_coarse_hybrid(): operation not fully initialized.");
    TEUCHOS_TEST_FOR_EXCEPT(X.getNumCols() != B.getNumCols());
    ReadyBufferHelper<Node> rbh(node_);   

    if (isEmpty_ == true) {
      // This makes no sense to try to call ...
      TEUCHOS_TEST_FOR_EXCEPT(true);
    }
    else if (begs1D_ != null) {
      Op1D wdp;
      rbh.begin();
      wdp.numRows = numRows_;
      wdp.numChunks = num_chunks;
      wdp.inds    = rbh.template addConstBuffer<Ordinal>(pbuf_inds1D_);
      wdp.vals    = rbh.template addConstBuffer<Scalar>(pbuf_vals1D_);
      wdp.x       = rbh.template addNonConstBuffer<Scalar>(X.getValuesNonConst());
      wdp.b       = rbh.template addConstBuffer<Scalar>(B.getValues());
      wdp.diag    = rbh.template addConstBuffer<Scalar>(diagonal_);
      wdp.damping_factor = dampingFactor_;
      wdp.xstride = X.getStride();
      wdp.bstride = B.getStride();
      rbh.end();
      const size_t numRHS = X.getNumCols();
      node_->template parallel_for<Op1D>(0,num_chunks*numRHS,wdp);
    }
    else {
      Op2D wdp;
      rbh.begin();
      wdp.numRows = numRows_;
      wdp.numChunks = num_chunks;
      wdp.numEntries = rbh.template addConstBuffer<size_t>(pbuf_numEntries_);
      wdp.inds_beg   = rbh.template addConstBuffer<const Ordinal *>(pbuf_inds2D_);
      wdp.vals_beg   = rbh.template addConstBuffer<const Scalar *>(pbuf_vals2D_);
      wdp.x       = rbh.template addNonConstBuffer<Scalar>(X.getValuesNonConst());
      wdp.b       = rbh.template addConstBuffer<Scalar>(B.getValues());
      wdp.diag    = rbh.template addConstBuffer<Scalar>(diagonal_);
      wdp.damping_factor = dampingFactor_;
      wdp.xstride = X.getStride();
      wdp.bstride = B.getStride();
      rbh.end();
      const size_t numRHS = X.getNumCols();
      node_->template parallel_for<Op2D>(0,num_chunks*numRHS,wdp);
    }
    return;
  } //sweep_coarse_hybrid()
#endif //ifdef ENABLE_ALL_OTHER_RELAXATION

  
  /********************************************************************/
  template <class Scalar, class Ordinal, class Node>
  void DefaultRelaxation<Scalar,Ordinal,Node>::sweep_jacobi(Scalar dampingFactor_,
                         MultiVector<Scalar,Node> &X, const MultiVector<Scalar,Node> &B) const{
    typedef DefaultJacobiOp1<Scalar,Ordinal>  Op1D;
    typedef DefaultJacobiOp2<Scalar,Ordinal>  Op2D;


    TEUCHOS_TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::sweep_jacobi(): operation not fully initialized.");
    TEUCHOS_TEST_FOR_EXCEPT(X.getNumCols() != B.getNumCols());
    ReadyBufferHelper<Node> rbh(node_);

    // Copy x over to the temp vector
    // NTS: The MultiVector copy constructor is a View. We need to do this the hard way.
    MultiVector<Scalar,Node> X0(X.getNode());
    UpdateJacobiTemp(X.getNumCols(),X.getNumRows());
    X0.initializeValues(X.getNumRows(),B.getNumCols(),tempJacobiVector_,B.getStride());
    DefaultArithmetic<MultiVector<Scalar,Node> >::Assign(X0,X);

    if (isEmpty_ == true) {
      // This makes no sense to try to call ...
      TEUCHOS_TEST_FOR_EXCEPT(true);
    }
    else if (begs1D_ != null) {
      Op1D wdp;
      rbh.begin();
      wdp.numRows = numRows_;
      wdp.begs    = rbh.template addConstBuffer<size_t>(begs1D_);
      wdp.ends    = rbh.template addConstBuffer<size_t>(ends1D_);
      wdp.inds    = rbh.template addConstBuffer<Ordinal>(pbuf_inds1D_);
      wdp.vals    = rbh.template addConstBuffer<Scalar>(pbuf_vals1D_);
      wdp.diag    = rbh.template addConstBuffer<Scalar>(diagonal_);
      wdp.x       = rbh.template addNonConstBuffer<Scalar>(X.getValuesNonConst());
      wdp.x0      = rbh.template addConstBuffer<Scalar>(X0.getValues());
      wdp.b       = rbh.template addConstBuffer<Scalar>(B.getValues());
      wdp.damping_factor = dampingFactor_;
      wdp.xstride = X.getStride();
      wdp.bstride = B.getStride();
      rbh.end();
      const size_t numRHS = X.getNumCols();
      node_->template parallel_for<Op1D>(0,numRows_*numRHS,wdp);
    }
    else {
      Op2D wdp;
      rbh.begin();
      wdp.numRows = numRows_;
      wdp.numEntries = rbh.template addConstBuffer<size_t>(pbuf_numEntries_);
      wdp.inds_beg   = rbh.template addConstBuffer<const Ordinal *>(indPtrs_);
      wdp.vals_beg   = rbh.template addConstBuffer<Scalar *>(valPtrs_);
      wdp.x       = rbh.template addNonConstBuffer<Scalar>(X.getValuesNonConst());
      wdp.x0      = rbh.template addConstBuffer<Scalar>(X0.getValues());
      wdp.b       = rbh.template addConstBuffer<Scalar>(B.getValues());
      wdp.diag    = rbh.template addConstBuffer<Scalar>(diagonal_);
      wdp.damping_factor = dampingFactor_;
      wdp.xstride = X.getStride();
      wdp.bstride = B.getStride();
      rbh.end();
      const size_t numRHS = X.getNumCols();
      node_->template parallel_for<Op2D>(0,numRows_*numRHS,wdp);
    }
    return;
  } //sweep_jacobi()

  /********************************************************************/
  template <class Scalar, class Ordinal, class Node>
  void DefaultRelaxation<Scalar,Ordinal,Node>::setup_chebyshev(const Scalar lambda_max, const Scalar lambda_min){
    lmax_=lambda_max;
    lmin_=lambda_min;

    Scalar alpha = lmin_;
    Scalar beta  = 1.1*lmax_;
    delta_ = 2.0 / (beta-alpha);
    Scalar theta_ = 0.5 * (beta+alpha);
    s1_    = theta_ * delta_;
    oneOverTheta_=1.0 / theta_;
    rho_=1.0/s1_;

    std::cout << std::endl;
    std::cout << "alpha  = " << alpha << std::endl;
    std::cout << "beta   = " << beta << std::endl;
    std::cout << "delta  = " << delta_ << std::endl;
    std::cout << "theta  = " << theta_ << std::endl;
    std::cout << "s1     = " << s1_ << std::endl;

    first_cheby_iteration_=true;
    cheby_setup_done_=true;
  } //setup_chebyshev()


  /********************************************************************/
  template <class Scalar, class Ordinal, class Node>
  void DefaultRelaxation<Scalar,Ordinal,Node>::sweep_chebyshev(MultiVector<Scalar,Node> &X, const MultiVector<Scalar,Node> &B) const{
    typedef DefaultChebyshevOp1<Scalar,Ordinal>  Op1D;
    //    typedef DefaultChebyOp2<Scalar,Ordinal>  Op2D;
    static int iterNumber=1;
    Scalar xnorm;

    if (first_cheby_iteration_) {
      xnorm=DefaultArithmetic<MultiVector<Scalar,Node> >::Norm2Squared(X);
      std::cout << "||x_0" << "|| = " << sqrt(xnorm) << std::endl;
    }
    
    TEUCHOS_TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
               Teuchos::typeName(*this) << "::sweep_chebyshev(): operation not fully initialized.");
    TEUCHOS_TEST_FOR_EXCEPT(X.getNumCols() != B.getNumCols());
    //size_t xstride = X.getStride();
    //size_t bstride = B.getStride();
    //TEUCHOS_TEST_FOR_EXCEPT(xstride != bstride); //TODO JJH don't think we want this b/c X is imported, B is local
    TEUCHOS_TEST_FOR_EXCEPT(cheby_setup_done_ == false);
    
    ReadyBufferHelper<Node> rbh(node_);
    
    // If the number of multivectors has changed, we need to run Cheby from scratch
    first_cheby_iteration_= UpdateChebyTemp(X.getNumCols(),X.getNumRows()) || first_cheby_iteration_;
    
    // Copy X to X0 so we get the matvec correct
    MultiVector<Scalar,Node> X0(X.getNode());
    X0.initializeValues(X.getNumRows(),X.getNumCols(),tempChebyVectorX_,X.getStride());
    DefaultArithmetic<MultiVector<Scalar,Node> >::Assign(X0,X);

    // Update Scalars 
    if(!first_cheby_iteration_){
      rho_new_ = 1.0 / (2.0*s1_ - rho_);
      dtemp1_=rho_*rho_new_;
      dtemp2_=2*rho_new_*delta_;
      rho_=rho_new_;
    }    
    //printf("    rho = %6.4e rho_new = %6.4e dtemp1=%6.4e dtemp2=%6.4e\n",rho_,rho_new_,dtemp1_,dtemp2_);


    
    if (isEmpty_ == true) {
      // This makes no sense to try to call ...
      TEUCHOS_TEST_FOR_EXCEPT(true);
    }
    else if (begs1D_ != null) {
      Op1D wdp;
      rbh.begin();
      wdp.numRows = numRows_;
      wdp.begs    = rbh.template addConstBuffer<size_t>(begs1D_);
      wdp.ends    = rbh.template addConstBuffer<size_t>(ends1D_);
      wdp.inds    = rbh.template addConstBuffer<Ordinal>(pbuf_inds1D_);
      wdp.vals    = rbh.template addConstBuffer<Scalar>(pbuf_vals1D_);
      wdp.x       = rbh.template addNonConstBuffer<Scalar>(X.getValuesNonConst());
      wdp.x0      = rbh.template addConstBuffer<Scalar>(tempChebyVectorX_);
      wdp.w       = rbh.template addNonConstBuffer<Scalar>(tempChebyVectorW_);
      wdp.b       = rbh.template addConstBuffer<Scalar>(B.getValues());
      wdp.diag    = rbh.template addConstBuffer<Scalar>(diagonal_);
      wdp.stride = X.getStride();
      wdp.first_step=first_cheby_iteration_;
      wdp.zero_initial_guess=false;//HAQ Fix me later
      wdp.oneOverTheta = oneOverTheta_;
      wdp.dtemp1 = dtemp1_;
      wdp.dtemp2 = dtemp2_;
      
      rbh.end();    
      const size_t numRHS = X.getNumCols();
      node_->template parallel_for<Op1D>(0,numRows_*numRHS,wdp);
    }
    else {
      /*      Op2D wdp;
          rbh.begin();
          wdp.numRows = numRows_;
          wdp.numEntries = rbh.template addConstBuffer<size_t>(pbuf_numEntries_);
          wdp.inds_beg   = rbh.template addConstBuffer<const Ordinal *>(pbuf_inds2D_);
          wdp.vals_beg   = rbh.template addConstBuffer<const Scalar *>(pbuf_vals2D_);
          wdp.x       = rbh.template addNonConstBuffer<Scalar>(X.getValuesNonConst());
          wdp.x0      = rbh.template addConstBuffer<Scalar>(X0.getValues());
          wdp.b       = rbh.template addConstBuffer<Scalar>(B.getValues());
          wdp.diag    = rbh.template addConstBuffer<Scalar>(diagonal_);
          wdp.damping_factor = dampingFactor_;
          wdp.xstride = X.getStride();
          wdp.bstride = B.getStride();
          rbh.end();
          const size_t numRHS = X.getNumCols();
          node_->template parallel_for<Op2D>(0,numRows_*numRHS,wdp);
      */
    }
    xnorm=DefaultArithmetic<MultiVector<Scalar,Node> >::Norm2Squared(X);
    std::cout << "||x_" << iterNumber++ << "|| = " << sqrt(xnorm) << std::endl;
    
    first_cheby_iteration_=false;

  } //sweep_chebyshev()

} // namespace Kokkos

#endif // KOKKOS_DEFAULTRELAXATION_HPP
