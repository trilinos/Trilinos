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

#ifndef KOKKOS_DEFAULTRELAXATION_HPP
#define KOKKOS_DEFAULTRELAXATION_HPP

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <stdexcept>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CrsMatrix.hpp" 
#include "Kokkos_CrsGraph.hpp" 
#include "Kokkos_MultiVector.hpp"
#include "Kokkos_NodeHelpers.hpp"
#include "Kokkos_DefaultArithmetic.hpp"
#include "Kokkos_DefaultRelaxationKernelOps.hpp"

#include <stdio.h>

namespace Kokkos {

  template <class Scalar, class Ordinal, class Node = DefaultNode::DefaultNodeType>
  class DefaultRelaxation {
  public:
    typedef Scalar  ScalarType;
    typedef Ordinal OrdinalType;
    typedef Node    NodeType;

    //! @name Constructors/Destructor

    //@{

    //! DefaultRelaxation constuctor 
    DefaultRelaxation(const Teuchos::RCP<Node> &node = DefaultNode::getDefaultNode());

    //! DefaultRelaxation Destructor
    ~DefaultRelaxation();

    //@}

    //! @name Accessor routines.
    //@{ 
    
    //! Node accessor.
    Teuchos::RCP<Node> getNode() const;

    //@}

    //! @name Initialization of structure

    //@{

    //! Initialize structure of matrix. NOT IMPLEMENTED!
    template <class GRAPH>
    Teuchos::DataAccess initializeStructure(GRAPH &graph, Teuchos::DataAccess cv);

    //! Initialize values of matrix.  NOT IMPLEMENTED!
    template <class MATRIX>
    Teuchos::DataAccess initializeValues(MATRIX &matrix, Teuchos::DataAccess cv);

    //! Initialize structure of the matrix, using Kokkos::CrsGraph
    Teuchos::DataAccess initializeStructure(CrsGraph<Ordinal,Node> &graph, Teuchos::DataAccess cv);

    //! Initialize values of the matrix, using Kokkos::CrsMatrix
    Teuchos::DataAccess initializeValues(CrsMatrix<Scalar,Node> &matrix, Teuchos::DataAccess cv);

    //! Sets the diagonal inverted for relaxation using a Kokkos::MultiVector
    void setDiagonal(MultiVector<Scalar,Node> & diag);    

    //! Clear all matrix structure and values.
    void clear();

    //! 

    //@}

    //! @name Computational methods

    //@{

    //! Applies a sweep of Jacobi
    void sweep_jacobi(Scalar dampingFactor_, MultiVector<Scalar,Node> &X, const MultiVector<Scalar,Node> &B) const;

    //! Applies a sweep of fine-grain Hybrid Gauss-Seidel
    void sweep_fine_hybrid(Scalar dampingFactor_, MultiVector<Scalar,Node> &X, const MultiVector<Scalar,Node> &B) const;

    //! Applies a sweep of coarse-grain Hybrid Gauss-Seidel
    void sweep_coarse_hybrid(Scalar dampingFactor_, size_t num_chunks, MultiVector<Scalar,Node> &X, const MultiVector<Scalar,Node> &B) const;

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

    // Update the temp vector size
    bool UpdateJacobiTemp(size_t num_vectors) const;
    bool UpdateChebyTemp(size_t num_vectors) const;

    // My node
    Teuchos::RCP<Node> node_;

    // we do this one of two ways: 
    // 1D/packed: array of offsets, pointer for ordinals, pointer for values. obviously the smallest footprint.
    Teuchos::ArrayRCP<const Ordinal> pbuf_inds1D_;
    Teuchos::ArrayRCP<const size_t>  pbuf_offsets1D_;
    Teuchos::ArrayRCP<const Scalar>  pbuf_vals1D_;
    // 2D: array of pointers
    Teuchos::ArrayRCP<const Ordinal *> pbuf_inds2D_;
    Teuchos::ArrayRCP<const Scalar  *> pbuf_vals2D_;
    Teuchos::ArrayRCP<size_t>          pbuf_numEntries_;
    
    // Array containing matrix diagonal for easy access
    Teuchos::ArrayRCP<Scalar> diagonal_;

    // Array containing a temp vector for Jacobi
    mutable Teuchos::ArrayRCP<Scalar> tempJacobiVector_;
    mutable size_t lastNumJacobiVectors_;

    // Arrays containing temp vectors for Chebyshev
    mutable Teuchos::ArrayRCP<Scalar> tempChebyVectorX_;
    mutable Teuchos::ArrayRCP<Scalar> tempChebyVectorW_;
    mutable size_t lastNumChebyVectors_;
    mutable bool cheby_setup_done_;
    mutable bool first_cheby_iteration_;

    // Constants for Chebyshev
    mutable Scalar lmin_,lmax_,delta_,s1_,oneOverTheta_,rho_,rho_new_,dtemp1_,dtemp2_;

    size_t numRows_;
    bool indsInit_, valsInit_, isPacked_, isEmpty_;
  };


  /**********************************************************************/
  template<class Scalar, class Ordinal, class Node>
  DefaultRelaxation<Scalar,Ordinal,Node>::DefaultRelaxation(const Teuchos::RCP<Node> &node)
  : node_(node)
  , lastNumJacobiVectors_(0)
  , lastNumChebyVectors_(0)
  , cheby_setup_done_(false)
  , first_cheby_iteration_(false)  
  , indsInit_(false)
  , valsInit_(false)
  , isPacked_(false) 
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
  template<class Scalar, class Ordinal, class Node>
  template <class GRAPH>
  Teuchos::DataAccess DefaultRelaxation<Scalar,Ordinal,Node>::initializeStructure(GRAPH &graph, Teuchos::DataAccess cv) {
    // not implemented for general sparse graphs
    TEST_FOR_EXCEPT(true);
  }

  /**********************************************************************/
  template<class Scalar, class Ordinal, class Node>
  template <class MATRIX>
  Teuchos::DataAccess DefaultRelaxation<Scalar,Ordinal,Node>::initializeValues(MATRIX &graph, Teuchos::DataAccess cv) {
    // not implemented for general sparse matrices
    TEST_FOR_EXCEPT(true);
  }

  /**********************************************************************/
  template <class Scalar, class Ordinal, class Node>
  Teuchos::DataAccess DefaultRelaxation<Scalar,Ordinal,Node>::initializeStructure(CrsGraph<Ordinal,Node> &graph, Teuchos::DataAccess cv) {
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(cv != Teuchos::View, std::runtime_error,
        Teuchos::typeName(*this) << "::initializeStructure(): requires View access.");
    TEST_FOR_EXCEPTION(indsInit_ == true, std::runtime_error, 
        Teuchos::typeName(*this) << "::initializeStructure(): structure already initialized.");
    numRows_ = graph.getNumRows();
    if (graph.isEmpty() || numRows_ == 0) {
      isEmpty_ = true;
    }
    else if (graph.isPacked()) {
      isEmpty_ = false;
      isPacked_ = true;
      pbuf_inds1D_    = graph.getPackedIndices();
      pbuf_offsets1D_ = graph.getPackedOffsets();
    }
    else {
      isEmpty_ = false;
      isPacked_ = false;
      pbuf_inds2D_     = node_->template allocBuffer<const Ordinal *>(numRows_);
      pbuf_numEntries_ = node_->template allocBuffer<size_t>(numRows_);
      ArrayRCP<const Ordinal *> inds2Dview = node_->template viewBufferNonConst<const Ordinal *>(WriteOnly, numRows_, pbuf_inds2D_);
      ArrayRCP<         size_t> numEntview = node_->template viewBufferNonConst<         size_t>(WriteOnly, numRows_, pbuf_numEntries_);
      for (size_t r=0; r < numRows_; ++r) {
        ArrayRCP<const Ordinal> rowinds = graph.get2DIndices(r);
        if (rowinds != Teuchos::null) {
          inds2Dview[r] = rowinds.getRawPtr();
          numEntview[r] = rowinds.size();
        }
        else {
          inds2Dview[r] = NULL;
          numEntview[r] = 0;
        }
      }
    }
    indsInit_ = true;
    return Teuchos::View;
  }

  /**********************************************************************/
  template <class Scalar, class Ordinal, class Node>
  Teuchos::DataAccess DefaultRelaxation<Scalar,Ordinal,Node>::initializeValues(CrsMatrix<Scalar,Node> &matrix, Teuchos::DataAccess cv) {
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(cv != Teuchos::View, std::runtime_error,
        Teuchos::typeName(*this) << "::initializeValues(): requires View access.");
    TEST_FOR_EXCEPTION(valsInit_ == true, std::runtime_error, 
        Teuchos::typeName(*this) << "::initializeValues(): values already initialized.");
    TEST_FOR_EXCEPTION(numRows_ != matrix.getNumRows() || (!isEmpty_ && isPacked_ != matrix.isPacked()), std::runtime_error,
        Teuchos::typeName(*this) << "::initializeValues(): matrix not compatible with previously supplied graph.");
    if (isEmpty_ || matrix.isEmpty() || numRows_ == 0) {
      isEmpty_ = true;
    }
    else if (matrix.isPacked()) {
      isEmpty_ = false;
      pbuf_vals1D_ = matrix.getPackedValues();
    }
    else {
      isEmpty_ = false;
      pbuf_vals2D_ = node_->template allocBuffer<const Scalar *>(numRows_);
      ArrayRCP<const Scalar *> vals2Dview = node_->template viewBufferNonConst<const Scalar *>(WriteOnly, numRows_, pbuf_vals2D_);
      for (size_t r=0; r < numRows_; ++r) {
        ArrayRCP<const Scalar> rowvals = matrix.get2DValues(r);
        if (rowvals != Teuchos::null) {
          vals2Dview[r] = rowvals.getRawPtr();
        }
        else {
          vals2Dview[r] = NULL;
        }
      }
    }
    
    valsInit_ = true;
    ExtractDiagonal();
    return Teuchos::View;
  }

  /**********************************************************************/
  template <class Scalar, class Ordinal, class Node>
  Teuchos::RCP<Node> DefaultRelaxation<Scalar,Ordinal,Node>::getNode() const { 
    return node_; 
  }

  /**********************************************************************/
  template <class Scalar, class Ordinal, class Node>
  void DefaultRelaxation<Scalar,Ordinal,Node>::clear() {
    pbuf_inds1D_      = Teuchos::null;
    pbuf_offsets1D_   = Teuchos::null;
    pbuf_vals1D_      = Teuchos::null;
    pbuf_inds2D_      = Teuchos::null;
    pbuf_vals2D_      = Teuchos::null;
    pbuf_numEntries_  = Teuchos::null;
    diagonal_         = Teuchos::null;
    tempJacobiVector_ = Teuchos::null;
    tempChebyVectorW_ = Teuchos::null;
    tempChebyVectorX_ = Teuchos::null;
    indsInit_ = false;
    valsInit_ = false;
    isPacked_ = false;
    isEmpty_  = false;
    lastNumJacobiVectors_=0;
    lastNumChebyVectors_=0;
  }

  /**********************************************************************/
  // Sets the diagonal inverted for relaxation using a Kokkos::MultiVector
  template <class Scalar, class Ordinal, class Node>
  void DefaultRelaxation<Scalar,Ordinal,Node>::setDiagonal(MultiVector<Scalar,Node> & diag){
    // Allocate space for diagonal
    diagonal_ = node_->template allocBuffer<Scalar>(numRows_);    

    // NTS: Copy diag over

    // Make it fail for now...
    TEST_FOR_EXCEPT(true);
  }
    
  /**********************************************************************/
  template <class Scalar, class Ordinal, class Node>
  void DefaultRelaxation<Scalar,Ordinal,Node>::ExtractDiagonal(){    
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(valsInit_ == false, std::runtime_error, 
        Teuchos::typeName(*this) << "::ExtractDiagonal(): initializeValues() hasn't been called.");

    // Allocate space for diagonal
    diagonal_ = node_->template allocBuffer<Scalar>(numRows_);    
      
    if (pbuf_vals1D_ != Teuchos::null){
      // Extract the diagonal for Type 1 storage
      typedef ExtractDiagonalOp1<Scalar,Ordinal>  Op1D;
      ReadyBufferHelper<Node> rbh(node_);
      Op1D wdp;
      rbh.begin();
      wdp.numRows = numRows_;
      wdp.offsets = rbh.template addConstBuffer<size_t>(pbuf_offsets1D_);
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
      wdp.inds_beg   = rbh.template addConstBuffer<const Ordinal *>(pbuf_inds2D_);
      wdp.vals_beg   = rbh.template addConstBuffer<const Scalar *>(pbuf_vals2D_);
      wdp.diag    = rbh.template addNonConstBuffer<Scalar>(diagonal_);
      rbh.end();
      rbh.end();
      node_->template parallel_for<Op2D>(0,numRows_,wdp);
      
    }
  }
  /**********************************************************************/
  template <class Scalar, class Ordinal, class Node>
  bool DefaultRelaxation<Scalar,Ordinal,Node>::UpdateJacobiTemp(size_t num_vectors) const{
    // Re-allocate memory if needed
    if(num_vectors > lastNumJacobiVectors_){
      tempJacobiVector_=Teuchos::null;
      lastNumJacobiVectors_=num_vectors;
      tempJacobiVector_ = node_->template allocBuffer<Scalar>(numRows_*lastNumJacobiVectors_); 
      return true;
    }
    return false;
  }

  /**********************************************************************/
  template <class Scalar, class Ordinal, class Node>
  bool DefaultRelaxation<Scalar,Ordinal,Node>::UpdateChebyTemp(size_t num_vectors) const{
    // Re-allocate memory if needed
    if(num_vectors > lastNumChebyVectors_){
      tempChebyVectorW_=Teuchos::null;
      tempChebyVectorX_=Teuchos::null;
      lastNumChebyVectors_=num_vectors;
      tempChebyVectorW_ = node_->template allocBuffer<Scalar>(numRows_*lastNumChebyVectors_); 
      tempChebyVectorX_ = node_->template allocBuffer<Scalar>(numRows_*lastNumChebyVectors_); 
      return true;
    }
    return false;
  }


  /**********************************************************************/
  template <class Scalar, class Ordinal, class Node>
  void DefaultRelaxation<Scalar,Ordinal,Node>::sweep_fine_hybrid(Scalar dampingFactor_,
						 MultiVector<Scalar,Node> &X, const MultiVector<Scalar,Node> &B) const{
    typedef DefaultFineGrainHybridGaussSeidelOp1<Scalar,Ordinal>  Op1D;
    typedef DefaultFineGrainHybridGaussSeidelOp2<Scalar,Ordinal>  Op2D;


    TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::sweep_fine_hybrid(): operation not fully initialized.");
    TEST_FOR_EXCEPT(X.getNumCols() != B.getNumCols());
    ReadyBufferHelper<Node> rbh(node_);
    if (isEmpty_ == true) {
      // This makes no sense to try to call ...
      TEST_FOR_EXCEPT(true);
    }
    else if (isPacked_ == true) {
      Op1D wdp;
      rbh.begin();
      wdp.numRows = numRows_;
      wdp.offsets = rbh.template addConstBuffer<size_t>(pbuf_offsets1D_);
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
  }

  /**********************************************************************/
  template <class Scalar, class Ordinal, class Node>
  void DefaultRelaxation<Scalar,Ordinal,Node>::sweep_coarse_hybrid(Scalar dampingFactor_,size_t num_chunks,
						 MultiVector<Scalar,Node> &X, const MultiVector<Scalar,Node> &B) const{
    typedef DefaultCoarseGrainHybridGaussSeidelOp1<Scalar,Ordinal>  Op1D;
    typedef DefaultCoarseGrainHybridGaussSeidelOp2<Scalar,Ordinal>  Op2D;


    TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::sweep_fine_hybrid(): operation not fully initialized.");
    TEST_FOR_EXCEPT(X.getNumCols() != B.getNumCols());
    ReadyBufferHelper<Node> rbh(node_);   

    if (isEmpty_ == true) {
      // This makes no sense to try to call ...
      TEST_FOR_EXCEPT(true);
    }
    else if (isPacked_ == true) {
      Op1D wdp;
      rbh.begin();
      wdp.numRows = numRows_;
      wdp.numChunks = num_chunks;
      wdp.offsets = rbh.template addConstBuffer<size_t>(pbuf_offsets1D_);
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
  }

  
  /********************************************************************/
  template <class Scalar, class Ordinal, class Node>
  void DefaultRelaxation<Scalar,Ordinal,Node>::sweep_jacobi(Scalar dampingFactor_,
						 MultiVector<Scalar,Node> &X, const MultiVector<Scalar,Node> &B) const{
    typedef DefaultJacobiOp1<Scalar,Ordinal>  Op1D;
    typedef DefaultJacobiOp2<Scalar,Ordinal>  Op2D;


    TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::sweep_jacobi(): operation not fully initialized.");
    TEST_FOR_EXCEPT(X.getNumCols() != B.getNumCols());
    ReadyBufferHelper<Node> rbh(node_);

    // Copy x over to the temp vector
    // NTS: The MultiVector copy constructor is a View. We need to do this the hard way.
    MultiVector<Scalar,Node> X0(X.getNode());
    UpdateJacobiTemp(X.getNumCols());
    X0.initializeValues(numRows_,B.getNumCols(),tempJacobiVector_,B.getStride());
    DefaultArithmetic<MultiVector<Scalar,Node> >::Assign(X0,X);

    if (isEmpty_ == true) {
      // This makes no sense to try to call ...
      TEST_FOR_EXCEPT(true);
    }
    else if (isPacked_ == true) {
      Op1D wdp;
      rbh.begin();
      wdp.numRows = numRows_;
      wdp.offsets = rbh.template addConstBuffer<size_t>(pbuf_offsets1D_);
      wdp.inds    = rbh.template addConstBuffer<Ordinal>(pbuf_inds1D_);
      wdp.vals    = rbh.template addConstBuffer<Scalar>(pbuf_vals1D_);
      wdp.x       = rbh.template addNonConstBuffer<Scalar>(X.getValuesNonConst());
      wdp.x0      = rbh.template addConstBuffer<Scalar>(X0.getValues());
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
  }

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

    first_cheby_iteration_=true;
    cheby_setup_done_=true;
  }


  /********************************************************************/
  template <class Scalar, class Ordinal, class Node>
  void DefaultRelaxation<Scalar,Ordinal,Node>::sweep_chebyshev(MultiVector<Scalar,Node> &X, const MultiVector<Scalar,Node> &B) const{
    typedef DefaultChebyshevOp1<Scalar,Ordinal>  Op1D;
    //    typedef DefaultChebyOp2<Scalar,Ordinal>  Op2D;
    
    TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
		       Teuchos::typeName(*this) << "::sweep_jacobi(): operation not fully initialized.");
    TEST_FOR_EXCEPT(X.getNumCols() != B.getNumCols());
    TEST_FOR_EXCEPT(X.getStride() != B.getStride());
    TEST_FOR_EXCEPT(cheby_setup_done_ == false);
    
    ReadyBufferHelper<Node> rbh(node_);
    
    // If the number of multivectors has changed, we need to run Cheby from scratch
    first_cheby_iteration_= UpdateChebyTemp(X.getNumCols()) || first_cheby_iteration_;
    
    // Copy X to X0 so we get the matvec correct
    MultiVector<Scalar,Node> X0(X.getNode());
    X0.initializeValues(numRows_,X.getNumCols(),tempChebyVectorX_,X.getStride());
    DefaultArithmetic<MultiVector<Scalar,Node> >::Assign(X0,X);

    // Update Scalars 
    if(!first_cheby_iteration_){
      rho_new_ = 1.0 / (2.0*s1_ - rho_);
      dtemp1_=rho_*rho_new_;
      dtemp2_=2*rho_new_*delta_;
      rho_=rho_new_;
    }    
    //    printf("    rho = %6.4e rho_new = %6.4e dtemp1=%6.4e dtemp2=%6.4e\n",rho_,rho_new_,dtemp1_,dtemp2_);


    
    if (isEmpty_ == true) {
      // This makes no sense to try to call ...
      TEST_FOR_EXCEPT(true);
    }
    else if (isPacked_ == true) {
      Op1D wdp;
      rbh.begin();
      wdp.numRows = numRows_;
      wdp.offsets = rbh.template addConstBuffer<size_t>(pbuf_offsets1D_);
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
    
    first_cheby_iteration_=false;
    return;
  }

} // namespace Kokkos

#endif /* KOKKOS_DEFAULTRELAXATION_HPP */
