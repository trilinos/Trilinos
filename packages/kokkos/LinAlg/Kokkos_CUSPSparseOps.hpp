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

#ifndef KOKKOS_CUSPSPARSEOPS_HPP
#define KOKKOS_CUSPSPARSEOPS_HPP

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_CompileTimeAssert.hpp>
#include <stdexcept>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CrsMatrix.hpp" 
#include "Kokkos_CrsGraph.hpp" 
#include "Kokkos_MultiVector.hpp"
#include "Kokkos_DefaultArithmetic.hpp"

#include <cusp/hyb_matrix.h>
#include <cusp/csr_matrix.h>
#include <cusp/device/detail/spmv.h>

namespace Kokkos {

  /** \brief Access to sparse matrix multiply and solve using the CUSP project.
      \ingroup kokkos_crs_ops
    */
  template <class Scalar, class Ordinal, class Node>
  class CUSPSparseOps {
  public:
    //@{ 
    //! @name Typedefs and structs

    //!
    typedef Scalar  ScalarType;
    //!
    typedef Ordinal OrdinalType;
    //!
    typedef Node    NodeType;

    /** \brief Rebind struct, for specifying type information for a different scalar.
          
        This specifies a CUSPSparseOps object, regardless of scalar type.
      */
    template <class S2>
    struct rebind {
      typedef CUSPSparseOps<S2,Ordinal,Node> other;
    };

    //@}
    //! @name Constructors/Destructor
    //@{

    //! CUSPSparseOps constuctor with variable number of indices per row.
    CUSPSparseOps(const RCP<Node> &node);

    //! CUSPSparseOps Destructor
    ~CUSPSparseOps();

    //@}
    //! @name Accessor routines.
    //@{ 
    
    //! Node accessor.
    RCP<Node> getNode() const;

    //@}
    //! @name Initialization of structure
    //@{

    //! Initialize structure of matrix, using CrsGraphHostCompute
    void initializeStructure(const CrsGraphHostCompute<Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> > &graph);

    //! Initialize values of matrix, using CrsMatrixHostCompute
    void initializeValues(const CrsMatrixHostCompute<Scalar,Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> > &matrix);

    //! Clear all matrix structure and values.
    void clear();

    //@}
    //! @name Computational methods
    //@{

    //! Applies the matrix to a MultiVector, overwriting Y.
    template <class DomainScalar, class RangeScalar>
    void multiply(Teuchos::ETransp trans, RangeScalar alpha, const MultiVector<DomainScalar,Node> &X, MultiVector<RangeScalar,Node> &Y) const;

    //! Applies the matrix to a MultiVector, accumulating into Y.
    template <class DomainScalar, class RangeScalar>
    void multiply(Teuchos::ETransp trans, 
                  RangeScalar alpha, const MultiVector<DomainScalar,Node> &X, RangeScalar beta, MultiVector<RangeScalar,Node> &Y) const;

    //! Solves the matrix for a given set of right-hand-sides.
    template <class DomainScalar, class RangeScalar>
    void solve(Teuchos::ETransp trans, Teuchos::EUplo uplo, Teuchos::EDiag diag, 
               const MultiVector<DomainScalar,Node> &Y, MultiVector<RangeScalar,Node> &X) const;

    //@}

  protected:
    //! Copy constructor (protected and unimplemented)
    CUSPSparseOps(const CUSPSparseOps& source);

    RCP<Node> node_;
    RCP<cusp::csr_matrix<Ordinal,Scalar,cusp::host_memory  > > hostCSR_;
    RCP<cusp::hyb_matrix<Ordinal,Scalar,cusp::device_memory> > devcHYB_;

    // we do this one of two ways: 
    // 1D/packed: arrays of offsets, array of ordinals, array of values.
    ArrayRCP<size_t>  begs1D_, ends1D_;
    ArrayRCP<Ordinal> inds1D_;
    // 2D: array of pointers
    ArrayRCP<ArrayRCP<Ordinal> > inds2D_;

    size_t numRows_;
    size_t numNZ_; 
    bool isEmpty_;
    bool indsInit_, valsInit_;
  };

  template<class Scalar, class Ordinal, class Node>
  CUSPSparseOps<Scalar,Ordinal,Node>::CUSPSparseOps(const RCP<Node> &node)
  : node_(node) 
  {
    clear();
  }

  template<class Scalar, class Ordinal, class Node>
  CUSPSparseOps<Scalar,Ordinal,Node>::~CUSPSparseOps() {
  }

  template <class Scalar, class Ordinal, class Node>
  RCP<Node> CUSPSparseOps<Scalar,Ordinal,Node>::getNode() const {
    return node_; 
  }

  template <class Scalar, class Ordinal, class Node>
  void CUSPSparseOps<Scalar,Ordinal,Node>::clear() {
    isEmpty_ = false;
    numRows_ = 0;
    numNZ_   = 0;
    indsInit_ = 0;
    valsInit_ = 0;
    //
    hostCSR_ = null;
    devcHYB_ = null;
    //
    begs1D_ = null;
    ends1D_ = null;
    inds1D_ = null;
    inds2D_ = null;
  }

  template <class Scalar, class Ordinal, class Node>
  void CUSPSparseOps<Scalar,Ordinal,Node>::initializeStructure(const CrsGraphHostCompute<Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> > &graph) {
    TEST_FOR_EXCEPTION(indsInit_ == true || valsInit_ == true, std::runtime_error, 
        Teuchos::typeName(*this) << "::initializeStructure(): structure already initialized.");
    numRows_ = graph.getNumRows();
    if (graph.isEmpty() || numRows_ == 0) {
      isEmpty_ = true;
    }
    else if (graph.is1DStructure()) {
      const_cast<CrsGraphHostviceCompute<Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> > &>(graph).getDeviceBuffers(inds, offs);
    }
    else {
      hostCSR_ = rcp(new cusp::csr_matrix<Ordinal,Scalar,cusp::host_memory>(numRows_,numCols_,numNZ_));
      if (graph.is2DStructure()) {
      }
      else {
        TEST_FOR_EXCEPT(true);
      }
    }
      ArrayRCP<Ordinal> inds;
      ArrayRCP<size_t > offs;
      const_cast<CrsGraphHostCompute<Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> > &>(graph).getDeviceBuffers(inds, offs);
      pbuf_inds1D_    = inds;
      pbuf_offsets1D_ = offs;
    }
    indsInit_ = true;
  }

  template <class Scalar, class Ordinal, class Node>
  void CUSPSparseOps<Scalar,Ordinal,Node>::initializeValues(const CrsMatrixHostCompute<Scalar,Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> > &matrix) {
    TEST_FOR_EXCEPTION(indsInit_ == false, std::runtime_error, 
        Teuchos::typeName(*this) << "::initializeValues(): must initialize values after graph.");
    TEST_FOR_EXCEPTION(numRows_ != matrix.getNumRows() || isEmpty_ != matrix.isEmpty(), std::runtime_error,
        Teuchos::typeName(*this) << "::initializeValues(): matrix not compatible with previously supplied graph.");
    if (!isEmpty_) {
      ArrayRCP<Scalar> vals;
      const_cast<CrsMatrixHostCompute<Scalar,Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> > &>(matrix).getDeviceBuffer(vals);
      pbuf_vals1D_ = vals;
    }
    valsInit_ = true;
  }


  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void CUSPSparseOps<Scalar,Ordinal,Node>::solve(
                      Teuchos::ETransp trans, Teuchos::EUplo uplo, Teuchos::EDiag diag, 
                      const MultiVector<DomainScalar,Node> &Y,
                            MultiVector<RangeScalar,Node> &X) const {
    TEST_FOR_EXCEPTION(true, std::logic_error, 
      Teuchos::typeName(*this) << "::solve(): this class does not provide support for transposed multipication. Consider manually transposing the matrix.");
    return;
  }


  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void CUSPSparseOps<Scalar,Ordinal,Node>::multiply(
                                Teuchos::ETransp trans, 
                                RangeScalar alpha,
                                const MultiVector<DomainScalar,Node> &X, 
                                MultiVector<RangeScalar,Node> &Y) const {
    // beta is not provided and the output multivector enjoys overwrite semantics
    TEST_FOR_EXCEPTION(trans != Teuchos::NO_TRANS, std::logic_error, 
      Teuchos::typeName(*this) << "::multiply(): this class does not provide support for transposed multipication. Consider manually transposing the matrix.");
    TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::multiply(): operation not fully initialized.");
    TEST_FOR_EXCEPT(X.getNumCols() != Y.getNumCols());
    ReadyBufferHelper<Node> rbh(node_);
    if (isEmpty_ == true) {
      // Y <= 0 * X
      //   <= 0
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Init(Y,Teuchos::ScalarTraits<RangeScalar>::zero());
    }
    else {
    }
    return;
  }

  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void CUSPSparseOps<Scalar,Ordinal,Node>::multiply(
                                Teuchos::ETransp trans, 
                                RangeScalar alpha, const MultiVector<DomainScalar,Node> &X, 
                                RangeScalar beta, MultiVector<RangeScalar,Node> &Y) const { 
    TEST_FOR_EXCEPTION(true, std::logic_error,
      Teuchos::typeName(*this) << "::multiply(): CUSP does not support multiple scalar types for sparse matrix-vector multiplication.");
  }

  template <class Scalar, class Ordinal, class Node>
  template <>
  void CUSPSparseOps<Scalar,Ordinal,Node>::multiply(
                                Teuchos::ETransp trans, 
                                Scalar alpha, const MultiVector<Scalar,Node> &X, 
                                Scalar beta, MultiVector<Scalar,Node> &Y) const {
    // beta is provided and the output multivector enjoys accumulate semantics
    TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::multiply(): operation not fully initialized.");
    TEST_FOR_EXCEPT(X.getNumCols() != Y.getNumCols());
    TEST_FOR_EXCEPTION(trans != Teuchos::NO_TRANS, std::logic_error, 
      Teuchos::typeName(*this) << "::multiply(): this class does not provide support for transposed multipication. Consider manually transposing the matrix.");
    const size_t numRHS = X.getNumCols(),
                 Xstride = X.getStride(),
                 Ystride = Y.getStride();
    ReadyBufferHelper<Node> rbh(node_);
    rbh.begin();
    const Scalar * X = rbh.template addConstBuffer<Scalar>(X.getValues());
    Scalar       * Y = rbh.template addNonConstBuffer<Scalar>(Y.getValuesNonConst());
    rbh.end();
    for (int v=0; v != numRHS; ++v) {
      cusp::detail::device<int,Scalar>(*devcHyb_, X, Y);
      X += Xstride;  
      Y += Ystride;  
    }
    return;
  }

  /** \example CrsMatrix_CUSP.cpp 
    * This is an example that unit tests and demonstrates the CUSPSparseOps interface.
    */

  /** \brief A partial specialization of CrsGraph for use with CUSPSparseOps.
      \ingroup kokkos_crs_ops
      
      This implementation is supported by CrsGraphHostCompute, even though the
      computation will occur on the GPU. The reason is that CUSP will be responsible
      for the organization of data and the implementation of the kernel; parent CrsGraphHostCompute 
      will handle of the interfaces needed by Kokkos.
    */
  template <class Ordinal, 
            class Node>
  class CrsGraph<Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> > : public CrsGraphHostCompute<Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> > {
  public:
    CrsGraph(size_t numRows, const Teuchos::RCP<Node> &node) : CrsGraphHostCompute<Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> >(numRows,node) {}
  private:
    CrsGraph(const CrsGraph<Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> > &graph); // not implemented
  };

  /** \brief A partial specialization of CrsMatrix for use with the CUSPSparseOps interface.
      \ingroup kokkos_crs_ops

      This implementation is supported by CrsMatrixHostCompute, even though the
      computation will occur on the GPU. The reason is that CUSP will be responsible
      for the organization of data and the implementation of the kernel; parent CrsMatrixHostCompute 
      will handle of the interfaces needed by Kokkos.
    */
  template <class Scalar,
            class Ordinal, 
            class Node>
  class CrsMatrix<Scalar,Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> > : public CrsMatrixHostCompute<Scalar,Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> > {
  public:
    CrsMatrix() : CrsMatrixHostCompute<Scalar,Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> >() {}
    CrsMatrix(CrsGraph<Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> > &graph) : CrsMatrixHostCompute<Scalar,Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> >(graph) {}
    CrsMatrix(const CrsGraph<Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> > &graph) : CrsMatrixHostCompute<Scalar,Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> >(graph) {}
    void setStaticGraph(const CrsGraph<Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> > &graph) {
      const CrsGraphHostCompute<Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> > &hgraph = dynamic_cast<const CrsGraphHostCompute<Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> > &>(graph);
      CrsMatrixHostCompute<Scalar,Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> >::setStaticGraph(hgraph);
    }
    void setOwnedGraph(CrsGraph<Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> > &graph) {
      CrsGraphHostCompute<Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> > &hgraph = dynamic_cast<CrsGraphHostCompute<Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> > &>(graph);
      CrsMatrixHostCompute<Scalar,Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> >::setOwnedGraph(hgraph);
    }
  private:
    CrsMatrix(const CrsMatrix<Scalar,Ordinal,Node,CUSPSparseOps<void,Ordinal,Node> > &mat); // not implemented
  };


} // namespace Kokkos

#endif /* KOKKOS_DEFAULTSPARSEOPS_HPP */

