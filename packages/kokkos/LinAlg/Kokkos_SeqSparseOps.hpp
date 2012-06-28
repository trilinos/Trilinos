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

#ifndef __Kokkos_MySparseOps_hpp
#define __Kokkos_MySparseOps_hpp

#include <Teuchos_DataAccess.hpp>
#include <Teuchos_CompileTimeAssert.hpp>
#include <iterator>
#include <stdexcept>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CrsMatrixBase.hpp"
#include "Kokkos_CrsGraphBase.hpp"

#include "Kokkos_MultiVector.hpp"
#include "Kokkos_DefaultArithmetic.hpp"

#include "Kokkos_Raw_SparseMatVec_decl.hpp"
#include "Kokkos_Raw_SparseMatVec_def.hpp"
#include "Kokkos_Raw_SparseTriangularSolve_decl.hpp"
#include "Kokkos_Raw_SparseTriangularSolve_def.hpp"


namespace Kokkos {

  /// \class SeqCrsGraph
  /// \brief Local sparse graph in compressed sparse row format;
  ///   suitable for host-based Kokkos Nodes.
  template <class Ordinal,class Node>
  class SeqCrsGraph {
  public:
    typedef Ordinal ordinal_type;
    typedef Node node_type;

    SeqCrsGraph (size_t numRows,
                 const Teuchos::RCP<Node>& node,
                 const Teuchos::RCP<Teuchos::ParameterList>& params);

    Teuchos::RCP<Node> getNode() const {
      return node_;
    }

    size_t getNumRows () const {
      return numRows_;
    }

    Teuchos::ArrayRCP<const size_t> getPointers() const {
      return ptr_;
    }

    Teuchos::ArrayRCP<const Ordinal> getIndices() const {
      return ind_;
    }

    bool isInitialized() const {
      return isInitialized_;
    }

    /// \brief Whether the matrix is empty.
    ///
    /// "Empty" means either that the matrix has zero rows, or zero
    /// stored entries.
    bool isEmpty() const {
      return isEmpty_;
    }

    void setStructure (const Teuchos::ArrayRCP<const size_t>& ptr,
                       const Teuchos::ArrayRCP<const Ordinal>& ind);
    void setMatDesc (Teuchos::EUplo uplo, Teuchos::EDiag diag);
    void getMatDesc (Teuchos::EUplo &uplo, Teuchos::EDiag &diag) const;

  private:
    Teuchos::RCP<Node> node_;
    size_t numRows_;
    //Teuchos::RCP<ParameterList> params_;
    bool isInitialized_;
    bool isEmpty_;
    Teuchos::EUplo tri_uplo_;
    Teuchos::EDiag unit_diag_;

    Teuchos::ArrayRCP<const size_t> ptr_;
    Teuchos::ArrayRCP<const Ordinal> ind_;
  };


  /// \class SeqCrsMatrix
  /// \brief Local sparse matrix in compressed sparse row format;
  ///   suitable for host-based Kokkos Nodes.
  ///
  /// \note Tied to a particular SeqCrsGraph instance that defines the
  ///   structure of the sparse matrix.
  template <class Scalar,
            class Ordinal,
            class Node>
  class SeqCrsMatrix {
  public:
    typedef Scalar scalar_type;
    typedef Ordinal ordinal_type;
    typedef Node node_type;
    typedef SeqCrsGraph<Ordinal,Node> graph_type;

    SeqCrsMatrix (const Teuchos::RCP<const SeqCrsGraph<Ordinal,Node> > &graph,
                  const Teuchos::RCP<Teuchos::ParameterList> &params);

    void setValues (const Teuchos::ArrayRCP<const Scalar>& val);

    Teuchos::ArrayRCP<const Scalar> getValues() const {
      return val_;
    }

    bool isInitialized() const {
      return isInitialized_;
    }

  private:
    Teuchos::RCP<const graph_type> graph_;
    Teuchos::ArrayRCP<const Scalar> val_;
    bool isInitialized_;
  };

  template <class Ordinal, class Node>
  SeqCrsGraph<Ordinal,Node>::
  SeqCrsGraph (size_t numRows,
               const Teuchos::RCP<Node> &node,
               const Teuchos::RCP<Teuchos::ParameterList> &params) :
    node_ (node),
    numRows_ (numRows),
    isInitialized_ (false),
    isEmpty_ (numRows == 0), // provisional; a matrix with numRows > 0
                             // may still have zero entries.
    tri_uplo_ (Teuchos::UNDEF_TRI),
    unit_diag_ (Teuchos::NON_UNIT_DIAG)
  {
    // Make sure that users only specialize for Kokkos Node types that
    // are host Nodes (vs. device Nodes, such as GPU Nodes).
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void) cta;

    // We don't use params currently.
    (void) params;
  }

  template <class Ordinal, class Node>
  void
  SeqCrsGraph<Ordinal,Node>::
  setStructure (const Teuchos::ArrayRCP<const size_t> &ptr,
                const Teuchos::ArrayRCP<const Ordinal> &ind)
  {
    std::string tfecfFuncName("setStructure(ptr,ind)");
    const size_t numRows = this->getNumRows();

    // mfh 19 June 2012: The tests expect std::runtime_error rather
    // than the arguably more appropriate std::invalid_argument, so
    // I'll throw std::runtime_error here.  Ditto for the other checks
    // below.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ptr.is_null (),
      std::runtime_error,
      ": The input array 'ptr' must be nonnull, even for a matrix with zero "
      "rows.");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      (size_t) ptr.size() != numRows+1,
      std::runtime_error,
      ": Graph input data are not coherent:\n"
      "-- ptr.size() = " << ptr.size() << " != numRows+1 = "
      << (numRows+1) << ".");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ptr[0] != 0,
      std::runtime_error,
      ": Graph input data are not coherent:\n"
      "-- ptr[0] = " << ptr[0] << " != 0.");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      (size_t) ind.size() != ptr[numRows],
      std::runtime_error,
      ": Graph input data are not coherent:\n"
      "-- ind.size() = " << ind.size() << " != ptr[numRows="
      << numRows << "] = " << ptr[numRows] << ".");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isInitialized_,
      std::runtime_error,
      ": Graph has already been initialized."
    )

    const size_t numEntries = ptr[numRows];
    if (numRows == 0 || numEntries == 0) {
      isEmpty_ = true;
    }
    else {
      isEmpty_ = false;
    }
    ptr_ = ptr;
    ind_ = ind;
    isInitialized_ = true;
  }

  template <class Ordinal, class Node>
  void
  SeqCrsGraph<Ordinal,Node>::
  setMatDesc (Teuchos::EUplo uplo, Teuchos::EDiag diag)
  {
    tri_uplo_ = uplo;
    unit_diag_ = diag;
  }

  template <class Ordinal, class Node>
  void
  SeqCrsGraph<Ordinal,Node>::
  getMatDesc (Teuchos::EUplo &uplo, Teuchos::EDiag &diag) const
  {
    uplo = tri_uplo_;
    diag = unit_diag_;
  }

  template <class Scalar, class Ordinal, class Node>
  SeqCrsMatrix<Scalar,Ordinal,Node>::
  SeqCrsMatrix (const Teuchos::RCP<const SeqCrsGraph<Ordinal,Node> >& graph,
                const Teuchos::RCP<Teuchos::ParameterList>& params) :
    graph_ (graph),
    isInitialized_ (false)
  {
    // Make sure that users only specialize for Kokkos Node types that
    // are host Nodes (vs. device Nodes, such as GPU Nodes).
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void) cta;

    // We don't use params currently.
    (void) params;
  }

  template <class Scalar, class Ordinal, class Node>
  void
  SeqCrsMatrix<Scalar,Ordinal,Node>::
  setValues (const Teuchos::ArrayRCP<const Scalar> &val)
  {
    std::string tfecfFuncName("setValues(val)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isInitialized_, std::runtime_error, ": The matrix is already initialized."
    );
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      graph_.is_null() && ! val.is_null(),
      std::runtime_error,
      ": The matrix has a null graph, but you're trying to give it a nonnull "
      "array of values."
    );
    val_ = val;
    if (val_.is_null ()) {
      isInitialized_ = false;
    }
    else {
      isInitialized_ = true;
    }
  }

  /// \class SeqSparseOps
  /// \brief Implementation of local sequential sparse matrix-vector multiply
  ///   and solve routines, for host-based Kokkos Node types.
  /// \ingroup kokkos_crs_ops
  ///
  /// \tparam Scalar The type of entries of the sparse matrix.
  /// \tparam Ordinal The type of (local) indices of the sparse matrix.
  /// \tparam Node The Kokkos Node type.
  template <class Scalar, class Ordinal, class Node>
  class SeqSparseOps {
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
    typedef SeqSparseOps<Scalar, Ordinal, Node> sparse_ops_type;

    //! Typedef for local graph class
    template <class O, class N>
    struct graph {
      typedef SeqCrsGraph<O,N> graph_type;
    };

    //! Typedef for local matrix class
    template <class S, class O, class N>
    struct matrix {
      typedef SeqCrsMatrix<S,O,N> matrix_type;
    };

    /// \brief Sparse operations type for a different scalar type.
    ///
    /// The bind_scalar struct defines the type responsible for sparse
    /// operations for a scalar type S2, which may be different from
    /// \c Scalar.
    ///
    /// Use by Tpetra CrsMatrix to bind a potentially "void" scalar
    /// type to the appropriate scalar.
    ///
    /// This always specifies a specialization of SeqSparseOps,
    /// regardless of the scalar type S2.
    ///
    /// \tparam S2 A scalar type possibly different from \c Scalar.
    template <class S2>
    struct bind_scalar {
      typedef SeqSparseOps<S2, Ordinal, Node> other_type;
    };

    //@}
    //! @name Constructors/Destructor
    //@{

    /// \brief Constructor.
    ///
    /// \param node [in/out] Kokkos Node instance.
    SeqSparseOps (const Teuchos::RCP<Node>& node);

    //! Destructor
    ~SeqSparseOps();

    //@}
    //! @name Accessor routines.
    //@{

    //! The Kokkos Node with which this object was instantiated.
    Teuchos::RCP<Node> getNode () const;

    //@}
    //! @name Initialization of graph and matrix
    //@{

    //! \brief Allocate and initialize the storage for the row pointers.
    static Teuchos::ArrayRCP<size_t>
    allocRowPtrs (const ArrayView<const size_t>& numEntriesPerRow);

    //! Allocate and initialize the storage for a sparse graph or matrix.
    template <class T>
    static Teuchos::ArrayRCP<T>
    allocStorage (const Teuchos::ArrayView<const size_t>& rowPtrs);

    //! Finalize a graph
    static void
    finalizeGraph (Teuchos::EUplo uplo,
                   Teuchos::EDiag diag,
                   SeqCrsGraph<Ordinal, Node>& graph,
                   const Teuchos::RCP<Teuchos::ParameterList> &params);

    //! Finalize the matrix of an already-finalized graph.
    static void
    finalizeMatrix (const SeqCrsGraph<Ordinal, Node>& graph,
                    SeqCrsMatrix<Scalar, Ordinal, Node>& matrix,
                    const Teuchos::RCP<Teuchos::ParameterList>& params);

    //! Finalize a graph and a matrix.
    static void
    finalizeGraphAndMatrix (Teuchos::EUplo uplo,
                            Teuchos::EDiag diag,
                            SeqCrsGraph<Ordinal, Node>& graph,
                            SeqCrsMatrix<Scalar, Ordinal, Node>& matrix,
                            const Teuchos::RCP<Teuchos::ParameterList>& params);

    //! Initialize sparse operations with a graph and matrix.
    void
    setGraphAndMatrix (const Teuchos::RCP<const SeqCrsGraph<Ordinal,Node> > &graph,
                       const Teuchos::RCP<const SeqCrsMatrix<Scalar,Ordinal,Node> > &matrix);

    //@}
    //! @name Computational methods
    //@{

    /// \brief Y := alpha * Op(A) * X.
    ///
    /// Apply the local sparse matrix A (or its transpose or conjugate
    /// transpose) to a multivector X, overwriting Y with the result.
    /// Op(A) means A, the transpose of A, or the conjugate transpose
    /// of A, depending on the \c trans argument.
    ///
    /// \tparam DomainScalar The type of entries in the input
    ///   multivector X.  This may differ from the type of entries in
    ///   A or in Y.
    ///
    /// \tparam RangeScalar The type of entries in the output
    ///   multivector Y.  This may differ from the type of entries in
    ///   A or in X.
    ///
    /// \param trans [in] Whether to apply the matrix, its transpose,
    ///   or its conjugate transpose (if applicable).
    ///
    /// \param alpha [in] Scalar constant \f$\alpha\f$ by which to
    ///   multiply the result of the sparse matrix-(multi)vector
    ///   multiply.
    ///
    /// \param X [in] Input multivector.
    ///
    /// \param Y [out] Result multivector. Contents will be overwritten.
    template <class DomainScalar, class RangeScalar>
    void
    multiply (Teuchos::ETransp trans,
              RangeScalar alpha,
              const MultiVector<DomainScalar,Node> &X,
              MultiVector<RangeScalar,Node> &Y) const;

    /// \brief Y := beta * Y + alpha * Op(A) * X.
    ///
    /// Apply the local sparse matrix A (or its transpose or conjugate
    /// transpose) to a multivector X, accumulating the result into Y.
    /// Op(A) means A, the transpose of A, or the conjugate transpose
    /// of A, depending on the \c trans argument.
    ///
    /// \tparam DomainScalar The type of entries in the input
    ///   multivector X.  This may differ from the type of entries in
    ///   A or in Y.
    ///
    /// \tparam RangeScalar The type of entries in the output
    ///   multivector Y.  This may differ from the type of entries in
    ///   A or in X.
    ///
    /// \param trans [in] Whether to apply the matrix, its transpose,
    ///   or its conjugate transpose (if applicable).
    ///
    /// \param alpha [in] Scalar constant \f$\alpha\f$ by which to
    ///   multiply the result of the sparse matrix-(multi)vector
    ///   multiply.
    ///
    /// \param X [in] Input multivector.
    ///
    /// \param beta [in] Scalar constant \f$\beta\f$ by which to
    ///   multiply Y when summing with the result of the sparse
    ///   matrix-(multi)vector multiply.
    ///
    /// \param Y [in/out] Result multivector.
    template <class DomainScalar, class RangeScalar>
    void
    multiply (Teuchos::ETransp trans,
              RangeScalar alpha,
              const MultiVector<DomainScalar,Node> &X,
              RangeScalar beta,
              MultiVector<RangeScalar,Node> &Y) const;

    /// \brief Solve Y = Op(A) X for X, where we assume A is triangular.
    ///
    /// Solve the (upper or lower) triangular system Y = Op(A) X.
    /// Op(A) means A, the transpose of A, or the conjugate transpose
    /// of A, depending on the \c trans argument.
    ///
    /// \tparam DomainScalar The type of entries in the input
    ///   multivector X.  This may differ from the type of entries in
    ///   A or in Y.
    ///
    /// \tparam RangeScalar The type of entries in the output
    ///   multivector Y.  This may differ from the type of entries in
    ///   A or in X.
    ///
    /// \param trans [in] Whether to solve with the matrix, its
    ///   transpose, or its conjugate transpose (if applicable).
    ///
    /// \param Y [in] Input multivector.
    ///
    /// \param X [out] Result multivector.
    template <class DomainScalar, class RangeScalar>
    void
    solve (Teuchos::ETransp trans,
           const MultiVector<DomainScalar,Node> &Y,
           MultiVector<RangeScalar,Node> &X) const;

    //@}

  private:
    //! Copy constructor (protected and unimplemented)
    SeqSparseOps (const SeqSparseOps& source);

    //! The Kokkos Node instance given to this object's constructor.
    Teuchos::RCP<Node> node_;

    // we do this one of two ways:
    // packed CRS: array of row pointers, array of indices, array of values.

    ArrayRCP<const Ordinal> ptr_;
    ArrayRCP<const Ordinal> ind_;
    ArrayRCP<const Scalar>  val_;

    Teuchos::EUplo  tri_uplo_;
    Teuchos::EDiag unit_diag_;

    size_t numRows_;
    bool isInitialized_;
    bool isEmpty_;
  };

  template <class Scalar, class Ordinal, class Node>
  void
  SeqSparseOps<Scalar,Ordinal,Node>::
  finalizeGraph (Teuchos::EUplo uplo,
                 Teuchos::EDiag diag,
                 SeqCrsGraph<Ordinal,Node>& graph,
                 const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! graph.isInitialized(), std::runtime_error,
      "Kokkos::SeqSparseOps::finalizeGraph: Graph has not yet been initialized."
    );
    graph.setMatDesc (uplo, diag);
  }

  template <class Scalar, class Ordinal, class Node>
  void
  SeqSparseOps<Scalar,Ordinal,Node>::
  finalizeMatrix (const SeqCrsGraph<Ordinal,Node> &graph,
                  SeqCrsMatrix<Scalar,Ordinal,Node> &matrix,
                  const Teuchos::RCP<Teuchos::ParameterList> &params)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! matrix.isInitialized(), std::runtime_error,
      "Kokkos::SeqSparseOps::finalizeMatrix(graph,matrix,params): "
      "Matrix has not yet been initialized."
    );
  }

  template <class Scalar, class Ordinal, class Node>
  void
  SeqSparseOps<Scalar,Ordinal,Node>::
  finalizeGraphAndMatrix (Teuchos::EUplo uplo,
                          Teuchos::EDiag diag,
                          SeqCrsGraph<Ordinal,Node>& graph,
                          SeqCrsMatrix<Scalar,Ordinal,Node>& matrix,
                          const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! graph.isInitialized(), std::runtime_error,
      "Kokkos::SeqSparseOps::finalizeGraphAndMatrix(graph,matrix,params): "
      "Graph has not yet been initialized."
    );
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! matrix.isInitialized(), std::runtime_error,
      "Kokkos::SeqSparseOps::finalizeGraphAndMatrix(graph,matrix,params): "
      "Matrix has not yet been initialized."
    );
    graph.setMatDesc (uplo, diag);
  }


  template<class Scalar, class Ordinal, class Node>
  SeqSparseOps<Scalar,Ordinal,Node>::
  SeqSparseOps (const Teuchos::RCP<Node>& node) :
    node_ (node),
    numRows_ (0),
    isInitialized_ (false),
    isEmpty_ (true) // Provisionally...
  {
    // Make sure that users only specialize SeqSparseOps for Kokkos
    // Node types that are host Nodes (vs. device Nodes, such as GPU
    // Nodes).
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;
  }

  template<class Scalar, class Ordinal, class Node>
  SeqSparseOps<Scalar,Ordinal,Node>::~SeqSparseOps() {
  }

  template <class Scalar, class Ordinal, class Node>
  RCP<Node> SeqSparseOps<Scalar,Ordinal,Node>::getNode() const {
    return node_;
  }

  template <class Scalar, class Ordinal, class Node>
  void
  SeqSparseOps<Scalar,Ordinal,Node>::
  setGraphAndMatrix (const Teuchos::RCP<const SeqCrsGraph<Ordinal,Node> > &opgraph,
                     const Teuchos::RCP<const SeqCrsMatrix<Scalar,Ordinal,Node> > &opmatrix)
  {
    using Teuchos::ArrayRCP;
    using Teuchos::arcp;
    using Teuchos::arcp_const_cast;
    using Teuchos::null;
    // using std::cerr;
    // using std::endl;

    std::string tfecfFuncName("setGraphAndMatrix(uplo,diag,graph,matrix)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isInitialized_, std::runtime_error, " operators already initialized."
    );

    ArrayRCP<const size_t> ptr = opgraph->getPointers ();
    ArrayRCP<const Ordinal> ind = opgraph->getIndices ();
    ArrayRCP<const Scalar> val = opmatrix->getValues ();
    const size_t numRows = opgraph->getNumRows ();

    // cerr << "SeqSparseOps::setGraphAndMatrix: on entry to routine:" << endl
    //      << "ptr = ";
    // std::copy (ptr.begin(), ptr.end(), std::ostream_iterator<size_t> (cerr, " "));
    // cerr << endl << "ind = ";
    // std::copy (ind.begin(), ind.end(), std::ostream_iterator<Ordinal> (cerr, " "));
    // cerr << endl << "val = ";
    // std::copy (val.begin(), val.end(), std::ostream_iterator<Scalar> (cerr, " "));
    // cerr << endl << "numRows = " << numRows << endl;

    // Verify the input data before setting internal state.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      (size_t) ptr.size() != numRows + 1,
      std::invalid_argument,
      ": ptr.size() = " << ptr.size() << " != numRows+1 = " << (numRows + 1) << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ind.size() != val.size(),
      std::invalid_argument,
      ": ind.size() = " << ind.size() << " != val.size() = " << val.size()
      << ", for ptr = opgraph->getPointers() and ind = opgraph->getIndices().");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ptr[numRows] != (size_t) ind.size(),
      std::invalid_argument,
      ": ptr[numRows = " << numRows << "] = " << ptr[numRows]
      << " != ind.size() = " << ind.size() << ", for ptr = "
      "opgraph->getPointers() and ind = opgraph->getIndices().");

    numRows_ = numRows;
    if (opgraph->isEmpty () || numRows_ == 0) {
      isEmpty_ = true;
      // We have to go through a little trouble because ptr_ is an
      // array of const Ordinal, but we need to set its first entry
      // here.
      ArrayRCP<Ordinal> myPtr = arcp<Ordinal> (1);
      myPtr[0] = 0;
      ptr_ = arcp_const_cast<const Ordinal> (myPtr);
      myPtr = null;

      // The matrix is empty, so ind_ and val_ have zero entries.
      ind_ = null;
      val_ = null;
    }
    else {
      isEmpty_ = false;
      // ptr is an array of size_t, ptr_ is an array of Ordinal.  We
      // compute with an array of Ordinal because Ordinal is generally
      // a 32-bit integer type, so it uses less bandwidth and storage.
      // We have to go through a little trouble because ptr_ is an
      // array of const Ordinal, but we need to set its entries here.
      //
      // FIXME (mfh 25 Jun 2012) In debug mode, we should check the
      // conversions from size_t to Ordinal for possible overflow.
      ArrayRCP<Ordinal> myPtr = arcp<Ordinal> (ptr.size ());
      std::copy (ptr.getRawPtr (),
                 ptr.getRawPtr () + ptr.size (),
                 myPtr.getRawPtr ());
      ptr_ = arcp_const_cast<const Ordinal> (myPtr);
      myPtr = null;

      // We can just use these arrays directly.
      ind_ = ind;
      val_ = val;
    }
    opgraph->getMatDesc (tri_uplo_, unit_diag_);
    isInitialized_ = true;

    // cerr << "SeqSparseOps::setGraphAndMatrix: on exit:" << endl
    //      << "ptr_ = ";
    // std::copy (ptr_.begin(), ptr_.end(), std::ostream_iterator<Ordinal> (cerr, " "));
    // cerr << endl << "ind_ = ";
    // std::copy (ind_.begin(), ind_.end(), std::ostream_iterator<Ordinal> (cerr, " "));
    // cerr << endl << "val_ = ";
    // std::copy (val_.begin(), val_.end(), std::ostream_iterator<Scalar> (cerr, " "));

    std::string triUplo;
    if (tri_uplo_ == Teuchos::UNDEF_TRI) {
      triUplo = "UNDEF_TRI";
    }
    else if (tri_uplo_ == Teuchos::LOWER_TRI) {
      triUplo = "LOWER_TRI";
    }
    else if (tri_uplo_ == Teuchos::UPPER_TRI) {
      triUplo = "UPPER_TRI";
    }
    std::string unitDiag;
    if (unit_diag_ == Teuchos::NON_UNIT_DIAG) {
      unitDiag = "NON_UNIT_DIAG";
    }
    else if (unit_diag_ == Teuchos::UNIT_DIAG) {
      unitDiag = "UNIT_DIAG";
    }
    // cerr << endl << "numRows_ = " << numRows_ << endl
    //      << "isEmpty_ = " << isEmpty_ << endl
    //      << "tri_uplo_ = " << triUplo << endl
    //      << "unit_diag_ = " << unitDiag << endl;
  }

  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void
  SeqSparseOps<Scalar,Ordinal,Node>::
  solve (Teuchos::ETransp trans,
         const MultiVector<DomainScalar, Node>& Y,
         MultiVector<RangeScalar, Node>& X) const
  {
    std::string tfecfFuncName("solve(trans,Y,X)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! isInitialized_,
      std::runtime_error,
      ": The solve was not fully initialized."
    );
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      X.getNumCols() != Y.getNumCols(),
      std::runtime_error,
      ": Input and output multivectors have different numbers of vectors (columns)."
    );
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      X.getNumRows() < numRows_,
      std::runtime_error,
      ": Output multivector X does not have enough rows.  X.getNumRows() == "
      << X.getNumRows() << ", but the matrix has " << numRows_ << " rows."
    );

    if (numRows_ == 0) {
      return; // Nothing to do
    }
    else if (isEmpty_) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        unit_diag_ != Teuchos::UNIT_DIAG,
        std::runtime_error,
        ": Solve with empty matrix is only valid for an implicit unit diagonal.");
      // solve I * X = Y for X = Y
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Assign (X, Y);
    }
    else {
      typedef ordinal_type OT;
      typedef scalar_type MST; // matrix scalar type
      typedef DomainScalar DST;
      typedef RangeScalar RST;

      RST* const X_raw = X.getValuesNonConst ().getRawPtr ();
      const Ordinal X_stride = (Ordinal) X.getStride ();
      const DST* const Y_raw = Y.getValues ().getRawPtr ();
      const Ordinal Y_stride = (Ordinal) Y.getStride ();

      const Ordinal* const ptr = ptr_.getRawPtr ();
      const Ordinal* const ind = ind_.getRawPtr ();
      const Scalar*  const val = val_.getRawPtr ();
      const Ordinal numRows = X.getNumRows ();
      const Ordinal numCols = Y.getNumRows ();
      const Ordinal numVecs = X.getNumCols ();

      if (trans == Teuchos::NO_TRANS) {
        if (tri_uplo_ == Teuchos::LOWER_TRI) {
          if (unit_diag_ == Teuchos::UNIT_DIAG) {
            using Kokkos::Raw::lowerTriSolveCsrColMajorUnitDiag;

            if (true) {
              lowerTriSolveCsrColMajorUnitDiag<OT, MST, DST, RST> (numRows, numCols,
                                                                   numVecs,
                                                                   X_raw, X_stride,
                                                                   ptr, ind, val,
                                                                   Y_raw, Y_stride);
            }
            else {
              for (Ordinal r = 0; r < numRows; ++r) {
                for (Ordinal j = 0; j < numVecs; ++j) {
                  X_raw[r + j*X_stride] = Y_raw[r + j*Y_stride];
                }
                for (Ordinal k = ptr[r]; k < ptr[r+1]; ++k) {
                  const MST A_rc = val[k];
                  const Ordinal c = ind[k];
                  TEUCHOS_TEST_FOR_EXCEPTION(c >= r, std::invalid_argument,
                                             "c = " << c << " >= r = " << r << ".");
                  for (Ordinal j = 0; j < numVecs; ++j) {
                    X_raw[r + j*X_stride] -= A_rc * X_raw[c + j*Y_stride];
                  }
                }
              }
            }
          }
          else { // non unit diagonal
            using Kokkos::Raw::lowerTriSolveCsrColMajor;
            lowerTriSolveCsrColMajor<OT, MST, DST, RST> (numRows, numCols,
                                                         numVecs,
                                                         X_raw, X_stride,
                                                         ptr, ind, val,
                                                         Y_raw, Y_stride);
          }
        }
        else { // upper triangular
          if (unit_diag_ == Teuchos::UNIT_DIAG) {
            using Kokkos::Raw::upperTriSolveCsrColMajorUnitDiag;
            upperTriSolveCsrColMajorUnitDiag<OT, MST, DST, RST> (numRows, numCols,
                                                                 numVecs,
                                                                 X_raw, X_stride,
                                                                 ptr, ind, val,
                                                                 Y_raw, Y_stride);
          }
          else { // non unit diagonal
            using Kokkos::Raw::upperTriSolveCsrColMajor;
            if (true) {
              upperTriSolveCsrColMajor<OT, MST, DST, RST> (numRows, numCols,
                                                           numVecs,
                                                           X_raw, X_stride,
                                                           ptr, ind, val,
                                                           Y_raw, Y_stride);
            }
            else {
              typedef Teuchos::ScalarTraits<MST> STS;

              for (Ordinal r = numRows-1; r >= 0; --r) {
                for (Ordinal j = 0; j < numVecs; ++j) {
                  X_raw[r + j*X_stride] = Y_raw[r + j*Y_stride];
                }
                // We assume there's an entry for the diagonal element.
                const MST A_rr = val[ptr[r]];
                for (Ordinal k = ptr[r]+1; k < ptr[r+1]; ++k) {
                  const MST A_rc = val[k];
                  const Ordinal c = ind[k];
                  for (Ordinal j = 0; j < numVecs; ++j) {
                    X_raw[r + j*X_stride] -= A_rc * X_raw[c + j*Y_stride];
                  }
                }
                for (Ordinal j = 0; j < numVecs; ++j) {
                  X_raw[r + j*X_stride] = X_raw[r + j*X_stride] / A_rr;
                }
              }
            }
          }
        }
      }
      else if (trans == Teuchos::TRANS) {
        if (tri_uplo_ == Teuchos::LOWER_TRI) {
          if (unit_diag_ == Teuchos::UNIT_DIAG) {
            using Kokkos::Raw::lowerTriSolveCscColMajorUnitDiag;
            // numRows resp. numCols come from the number of rows in Y
            // resp. X, so they still appear in the same order as
            // in the not transposed cases above.
            lowerTriSolveCscColMajorUnitDiag<OT, MST, DST, RST> (numRows, numCols,
                                                                 numVecs,
                                                                 X_raw, X_stride,
                                                                 ptr, ind, val,
                                                                 Y_raw, Y_stride);
          }
          else {
            using Kokkos::Raw::lowerTriSolveCscColMajor;
            lowerTriSolveCscColMajor<OT, MST, DST, RST> (numRows, numCols,
                                                         numVecs,
                                                         X_raw, X_stride,
                                                         ptr, ind, val,
                                                         Y_raw, Y_stride);
          }
        }
        else { // upper triangular
          if (unit_diag_ == Teuchos::UNIT_DIAG) {
            using Kokkos::Raw::upperTriSolveCscColMajorUnitDiag;
            upperTriSolveCscColMajorUnitDiag<OT, MST, DST, RST> (numRows, numCols,
                                                                 numVecs,
                                                                 X_raw, X_stride,
                                                                 ptr, ind, val,
                                                                 Y_raw, Y_stride);
          }
          else {
            using Kokkos::Raw::upperTriSolveCscColMajor;
            upperTriSolveCscColMajor<OT, MST, DST, RST> (numRows, numCols,
                                                         numVecs,
                                                         X_raw, X_stride,
                                                         ptr, ind, val,
                                                         Y_raw, Y_stride);
          }
        }
      }
      else if (trans == Teuchos::CONJ_TRANS) {
        if (tri_uplo_ == Teuchos::LOWER_TRI) {
          if (unit_diag_ == Teuchos::UNIT_DIAG) {
            using Kokkos::Raw::lowerTriSolveCscColMajorUnitDiagConj;
            lowerTriSolveCscColMajorUnitDiagConj<OT, MST, DST, RST> (numRows, numCols,
                                                                     numVecs,
                                                                     X_raw, X_stride,
                                                                     ptr, ind, val,
                                                                     Y_raw, Y_stride);
          }
          else {
            using Kokkos::Raw::lowerTriSolveCscColMajorConj;
            lowerTriSolveCscColMajorConj<OT, MST, DST, RST> (numRows, numCols,
                                                             numVecs,
                                                             X_raw, X_stride,
                                                             ptr, ind, val,
                                                             Y_raw, Y_stride);
          }
        }
        else { // upper triangular
          if (unit_diag_ == Teuchos::UNIT_DIAG) {
            using Kokkos::Raw::upperTriSolveCscColMajorUnitDiagConj;
            upperTriSolveCscColMajorUnitDiagConj<OT, MST, DST, RST> (numRows, numCols,
                                                                     numVecs,
                                                                     X_raw, X_stride,
                                                                     ptr, ind, val,
                                                                     Y_raw, Y_stride);
          }
          else {
            using Kokkos::Raw::upperTriSolveCscColMajorConj;
            upperTriSolveCscColMajorConj<OT, MST, DST, RST> (numRows, numCols,
                                                             numVecs,
                                                             X_raw, X_stride,
                                                             ptr, ind, val,
                                                             Y_raw, Y_stride);
          }
        }
      }
    }
  }

  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void
  SeqSparseOps<Scalar,Ordinal,Node>::
  multiply (Teuchos::ETransp trans,
            RangeScalar alpha,
            const MultiVector<DomainScalar, Node>& X,
            MultiVector<RangeScalar, Node>& Y) const
  {
    typedef DomainScalar DST;
    typedef RangeScalar RST;
    const RST beta = Teuchos::ScalarTraits<RST>::zero ();
    this->template multiply<DST, RST> (trans, alpha, X, beta, Y);
  }

  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void
  SeqSparseOps<Scalar,Ordinal,Node>::
  multiply (Teuchos::ETransp trans,
            RangeScalar alpha,
            const MultiVector<DomainScalar, Node> &X,
            RangeScalar beta,
            MultiVector<RangeScalar, Node> &Y) const
  {
    using std::cerr;
    using std::endl;

    std::string tfecfFuncName("multiply(trans,alpha,X,beta,Y)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! isInitialized_,
      std::runtime_error,
      ": Sparse ops not initialized.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      X.getNumCols() != Y.getNumCols(),
      std::runtime_error,
      ": X and Y do not have the same number of columns.");

    typedef ordinal_type OT;
    typedef scalar_type MST; // matrix scalar type
    typedef DomainScalar DST;
    typedef RangeScalar RST;

    // These dimensions come from the input and output multivectors,
    // so they apply for the transpose case as well.
    const Ordinal numRows = Y.getNumRows ();
    const Ordinal numCols = X.getNumRows ();
    const Ordinal numVecs = X.getNumCols ();
    RST* const Y_raw = Y.getValuesNonConst ().getRawPtr ();
    const Ordinal Y_stride = Y.getStride ();
    const DST* const X_raw = X.getValues ().getRawPtr ();
    const Ordinal X_stride = X.getStride ();

    const Ordinal* const ptr = ptr_.getRawPtr ();
    const Ordinal* const ind = ind_.getRawPtr ();
    const Scalar*  const val = val_.getRawPtr ();

    // FIXME (mfh 25 Jun 2012) Think about whether the first argument
    // should be always numRows, or should be numRows for CSR and
    // numCols for CSC.
    if (trans == Teuchos::NO_TRANS) {
      using Kokkos::Raw::matVecCsrColMajor;

      if (true) { // One 'for' loop with 'while' loop to advance row index
        if (true) {
          matVecCsrColMajor<OT, MST, DST, RST> (numRows, numCols, numVecs,
                                                beta, Y_raw, Y_stride,
                                                alpha, ptr, ind, val,
                                                X_raw, X_stride);
        }
        else {
          for (Ordinal c = 0; c < numVecs; ++c) {
            Y_raw[c*Y_stride] = beta * Y_raw[c*Y_stride];
          }
          Ordinal i = 0;
          const Ordinal nnz = ptr[numRows];
          for (Ordinal k = 0; k < nnz; ++k) {
            const MST A_ij = val[k];
            const Ordinal j = ind[k];
            while (k >= ptr[i+1]) {
              ++i;
              RST* const Y_i = &Y_raw[i];
              for (Ordinal c = 0; c < numVecs; ++c) {
                Y_i[c*Y_stride] = beta * Y_i[c*Y_stride];
              }
            }
            const DST* const X_j = &X_raw[j];
            RST* const Y_i = &Y_raw[i];
            for (Ordinal c = 0; c < numVecs; ++c) {
              Y_i[c*Y_stride] += alpha * A_ij * X_j[c*X_stride];
            }
          }
        }
      }
      else { // Two nested 'for' loops
        for (Ordinal i = 0; i < numRows; ++i) {
          RST* const Y_i = &Y_raw[i];
          for (Ordinal c = 0; c < numVecs; ++c) {
            Y_i[c*Y_stride] = beta * Y_i[c*Y_stride];
          }
          for (Ordinal k = ptr[i]; k < ptr[i+1]; ++k) {
            const MST A_ij = val[k];
            const Ordinal j = ind[k];
            const DST* const X_j = &X_raw[j];
            for (Ordinal c = 0; c < numVecs; ++c) {
              Y_i[c*Y_stride] += alpha * A_ij * X_j[c*X_stride];
            }
          }
        }
      }
    }
    else if (trans == Teuchos::TRANS) {
      using Kokkos::Raw::matVecCscColMajor;
      if (true) {
        matVecCscColMajor<OT, MST, DST, RST> (numRows, numCols, numVecs,
                                              beta, Y_raw, Y_stride,
                                              alpha, ptr, ind, val,
                                              X_raw, X_stride);
      }
      else {
        // Prescale Y.
        for (Ordinal i = 0; i < numRows; ++i) {
          RST* const Y_i = &Y_raw[i];
          for (Ordinal c = 0; c < numVecs; ++c) {
            Y_i[c*Y_stride] = beta * Y_i[c*Y_stride];
          }
        }
        // It's still numRows, even though we're iterating over columns
        // of the matrix.
        for (Ordinal j = 0; j < numRows; ++j) {
          const DST* const X_j = &X_raw[j];
          for (Ordinal k = ptr[j]; k < ptr[j+1]; ++k) {
            const MST A_ij = val[k];
            const Ordinal i = ind[k];
            RST* const Y_i = &Y_raw[i];
            for (Ordinal c = 0; c < numVecs; ++c) {
              Y_i[c*Y_stride] += alpha * A_ij * X_j[c*X_stride];
            }
          }
        }
      }
    }
    else if (trans == Teuchos::CONJ_TRANS) {
      using Kokkos::Raw::matVecCscColMajorConj;
      matVecCscColMajorConj<OT, MST, DST, RST> (numRows, numCols, numVecs,
                                                beta, Y_raw, Y_stride,
                                                alpha, ptr, ind, val,
                                                X_raw, X_stride);
    }
  }

  template <class Scalar, class Ordinal, class Node>
  Teuchos::ArrayRCP<size_t>
  SeqSparseOps<Scalar, Ordinal, Node>::
  allocRowPtrs (const Teuchos::ArrayView<const size_t>& numEntriesPerRow)
  {
    using Teuchos::ArrayRCP;
    using Teuchos::arcp;

    ArrayRCP<size_t> ptr = arcp<size_t> (numEntriesPerRow.size() + 1);
    ptr[0] = 0;
    std::partial_sum (numEntriesPerRow.getRawPtr (),
                      numEntriesPerRow.getRawPtr () + numEntriesPerRow.size (),
                      ptr.getRawPtr () + 1);
    return ptr;
  }

  template <class Scalar, class Ordinal, class Node>
  template <class T>
  Teuchos::ArrayRCP<T>
  SeqSparseOps<Scalar, Ordinal, Node>::
  allocStorage (const Teuchos::ArrayView<const size_t>& rowPtrs)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(rowPtrs.size() == 0, std::invalid_argument,
      "SeqSparseOps::allocStorage: The input rowPtrs array must have length at "
      "least one, but rowPtrs.size() = " << rowPtrs.size() << ".");

    const size_t totalNumEntries = rowPtrs[rowPtrs.size() - 1];
    // TEUCHOS_TEST_FOR_EXCEPTION(totalNumEntries < Ordinal(0), std::invalid_argument,
    //   "SeqSparseOps::allocStorage: The last element of the input rowPtrs array, "
    //   "representing the total number of stored entries in the matrix, is negative: "
    //   "rowPtrs[" << (rowPtrs.size() - 1) << "] = " << totalNumEntries << " < 0.");
    Teuchos::ArrayRCP<T> val = Teuchos::arcp<T> (totalNumEntries);
    std::fill (val.getRawPtr (),
               val.getRawPtr () + totalNumEntries,
               Teuchos::ScalarTraits<T>::zero ());
    return val;
  }

} // namespace Kokkos

#endif // #ifndef __Kokkos_MySparseOps_hpp

