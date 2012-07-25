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

#include <Teuchos_CompileTimeAssert.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
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

    SeqCrsGraph (Ordinal numRows, Ordinal numCols,
                 const Teuchos::RCP<Node>& node,
                 const Teuchos::RCP<Teuchos::ParameterList>& params);

    Teuchos::RCP<Node> getNode() const {
      return node_;
    }

    Ordinal getNumRows () const {
      return numRows_;
    }

    Ordinal getNumCols () const {
      return numCols_;
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

    /// \brief Whether the graph is empty.
    ///
    /// "Empty" means either that the graph has no rows (the number of
    /// rows is zero), or that the graph has no stored entries.
    bool isEmpty() const {
      return isEmpty_;
    }

    /// \brief Whether the graph has empty rows.
    ///
    /// An empty graph (see isEmpty()) trivially has empty rows.
    /// Otherwise, the graph has empty rows if one or more rows
    /// contains no stored entries.
    bool hasEmptyRows() const {
      return hasEmptyRows_;
    }

    void setStructure (const Teuchos::ArrayRCP<const size_t>& ptr,
                       const Teuchos::ArrayRCP<const Ordinal>& ind);
    void setMatDesc (Teuchos::EUplo uplo, Teuchos::EDiag diag);
    void getMatDesc (Teuchos::EUplo &uplo, Teuchos::EDiag &diag) const;

  private:
    Teuchos::RCP<Node> node_;
    Ordinal numRows_, numCols_;
    //Teuchos::RCP<ParameterList> params_;
    bool isInitialized_;
    bool isEmpty_;
    bool hasEmptyRows_;
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
  SeqCrsGraph (Ordinal numRows, Ordinal numCols,
               const Teuchos::RCP<Node> &node,
               const Teuchos::RCP<Teuchos::ParameterList> &params) :
    node_ (node),
    numRows_ (numRows),
    numCols_ (numCols),
    isInitialized_ (false),
    isEmpty_ (numRows == 0), // provisional; a matrix with numRows > 0
                             // may still have zero entries.
    hasEmptyRows_ (true), // provisional
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
    const Ordinal numRows = this->getNumRows();

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
      (size_t) ptr.size() != (size_t) numRows+1,
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
      (size_t) ind.size() != (size_t) ptr[numRows],
      std::runtime_error,
      ": Graph input data are not coherent:\n"
      "-- ind.size() = " << ind.size() << " != ptr[numRows="
      << numRows << "] = " << ptr[numRows] << ".");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isInitialized_,
      std::runtime_error,
      ": Graph has already been initialized."
    )

    const Ordinal numEntries = ptr[numRows];
    if (numRows == 0 || numEntries == 0) {
      isEmpty_ = true;
      hasEmptyRows_ = true; // trivially
    }
    else {
      isEmpty_ = false;
      // Check whether the graph has any empty rows.
      bool emptyRows = false;
      for (Ordinal i = 0; i < numRows; ++i) {
        if (ptr[i] == ptr[i+1]) {
          emptyRows = true;
          break;
        }
      }
      hasEmptyRows_ = emptyRows;
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
  class SeqSparseOps : public Teuchos::Describable {
  public:
    //! \name Typedefs and structs
    //@{

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
    //! \name Constructors and destructor
    //@{

    /// \brief Constructor, with default parameters.
    ///
    /// We syntactically forbid setting parameters after construction,
    /// since setting parameters after calling setGraphAndMatrix()
    /// would require reorganizing the already optimized sparse matrix
    /// storage.  If you want to set nondefault values of parameters,
    /// you must use the constructor that takes a ParameterList.
    ///
    /// \param node [in/out] Kokkos Node instance.
    SeqSparseOps (const Teuchos::RCP<Node>& node);

    /// \brief Constructor, with custom parameters.
    ///
    /// Both this constructor and finalizeGraphAndMatrix() accept a
    /// ParameterList.  However, those sets of parameters are
    /// different.  The constructor's parameters concern the
    /// algorithm, and the parameters for finalizeGraphAndMatrix()
    /// concern the data structure.  It's possible to use different
    /// algorithms with the same data structure.
    ///
    /// \param node [in/out] Kokkos Node instance.
    ///
    /// \param params [in/out] Parameters for the solve.  We fill in
    ///   the given ParameterList with its default values, but we
    ///   don't keep it around.  (This saves a bit of memory.)
    SeqSparseOps (const Teuchos::RCP<Node>& node,
                  Teuchos::ParameterList& plist);

    //! Destructor
    ~SeqSparseOps();

    /// \brief Get a default ParameterList.
    ///
    /// The returned ParameterList has all accepted parameters, their
    /// default values, documentation, and validators (if applicable).
    ///
    /// This is a class (static) method so that you can get the
    /// default ParameterList (with built-in documentation) before
    /// constructing a SeqSparseOps instance.
    static Teuchos::RCP<const Teuchos::ParameterList>
    getDefaultParameters ()
    {
      using Teuchos::ParameterList;
      using Teuchos::parameterList;
      using Teuchos::RCP;
      using Teuchos::rcp_const_cast;

      RCP<ParameterList> plist = parameterList ("SeqSparseOps");
      setDefaultParameters (*plist);
      return rcp_const_cast<const ParameterList> (plist);
    }

    //@}
    //! \name Implementation of Teuchos::Describable
    //@{

    //! One-line description of this instance.
    std::string description () const {
      using Teuchos::TypeNameTraits;
      std::ostringstream os;
      os <<  "Kokkos::SeqSparseOps<"
         << "Scalar=" << TypeNameTraits<Scalar>::name()
         << ", Ordinal=" << TypeNameTraits<Ordinal>::name()
         << ", Node=" << TypeNameTraits<Node>::name()
         << ">";
      return os.str();
    }

    //! Write a possibly more verbose description of this instance to out.
    void
    describe (Teuchos::FancyOStream& out,
              const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const
    {
      using Teuchos::EVerbosityLevel;
      using Teuchos::includesVerbLevel;
      using Teuchos::OSTab;
      using Teuchos::rcpFromRef;
      using Teuchos::VERB_DEFAULT;
      using Teuchos::VERB_NONE;
      using Teuchos::VERB_LOW;
      using Teuchos::VERB_MEDIUM;
      using Teuchos::VERB_HIGH;
      using Teuchos::VERB_EXTREME;
      using std::endl;

      // Interpret the default verbosity level as VERB_LOW.
      const EVerbosityLevel vl =
        (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;

      if (vl == VERB_NONE) {
        return;
      }
      else if (includesVerbLevel (VERB_LOW, vl)) { // VERB_LOW >= vl
        out << this->description();

        if (includesVerbLevel (vl, VERB_MEDIUM)) { // vl >= VERB_MEDIUM
          out << ":" << endl;
          OSTab tab1 (rcpFromRef (out));

          out << "matVecVariant_ = " << matVecVariant_ << endl
              << "unroll_ = " << unroll_ << endl
              << "isInitialized_ = " << isInitialized_ << endl;
          if (isInitialized_) {
            out << "numRows_ = " << numRows_ << endl
                << "isEmpty_ = " << isEmpty_ << endl
                << "hasEmptyRows_ = " << hasEmptyRows_ << endl;
            if (ptr_.size() > 0) {
              out << "numEntries = " << ptr_[ptr_.size()-1] << endl;
            }
            else {
              out << "numEntries = 0" << endl;
            }
          }
        }
      }
    }

    //@}
    //! \name Accessor routines
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

    /// \brief Finalize a graph and a matrix.
    ///
    /// Both the constructor and this method accept a ParameterList.
    /// However, those sets of parameters are different.  The
    /// constructor's parameters concern the algorithm, and the
    /// parameters for this method concern the data structure.  It's
    /// possible to use different algorithms with the same data
    /// structure.
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
    /// \brief Which algorithm variant to use for sparse matrix-vector multiply.
    ///
    /// The textbook compressed sparse row (CSR) and compressed sparse
    /// column (CSC) sparse matrix-vector multiply algorithms have two
    /// nested 'for' loops.  The outer for loop is for the rows (for
    /// CSR; columns for CSC), and the inner for loop is for the
    /// entries within a row (for CSR; column for CSC).  We call this
    /// the 'for-for' variant.
    ///
    /// We also make available two variants that use a single 'for'
    /// loop over all the entries in the sparse matrix.  The first,
    /// which we call 'for-while', has an inner whlie loop for
    /// incrementing the current row (for CSR; column for CSC) index.
    /// The second, which we call 'for-if', replaces the while loop in
    /// 'for-while' with a single if statement.  The 'for-if' variant
    /// is only correct if the sparse matrix contains no empty rows
    /// (for CSR; columns for CSC).  If you specify the for-if
    /// variant, we check for empty rows, and use the for-while
    /// variant if there are any empty rows.
    enum EMatVecVariant {
      FOR_FOR,
      FOR_WHILE,
      FOR_IF
    };

    //! Fill the given ParameterList with defaults and validators.
    static void
    setDefaultParameters (Teuchos::ParameterList& plist);

    /// \brief Set the default mat-vec algorithm variant parameter.
    ///
    /// Use this to construct a ParameterList with default values and
    /// validators.
    static void
    setDefaultMatVecVariantParameter (Teuchos::ParameterList& plist);

    //! Set the default multivector unroll parameter.
    static void
    setDefaultUnrollParameter (Teuchos::ParameterList& plist);

    //! Copy constructor (protected and unimplemented)
    SeqSparseOps (const SeqSparseOps& source);

    //! The Kokkos Node instance given to this object's constructor.
    Teuchos::RCP<Node> node_;

    // we do this one of two ways:
    // packed CRS: array of row pointers, array of indices, array of values.

    ArrayRCP<const size_t>  ptr_;
    ArrayRCP<const Ordinal> ind_;
    ArrayRCP<const Scalar>  val_;

    Teuchos::EUplo  tri_uplo_;
    Teuchos::EDiag unit_diag_;
    EMatVecVariant matVecVariant_;
    bool unroll_;

    Ordinal numRows_;
    bool isInitialized_;
    bool isEmpty_;
    bool hasEmptyRows_;
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
    tri_uplo_ (Teuchos::UNDEF_TRI),      // Provisionally
    unit_diag_ (Teuchos::NON_UNIT_DIAG), // Provisionally
    matVecVariant_ (SeqSparseOps<Scalar,Ordinal,Node>::FOR_FOR),
    unroll_ (true),
    numRows_ (0),                        // Provisionally
    isInitialized_ (false),
    isEmpty_ (true),                     // Provisionally
    hasEmptyRows_ (true)                 // Provisionally
  {
    // Make sure that users only specialize SeqSparseOps for Kokkos
    // Node types that are host Nodes (vs. device Nodes, such as GPU
    // Nodes).
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;
  }

  template<class Scalar, class Ordinal, class Node>
  SeqSparseOps<Scalar,Ordinal,Node>::
  SeqSparseOps (const Teuchos::RCP<Node>& node,
                Teuchos::ParameterList& params) :
    node_ (node),
    tri_uplo_ (Teuchos::UNDEF_TRI),      // Provisionally
    unit_diag_ (Teuchos::NON_UNIT_DIAG), // Provisionally
    matVecVariant_ (SeqSparseOps<Scalar,Ordinal,Node>::FOR_FOR),
    unroll_ (true),
    numRows_ (0),                        // Provisionally
    isInitialized_ (false),
    isEmpty_ (true),                     // Provisionally
    hasEmptyRows_ (true)                 // Provisionally
  {
    using Teuchos::getIntegralValue;
    using Teuchos::ParameterList;
    using Teuchos::RCP;

    // Make sure that users only specialize SeqSparseOps for Kokkos
    // Node types that are host Nodes (vs. device Nodes, such as GPU
    // Nodes).
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;

    params.validateParametersAndSetDefaults (*getDefaultParameters ());
    matVecVariant_ = getIntegralValue<EMatVecVariant> (params, "Sparse matrix-vector multiply variant");
    unroll_ = params.get<bool> ("Unroll across multivectors");
  }

  template<class Scalar, class Ordinal, class Node>
  SeqSparseOps<Scalar,Ordinal,Node>::~SeqSparseOps() {
  }

  template<class Scalar, class Ordinal, class Node>
  void
  SeqSparseOps<Scalar,Ordinal,Node>::
  setDefaultParameters (Teuchos::ParameterList& plist)
  {
    setDefaultMatVecVariantParameter (plist);
    setDefaultUnrollParameter (plist);
  }

  template<class Scalar, class Ordinal, class Node>
  void
  SeqSparseOps<Scalar,Ordinal,Node>::
  setDefaultUnrollParameter (Teuchos::ParameterList& plist)
  {
    const bool unroll = true;
    plist.set ("Unroll across multivectors", unroll, "Whether to unroll reads "
               "and writes of multivectors across columns of the input and "
               "ouput multivectors");
  }

  template<class Scalar, class Ordinal, class Node>
  void
  SeqSparseOps<Scalar,Ordinal,Node>::
  setDefaultMatVecVariantParameter (Teuchos::ParameterList& plist)
  {
    using Teuchos::Array;
    using Teuchos::ParameterList;
    using Teuchos::stringToIntegralParameterEntryValidator;

    Array<std::string> strs (3);
    strs[0] = "for-for";
    strs[1] = "for-while";
    strs[2] = "for-if";

    Array<std::string> docs (3);
    docs[0] = "Two nested for loops (textbook algorithm)";
    docs[1] = "Outer for loop, inner while loop";
    docs[2] = "Outer for loop, inner if statement";

    Array<EMatVecVariant> vals (3);
    vals[0] = FOR_FOR;
    vals[1] = FOR_WHILE;
    vals[2] = FOR_IF;

    const std::string paramName ("Sparse matrix-vector multiply variant");
    const std::string paramDoc ("Which algorithm variant to use for sparse matrix-vector multiply");
    plist.set (paramName, strs[0], paramDoc,
               stringToIntegralParameterEntryValidator<EMatVecVariant> (strs(), docs(), vals(), strs[0]));
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

    ArrayRCP<const  size_t> ptr = opgraph->getPointers ();
    ArrayRCP<const Ordinal> ind = opgraph->getIndices ();
    ArrayRCP<const Scalar> val = opmatrix->getValues ();
    const Ordinal numRows = opgraph->getNumRows ();

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
      (size_t) ptr.size() != (size_t) numRows + 1,
      std::invalid_argument,
      ": ptr.size() = " << ptr.size() << " != numRows+1 = " << (numRows + 1) << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ind.size() != val.size(),
      std::invalid_argument,
      ": ind.size() = " << ind.size() << " != val.size() = " << val.size()
      << ", for ptr = opgraph->getPointers() and ind = opgraph->getIndices().");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      (size_t) ptr[numRows] != (size_t) ind.size(),
      std::invalid_argument,
      ": ptr[numRows = " << numRows << "] = " << ptr[numRows]
      << " != ind.size() = " << ind.size() << ", for ptr = "
      "opgraph->getPointers() and ind = opgraph->getIndices().");

    numRows_ = numRows;
    hasEmptyRows_ = opgraph->hasEmptyRows ();

    if (opgraph->isEmpty () || numRows_ == 0) {
      isEmpty_ = true;
      // We have to go through a little trouble because ptr_ is an
      // array of const size_t, but we need to set its first entry
      // here.
      ArrayRCP<size_t> myPtr = arcp<size_t> (1);
      myPtr[0] = 0;
      ptr_ = arcp_const_cast<const size_t> (myPtr);
      myPtr = null;

      // The matrix is empty, so ind_ and val_ have zero entries.
      ind_ = null;
      val_ = null;
    }
    else {
      isEmpty_ = false;
      // We can just use these arrays directly.
      ptr_ = ptr;
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

//     std::string triUplo;
//     if (tri_uplo_ == Teuchos::UNDEF_TRI) {
//       triUplo = "UNDEF_TRI";
//     }
//     else if (tri_uplo_ == Teuchos::LOWER_TRI) {
//       triUplo = "LOWER_TRI";
//     }
//     else if (tri_uplo_ == Teuchos::UPPER_TRI) {
//       triUplo = "UPPER_TRI";
//     }
//     std::string unitDiag;
//     if (unit_diag_ == Teuchos::NON_UNIT_DIAG) {
//       unitDiag = "NON_UNIT_DIAG";
//     }
//     else if (unit_diag_ == Teuchos::UNIT_DIAG) {
//       unitDiag = "UNIT_DIAG";
//     }
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
      static_cast<size_t> (X.getNumRows()) < static_cast<size_t> (numRows_),
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

      const  size_t* const ptr = ptr_.getRawPtr ();
      const Ordinal* const ind = ind_.getRawPtr ();
      const Scalar*  const val = val_.getRawPtr ();
      const Ordinal numRows = X.getNumRows ();
      const Ordinal numCols = Y.getNumRows ();
      const Ordinal numVecs = X.getNumCols ();

      if (trans == Teuchos::NO_TRANS) {
        if (tri_uplo_ == Teuchos::LOWER_TRI) {
          if (unit_diag_ == Teuchos::UNIT_DIAG) {
            using Kokkos::Raw::lowerTriSolveCsrColMajorUnitDiag;
            lowerTriSolveCsrColMajorUnitDiag<OT, MST, DST, RST> (numRows, numCols,
                                                                 numVecs,
                                                                 X_raw, X_stride,
                                                                 ptr, ind, val,
                                                                 Y_raw, Y_stride);
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
            upperTriSolveCsrColMajor<OT, MST, DST, RST> (numRows, numCols,
                                                         numVecs,
                                                         X_raw, X_stride,
                                                         ptr, ind, val,
                                                         Y_raw, Y_stride);
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

    const  size_t* const ptr = ptr_.getRawPtr ();
    const Ordinal* const ind = ind_.getRawPtr ();
    const Scalar*  const val = val_.getRawPtr ();

    // Pointer to the sparse matrix-vector multiply routine to use.
    void (*matVec) (const OT, const OT, const OT,
                    const RST&, RST[], const OT,
                    const RST&, const size_t[], const OT[], const MST[],
                    const DST[], const OT);

    // The following very long switch statement selects one of the
    // hard-coded-numVecs routines for certain values of numVecs.
    // (Hard-coding the number of columns in the multivectors avoids
    // two branches and an integer addition.)  Otherwise, it picks a
    // general routine.  Here is also where we use the parameters
    // given to the constructor to pick the algorithm variant.
    //
    // Note that we're taking numRows and numCols from Y resp. X.
    // Assuming that the dimensions of X and Y are correct, then
    // whether or not we're applying the transpose, the (transposed,
    // if applicable) matrix has dimensions numRows by numCols.
    // That's why we don't switch the order of numRows, numCols in the
    // invocations below.
    switch (numVecs) {
    case 1:
      if (trans == Teuchos::NO_TRANS) {
        if (matVecVariant_ == FOR_FOR) {
          matVec = &Kokkos::Raw::matVecCsrColMajorForfor1Vec;
        }
        else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
          matVec = &Kokkos::Raw::matVecCsrColMajorForif1Vec;
        }
        else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
          matVec = &Kokkos::Raw::matVecCsrColMajorForwhile1Vec;
        }
      }
      else if (trans == Teuchos::TRANS) {
        if (matVecVariant_ == FOR_FOR) {
          matVec = &Kokkos::Raw::matVecCscColMajorForfor1Vec;
        }
        else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
          matVec = &Kokkos::Raw::matVecCscColMajorForif1Vec;
        }
        else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
          matVec = &Kokkos::Raw::matVecCscColMajorForwhile1Vec;
        }
      }
      else { // if (trans == Teuchos::CONJ_TRANS) {
        if (matVecVariant_ == FOR_FOR) {
          matVec = &Kokkos::Raw::matVecCscColMajorForforConj1Vec;
        }
        else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
          matVec = &Kokkos::Raw::matVecCscColMajorForifConj1Vec;
        }
        else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
          matVec = &Kokkos::Raw::matVecCscColMajorForwhileConj1Vec;
        }
      }
      break;
    case 2:
      if (trans == Teuchos::NO_TRANS) {
        if (matVecVariant_ == FOR_FOR) {
          matVec = &Kokkos::Raw::matVecCsrColMajorForfor2Vec;
        }
        else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
          matVec = &Kokkos::Raw::matVecCsrColMajorForif2Vec;
        }
        else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
          matVec = &Kokkos::Raw::matVecCsrColMajorForwhile2Vec;
        }
      }
      else if (trans == Teuchos::TRANS) {
        if (matVecVariant_ == FOR_FOR) {
          matVec = &Kokkos::Raw::matVecCscColMajorForfor2Vec;
        }
        else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
          matVec = &Kokkos::Raw::matVecCscColMajorForif2Vec;
        }
        else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
          matVec = &Kokkos::Raw::matVecCscColMajorForwhile2Vec;
        }
      }
      else { // if (trans == Teuchos::CONJ_TRANS) {
        if (matVecVariant_ == FOR_FOR) {
          matVec = &Kokkos::Raw::matVecCscColMajorForforConj2Vec;
        }
        else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
          matVec = &Kokkos::Raw::matVecCscColMajorForifConj2Vec;
        }
        else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
          matVec = &Kokkos::Raw::matVecCscColMajorForwhileConj2Vec;
        }
      }
      break;
    case 3:
      if (trans == Teuchos::NO_TRANS) {
        if (matVecVariant_ == FOR_FOR) {
          matVec = &Kokkos::Raw::matVecCsrColMajorForfor3Vec;
        }
        else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
          matVec = &Kokkos::Raw::matVecCsrColMajorForif3Vec;
        }
        else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
          matVec = &Kokkos::Raw::matVecCsrColMajorForwhile3Vec;
        }
      }
      else if (trans == Teuchos::TRANS) {
        if (matVecVariant_ == FOR_FOR) {
          matVec = &Kokkos::Raw::matVecCscColMajorForfor3Vec;
        }
        else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
          matVec = &Kokkos::Raw::matVecCscColMajorForif3Vec;
        }
        else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
          matVec = &Kokkos::Raw::matVecCscColMajorForwhile3Vec;
        }
      }
      else { // if (trans == Teuchos::CONJ_TRANS) {
        if (matVecVariant_ == FOR_FOR) {
          matVec = &Kokkos::Raw::matVecCscColMajorForforConj3Vec;
        }
        else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
          matVec = &Kokkos::Raw::matVecCscColMajorForifConj3Vec;
        }
        else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
          matVec = &Kokkos::Raw::matVecCscColMajorForwhileConj3Vec;
        }
      }
      break;
    case 4:
      if (trans == Teuchos::NO_TRANS) {
        if (matVecVariant_ == FOR_FOR) {
          matVec = &Kokkos::Raw::matVecCsrColMajorForfor4Vec;
        }
        else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
          matVec = &Kokkos::Raw::matVecCsrColMajorForif4Vec;
        }
        else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
          matVec = &Kokkos::Raw::matVecCsrColMajorForwhile4Vec;
        }
      }
      else if (trans == Teuchos::TRANS) {
        if (matVecVariant_ == FOR_FOR) {
          matVec = &Kokkos::Raw::matVecCscColMajorForfor4Vec;
        }
        else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
          matVec = &Kokkos::Raw::matVecCscColMajorForif4Vec;
        }
        else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
          matVec = &Kokkos::Raw::matVecCscColMajorForwhile4Vec;
        }
      }
      else { // if (trans == Teuchos::CONJ_TRANS) {
        if (matVecVariant_ == FOR_FOR) {
          matVec = &Kokkos::Raw::matVecCscColMajorForforConj4Vec;
        }
        else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
          matVec = &Kokkos::Raw::matVecCscColMajorForifConj4Vec;
        }
        else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
          matVec = &Kokkos::Raw::matVecCscColMajorForwhileConj4Vec;
        }
      }
      break;
    default: // The "general case"
      if (unroll_) {
        if (trans == Teuchos::NO_TRANS) {
          if (matVecVariant_ == FOR_FOR) {
            matVec = &Kokkos::Raw::matVecCsrColMajorForfor4Unrolled;
          }
          else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
            matVec = &Kokkos::Raw::matVecCsrColMajorForif4Unrolled;
          }
          else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
            matVec = &Kokkos::Raw::matVecCsrColMajorForwhile4Unrolled;
          }
        }
        else if (trans == Teuchos::TRANS) {
          if (matVecVariant_ == FOR_FOR) {
            matVec = &Kokkos::Raw::matVecCscColMajorForfor4Unrolled;
          }
          else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
            matVec = &Kokkos::Raw::matVecCscColMajorForif4Unrolled;
          }
          else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
            matVec = &Kokkos::Raw::matVecCscColMajorForwhile4Unrolled;
          }
        }
        else { // if (trans == Teuchos::CONJ_TRANS) {
          if (matVecVariant_ == FOR_FOR) {
            matVec = &Kokkos::Raw::matVecCscColMajorForforConj4Unrolled;
          }
          else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
            matVec = &Kokkos::Raw::matVecCscColMajorForifConj4Unrolled;
          }
          else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
            matVec = &Kokkos::Raw::matVecCscColMajorForwhileConj4Unrolled;
          }
        }
      }
      else { // Don't unroll across multivector columns
        if (trans == Teuchos::NO_TRANS) {
          if (matVecVariant_ == FOR_FOR) {
            matVec = &Kokkos::Raw::matVecCsrColMajorForfor;
          }
          else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
            matVec = &Kokkos::Raw::matVecCsrColMajorForif;
          }
          else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
            matVec = &Kokkos::Raw::matVecCsrColMajorForwhile;
          }
        }
        else if (trans == Teuchos::TRANS) {
          if (matVecVariant_ == FOR_FOR) {
            matVec = &Kokkos::Raw::matVecCscColMajorForfor;
          }
          else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
            matVec = &Kokkos::Raw::matVecCscColMajorForif;
          }
          else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
            matVec = &Kokkos::Raw::matVecCscColMajorForwhile;
          }
        }
        else { // if (trans == Teuchos::CONJ_TRANS) {
          if (matVecVariant_ == FOR_FOR) {
            matVec = &Kokkos::Raw::matVecCscColMajorForforConj;
          }
          else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
            matVec = &Kokkos::Raw::matVecCscColMajorForifConj;
          }
          else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
            matVec = &Kokkos::Raw::matVecCscColMajorForwhileConj;
          }
        }
      }
    }

    // Now we know what mat-vec routine to call, so call it.
    matVec (numRows, numCols, numVecs, beta, Y_raw, Y_stride,
            alpha, ptr, ind, val, X_raw, X_stride);
  }

  template <class Scalar, class Ordinal, class Node>
  Teuchos::ArrayRCP< size_t>
  SeqSparseOps<Scalar, Ordinal, Node>::
  allocRowPtrs (const Teuchos::ArrayView<const  size_t>& numEntriesPerRow)
  {
    Teuchos::ArrayRCP< size_t> ptr (numEntriesPerRow.size() + 1);
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
  allocStorage (const Teuchos::ArrayView<const  size_t>& rowPtrs)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(rowPtrs.size() == 0, std::invalid_argument,
      "SeqSparseOps::allocStorage: The input rowPtrs array must have length at "
      "least one, but rowPtrs.size() = " << rowPtrs.size() << ".");

    const Ordinal totalNumEntries = rowPtrs[rowPtrs.size() - 1];
    Teuchos::ArrayRCP<T> val (totalNumEntries);
    std::fill (val.getRawPtr (),
               val.getRawPtr () + totalNumEntries,
               Teuchos::ScalarTraits<T>::zero ());
    return val;
  }

} // namespace Kokkos

#endif // #ifndef __Kokkos_MySparseOps_hpp

