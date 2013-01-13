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

#ifndef __Kokkos_MklSparseOps_hpp
#define __Kokkos_MklSparseOps_hpp

/// \file Kokkos_MklSparseOps.hpp
/// \brief MklSparseOps: Local sparse kernels using Intel MKL.
/// \ingroup kokkos_crs_ops
///
/// The main point of this header is the declaration and definition of
/// MklSparseOps.  This is an implementation of the "LocalMatOps"
/// template parameter of Tpetra::CrsGraph and Tpetra::CrsMatrix, that
/// works by calling the sparse kernels in Intel's Math Kernel Library
/// (MKL).

// The allocators, graph, and matrix implementations used by
// DefaultHostSparseOps are perfectly fine for MKL's sparse kernels as
// well.  MKL uses a different parallelization scheme which breaks
// first touch, but we can't do anything about that.
#include "Kokkos_DefaultSparseOps.hpp"
#include "Kokkos_Mkl_MatrixDescriptor.hpp"
#include "Kokkos_Mkl_RawSparseKernels.hpp"
#include "Teuchos_ConfigDefs.hpp" // HAVE_TEUCHOS_COMPLEX

#ifdef HAVE_KOKKOSCLASSIC_MKL
// Unfortunately, we have to include the header file in order to have
// a definition of MKL_INT.
#  include <mkl.h>
#else
#  error "Using this header file requires having built Trilinos with the Intel MKL third-party library."
#  ifdef MKL_INT
#    undef MKL_INT
#  endif // MKL_INT
#  define MKL_INT int
#endif // HAVE_KOKKOSCLASSIC_MKL


namespace Kokkos {

  // Forward declaration, for MklBindScalarAndOrdinal below.
  template <class Scalar, class Ordinal, class Node, class Allocator>
  struct MklSparseOps;

  namespace { // (anonymous)

    /// \struct MklBindScalarAndOrdinal
    /// \brief Given Scalar and Ordinal types, decide whether to use
    ///   MKL or a fall-back for local sparse ops.
    ///
    /// \note This class is _not_ meant for users.  It is an
    ///   implementation detail of MklSparseOps (which is defined
    ///   below in this header file).
    ///
    /// \tparam Scalar The type of entries in the sparse matrix; same
    ///   as the Scalar template parameter of MklSparseOps.
    /// \tparam Ordinal The type of column indices in the sparse
    ///   matrix; same as the Ordinal template parameter of
    ///   MklSparseOps.
    /// \tparam Node The Kokkos Node type; same as the Node template
    ///   parameter of MklSparseOps.
    /// \tparam Allocator The type that implements allocation of
    ///   sparse matrix arrays; same as the Allocator template
    ///   parameter of MklSparseOps.
    ///
    /// Intel's Math Kernel Library (MKL) only provides local sparse
    /// kernels (matrix-(multi)vector multiply and triangular solve)
    /// for certain Scalar and Ordinal types.  In turn, our MKL
    /// wrapper, MklSparseOps (see below in this header file), only
    /// has valid specializations for those types.  This means that if
    /// the user's combination of Scalar and Ordinal types is not
    /// supported by MKL, we need a fall-back.
    ///
    /// This struct specifies the fall-back, given a Scalar and an
    /// Ordinal type.  If MKL provides kernels for the given Scalar
    /// and Ordinal types, this struct's \c other_type typedef is a
    /// partial specialization of MklSparseOps.  Otherwise, this
    /// struct's \c other_type typedef is a specialization of
    /// DefaultHostSparseOps.  We say that this struct "binds" the
    /// given Scalar and Ordinal types to the appropriate local sparse
    /// ops implementation, hence the name ("bind Scalar and
    /// Ordinal").
    ///
    /// \section Notes for developers
    ///
    /// \subsection Redefining the fall-back kernels
    ///
    /// The place to define the fall-back implementation of sparse
    /// kernels is in the generic definition of the other_type typedef
    /// below.  For now, we use DefaultHostSparseOps as the generic
    /// fall-back.  Also, you may wish to provide other partial
    /// specializations as fall-backs for Scalar or Ordinal types that
    /// MKL does not support, in case a different implementation
    /// performs better for those Scalar or Ordinal types than
    /// DefaultHostSparseOps.
    ///
    /// \subsection Why this class exists
    ///
    /// MklBindScalarAndOrdinal exists because it was easier to define
    /// partial specializations of a traits class outside of
    /// MklSparseOps, than it was to try to define partial
    /// specializations of the bind_scalar and bind_ordinal traits
    /// classes inside of MklSparseOps.  Thus, bind_scalar and
    /// bind_ordinal in MklSparseOps get their typedefs from
    /// MklBindScalarAndOrdinal.
    template <class Scalar, class Ordinal, class Node, class Allocator>
    struct MklBindScalarAndOrdinal {
      /// \typedef other_type
      /// \brief Implementation of sparse kernels for the given Scalar and Ordinal types.
      typedef DefaultHostSparseOps<Scalar, Ordinal, Node, Allocator> other_type;
    };

    // Partial specialization for Scalar=float and Ordinal=MKL_INT.
    template <class Node, class Allocator>
    struct MklBindScalarAndOrdinal<float, MKL_INT, Node, Allocator> {
      typedef MklSparseOps<float, MKL_INT, Node, Allocator> other_type;
    };

    // Partial specialization for Scalar=double and Ordinal=MKL_INT.
    template <class Node, class Allocator>
    struct MklBindScalarAndOrdinal<double, MKL_INT, Node, Allocator> {
      typedef MklSparseOps<double, MKL_INT, Node, Allocator> other_type;
    };

#ifdef HAVE_TEUCHOS_COMPLEX
    // Partial specialization for Scalar=std::complex<float> and Ordinal=MKL_INT.
    template <class Node, class Allocator>
    struct MklBindScalarAndOrdinal<std::complex<float>, MKL_INT, Node, Allocator> {
      typedef MklSparseOps<std::complex<float>, MKL_INT, Node, Allocator> other_type;
    };

    // Partial specialization for Scalar=std::complex<double> and Ordinal=MKL_INT.
    template <class Node, class Allocator>
    struct MklBindScalarAndOrdinal<std::complex<double>, MKL_INT, Node, Allocator> {
      typedef MklSparseOps<std::complex<double>, MKL_INT, Node, Allocator> other_type;
    };
#endif // HAVE_TEUCHOS_COMPLEX

    /// \class ConvertRowPtrs
    /// \brief Kokkos kernel for converting sparse matrix row offsets
    ///   between 0-based and 1-based.
    ///
    /// \note This class is _not_ meant for users.  It is an
    ///   implementation detail of MklSparseOps.
    class ConvertRowPtrs {
    private:
      MKL_INT* const outRowPtr_;
      const MKL_INT* const inRowPtr_;
      const int newIndexBase_;

    public:
      /// Constructor.
      ///
      /// \param outRowPtr [out] The output row offsets array, that
      ///   uses the new index base.
      /// \param inRowPtr [in] The input row offsets array, that uses
      ///   the old index base.
      /// \param newIndexBase [in] The new index base.  It must be
      ///   either 0 or 1.  If 0, the old index base must have been 1;
      ///   if 1, the old index base must have been 0.
      ConvertRowPtrs (MKL_INT* const outRowPtr,
                      const MKL_INT* const inRowPtr,
                      const int newIndexBase) :
        outRowPtr_ (outRowPtr),
        inRowPtr_ (inRowPtr),
        newIndexBase_ (newIndexBase)
      {}

      /// \brief Method for the Kokkos Node's parallel_for to execute
      ///   in parallel.
      ///
      /// \param i [in] Current row index (0-based), over which
      ///   parallel_for parallelizes.  0 <= i <= M, where M is the
      ///   number of rows.  (Yes, the case i == M is included.
      ///   Remember that the row offsets array has one more entry
      ///   than the number of rows in the matrix.)
      void execute (const MKL_INT i) {
        // base = 0: subtract 1 (converting from base = 1)
        // base = 1: add 1 (converting from base = 0)
        outRowPtr_[i] = inRowPtr_[i] + (2*newIndexBase_ - 1);
      }
    };

    /// \class ConvertColInds
    /// \brief Kokkos kernel for converting sparse matrix column
    ///   indices between 0-based and 1-based.
    ///
    /// \note This class is _not_ meant for users.  It is an
    ///   implementation detail of MklSparseOps.
    class ConvertColInds {
    private:
      MKL_INT* const newColInd_; // 1-based, on output
      const MKL_INT* const newRowPtr_; // 1-based
      const MKL_INT* const oldColInd_; // 0-based
      const int newIndexBase_;

    public:
      /// Constructor.
      ///
      /// \param outColInd [out] The output column indices array, that
      ///   uses the new index base.
      /// \param inColInd [in] The input column indices array, that
      ///   uses the old index base.
      /// \param newIndexBase [in] The new index base.  It must be
      ///   either 0 or 1.  If 0, the old index base must have been 1;
      ///   if 1, the old index base must have been 0.
      ConvertColInds (MKL_INT* const newColInd, // new index base, on output
                      const MKL_INT* const newRowPtr, // new index base
                      const MKL_INT* const oldColInd, // old index base
                      const int newIndexBase) :
        newColInd_ (newColInd),
        newRowPtr_ (newRowPtr),
        oldColInd_ (oldColInd),
        newIndexBase_ (newIndexBase)
      {}

      /// \brief Method for the Kokkos Node's parallel_for to execute
      ///   in parallel.
      ///
      /// \param i [in] Current row index (0-based), over which
      ///   parallel_for parallelizes.  0 <= i < M, where M is the
      ///   number of rows.  (Here, the case i == M is _not_
      ///   included.)
      void execute (const MKL_INT r) {
        const int baseShift = 2*newIndexBase_ - 1;
        const MKL_INT start = newRowPtr_[r] - baseShift;
        const MKL_INT end = newRowPtr_[r+1] - baseShift;

        for (MKL_INT k = start; k < end; ++k) {
          newColInd_[k] = oldColInd_[k] + baseShift;
        }
      }
    };

    /// \class CopyRowPtrs
    /// \brief Kokkos kernel for copying a sparse matrix's row offsets
    ///   array from an array of size_t (the input format) to an array
    ///   of MKL_INT (the internal storage format).
    ///
    /// This kernel makes no attempt to check for overflow when
    /// converting from size_t to MKL_INT.  However, if the input row
    /// offsets array is correct, this is easy to do: just check the
    /// last entry (since the entries are nondecreasing).
    ///
    /// \note This class is _not_ meant for users.  It is an
    ///   implementation detail of MklSparseOps.
    class CopyRowPtrs {
    private:
      MKL_INT* const outRowPtr_;
      const size_t* const inRowPtr_;

    public:
      CopyRowPtrs (MKL_INT* const outRowPtr,
                   const size_t* const inRowPtr) :
        outRowPtr_ (outRowPtr),
        inRowPtr_ (inRowPtr)
      {}

      /// \brief Method for the Kokkos Node's parallel_for to execute
      ///   in parallel.
      ///
      /// \param i [in] Current row index (0-based), over which
      ///   parallel_for parallelizes.  0 <= i < M, where M is the
      ///   number of rows.  (Here, the case i == M is _not_
      ///   included.)
      void execute (const MKL_INT i) {
        outRowPtr_[i] = inRowPtr_[i];
      }
    };
  } // namespace (anonymous)

  /// \class MklSparseOps
  /// \brief Default implementation of sparse matrix-vector multiply
  ///   and solve routines, for host-based Kokkos Node types, using MKL.
  /// \ingroup kokkos_crs_ops
  ///
  /// \tparam Scalar The type of entries of the sparse matrix.
  /// \tparam Ordinal The type of (local) column indices of the sparse
  ///   matrix.
  /// \tparam Node The Kokkos Node type.
  /// \tparam Allocator The allocator to use when allocating sparse
  ///   matrix data.  Depending on the particular Allocator, this may
  ///   or may not do first-touch initialization.
  ///
  /// This is an implementation of the "LocalMatOps" template
  /// parameter of Tpetra::CrsGraph and Tpetra::CrsMatrix, that works
  /// by calling the sparse kernels in Intel's Math Kernel Library
  /// (MKL).  Like all such implementations, it provides local sparse
  /// matrix-(multi)vector multiply and local sparse triangular solve.
  /// It is intended only for host-based Kokkos Nodes, just like
  /// DefaultHostSparseOps.
  ///
  /// MKL only defines sparse matrix kernels for a limited set of
  /// Scalar types, and only for one Ordinal type.  The Ordinal type
  /// for which MKL is defined depends on the MKL library with which
  /// you link.  If MKL_ILP64 is #defined, then the only valid Ordinal
  /// type is a 64-bit signed integer type (int64_t); otherwise,
  /// Ordinal is a 32-bit signed integer type (int, or int32_t).
  ///
  /// We provide specializations of MklSparseOps for all Scalar and
  /// Ordinal types.  For those which MKL do not support, we default
  /// to use our own generic kernels.  As a result, you should always
  /// access MklSparseOps type via the following typedef:
  ///
  /// \code
  /// typedef typename MklSparseOps<Scalar, Ordinal, Node, Allocator>::bind_scalar<Scalar>::other_type::typename bind_ordinal<Ordinal>::other_type sparse_ops_type;
  /// \endcode
  ///
  /// The resulting sparse_ops_type typedef will always be a
  /// functioning implementation of local sparse ops, whether or not
  /// the "original" implementation of local sparse ops (in this case,
  /// MklSparseOps) supports your given combination of Scalar and
  /// Ordinal types.  If MKL supports your Scalar and Ordinal types,
  /// then sparse_ops_type (in the typedef above) will be a
  /// specialization of MklSparseOps.  Otherwise, it will be a
  /// specialization of DefaultHostSparseOps.
  ///
  /// The MKL library itself does not use the Kokkos Node API.  If you
  /// link with the right combination of MKL libraries and specify the
  /// right environment variables, MKL will use OpenMP internally to
  /// parallelize its sparse kernels.  However, depending on the
  /// Allocator template parameter, we may use the Kokkos Node API to
  /// do first-touch allocation in parallel of the sparse graph and
  /// matrix data structures.  It would be best in that case to use
  /// Node=OpenMPNode, so that the first-touch allocation most closely
  /// matches the parallelization scheme.
  ///
  /// We cannot and should not control MKL's parallelization scheme.
  /// We cannot, because this depends on the specific MKL libraries
  /// you are using, and on environment variables you specify at run
  /// time.  We should not, because your application which uses
  /// Trilinos may use MKL for its own purposes, so it would be rude
  /// for Trilinos to change MKL's behavior.  This means that you are
  /// responsible for configuring MKL correctly so as to achieve the
  /// desired parallelization scheme.
  template <class Scalar,
            class Ordinal,
            class Node,
            class Allocator = details::FirstTouchCRSAllocator>
  class MklSparseOps : public Teuchos::Describable {
  public:
    //! \name Typedefs and structs
    //@{

    //! The type of the individual entries of the sparse matrix.
    typedef Scalar  scalar_type;
    //! The type of the (local) column indices that describe the structure of the sparse matrix.
    typedef Ordinal ordinal_type;
    //! The Kokkos Node type.
    typedef Node    node_type;
    //! The allocator type.
    typedef Allocator allocator_type;
    //! The type of this object, the sparse operator object.
    typedef MklSparseOps<Scalar, Ordinal, Node, Allocator> sparse_ops_type;

    /// \brief Typedef for local graph class.
    ///
    /// This is the graph type that MklSparseOps uses as its input
    /// format.  MKL takes the same input graph format as
    /// DefaultHostSparseOps, so we can use the same graph data
    /// structure.
    template <class O, class N>
    struct graph {
      typedef DefaultCrsGraph<O,N> graph_type;
    };

    /// \brief Typedef for local matrix class.
    ///
    /// This is the matrix type that MklSparseOps uses as its input
    /// format.  MKL takes the same input matrix format as
    /// DefaultHostSparseOps, so we can use the same matrix data
    /// structure.
    template <class S, class O, class N>
    struct matrix {
      typedef DefaultCrsMatrix<S,O,N> matrix_type;
    };

    /// \brief Local sparse operations type for a different scalar type.
    ///
    /// The bind_scalar struct defines the type responsible for sparse
    /// operations for a scalar type S2, which may be different from
    /// \c Scalar.
    ///
    /// This class' typedef is used by Tpetra::CrsMatrix to bind a
    /// potentially "void" scalar type to the appropriate scalar.  The
    /// other_type typedef tells Tpetra::CrsMatrix which local sparse
    /// ops type to use, as a function of Tpetra's Scalar template
    /// parameter.
    ///
    /// For MklSparseOps, the other_type typedef is a specialization
    /// of MklSparseOps when S2 and Ordinal are scalar resp. index
    /// types that MKL supports.  Otherwise, it is a specialization of
    /// DefaultHostSparseOps.  Tpetra developers should always get
    /// their local sparse ops type from the other_type typedef.
    ///
    /// \tparam S2 A scalar type possibly different from \c Scalar.
    template <class S2>
    struct bind_scalar {
      typedef typename MklBindScalarAndOrdinal<S2, Ordinal, Node, Allocator>::other_type other_type;
    };

    /// \brief Local sparse operations type for a different ordinal type.
    ///
    /// The bind_ordinal struct defines the type responsible for
    /// sparse operations for an ordinal type O2, which may be
    /// different from \c Ordinal.
    ///
    /// This is used by Tpetra::CrsMatrix to bind the local sparse ops
    /// type to use its own (Local)Ordinal type.  The other_type
    /// typedef is a specialization of MklSparseOps when Scalar and O2
    /// are scalar resp. ordinal types that MKL supports.  Otherwise,
    /// it is a specialization of DefaultHostSparseOps.  Tpetra
    /// developers should always get their local sparse ops type from
    /// the other_type typedef.
    ///
    /// \tparam O2 An ordinal type possibly different from \c Ordinal.
    template <class O2>
    struct bind_ordinal {
      typedef typename MklBindScalarAndOrdinal<Scalar, O2, Node, Allocator>::other_type other_type;
    };

    //@}
    //! \name Constructors/Destructor
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
    MklSparseOps(const RCP<Node> &node);

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
    MklSparseOps(const RCP<Node> &node, Teuchos::ParameterList& params);

    //! Destructor
    ~MklSparseOps();

    //@}
    //! \name Implementation of Teuchos::Describable
    //@{

    //! One-line description of this instance.
    std::string description () const;

    //! Write a possibly more verbose description of this instance to out.
    void
    describe (Teuchos::FancyOStream& out,
              const Teuchos::EVerbosityLevel verbLevel =
              Teuchos::Describable::verbLevel_default) const;

    /// \brief Convert to dense matrix and return.
    ///
    /// \warning This method is for debugging only.  It uses a lot of
    ///   memory.  Users should never call this method.  Do not rely
    ///   on this method continuing to exist in future releases.
    ///
    /// \warning SerialDenseMatrix currently requires Ordinal=int
    ///   indices for BLAS and LAPACK compatibility.  We make no
    ///   attempt to check whether Ordinal -> int conversions
    ///   overflow.
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int, scalar_type> >
    asDenseMatrix () const {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
    }

    //@}
    //! \name Accessor routines.
    //@{

    //! The Kokkos Node with which this object was instantiated.
    RCP<Node> getNode() const {
      return node_;
    }

    //@}
    //! @name Initialization of graph and matrix
    //@{

    /// \brief Allocate and initialize the storage for the row pointers.
    ///
    /// \param node [in/out] Kokkos Node instance.
    /// \param numEntriesPerRow [in] numEntriesPerRow[i] is the number
    ///   of entries in row i, for all rows of the local sparse matrix.
    static ArrayRCP<size_t>
    allocRowPtrs (const Teuchos::RCP<Node>& node,
                  const Teuchos::ArrayView<const size_t>& numEntriesPerRow)
    {
      return Allocator::allocRowPtrs (node, numEntriesPerRow);
    }

    /// \brief Allocate storage for column indices or matrix values.
    ///
    /// \param node [in/out] Kokkos Node instance.
    /// \param rowPtrs [in] The array of row offsets; the 'ptr' array
    ///   in the compressed sparse row storage format.  rowPtrs.size()
    ///   is one plus the number of rows in the local sparse matrix.
    template <class T>
    static ArrayRCP<T>
    allocStorage (const Teuchos::RCP<Node>& node,
                  const Teuchos::ArrayView<const size_t>& rowPtrs)
    {
      return Allocator::template allocStorage<T,size_t>(node,rowPtrs);
    }

    /// \brief Finalize a graph.
    ///
    /// \param uplo [in] Whether the structure of the graph is lower
    ///   triangular (Teuchos::LOWER_TRI), upper triangular
    ///   (Teuchos::UPPER_TRI), or neither (Teuchos::UNDEF_TRI).
    /// \param diag [in] Whether the graph has an implicitly stored
    ///   diagonal (Teuchos::UNIT_DIAG) or does not
    ///   (Teuchos::NON_UNIT_DIAG).  This currently only affects
    ///   sparse triangular solve.
    /// \param graph [in/out] The graph to finalize.
    /// \param params [in/out] Parameters for finalization.
    static void
    finalizeGraph (Teuchos::EUplo uplo,
                   Teuchos::EDiag diag,
                   DefaultCrsGraph<Ordinal,Node> &graph,
                   const Teuchos::RCP<Teuchos::ParameterList> &params);

    /// \brief Finalize the matrix of an already-finalized graph.
    ///
    /// \param graph [in] The graph, which must have already been
    ///   finalized using finalizeGraph().
    /// \param matrix [in/out] The matrix to finalize.  It must have
    ///   been created with the above graph as its graph.
    /// \params [in/out] Parameters for finalization.
    static void
    finalizeMatrix (const DefaultCrsGraph<Ordinal,Node> &graph,
                    DefaultCrsMatrix<Scalar,Ordinal,Node> &matrix,
                    const Teuchos::RCP<Teuchos::ParameterList> &params);

    /// \brief Finalize a graph and a matrix.
    ///
    /// \param uplo [in] Whether the structure of the graph and matrix
    ///   is lower triangular (Teuchos::LOWER_TRI), upper triangular
    ///   (Teuchos::UPPER_TRI), or neither (Teuchos::UNDEF_TRI).
    /// \param diag [in] Whether the matrix has an implicitly stored
    ///   unit diagonal (Teuchos::UNIT_DIAG) or does not
    ///   (Teuchos::NON_UNIT_DIAG).  This currently only affects
    ///   sparse triangular solve.
    /// \param graph [in/out] The graph to finalize.
    /// \param matrix [in/out] The matrix to finalize.  It must have
    ///   been created with the above graph as its graph.
    /// \param params [in/out] Parameters for finalization.
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
                            DefaultCrsGraph<Ordinal,Node> &graph,
                            DefaultCrsMatrix<Scalar,Ordinal,Node> &matrix,
                            const Teuchos::RCP<Teuchos::ParameterList> &params);

    /// \brief Initialize sparse operations with a graph and matrix.
    ///
    /// This is the point at which this object initializes its
    /// internal representation of the local sparse matrix.  In some
    /// cases, this merely involves asking the graph and matrix for
    /// pointers to their data.  In other cases, this involves copying
    /// the data into a completely different storage format.  Whatever
    /// happens is an implementation detail of this object.
    void
    setGraphAndMatrix (const Teuchos::RCP<const DefaultCrsGraph<Ordinal,Node> > &graph,
                       const Teuchos::RCP<const DefaultCrsMatrix<Scalar,Ordinal,Node> > &matrix);

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

    template <class DomainScalar, class RangeScalar>
    void 
    gaussSeidel (const MultiVector<DomainScalar,Node> &B,
		 MultiVector< RangeScalar,Node> &X,
		 const MultiVector<Scalar,Node> &D,
		 const RangeScalar& dampingFactor,
		 const ESweepDirection direction) const
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
				 "MklSparseOps: gaussSeidel not implemented");
    }
    //@}

  private:
    //! Copy constructor (protected and unimplemented)
    MklSparseOps (const MklSparseOps& source);

    /// Convert between 0-based indices and 1-based indices.
    ///
    /// \param newIndexBase [in] 0 or 1.  If the new index base is the
    ///   same as the old one (see getIndexBase()), nothing happens.
    void convertIndexBase (const int newIndexBase);

    //! Get the current index base (0 or 1).
    int getIndexBase () const {
      return matrixDescriptor_.getIndexBase ();
    }

    template <class DomainScalar, class RangeScalar, class OffsetType>
    void solvePrivate(Teuchos::ETransp trans,
                      const RangeScalar& alpha,
                      const MultiVector<DomainScalar,Node> &Y,
                            MultiVector< RangeScalar,Node> &X) const;

    template <class DomainScalar, class RangeScalar, class OffsetType>
    void multiplyPrivate(Teuchos::ETransp trans,
                         RangeScalar alpha,
                         const MultiVector<DomainScalar,Node> &X,
                               MultiVector<RangeScalar,Node> &Y) const;

    template <class DomainScalar, class RangeScalar, class OffsetType>
    void multiplyPrivate(Teuchos::ETransp trans,
                         RangeScalar alpha,
                         const MultiVector<DomainScalar,Node> &X,
                         RangeScalar beta,
                         MultiVector<RangeScalar,Node> &Y) const;

    //! The Kokkos Node instance given to this object's constructor.
    RCP<Node> node_;

    /// \name Packed compressed sparse row (CSR) storage
    ///
    /// The indices and values of row i of the sparse matrix live in
    /// ind_[k] resp. val_[k] for k = ptr_[i] .. ptr_[i+1]-1
    /// (inclusive).
    ///
    /// We use the conventional three-array storage format: an array
    /// of row offsets ("pointers," ptr_), an array of column indices
    /// (ind_), and an array of values (val_).
    ///
    /// MKL requires that the row offsets and column indices have the
    /// same type, which must be MKL_INT.  We check this at run time
    /// (in setGraphAndMatrix()) and throw an exception if the max row
    /// offset does not fit in MKL_INT.
    //@{
    Teuchos::ArrayRCP<const MKL_INT> ptr_; //!< Array of row offsets
    Teuchos::ArrayRCP<const MKL_INT> ind_; //!< Array of column indices
    Teuchos::ArrayRCP<const Scalar>  val_; //!< Array of sparse matrix values
    //@}

    //! Encapsulation of MKL's "matrix descriptor" character array.
    Mkl::MatrixDescriptor matrixDescriptor_;
    MKL_INT numRows_; //!< Number of rows in the matrix.
    MKL_INT numCols_; //!< Number of columns in the matrix.
    //! Whether the matrix was initialized (via setGraphAndMatrix()).
    bool isInitialized_;
    bool isEmpty_; //!< Whether the matrix is empty.

    Teuchos::EUplo tri_uplo_;
    Teuchos::EDiag unit_diag_;
  };

  // // Specializations of MklSparseOps for the Scalar types that it
  // // supports.
  // template <class Scalar, class Ordinal, class Node, class Allocator>
  // template <>
  // struct MklSparseOps<Scalar, Ordinal, Node, Allocator>::bind_scalar<float> {
  //   typedef MklSparseOps<float, Ordinal, Node, Allocator> other_type;
  // };

  // template <class Scalar, class Ordinal, class Node, class Allocator>
  // template <>
  // struct MklSparseOps<Scalar, Ordinal, Node, Allocator>::bind_scalar<double> {
  //   typedef MklSparseOps<double, Ordinal, Node, Allocator> other_type;
  // };

  // template <class Scalar, class Ordinal, class Node, class Allocator>
  // template <>
  // struct MklSparseOps<Scalar, Ordinal, Node, Allocator>::bind_scalar<std::complex<float> > {
  //   typedef MklSparseOps<std::complex<float>, Ordinal, Node, Allocator> other_type;
  // };

  // template <class Scalar, class Ordinal, class Node, class Allocator>
  // template <>
  // struct MklSparseOps<Scalar, Ordinal, Node, Allocator>::bind_scalar<std::complex<double> > {
  //   typedef MklSparseOps<std::complex<double>, Ordinal, Node, Allocator> other_type;
  // };

  // template <class Scalar, class Ordinal, class Node, class Allocator>
  // template <>
  // struct MklSparseOps<Scalar, Ordinal, Node, Allocator>::bind_ordinal<MKL_INT> {
  //   typedef MklSparseOps<Scalar, MKL_INT, Node, Allocator> other_type;
  // };


  template <class Scalar, class Ordinal, class Node, class Allocator>
  void
  MklSparseOps<Scalar,Ordinal,Node,Allocator>::
  finalizeGraph (Teuchos::EUplo uplo,
                 Teuchos::EDiag diag,
                 DefaultCrsGraph<Ordinal,Node> &graph,
                 const RCP<ParameterList> &params)
  {
    using Teuchos::as;

    TEUCHOS_TEST_FOR_EXCEPTION(
      graph.isInitialized() == false,
      std::runtime_error,
      "Kokkos::MklSparseOps::finalizeGraph: Graph has not yet been initialized.");
    RCP<Node> node = graph.getNode ();

    // Determine how many entries are in the graph, so that we can
    // decide whether the row offsets will fit in MKL_INT.
    ArrayRCP<const size_t> bigPtrs = graph.getPointers();
    TEUCHOS_TEST_FOR_EXCEPTION(
      bigPtrs.size() == 0, std::logic_error, "Kokkos::MklSparseOps::finalize"
      "Graph: DefaultCrsGraph's getPointers() returned a row offsets array of "
      "length 0.  This is not allowed even for a graph of zero rows.  Please "
      "report this bug to the Kokkos developers.");

    const size_t numRows = bigPtrs.size() - 1;
    const size_t numnz = bigPtrs[numRows];
    if (numnz < as<size_t> (Teuchos::OrdinalTraits<MKL_INT>::max())) {
      // We may safely convert the row offsets array to MKL_INT.
      // The above test assumes that any Ordinal fits in a size_t.
      ArrayRCP<MKL_INT> smallPtrs;
      {
        // mfh 03 Aug 2012: ArrayRCP isn't thread safe (yet), so use
        // raw pointers in the parallel for.
        const size_t* const ptrIn = bigPtrs.getRawPtr ();
        MKL_INT* const ptrOut = new MKL_INT [bigPtrs.size ()];

        // Copy row offsets in parallel.
        node->parallel_for (0, numRows+1, CopyRowPtrs (ptrOut, ptrIn));
        smallPtrs = arcp<MKL_INT> (ptrOut, 0, numRows+1, true);
      }
      graph.setSmallPointers (smallPtrs);

      // mfh 23 Aug 2012: Don't convert the column indices to 1-based
      // here.  (We can't do this with setStructure(), because that
      // throws an exception if the graph has already been
      // initialized.)  We'll convert to 1-based indices on demand,
      // only if necessary.
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument, "Kokkos::MklSparseOps::finalizeGraph: "
        "The given graph has more entries than can fit in MKL's index type "
        "MKL_INT=" << Teuchos::TypeNameTraits<MKL_INT>::name() << ".  MKL "
        "requires that the row offsets array, and therefore the number of "
        "entries in the graph, has the same type as the column indices.");
    }
    graph.setMatDesc (uplo, diag);
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  void MklSparseOps<Scalar,Ordinal,Node,Allocator>::
  finalizeMatrix (const DefaultCrsGraph<Ordinal,Node> &graph,
                  DefaultCrsMatrix<Scalar,Ordinal,Node> &matrix,
                  const RCP<ParameterList> &params)
  {
    std::string FuncName("Kokkos::MklSparseOps::finalizeMatrix(graph,matrix,params)");
    (void) params; // not using these for now
    // nothing much to do here
    TEUCHOS_TEST_FOR_EXCEPTION(
      matrix.isInitialized() == false,
      std::runtime_error, FuncName << ": matrix has not yet been initialized."
    )
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  void MklSparseOps<Scalar,Ordinal,Node,Allocator>::
  finalizeGraphAndMatrix (Teuchos::EUplo uplo,
                          Teuchos::EDiag diag,
                          DefaultCrsGraph<Ordinal,Node> &graph,
                          DefaultCrsMatrix<Scalar,Ordinal,Node> &matrix,
                          const RCP<ParameterList> &params)
  {
    // finalize them individually; no benefit to doing them together
    finalizeGraph (uplo, diag, graph, params);
    finalizeMatrix (graph, matrix, params);
  }

  template<class Scalar, class Ordinal, class Node, class Allocator>
  MklSparseOps<Scalar,Ordinal,Node,Allocator>::MklSparseOps(const RCP<Node> &node)
  : node_(node)
  , matrixDescriptor_ ("XXXXXX") // invalid; will fill below
  , numRows_(0)
  , isInitialized_(false)
  , isEmpty_(false)
  {
    // Make sure that users only specialize MklSparseOps for
    // Kokkos Node types that are host Nodes (vs. device Nodes, such
    // as GPU Nodes).
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;

    // mfh 04 Aug 2012: We don't currently support symmetric (and
    // related) or diagonal storage formats, but we allow them in the
    // matrix descriptor for future expansion.
    //
    // mfh 23 Aug 2012: We start with 0-based indices.  This will
    // change post-initialization on demand, if multiply() or solve()
    // are given multiple right-hand sides.
    const int indexBase = 0;
    matrixDescriptor_ =
      Mkl::MatrixDescriptor (false, false, false, false, false,
                             Teuchos::UNDEF_TRI, Teuchos::NON_UNIT_DIAG, indexBase);
  }

  template<class Scalar, class Ordinal, class Node, class Allocator>
  MklSparseOps<Scalar,Ordinal,Node,Allocator>::
  MklSparseOps (const RCP<Node> &node, Teuchos::ParameterList& params)
  : node_(node)
  , numRows_(0)
  , isInitialized_(false)
  , isEmpty_(false)
  {
    (void) params; // Not using this yet.

    // Make sure that users only specialize MklSparseOps for
    // Kokkos Node types that are host Nodes (vs. device Nodes, such
    // as GPU Nodes).
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;
  }

  template<class Scalar, class Ordinal, class Node, class Allocator>
  MklSparseOps<Scalar,Ordinal,Node,Allocator>::~MklSparseOps() {
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  std::string
  MklSparseOps<Scalar,Ordinal,Node,Allocator>::description() const
  {
    using Teuchos::TypeNameTraits;
    std::ostringstream os;
    os <<  "Kokkos::MklSparseOps<"
       << "Scalar=" << TypeNameTraits<Scalar>::name()
       << ", Ordinal=" << TypeNameTraits<Ordinal>::name()
       << ", Node=" << TypeNameTraits<Node>::name()
       << ", Allocator=" << TypeNameTraits<Allocator>::name()
       << ">";
    return os.str();
  }

  // Partial specialization for Scalar=void avoids TypeNameTraits<void>.
  // template <class Ordinal, class Node, class Allocator>
  // std::string
  // MklSparseOps<void,Ordinal,Node,Allocator>::description() const
  // {
  //   using Teuchos::TypeNameTraits;
  //   std::ostringstream os;
  //   os <<  "Kokkos::MklSparseOps<"
  //      << "Scalar=void"
  //      << ", Ordinal=" << TypeNameTraits<Ordinal>::name()
  //      << ", Node=" << TypeNameTraits<Node>::name()
  //      << ", Allocator=" << TypeNameTraits<Allocator>::name()
  //      << ">";
  //   return os.str();
  // }


  template <class Scalar, class Ordinal, class Node, class Allocator>
  void
  MklSparseOps<Scalar,Ordinal,Node,Allocator>::
  describe (Teuchos::FancyOStream& out,
            const Teuchos::EVerbosityLevel verbLevel) const
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

    // Interpret the default verbosity level as VERB_MEDIUM.
    const EVerbosityLevel vl =
      (verbLevel == VERB_DEFAULT) ? VERB_MEDIUM : verbLevel;

    if (vl == VERB_NONE) {
      return;
    }
    else if (includesVerbLevel (vl, VERB_LOW)) { // vl >= VERB_LOW
      out << this->description();

      if (! includesVerbLevel (vl, VERB_MEDIUM)) {
        out << endl;
      }
      else if (includesVerbLevel (vl, VERB_MEDIUM)) { // vl >= VERB_MEDIUM
        out << ":" << endl;
        OSTab tab1 (rcpFromRef (out));

        out << "isInitialized_: " << (isInitialized_ ? "true" : "false") << endl;
        if (isInitialized_) {
          out << "numRows_: " << numRows_ << endl
              << "isEmpty_: " << (isEmpty_ ? "true" : "false") << endl;
          {
            size_t numEntries = 0;
            if (ptr_.size() > 0) {
              numEntries = ptr_[ptr_.size()-1];
            }
            else {
              numEntries = 0;
            }
            out << "numEntries: " << numEntries << endl;
          }
          matrixDescriptor_.print (out);

          if (includesVerbLevel (vl, VERB_EXTREME)) { // vl >= VERB_EXTREME
            // Only print out all the sparse matrix's data in
            // extreme verbosity mode.
            out << "ptrs: [";
            if (ptr_.size() > 0) {
              std::copy (ptr_.begin(), ptr_.end(),
                         std::ostream_iterator<Ordinal> (out, ", "));
            }
            out << "]" << endl << "ind_: [";
            std::copy (ind_.begin(), ind_.end(),
                       std::ostream_iterator<Ordinal> (out, ", "));
            out << "]" << endl << "val_: [";
            std::copy (val_.begin(), val_.end(),
                       std::ostream_iterator<Scalar> (out, ", "));
            out << "]" << endl;
          } // vl >= VERB_EXTREME
        } // if is initialized
      } // vl >= VERB_MEDIUM
    } // vl >= VERB_LOW
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  void MklSparseOps<Scalar,Ordinal,Node,Allocator>::
  setGraphAndMatrix (const RCP<const DefaultCrsGraph<Ordinal,Node> > &opgraph,
                     const RCP<const DefaultCrsMatrix<Scalar,Ordinal,Node> > &opmatrix)
  {
    using Teuchos::arcp;
    using Teuchos::as;

    std::string tfecfFuncName("setGraphAndMatrix(uplo,diag,graph,matrix)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isInitialized_ == true, std::runtime_error,
      " The sparse kernels object is already initialized.");

    ArrayRCP<const MKL_INT> ptr = opgraph->getSmallPointers ();
    ArrayRCP<const MKL_INT> ind = opgraph->getIndices ();
    ArrayRCP<const Scalar> val = opmatrix->getValues ();
    const Ordinal numRows = opgraph->getNumRows ();
    const Ordinal numCols = opgraph->getNumCols ();

    // Verify the input data before setting internal state.

    TEUCHOS_TEST_FOR_EXCEPTION(
      ptr.is_null () && ! opgraph->getPointers ().is_null (),
      std::invalid_argument,
      "Kokkos::MklSparseOps::setGraphAndMatrix: "
      "The given graph and matrix have more entries than can fit in the MKL_INT="
      << Teuchos::TypeNameTraits<MKL_INT>::name() << " type.  MKL requires "
      "that the row offsets array, and therefore the number of entries in the "
      "graph, has the same type as the column indices.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ptr.is_null (), std::invalid_argument,
      ": The input graph's row offsets array (the return value of its "
      "getSmallPointers() method) is null.  This probably means that has not "
      "yet been finalized.  You can finalize the graph by calling "
      "finalizeGraph() or finalizeGraphAndMatrix().");
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
    numCols_ = numCols;

    if (opgraph->isEmpty() || numRows_ == 0) {
      isEmpty_ = true;
      ind_ = null;
      val_ = null;
    }
    else {
      isEmpty_ = false;
      // MKL requires one-based row offsets and column indices when
      // using column-major multivectors with > 1 column.  However,
      // that is not a common case.  We optimize for the common case
      // by preserving zero-based indices for now, and only copying if
      // multiply() or solve() encounters input with multiple columns.
      ptr_ = ptr;
      ind_ = ind;
      val_ = val;
    }
    opgraph->getMatDesc (tri_uplo_, unit_diag_);

    const bool triangular = (tri_uplo_ != Teuchos::UNDEF_TRI);
    // Index base is 0 for now.  We will convert to 1-based indices
    // (which MKL needs when using column-major multivectors, as we
    // do) on demand, if users invoke multiply() or solve() with
    // multiple right-hand sides.
    const int indexBase = 0;
    matrixDescriptor_ =
      Mkl::MatrixDescriptor (false, false, triangular, false, false,
                             tri_uplo_, unit_diag_, indexBase);
    isInitialized_ = true;
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  void MklSparseOps<Scalar,Ordinal,Node,Allocator>::
  convertIndexBase (const int newIndexBase)
  {
    // This is a logic error because users are not allowed to invoke
    // convertIndexBase(), and the only index bases we should be using
    // are 0 or 1.
    TEUCHOS_TEST_FOR_EXCEPTION(newIndexBase != 0 && newIndexBase != 1,
      std::logic_error, "convertIndexBase only allows an index base of 0 or 1, "
      "but you provided an index base of " << newIndexBase << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(! isInitialized_, std::logic_error, "Converting "
      "the index base is only valid if the local sparse ops have been "
      "initialized by calling setGraphAndMatrix().");
    TEUCHOS_TEST_FOR_EXCEPTION(ptr_.is_null (), std::logic_error, "The local "
      "sparse ops say that they have been initialized, but the row offsets "
      "array is null.  This is invalid even for an empty matrix.");

    const int oldIndexBase = matrixDescriptor_.getIndexBase ();
    if (newIndexBase != oldIndexBase) {
      const MKL_INT nnz = ptr_[numRows_] - oldIndexBase;

      TEUCHOS_TEST_FOR_EXCEPTION(nnz < 0, std::logic_error, "The local sparse "
        "matrix reports that its number of entries nnz = " << nnz << " is "
        "negative.");

      const MKL_INT* const rawPtrIn = ptr_.getRawPtr ();
      MKL_INT* rawPtrOut = NULL;
      try {
        rawPtrOut = new int [numRows_+1];
        node_->parallel_for (0, numRows_+1, ConvertRowPtrs (rawPtrOut, rawPtrIn, newIndexBase));
        ptr_ = null; // We don't need the original row offsets anymore.
      }
      catch (...) {
        if (rawPtrOut != NULL) {
          delete [] rawPtrOut;
        }
        throw;
      }
      ptr_ = arcp<const MKL_INT> (const_cast<const MKL_INT*> (rawPtrOut),
                                  0, numRows_+1, true);

      const MKL_INT* const rawIndIn = ind_.getRawPtr();
      MKL_INT* rawIndOut = NULL;
      try {
        rawIndOut = new int [nnz];
        ConvertColInds wdp (rawIndOut, rawPtrOut, rawIndIn, newIndexBase);
        node_->parallel_for (0, numRows_, wdp);
        ind_ = null; // We don't need the original column indices anymore.
      }
      catch (...) {
        if (rawIndOut != NULL) {
          delete [] rawIndOut;
        }
        throw;
      }
      ind_ = arcp<const MKL_INT> (const_cast<const MKL_INT*> (rawIndOut),
                                  0, nnz, true);
      const bool triangular = (tri_uplo_ != Teuchos::UNDEF_TRI);
      matrixDescriptor_ =
        Mkl::MatrixDescriptor (false, false, triangular, false, false,
                               tri_uplo_, unit_diag_, newIndexBase);
    }
  }

  // MKL only implements sparse triangular solve for DomainScalar ==
  // RangeScalar == Scalar, and for OffsetType == Ordinal.  Here, we
  // specialize for the general (types-not-equal) case by throwing a
  // run-time exception.  That's probably less confusing for the end
  // user than letting the compiler spew template errors all over
  // their screen.  Also, we can take the opportunity to teach users
  // how to access fall-back sparse kernels.
  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class DomainScalar, class RangeScalar, class OffsetType>
  void MklSparseOps<Scalar,Ordinal,Node,Allocator>::
  solvePrivate (Teuchos::ETransp trans,
                const RangeScalar& alpha,
                const MultiVector<DomainScalar,Node> &Y,
                MultiVector< RangeScalar,Node> &X) const
  {
    // using Teuchos::TypeNameTraits;
    // std::string tfecfFuncName("solve(trans,Y,X)");
    // TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    //   true, std::logic_error, " Intel's Math Kernel Library (MKL) does not "
    //   "define sparse triangular solve for this combination of (Matrix)Scalar="
    //   << TypeNameTraits<Scalar>::name() << ", DomainScalar="
    //   << TypeNameTraits<DomainScalar>::name() << ", RangeScalar="
    //   << TypeNameTraits<RangeScalar>::name() << ", Ordinal="
    //   << TypeNameTraits<Ordinal>::name() << ", OffsetType="
    //   << TypeNameTraits<OffsetType>::name() << ".  Please refer to the "
    //   "documentation of the Kokkos::MklSparseOps class to learn how to access a "
    //   "fall-back implementation of sparse kernels in that case.");

    // We know at this point that DomainScalar == RangeScalar ==
    // Scalar, so we only need to use Scalar inside this method.
    // Likewise, we know that OffsetType == Ordinal.
    using Teuchos::TypeNameTraits;
    typedef Mkl::RawSparseKernels<Scalar, Ordinal> kernels_type;
    std::string tfecfFuncName ("solve(trans,Y,X)");

    const Teuchos::EDiag diag = matrixDescriptor_.getDiag ();
    if (numRows_ == 0) {
      // null op
    }
    else if (isEmpty_) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        diag != Teuchos::UNIT_DIAG, std::runtime_error,
        ": Solve with an empty matrix is only valid with an implicit unit diagonal.");
      // solve I * X = Y for X = Y
      DefaultArithmetic<MultiVector<Scalar, Node> >::Assign (X, Y);
    }
    else {
      const char* const transa =
        trans == Teuchos::NO_TRANS ? "N" :
        (trans == Teuchos::TRANS ? "T" : "C");
      const Ordinal M = numRows_;
      const Ordinal N = Y.getNumCols ();

      // Convert the sparse matrix's index base from 0 to 1 in order
      // to support MKL's requirement of 1-based indices when using
      // column-major multivectors.
      if (N != 1 && getIndexBase () != 1) {
        typedef MklSparseOps<Scalar, Ordinal, Node, Allocator> this_type;
        // mfh 23 Aug 2012: Don't try this at home, kids.
        const_cast<this_type*> (this)->convertIndexBase (1);
      }

      // Don't keep this pointer around.  It's not valid if you change
      // the matrix descriptor.
      const char* const descr = matrixDescriptor_.getRawPtr ();
      const Scalar* const val = val_.getRawPtr ();
      const Ordinal* const ind = ind_.getRawPtr ();
      const Ordinal* const ptr = ptr_.getRawPtr ();
      Scalar* const X_raw = X.getValuesNonConst ().getRawPtr ();
      const Scalar* const Y_raw = Y.getValues ().getRawPtr ();
      const Ordinal LDX = X.getStride ();
      const Ordinal LDY = Y.getStride ();

      if (N == 1) {
        kernels_type::csrsv (transa, M, alpha, descr, val, ind, ptr, &ptr[1],
                             Y_raw, X_raw);
      }
      else {
        kernels_type::csrsm (transa, M, N, alpha, descr, val, ind, ptr, &ptr[1],
                             Y_raw, LDY, X_raw, LDX);
      }
    }
  }

  // MKL only implements sparse triangular solve for DomainScalar ==
  // RangeScalar == Scalar, and for OffsetType == Ordinal.  Here, we
  // specialize for that valid case.
#if 0
  template <class Node, class Allocator>
  template <>
  void MklSparseOps<double, MKL_INT, Node, Allocator>::
  solvePrivate<double, double, MKL_INT> (Teuchos::ETransp trans,
                                         const double& alpha,
                                         const MultiVector<double, Node> &Y,
                                         MultiVector<double, Node> &X) const
  {
    // We know at this point that DomainScalar == RangeScalar ==
    // Scalar, so we only need to use Scalar inside this method.
    // Likewise, we know that OffsetType == Ordinal.
    using Teuchos::TypeNameTraits;
    typedef Mkl::RawSparseKernels<double, MKL_INT> kernels_type;
    std::string tfecfFuncName ("solve(trans,Y,X)");

    // Get the array of row offsets.
    ArrayRCP<const MKL_INT> ptrs = ptr_;

    const Teuchos::EDiag diag = matrixDescriptor_.getDiag ();
    if (numRows_ == 0) {
      // null op
    }
    else if (isEmpty_) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        diag != Teuchos::UNIT_DIAG, std::runtime_error,
        ": Solve with an empty matrix is only valid with an implicit unit diagonal.");
      // solve I * X = Y for X = Y
      DefaultArithmetic<MultiVector<Scalar, Node> >::Assign (X, Y);
    }
    else {
      const char* const transa =
        trans == Teuchos::NO_TRANS ? "N" :
        (trans == Teuchos::TRANS ? "T" : "C");
      const MKL_INT M = numRows_;
      const MKL_INT N = Y.getNumCols ();

      // Convert the sparse matrix's index base from 0 to 1 in order
      // to support MKL's requirement of 1-based indices when using
      // column-major multivectors.
      if (N != 1 && getIndexBase () != 1) {
        typedef MklSparseOps<double, MKL_INT, Node, Allocator> this_type;
        // mfh 23 Aug 2012: Don't try this at home, kids.
        const_cast<this_type*> (this)->convertIndexBase (1);
      }

      // Don't keep this pointer around.  It's not valid if you change
      // the matrix descriptor.
      const char* const descr = matrixDescriptor_.getRawPtr ();
      const double* const val = vals.getRawPtr ();
      const MKL_INT* const ind = ind_.getRawPtr ();
      const MKL_INT* const ptr = ptrs.getRawPtr ();
      const double* const X_raw = X.getValues ().getRawPtr ();
      double* const Y_raw = Y.getValues ().getRawPtr ();
      const MKL_INT LDX = X.getStride ();
      const MKL_INT LDY = Y.getStride ();

      if (N == 1) {
        kernels_type::csrsv (transa, M, alpha, descr, val, ind, ptr, &ptr[1],
                             X_raw, Y_raw);
      }
      else {
        kernels_type::csrsv (transa, M, N, alpha, descr, val, ind, ptr, &ptr[1],
                             X_raw, LDX, Y_raw, LDY);
      }
    }
  }
#endif // 0

  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class DomainScalar, class RangeScalar>
  void
  MklSparseOps<Scalar,Ordinal,Node,Allocator>::
  solve (Teuchos::ETransp trans,
         const MultiVector<DomainScalar,Node> &Y,
         MultiVector< RangeScalar,Node> &X) const
  {
    std::string tfecfFuncName("solve(trans,Y,X)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == false,
        std::runtime_error, " The local sparse kernels object was not fully initialized."
    );
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        (size_t)X.getNumCols() != (size_t)Y.getNumCols(),
        std::runtime_error, " Left hand side and right hand side multivectors have differing numbers of vectors."
    );
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        (size_t)X.getNumRows() < (size_t)numRows_,
        std::runtime_error, " Left-hand-side multivector does not have enough rows. "
                            "Likely cause is that the column map was not provided to "
                            "the Tpetra::CrsMatrix in the case of an implicit unit diagonal."
    );

    const RangeScalar alpha = Teuchos::ScalarTraits<RangeScalar>::one();
    solvePrivate<DomainScalar, RangeScalar, Ordinal> (trans, alpha, Y, X);
  }


  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class DomainScalar, class RangeScalar>
  void
  MklSparseOps<Scalar,Ordinal,Node,Allocator>::
  multiply (Teuchos::ETransp trans,
            RangeScalar alpha,
            const MultiVector<DomainScalar,Node> &X,
            MultiVector<RangeScalar ,Node> &Y) const
  {
    const RangeScalar beta = Teuchos::ScalarTraits<RangeScalar>::zero ();
    this->multiply<DomainScalar, RangeScalar> (trans, alpha, X, beta, Y);
  }


  // MKL only implements sparse matrix-(multi)vector multiply for
  // DomainScalar == RangeScalar == Scalar, and for OffsetType ==
  // Ordinal.  Here, we specialize for the general (types-not-equal)
  // case by throwing a run-time exception.  That's probably less
  // confusing for the end user than letting the compiler spew
  // template errors all over their screen.  Also, we can take the
  // opportunity to teach users how to access fall-back sparse
  // kernels.
  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class DomainScalar, class RangeScalar, class OffsetType>
  void MklSparseOps<Scalar,Ordinal,Node,Allocator>::
  multiplyPrivate (Teuchos::ETransp trans,
                   RangeScalar alpha,
                   const MultiVector<DomainScalar,Node> &X,
                   RangeScalar beta,
                   MultiVector<RangeScalar,Node> &Y) const
  {
    // using Teuchos::TypeNameTraits;
    // std::string tfecfFuncName("solve(trans,Y,X)");
    // TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    //   true, std::logic_error, " Intel's Math Kernel Library (MKL) does not "
    //   "define sparse matrix-(multi)vector multiply for this combination of "
    //   "(Matrix)Scalar=" << TypeNameTraits<Scalar>::name() << ", DomainScalar="
    //   << TypeNameTraits<DomainScalar>::name() << ", RangeScalar="
    //   << TypeNameTraits<RangeScalar>::name() << ", Ordinal="
    //   << TypeNameTraits<Ordinal>::name() << ", OffsetType="
    //   << TypeNameTraits<OffsetType>::name() << ".  Please refer to the "
    //   "documentation of the Kokkos::MklSparseOps class to learn how to access a "
    //   "fall-back implementation of sparse kernels in that case.");

    // We know at this point that DomainScalar == RangeScalar ==
    // Scalar, so we only need to use Scalar inside this method.
    // Likewise, we know that OffsetType == Ordinal.
    using Teuchos::TypeNameTraits;
    typedef Mkl::RawSparseKernels<Scalar, Ordinal> kernels_type;
    std::string tfecfFuncName("multiply(trans,alpha,X,beta,Y)");

    if (numRows_ == 0) {
      // null op
    }
    else if (isEmpty_) {
      // Y := alpha * 0 * X + beta * Y = beta * Y.
      DefaultArithmetic<MultiVector<Scalar, Node> >::Scale (Y, beta);
    }
    else {
      const char* const transa =
        trans == Teuchos::TRANS ? "T" :
        (trans == Teuchos::CONJ_TRANS ? "C" : "N");
      const Ordinal M = numRows_;
      const Ordinal N = Y.getNumCols ();
      const Ordinal K = X.getNumRows ();

      // Convert the sparse matrix's index base from 0 to 1 in order
      // to support MKL's requirement of 1-based indices when using
      // column-major multivectors.
      if (N != 1 && matrixDescriptor_.getIndexBase () != 1) {
        typedef MklSparseOps<Scalar, Ordinal, Node, Allocator> this_type;
        // mfh 23 Aug 2012: Don't try this at home, kids.
        const_cast<this_type*> (this)->convertIndexBase (1);
      }

      // Don't keep this pointer around.  It's not valid if you change
      // the matrix descriptor.
      const char* const descr = matrixDescriptor_.getRawPtr ();
      const Scalar* const val = val_.getRawPtr ();
      const Ordinal* const ind = ind_.getRawPtr ();
      const Ordinal* const ptr = ptr_.getRawPtr ();
      const Scalar* const X_raw = X.getValues ().getRawPtr ();
      Scalar* const Y_raw = Y.getValuesNonConst ().getRawPtr ();
      const Ordinal LDX = X.getStride ();
      const Ordinal LDY = Y.getStride ();

      if (N == 1) {
        kernels_type::csrmv (transa, M, K, alpha, descr,
                             val, ind, ptr, &ptr[1],
                             X_raw, beta, Y_raw);
      }
      else {
        kernels_type::csrmm (transa, M, N, K, alpha, descr,
                             val, ind, ptr, &ptr[1],
                             X_raw, LDX, beta, Y_raw, LDY);
      }
    }
  }

#if 0
  // MKL only implements sparse matrix-(multi)vector multiply for
  // DomainScalar == RangeScalar == Scalar, and for OffsetType ==
  // Ordinal.  Here, we specialize for that valid case.
  template <class Node, class Allocator>
  template <>
  void MklSparseOps<double, MKL_INT, Node, Allocator>::
  multiplyPrivate<double, double, MKL_INT> (Teuchos::ETransp trans,
                                            double alpha,
                                            const MultiVector<double, Node> &X,
                                            double beta,
                                            MultiVector<double, Node> &Y) const
  {
    // We know at this point that DomainScalar == RangeScalar ==
    // Scalar, so we only need to use Scalar inside this method.
    // Likewise, we know that OffsetType == Ordinal.
    using Teuchos::TypeNameTraits;
    typedef Mkl::RawSparseKernels<double, MKL_INT> kernels_type;
    std::string tfecfFuncName("multiply(trans,alpha,X,beta,Y)");

    // Get the array of row offsets.
    ArrayRCP<const MKL_INT> ptrs = ptr_;

    if (numRows_ == 0) {
      // null op
    }
    else if (isEmpty_) {
      // Y := alpha * 0 * X + beta * Y = beta * Y.
      DefaultArithmetic<MultiVector<double, Node> >::Scale (Y, beta);
    }
    else {
      const char* const transa =
        trans == Teuchos::TRANS ? "T" :
        (trans == Teuchos::CONJ_TRANS ? "C" : "N");
      const MKL_INT M = numRows_;
      const MKL_INT N = Y.getNumCols ();
      const MKL_INT K = X.getNumRows ();

      // Convert the sparse matrix's index base from 0 to 1 in order
      // to support MKL's requirement of 1-based indices when using
      // column-major multivectors.
      if (N != 1 && matrixDescriptor_.getIndexBase () != 1) {
        typedef MklSparseOps<double, MKL_INT, Node, Allocator> this_type;
        // mfh 23 Aug 2012: Don't try this at home, kids.
        const_cast<this_type*> (this)->convertIndexBase (1);
      }

      // Don't keep this pointer around.  It's not valid if you change
      // the matrix descriptor.
      const char* const descr = getMatrixDescriptor ();
      const double* const val = vals.getRawPtr ();
      const MKL_INT* const ind = ind_.getRawPtr ();
      const MKL_INT* const ptr = ptrs.getRawPtr ();
      const double* const X_raw = X.getValues ().getRawPtr ();
      double* const Y_raw = Y.getValues ().getRawPtr ();
      const MKL_INT LDX = X.getStride ();
      const MKL_INT LDY = Y.getStride ();

      if (N == 1) {
        kernels_type::csrmv (transa, M, K, alpha, descr,
                             val, ind, ptr, &ptr[1],
                             X_raw, beta, Y_raw);
      }
      else {
        kernels_type::csrmm (transa, M, N, K, alpha, descr,
                             val, ind, ptr, &ptr[1],
                             X_raw, LDX, beta, Y_raw, LDY);
      }
    }
  }
#endif // 0

  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class DomainScalar, class RangeScalar>
  void MklSparseOps<Scalar,Ordinal,Node,Allocator>::
  multiply (Teuchos::ETransp trans,
            RangeScalar alpha,
            const MultiVector<DomainScalar,Node> &X,
            RangeScalar beta,
            MultiVector<RangeScalar,Node> &Y) const
  {
    std::string tfecfFuncName("multiply(trans,alpha,X,beta,Y)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isInitialized_ == false,
      std::runtime_error,
      " Sparse kernels have not yet been initialized.  Call the MklSparseOps' "
      "setGraphAndMatrix() method with a valid graph and matrix in order to "
      "initialize correctly, before calling multiply() or solve().");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      X.getNumCols() != Y.getNumCols(),
      std::invalid_argument,
      " X and Y do not have the same number of columns.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      Y.getNumRows() != Teuchos::as<size_t> (numRows_),
      std::invalid_argument,
      " Y has " << Y.getNumRows() << " rows, but the matrix has "
      << numRows_ << " rows.");

    multiplyPrivate<DomainScalar, RangeScalar, Ordinal> (trans, alpha, X, beta, Y);
  }

} // namespace Kokkos

#endif /* __Kokkos_MklSparseOps_hpp */

