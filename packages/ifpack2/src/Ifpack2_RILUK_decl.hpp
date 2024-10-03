// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Ifpack2_RILUK_decl.hpp
/// \brief Declaration of RILUK interface

#ifndef IFPACK2_CRSRILUK_DECL_HPP
#define IFPACK2_CRSRILUK_DECL_HPP

#include "KokkosSparse_spiluk.hpp"

#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Details_CanChangeMatrix.hpp"
#include "Tpetra_CrsMatrix_decl.hpp"
#include "Ifpack2_ScalingType.hpp"
#include "Ifpack2_IlukGraph.hpp"
#include "Ifpack2_LocalSparseTriangularSolver_decl.hpp"

#include <memory>
#include <type_traits>

namespace Teuchos {
  class ParameterList; // forward declaration
}

namespace Ifpack2 {

/** \class RILUK
\brief ILU(k) factorization of a given Tpetra::RowMatrix.
\tparam MatrixType A specialization of Tpetra::RowMatrix.

This class implements a "relaxed" incomplete ILU (ILU) factorization
with level k fill.  It is based upon the ILU algorithms outlined in
Yousef Saad's "Iterative Methods for Sparse Linear Systems", 2nd
edition, Chapter 10.

\section Ifpack2_RILUK_Parameters Parameters

For a complete list of valid parameters, see the documentation of setParameters().

The computed factorization is a function of several parameters:
<ul>
<li>
The graph structure (sparsity pattern) of the matrix: All fill is
derived from the original matrix nonzero structure.  Level zero fill
is defined as the original matrix pattern (nonzero structure), even if
the matrix value at an entry is stored as a zero. (Thus it is possible
to add entries to the ILU factors by adding zero entries to the
original matrix.)
</li>

<li>
Level of fill: Starting with the original matrix pattern as level
fill of zero, the next level of fill is determined by analyzing the
graph of the previous level and determining nonzero fill that is a
result of combining entries that were from previous level only (not
the current level).  This rule limits fill to entries that are direct
decendents from the previous level graph.  Fill for level k is
determined by applying this rule recursively.  For sufficiently large
values of k, the fill would eventually be complete and an exact LU
factorization would be computed.
</li>

<li>
Fraction of relaxation: Ifpack2::RILUK computes the ILU factorization
row-by-row.  As entries at a given row are computed, some number of
them will be dropped because they do match the prescribed sparsity
pattern.  The relaxation factor determines how these dropped values
will be handled.  If the factor is zero, then these extra entries will
by dropped.  This is a classical ILU approach.  If the RelaxValue is
1, then the sum of the extra entries will be added to the diagonal.
This is a classical Modified ILU (MILU) approach.  If RelaxValue is
between 0 and 1, then the factor times the sum of extra entries will
be added to the diagonal.

For most situations, the relaxation factor should be set to zero.  For
certain kinds of problems, e.g., reservoir modeling, there is a
conservation principle involved such that any operator should obey a
zero row-sum property.  MILU was designed for these cases and you
should set the relaxation factor to 1.  For other situations, setting
RelaxValue to some nonzero value may improve the stability of
factorization, and can be used if the computed ILU factors are poorly
conditioned.
</li>

<li>
Diagonal perturbation: Prior to computing the factorization, it is
possible to modify the diagonal entries of the matrix for which the
factorization will be computing.  If the absolute and relative
perturbation values are zero and one, respectively, the factorization
will be compute for the original user matrix A.  Otherwise, the
factorization will computed for a matrix that differs from the
original user matrix in the diagonal values only.  Below we discuss
the details of diagonal perturbations.
</li>

</ul>

\section Ifpack2_RILUK_GlobalOrdering An important note about ordering

Note that the factorization is calculated based upon local ordering.   This means
that the ordering of the GIDs in the row map is ignored.
Initial entries in \f$L\f$, the strictly lower triangular part of A, and \f$U\f$, the strictly upper
triangular part of A, are given by

\f$L(i,j) = A(i,j)\f$ if \f$j < i\f$, for local IDs \f$i\f$ and \f$j\f$, even if GID\f$(j)\f$ \f$>\f$ GID\f$(i)\f$,

and

\f$U(i,j) = A(i,j)\f$ if \f$i < j\f$, for local IDs \f$i\f$ and \f$j\f$, even if GID\f$(j)\f$ \f$<\f$ GID\f$(i)\f$.

In particular, if the row map GIDs are not in ascending
order on processor, then the incomplete factors will be different than those produced by ILU(k) using global IDs.
If the row map GIDs are in ascending order, then the factors produced based on LID and GID ordering are the same.

\section Ifpack2_RILUK_CondEst Estimating preconditioner condition numbers

For ill-conditioned matrices, we often have difficulty computing
usable incomplete factorizations.  The most common source of problems
is that the factorization may encounter a small or zero pivot.  In
that case, the factorization may fail.  Even if the factorization
succeeds, the factors may be so poorly conditioned that use of them in
the iterative phase produces meaningless results.  Before we can fix
this problem, we must be able to detect it.  To this end, we use a
simple but effective condition number estimate for \f$(LU)^{-1}\f$.

The condition number of a matrix \f$B\f$, called \f$cond_p(B)\f$, is
defined as \f$cond_p(B) = \|B\|_p\|B^{-1}\|_p\f$ in some appropriate
norm \f$p\f$.  \f$cond_p(B)\f$ gives some indication of how many
accurate floating point digits can be expected from operations
involving the matrix and its inverse.  A condition number approaching
the accuracy of a given floating point number system, about 15 decimal
digits in IEEE double precision, means that any results involving
\f$B\f$ or \f$B^{-1}\f$ may be meaningless.

The \f$\infty\f$-norm of a vector \f$y\f$ is defined as the maximum of
the absolute values of the vector entries, and the \f$\infty\f$-norm
of a matrix C is defined as \f$\|C\|_\infty = \max_{\|y\|_\infty = 1}
\|Cy\|_\infty\f$.  A crude lower bound for the \f$cond_\infty(C)\f$ is
\f$\|C^{-1}e\|_\infty\f$ where \f$e = (1, 1, \ldots, 1)^T\f$.  It is a
lower bound because \f$cond_\infty(C) = \|C\|_\infty\|C^{-1}\|_\infty
\ge \|C^{-1}\|_\infty \ge |C^{-1}e\|_\infty\f$.

For our purposes, we want to estimate \f$cond_\infty(LU)\f$, where
\f$L\f$ and \f$U\f$ are our incomplete factors.  Edmond in his
Ph.D. thesis demonstrates that \f$\|(LU)^{-1}e\|_\infty\f$ provides an
effective estimate for \f$cond_\infty(LU)\f$.  Furthermore, since
finding \f$z\f$ such that \f$LUz = y\f$ is a basic kernel for applying
the preconditioner, computing this estimate of \f$cond_\infty(LU)\f$
is performed by setting \f$y = e\f$, calling the solve kernel to
compute \f$z\f$ and then computing \f$\|z\|_\infty\f$.

\section Ifpack2_RILUK_DiagPerturb A priori diagonal perturbations

If we detect using the above method that our factorization is too
ill-conditioned, we can improve the conditioning by perturbing the
matrix diagonal and restarting the factorization using this more
diagonally dominant matrix.  In order to apply perturbation, prior to
starting the factorization, we compute a diagonal perturbation of our
matrix \f$A\f$ and perform the factorization on this perturbed matrix.
The overhead cost of perturbing the diagonal is minimal since the
first step in computing the incomplete factors is to copy the matrix
\f$A\f$ into the memory space for the incomplete factors.  We simply
compute the perturbed diagonal at this point.

The actual perturbation values we use are the diagonal values \f$(d_1,
d_2, \ldots, d_n)\f$ with \f$d_i = sgn(d_i)\alpha + d_i\rho\f$,
\f$i=1, 2, \ldots, n\f$, where \f$n\f$ is the matrix dimension and
\f$sgn(d_i)\f$ returns the sign of the diagonal entry.  This has the
effect of forcing the diagonal values to have minimal magnitude of
\f$\alpha\f$ and to increase each by an amount proportional to
\f$\rho\f$, and still keep the sign of the original diagonal entry.

\section Ifpack2_RILUK_Phases Phases of computation

Every Ifpack2 preconditioner has the following phases of computation:
<ol>
  <li> initialize() </li>
  <li> compute() </li>
  <li> apply() </li>
</ol>

RILUK constructs the symbolic incomplete factorization (that is, the
structure of the incomplete factors) in the initialize() phase.  It
computes the numerical incomplete factorization (that is, it fills in
the factors' entries with their correct values) in the compute()
phase.  The apply() phase applies the incomplete factorization to a
given multivector using two triangular solves.

\section Ifpack2_RILUK_Measuring Measuring performance

Each RILUK object keeps track of both the time required for various
operations, and the number of times those operations have been applied
for that object.  The operations tracked include:
  - initialize() (via getNumInitialize() and getInitializeTime())
  - compute() (via getNumCompute() and getComputeTime())
  - apply() (via getNumApply() and getApplyTime())

The <tt>getNum*</tt> methods return the number of times that operation
was called.  The <tt>get*Time</tt> methods return the total number of
seconds spent in <i>all</i> invocations of that operation.  For
example, getApplyTime() returns the number of seconds spent in all
apply() calls.  For an average time per apply() call, divide by
getNumApply(), the total number of calls to apply().
*/
template<class MatrixType>
class RILUK:
    virtual public Ifpack2::Preconditioner<typename MatrixType::scalar_type,
                                           typename MatrixType::local_ordinal_type,
                                           typename MatrixType::global_ordinal_type,
                                           typename MatrixType::node_type>,
    virtual public Ifpack2::Details::CanChangeMatrix<Tpetra::RowMatrix<typename MatrixType::scalar_type,
                                                                       typename MatrixType::local_ordinal_type,
                                                                       typename MatrixType::global_ordinal_type,
                                                                       typename MatrixType::node_type> >
{
 public:
  //! The type of the entries of the input MatrixType.
  typedef typename MatrixType::scalar_type scalar_type;

  //! The type of local indices in the input MatrixType.
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;

  //! The type of global indices in the input MatrixType.
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;

  //! The Node type used by the input MatrixType.
  typedef typename MatrixType::node_type node_type;

  //! The type of the magnitude (absolute value) of a matrix entry.
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  //! The Kokkos device type of the input MatrixType.
  typedef typename node_type::device_type device_type;

  //! The Kokkos execution space of the input MatrixType.
  typedef typename node_type::execution_space execution_space;

  //! Tpetra::RowMatrix specialization used by this class.
  typedef Tpetra::RowMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> row_matrix_type;


  static_assert(std::is_same<MatrixType, row_matrix_type>::value, "Ifpack2::RILUK: The template parameter MatrixType must be a Tpetra::RowMatrix specialization.  Please don't use Tpetra::CrsMatrix (a subclass of Tpetra::RowMatrix) here anymore.");

  //! Tpetra::CrsMatrix specialization used by this class for representing L and U.
  typedef Tpetra::CrsMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> crs_matrix_type;

  //! Scalar type stored in Kokkos::Views (CrsMatrix and MultiVector)
  typedef typename crs_matrix_type::impl_scalar_type impl_scalar_type;

  template <class NewMatrixType> friend class RILUK;

  typedef typename crs_matrix_type::global_inds_host_view_type global_inds_host_view_type;
  typedef typename crs_matrix_type::local_inds_host_view_type local_inds_host_view_type;
  typedef typename crs_matrix_type::values_host_view_type values_host_view_type;


  typedef typename crs_matrix_type::nonconst_global_inds_host_view_type nonconst_global_inds_host_view_type;
  typedef typename crs_matrix_type::nonconst_local_inds_host_view_type nonconst_local_inds_host_view_type;
  typedef typename crs_matrix_type::nonconst_values_host_view_type nonconst_values_host_view_type;


  //@}
  //! \name Implementation of Kokkos Kernels ILU(k).
  //@{

  typedef typename crs_matrix_type::local_matrix_device_type local_matrix_device_type;
  typedef typename local_matrix_device_type::StaticCrsGraphType::row_map_type lno_row_view_t;
  typedef typename local_matrix_device_type::StaticCrsGraphType::entries_type lno_nonzero_view_t;
  typedef typename local_matrix_device_type::values_type scalar_nonzero_view_t;
  typedef typename local_matrix_device_type::StaticCrsGraphType::device_type::memory_space TemporaryMemorySpace;
  typedef typename local_matrix_device_type::StaticCrsGraphType::device_type::memory_space PersistentMemorySpace;
  typedef typename local_matrix_device_type::StaticCrsGraphType::device_type::execution_space HandleExecSpace;
  typedef typename KokkosKernels::Experimental::KokkosKernelsHandle
    <typename lno_row_view_t::const_value_type, typename lno_nonzero_view_t::const_value_type, typename scalar_nonzero_view_t::value_type,
    HandleExecSpace, TemporaryMemorySpace,PersistentMemorySpace > kk_handle_type;
  typedef Ifpack2::IlukGraph<Tpetra::CrsGraph<local_ordinal_type, global_ordinal_type, node_type>, kk_handle_type> iluk_graph_type;

  /// \brief Constructor that takes a Tpetra::RowMatrix.
  ///
  /// \param A_in [in] The input matrix.
  RILUK (const Teuchos::RCP<const row_matrix_type>& A_in);

  /// \brief Constructor that takes a Tpetra::CrsMatrix.
  ///
  /// This constructor exists to avoid "ambiguous constructor"
  /// warnings.  It does the same thing as the constructor that takes
  /// a Tpetra::RowMatrix.
  ///
  /// \param A_in [in] The input matrix.
  RILUK (const Teuchos::RCP<const crs_matrix_type>& A_in);

 private:
  /// \brief Copy constructor: declared private but not defined, so
  ///   that calling it is syntactically forbidden.
  RILUK (const RILUK<MatrixType> & src);

 public:
  //! Destructor (declared virtual for memory safety).
  virtual ~RILUK ();

  /// Set parameters for the incomplete factorization.
  ///
  /// This preconditioner supports the following parameters:
  ///   - "fact: iluk level-of-fill" (int)
  ///   - "fact: absolute threshold" (magnitude_type)
  ///   - "fact: relative threshold" (magnitude_type)
  ///   - "fact: relax value" (magnitude_type)
  ///   - "fact: iluk overalloc" (double)
  void setParameters (const Teuchos::ParameterList& params);

  //! Initialize by computing the symbolic incomplete factorization.
  void initialize ();

  /// \brief Compute the (numeric) incomplete factorization.
  ///
  /// This function computes the RILU(k) factors L and U using the current:
  /// - Ifpack2_IlukGraph specifying the structure of L and U.
  /// - Value for the RILU(k) relaxation parameter.
  /// - Value for the a priori diagonal threshold values.
  ///
  /// initialize() must be called first, before this method may be called.
  void compute ();

  //! Whether initialize() has been called on this object.
  bool isInitialized () const {
    return isInitialized_;
  }
  //! Whether compute() has been called on this object.
  bool isComputed () const {
    return isComputed_;
  }

  //! Number of successful initialize() calls for this object.
  int getNumInitialize () const {
    return numInitialize_;
  }
  //! Number of successful compute() calls for this object.
  int getNumCompute () const {
    return numCompute_;
  }
  //! Number of successful apply() calls for this object.
  int getNumApply () const {
    return numApply_;
  }

  //! Total time in seconds taken by all successful initialize() calls for this object.
  double getInitializeTime () const {
    return initializeTime_;
  }
  //! Total time in seconds taken by all successful compute() calls for this object.
  double getComputeTime () const {
    return computeTime_;
  }
  //! Total time in seconds taken by all successful apply() calls for this object.
  double getApplyTime () const {
    return applyTime_;
  }

  //! Get a rough estimate of cost per iteration
  size_t getNodeSmootherComplexity() const;



  //! \name Implementation of Ifpack2::Details::CanChangeMatrix
  //@{

  /// \brief Change the matrix to be preconditioned.
  ///
  /// \param A [in] The new matrix.
  ///
  /// \post <tt>! isInitialized ()</tt>
  /// \post <tt>! isComputed ()</tt>
  ///
  /// Calling this method resets the preconditioner's state.  After
  /// calling this method with a nonnull input, you must first call
  /// initialize() and compute() (in that order) before you may call
  /// apply().
  ///
  /// You may call this method with a null input.  If A is null, then
  /// you may not call initialize() or compute() until you first call
  /// this method again with a nonnull input.  This method invalidates
  /// any previous factorization whether or not A is null, so calling
  /// setMatrix() with a null input is one way to clear the
  /// preconditioner's state (and free any memory that it may be
  /// using).
  ///
  /// The new matrix A need not necessarily have the same Maps or even
  /// the same communicator as the original matrix.
  virtual void
  setMatrix (const Teuchos::RCP<const row_matrix_type>& A);

  //@}
  //! @name Implementation of Teuchos::Describable interface
  //@{

  //! A one-line description of this object.
  std::string description () const;

  //@}
  //! \name Implementation of Tpetra::Operator
  //@{

  //! Returns the Tpetra::Map object associated with the domain of this operator.
  Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> >
  getDomainMap () const;

  //! Returns the Tpetra::Map object associated with the range of this operator.
  Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> >
  getRangeMap () const;

  /// \brief Apply the (inverse of the) incomplete factorization to X, resulting in Y.
  ///
  /// For an incomplete factorization \f$A \approx LDU\f$, this method
  /// computes the following, depending on the value of \c mode:
  /// <ul>
  /// <li> If mode = Teuchos::NO_TRANS, it computes
  ///      <tt>Y = beta*Y + alpha*(U \ (D \ (L \ X)))</tt> </li>
  /// <li> If mode = Teuchos::TRANS, it computes
  ///      <tt>Y = beta*Y + alpha*(L^T \ (D^T \ (U^T \ X)))</tt> </li>
  /// <li> If mode = Teuchos::CONJ_TRANS, it computes
  ///      <tt>Y = beta*Y + alpha*(L^* \ (D^* \ (U^* \ X)))</tt>,
  ///      where the asterisk indicates the conjugate transpose. </li>
  /// </ul>
  /// If alpha is zero, then the result of applying the operator to a
  /// vector is ignored.  This matters because zero times NaN (not a
  /// number) is NaN, not zero.  Analogously, if beta is zero, then
  /// any values in Y on input are ignored.
  ///
  /// \param X [in] The input multivector.
  ///
  /// \param Y [in/out] The output multivector.
  ///
  /// \param mode [in] If Teuchos::TRANS resp. Teuchos::CONJ_TRANS,
  ///   apply the transpose resp. conjugate transpose of the incomplete
  ///   factorization.  Otherwise, don't apply the tranpose.
  ///
  /// \param alpha [in] Scaling factor for the result of applying the preconditioner.
  ///
  /// \param beta [in] Scaling factor for the initial value of Y.
  void
  apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
         Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one (),
         scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero ()) const;
  //@}

private:
  /// \brief Apply the incomplete factorization (as a product) to X, resulting in Y.
  ///
  /// Given an incomplete factorization is \f$A \approx LDU\f$, this
  /// method computes the following, depending on the value of \c mode:
  ///
  ///   - If mode = Teuchos::NO_TRANS, it computes
  ///     <tt>Y = beta*Y + alpha*(L \ (D \ (U \ X)))</tt>
  ///   - If mode = Teuchos::TRANS, it computes
  ///     <tt>Y = beta*Y + alpha*(U^T \ (D^T \ (L^T \ X)))</tt>
  ///   - If mode = Teuchos::CONJ_TRANS, it computes
  ///     <tt>Y = beta*Y + alpha*(U^* \ (D^* \ (L^* \ X)))</tt>,
  ///     where the asterisk indicates the conjugate transpose.
  ///
  /// \param X [in] The input multivector.
  ///
  /// \param Y [in/out] The output multivector.
  ///
  /// \param mode [in] If Teuchos::TRANS resp. Teuchos::CONJ_TRANS,
  ///   apply the transpose resp. conjugate transpose of the
  ///   incomplete factorization.  Otherwise, don't apply the
  ///   transpose.
  void
  multiply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
            Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
            const Teuchos::ETransp mode = Teuchos::NO_TRANS) const;
public:

  //! Get the input matrix.
  Teuchos::RCP<const row_matrix_type> getMatrix () const;

  // Attribute access functions

  //! Get RILU(k) relaxation parameter
  magnitude_type getRelaxValue () const { return RelaxValue_; }

  //! Get absolute threshold value
  magnitude_type getAbsoluteThreshold () const { return Athresh_; }

  //! Get relative threshold value
  magnitude_type getRelativeThreshold () const {return Rthresh_;}

  //! Get level of fill (the "k" in ILU(k)).
  int getLevelOfFill () const { return LevelOfFill_; }

  //! Get overlap mode type
  Tpetra::CombineMode getOverlapMode () {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Ifpack2::RILUK::SetOverlapMode: "
      "RILUK no longer implements overlap on its own.  "
      "Use RILUK with AdditiveSchwarz if you want overlap.");
  }

  //! Returns the number of nonzero entries in the global graph.
  Tpetra::global_size_t getGlobalNumEntries () const {
    return getL ().getGlobalNumEntries () + getU ().getGlobalNumEntries ();
  }

  //! Return the Ifpack2::IlukGraph associated with this factored matrix.
  Teuchos::RCP<Ifpack2::IlukGraph<Tpetra::CrsGraph<local_ordinal_type,
                                                   global_ordinal_type,
                                                   node_type>, kk_handle_type> > getGraph () const {
    return Graph_;
  }

  //! Return the L factor of the ILU factorization.
  const crs_matrix_type& getL () const;

  //! Return the diagonal entries of the ILU factorization.
  const Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>&
  getD () const;

  //! Return the U factor of the ILU factorization.
  const crs_matrix_type& getU () const;

  //! Return the input matrix A as a Tpetra::CrsMatrix, if possible; else throws.
  Teuchos::RCP<const crs_matrix_type> getCrsMatrix () const;

  /// \brief Return A, wrapped in a LocalFilter, if necessary.
  ///
  /// "If necessary" means that if A is already a LocalFilter, or if
  /// its communicator only has one process, then we don't need to
  /// wrap it, so we just return A.
  static Teuchos::RCP<const row_matrix_type>
  makeLocalFilter (const Teuchos::RCP<const row_matrix_type>& A);

private:
  typedef Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> MV;
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;

  void allocateSolvers ();
  void allocate_L_and_U ();
  static void checkOrderingConsistency (const row_matrix_type& A);
  void initAllValues (const row_matrix_type& A);

  void compute_serial();
  void compute_kkspiluk();
// Workaround Cuda limitation of KOKKOS_LAMBDA in private/protected member functions
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_kkspiluk_stream();
#ifdef KOKKOS_ENABLE_CUDA
private:
#endif

protected:
  typedef Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> vec_type;

  //! The (original) input matrix for which to compute ILU(k).
  Teuchos::RCP<const row_matrix_type> A_;

  //! The ILU(k) graph.
  Teuchos::RCP<iluk_graph_type> Graph_;
  std::vector< Teuchos::RCP<iluk_graph_type> > Graph_v_;
  /// \brief The matrix whos numbers are used to to compute ILU(k). The graph
  /// may be computed using a crs_matrix_type that initialize() constructs
  /// temporarily.
  Teuchos::RCP<const row_matrix_type> A_local_;
  Teuchos::RCP<const crs_matrix_type> A_local_crs_;
  Teuchos::RCP<crs_matrix_type> A_local_crs_nc_;
  std::vector<local_matrix_device_type> A_local_diagblks;
  std::vector< lno_row_view_t > A_local_diagblks_rowmap_v_;
  std::vector< lno_nonzero_view_t > A_local_diagblks_entries_v_;
  std::vector< scalar_nonzero_view_t > A_local_diagblks_values_v_;

  //! The L (lower triangular) factor of ILU(k).
  Teuchos::RCP<crs_matrix_type> L_;
  std::vector< Teuchos::RCP<crs_matrix_type> > L_v_;
  //! Sparse triangular solver for L
  Teuchos::RCP<LocalSparseTriangularSolver<row_matrix_type> > L_solver_;
  //! The U (upper triangular) factor of ILU(k).
  Teuchos::RCP<crs_matrix_type> U_;
  std::vector< Teuchos::RCP<crs_matrix_type> > U_v_;
  //! Sparse triangular solver for U
  Teuchos::RCP<LocalSparseTriangularSolver<row_matrix_type> > U_solver_;

  //! The diagonal entries of the ILU(k) factorization.
  Teuchos::RCP<vec_type> D_;

  int LevelOfFill_;
  double Overalloc_;

  bool isAllocated_;
  bool isInitialized_;
  bool isComputed_;

  int numInitialize_;
  int numCompute_;
  mutable int numApply_;

  double initializeTime_;
  double computeTime_;
  mutable double applyTime_;

  magnitude_type RelaxValue_;
  magnitude_type Athresh_;
  magnitude_type Rthresh_;

  //! Optional KokkosKernels implementation.
  bool isKokkosKernelsSpiluk_;
  Teuchos::RCP<kk_handle_type> KernelHandle_;
  std::vector< Teuchos::RCP<kk_handle_type> > KernelHandle_v_;
  bool isKokkosKernelsStream_;
  int num_streams_;
  std::vector<execution_space> exec_space_instances_;
  bool hasStreamReordered_;
  std::vector<typename lno_nonzero_view_t::non_const_type> perm_v_;
  std::vector<typename lno_nonzero_view_t::non_const_type> reverse_perm_v_;
  mutable std::unique_ptr<MV> Y_tmp_;
  mutable std::unique_ptr<MV> reordered_x_;
  mutable std::unique_ptr<MV> reordered_y_;
};

// NOTE (mfh 11 Feb 2015) This used to exist in order to deal with
// different behavior of Tpetra::Crs{Graph,Matrix} for
// KokkosClassic::ThrustGPUNode.  In particular, fillComplete on a
// CrsMatrix used to make the graph go away by default, so we had to
// pass in a parameter to keep a host copy of the graph.  With the new
// (Kokkos refactor) version of Tpetra, this problem has gone away.
namespace detail {
  template<class MatrixType, class NodeType>
  struct setLocalSolveParams{
    static Teuchos::RCP<Teuchos::ParameterList>
    setParams (const Teuchos::RCP<Teuchos::ParameterList>& param) {
      return param;
    }
  };
} // namespace detail

} // namespace Ifpack2

#endif /* IFPACK2_CRSRILUK_DECL_HPP */
