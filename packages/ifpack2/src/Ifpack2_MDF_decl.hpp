// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Ifpack2_RILUK_decl.hpp
/// \brief Declaration of MDF interface

#ifndef IFPACK2_MDF_DECL_HPP
#define IFPACK2_MDF_DECL_HPP

#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Details_CanChangeMatrix.hpp"
#include "Tpetra_CrsMatrix_decl.hpp"
#include "Ifpack2_ScalingType.hpp"
#include "Ifpack2_IlukGraph.hpp"
#include "Ifpack2_LocalSparseTriangularSolver_decl.hpp"
#include "KokkosSparse_mdf.hpp"

#include <type_traits>

namespace Teuchos {
  class ParameterList; // forward declaration
}
namespace Ifpack2 {

/// \class MDF
/// \brief MDF (incomplete LU factorization with minimum discarded fill reordering) of a
///   Tpetra sparse matrix
/// \tparam A specialization of Tpetra::RowMatrix.
///
/// This class computes a sparse MDF (incomplete LU) factorization
/// with a reordering that minimizes the discarded fill of the local part of a
/// given sparse matrix represented as a Tpetra::RowMatrix or
/// Tpetra::CrsMatrix.  The "local part" is the square diagonal block
/// of the matrix owned by the calling process.  Thus, if the input
/// matrix is distributed over multiple MPI processes, this
/// preconditioner is equivalent to nonoverlapping additive Schwarz
/// domain decomposition over the MPI processes, with MDF as the
/// subdomain solver on each process.
///
/// @remark See the documentation of setParameters() for a list of valid
/// parameters.
///
template<class MatrixType>
class MDF:
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


  static_assert(std::is_same<MatrixType, row_matrix_type>::value, "Ifpack2::MDF: The template parameter MatrixType must be a Tpetra::RowMatrix specialization.  Please don't use Tpetra::CrsMatrix (a subclass of Tpetra::RowMatrix) here anymore.");

  //! Tpetra::CrsMatrix specialization used by this class for representing L and U.
  typedef Tpetra::CrsMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> crs_matrix_type;

  //! Scalar type stored in Kokkos::Views (CrsMatrix and MultiVector)
  typedef typename crs_matrix_type::impl_scalar_type impl_scalar_type;

  template <class NewMatrixType> friend class MDF;

  typedef typename crs_matrix_type::global_inds_host_view_type global_inds_host_view_type;
  typedef typename crs_matrix_type::local_inds_host_view_type local_inds_host_view_type;
  typedef typename crs_matrix_type::values_host_view_type values_host_view_type;


  typedef typename crs_matrix_type::nonconst_global_inds_host_view_type nonconst_global_inds_host_view_type;
  typedef typename crs_matrix_type::nonconst_local_inds_host_view_type nonconst_local_inds_host_view_type;
  typedef typename crs_matrix_type::nonconst_values_host_view_type nonconst_values_host_view_type;


  //@}
  //! \name Implementation of Kokkos Kernels MDF.
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
  
  /// \brief Constructor that takes a Tpetra::RowMatrix.
  ///
  /// \param A_in [in] The input matrix.
  MDF (const Teuchos::RCP<const row_matrix_type>& A_in);

  /// \brief Constructor that takes a Tpetra::CrsMatrix.
  ///
  /// This constructor exists to avoid "ambiguous constructor"
  /// warnings.  It does the same thing as the constructor that takes
  /// a Tpetra::RowMatrix.
  ///
  /// \param A_in [in] The input matrix.
  MDF (const Teuchos::RCP<const crs_matrix_type>& A_in);

 private:
  /// \brief Copy constructor: declared private but not defined, so
  ///   that calling it is syntactically forbidden.
  MDF (const MDF<MatrixType> & src);

 public:
  //! Destructor (declared virtual for memory safety).
  virtual ~MDF () = default;

  /// Set parameters for the incomplete factorization.
  ///
  /// This preconditioner supports the following parameters:
  ///   - "fact: mdf level-of-fill" (int)
  ///   - "fact: relax value" (magnitude_type)
  ///   - "fact: mdf overalloc" (double)
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

  // Split off to a different impl call so that nested apply calls don't mess up apply counts/timers
  void apply_impl (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
         Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one (),
         scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero ()) const;

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
  using MDF_handle_device_type = KokkosSparse::Experimental::MDF_handle<local_matrix_device_type>;
  using permutations_type = Teuchos::ArrayRCP<local_ordinal_type>;

  //! Get the input matrix.
  Teuchos::RCP<const row_matrix_type> getMatrix () const;

  //! Get level of fill (the "k" in ILU(k)).
  int getLevelOfFill () const { return LevelOfFill_; }

  //! Get overlap mode type
  Tpetra::CombineMode getOverlapMode () {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Ifpack2::MDF::SetOverlapMode: "
      "MDF no longer implements overlap on its own.  "
      "Use MDF with AdditiveSchwarz if you want overlap.");
  }

  //! Returns the number of nonzero entries in the global graph.
  Tpetra::global_size_t getGlobalNumEntries () const {
    return getL ().getGlobalNumEntries () + getU ().getGlobalNumEntries ();
  }

  //! Return the L factor of the MDF factorization.
  const crs_matrix_type& getL () const;

  //! Return the U factor of the MDF factorization.
  const crs_matrix_type& getU () const;

  //! Return the permutations of the MDF factorization
  permutations_type & getPermutations() const;

  //! Return the reverse permutations of the MDF factorization
  permutations_type & getReversePermutations() const;

  //! Return the input matrix A as a Tpetra::CrsMatrix, if possible; else throws.
  Teuchos::RCP<const crs_matrix_type> getCrsMatrix () const;

private:
  typedef Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> MV;
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;

  void allocateSolvers ();
  void allocatePermutations(bool force = false);
  static void checkOrderingConsistency (const row_matrix_type& A);
  // void initAllValues (const row_matrix_type& A);

  /// \brief Return A, wrapped in a LocalFilter, if necessary.
  ///
  /// "If necessary" means that if A is already a LocalFilter, or if
  /// its communicator only has one process, then we don't need to
  /// wrap it, so we just return A.
  static Teuchos::RCP<const row_matrix_type>
  makeLocalFilter (const Teuchos::RCP<const row_matrix_type>& A);

protected:
  typedef Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> vec_type;

  //! The (original) input matrix for which to compute ILU(k).
  Teuchos::RCP<const row_matrix_type> A_;

  // The MDF handle
  Teuchos::RCP<MDF_handle_device_type> MDF_handle_;

  /// \brief The matrix whos numbers are used to to compute ILU(k). The graph
  /// may be computed using a crs_matrix_type that initialize() constructs
  /// temporarily.
  Teuchos::RCP<const row_matrix_type> A_local_;
  lno_row_view_t A_local_rowmap_; 
  lno_nonzero_view_t A_local_entries_; 
  scalar_nonzero_view_t A_local_values_;

  //! The L (lower triangular) factor of ILU(k).
  Teuchos::RCP<crs_matrix_type> L_;
  //! Sparse triangular solver for L
  Teuchos::RCP<LocalSparseTriangularSolver<row_matrix_type> > L_solver_;
  //! The U (upper triangular) factor of ILU(k).
  Teuchos::RCP<crs_matrix_type> U_;
  //! Sparse triangular solver for U
  Teuchos::RCP<LocalSparseTriangularSolver<row_matrix_type> > U_solver_;

  //! The computed permuations from MDF factorization
  permutations_type permutations_;

  //! The reverse permuations from MDF factorization
  permutations_type reversePermutations_;

  int Verbosity_;

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
};

} // namespace Ifpack2

#endif /* IFPACK2_MDF_DECL_HPP */
