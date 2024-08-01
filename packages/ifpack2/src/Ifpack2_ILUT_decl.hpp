// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Ifpack2_ILUT_decl.hpp
/// \brief Declaration of ILUT preconditioner

#ifndef IFPACK2_ILUT_DECL_HPP
#define IFPACK2_ILUT_DECL_HPP

#include "KokkosSparse_par_ilut.hpp"

#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Details_CanChangeMatrix.hpp"
#include "Tpetra_CrsMatrix_decl.hpp"
#include "Ifpack2_LocalSparseTriangularSolver_decl.hpp"

#include <string>
#include <sstream>
#include <iostream>
#include <cmath>
#include <type_traits>

namespace Teuchos {
  class ParameterList; // forward declaration
}

namespace Ifpack2 {

/// \class ILUT
/// \brief ILUT (incomplete LU factorization with threshold) of a
///   Tpetra sparse matrix
/// \tparam A specialization of Tpetra::RowMatrix.
///
/// This class computes a sparse ILUT (incomplete LU) factorization
/// with specified fill and drop tolerance, of the local part of a
/// given sparse matrix represented as a Tpetra::RowMatrix or
/// Tpetra::CrsMatrix.  The "local part" is the square diagonal block
/// of the matrix owned by the calling process.  Thus, if the input
/// matrix is distributed over multiple MPI processes, this
/// preconditioner is equivalent to nonoverlapping additive Schwarz
/// domain decomposition over the MPI processes, with ILUT as the
/// subdomain solver on each process.
///
/// @remark See the documentation of setParameters() for a list of valid
/// parameters.
///
/// @remark This version of ILUT is a translation of Aztec's ILUT
/// implementation, which was written by Ray Tuminaro.
///
/// @remark There is an important difference between this implementation and the version
/// described in Saad's paper.  See setParameters() for details.
///
template<class MatrixType>
class ILUT :
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
  //! \name Typedefs
  //@{

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

  //! Type of the Tpetra::RowMatrix specialization that this class uses.
  typedef Tpetra::RowMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> row_matrix_type;

  typedef typename row_matrix_type::global_inds_host_view_type global_inds_host_view_type;
  typedef typename row_matrix_type::local_inds_host_view_type local_inds_host_view_type;
  typedef typename row_matrix_type::values_host_view_type values_host_view_type;

  typedef typename row_matrix_type::nonconst_global_inds_host_view_type nonconst_global_inds_host_view_type;
  typedef typename row_matrix_type::nonconst_local_inds_host_view_type nonconst_local_inds_host_view_type;
  typedef typename row_matrix_type::nonconst_values_host_view_type nonconst_values_host_view_type;

  static_assert(std::is_same<MatrixType, row_matrix_type>::value, "Ifpack2::ILUT: The template parameter MatrixType must be a Tpetra::RowMatrix specialization.  Please don't use Tpetra::CrsMatrix (a subclass of Tpetra::RowMatrix) here anymore.  The constructor can take either a RowMatrix or a CrsMatrix just fine.");

  //! Type of the Tpetra::CrsMatrix specialization that this class uses for the L and U factors.
  typedef Tpetra::CrsMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> crs_matrix_type;

  //! \name For implementation of Kokkos Kernels parallel ILUt (thresholded ILU)
  typedef typename crs_matrix_type::local_matrix_device_type local_matrix_device_type;
  typedef typename local_matrix_device_type::StaticCrsGraphType::row_map_type lno_row_view_t;
  ////////////////////////////////////
  typedef Tpetra::CrsGraph<local_ordinal_type, global_ordinal_type, node_type> crs_graph_type;
  typedef typename crs_graph_type::local_graph_device_type local_graph_device_type;
  typedef typename local_graph_device_type::array_layout   array_layout;
  typedef typename local_graph_device_type::device_type    device_type;
  typedef typename local_graph_device_type::size_type      usize_type;
  //KokkosKernels requires unsigned
  //typedef typename Kokkos::View<size_type*, array_layout, device_type> lno_row_view_t;
  typedef typename Kokkos::View<usize_type*, array_layout, device_type> lno_urow_view_t;
  ////////////////////////////////////
  typedef typename local_matrix_device_type::StaticCrsGraphType::entries_type lno_nonzero_view_t;
  //typedef typename Kokkos::View<lno_nonzero_view_t*, array_layout, device_type> static_graph_entries_t;
  typedef typename Kokkos::View<typename local_matrix_device_type::non_const_ordinal_type*, array_layout, device_type> static_graph_entries_t;
  typedef typename local_matrix_device_type::values_type scalar_nonzero_view_t;
  //typedef typename Kokkos::View<scalar_nonzero_view_t*, array_layout, device_type> local_matrix_values_t;
  typedef typename Kokkos::View<typename local_matrix_device_type::non_const_value_type*, array_layout, device_type> local_matrix_values_t;
  typedef typename local_matrix_device_type::StaticCrsGraphType::device_type::memory_space TemporaryMemorySpace;
  typedef typename local_matrix_device_type::StaticCrsGraphType::device_type::memory_space PersistentMemorySpace;
  typedef typename local_matrix_device_type::StaticCrsGraphType::device_type::execution_space HandleExecSpace;
  typedef typename KokkosKernels::Experimental::KokkosKernelsHandle
    <typename lno_row_view_t::const_value_type, typename lno_nonzero_view_t::const_value_type, typename scalar_nonzero_view_t::value_type,
    HandleExecSpace, TemporaryMemorySpace,PersistentMemorySpace > kk_handle_type;

  //! Type of the Tpetra::Map specialization that this class uses.
  typedef Tpetra::Map<local_ordinal_type,
                      global_ordinal_type,
                      node_type> map_type;
  //@}
  //! \name Constructors and Destructors
  //@{

  /// \brief Constructor
  ///
  /// \param A [in] The sparse matrix to factor, as a
  ///   Tpetra::RowMatrix.  (Tpetra::CrsMatrix inherits from this, so
  ///   you may use a Tpetra::CrsMatrix here instead.)
  ///
  /// The factorization will <i>not</i> modify the input matrix.  It
  /// stores the L and U factors in the incomplete factorization
  /// separately.
  explicit ILUT (const Teuchos::RCP<const row_matrix_type>& A);

  //! Destructor
  virtual ~ILUT () = default;

  //@}
  //! \name Methods for setting up and computing the incomplete factorization
  //@{

  /// \brief Set preconditioner parameters.
  ///
  /// ILUT implements the following parameters:
  /// <ul>
  /// <li> "fact: ilut level-of-fill" (\c double)
  /// <li> "fact: drop tolerance" (\c magnitude_type)
  /// <li> "fact: absolute threshold" (\c magnitude_type)
  /// <li> "fact: relative threshold" (\c magnitude_type)
  /// <li> "fact: relax value" (\c magnitude_type)
  /// </ul>
  /// "fact: drop tolerance" is the magnitude threshold for dropping
  /// entries.  It corresponds to the \f$\tau\f$ parameter in Saad's
  /// original description of ILUT.  "fact: ilut level-of-fill" controls the
  /// number of entries to keep in the strict upper triangle of the
  /// current row, and in the strict lower triangle of the current
  /// row.  It does <B>not</B> correspond to the \f$p\f$ parameter in Saad's original
  /// description. This parameter represents a maximum fill fraction.
  /// In this implementation, the L and U factors always contains nonzeros corresponding
  /// to the original sparsity pattern of A, so this value should be >= 1.0.
  /// Letting \f$fill = \frac{(level-of-fill - 1)*nnz(A)}{2*N}\f$,
  /// each row of the computed L and U factors contains at most \f$fill\f$
  /// nonzero elements in addition to those from the sparsity pattern of A.
  /// ILUT always keeps the diagonal entry in the
  /// current row, regardless of the drop tolerance or fill level.
  ///
  /// The absolute and relative threshold parameters affect how this
  /// code modifies the diagonal entry of the output factor.  These
  /// parameters are not part of the original ILUT algorithm, but we
  /// include them for consistency with other Ifpack2 preconditioners.
  ///
  /// The "fact: relax value" parameter currently has no effect.
  void setParameters (const Teuchos::ParameterList& params);

  /// \brief Clear any previously computed factors, and potentially
  ///        compute sparsity patterns of factors.
  ///
  /// You may call this before calling compute().  The compute()
  /// method will call this automatically if it has not yet been
  /// called.  If you call this after calling compute(), you must
  /// recompute the factorization (by calling compute() again) before
  /// you may call apply().
  ///
  /// If your are using Par_ILUT from Kokkos Kernels, initialize()
  /// will also perform a symbolic factorization (i.e., compute sparsity patterns of factors).
  void initialize ();

  //! Returns \c true if the preconditioner has been successfully initialized.
  inline bool isInitialized() const {
    return IsInitialized_;
  }

  //! Compute factors L and U using the specified diagonal perturbation thresholds and relaxation parameters.
  /*! This function computes the ILUT factors L and U using the current:
    <ol>
    <li> Value for the drop tolerance
    <li> Value for the level of fill
    <li> Value for the \e a \e priori diagonal threshold values.
    </ol>
   */
  void compute();

  //! If compute() is completed, this query returns true, otherwise it returns false.
  inline bool isComputed() const {
    return IsComputed_;
  }

  //@}
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
  //! \name Implementation of Tpetra::Operator
  //@{

  /// \brief Apply the ILUT preconditioner to X, resulting in Y.
  ///
  /// \param X [in] Input multivector; "right-hand side" of the solve.
  /// \param Y [out] Output multivector; result of the solve.
  void
  apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
         Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
         scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //! Tpetra::Map representing the domain of this operator.
  Teuchos::RCP<const map_type> getDomainMap() const;

  //! Tpetra::Map representing the range of this operator.
  Teuchos::RCP<const map_type> getRangeMap() const;

  //! Whether this object's apply() method can apply the transpose (or conjugate transpose, if applicable).
  bool hasTransposeApply() const;

  //@}
  //! \name Mathematical functions
  //@{

  //! Returns the input matrix's communicator.
  Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

  //! Returns a reference to the matrix to be preconditioned.
  Teuchos::RCP<const row_matrix_type> getMatrix () const;

  //! Returns a reference to the L factor.
  Teuchos::RCP<const crs_matrix_type> getL () const { return L_; }

  //! Returns a reference to the U factor.
  Teuchos::RCP<const crs_matrix_type> getU () const { return U_; }

  //! Returns the number of calls to Initialize().
  int getNumInitialize() const;

  //! Returns the number of calls to Compute().
  int getNumCompute() const;

  //! Returns the number of calls to apply().
  int getNumApply() const;

  //! Returns the time spent in Initialize().
  double getInitializeTime() const;

  //! Returns the time spent in Compute().
  double getComputeTime() const;

  //! Returns the time spent in apply().
  double getApplyTime() const;

  //! Get a rough estimate of cost per iteration
  size_t getNodeSmootherComplexity() const;


  /// \brief The level of fill.
  ///
  /// For ILUT, this means the maximum number of entries in each row
  /// of the resulting L and U factors (each considered separately),
  /// not including the diagonal entry in that row (which is always
  /// part of U).  This has a different meaning for ILUT than it does
  /// for ILU(k).
  inline double getLevelOfFill() const {
    return LevelOfFill_;
  }

  //! Get absolute threshold value
  inline magnitude_type getAbsoluteThreshold() const {
    return(Athresh_);
  }

  //! Get relative threshold value
  inline magnitude_type getRelativeThreshold() const {
    return(Rthresh_);
  }

  //! Get the relax value
  inline magnitude_type getRelaxValue() const {
    return(RelaxValue_);
  }

  //! Gets the dropping tolerance
  inline magnitude_type getDropTolerance() const {
    return(DropTolerance_);
  }

  //! Returns the number of nonzero entries in the global graph.
  global_size_t getGlobalNumEntries() const;

  //! Returns the number of nonzero entries in the local graph.
  size_t getLocalNumEntries() const;

  //@}
  //! \name Implementation of Teuchos::Describable
  //@{

  /** \brief Return a simple one-line description of this object. */
  std::string description() const;

  /** \brief Print the object with some verbosity level to an FancyOStream object. */
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

  //@}

private:
  typedef Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> MV;
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;
  typedef typename Teuchos::Array<local_ordinal_type>::size_type size_type;

  //! Copy constructor (declared private and undefined; may not be used)
  ILUT (const ILUT<MatrixType>& RHS);

  void allocateSolvers ();

  //! operator= (declared private and undefined; may not be used)
  ILUT<MatrixType>& operator= (const ILUT<MatrixType>& RHS);

  /// \brief Wrap the given matrix in a "local filter," if necessary.
  ///
  /// A "local filter" excludes rows and columns that do not belong to
  /// the calling process.  It also uses a "serial" communicator
  /// (equivalent to MPI_COMM_SELF) rather than the matrix's original
  /// communicator.
  ///
  /// If the matrix's communicator only contains one process, then the
  /// matrix is already "local," so this function just returns the
  /// original input.
  static Teuchos::RCP<const row_matrix_type>
  makeLocalFilter (const Teuchos::RCP<const row_matrix_type>& A);

  // \name The matrix and its incomplete LU factors
  //@{

  //! The matrix to be preconditioned.
  Teuchos::RCP<const row_matrix_type> A_;
  //! "Local filter" version of A_.
  Teuchos::RCP<const row_matrix_type> A_local_;
  lno_row_view_t A_local_rowmap_;
  lno_nonzero_view_t A_local_entries_;
  scalar_nonzero_view_t A_local_values_;
  //! L factor of the incomplete LU factorization of A_local_.
  Teuchos::RCP<crs_matrix_type> L_;
  lno_urow_view_t     L_rowmap_;
  lno_urow_view_t     U_rowmap_;
  lno_urow_view_t     L_rowmap_orig_;
  lno_urow_view_t     U_rowmap_orig_;
  //! Sparse triangular solver for L
  Teuchos::RCP<LocalSparseTriangularSolver<row_matrix_type> > L_solver_;
  //! U factor of the incomplete LU factorization of A_local_.
  Teuchos::RCP<crs_matrix_type> U_;
  //! Sparse triangular solver for U
  Teuchos::RCP<LocalSparseTriangularSolver<row_matrix_type> > U_solver_;

  //@}
  // \name Parameters (set by setParameters())
  //@{

  magnitude_type Athresh_; //!< Absolute threshold
  magnitude_type Rthresh_; //!< Relative threshold
  magnitude_type RelaxValue_; //!< Relax value
  double LevelOfFill_; //!< Max fill level
  //! Discard all elements below this tolerance
  magnitude_type DropTolerance_;
  // See https://kokkos-kernels.readthedocs.io/en/latest/developer/apidocs/sparse.html#par-ilut
  // for more information on the following options.
  mutable struct par_ilut_option_struct {
    int max_iter;
    magnitude_type residual_norm_delta_stop;
    int team_size;
    int vector_size;
    double fill_in_limit; //Note: par_ilut declares this as float
    bool verbose;
  } par_ilut_options_;

  //@}
  // \name Other internal data
  //@{

  //! Total time in seconds for all successful calls to initialize().
  double InitializeTime_;
  //! Total time in seconds for all successful calls to compute().
  double ComputeTime_;
  //! Total time in seconds for all successful calls to apply().
  mutable double ApplyTime_;
  //! The number of successful calls to initialize().
  int NumInitialize_;
  //! The number of successful call to compute().
  int NumCompute_;
  //! The number of successful call to apply().
  mutable int NumApply_;
  //! \c true if \c this object has been initialized
  bool IsInitialized_;
  //! \c true if \c this object has been computed
  bool IsComputed_;
  //@}

  //! Optional KokkosKernels implementation.
  bool useKokkosKernelsParILUT_;
  Teuchos::RCP<kk_handle_type> KernelHandle_;

}; // class ILUT

} // namespace Ifpack2

#endif /* IFPACK2_ILUT_HPP */
