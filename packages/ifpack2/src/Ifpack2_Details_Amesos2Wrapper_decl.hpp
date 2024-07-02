// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Ifpack2_Details_Amesos2Wrapper_decl.hpp
/// \brief Declaration of wrapper for Amesos2 solvers
///
/// \warning This header file and its contents is an implementation
///   detail of Ifpack2.  Users of Ifpack2 must NOT depend on this
///   header file continuing to exist, on its name, or on any of its
///   contents.

#ifndef IFPACK2_DETAILS_AMESOS2WRAPPER_DECL_HPP
#define IFPACK2_DETAILS_AMESOS2WRAPPER_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"

#ifdef HAVE_IFPACK2_AMESOS2

#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Details_CanChangeMatrix.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include <type_traits>

namespace Teuchos {
  // forward declaration
  class ParameterList;
}

namespace Trilinos {
namespace Details {
  template<class MV, class OP, class NormType>
  class LinearSolver; // forward declaration
} // namespace Details
} // namespace Trilinos


namespace Ifpack2 {
namespace Details {
/// \class Amesos2Wrapper
/// \brief Wrapper class for direct solvers in Amesos2.
/// \tparam MatrixType A specialization of Tpetra::RowMatrix.
///
/// This class computes a sparse factorization of the input
/// matrix A using Amesos2.  The apply() method solves linear
/// system(s) using that factorization.  As with all Ifpack2
/// preconditioners, initialize() computes the symbolic factorization,
/// and compute() computes the numeric factorization.
///
/// \warning This class is an implementation detail of Ifpack2.  Users
///   must not rely on this class.  It may go away or its interface
///   may change at any time.
///
/// \warning This class creates a local filter.  In particular, if the
///   matrix is not a Tpetra::CrsMatrix instance (this class will
///   check this at run time using a dynamic cast), then this class
///   will perform a deep copy to produce a CrsMatrix.  This will
///   happen, for example, if you are doing additive Schwarz with
///   nonzero overlap, and apply Amesos2 as the subdomain solve.  This
///   deep copy is required by Amesos2, and is in addition to any data
///   copying that Amesos2 may do internally to satisfy TPL storage
///   formats.
template<class MatrixType>
class Amesos2Wrapper :
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

  static_assert(std::is_same<MatrixType, row_matrix_type>::value,
                "Ifpack2::Details::Amesos2Wrapper: Please use MatrixType = Tpetra::RowMatrix.");

  //! Type of the Tpetra::Map specialization that this class uses.
  typedef Tpetra::Map<local_ordinal_type,
                      global_ordinal_type,
                      node_type> map_type;

  //! Type of the Tpetra::CrsMatrix specialization that this class uses.
  typedef Tpetra::CrsMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> crs_matrix_type;
  //@}
  //! \name Constructors and destructor
  //@{

  /// \brief Constructor
  ///
  /// \param A [in] The sparse matrix to factor, as a
  ///   Tpetra::RowMatrix.  (Tpetra::CrsMatrix inherits from this, so
  ///   you may use a Tpetra::CrsMatrix here instead.)
  explicit Amesos2Wrapper (const Teuchos::RCP<const row_matrix_type>& A);

  //! Destructor
  virtual ~Amesos2Wrapper();

  //@}
  //! \name Methods for setting up and computing the factorization
  //@{

  /// \brief Set parameters.
  void setParameters (const Teuchos::ParameterList& params);

  /// \brief Compute the preordering and symbolic factorization of the matrix.
  ///
  /// You must call this method before you may call compute() or
  /// apply(), under the following conditions:
  /// <ol>
  /// <li> If you have not yet called initialize() before on this instance </li>
  /// <li> If you have just called setMatrix() with a nonnull matrix </li>
  /// <li> If the structure of the sparse matrix has changed </li>
  /// </ol>
  /// Please also see the documentation of compute() for the
  /// conditions under which you must call compute() before calling
  /// apply().
  void initialize ();

  //! Returns \c true if the preconditioner has been successfully initialized.
  inline bool isInitialized() const {
    return IsInitialized_;
  }

  /// \brief Compute the numeric factorization of the matrix.
  ///
  /// You must call this method before you may call apply(), under the
  /// following conditions:
  /// <ol>
  /// <li> If you have not yet called compute() before on this instance </li>
  /// <li> If the values in the sparse matrix has changed </li>
  /// </ol>
  /// Please also see the documentation of initialize() for the
  /// conditions under which you must call initialize() before calling
  /// compute() or apply().
  void compute();

  //! True if compute() completed successfully, else false.
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

  /// \brief Apply the preconditioner to X, resulting in Y.
  ///
  /// If \f$M^{-1}\f$ represents the action of the preconditioner,
  /// then this routine computes \f$Y := beta \cdot Y + alpha \cdot
  /// Op(M^{-1}) X\f$, where \f$Op(M^{-1})\f$ can be either
  /// \f$M^{-1}\f$, \f$M^{-T}\f$, or \f$M^{-H}\f$, depending on \c
  /// mode.
  ///
  /// \param X [in] Input multivector; "right-hand side" of the solve.
  /// \param Y [out] Output multivector; result of the solve.
  /// \param mode [in] Whether to solve with the original matrix, its
  ///   transpose, or its conjugate transpose.
  /// \param alpha [in] Scaling factor for the result of applying the
  ///   preconditioner.
  /// \param alpha [in] Scaling factor for Y.
  void
  apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
         Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
         scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //! Tpetra::Map representing the domain of this operator.
  Teuchos::RCP<const map_type> getDomainMap () const;

  //! Tpetra::Map representing the range of this operator.
  Teuchos::RCP<const map_type> getRangeMap () const;

  //! Whether this object's apply() method can apply the transpose (or conjugate transpose, if applicable).
  bool hasTransposeApply () const;

  //@}
  //! \name Mathematical functions
  //@{

  //! The input matrix's communicator.
  Teuchos::RCP<const Teuchos::Comm<int> > getComm () const;

  //! The input matrix; the matrix to be preconditioned.
  Teuchos::RCP<const row_matrix_type> getMatrix () const;

  //! The total number of successful calls to initialize().
  int getNumInitialize () const;

  //! The total number of successful calls to compute().
  int getNumCompute () const;

  //! The total number of successful calls to apply().
  int getNumApply () const;

  //! The total time in seconds spent in successful calls to initialize().
  double getInitializeTime () const;

  //! The total time in seconds spent in successful calls to compute().
  double getComputeTime () const;

  //! The total time in seconds spent in successful calls to apply().
  double getApplyTime () const;

  //@}
  //! \name Implementation of Teuchos::Describable
  //@{

  //! A one-line description of this object.
  std::string description () const;

  //! Print the object with some verbosity level to the given output stream.
  void
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel =
            Teuchos::Describable::verbLevel_default) const;
  //@}

private:
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;
  typedef typename Teuchos::Array<local_ordinal_type>::size_type size_type;

  //! Type of the Tpetra::MultiVector specialization that this class uses.
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> MV;
  //! Type of the Tpetra::Operator specialization that this class uses.
  typedef Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type, node_type> OP;

  //! Copy constructor (declared private and undefined; may not be used)
  Amesos2Wrapper (const Amesos2Wrapper<MatrixType>& RHS);

  //! operator= (declared private and undefined; may not be used)
  Amesos2Wrapper<MatrixType>& operator= (const Amesos2Wrapper<MatrixType>& RHS);

  //! Amesos2 solver; it contains the factorization of the matrix A_.
  Teuchos::RCP<Trilinos::Details::LinearSolver<MV, OP, typename MV::mag_type> > solver_;

  /// \brief Return A, wrapped in a LocalFilter, if necessary.
  ///
  /// "If necessary" means that if A is already a LocalFilter, or if
  /// its communicator only has one process, then we don't need to
  /// wrap it, so we just return A.
  static Teuchos::RCP<const row_matrix_type>
  makeLocalFilter (const Teuchos::RCP<const row_matrix_type>& A);

  //! The (original) input matrix to be preconditioned.
  //Teuchos::RCP<const MatrixType> A_;
  Teuchos::RCP<const row_matrix_type> A_;

  /// \brief The matrix used to compute the Amesos2 preconditioner.
  ///
  /// If A_local (the local filter of the original input matrix) is a
  /// Tpetra::CrsMatrix, then this is just A_local.  Otherwise, this
  /// class reserves the right for A_local_crs_ to be a copy of
  /// A_local.  This is because the current adapters in Amesos2
  /// only accept a Tpetra::CrsMatrix.  That may change
  /// in the future.
  Teuchos::RCP<const crs_matrix_type> A_local_crs_;

  //@}
  // \name Parameters (set by setParameters())
  //@{

  //! Caches parameters passed into setParameters
  //! if the concrete inner preconditioner doesn't exist yet.
  Teuchos::RCP<const Teuchos::ParameterList> parameterList_;

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
  //! @brief The name of the solver in Amesos2 to use.
  //! See the Amesos2 documentation for valid names.
  std::string SolverName_;
  //@}
}; // class Amesos2Wrapper

} // namespace Details
} // namespace Ifpack2

#endif // HAVE_IFPACK2_AMESOS2

#endif // IFPACK2_DETAILS_AMESOS2WRAPPER_DECL_HPP
