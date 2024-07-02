// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Ifpack2_SupportGraph_decl.hpp
/// \brief Declaration of SupportGraph preconditioner

#ifndef IFPACK2_SUPPORTGRAPH_DECL_HPP
#define IFPACK2_SUPPORTGRAPH_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Details_CanChangeMatrix.hpp"

// Ifpack2's CMake system should (and does) prevent Trilinos from
// attempting to build or install this class, if Amesos2 is not
// enabled.  We check for this case regardless, in order to catch any
// bugs that future development might introduce in the CMake scripts.

#ifdef HAVE_IFPACK2_AMESOS2
#  include <Amesos2.hpp>
#  include <Amesos2_Version.hpp>
#else
#  error "Ifpack2::SupportGraph requires that Trilinos be built with the Amesos2 package enabled."
#endif // HAVE_IFPACK2_AMESOS2


namespace Ifpack2 {

/// \class SupportGraph
/// \brief SupportGraph of a Tpetra sparse matrix.
/// \tparam Specialization of Tpetra::RowMatrix.
///
/// This class computes a maximum weight spanning tree
/// or multiple trees (forest), of a given sparse matrix
/// represented as a Tpetra::RowMatrix.
///
/// \warning Do not attempt to use this class unless Trilinos was
///   built with support for the Lemon and Cholmod third-party
///   libraries.
///
/// \warning This class will not be installed unless the CMake option
///   Ifpack2_ENABLE_Experimental was set to ON when configuring
///   Trilinos.
///
/// \warning If the matrix is distributed over multiple MPI processes,
///   this class will not work correctly by itself.  You must use it
///   as a subdomain solver inside of a domain decomposition method
///   like AdditiveSchwarz (which see).  If you use Factory to create
///   an SupportGraph preconditioner, the Factory will automatically
///   wrap SupportGraph in AdditiveSchwarz for you, if the matrix's
///   communicator contains multiple processes.
///
/// See the documentation of setParameters() for a list of valid
/// parameters.
template<class MatrixType>
class SupportGraph :
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

  //! Type of the Tpetra::Map specialization that this class uses.
  typedef Tpetra::Map<local_ordinal_type,
                      global_ordinal_type,
                      node_type> map_type;

  //@}
  //! \name Constructors and destructor
  //@{

  /// \brief Constructor
  ///
  /// \param A [in] The sparse matrix to factor, as a
  ///   Tpetra::RowMatrix.  (Tpetra::CrsMatrix inherits from this, so
  ///   you may use a Tpetra::CrsMatrix here instead.)
  ///
  /// This will <i>not</i> modify the input matrix.  It
  /// stores the support graph separately.
  explicit SupportGraph (const Teuchos::RCP<const row_matrix_type>& A);

  //! Destructor
  virtual ~SupportGraph();

  //@}
  //! \name Methods for setting up and computing support graph
  //@{

  /// \brief Set preconditioner parameters.
  ///
  /// SupportGraph implements the following parameters:
  /// <ul>
  /// <li> "fact: absolute threshold" (\c magnitude_type)
  /// <li> "fact: relative threshold" (\c magnitude_type)
  /// </ul>
  /// The absolute and relative threshold parameters affect how this
  /// code modifies the diagonal entry of the output factor.
  ///
  void setParameters (const Teuchos::ParameterList& params);

  /// \brief Clear any previously computed support graph.
  ///
  /// You may call this before calling compute().  The compute()
  /// method will call this automatically if it has not yet been
  /// called.  If you call this after calling compute(), you must
  /// recompute the factorization (by calling compute() again) before
  /// you may call apply().
  void initialize ();

  //! Returns \c true if the preconditioner has been successfully initialized.
  inline bool isInitialized() const {
    return IsInitialized_;
  }

  //! Compute factor by having Cholmod do a complete factorization.
  /*! This function computes the SupportGraph factor using the current:
    <ol>
    <li> Value for the \e a \e priori diagonal threshold values.
    </ol>
   */
  void compute();

  //! If compute() is completed, this query returns true, otherwise it returns false.
  inline bool isComputed() const {
    return IsComputed_;
  }

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

  /// \brief Apply the support graph preconditioner to X, resulting in Y.
  ///
  /// \param X [in] Input multivector; "right-hand side" of the solve.
  /// \param Y [out] Output multivector; result of the solve.
  void
  apply (const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type,node_type>& X,
         Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one (),
         scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero ()) const;

  //! Tpetra::Map representing the domain of this operator.
  Teuchos::RCP<const map_type> getDomainMap () const;

  //! Tpetra::Map representing the range of this operator.
  Teuchos::RCP<const map_type> getRangeMap () const;

  //! Whether this object's apply() method can apply the transpose
  //! (or conjugate transpose, if applicable).
  bool hasTransposeApply () const;

  //@}
  //! \name Mathematical functions
  //@{

  //! Return the operator's communicator.
  Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

  //! The matrix to be preconditioned.
  Teuchos::RCP<const row_matrix_type> getMatrix () const;

  //! Returns the number of calls to initialize().
  int getNumInitialize() const;

  //! Returns the number of calls to compute().
  int getNumCompute() const;

  //! Returns the number of calls to apply().
  int getNumApply() const;

  //! The total time in seconds spent in all initialize() calls.
  double getInitializeTime() const;

  //! The total time in seconds spent in all compute() calls.
  double getComputeTime() const;

  //! The total time in seconds spent in all apply() calls.
  double getApplyTime() const;

  //! Get absolute threshold value
  inline magnitude_type getAbsoluteThreshold() const {
    return(Athresh_);
  }

  //! Get relative threshold value
  inline magnitude_type getRelativeThreshold() const {
    return(Rthresh_);
  }

  // @}
  //! \name Implementation of Teuchos::Describable
  //@{

  /** \brief Return a simple one-line description of this object. */
  std::string description() const;

  //! Print the object with some verbosity level to a Teuchos::nFancyOStream.
  void
  describe (Teuchos::FancyOStream& out,
            const Teuchos::EVerbosityLevel verbLevel =
            Teuchos::Describable::verbLevel_default) const;

  //@}

private:
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;
  typedef typename Teuchos::Array<local_ordinal_type>::size_type size_type;
  typedef Tpetra::MultiVector<scalar_type,
                              local_ordinal_type,
                              global_ordinal_type,
                              node_type> MV;
  typedef Tpetra::CrsMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> crs_matrix_type;

  //! Find the maximum weight spanning tree and construct a CrsMatrix with it
  void findSupport ();

  //! Copy constructor (declared private and undefined; may not be used)
  SupportGraph (const SupportGraph<MatrixType>& RHS);

  //! operator= (declared private and undefined; may not be used)
  SupportGraph<MatrixType>& operator= (const SupportGraph<MatrixType>& RHS);

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

  //@}
  // \name The matrix, its support graph, and the (local) solver
  //@{

  //! The matrix to be preconditioned; input to the constructor or setMatrix().
  Teuchos::RCP<const row_matrix_type> A_;

  /// \brief "Local filter" version of A_.
  ///
  /// SupportGraph only knows how to precondition a square matrix on a
  /// single process.  Thus, if A_ has multiple processes in its
  /// communicator, we instead apply the preconditioner to the "local
  /// filter" version of A.  See the documentation of LocalFilter for
  /// an explanation.
  Teuchos::RCP<const row_matrix_type> A_local_;

  Teuchos::RCP<crs_matrix_type> Support_;

  Teuchos::RCP<Amesos2::Solver<crs_matrix_type,
                               Tpetra::MultiVector<scalar_type,
                                                   local_ordinal_type,
                                                   global_ordinal_type,
                                                   node_type> > > solver_;

  //@}
  // \name Parameters (set by the input ParameterList)
  //@{

  magnitude_type Athresh_; //!< Absolute threshold
  magnitude_type Rthresh_; //!< Relative threshold
  int Randomize_;
  int NumForests_;
  double KeepDiag_;
  //@}
  // \name Other internal data
  //@{


  //! Total time in seconds for all successful calls to initialize().
  double InitializeTime_;
  //! Total time in seconds for all successful calls to compute().
  double ComputeTime_;
  //! Total timer in seconds for all successful calls to apply().
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
}; // class SupportGraph

} // namespace Ifpack2

#endif /* IFPACK2_SUPPORTGRAPH_HPP */
