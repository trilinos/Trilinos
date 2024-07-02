// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Ifpack2_DatabaseSchwarz_decl.hpp
/// \brief Declaration of DatabaseSchwarz class

#ifndef IFPACK2_DATABASESCHWARZ_DECL_HPP
#define IFPACK2_DATABASESCHWARZ_DECL_HPP

#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Details_CanChangeMatrix.hpp"

// We only need the declaration here, and only for the method
// getCrsMatrix.  I would very much prefer that this method not exist.
// Furthermore, it is both unsafe (MatrixType need not be CrsMatrix)
// and completely redundant (just call getMatrix() and do the
// dynamic_cast yourself).
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Tpetra_CrsMatrix_decl.hpp"

#include <type_traits>

namespace Ifpack2 {

/// \class DatabaseSchwarz
/// \brief Overlapping Schwarz where redundant patches are not stored explicitly.
/// \tparam MatrixType The type of matrix to use.
///
/// This class implements an overlapping Schwarz assuming it has already
/// been provided overlapping data containers. Additionally, it is assumed
/// there is at least one row with PatchSize_ nonzeros for each patch.
/// This means that a matrix obtained from piecewise linear FEMs in 2d,
/// for example, will not have rows with 4 nonzeros for each patch.
/// However, piecewise quadratic FEMs in 2d will have one row with 9 nonzeros
/// for each cell in the mesh, meaning all the patches will be detected.
/// For a matrix corresponding to a Taylor-Hood discretization in 2d
/// (quadratic velocities and linear pressures), the patch size would be 22.
///
/// The general algorithm proceeds as follows:
/// <ol>
/// <li> The rows of A are analyzed sequentially for any row with num_entries == PatchSize_ </li>
/// <li> If the current row corresponds to a DOF that has already been "visited", it is skipped </li>
/// <li> Then, all nonzero indices belonging to the row are marked as "visited" </li>
/// <li> The local patch matrix is formed and compared to a database of previous patch matrices.
///      If any patch matrix in the database has an l1 distance to the current patch matrix less than tol,
///      the current patch is not matrix is not stored, but an index pointing it to a replacement is stored instead </li>
/// <li> Finally, if this patch matrix is not sufficiently close to one that has already been seen, it is added to the database </li>
/// <li> The compute phase then inverts only the patch matcies in the database,
///      and the apply phase loops over all patches and applies the inverse of each appropriate patch matrix </li>
/// </ol>
///
/// In general, there is a noticeable speedup when using this method
/// compared to a typical method. This speedup may be further improved
/// by using more advanced linear algebra interfaces such as batched
/// Kokkos solves instead of the current LAPACK approach.
template<class MatrixType>
class DatabaseSchwarz :
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

  //! The template parameter of this class.
  typedef MatrixType matrix_type;

  //! The type of the entries of the input MatrixType.
  typedef typename MatrixType::scalar_type scalar_type;

  //! The type of local indices in the input MatrixType.
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;

  //! The type of global indices in the input MatrixType.
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;

  //! The Kokkos::Device specialization used by the input MatrixType.
  typedef typename MatrixType::node_type::device_type device_type;

  //! The Node type used by the input MatrixType.
  typedef typename MatrixType::node_type node_type;

  //! The type of the magnitude (absolute value) of a matrix entry.
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  /// \brief The Tpetra::RowMatrix specialization matching MatrixType.
  ///
  /// MatrixType must be a Tpetra::RowMatrix specialization.  This
  /// typedef will always be a Tpetra::RowMatrix specialization.
  typedef Tpetra::RowMatrix<scalar_type, local_ordinal_type,
                            global_ordinal_type, node_type> row_matrix_type;

  static_assert (std::is_same<MatrixType, row_matrix_type>::value,
                 "Ifpack2::DatabaseSchwarz: MatrixType must be a Tpetra::RowMatrix "
                 "specialization.  Don't use Tpetra::CrsMatrix here.");

  //! The Tpetra::Map specialization matching MatrixType.
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

  /// \brief The Tpetra::Vector specialization matching MatrixType.
  ///
  /// If you wish to supply setParameters() a precomputed vector of
  /// diagonal entries of the matrix, use a pointer to an object of
  /// this type.
  typedef Tpetra::Vector<scalar_type, local_ordinal_type,
                         global_ordinal_type, node_type> vector_type;

  //@}
  // \name Constructors and destructors
  //@{

  /// \brief Constructor.
  ///
  /// \param[in] A The sparse matrix to which to apply DatabaseSchwarz
  ///   iteration.  The matrix A must be square, and its domain Map
  ///   and range Map must be the same.  The latter means that the
  ///   vectors x and y in the sparse matrix-vector product y = A*x
  ///   must both have the same distribution over process(es).
  ///
  explicit DatabaseSchwarz(const Teuchos::RCP<const row_matrix_type>& A);

  /// \brief Constructor.
  ///
  /// \param[in] A The sparse matrix to which to apply DatabaseSchwarz
  ///   iteration.  The matrix A must be square, and its domain Map
  ///   and range Map must be the same.  The latter means that the
  ///   vectors x and y in the sparse matrix-vector product y = A*x
  ///   must both have the same distribution over process(es).
  /// \param[in] params The parameterlist containing settings for the
  ///   object, such as the patch size to search for.
  ///
  DatabaseSchwarz (const Teuchos::RCP<const row_matrix_type>& A,
                   Teuchos::ParameterList& params);

  //! Destructor.
  virtual ~DatabaseSchwarz();

  //@}
  //! \name Preconditioner computation methods
  //@{

  /// \brief Set (or reset) parameters.
  void setParameters(const Teuchos::ParameterList& params);

  bool supportsZeroStartingSolution() { return true; }

  void setZeroStartingSolution(bool zeroStartingSolution);

  /// \brief Initialize the preconditioner.
  ///
  /// The compute() method will call initialize() automatically if it
  /// has not yet been called, so you do not normally need to call
  /// this.  However, it is correct to call initialize() yourself, and
  /// compute() will not call it again if it already has been called.
  void initialize();

  /// Whether the preconditioner has been successfully initialized
  /// (by calling initialize()).
  inline bool isInitialized() const {
    return IsInitialized_;
  }

  /// \brief (Re)compute the left scaling, and (if applicable)
  ///   estimate max and min eigenvalues of D_inv * A.
  void compute();

  /// Whether compute() has been called at least once.
  inline bool isComputed() const {
    return IsComputed_;
  }

  //@}
  //! \name Implementation of Ifpack2::Details::CanChangeMatrix
  //@{

  /// \brief Change the matrix to be preconditioned.
  virtual void
  setMatrix(const Teuchos::RCP<const row_matrix_type>& A);

  //@}
  //! \name Implementation of Tpetra::Operator
  //@{

  /// \brief Apply the preconditioner to X, returning the result in Y.
  /// Y = alpha*Op(A)*X + beta*Y
  void
  apply(const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
        Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
        Teuchos::ETransp mode = Teuchos::NO_TRANS,
        scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
        scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //! The Tpetra::Map representing the domain of this operator.
  Teuchos::RCP<const map_type> getDomainMap() const;

  //! The Tpetra::Map representing the range of this operator.
  Teuchos::RCP<const map_type> getRangeMap() const;

  //! Whether it's possible to apply the transpose of this operator.
  bool hasTransposeApply() const;

  /// \brief Compute Y = Op(A)*X, where Op(A) is either A, \f$A^T\f$, or \f$A^H\f$.
  void
  applyMat(const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
           Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
           Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

  //@}
  //! \name Attribute accessor methods
  //@{

  //! The communicator over which the matrix is distributed.
  Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

  //! The matrix for which this is a preconditioner.
  Teuchos::RCP<const row_matrix_type> getMatrix() const;

  //! The matrix for which this is a preconditioner.
  Teuchos::RCP<const row_matrix_type> A_;

  /// \brief Attempt to return the matrix A as a Tpetra::CrsMatrix.
  ///
  /// This class does not require that A be a Tpetra::CrsMatrix.
  /// If it is NOT, this method will return Teuchos::null.
  Teuchos::RCP<const Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> >
  getCrsMatrix() const;

  //! The total number of floating-point operations taken by all calls to compute().
  double getComputeFlops() const;

  //! The total number of floating-point operations taken by all calls to apply().
  double getApplyFlops() const;

  //! The total number of successful calls to initialize().
  int getNumInitialize() const;

  //! The total number of successful calls to compute().
  int getNumCompute() const;

  //! The total number of successful calls to apply().
  int getNumApply() const;

  //! The total time spent in all calls to initialize().
  double getInitializeTime() const;

  //! The total time spent in all calls to compute().
  double getComputeTime() const;

  //! The total time spent in all calls to apply().
  double getApplyTime() const;

  //! Get a rough estimate of cost per iteration
  size_t getNodeSmootherComplexity() const;  

  //@}
  //! @name Implementation of Teuchos::Describable
  //@{

  //! A simple one-line description of this object.
  std::string description() const;

  //! Print the object with some verbosity level to a Teuchos::FancyOStream.
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

  //@}

private:

  //! Implementation of parameter setting
  void setParametersImpl(Teuchos::ParameterList& params);

  //! Abbreviation for the Teuchos::ScalarTraits specialization for scalar_type.
  typedef Teuchos::ScalarTraits<typename MatrixType::scalar_type> STS;

  //! Abbreviation for the Tpetra::MultiVector specialization used in methods like apply().
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> MV;

  //! Copy constructor (use is syntactically forbidden)
  DatabaseSchwarz(const DatabaseSchwarz<MatrixType>&);

  //! Assignment operator (use is syntactically forbidded)
  DatabaseSchwarz<MatrixType>& operator= (const DatabaseSchwarz<MatrixType>&);

  //! \name Internal state
  //@{

  //! If \c true, initialize() has completed successfully.
  bool IsInitialized_;
  //! If \c true, compute() has completed successfully.
  bool IsComputed_;
  //! The total number of successful calls to initialize().
  int NumInitialize_;
  //! The total number of successful calls to compute().
  int NumCompute_;
  /// \brief The total number of successful calls to apply().
  ///
  /// This is "mutable" because apply() is a const method; apply() is
  /// const because it is declared this way in Tpetra::Operator.
  mutable int NumApply_;
  //! The total time in seconds over all calls to initialize().
  double InitializeTime_;
  //! The total time in seconds over all calls to compute().
  double ComputeTime_;
  /// \brief The total time in seconds over all calls to apply().
  ///
  /// This is "mutable" because apply() is a const method; apply() is
  /// const because it is declared this way in Tpetra::Operator.
  mutable double ApplyTime_;
  //! The total number of floating-point operations over all calls to compute().
  double ComputeFlops_;
  /// \brief The total number of floating-point operations over all calls to apply().
  ///
  /// This is "mutable" because apply() is a const method; apply() is
  /// const because it is declared this way in Tpetra::Operator.
  mutable double ApplyFlops_;

  /// Size of inner patches to search for
  local_ordinal_type PatchSize_;

  /// Number of found patches
  mutable size_t NumPatches_;

  /// Tolerance at which to consider two patches "equal"
  double PatchTolerance_;

  /// Boolean to indicate we'd like to skip all database comparisons and instead invert all patches
  bool SkipDatabase_;

  /// Boolean to indicate we'd like to skip all database comparisons and instead invert all patches
  bool Verbose_;

  /// A vector of vectors where each row of PatchIndices_ is the indices corresponding to a patch
  mutable std::vector<std::vector<typename row_matrix_type::local_ordinal_type> > PatchIndices_;

  /// Size of database
  mutable size_t DatabaseSize_;

  /// Database patches
  mutable std::vector<Teuchos::RCP<typename Teuchos::SerialDenseMatrix<typename row_matrix_type::local_ordinal_type,typename row_matrix_type::scalar_type> > > DatabaseMatrices_;

  /// Database indices
  std::vector<int> DatabaseIndices_;

  /// Patch weights
  std::vector<magnitude_type> Weights_;

  /// Pivots used for LAPACK
  mutable Teuchos::Array<int> ipiv_;

  //@}
}; // class DatabaseSchwarz

} // namespace Ifpack2

#endif // IFPACK2_DATABASESCHWARZ_DECL_HPP

