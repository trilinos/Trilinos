// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_DIAGONAL_DECL_HPP
#define IFPACK2_DIAGONAL_DECL_HPP

#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Details_CanChangeMatrix.hpp"
#include "Tpetra_CrsMatrix_decl.hpp"
#include <type_traits>

namespace Ifpack2 {

/// \class Diagonal
/// \brief Diagonal preconditioner.
/// \tparam MatrixType A specialization of Tpetra::RowMatrix.
///
/// This class wraps a Tpetra::Vector as a diagonal preconditioner.
/// The preconditioner is defined as
/// \f[
/// z_i = D_i r_i,
/// \f]
/// where \f$r\f$ is the Vector to be preconditioned, \f$z\f$ is the
/// Vector result of applying the preconditioner to \f$r\f$, and
/// \f$D_i\f$ is the i-th element of the scaling vector.
///
/// When Diagonal is constructed with a matrix, \f$D\f$ contains the
/// inverted diagonal elements from the matrix.  When Diagonal is
/// constructed with a Tpetra::Vector, \f$D\f$ is the caller-supplied
/// vector.
template<class MatrixType>
class Diagonal :
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
  typedef typename MatrixType::scalar_type scalar_type;
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;
  typedef typename MatrixType::node_type node_type;
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  //! Tpetra::RowMatrix specialization used by this class.
  typedef Tpetra::RowMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> row_matrix_type;

  static_assert(std::is_same<MatrixType, row_matrix_type>::value, "Ifpack2::Diagonal: The template parameter MatrixType must be a Tpetra::RowMatrix specialization.  Please don't use Tpetra::CrsMatrix (a subclass of Tpetra::RowMatrix) here anymore.  The constructor can take either a RowMatrix or a CrsMatrix just fine.");

  //! Tpetra::CrsMatrix specialization used by this class.
  typedef Tpetra::CrsMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> crs_matrix_type;
  //! Tpetra::Vector specialization used by this class.
  typedef Tpetra::Vector<scalar_type,
                         local_ordinal_type,
                         global_ordinal_type,
                         node_type> vector_type;
  //! Tpetra::Map specialization used by this class.
  typedef Tpetra::Map<local_ordinal_type,
                      global_ordinal_type,
                      node_type> map_type;

  /// \brief Constructor that takes a Tpetra::RowMatrix.
  ///
  /// \param A_in [in] The input matrix.
  Diagonal (const Teuchos::RCP<const row_matrix_type>& A);

  /// \brief Constructor that takes a Tpetra::CrsMatrix.
  ///
  /// This constructor exists to avoid "ambiguous constructor"
  /// warnings.  It does the same thing as the constructor that takes
  /// a Tpetra::RowMatrix.
  ///
  /// \param A_in [in] The input matrix.
  Diagonal (const Teuchos::RCP<const crs_matrix_type>& A_in);

  /// \brief Constructor that accepts a Tpetra::Vector of inverse diagonal entries.
  ///
  /// \param diag [in] Vector of inverse diagonal entries.
  ///
  /// If your compiler complains about this constructor being ambigous
  /// with the other constructor overload, instead call the
  /// free-standing function Ifpack2::createDiagonalPreconditioner
  /// which is located at the bottom of this header file.  (This issue
  /// may arise if this constructor is called with a
  /// <tt>Teuchos::RCP<Tpetra::Vector></tt> that isn't const-qualified
  /// exactly as declared here.)
  Diagonal (const Teuchos::RCP<const vector_type>& diag);

  //! Destructor
  virtual ~Diagonal ();

  //! Sets parameters on this object.
  /**
    Currently doesn't need any parameters.
  */
  void setParameters (const Teuchos::ParameterList& params);

  //! Initialize
  void initialize ();

  //! Returns \c true if the preconditioner has been successfully initialized.
  bool isInitialized () const {
    return isInitialized_;
  }

  //! Compute the preconditioner.
  void compute ();

  //! Return true if compute() has been called.
  bool isComputed () const {
    return isComputed_;
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
  /// Calling this method with a matrix different than the current
  /// matrix resets the preconditioner's state.  After calling this
  /// method with a nonnull input, you must first call initialize()
  /// and compute() (in that order) before you may call apply().
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

  /// \brief Apply the preconditioner to X, putting the result in Y.
  ///
  /// If the result of applying this preconditioner to a vector X is
  /// \f$F \cdot X$, then this method computes \f$\beta Y + \alpha F \cdot X\f$.
  /// The typical case is \f$\beta = 0\f$ and \f$\alpha = 1\f$.
  void
  apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
         Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
         scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //! The Tpetra::Map representing this operator's domain.
  Teuchos::RCP<const map_type> getDomainMap () const;

  //! The Tpetra::Map representing this operator's range.
  Teuchos::RCP<const map_type> getRangeMap () const;

  //@}
  //! \name Attribute accessor methods
  //@{

  //! Return the communicator associated with this matrix operator.
  //Teuchos::RCP<const Teuchos::Comm<int> > getComm () const;

  /// \brief The original input matrix to be preconditioned.
  ///
  /// This could be null, for example if the user created this object
  /// using the constructor that takes a Tpetra::Vector, or if the
  /// user called setMatrix() with a null input.
  Teuchos::RCP<const row_matrix_type> getMatrix () const {
    return matrix_;
  }

  //! Return the number of flops in the computation phase.
  double getComputeFlops() const;

  //! Return the number of flops for the application of the preconditioner.
  double getApplyFlops() const;

  //! Return the number of calls to initialize().
  int getNumInitialize() const;

  //! Return the number of calls to compute().
  int getNumCompute() const;

  //! Return the number of calls to apply().
  int getNumApply() const;

  //! Return the time spent in initialize().
  double getInitializeTime() const;

  //! Return the time spent in compute().
  double getComputeTime() const;

  //! Return the time spent in apply().
  double getApplyTime() const;

  //@}
  //! @name Implementation of Teuchos::Describable
  //@{

  //! Return a one-line description of this object.
  std::string description() const;

  //! Print the object with some verbosity level to an FancyOStream object.
  void
  describe (Teuchos::FancyOStream& out,
            const Teuchos::EVerbosityLevel verbLevel =
            Teuchos::Describable::verbLevel_default) const;
  //@}

private:
  //! Reset the preconditioner's state.  Call in setMatrix().
  void reset ();

  //! The input matrix provided by the user.
  Teuchos::RCP<const row_matrix_type> matrix_;

  /// \brief Inverse diagonal provided by the user.
  ///
  /// This is only nonnull if this Diagonal instance was created using
  /// the constructor that takes a pointer to a Tpetra::Vector.
  Teuchos::RCP<const vector_type> userInverseDiag_;

  //! The vector of inverse diagonal entries to use in apply().
  Teuchos::RCP<const vector_type> inverseDiag_;

  typedef Kokkos::View<size_t*, typename node_type::device_type> offsets_type;
  offsets_type offsets_;

  double initializeTime_;
  double computeTime_;
  mutable double applyTime_;

  int numInitialize_;
  int numCompute_;
  mutable int numApply_;

  bool isInitialized_;
  bool isComputed_;
};

/** Function to construct a Diagonal preconditioner with Vector input.
* This is just a nonmember function version of Diagonal's constructor.
*
* Example usage:
*
* \code
* using Teuchos::RCP;
* typedef Tpetra::RowMatrix<> row_matrix_type;
* typedef Tpetra::Vector<> vec_type;
* typedef Tpetra::Preconditioner<> prec_type;
*
* RCP<vec_type> D = ...
* RCP<prec_type> P =
*   Ifpack2::createDiagonalPreconditioner<row_matrix_type, vec_type> (D);
* \endcode
*/
template<class MatrixType, class VectorType>
Teuchos::RCP<Ifpack2::Diagonal<Tpetra::RowMatrix<typename MatrixType::scalar_type,
                                                 typename MatrixType::local_ordinal_type,
                                                 typename MatrixType::global_ordinal_type,
                                                 typename MatrixType::node_type> > >
createDiagonalPreconditioner (const Teuchos::RCP<const VectorType>& invdiag)
{
  typedef Tpetra::RowMatrix<typename MatrixType::scalar_type,
    typename MatrixType::local_ordinal_type,
    typename MatrixType::global_ordinal_type,
    typename MatrixType::node_type> row_matrix_type;

  return Teuchos::rcp (new Ifpack2::Diagonal<row_matrix_type> (invdiag));
}

}//namespace Ifpack2

#endif
