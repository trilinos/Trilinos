/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_LOCALSPARSETRIANGULARSOLVER_DECL_HPP
#define IFPACK2_LOCALSPARSETRIANGULARSOLVER_DECL_HPP

#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Details_CanChangeMatrix.hpp"
#include "Tpetra_CrsMatrix_fwd.hpp"
#include "Teuchos_FancyOStream.hpp"
#include <type_traits>

namespace Ifpack2 {

/// \brief "Preconditioner" that solves local sparse triangular systems.
/// \tparam MatrixType Specialization of Tpetra::RowMatrix.
///
/// This class solves local sparse triangular systems.  "Local" means
/// "per MPI process."  The matrix itself may be distributed across
/// multiple MPI processes, but this class works on each MPI process'
/// part of the matrix, and the input and output multivectors,
/// separately.  (See this class' constructor for details.)
///
/// This effectively assumes that the global matrix is block diagonal.
/// Triangular solves usually imply that these blocks are square.  If
/// a particular triangular solver knows how to deal with nonsquare
/// blocks, though, this is allowed.
///
/// The implementation currently requires that the input
/// Tpetra::RowMatrix actually be a Tpetra::CrsMatrix.  This lets us
/// optimize without necessarily copying data structures.  We may
/// relax this restriction in the future.
///
/// If you are writing a new Ifpack2 class that needs to solve local
/// sparse triangular systems stored as Tpetra::CrsMatrix, use this
/// class <i>only</i>.
template<class MatrixType>
class LocalSparseTriangularSolver :
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
  //! Type of the entries of the input matrix.
  typedef typename MatrixType::scalar_type scalar_type;
  //! Type of the local indices of the input matrix.
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;
  //! Type of the global indices of the input matrix.
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;
  //! Node type of the input matrix.
  typedef typename MatrixType::node_type node_type;

  //! Type of the absolute value (magnitude) of a \c scalar_type value.
  typedef typename MatrixType::mag_type magnitude_type;
  //! Specialization of Tpetra::Map used by this class.
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;
  //! Specialization of Tpetra::RowMatrix used by this class.
  typedef Tpetra::RowMatrix<scalar_type, local_ordinal_type,
                            global_ordinal_type, node_type> row_matrix_type;
  //! Specialization of Tpetra::CrsMatrix used by this class.
  typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type,
                            global_ordinal_type, node_type> crs_matrix_type;

  static_assert (std::is_same<MatrixType, row_matrix_type>::value,
                 "Ifpack2::LocalSparseTriangularSolver: The template parameter "
                 "MatrixType must be a Tpetra::RowMatrix specialization.  "
                 "Please don't use Tpetra::CrsMatrix (a subclass of "
                 "Tpetra::RowMatrix) here anymore.  The constructor can take "
                 "either a RowMatrix or a CrsMatrix just fine.");

  /// \brief Constructor
  ///
  /// \param A [in] The input sparse matrix.  Though its type is
  ///   Tpetra::RowMatrix for consistency with other Ifpack2 solvers,
  ///   this must be a Tpetra::CrsMatrix specialization.
  ///
  /// The input matrix A may be distributed across multiple MPI
  /// processes.  This class' apply() method will use A's Import
  /// object, if it exists, to Import the input MultiVector from the
  /// domain Map to the column Map.  It will also use A's Export
  /// object, if it exists, to Export the output MultiVector from the
  /// row Map to the range Map.  Thus, to avoid MPI communication and
  /// local permutations, construct A so that the row, column, range,
  /// and domain Maps are all identical.
  ///
  /// On the other hand, you may encode local permutations in the
  /// matrix's Maps, and let Import and/or Export execute them for
  /// you.
  ///
  /// The input matrix must have local properties corresponding to the
  /// way in which one wants to solve.  ("Local" means "to each MPI
  /// process.")  For example, if one wants to solve lower triangular
  /// systems with an implicit unit diagonal, the matrix A must have
  /// these properties.  If the matrix does not know whether it has
  /// these properties and the user does not specify them, then this
  /// class is responsible for figuring out whether the matrix has
  /// those properties.
  LocalSparseTriangularSolver (const Teuchos::RCP<const row_matrix_type>& A);

  /// \brief Constructor that takes an optional debug output stream.
  ///
  /// \param A [in] The input sparse matrix.  Though its type is
  ///   Tpetra::RowMatrix for consistency with other Ifpack2 solvers,
  ///   this must be a Tpetra::CrsMatrix specialization.
  ///
  /// \param out [in/out] Optional debug output stream.  If nonnull,
  ///   this solver will print copious debug output to the stream.
  LocalSparseTriangularSolver (const Teuchos::RCP<const row_matrix_type>& A,
                               const Teuchos::RCP<Teuchos::FancyOStream>& out);

  /// \brief Constructor that takes no input matrix.
  ///
  /// Call setMatrix to initialize the matrix.
  LocalSparseTriangularSolver ();

  /// \brief Constructor that takes no input matrix.
  ///
  /// Call setMatrix to initialize the matrix.
  ///
  /// \param unused [in] Disambiguate the constructor so that the ctor taking A
  ///   can cast A implicitly. This ctor is meant for debugging, so I prefer to
  ///   make it a bit awkward than require careful explicit casting of A in the
  ///   other ctor.
  ///
  /// \param out [in/out] Optional debug output stream.  If nonnull,
  ///   this solver will print copious debug output to the stream.
  LocalSparseTriangularSolver (const bool /* unused */, const Teuchos::RCP<Teuchos::FancyOStream>& out);

  //! Destructor (virtual for memory safety).
  virtual ~LocalSparseTriangularSolver ();

  /// \brief Set this object's parameters.
  ///
  /// \param plist [in] List of parameters.
  ///
  /// - "trisolver: reverse U" (\c bool): reverse storage for upper triangular matrices
  ///   to be more cache-efficient
  ///
  /// If Trilinos_ENABLE_ShyLU_NodeHTS=TRUE, then these parameters are available:
  ///   - "trisolver: type" (\c string): One of {"Internal" (default), "HTS"}.
  ///   - "trisolver: block size" (\c int): The triangular matrix can usefully be
  ///     thought of as being blocked int little blocks of this size. Default
  ///     is 1.
  void setParameters (const Teuchos::ParameterList& params);

  /// \brief "Symbolic" phase of setup
  ///
  /// Call this before calling compute() or apply() if the matrix
  /// object itself changes, or if the matrix's graph structure may
  /// have changed.
  void initialize ();

  //! Return \c true if the preconditioner has been successfully initialized.
  inline bool isInitialized () const {
    return isInitialized_;
  }

  /// \brief "Numeric" phase of setup
  ///
  /// Call this before calling apply() if the values in the matrix may
  /// have changed.
  void compute ();

  //! Return true if compute() has been called.
  inline bool isComputed () const {
    return isComputed_;
  }

  //! @name Implementation of Tpetra::Operator
  //@{

  /// \brief Apply the preconditioner to X, and put the result in Y.
  ///
  /// If this preconditioner is an operator M, this method computes
  /// <tt> Y := beta * Y + alpha * (M * X) </tt>.
  ///
  /// \param X [in] MultiVector to which to apply the preconditioner.
  /// \param Y [in/out] On input: Initial guess, if applicable.
  ///   On output: Result of applying the preconditioner. Y may alias X.
  /// \param mode [in] Whether to apply the transpose (Teuchos::TRANS)
  ///   or conjugate transpose (Teuchos::CONJ_TRANS).  The default
  ///   (Teuchos::NO_TRANS) is not to apply the transpose or conjugate
  ///   transpose.
  /// \param alpha [in] Scalar factor by which to multiply the result
  ///   of applying this operator to X.
  /// \param beta [in] Scalar factor for Y.
  void
  apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
         Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one (),
         scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero ()) const;

  //! The domain of this operator.
  Teuchos::RCP<const map_type> getDomainMap () const;

  //! The range of this operator.
  Teuchos::RCP<const map_type> getRangeMap () const;

  /// \brief Apply the original input matrix.
  ///
  /// \param X [in] MultiVector input.
  /// \param Y [in/out] Result of applying Op(A) to X, where Op(A) is
  ///   A (mode == Teuchos::NO_TRANS), the transpose of A (<tt>mode ==
  ///   Teuchos::TRANS</tt>), or the conjugate transpose of A
  ///   (<tt>mode == Teuchos::CONJ_TRANS</tt>)
  /// \param mode [in] See above.
  void
  applyMat (const Tpetra::MultiVector<scalar_type, local_ordinal_type,
              global_ordinal_type, node_type>& X,
            Tpetra::MultiVector<scalar_type, local_ordinal_type,
              global_ordinal_type, node_type>& Y,
            Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

  //! This operator's communicator.
  Teuchos::RCP<const Teuchos::Comm<int> > getComm () const;

  //! The original input matrix.
  Teuchos::RCP<const row_matrix_type> getMatrix () const {
    return A_;
  }

  //! Return the number of flops in the computation phase.
  double getComputeFlops () const;

  //! Return the number of flops for the application of the preconditioner.
  double getApplyFlops () const;

  //! Return the number of calls to initialize().
  int getNumInitialize () const;

  //! Return the number of calls to compute().
  int getNumCompute () const;

  //! Return the number of calls to apply().
  int getNumApply () const;

  //! Return the time spent in initialize().
  double getInitializeTime () const;

  //! Return the time spent in compute().
  double getComputeTime () const;

  //! Return the time spent in apply().
  double getApplyTime () const;

  //@}
  //! @name Implementation of Teuchos::Describable
  //@{

  //! A one-line description of this object.
  std::string description() const;

  /// \brief Print this object with given verbosity to the given output stream.
  ///
  /// \param out [out] Output stream to which to print
  /// \param verbLevel [in] Verbosity level
  ///
  /// You may create a Teuchos::FancyOStream from any std::ostream.
  /// For example, to wrap std::cout in a FancyOStream, do this:
  /// \code
  /// Teuchos::RCP<Teuchos::FancyOStream> out =
  ///   Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cout));
  /// \endcode
  /// To wrap a new std::ostringstream in a FancyOStream, do this:
  /// \code
  /// auto osPtr = Teuchos::rcp (new std::ostringstream ());
  /// Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::getFancyOStream (osPtr);
  ///
  /// // ... use out ...
  ///
  /// // Print to std::cout whatever the std::ostringstream got.
  /// std::cout << osPtr->str () << std::endl;
  /// \endcode
  void
  describe (Teuchos::FancyOStream& out,
            const Teuchos::EVerbosityLevel verbLevel =
            Teuchos::Describable::verbLevel_default) const;

  /// \brief Set this preconditioner's matrix.
  ///
  /// After calling this method, you must call first initialize(),
  /// then compute(), before you may call apply().
  virtual void setMatrix (const Teuchos::RCP<const row_matrix_type>& A);

  //@}

private:
  //! The original input matrix.
  Teuchos::RCP<const row_matrix_type> A_;
  //! Debug output stream; may be null (not used in that case)
  Teuchos::RCP<Teuchos::FancyOStream> out_;
  //! The original input matrix, as a Tpetra::CrsMatrix.
  Teuchos::RCP<const crs_matrix_type> A_crs_;

  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> MV;
  mutable Teuchos::RCP<MV> X_colMap_;
  mutable Teuchos::RCP<MV> Y_rowMap_;

  bool isInitialized_;
  bool isComputed_;
  /// \brief True if and only if this class' internal storage
  ///   representation of the matrix is not the same as A_.
  ///
  /// If true, then one of two things has happened:
  ///
  /// <ol>
  /// <li> A_crs_ is actually a copy, with permuted / reversed storage. </li>
  /// <li> htsImpl_->initialize(*A_crs_) has been called. </li>
  /// </ol>
  bool isInternallyChanged_;
  bool reverseStorage_;

  mutable int numInitialize_;
  mutable int numCompute_;
  mutable int numApply_;

  double initializeTime_;
  double computeTime_;
  double applyTime_;

  //! Optional HTS implementation.
  class HtsImpl;
  Teuchos::RCP<HtsImpl> htsImpl_;

  /// \brief "L" if the matrix is locally lower triangular, "U" if the
  ///   matrix is locally upper triangular, or "N" if unknown or
  ///   otherwise.
  std::string uplo_;
  /// \brief "U" if the matrix is known to have an implicitly stored
  ///   unit diagonal, else "N".
  std::string diag_;

  /// \brief The purely local part of apply().
  ///
  /// This is where all implementation effort should go.  If you want
  /// to plug in a new triangular solver, put it here.  No MPI
  /// communication (use of Import or Export) happens here.
  ///
  /// \param X [in] Input MultiVector; distributed according to the
  ///   input matrix's column Map.
  /// \param Y [in/out] Output MultiVector; distributed according to
  ///   the input matrix's row Map.  On input: Initial guess, if
  ///   applicable.
  /// \param mode [in] Whether to apply the transpose (Teuchos::TRANS)
  ///   or conjugate transpose (Teuchos::CONJ_TRANS).  The default
  ///   (Teuchos::NO_TRANS) is not to apply the transpose or conjugate
  ///   transpose.
  /// \param alpha [in] Scalar factor by which to multiply the result
  ///   of applying this operator to X.
  /// \param beta [in] Scalar factor for Y.
  void
  localApply (const MV& X,
              MV& Y,
              const Teuchos::ETransp mode,
              const scalar_type& alpha,
              const scalar_type& beta) const;

  //! Replacement for Tpetra::CrsMatrix::localSolve.
  void
  localTriangularSolve (const MV& Y,
                        MV& X,
                        const Teuchos::ETransp mode) const;

  void initializeState();
};

} // namespace Ifpack2

#endif // IFPACK2_LOCALSPARSETRIANGULARSOLVER_DECL_HPP
