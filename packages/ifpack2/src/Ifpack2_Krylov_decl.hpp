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

/// \file Ifpack2_Krylov_decl.hpp
/// \brief Declaration of Ifpack2::DeprecatedAndMayDisappearAtAnyTime::Krylov class.
/// \author Paul Tsuji
///
/// \warning This file and its contents are DEPRECATED.
///   DO NOT USE THEM ANY MORE.  THEY MAY DISAPPEAR AT ANY TIME.

#ifndef IFPACK2_KRYLOV_DECL_HPP
#define IFPACK2_KRYLOV_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#ifdef HAVE_IFPACK2_DEPRECATED_CODE

#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Details_CanChangeMatrix.hpp"

#include "BelosSolverManager.hpp"
#include "BelosTpetraAdapter.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include <type_traits>

namespace Teuchos {
  class ParameterList; // forward declaration
}

namespace Ifpack2 {
namespace DeprecatedAndMayDisappearAtAnyTime {

  /// \struct BelosScalarType
  /// \brief Traits class for determining the scalar type to use for Belos
  ///
  /// \warning This is an implementation detail of Ifpack2.  Users
  ///   must not rely on this struct existing or on any details of its
  ///   interface.  It may go away at any time.
  ///
  /// This exists to allow a ScalarType to override what scalar type
  /// is used for Belos, if those are not equal, by specializing this
  /// class.
  template <typename ScalarType>
  struct BelosScalarType {
    typedef ScalarType type;
  };

  /// \class Krylov
  /// \brief Wrapper for iterative linear solvers.  DO NOT USE THIS
  ///   CLASS.  IT IS DEPRECATED AND MAY GO AWAY AT ANY TIME.
  /// \tparam MatrixType A specialization of Tpetra::RowMatrix.
  ///
  /// \warning This class has been DEPRECATED.  IT HAS BEEN DEPRECATED
  ///   SINCE APRIL 2015 SO PLEASE DON'T USE IT.  IT MAY GO AWAY AT
  ///   ANY TIME.  IT WON'T EVEN BUILD BY DEFAULT.
  ///
  /// Please do not use this class.  If you insist on enabling this
  /// long deprecated class, you must first set the CMake option
  /// Ifpack2_ENABLE_DEPRECATED_CODE:BOOL=ON, then rerun CMake, then
  /// rebuild and reinstall Trilinos.  Otherwise, this class does not
  /// exist and cannot be used.
  template<class MatrixType>
  class IFPACK2_DEPRECATED Krylov :
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
    // \name Public typedefs
    //@{

    //! The type of the entries of the input MatrixType.
    typedef typename MatrixType::scalar_type scalar_type;

    //! The scalar type used by Belos (which may not be scalar_type)
    typedef typename BelosScalarType<scalar_type>::type belos_scalar_type;

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

    static_assert(std::is_same<MatrixType, row_matrix_type>::value, "Ifpack2::Krylov: The template parameter MatrixType must be a Tpetra::RowMatrix specialization.  Please don't use Tpetra::CrsMatrix (a subclass of Tpetra::RowMatrix) here anymore.  The constructor can take either a RowMatrix or a CrsMatrix just fine.");

    //! Type of the Ifpack2::Preconditioner specialization from which this class inherits.
    typedef Ifpack2::Preconditioner<scalar_type,
                                    local_ordinal_type,
                                    global_ordinal_type,
                                    node_type> prec_type;

    //@}
    // \name Constructors and Destructors
    //@{

    //! Constructor that takes a Tpetra::RowMatrix.
    explicit Krylov (const Teuchos::RCP<const row_matrix_type>& A);

    //! Destructor
    virtual ~Krylov ();

    //@}
    //! \name Methods for setting parameters and initialization
    //@{

    /// \brief Set the preconditioner's parameters.
    ///
    /// This preconditioner accepts the following parameters:
    ///   - "krylov: iteration type" (\c std::string)
    ///   - "krylov: number of iterations" (\c int)
    ///   - "krylov: residual tolerance" (\c magnitude_type)
    ///
    /// The "krylov: iteration type" parameter specifies the name of
    /// the iterative linear solver to use.
    ///
    /// Note: Because some of the iterative solvers in Belos
    /// are not currently supported for complex types, the
    /// BelosSolverFactory is not used; Ifpack2::Krylov
    /// either uses Block GMRES or Block CG.
    void setParameters (const Teuchos::ParameterList& params);

    //! Do any initialization that depends on the input matrix's structure.
    void initialize ();

    //! Return \c true if initialize() completed successfully, else \c false.
    inline bool isInitialized () const {
      return IsInitialized_;
    }

    //! Do any initialization that depends on the input matrix's values.
    void compute ();

    //! Return \c true if compute() completed successfully, else \c false.
    inline bool isComputed() const {
      return IsComputed_;
    }

    //@}
    //! \name Implementation of Ifpack2::Details::CanChangeMatrix
    //@{

    /// \brief Change the matrix to be preconditioned.
    ///
    /// \param[in] A The new matrix.
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
    virtual void setMatrix (const Teuchos::RCP<const row_matrix_type>& A);

    //@}
    //! @name Implementation of Tpetra::Operator
    //@{

    //! Apply the preconditioner to X, putting the result in Y.
    void
    apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
           Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
           Teuchos::ETransp mode = Teuchos::NO_TRANS,
           scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
           scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

    //! Tpetra::Map representing the domain of this operator.
    Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> > getDomainMap() const;

    //! Tpetra::Map representing the range of this operator.
    Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> > getRangeMap() const;

    //! Whether this object's apply() method can apply the transpose (or conjugate transpose, if applicable).
    bool hasTransposeApply() const;

    //@}
    //! \name Mathematical functions.
    //@{

    //! Return the operator's communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

    //! Returns a reference to the matrix to be preconditioned.
    Teuchos::RCP<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > getMatrix() const;

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

    //! @name Overridden from Teuchos::Describable
    //@{

    /** \brief Return a simple one-line description of this object. */
    std::string description() const;

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

    //@}

  private:
    typedef Teuchos::ScalarTraits<scalar_type> STS;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;

    //! Copy constructor (should never be used)
    Krylov (const Krylov<MatrixType>& RHS);

    //! operator= (should never be used)
    Krylov<MatrixType>& operator= (const Krylov<MatrixType>& RHS);

    /// \brief The input matrix to be preconditioned.
    ///
    /// This may be null.  If this is null, initialize(), compute(),
    /// and apply() may not be called.
    Teuchos::RCP<const row_matrix_type> A_;

    //! General CG/GMRES parameters

    /// \brief Belos iterative linear solver name.
    ///
    /// Default is "GMRES".
    std::string iterationType_;

    //! Number of iterations
    int numIters_;

    //! Residual Tolerance
    magnitude_type resTol_;

    //! Block size
    int BlockSize_;

    //! If true, the starting solution is always the zero vector.
    bool ZeroStartingSolution_;

    //! Preconditioner Type
    // 1 for relaxation
    // 2 for ILUT
    // 3 for RILUK
    // 4 for Chebyshev
    int PreconditionerType_;

    /// \brief Inner preconditioner parameters.
    ///
    /// The "inner preconditioner" for Krylov means the preconditioner
    /// for the Krylov subspace method.
    Teuchos::ParameterList precParams_;

    //! \c true if \c this object has been initialized
    bool IsInitialized_;
    //! \c true if \c this object has been computed
    bool IsComputed_;
    //! Contains the number of successful calls to Initialize().
    int NumInitialize_;
    //! Contains the number of successful call to Compute().
    int NumCompute_;
    //! Contains the number of successful call to apply().
    mutable int NumApply_;
    //! Contains the time for all successful calls to Initialize().
    double InitializeTime_;
    //! Contains the time for all successful calls to Compute().
    double ComputeTime_;
    //! Contains the time for all successful calls to apply().
    mutable double ApplyTime_;

    //! Belos' encapsulation of the linear problem to solve.
    Teuchos::RCP<Belos::LinearProblem<belos_scalar_type,
                                      Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>,
                                      Tpetra::Operator<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > > belosProblem_;

    //! The Belos solver (implementation of the Krylov method).
    Teuchos::RCP<Belos::SolverManager<belos_scalar_type,
                                      Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>,
                                      Tpetra::Operator<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > > belosSolver_;

    //! The inner preconditioner (the preconditioner for the Krylov method).
    Teuchos::RCP<prec_type> ifpack2_prec_;
  };

} // namespace DeprecatedAndMayDisappearAtAnyTime
} // namespace Ifpack2

#endif // HAVE_IFPACK2_DEPRECATED_CODE
#endif // IFPACK2_KRYLOV_DECL_HPP
