/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
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

/// \file Ifpack2_ILUT_decl.hpp
/// \brief Declaration of ILUT preconditioner

#ifndef IFPACK2_ILUT_DECL_HPP
#define IFPACK2_ILUT_DECL_HPP

#include <Ifpack2_ConfigDefs.hpp>
#include <Ifpack2_Preconditioner.hpp>
#include <Ifpack2_Details_CanChangeMatrix.hpp>
#include <Tpetra_CrsMatrix_decl.hpp>

#include <string>
#include <sstream>
#include <iostream>
#include <cmath>


namespace Teuchos {
  class ParameterList; // forward declaration
}

namespace Ifpack2 {

/// \class ILUT
/// \brief ILUT (incomplete LU factorization with threshold) of a
///   Tpetra sparse matrix
/// \tparam Specialization of Tpetra::RowMatrix (preferred) or
///   Tpetra::CrsMatrix
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

  //! Preserved only for backwards compatibility.  Please use "scalar_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::scalar_type Scalar;


  //! The type of local indices in the input MatrixType.
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;

  //! Preserved only for backwards compatibility.  Please use "local_ordinal_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::local_ordinal_type LocalOrdinal;


  //! The type of global indices in the input MatrixType.
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;

  //! Preserved only for backwards compatibility.  Please use "global_ordinal_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::global_ordinal_type GlobalOrdinal;


  //! The type of the Kokkos Node used by the input MatrixType.
  typedef typename MatrixType::node_type node_type;

  //! Preserved only for backwards compatibility.  Please use "node_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::node_type Node;


  //! The type of the magnitude (absolute value) of a matrix entry.
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  //! Preserved only for backwards compatibility.  Please use "magnitude_type".
  TEUCHOS_DEPRECATED typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitudeType;

  //! Type of the Tpetra::RowMatrix specialization that this class uses.
  typedef Tpetra::RowMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> row_matrix_type;

  //! Type of the Tpetra::CrsMatrix specialization that this class uses for the L and U factors.
  typedef Tpetra::CrsMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> crs_matrix_type;

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
  virtual ~ILUT();

  //@}
  //! \name Methods for setting up and computing the incomplete factorization
  //@{

  /// \brief Set preconditioner parameters.
  ///
  /// ILUT implements the following parameters:
  /// <ul>
  /// <li> "fact: ilut level-of-fill" (\c int)
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
  /// description.
  /// Each row has at most \f$level-of-fill + nnz(A(i; 1 : i))\f$
  /// nonzero elements.
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

  /// \brief Clear any previously computed factors.
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

  /// \brief Compute the condition number estimate and return its value.
  ///
  /// \warning This method is DEPRECATED.  It was inherited from
  ///   Ifpack, and Ifpack never clearly stated what this method
  ///   computes.  Furthermore, Ifpack's method just estimates the
  ///   condition number of the matrix A, and ignores the
  ///   preconditioner -- which is probably not what users thought it
  ///   did.  If there is sufficient interest, we might reintroduce
  ///   this method with a different meaning and a better algorithm.
  virtual magnitude_type TEUCHOS_DEPRECATED
  computeCondEst (CondestType CT = Cheap,
                  local_ordinal_type MaxIters = 1550,
                  magnitude_type Tol = 1e-9,
                  const Teuchos::Ptr<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > &Matrix_in = Teuchos::null);

  /// \brief Return the computed condition number estimate, or -1 if not computed.
  ///
  /// \warning This method is DEPRECATED.  See warning for computeCondEst().
  virtual magnitude_type TEUCHOS_DEPRECATED getCondEst() const {
    return Condest_;
  }

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

  /// \brief The level of fill.
  ///
  /// For ILUT, this means the maximum number of entries in each row
  /// of the resulting L and U factors (each considered separately),
  /// not including the diagonal entry in that row (which is always
  /// part of U).  This has a different meaning for ILUT than it does
  /// for ILU(k).
  inline int getLevelOfFill() const {
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
  size_t getNodeNumEntries() const;

  //@}
  //! \name Implementation of Teuchos::Describable
  //@{

  /** \brief Return a simple one-line description of this object. */
  std::string description() const;

  /** \brief Print the object with some verbosity level to an FancyOStream object. */
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

  //@}

private:
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;
  typedef typename Teuchos::Array<local_ordinal_type>::size_type size_type;

  //! Copy constructor (declared private and undefined; may not be used)
  ILUT (const ILUT<MatrixType>& RHS);

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
  //! L factor of the incomplete LU factorization of A_local_.
  Teuchos::RCP<crs_matrix_type> L_;
  //! U factor of the incomplete LU factorization of A_local_.
  Teuchos::RCP<crs_matrix_type> U_;

  //@}
  // \name Parameters (set by setParameters())
  //@{

  magnitude_type Athresh_; //!< Absolute threshold
  magnitude_type Rthresh_; //!< Relative threshold
  magnitude_type RelaxValue_; //!< Relax value
  int LevelOfFill_; //!< Max fill level
  //! Discard all elements below this tolerance
  magnitude_type DropTolerance_;
  //! Condition number estimate
  magnitude_type Condest_;

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
}; // class ILUT

} // namespace Ifpack2

#endif /* IFPACK2_ILUT_HPP */
