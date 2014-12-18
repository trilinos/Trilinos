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

#ifndef IFPACK2_IDENTITY_SOLVER_DECL_HPP
#define IFPACK2_IDENTITY_SOLVER_DECL_HPP

#include <Ifpack2_Preconditioner.hpp>
#include <Ifpack2_Details_CanChangeMatrix.hpp>


namespace Ifpack2 {

/// \brief "Identity" preconditioner.
/// \tparam MatrixType Specialization of Tpetra::RowMatrix.
///
/// This class is mainly useful for testing.  It implements the
/// identity operator; its apply() method just copies the input to the
/// output.
template<class MatrixType>
class IdentitySolver :
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
  typedef typename MatrixType::scalar_type                                scalar_type;
  //! Type of the local indices of the input matrix.
  typedef typename MatrixType::local_ordinal_type                         local_ordinal_type;
  //! Type of the global indices of the input matrix.
  typedef typename MatrixType::global_ordinal_type                        global_ordinal_type;
  //! Node type of the input matrix.
  typedef typename MatrixType::node_type                                  node_type;

  //! Type of the absolute value (magnitude) of a \c scalar_type value.
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType      magnitude_type;
  //! Specialization of Tpetra::Map used by this class.
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;
  //! Specialization of Tpetra::RowMatrix used by this class.
  typedef Tpetra::RowMatrix<scalar_type, local_ordinal_type,
                            global_ordinal_type, node_type> row_matrix_type;

  //! Constructor: Takes the matrix to precondition.
  IdentitySolver (const Teuchos::RCP<const row_matrix_type>& A);

  //! Destructor.
  virtual ~IdentitySolver();

  /// \brief Set this object's parameters.
  ///
  /// This object does not currently take any parameters.
  void setParameters (const Teuchos::ParameterList& params);

  //! Initialize
  void initialize();

  //! Return \c true if the preconditioner has been successfully initialized.
  inline bool isInitialized() const {
    return(isInitialized_);
  }

  //! Compute the preconditioner
  void compute();

  //! Return true if compute() has been called.
  inline bool isComputed() const {
    return(isComputed_);
  }

  //! @name Implementation of Tpetra::Operator
  //@{

  /// \brief Apply the preconditioner to X, and put the result in Y.
  ///
  /// \param X [in] MultiVector to which to apply the preconditioner.
  ///
  /// \param Y [in/out] On input: Initial guess, if applicable.
  ///   On poutput: Result of applying the preconditioner.
  void
  apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
         Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
         scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //! Return the Tpetra::Map object associated with the domain of this operator.
  Teuchos::RCP<const map_type> getDomainMap () const;

  //! Return the Tpetra::Map object associated with the range of this operator.
  Teuchos::RCP<const map_type> getRangeMap () const;

  /// \brief Apply the original input matrix.
  ///
  /// \param X [in] MultiVector input.
  /// \param Y [in/out] Result of applying the matrix A to X.
  void
  applyMat (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
            Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
            Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

  //@}
  //! \name Mathematical functions.
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
                  const Teuchos::Ptr<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > &matrix = Teuchos::null);

  //@}
  //! \name Attribute accessor methods
  //@{

  /// \brief Return the computed condition number estimate, or -1 if not computed.
  ///
  /// \warning This method is DEPRECATED.  See warning for computeCondEst().
  virtual magnitude_type TEUCHOS_DEPRECATED getCondEst() const {
    return condEst_;
  }

  //! Return the communicator associated with this matrix operator.
  Teuchos::RCP<const Teuchos::Comm<int> > getComm () const;

  //! Return a reference to the matrix to be preconditioned.
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

  /** \brief Return a simple one-line description of this object. */
  std::string description() const;

  /** \brief Print the object with some verbosity level to an FancyOStream object. */
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

  /// \brief Set this preconditioner's matrix.
  ///
  /// After calling this method, you must call first initialize(),
  /// then compute(), before you may call apply().
  virtual void setMatrix (const Teuchos::RCP<const row_matrix_type>& A);

  //@}

private:
  //! Specialization of Tpetra::Export used by this preconditioner.
  typedef Tpetra::Export<local_ordinal_type, global_ordinal_type, node_type> export_type;

  Teuchos::RCP<const row_matrix_type> matrix_;

  /// \brief Cached communication pattern from domain to range Map of \c matrix_.
  ///
  /// This is only nonnull, and is only used, when the domain and
  /// range Maps are not the same.  They must be compatible, but they
  /// need not necessarily be the same.  We need this because if the
  /// Maps are compatible but not the same, just copying X to Y in
  /// apply() would be a permutation, not the identity operation.
  Teuchos::RCP<const export_type> export_;

  bool isInitialized_;
  bool isComputed_;

  mutable int numInitialize_;
  mutable int numCompute_;
  mutable int numApply_;

  double initializeTime_;
  double computeTime_;
  double applyTime_;

  magnitude_type condEst_;
};

}//namespace Ifpack2

#endif
