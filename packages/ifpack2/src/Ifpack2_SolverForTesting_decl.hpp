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

#ifndef IFPACK2_SOLVER_FOR_TESTING_DECL_HPP
#define IFPACK2_SOLVER_FOR_TESTING_DECL_HPP

#include "Ifpack2_Preconditioner.hpp"

namespace Ifpack2 {

//! A class that does nothing.
template<class MatrixType>
class SolverForTesting :
    virtual public Ifpack2::Preconditioner<typename MatrixType::scalar_type,
                                           typename MatrixType::local_ordinal_type,
                                           typename MatrixType::global_ordinal_type,
                                           typename MatrixType::node_type> {
public:
  typedef typename MatrixType::scalar_type Scalar;
  typedef typename MatrixType::local_ordinal_type LocalOrdinal;
  typedef typename MatrixType::global_ordinal_type GlobalOrdinal;
  typedef typename MatrixType::node_type Node;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;

  SolverForTesting (const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A);

  //! Destructor
  virtual ~SolverForTesting();

  //! Sets parameters on this object.
  /**
    Currently doesn't need any parameters.
  */
  void setParameters (const Teuchos::ParameterList& params);

  //! Initialize
  void initialize();

  //! Returns \c true if the preconditioner has been successfully initialized.
  inline bool isInitialized() const {
    return(isInitialized_);
  }

  //! compute the preconditioner
  void compute();

  //! Return true if compute() has been called.
  inline bool isComputed() const {
    return(isComputed_);
  }

  //! @name Methods implementing a Tpetra::Operator interface.
  //@{

  //! Applies the preconditioner to X, returns the result in Y.
  /*!
    \param
    X - (In) A Tpetra::MultiVector of dimension NumVectors to be preconditioned.
    \param
    Y - (InOut) A Tpetra::MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning This routine is NOT AztecOO compliant.
  */
  void
  apply (const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
         Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
         Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const;

  //! Returns the Tpetra::Map object associated with the domain of this operator.
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
  getDomainMap () const {
    return domainMap_;
  }

  //! Returns the Tpetra::Map object associated with the range of this operator.
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
  getRangeMap() const {
    return rangeMap_;
  }

  //! Applies the matrix to a Tpetra::MultiVector.
  /*!
    \param
    X - (In) A Tpetra::MultiVector of dimension NumVectors to multiply with matrix.
    \param
    Y - (Out) A Tpetra::MultiVector of dimension NumVectors containing the result.
    */
  void applyMat(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

  //@}
  //! \name Mathematical functions.
  //@{

  //! Computes the estimated condition number and returns the value.
  magnitudeType
  computeCondEst (CondestType CT = Cheap,
                  LocalOrdinal MaxIters = 1550,
                  magnitudeType Tol = 1e-9,
                  const Teuchos::Ptr<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &matrix = Teuchos::null);

  //@}
  //! \name Attribute accessor methods
  //@{

  //! Return the computed estimated condition number, or -1.0 if no computed.
  magnitudeType getCondEst() const
  { return condEst_; }

  //! Return the communicator associated with this matrix operator.
  Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

  //! Return a reference to the matrix to be preconditioned.
  Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getMatrix () const {
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

  //@}

  private:
    bool isInitialized_;
    bool isComputed_;
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > domainMap_;
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rangeMap_;
    //Teuchos::RCP<const MatrixType> matrix_;
    const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > matrix_;

    mutable int numInitialize_;
    mutable int numCompute_;
    mutable int numApply_;

    double initializeTime_;
    double computeTime_;
    double applyTime_;

    magnitudeType condEst_;
};

}//namespace Ifpack2

#endif
