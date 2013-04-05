/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_DIAGONAL_DECL_HPP
#define IFPACK2_DIAGONAL_DECL_HPP

#include "Ifpack2_Preconditioner.hpp"

namespace Ifpack2 {

//! A class for diagonal preconditioning.
/**
Ifpack2::Diagonal wraps a vector as a diagonal preconditioner.
The preconditioner is simply defined by
\f[
z_i = D_i r_i,
\f]
where \f$r,z\f$ are the vector to be preconditioned and the preconditioned vector, and \f$D_i\f$ is the i-th element of the scaling vector.

When Ifpack2::Diagonal is constructed with a matrix, \f$D\f$ contains the inverted diagonal elements from the matrix.

When Ifpack2::Diagonal is constructed with a vector, \f$D\f$ is the caller-supplied vector.

\author Michael Heroux, ETHZ/D-INFK

\date Ifpack2 conversion (from Ifpack code) 31-Mar-2010
 */
template<class MatrixType>
class Diagonal : virtual public Ifpack2::Preconditioner<typename MatrixType::scalar_type,typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> {

public:
  typedef typename MatrixType::scalar_type Scalar;
  typedef typename MatrixType::local_ordinal_type LocalOrdinal;
  typedef typename MatrixType::global_ordinal_type GlobalOrdinal;
  typedef typename MatrixType::node_type Node;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;

  //! Constructor to create a Diagonal preconditioner using a Tpetra::CrsMatrix.
  Diagonal(const Teuchos::RCP<const MatrixType>& A);

  //! Constructor to create a Diagonal preconditioner using a Tpetra::Vector.
  /**
   * If your compiler complains about this constructor being ambigous with the
   * other constructor overload, instead call the free-standing function
   * Ifpack2::createDiagonalPreconditioner which is located at the bottom
   * of this header file.
  * (This issue arises if this constructor is called with a RCP<Tpetra::Vector>
  * that isn't const-qualified exactly as declared here.)
  */
  Diagonal(const Teuchos::RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& diag);

public:
  //! Destructor
  virtual ~Diagonal();

  //! Sets parameters on this object.
  /**
    Currently doesn't need any parameters.
  */
  void setParameters(const Teuchos::ParameterList& params);

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
  void apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
             Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
             Teuchos::ETransp mode = Teuchos::NO_TRANS,
                 Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
                 Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const;

  //! Returns the Tpetra::Map object associated with the domain of this operator.
  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& getDomainMap() const
  { return domainMap_; }

  //! Returns the Tpetra::Map object associated with the range of this operator.
  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& getRangeMap() const
  { return rangeMap_; }

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


  //@{
  //! \name Mathematical functions.

  //! Applies the preconditioner to X, returns the result in Y.
  /*! 
    \param
    X - (In) A Tpetra::MultiVector of dimension NumVectors to be preconditioned.
    \param
    Y - (InOut) A Tpetra::MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning This routine is NOT AztecOO compliant.
  */
  template <class DomainScalar, class RangeScalar>
  void applyTempl(const Tpetra::MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node>& X,
             Tpetra::MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
             Teuchos::ETransp mode = Teuchos::NO_TRANS,
                 RangeScalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
                 RangeScalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const;

  //! Computes the estimated condition number and returns the value.
  magnitudeType computeCondEst(CondestType CT = Cheap,
                               LocalOrdinal MaxIters = 1550,
                               magnitudeType Tol = 1e-9,
                               const Teuchos::Ptr<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &matrix = Teuchos::null);

  //@}

  //@{ 
  //! \name Attribute accessor methods

  //! Returns the computed estimated condition number, or -1.0 if no computed.
  magnitudeType getCondEst() const
  { return condEst_; }

  //! Returns the Tpetra::BlockMap object associated with the range of this matrix operator.
  const Teuchos::RCP<const Teuchos::Comm<int> > & getComm() const;

  //! Returns a reference to the matrix to be preconditioned.
  Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > getMatrix() const
  { return matrix_; }

  //! Returns the number of flops in the computation phase.
  double getComputeFlops() const;

  //! Returns the number of flops for the application of the preconditioner.
  double getApplyFlops() const;

  //! Returns the number of calls to initialize().
  int getNumInitialize() const;

  //! Returns the number of calls to compute().
  int getNumCompute() const;

  //! Returns the number of calls to apply().
  int getNumApply() const;

  //! Returns the time spent in initialize().
  double getInitializeTime() const;

  //! Returns the time spent in compute().
  double getComputeTime() const;

  //! Returns the time spent in apply().
  double getApplyTime() const;

  //@}

  //! @name Overridden from Teuchos::Describable 
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
    Teuchos::RCP<const MatrixType> matrix_;
    Teuchos::RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > inversediag_;

    mutable int numInitialize_;
    mutable int numCompute_;
    mutable int numApply_;

    double initializeTime_;
    double computeTime_;
    double applyTime_;

    magnitudeType condEst_;
};

/** Function to construct a Diagonal preconditioner with vector input.
* The input vector is assumed to contain the equivalent of the inverted
* diagonal of a matrix.
*
* Example usage:<br>
* typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> TCrsMatrix;<br>
* typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TVector;<br>
* typedef Tpetra::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> TPrec;
*
* Teuchos::RCP<TVector> myvec = ...
*
* Teuchos::RCP<TPrec> prec = Ifpack2::createDiagonalPreconditioner<TCrsMatrix,TVector>(myvec);
*
*/
template<class MatrixType,class VectorType>
Teuchos::RCP<Ifpack2::Diagonal<MatrixType> >
createDiagonalPreconditioner(const Teuchos::RCP<const VectorType>& invdiag)
{
  return Teuchos::rcp(new Ifpack2::Diagonal<MatrixType>(invdiag));
}

}//namespace Ifpack2

#endif
