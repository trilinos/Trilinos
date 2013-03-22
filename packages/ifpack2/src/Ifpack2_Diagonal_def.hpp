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

#ifndef IFPACK2_DIAGONAL_DEF_HPP
#define IFPACK2_DIAGONAL_DEF_HPP

#include "Ifpack2_Diagonal_decl.hpp"
#include "Ifpack2_Condest.hpp"

namespace Ifpack2 {

template<class MatrixType>
Diagonal<MatrixType>::Diagonal(const Teuchos::RCP<const MatrixType>& A)
 : isInitialized_(false),
   isComputed_(false),
   domainMap_(A->getDomainMap()),
   rangeMap_(A->getRangeMap()),
   matrix_(A),
   inversediag_(),
   numInitialize_(0),
   numCompute_(0),
   numApply_(0),
   condEst_ (-Teuchos::ScalarTraits<magnitudeType>::one ())
{
}

template<class MatrixType>
Diagonal<MatrixType>::Diagonal(const Teuchos::RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& diag)
 : isInitialized_(false),
   isComputed_(false),
   domainMap_(),
   rangeMap_(),
   matrix_(),
   inversediag_(diag),
   numInitialize_(0),
   numCompute_(0),
   numApply_(0),
   condEst_ (-Teuchos::ScalarTraits<magnitudeType>::one ())
{
}

template<class MatrixType>
Diagonal<MatrixType>::~Diagonal()
{
}

template<class MatrixType>
void Diagonal<MatrixType>::setParameters(const Teuchos::ParameterList& /*params*/)
{
}

template<class MatrixType>
void Diagonal<MatrixType>::initialize()
{
  if (isInitialized_ == true) return;
  isInitialized_ = true;
  ++numInitialize_;
  //nothing to do
}

template<class MatrixType>
void Diagonal<MatrixType>::compute()
{
  initialize();
  ++numCompute_;

  if (isComputed_ == true) return;

  isComputed_ = true;

  if (matrix_ == Teuchos::null) return;

  Teuchos::RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tmp_vec = Teuchos::rcp(new Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(matrix_->getRowMap()));

  matrix_->getLocalDiagCopy(*tmp_vec);
  tmp_vec->reciprocal(*tmp_vec);

  inversediag_ = tmp_vec;
}

template<class MatrixType>
void Diagonal<MatrixType>::apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
             Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
             Teuchos::ETransp /*mode*/,
                 Scalar alpha,
                 Scalar beta) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isComputed(), std::runtime_error,
    "Ifpack2::Diagonal::apply() ERROR, compute() hasn't been called yet.");

  ++numApply_;
  Y.elementWiseMultiply(alpha, *inversediag_, X, beta);
}

template<class MatrixType>
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType
Diagonal<MatrixType>::computeCondEst(
                     CondestType CT,
                     LocalOrdinal MaxIters,
                     magnitudeType Tol,
                     const Teuchos::Ptr<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &matrix)
{
  const magnitudeType minusOne = Teuchos::ScalarTraits<magnitudeType>::one ();

  if (!isComputed()) { // cannot compute right now
    return minusOne;
  }
  // NOTE: this is computing the *local* condest
  if (condEst_ == minusOne) {
    condEst_ = Ifpack2::Condest(*this, CT, MaxIters, Tol, matrix);
  }
  return(condEst_);
}

template <class MatrixType>
int Diagonal<MatrixType>::getNumInitialize() const {
  return(numInitialize_);
}

template <class MatrixType>
int Diagonal<MatrixType>::getNumCompute() const {
  return(numCompute_);
}

template <class MatrixType>
int Diagonal<MatrixType>::getNumApply() const {
  return(numApply_);
}

template <class MatrixType>
double Diagonal<MatrixType>::getInitializeTime() const {
  return(initializeTime_);
}

template<class MatrixType>
double Diagonal<MatrixType>::getComputeTime() const {
  return(computeTime_);
}

template<class MatrixType>
double Diagonal<MatrixType>::getApplyTime() const {
  return(applyTime_);
}

template <class MatrixType>
std::string Diagonal<MatrixType>::description() const
{
  return std::string("Ifpack2::Diagonal");
}

template <class MatrixType>
void Diagonal<MatrixType>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const
{
  if (verbLevel != Teuchos::VERB_NONE) {
    out << this->description() << std::endl;
    out << "  numApply: " << numApply_ << std::endl;
  }
}

}//namespace Ifpack2

#endif
