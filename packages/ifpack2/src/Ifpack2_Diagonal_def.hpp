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

#ifndef IFPACK2_DIAGONAL_DEF_HPP
#define IFPACK2_DIAGONAL_DEF_HPP

#include "Ifpack2_Diagonal_decl.hpp"
#include "Tpetra_CrsMatrix_def.hpp"
#include "Ifpack2_Condest.hpp"

namespace Ifpack2 {

template<class MatrixType>
Diagonal<MatrixType>::Diagonal (const Teuchos::RCP<const MatrixType>& A)
 : isInitialized_ (false),
   isComputed_ (false),
   domainMap_ (A->getDomainMap ()),
   rangeMap_ (A->getRangeMap ()),
   matrix_ (A),
   numInitialize_ (0),
   numCompute_ (0),
   numApply_ (0),
   condEst_ (-Teuchos::ScalarTraits<magnitude_type>::one ())
{}

template<class MatrixType>
Diagonal<MatrixType>::
Diagonal (const Teuchos::RCP<const Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> >& diag)
 : isInitialized_ (false),
   isComputed_ (false),
   inversediag_ (diag),
   numInitialize_ (0),
   numCompute_ (0),
   numApply_ (0),
   condEst_ (-Teuchos::ScalarTraits<magnitude_type>::one ())
{}

template<class MatrixType>
Diagonal<MatrixType>::~Diagonal()
{}

template<class MatrixType>
void Diagonal<MatrixType>::setParameters(const Teuchos::ParameterList& /*params*/)
{}

template<class MatrixType>
void Diagonal<MatrixType>::initialize()
{
  // mfh 13 Dec 2013: If you call initialize(), it means that you are
  // asserting that the structure of the input sparse matrix may have
  // changed.  This means we should definitely recompute the diagonal
  // offsets.

  // if (isInitialized_) {
  //   return;
  // }

  // Precompute diagonal offsets so we don't have to search for them
  // later.  Only do this if the input matrix is nonnull and is a
  // Tpetra::CrsMatrix.
  if (matrix_.is_null ()) {
    offsets_ = Teuchos::null; // offsets are no longer valid
  }
  else { // matrix is not null
    Teuchos::RCP<const crs_matrix_type> A_crs =
      Teuchos::rcp_dynamic_cast<const crs_matrix_type> (matrix_);
    if (A_crs.is_null ()) {
      offsets_ = Teuchos::null; // offsets are no longer valid
    }
    else {
      A_crs->getLocalDiagOffsets (offsets_);
    }
  }

  isInitialized_ = true;
  ++numInitialize_;
}

template<class MatrixType>
void Diagonal<MatrixType>::compute()
{
  if (! isInitialized_) {
    initialize ();
  }
  isComputed_ = false;
  if (matrix_.is_null ()) {
    isComputed_ = true;
    return;
  }

  Teuchos::RCP<vector_type> tmp_vec (new vector_type (matrix_->getRowMap ()));
  Teuchos::RCP<const crs_matrix_type> A_crs =
    Teuchos::rcp_dynamic_cast<const crs_matrix_type> (matrix_);
  if (A_crs.is_null ()) {
    // Get the diagonal entries from the Tpetra::RowMatrix.
    matrix_->getLocalDiagCopy (*tmp_vec);
  }
  else {
    // Get the diagonal entries from the Tpetra::CrsMatrix using the
    // precomputed offsets.
    A_crs->getLocalDiagCopy (*tmp_vec, offsets_ ());
  }
  // Invert the diagonal entries.
  tmp_vec->reciprocal (*tmp_vec);
  inversediag_ = tmp_vec;

  isComputed_ = true;
  ++numCompute_;
}

template<class MatrixType>
void Diagonal<MatrixType>::
apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
       Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
       Teuchos::ETransp /*mode*/,
       scalar_type alpha,
       scalar_type beta) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isComputed(), std::runtime_error,
    "Ifpack2::Diagonal::apply() ERROR, compute() hasn't been called yet.");

  ++numApply_;
  Y.elementWiseMultiply (alpha, *inversediag_, X, beta);
}

template<class MatrixType>
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType
Diagonal<MatrixType>::
computeCondEst (CondestType CT,
                local_ordinal_type MaxIters,
                magnitude_type Tol,
                const Teuchos::Ptr<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > &matrix)
{
  const magnitude_type minusOne = -Teuchos::ScalarTraits<magnitude_type>::one ();

  if (! isComputed ()) { // cannot compute right now
    return minusOne;
  }
  // NOTE: this is computing the *local* condest
  if (condEst_ == minusOne) {
    condEst_ = Ifpack2::Condest (*this, CT, MaxIters, Tol, matrix);
  }
  return condEst_;
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
