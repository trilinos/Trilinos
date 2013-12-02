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

#ifndef IFPACK2_AMESOS2SOLVER_DEF_HPP
#define IFPACK2_AMESOS2SOLVER_DEF_HPP

// disable clang warnings
#ifdef __clang__
#pragma clang system_header
#endif

#include <Ifpack2_Heap.hpp>
#include <Ifpack2_Condest.hpp>
#include <Ifpack2_Parameters.hpp>
#include <Amesos2.hpp>

#include <Teuchos_Time.hpp>
#include <Teuchos_TypeNameTraits.hpp>


namespace Ifpack2 {

template <class MatrixType>
Amesos2solver<MatrixType>::Amesos2solver(const Teuchos::RCP<const row_matrix_type>& A) :
  A_ (A),
  Condest_ (-Teuchos::ScalarTraits<magnitude_type>::one ()),
  InitializeTime_ (0.0),
  ComputeTime_ (0.0),
  ApplyTime_ (0.0),
  NumInitialize_ (0),
  NumCompute_ (0),
  NumApply_ (0),
  IsInitialized_ (false),
  IsComputed_ (false)
{}

template <class MatrixType>
Amesos2solver<MatrixType>::~Amesos2solver()
{}

template <class MatrixType>
void Amesos2solver<MatrixType>::setParameters (const Teuchos::ParameterList& params)
{
  // no parameters set right now
}


template <class MatrixType>
Teuchos::RCP<const Teuchos::Comm<int> >
Amesos2solver<MatrixType>::getComm () const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Amesos2solver::getComm: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");
  return A_->getComm ();
}


template <class MatrixType>
Teuchos::RCP<const typename Amesos2solver<MatrixType>::row_matrix_type>
Amesos2solver<MatrixType>::getMatrix () const {
  return A_;
}


template <class MatrixType>
Teuchos::RCP<const typename Amesos2solver<MatrixType>::map_type>
Amesos2solver<MatrixType>::getDomainMap () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Amesos2solver::getDomainMap: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");
  return A_->getDomainMap ();
}


template <class MatrixType>
Teuchos::RCP<const typename Amesos2solver<MatrixType>::map_type>
Amesos2solver<MatrixType>::getRangeMap () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Amesos2solver::getRangeMap: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");
  return A_->getRangeMap ();
}


template <class MatrixType>
bool Amesos2solver<MatrixType>::hasTransposeApply () const {
  return true;
}


template <class MatrixType>
int Amesos2solver<MatrixType>::getNumInitialize () const {
  return NumInitialize_;
}


template <class MatrixType>
int Amesos2solver<MatrixType>::getNumCompute () const {
  return NumCompute_;
}


template <class MatrixType>
int Amesos2solver<MatrixType>::getNumApply () const {
  return NumApply_;
}


template <class MatrixType>
double Amesos2solver<MatrixType>::getInitializeTime () const {
  return InitializeTime_;
}


template<class MatrixType>
double Amesos2solver<MatrixType>::getComputeTime () const {
  return ComputeTime_;
}


template<class MatrixType>
double Amesos2solver<MatrixType>::getApplyTime () const {
  return ApplyTime_;
}

template<class MatrixType>
typename Amesos2solver<MatrixType>::magnitude_type
Amesos2solver<MatrixType>::
computeCondEst (CondestType CT,
                local_ordinal_type MaxIters,
                magnitude_type Tol,
                const Teuchos::Ptr<const row_matrix_type>& matrix)
{
  if (! isComputed ()) {
    return -STM::one ();
  }
  // NOTE: this is computing the *local* condest
  if (Condest_ == -STM::one ()) {
    Condest_ = Ifpack2::Condest(*this, CT, MaxIters, Tol, matrix);
  }
  return Condest_;
}


template<class MatrixType>
void Amesos2solver<MatrixType>::setMatrix (const Teuchos::RCP<const row_matrix_type>& A)
{
  // Check in serial or one-process mode if the matrix is square.
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! A.is_null () && A->getComm ()->getSize () == 1 &&
    A->getNodeNumRows () != A->getNodeNumCols (),
    std::runtime_error, "Ifpack2::Amesos2solver::setMatrix: If A's communicator only "
    "contains one process, then A must be square.  Instead, you provided a "
    "matrix A with " << A->getNodeNumRows () << " rows and "
    << A->getNodeNumCols () << " columns.");

  // It's legal for A to be null; in that case, you may not call
  // initialize() until calling setMatrix() with a nonnull input.
  // Regardless, setting the matrix invalidates any previous
  // factorization.
  IsInitialized_ = false;
  IsComputed_ = false;
  A_local_ = Teuchos::null;
  A_ = A;
}


template<class MatrixType>
void Amesos2solver<MatrixType>::initialize ()
{
  Teuchos::Time timer ("Amesos2solver::initialize");
  {
    Teuchos::TimeMonitor timeMon (timer);

    // Check that the matrix is nonnull.
    TEUCHOS_TEST_FOR_EXCEPTION(
      A_.is_null (), std::runtime_error, "Ifpack2::Amesos2solver::initialize: "
      "The matrix to precondition is null.  Please call setMatrix() with a "
      "nonnull input before calling this method.");

    // Clear any previous computations.
    IsInitialized_ = false;
    IsComputed_ = false;
    A_local_ = Teuchos::null;
    A_local_ = makeLocalMatrix(A_); // Construct the local matrix.

    std::string solver_type;

#if defined(HAVE_AMESOS2_SUPERLU)
    solver_type = "superlu";
#elif defined(HAVE_AMESOS2_KLU2)
    solver_type = "klu";
#elif defined(HAVE_AMESOS2_SUPERLUDIST)
    solver_type = "superludist";
#elif defined(HAVE_AMESOS2_LAPACK)
    solver_type = "lapack";
#endif
    
    amesos2solver_ = Amesos2::create<MatrixType,multivector_type>(solver_type, A_local_);
    amesos2solver_ -> preOrdering();

    IsInitialized_ = true;
    ++NumInitialize_;
  }
  InitializeTime_ += timer.totalElapsedTime ();
}

template<class MatrixType>
void Amesos2solver<MatrixType>::compute ()
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::rcp;
  using Teuchos::reduceAll;

  // Don't count initialization in the compute() time.
  if (! isInitialized ()) {
    initialize ();
  }

  Teuchos::Time timer ("Amesos2solver::compute");
  {
    // Timer scope for timing compute()
    Teuchos::TimeMonitor timeMon (timer, true);
    amesos2solver_ -> symbolicFactorization();
    amesos2solver_ -> numericFactorization();
  }
  ComputeTime_ += timer.totalElapsedTime ();
  IsComputed_ = true;
  ++NumCompute_;
}


template <class MatrixType>
void Amesos2solver<MatrixType>::
apply (const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
       Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> MV;

  Teuchos::Time timer ("Amesos2solver::apply");
  { 
    // Timer scope for timing apply()
    Teuchos::TimeMonitor timeMon (timer, true);

    TEUCHOS_TEST_FOR_EXCEPTION(
      X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
      "Ifpack2::Amesos2solver::apply: X and Y must have the same number of columns.  "
      "X has " << X.getNumVectors () << " columns, but Y has "
      << Y.getNumVectors () << " columns.");

    // If beta != 0, create a temporary multivector Y_temp to hold the
    // contents of alpha*M^{-1}*X.  Otherwise, alias Y_temp to Y.
    RCP<MV> Y_temp;
    Y_temp = rcp (new MV (Y.getMap (), Y.getNumVectors ()));

    // If X and Y are pointing to the same memory location, create an
    // auxiliary vector, X_temp, so that we don't clobber the input
    // when computing the output.  Otherwise, alias X_temp to X.
    RCP<const MV> X_temp;
    if (X.getLocalMV ().getValues () == Y.getLocalMV ().getValues ()) {
      X_temp = rcp (new MV (X));
    } else {
      X_temp = rcpFromRef (X);
    }

    // construct local vectors    
    size_t numvecs = X_temp->getNumVectors();
    Teuchos::RCP<const map_type> globalRowMap = A_->getRowMap();
    Teuchos::ArrayView<const global_ordinal_type> globalRows = globalRowMap -> getNodeElementList();
    local_ordinal_type numRows = globalRowMap->getNodeNumElements();
    Teuchos::RCP<MV> localY = Teuchos::rcp( new MV(localRowMap_,numvecs) );
    Teuchos::RCP<MV> localX = Teuchos::rcp( new MV(localRowMap_,numvecs) );
    // extract values
    for(size_t j=0; j<numvecs; j++) {
      Teuchos::ArrayRCP<const scalar_type> vecj = X_temp -> getData(j);
      for(local_ordinal_type i = 0; i < numRows; i++) {
	localX->replaceLocalValue(i,j,vecj[i]);
      }
    }

    // solve
    amesos2solver_->setX(localY);
    amesos2solver_->setB(localX);
    amesos2solver_->solve();

    // extract to global vector
    for(size_t j=0; j<localY->getNumVectors(); j++) {
      Teuchos::ArrayRCP<const scalar_type> localview = localY->getData(j);
      for(unsigned int i=0; i<globalRowMap->getNodeNumElements(); i++) {
	Y_temp->replaceLocalValue(i,j,localview[i]);
      }
    }

    // update Y
    Y.update(alpha, *Y_temp, beta);

  }
  ++NumApply_;
  ApplyTime_ += timer.totalElapsedTime ();
}


template <class MatrixType>
std::string Amesos2solver<MatrixType>::description() const {
  using Teuchos::TypeNameTraits;
  std::ostringstream os;

  os << "Ifpack2::Amesos2solver<" << TypeNameTraits<MatrixType>::name ()
     << ">: {";
  if (this->getObjectLabel () != "") {
    os << "Label: \"" << this->getObjectLabel () << "\", ";
  }
  os << "Initialized: " << (isInitialized () ? "true" : "false")
     << ", "
     << "Computed: " << (isComputed () ? "true" : "false")
     << ", "
     << "Number of rows: " << A_->getGlobalNumRows ()
     << ", "
     << "Number of columns: " << A_->getGlobalNumCols ()
     << "}";
  return os.str();
}


template <class MatrixType>
void
Amesos2solver<MatrixType>::
describe (Teuchos::FancyOStream& out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using Teuchos::Comm;
  using Teuchos::OSTab;
  using Teuchos::RCP;
  using Teuchos::TypeNameTraits;
  using std::endl;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;

  const Teuchos::EVerbosityLevel vl = (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;
  OSTab tab0 (out);

  if (vl > VERB_NONE) {

    out << "Ifpack2::Amesos2solver:" << endl;
    OSTab tab1 (out);
    out << "MatrixType: " << TypeNameTraits<MatrixType>::name () << endl;

    if (this->getObjectLabel () != "") {
      out << "Label: \"" << this->getObjectLabel () << "\"" << endl;
    }

    out << "Initialized: " << (isInitialized () ? "true" : "false") << endl;
    out << "Computed: " << (isComputed () ? "true" : "false") << endl;

    out << "Number of initialize calls: " << getNumInitialize () << endl;
    out << "Number of compute calls: " << getNumCompute () << endl;
    out << "Number of apply calls: " << getNumApply () << endl;
    out << "Total time in seconds for initialize: " << getInitializeTime () << endl;
    out << "Total time in seconds for compute: " << getComputeTime () << endl;
    out << "Total time in seconds for apply: " << getApplyTime () << endl;
    out << "Local matrix:" << endl;
    A_local_->describe (out, vl);
  }
}

template <class MatrixType>
Teuchos::RCP<MatrixType>
Amesos2solver<MatrixType>::makeLocalMatrix(Teuchos::RCP<const row_matrix_type> A)
{
  // local communicator
  Teuchos::RCP<const Teuchos::Comm<int> > localComm;
#ifdef HAVE_MPI
  localComm = Teuchos::rcp(new Teuchos::MpiComm<int> (MPI_COMM_SELF));
#else
  localComm = Teuchos::rcp(new Teuchos::SerialComm<int> ());
#endif

  // get row map and setup local matrix
  Teuchos::RCP<const map_type> globalRowMap = A->getRowMap();
  Teuchos::RCP<const map_type> globalColMap = A->getColMap();
  local_ordinal_type numRows = globalRowMap -> getNodeNumElements();
  localRowMap_ = Teuchos::rcp(new map_type(numRows, 0, localComm,
					   Tpetra::GloballyDistributed, A->getNode()));
  Teuchos::RCP<MatrixType> Alocal = Teuchos::rcp(new MatrixType(localRowMap_,localRowMap_,100));

  // extract rows
  for(local_ordinal_type i = 0; i < numRows; i++) {
    Teuchos::ArrayView<const local_ordinal_type> indices;
    Teuchos::ArrayView<const scalar_type> values;
    A -> getLocalRowView(i,indices,values);
    std::vector<local_ordinal_type> indices_vec;
    std::vector<scalar_type> values_vec;
    indices_vec.resize(0);
    values_vec.resize(0);
    for(unsigned int j=0; j < indices.size(); j++ ) {
      local_ordinal_type local_col = indices[j];
      global_ordinal_type global_col = globalColMap->getGlobalElement(local_col);
      if(globalRowMap->isNodeGlobalElement(global_col)) {
	local_col = globalRowMap->getLocalElement(global_col);
	indices_vec.push_back(local_col);
	values_vec.push_back(values[j]);
      }
    }
    Alocal->insertLocalValues(i,
			      Teuchos::ArrayView<local_ordinal_type>(indices_vec),
			      Teuchos::ArrayView<scalar_type>(values_vec));
  }
  Alocal->fillComplete();
  return Alocal;
}

}//namespace Ifpack2

#endif /* IFPACK2_AMESOS2SOLVER_DEF_HPP */

