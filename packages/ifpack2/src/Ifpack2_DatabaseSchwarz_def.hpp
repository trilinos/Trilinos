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

#ifndef IFPACK2_DATABASESCHWARZ_DEF_HPP
#define IFPACK2_DATABASESCHWARZ_DEF_HPP

#include "Ifpack2_Parameters.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_LAPACK.hpp"
#include <iostream>
#include <sstream>


namespace Ifpack2 {

template<class MatrixType>
DatabaseSchwarz<MatrixType>::
DatabaseSchwarz (const Teuchos::RCP<const row_matrix_type>& A)
  : A_(A),
    IsInitialized_(false),
    IsComputed_(false),
    NumInitialize_(0),
    NumCompute_(0),
    NumApply_(0),
    InitializeTime_(0.0),
    ComputeTime_(0.0),
    ApplyTime_(0.0),
    ComputeFlops_(0.0),
    ApplyFlops_(0.0)
{
  this->setObjectLabel("Ifpack2::DatabaseSchwarz");
}


template<class MatrixType>
DatabaseSchwarz<MatrixType>::~DatabaseSchwarz() {
}


template<class MatrixType>
void DatabaseSchwarz<MatrixType>::setMatrix(const Teuchos::RCP<const row_matrix_type>& A)
{
  // ASSERT NON-NULL INPUT
  if (A.getRawPtr() != A_.getRawPtr()) {
    IsInitialized_ = false;
    IsComputed_ = false;
    A_ = A;
  }
}

template<class MatrixType>
void
DatabaseSchwarz<MatrixType>::setParameters (const Teuchos::ParameterList& List)
{
  (void) List;
}


template<class MatrixType>
void
DatabaseSchwarz<MatrixType>::setZeroStartingSolution (bool zeroStartingSolution)
{
  (void) zeroStartingSolution;
}

template<class MatrixType>
Teuchos::RCP<const Teuchos::Comm<int> >
DatabaseSchwarz<MatrixType>::getComm() const
{
  Teuchos::RCP<const row_matrix_type> A = getMatrix();
  TEUCHOS_TEST_FOR_EXCEPTION(
    A.is_null(), std::runtime_error, "Ifpack2::DatabaseSchwarz::getComm: The input "
    "matrix A is null.  Please call setMatrix() with a nonnull input matrix "
    "before calling this method.");
  return A->getRowMap()->getComm();
}


template<class MatrixType>
Teuchos::RCP<const typename DatabaseSchwarz<MatrixType>::row_matrix_type>
DatabaseSchwarz<MatrixType>::getMatrix() const
{
  return getMatrix();
}


template<class MatrixType>
Teuchos::RCP<const Tpetra::CrsMatrix<typename MatrixType::scalar_type,
                                     typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >
DatabaseSchwarz<MatrixType>::
getCrsMatrix() const {
  typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type,
    global_ordinal_type, node_type> crs_matrix_type;
  return Teuchos::rcp_dynamic_cast<const crs_matrix_type> (getMatrix());
}


template<class MatrixType>
Teuchos::RCP<const typename DatabaseSchwarz<MatrixType>::map_type>
DatabaseSchwarz<MatrixType>::
getDomainMap() const
{
  Teuchos::RCP<const row_matrix_type> A = getMatrix();
  TEUCHOS_TEST_FOR_EXCEPTION(
    A.is_null(), std::runtime_error, "Ifpack2::DatabaseSchwarz::getDomainMap: The "
    "input matrix A is null.  Please call setMatrix() with a nonnull input "
    "matrix before calling this method.");
  return A->getDomainMap();
}


template<class MatrixType>
Teuchos::RCP<const typename DatabaseSchwarz<MatrixType>::map_type>
DatabaseSchwarz<MatrixType>::
getRangeMap() const
{
  Teuchos::RCP<const row_matrix_type> A = getMatrix();
  TEUCHOS_TEST_FOR_EXCEPTION(
    A.is_null(), std::runtime_error, "Ifpack2::DatabaseSchwarz::getRangeMap: The "
    "input matrix A is null.  Please call setMatrix() with a nonnull input "
    "matrix before calling this method.");
  return A->getRangeMap();
}


template<class MatrixType>
bool DatabaseSchwarz<MatrixType>::hasTransposeApply() const {
  return true; // GH: redo
}


template<class MatrixType>
int DatabaseSchwarz<MatrixType>::getNumInitialize() const {
  return NumInitialize_;
}


template<class MatrixType>
int DatabaseSchwarz<MatrixType>::getNumCompute() const {
  return NumCompute_;
}


template<class MatrixType>
int DatabaseSchwarz<MatrixType>::getNumApply() const {
  return NumApply_;
}


template<class MatrixType>
double DatabaseSchwarz<MatrixType>::getInitializeTime() const {
  return InitializeTime_;
}


template<class MatrixType>
double DatabaseSchwarz<MatrixType>::getComputeTime() const {
  return ComputeTime_;
}


template<class MatrixType>
double DatabaseSchwarz<MatrixType>::getApplyTime() const {
  return ApplyTime_;
}


template<class MatrixType>
double DatabaseSchwarz<MatrixType>::getComputeFlops() const {
  return ComputeFlops_;
}


template<class MatrixType>
double DatabaseSchwarz<MatrixType>::getApplyFlops() const {
  return ApplyFlops_;
}

template<class MatrixType>
size_t DatabaseSchwarz<MatrixType>::getNodeSmootherComplexity() const {
  Teuchos::RCP<const row_matrix_type> A = getMatrix();
  TEUCHOS_TEST_FOR_EXCEPTION(
    A.is_null(), std::runtime_error, "Ifpack2::DatabaseSchwarz::getNodeSmootherComplexity: "
    "The input matrix A is null.  Please call setMatrix() with a nonnull "
    "input matrix, then call compute(), before calling this method.");
  // DatabaseSchwarz costs roughly one apply + one diagonal inverse per iteration
  return A->getLocalNumRows() + A->getLocalNumEntries();
}



template<class MatrixType>
void
DatabaseSchwarz<MatrixType>::
apply(const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
      Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& Y,
      Teuchos::ETransp mode,
      scalar_type alpha,
      scalar_type beta) const
{
  const std::string timerName ("Ifpack2::DatabaseSchwarz::apply");
  Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::lookupCounter (timerName);
  if (timer.is_null()) {
    timer = Teuchos::TimeMonitor::getNewCounter (timerName);
  }

  double startTime = timer->wallTime();

  // Start timing here.
  {
    Teuchos::TimeMonitor timeMon (*timer);

    // compute() calls initialize() if it hasn't already been called.
    // Thus, we only need to check isComputed().
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! isComputed(), std::runtime_error,
      "Ifpack2::DatabaseSchwarz::apply(): You must call the compute() method before "
      "you may call apply().");
    TEUCHOS_TEST_FOR_EXCEPTION(
      X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
      "Ifpack2::DatabaseSchwarz::apply(): X and Y must have the same number of "
      "columns.  X.getNumVectors() = " << X.getNumVectors() << " != "
      << "Y.getNumVectors() = " << Y.getNumVectors() << ".");

    // 1. Compute beta*Y
    Y.scale(beta);

    // 2. Solve prec on X
    // 2a. Split X into Xk on each local patch

    // 2b. Solve each Xk using Lapack::GETRS(const char & TRANS, const OrdinalType & n, const OrdinalType & nrhs, const ScalarType * A, const OrdinalType & lda, const OrdinalType * IPIV, ScalarType * B, const OrdinalType & ldb, OrdinalType * info)
    
    // 2c. Add alpha*weights*Xk into Y for each Xk
    
  }
  ++NumApply_;
  ApplyTime_ += (timer->wallTime() - startTime);
}


template<class MatrixType>
void
DatabaseSchwarz<MatrixType>::
applyMat (const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
          Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& Y,
          Teuchos::ETransp mode) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.getNumVectors() != Y.getNumVectors(), std::invalid_argument,
    "Ifpack2::DatabaseSchwarz::applyMat: X.getNumVectors() != Y.getNumVectors().");

  Teuchos::RCP<const row_matrix_type> A = getMatrix();
  TEUCHOS_TEST_FOR_EXCEPTION(
    A.is_null(), std::runtime_error, "Ifpack2::DatabaseSchwarz::applyMat: The input "
    "matrix A is null.  Please call setMatrix() with a nonnull input matrix "
    "before calling this method.");

  A->apply (X, Y, mode);
}


template<class MatrixType>
void DatabaseSchwarz<MatrixType>::initialize() {
  // We create the timer, but this method doesn't do anything, so
  // there is no need to start the timer.  The resulting total time
  // will always be zero.
  const std::string timerName ("Ifpack2::DatabaseSchwarz::initialize");
  Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::lookupCounter (timerName);
  if (timer.is_null()) {
    timer = Teuchos::TimeMonitor::getNewCounter (timerName);
  }
  IsInitialized_ = true;
  ++NumInitialize_;
}


template<class MatrixType>
void DatabaseSchwarz<MatrixType>::compute()
{
  const std::string timerName ("Ifpack2::DatabaseSchwarz::compute");
  Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::lookupCounter (timerName);
  if (timer.is_null()) {
    timer = Teuchos::TimeMonitor::getNewCounter (timerName);
  }

  double startTime = timer->wallTime();

  // Start timing here.
  {
    Teuchos::TimeMonitor timeMon(*timer);
    if (!isInitialized()) {
      initialize();
    }
    IsComputed_ = false;
    

    // we do not change the local row, so we avoid doing a copy
    // for(unsigned int i=0; i<A_->getLocalNumRows(); ++i) {
    //   typename row_matrix_type::local_inds_host_view_type indices;
    //   typename row_matrix_type::values_host_view_type values;
    //   size_t numEntries = 0;
    //   A_->getLocalRowView(i, indices, values); // LocalRowView complains about temporaries and nonconst lvalues, we'll revisit
    //   std::cout << indices.extent(0) << std::endl;
    //   std::cout << values.extent(0) << std::endl;
    //   std::cout << numEntries << std::endl;
    // }

    // Phase 1: Loop over rows of A_, construct patch indices, and construct patch matrices
    std::vector<std::vector<typename row_matrix_type::local_ordinal_type> > patchIndices;
    NumPatches_ = 0;
    // loop over potential candidates by checking rows of A_
    for(unsigned int irow=0; irow<A_->getLocalNumRows(); ++irow) {
      // grab row irow of A_
      typename row_matrix_type::nonconst_local_inds_host_view_type rowIndices("row indices", 10);
      typename row_matrix_type::nonconst_values_host_view_type rowValues("row values", 10);
      size_t numEntries;
      A_->getLocalRowCopy(irow, rowIndices, rowValues, numEntries);

      // if irow is a potential patch candidate
      if(numEntries == PatchSize_) {
        // check if we've added DOF irow before
        bool isNewPatch = true;
        for(size_t ipatch=0; ipatch<patchIndices.size(); ++ipatch) {
          for(size_t i=0; i<patchIndices[ipatch].size(); ++i) {
            if(patchIndices[ipatch][i]==irow) {
              isNewPatch = false;
              // likely the ugliest way to break out other than using goto TheEnd:
              ipatch=patchIndices.size();
              break;
            }
          }
        }

        // if this patch is new, append the indices and then compute the local Ak
        if(isNewPatch) {
          // append indices
          Teuchos::ArrayView<typename row_matrix_type::local_ordinal_type> indicesArrayView(rowIndices.data(),numEntries);
          std::vector<typename row_matrix_type::local_ordinal_type> indicesVector = Teuchos::createVector(indicesArrayView);
          patchIndices.push_back(indicesVector);
          NumPatches_++;
        }
      }
    }
    // sanity check everything
    std::cout << "Found " << NumPatches_ << " patches!" << std::endl;
    for(size_t ipatch=0; ipatch<patchIndices.size(); ++ipatch) {
      std::cout << "patch " << ipatch << " = [";
      for(size_t i=0; i<patchIndices[ipatch].size(); ++i)
        std::cout << patchIndices[ipatch][i] << " ";
      std::cout << "]" << std::endl;
    }

    // Phase 2: construct the list of local patch matrices
    // std::vector<Teuchos::SerialDenseMatrix<row_matrix_type::local_ordinal_type,row_matrix_type::scalar_type> > patchMatrices;
    // // compute the local patch matrix A_k via Vk*A_*Vk^T
    // std::vector<typename row_matrix_type::scalar_type> Ak_values(PatchSize_*PatchSize_,0);
    // for(size_t i=0; i<PatchSize_; ++i) {
    //   A_->getLocalRowCopy(indicesVector[i], rowIndices, rowValues, numEntries);
    //   for(size_t j=0; j<PatchSize_; ++j) {
    //     Ak_values.push_back(rowValues(j));
    //   }
    // }
    // Teuchos::SerialDenseMatrix<int,double> patchMatrix(Teuchos::DataAccess::Copy, Ak_values.data(), 1, PatchSize_, PatchSize_);

    // // check if the local patch matrix has been seen before
    // // ignore for now

    // //Teuchos::ArrayView<double> densematvalues[PatchSize_*PatchSize_];
    // //Teuchos::SerialDenseMatrix<int,double> patchmatrix(PatchSize_, PatchSize_);
    // //Teuchos::SerialDenseMatrix<int,double> patchmatrix(Teuchos::DataAccess::Copy, densematvalues, 1, PatchSize_, PatchSize_);
    // //Teuchos::RCP<Teuchos::SerialDenseMatrix<int,double> > patchmatrix = Teuchos::rcp(new Teuchos::SerialDenseMatrix<int,double>(Teuchos::DataAccess::Copy, densematvalues, 1, PatchSize_, PatchSize_));

    // // add the matrix to the database
    // patchMatrices.push_back(patchMatrix);

    // Phase 3: factor the patches using LAPACK (GETRF for factorization)
    Teuchos::LAPACK<int, double> lapack;
    int INFO = 0;
    // int* blockIpiv = &ipiv_[this->blockOffsets_[i] * this->scalarsPerRow_];
    // lapack.GETRF(diagBlocks_[i].numRows(),
    //             diagBlocks_[i].numCols(),
    //             diagBlocks_[i].values(),
    //             diagBlocks_[i].stride(),
    //             blockIpiv, &INFO);
    // INFO < 0 is a bug.
    TEUCHOS_TEST_FOR_EXCEPTION(
      INFO < 0, std::logic_error, "Ifpack2::DenseContainer::factor: "
      "LAPACK's _GETRF (LU factorization with partial pivoting) was called "
      "incorrectly.  INFO = " << INFO << " < 0.  "
      "Please report this bug to the Ifpack2 developers.");
    // INFO > 0 means the matrix is singular.  This is probably an issue
    // either with the choice of rows the rows we extracted, or with the
    // input matrix itself.
    TEUCHOS_TEST_FOR_EXCEPTION(
      INFO > 0, std::runtime_error, "Ifpack2::DenseContainer::factor: "
      "LAPACK's _GETRF (LU factorization with partial pivoting) reports that the "
      "computed U factor is exactly singular.  U(" << INFO << "," << INFO << ") "
      "(one-based index i) is exactly zero.  This probably means that the input "
      "patch is singular.");


  }
  IsComputed_ = true;
  ++NumCompute_;

  ComputeTime_ += (timer->wallTime() - startTime);
}


template <class MatrixType>
std::string DatabaseSchwarz<MatrixType>::description() const {
  std::ostringstream out;
  std::cout << "DatabaseSchwarz description 1..." << std::endl;
  // Output is a valid YAML dictionary in flow style.  If you don't
  // like everything on a single line, you should call describe()
  // instead.
  out << "\"Ifpack2::DatabaseSchwarz\": {";
  out << "Initialized: " << (isInitialized() ? "true" : "false") << ", "
      << "Computed: " << (isComputed() ? "true" : "false") << ", ";
  std::cout << "DatabaseSchwarz description 2..." << std::endl;
  if (getMatrix().is_null()) {
    out << "Matrix: null";
  }
  else {
    out << "Global matrix dimensions: ["
        << getMatrix()->getGlobalNumRows() << ", "
        << getMatrix()->getGlobalNumCols() << "]"
        << ", Global nnz: " << getMatrix()->getGlobalNumEntries();
  }
  std::cout << "DatabaseSchwarz description 3..." << std::endl;

  out << "}";
  std::cout << "DatabaseSchwarz description done!" << std::endl;
  return out.str();
}


template <class MatrixType>
void DatabaseSchwarz<MatrixType>::
describe (Teuchos::FancyOStream& out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using Teuchos::TypeNameTraits;
  using std::endl;

  // Default verbosity level is VERB_LOW
  const Teuchos::EVerbosityLevel vl =
    (verbLevel == Teuchos::VERB_DEFAULT) ? Teuchos::VERB_LOW : verbLevel;

  if (vl == Teuchos::VERB_NONE) {
    return; // print NOTHING, not even the class name
  }

  // By convention, describe() starts with a tab.
  //
  // This does affect all processes on which it's valid to print to
  // 'out'.  However, it does not actually print spaces to 'out'
  // unless operator<< gets called, so it's safe to use on all
  // processes.
  Teuchos::OSTab tab0 (out);
  const int myRank = this->getComm()->getRank();
  if (myRank == 0) {
    // Output is a valid YAML dictionary.
    // In particular, we quote keys with colons in them.
    out << "\"Ifpack2::DatabaseSchwarz\":" << endl;
  }

  Teuchos::OSTab tab1 (out);
  if (vl >= Teuchos::VERB_LOW && myRank == 0) {
    out << "Template parameters:" << endl;
    {
      Teuchos::OSTab tab2 (out);
      out << "Scalar: " << TypeNameTraits<scalar_type>::name() << endl
          << "LocalOrdinal: " << TypeNameTraits<local_ordinal_type>::name() << endl
          << "GlobalOrdinal: " << TypeNameTraits<global_ordinal_type>::name() << endl
          << "Device: " << TypeNameTraits<device_type>::name() << endl;
    }
    out << "Initialized: " << (isInitialized() ? "true" : "false") << endl
        << "Computed: " << (isComputed() ? "true" : "false") << endl;

    if (getMatrix().is_null()) {
      out << "Matrix: null" << endl;
    }
    else {
      out << "Global matrix dimensions: ["
          << getMatrix()->getGlobalNumRows() << ", "
          << getMatrix()->getGlobalNumCols() << "]" << endl
          << "Global nnz: " << getMatrix()->getGlobalNumEntries() << endl;
    }
  }
}



}//namespace Ifpack2

#define IFPACK2_DATABASESCHWARZ_INSTANT(S,LO,GO,N)                            \
  template class Ifpack2::DatabaseSchwarz< Tpetra::RowMatrix<S, LO, GO, N> >;

#endif // IFPACK2_DATABASESCHWARZ_DEF_HPP
