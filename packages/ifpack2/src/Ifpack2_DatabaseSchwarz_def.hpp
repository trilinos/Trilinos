// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_DATABASESCHWARZ_DEF_HPP
#define IFPACK2_DATABASESCHWARZ_DEF_HPP

#include "Ifpack2_Parameters.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_oblackholestream.hpp"
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
    ApplyFlops_(0.0),
    PatchSize_(9),
    PatchTolerance_(1e-3),
    SkipDatabase_(false),
    Verbose_(false)
{
  this->setObjectLabel("Ifpack2::DatabaseSchwarz");
}


template<class MatrixType>
DatabaseSchwarz<MatrixType>::
DatabaseSchwarz (const Teuchos::RCP<const row_matrix_type>& A,
                 Teuchos::ParameterList& params)
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
    ApplyFlops_(0.0),
    PatchSize_(9),
    PatchTolerance_(1e-3),
    SkipDatabase_(false),
    Verbose_(false)
{
  this->setObjectLabel("Ifpack2::DatabaseSchwarz");
  this->setParameters(params);
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
DatabaseSchwarz<MatrixType>::setParameters(const Teuchos::ParameterList& params)
{
  // GH: Copied from CAG and others. Yes, const_cast bad.
  this->setParametersImpl(const_cast<Teuchos::ParameterList&>(params));
}


template<class MatrixType>
void
DatabaseSchwarz<MatrixType>::setParametersImpl(Teuchos::ParameterList& params)
{
  // GH: since patch size varies dramatically, it doesn't make sense to force a default size
  // but I don't know if it's any better to throw if the user doesn't provide a patch size
  const int defaultPatchSize = 9;
  const double defaultPatchTolerance = 1e-3;
  const bool defaultSkipDatabase = false;
  const bool defaultVerbose = false;

  // the size of patch to search for
  PatchSize_ = params.get<int>("database schwarz: patch size",defaultPatchSize);
  
  TEUCHOS_TEST_FOR_EXCEPTION(
    PatchSize_ < 0, std::invalid_argument,
    "Ifpack2::DatabaseSchwarz::setParameters: \"database schwarz: patch size\" parameter "
    "must be a nonnegative integer.  You gave a value of " << PatchSize_ << ".");

  // the tolerance at which two patch matrices are considered "equal"
  PatchTolerance_ = params.get("database schwarz: patch tolerance", defaultPatchTolerance);

  TEUCHOS_TEST_FOR_EXCEPTION(
    PatchTolerance_ <= 0, std::invalid_argument,
    "Ifpack2::DatabaseSchwarz::setParameters: \"database schwarz: patch tolerance\" parameter "
    "must be a positive double.  You gave a value of " << PatchTolerance_ << ".");

  // whether to skip the database computation and invert all patches or not
  SkipDatabase_ = params.get<bool>("database schwarz: skip database",defaultSkipDatabase);

  // verbosity: controls whether to print a database summary at the end of the compute phase
  Verbose_ = params.get<bool>("database schwarz: print database summary",defaultVerbose);
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
  return A_;
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
  return false;
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
    auto X_view = X.getLocalViewHost(Tpetra::Access::ReadOnly);
    auto Y_view = Y.getLocalViewHost(Tpetra::Access::ReadWrite);

    Teuchos::LAPACK<int, typename row_matrix_type::scalar_type> lapack;
    int INFO = 0;
    for(unsigned int ipatch=0; ipatch<NumPatches_; ipatch++) {
      int idatabase = DatabaseIndices_[ipatch];

      // 2a. Split X into Xk on each local patch
      Teuchos::Array<typename row_matrix_type::scalar_type> x_patch(PatchSize_);
      for(unsigned int c=0; c<X_view.extent(1); ++c) {
        for(unsigned int i=0; i<x_patch.size(); ++i) {
          x_patch[i] = X_view(PatchIndices_[ipatch][i],c);
        }
      }

      // 2b. Solve each using Lapack::GETRS
      // GH: TODO: can improve this by grouping all patches such that DatabaseIndices_[ipatch] is equal
      // in the compute phase and then utilizing the multiple RHS capability here.
      int numRhs = 1;
      
      int* ipiv = &ipiv_[idatabase*PatchSize_];

      lapack.GETRS('N', DatabaseMatrices_[idatabase]->numRows(), numRhs,
                   DatabaseMatrices_[idatabase]->values(), DatabaseMatrices_[idatabase]->numRows(),
                   ipiv, x_patch.getRawPtr(),
                   DatabaseMatrices_[idatabase]->numRows(), &INFO);

      // INFO < 0 is a bug.
      TEUCHOS_TEST_FOR_EXCEPTION(
        INFO < 0, std::logic_error, "Ifpack2::DatabaseSchwarz::compute: "
        "LAPACK's _GETRF (LU factorization with partial pivoting) was called "
        "incorrectly.  INFO = " << INFO << " < 0.  "
        "Please report this bug to the Ifpack2 developers.");
      // INFO > 0 means the matrix is singular.  This is probably an issue
      // either with the choice of rows the rows we extracted, or with the
      // input matrix itself.
      TEUCHOS_TEST_FOR_EXCEPTION(
        INFO > 0, std::runtime_error, "Ifpack2::DatabaseSchwarz::compute: "
        "LAPACK's _GETRF (LU factorization with partial pivoting) reports that the "
        "computed U factor is exactly singular.  U(" << INFO << "," << INFO << ") "
        "(one-based index i) is exactly zero.  This probably means that the input "
        "patch is singular.");

      // 2c. Add alpha*weights*Xk into Y for each Xk
      for(unsigned int c=0; c<Y_view.extent(1); ++c) {
        for(unsigned int i=0; i<x_patch.size(); ++i) {
          Y_view(PatchIndices_[ipatch][i],c) += alpha*Weights_[PatchIndices_[ipatch][i]]*x_patch[i];
        }
      }
    }
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
    const int maxNnzPerRow = A_->getGlobalMaxNumRowEntries();

    // Phase 1: Loop over rows of A_, construct patch indices, and construct patch matrices
    PatchIndices_.resize(0);
    NumPatches_ = 0;
    // loop over potential candidates by checking rows of A_
    for(local_ordinal_type irow=0; irow < (local_ordinal_type) A_->getLocalNumRows(); ++irow) {
      size_t num_entries = A_->getNumEntriesInLocalRow(irow);

      // if irow is a potential patch candidate
      if((local_ordinal_type) num_entries == PatchSize_) {
        // grab row irow of A_
        typename row_matrix_type::nonconst_local_inds_host_view_type row_inds("row indices", maxNnzPerRow);
        typename row_matrix_type::nonconst_values_host_view_type row_vals("row values", maxNnzPerRow);
        A_->getLocalRowCopy(irow, row_inds, row_vals, num_entries);

        // check if we've added DOF irow before
        bool is_new_patch = true;
        for(size_t ipatch=0; ipatch<PatchIndices_.size(); ++ipatch) {
          for(size_t i=0; i<PatchIndices_[ipatch].size(); ++i) {
            if(PatchIndices_[ipatch][i] == irow) {
              is_new_patch = false;
              ipatch=PatchIndices_.size(); // likely the ugliest way to break out other than using goto TheEnd:
              break;
            }
          }
        }

        // if this patch is new, append the indices
        if(is_new_patch) {
          Teuchos::ArrayView<typename row_matrix_type::local_ordinal_type> indices_array_view(row_inds.data(),num_entries);
          std::vector<typename row_matrix_type::local_ordinal_type> indices_vector = Teuchos::createVector(indices_array_view);
          PatchIndices_.push_back(indices_vector);
          NumPatches_++;
        }
      }
    }
    
    // Phase 2: construct the list of local patch matrices
    typedef typename Teuchos::SerialDenseMatrix<typename row_matrix_type::local_ordinal_type,typename row_matrix_type::scalar_type> DenseMatType;
    typedef Teuchos::RCP<DenseMatType> DenseMatRCP;
    DatabaseIndices_.resize(NumPatches_,-1);
    Weights_.resize(A_->getLocalNumRows(),0);
    std::vector<double> index_count(A_->getLocalNumRows(),0);
    // compute the local patch matrix A_k by grabbing values from the rows of A_
    for(unsigned int ipatch=0; ipatch< NumPatches_; ++ipatch) {
      // form a local patch matrix and grab the indices for its rows/columns
      DenseMatRCP patch_matrix = Teuchos::rcp(new DenseMatType(PatchSize_, PatchSize_));
      auto indices_vector = PatchIndices_[ipatch];

      for(local_ordinal_type i=0; i<PatchSize_; ++i) {
        index_count[indices_vector[i]]++;
        // grab each row from A_ and throw them into patch_matrix
        typename row_matrix_type::nonconst_local_inds_host_view_type row_inds("row indices", maxNnzPerRow);
        typename row_matrix_type::nonconst_values_host_view_type row_vals("row values", maxNnzPerRow);
        size_t num_entries;
        A_->getLocalRowCopy(indices_vector[i], row_inds, row_vals, num_entries);
        for(local_ordinal_type j=0; j<PatchSize_; ++j) {
          for(size_t k=0; k<num_entries; ++k) {
            if(row_inds(k) == indices_vector[j]) {
              (*patch_matrix)(i,j) = row_vals(k);
            }
          }
        }
      }
      
      // check if the local patch matrix has been seen before
      // this is skipped and the patch is stored anyway if SkipDatabase_ is true
      bool found_match = false;
      if(!SkipDatabase_) {
        for(size_t idatabase=0; idatabase<DatabaseMatrices_.size(); ++idatabase) {

          // sum errors
          typename Teuchos::ScalarTraits<typename row_matrix_type::scalar_type>::magnitudeType abserror = 0.0;
          for(local_ordinal_type irow=0; irow<PatchSize_; irow++) {
            for(local_ordinal_type icol=0; icol<PatchSize_; ++icol) {
              DenseMatRCP database_candidate = DatabaseMatrices_[idatabase];
              abserror += Teuchos::ScalarTraits<typename row_matrix_type::scalar_type>::magnitude((*patch_matrix)(irow,icol)-(*database_candidate)(irow,icol));
            }
            // break out early if we finish a row and the error is already too high
            if(abserror > Teuchos::as<typename Teuchos::ScalarTraits<typename row_matrix_type::scalar_type>::magnitudeType>(PatchTolerance_))
              break;
          }

          // check if this error is acceptable; if so, mark the match and break
          if(abserror < Teuchos::as<typename Teuchos::ScalarTraits<typename row_matrix_type::scalar_type>::magnitudeType>(PatchTolerance_)) {
            DatabaseIndices_[ipatch] = idatabase;
            found_match = true;
            break;
          }
        }
      }

      // if no match was found, append patch_matrix to the database
      if(!found_match) {
        DatabaseMatrices_.push_back(patch_matrix);
        DatabaseIndices_[ipatch] = DatabaseMatrices_.size()-1;
        TEUCHOS_TEST_FOR_EXCEPTION(DatabaseMatrices_[DatabaseMatrices_.size()-1].is_null(), std::logic_error,
          "Ifpack2::DatabaseSchwarz::compute: A matrix was added to the database, but appears to be null!"
          "Please report this bug to the Ifpack2 developers.");
      }
    }

    // compute proc-local overlap weights
    for(unsigned int i=0; i<index_count.size(); ++i) {
      TEUCHOS_TEST_FOR_EXCEPTION(index_count[i] == 0.0, std::logic_error,
        "Ifpack2::DatabaseSchwarz::compute: DOF " << i << " has no corresponding patch! "
        "All DOFs must be able to locate in a patch to use this method! "
        "Have you considered adjusting the patch size? Currently it is " << PatchSize_ << ".");
      Weights_[i] = 1./index_count[i];
    }
    DatabaseSize_ = DatabaseMatrices_.size();
    
    // compute how many patches refer to a given database entry
    std::vector<int> database_counts(DatabaseSize_,0);
    for(unsigned int ipatch=0; ipatch<NumPatches_; ++ipatch) {
      database_counts[DatabaseIndices_[ipatch]]++;
    }

    // Phase 3: factor the patches using LAPACK (GETRF for factorization)
    Teuchos::LAPACK<int, typename row_matrix_type::scalar_type> lapack;
    int INFO = 0;
    ipiv_.resize(DatabaseSize_*PatchSize_);
    std::fill(ipiv_.begin (), ipiv_.end (), 0);
    for(unsigned int idatabase=0; idatabase<DatabaseSize_; idatabase++) {
      int* ipiv = &ipiv_[idatabase*PatchSize_];

      lapack.GETRF(DatabaseMatrices_[idatabase]->numRows(),
                   DatabaseMatrices_[idatabase]->numCols(),
                   DatabaseMatrices_[idatabase]->values(),
                   DatabaseMatrices_[idatabase]->stride(),
                   ipiv, &INFO);

      // INFO < 0 is a bug.
      TEUCHOS_TEST_FOR_EXCEPTION(
        INFO < 0, std::logic_error, "Ifpack2::DatabaseSchwarz::compute: "
        "LAPACK's _GETRF (LU factorization with partial pivoting) was called "
        "incorrectly.  INFO = " << INFO << " < 0.  "
        "Please report this bug to the Ifpack2 developers.");
      // INFO > 0 means the matrix is singular.  This is probably an issue
      // either with the choice of rows the rows we extracted, or with the
      // input matrix itself.
      if(INFO > 0) {
        std::cout << "SINGULAR LOCAL MATRIX, COUNT=" << database_counts[idatabase] << std::endl;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(
        INFO > 0, std::runtime_error, "Ifpack2::DatabaseSchwarz::compute: "
        "LAPACK's _GETRF (LU factorization with partial pivoting) reports that the "
        "computed U factor is exactly singular.  U(" << INFO << "," << INFO << ") "
        "(one-based index i) is exactly zero.  This probably means that the input "
        "patch is singular.");
    }
  }
  IsComputed_ = true;
  ++NumCompute_;

  ComputeTime_ += (timer->wallTime() - startTime);

  // print a summary after compute finishes if Verbose_ is true (TODO: fancyostream)
  if(Verbose_) {
    std::cout << "Ifpack2::DatabaseSchwarz()::Compute() summary\n";
    std::cout << "Found " << NumPatches_ << " patches of size " << PatchSize_ << " in matrix A\n";
    std::cout << "Database tol = " << PatchTolerance_ << "\n";
    std::cout << "Database size = " << DatabaseSize_ << " patches\n";
  }
}


template <class MatrixType>
std::string DatabaseSchwarz<MatrixType>::description() const {
  std::ostringstream out;

  // Output is a valid YAML dictionary in flow style.  If you don't
  // like everything on a single line, you should call describe()
  // instead.
  out << "\"Ifpack2::DatabaseSchwarz\": {";
  out << "Initialized: " << (isInitialized() ? "true" : "false")
      << ", Computed: " << (isComputed() ? "true" : "false")
      << ", patch size: " << PatchSize_
      << ", patch tolerance: " << PatchTolerance_
      << ", skip database: " << (SkipDatabase_ ? "true" : "false")
      << ", print database summary: " << (Verbose_ ? "true" : "false");

  if (getMatrix().is_null()) {
    out << "Matrix: null";
  }
  else {
    out << "Global matrix dimensions: ["
        << getMatrix()->getGlobalNumRows() << ", "
        << getMatrix()->getGlobalNumCols() << "]"
        << ", Global nnz: " << getMatrix()->getGlobalNumEntries();
  }
  
  out << "}";
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
