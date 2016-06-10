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

#ifndef IFPACK2_BLOCKRELAXATION_DEF_HPP
#define IFPACK2_BLOCKRELAXATION_DEF_HPP

#include "Ifpack2_BlockRelaxation_decl.hpp"
#include "Ifpack2_LinearPartitioner.hpp"
#include "Ifpack2_LinePartitioner.hpp"
#include "Ifpack2_Details_UserPartitioner_decl.hpp"
#include "Ifpack2_Details_UserPartitioner_def.hpp"
#include <Ifpack2_Parameters.hpp>

namespace Ifpack2 {

template<class MatrixType,class ContainerType>
void BlockRelaxation<MatrixType,ContainerType>::
setMatrix (const Teuchos::RCP<const row_matrix_type>& A)
{
  if (A.getRawPtr () != A_.getRawPtr ()) { // it's a different matrix
    IsInitialized_ = false;
    IsComputed_ = false;
    Partitioner_ = Teuchos::null;
    Importer_ = Teuchos::null;
    W_ = Teuchos::null;
    hasBlockCrsMatrix_ = false;

    if (! A.is_null ()) {
      IsParallel_ = (A->getRowMap ()->getComm ()->getSize () > 1);
    }

    std::vector<Teuchos::RCP<ContainerType> > emptyVec;
    std::swap (Containers_, emptyVec);
    NumLocalBlocks_ = 0;

    A_ = A;
  }
}

template<class MatrixType,class ContainerType>
BlockRelaxation<MatrixType,ContainerType>::
BlockRelaxation (const Teuchos::RCP<const row_matrix_type>& A)
: A_ (A),
  Time_ (Teuchos::rcp (new Teuchos::Time ("Ifpack2::BlockRelaxation"))),
  OverlapLevel_ (0),
  PartitionerType_ ("linear"),
  NumSweeps_ (1),
  NumLocalBlocks_(0),
  PrecType_ (Ifpack2::Details::JACOBI),
  DampingFactor_ (STS::one ()),
  IsParallel_ (false),
  ZeroStartingSolution_ (true),
  DoBackwardGS_ (false),
  IsInitialized_ (false),
  IsComputed_ (false),
  NumInitialize_ (0),
  NumCompute_ (0),
  NumApply_ (0),
  InitializeTime_ (0.0),
  ComputeTime_ (0.0),
  ApplyTime_ (0.0),
  ComputeFlops_ (0.0),
  ApplyFlops_ (0.0),
  NumMyRows_ (0),
  NumGlobalRows_ (0),
  NumGlobalNonzeros_ (0),
  hasBlockCrsMatrix_ (false)
{}

template<class MatrixType,class ContainerType>
BlockRelaxation<MatrixType,ContainerType>::
~BlockRelaxation ()
{}

template<class MatrixType,class ContainerType>
void
BlockRelaxation<MatrixType,ContainerType>::
setParameters (const Teuchos::ParameterList& List)
{
  Teuchos::ParameterList validparams;
  Ifpack2::getValidParameters (validparams);
  List.validateParameters (validparams);

  std::string PT;
  if (PrecType_ == Ifpack2::Details::JACOBI) {
    PT = "Jacobi";
  } else if (PrecType_ == Ifpack2::Details::GS) {
    PT = "Gauss-Seidel";
  } else if (PrecType_ == Ifpack2::Details::SGS) {
    PT = "Symmetric Gauss-Seidel";
  }

  Ifpack2::getParameter (List, "relaxation: type", PT);

  if (PT == "Jacobi") {
    PrecType_ = Ifpack2::Details::JACOBI;
  }
  else if (PT == "Gauss-Seidel") {
    PrecType_ = Ifpack2::Details::GS;
  }
  else if (PT == "Symmetric Gauss-Seidel") {
    PrecType_ = Ifpack2::Details::SGS;
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::invalid_argument, "Ifpack2::BlockRelaxation::setParameters: "
      "Invalid parameter value \"" << PT << "\" for parameter \"relaxation: "
      "type\".");
  }

  Ifpack2::getParameter (List, "relaxation: sweeps",NumSweeps_);
  Ifpack2::getParameter (List, "relaxation: damping factor", DampingFactor_);
  Ifpack2::getParameter (List, "relaxation: zero starting solution", ZeroStartingSolution_);
  Ifpack2::getParameter (List, "relaxation: backward mode",DoBackwardGS_);
  Ifpack2::getParameter (List, "partitioner: type",PartitionerType_);
  Ifpack2::getParameter (List, "partitioner: local parts",NumLocalBlocks_);
  Ifpack2::getParameter (List, "partitioner: overlap",OverlapLevel_);

  // check parameters
  if (PrecType_ != Ifpack2::Details::JACOBI) {
    OverlapLevel_ = 0;
  }
  if (NumLocalBlocks_ < 0) {
    NumLocalBlocks_ = A_->getNodeNumRows() / (-NumLocalBlocks_);
  }
  // other checks are performed in Partitioner_

  // NTS: Sanity check to be removed at a later date when Backward mode is enabled
  TEUCHOS_TEST_FOR_EXCEPTION(
    DoBackwardGS_, std::runtime_error,
    "Ifpack2::BlockRelaxation:setParameters: Setting the \"relaxation: "
    "backward mode\" parameter to true is not yet supported.");

  // copy the list as each subblock's constructor will
  // require it later
  List_ = List;
}

template<class MatrixType,class ContainerType>
Teuchos::RCP<const Teuchos::Comm<int> >
BlockRelaxation<MatrixType,ContainerType>::getComm () const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error, "Ifpack2::BlockRelaxation::getComm: "
     "The matrix is null.  You must call setMatrix() with a nonnull matrix "
     "before you may call this method.");
  return A_->getComm ();
}

template<class MatrixType,class ContainerType>
Teuchos::RCP<const Tpetra::RowMatrix<typename MatrixType::scalar_type,
                                     typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >
BlockRelaxation<MatrixType,ContainerType>::getMatrix () const {
  return A_;
}

template<class MatrixType,class ContainerType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                               typename MatrixType::global_ordinal_type,
                               typename MatrixType::node_type> >
BlockRelaxation<MatrixType,ContainerType>::
getDomainMap () const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error, "Ifpack2::BlockRelaxation::"
     "getDomainMap: The matrix is null.  You must call setMatrix() with a "
     "nonnull matrix before you may call this method.");
  return A_->getDomainMap ();
}

template<class MatrixType,class ContainerType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                               typename MatrixType::global_ordinal_type,
                               typename MatrixType::node_type> >
BlockRelaxation<MatrixType,ContainerType>::
getRangeMap () const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error, "Ifpack2::BlockRelaxation::"
     "getRangeMap: The matrix is null.  You must call setMatrix() with a "
     "nonnull matrix before you may call this method.");
  return A_->getRangeMap ();
}

template<class MatrixType,class ContainerType>
bool
BlockRelaxation<MatrixType,ContainerType>::
hasTransposeApply () const {
  return true;
}

template<class MatrixType,class ContainerType>
int
BlockRelaxation<MatrixType,ContainerType>::
getNumInitialize () const {
  return NumInitialize_;
}

template<class MatrixType,class ContainerType>
int
BlockRelaxation<MatrixType,ContainerType>::
getNumCompute () const
{
  return NumCompute_;
}

template<class MatrixType,class ContainerType>
int
BlockRelaxation<MatrixType,ContainerType>::
getNumApply () const
{
  return NumApply_;
}

template<class MatrixType,class ContainerType>
double
BlockRelaxation<MatrixType,ContainerType>::
getInitializeTime () const
{
  return InitializeTime_;
}

template<class MatrixType,class ContainerType>
double
BlockRelaxation<MatrixType,ContainerType>::
getComputeTime () const
{
  return ComputeTime_;
}

template<class MatrixType,class ContainerType>
double
BlockRelaxation<MatrixType,ContainerType>::
getApplyTime () const
{
  return ApplyTime_;
}

template<class MatrixType,class ContainerType>
double
BlockRelaxation<MatrixType,ContainerType>::
getComputeFlops () const
{
  return ComputeFlops_;
}

template<class MatrixType,class ContainerType>
double
BlockRelaxation<MatrixType,ContainerType>::
getApplyFlops () const
{
  return ApplyFlops_;
}

template<class MatrixType,class ContainerType>
void
BlockRelaxation<MatrixType,ContainerType>::
apply (const Tpetra::MultiVector<typename MatrixType::scalar_type,
                                 typename MatrixType::local_ordinal_type,
                                 typename MatrixType::global_ordinal_type,
                                 typename MatrixType::node_type>& X,
       Tpetra::MultiVector<typename MatrixType::scalar_type,
                           typename MatrixType::local_ordinal_type,
                           typename MatrixType::global_ordinal_type,
                           typename MatrixType::node_type>& Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error, "Ifpack2::BlockRelaxation::apply: "
     "The matrix is null.  You must call setMatrix() with a nonnull matrix, "
     "then call initialize() and compute() (in that order), before you may "
     "call this method.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! isComputed (), std::runtime_error, "Ifpack2::BlockRelaxation::apply: "
    "isComputed() must be true prior to calling apply.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.getNumVectors () != Y.getNumVectors (), std::invalid_argument,
    "Ifpack2::BlockRelaxation::apply: X.getNumVectors() = "
    << X.getNumVectors () << " != Y.getNumVectors() = "
    << Y.getNumVectors () << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(
    mode != Teuchos::NO_TRANS, std::logic_error, "Ifpack2::BlockRelaxation::"
    "apply: This method currently only implements the case mode == Teuchos::"
    "NO_TRANS.  Transposed modes are not currently supported.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    alpha != Teuchos::ScalarTraits<scalar_type>::one (), std::logic_error,
    "Ifpack2::BlockRelaxation::apply: This method currently only implements "
    "the case alpha == 1.  You specified alpha = " << alpha << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(
    beta != Teuchos::ScalarTraits<scalar_type>::zero (), std::logic_error,
    "Ifpack2::BlockRelaxation::apply: This method currently only implements "
    "the case beta == 0.  You specified beta = " << beta << ".");

  Time_->start(true);

  // If X and Y are pointing to the same memory location,
  // we need to create an auxiliary vector, Xcopy
  Teuchos::RCP<const MV> X_copy;
  {
    auto X_lcl_host = X.template getLocalView<Kokkos::HostSpace> ();
    auto Y_lcl_host = Y.template getLocalView<Kokkos::HostSpace> ();
    if (X_lcl_host.ptr_on_device () == Y_lcl_host.ptr_on_device ()) {
      X_copy = rcp (new MV (X, Teuchos::Copy));
    } else {
      X_copy = rcpFromRef (X);
    }
  }

  if (ZeroStartingSolution_) {
    Y.putScalar (STS::zero ());
  }

  // Flops are updated in each of the following.
  switch (PrecType_) {
  case Ifpack2::Details::JACOBI:
    ApplyInverseJacobi(*X_copy,Y);
    break;
  case Ifpack2::Details::GS:
    ApplyInverseGS(*X_copy,Y);
    break;
  case Ifpack2::Details::SGS:
    ApplyInverseSGS(*X_copy,Y);
    break;
  default:
    TEUCHOS_TEST_FOR_EXCEPTION
      (true, std::logic_error, "Ifpack2::BlockRelaxation::apply: Invalid "
       "PrecType_ enum value " << PrecType_ << ".  Valid values are Ifpack2::"
       "Details::JACOBI = " << Ifpack2::Details::JACOBI << ", Ifpack2::Details"
       "::GS = " << Ifpack2::Details::GS << ", and Ifpack2::Details::SGS = "
       << Ifpack2::Details::SGS << ".  Please report this bug to the Ifpack2 "
       "developers.");
  }

  ++NumApply_;
  Time_->stop();
  ApplyTime_ += Time_->totalElapsedTime();
}

template<class MatrixType,class ContainerType>
void
BlockRelaxation<MatrixType,ContainerType>::
applyMat (const Tpetra::MultiVector<typename MatrixType::scalar_type,
                                    typename MatrixType::local_ordinal_type,
                                    typename MatrixType::global_ordinal_type,
                                    typename MatrixType::node_type>& X,
          Tpetra::MultiVector<typename MatrixType::scalar_type,
                             typename MatrixType::local_ordinal_type,
                             typename MatrixType::global_ordinal_type,
                             typename MatrixType::node_type>& Y,
          Teuchos::ETransp mode) const
{
  A_->apply (X, Y, mode);
}

template<class MatrixType,class ContainerType>
void
BlockRelaxation<MatrixType,ContainerType>::
initialize ()
{
  using Teuchos::rcp;
  typedef Tpetra::RowGraph<local_ordinal_type, global_ordinal_type, node_type>
    row_graph_type;

  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error, "Ifpack2::BlockRelaxation::initialize: "
     "The matrix is null.  You must call setMatrix() with a nonnull matrix "
     "before you may call this method.");

  // Check whether we have a BlockCrsMatrix
  Teuchos::RCP<const block_crs_matrix_type> A_bcrs =
    Teuchos::rcp_dynamic_cast<const block_crs_matrix_type> (A_);
  if (A_bcrs.is_null ()) {
    hasBlockCrsMatrix_ = false;
  }
  else {
    hasBlockCrsMatrix_ = true;
  }

  IsInitialized_ = false;
  Time_->start (true);

  NumMyRows_         = A_->getNodeNumRows ();
  NumGlobalRows_     = A_->getGlobalNumRows ();
  NumGlobalNonzeros_ = A_->getGlobalNumEntries ();

  // NTS: Will need to add support for Zoltan2 partitions later Also,
  // will need a better way of handling the Graph typing issue.
  // Especially with ordinal types w.r.t the container.

  if (PartitionerType_ == "linear") {
    Partitioner_ =
      rcp (new Ifpack2::LinearPartitioner<row_graph_type> (A_->getGraph ()));
  } else if (PartitionerType_ == "line") {
    Partitioner_ =
      rcp (new Ifpack2::LinePartitioner<row_graph_type,typename MatrixType::scalar_type> (A_->getGraph ()));
  } else if (PartitionerType_ == "user") {
    Partitioner_ =
      rcp (new Ifpack2::Details::UserPartitioner<row_graph_type> (A_->getGraph () ) );
  } else {
    // We should have checked for this in setParameters(), so it's a
    // logic_error, not an invalid_argument or runtime_error.
    TEUCHOS_TEST_FOR_EXCEPTION
      (true, std::logic_error, "Ifpack2::BlockRelaxation::initialize: Unknown "
       "partitioner type " << PartitionerType_ << ".  Valid values are "
       "\"linear\", \"line\", and \"user\".");
  }

  // need to partition the graph of A
  Partitioner_->setParameters (List_);
  Partitioner_->compute ();

  // get actual number of partitions
  NumLocalBlocks_ = Partitioner_->numLocalParts ();

  // Note: Unlike Ifpack, we'll punt on computing W_ until compute(), which is where
  // we assume that the type of relaxation has been chosen.

  if (A_->getComm()->getSize() != 1) {
    IsParallel_ = true;
  } else {
    IsParallel_ = false;
  }

  ++NumInitialize_;
  Time_->stop ();
  InitializeTime_ += Time_->totalElapsedTime ();
  IsInitialized_ = true;
}


template<class MatrixType,class ContainerType>
void
BlockRelaxation<MatrixType,ContainerType>::
computeBlockCrs ()
{
  typedef Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type>     map_type;
  typedef Tpetra::Import<local_ordinal_type,global_ordinal_type, node_type> import_type;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Array;
  using Teuchos::ArrayView;

  Time_->start (true);

  // reset values
  IsComputed_ = false;

  // Extract the submatrices
  ExtractSubmatrices ();

  // Compute the weight vector if we're doing overlapped Jacobi (and
  // only if we're doing overlapped Jacobi).
  TEUCHOS_TEST_FOR_EXCEPTION
    (PrecType_ == Ifpack2::Details::JACOBI && OverlapLevel_ > 0, std::runtime_error,
     "Ifpack2::BlockRelaxation::computeBlockCrs: "
     "We do not support overlapped Jacobi yet for Tpetra::BlockCrsMatrix.  Sorry!");

  // We need to import data from external processors. Here I create a
  // Tpetra::Import object if needed (stealing from A_ if possible)
  // Marzio's comment:
  // Note that I am doing some strange stuff to set the components of Y
  // from Y2 (to save some time).
  //
  if (IsParallel_ && (PrecType_ == Ifpack2::Details::GS ||
                      PrecType_ == Ifpack2::Details::SGS)) {
    // Get the block size
    RCP<const block_crs_matrix_type> A_bcrs =
      Teuchos::rcp_dynamic_cast<const block_crs_matrix_type> (A_);
    int bs = A_bcrs->getBlockSize();

    // Get the maps describing the block data
    // Each "element" is actually a block of numbers
    RCP<const map_type> oldDomainMap = A_->getDomainMap();
    RCP<const map_type> oldColMap = A_->getColMap();

    // Because A is a block CRS matrix, import will not do what you think it does
    // We have to construct the correct maps for it
    global_size_t numGlobalElements = oldColMap->getGlobalNumElements()*bs;
    global_ordinal_type indexBase = oldColMap->getIndexBase();
    RCP<const Teuchos::Comm<int> >comm = oldColMap->getComm();
    ArrayView<const global_ordinal_type> oldColElements = oldColMap->getNodeElementList();
    Array<global_ordinal_type> newColElements(bs*oldColElements.size());

    for(int i=0; i<oldColElements.size(); i++) {
      for(int j=0; j<bs; j++) {
        newColElements[i*bs+j] = oldColElements[i]*bs+j;
      }
    }
    RCP<map_type> colMap = rcp(new map_type(numGlobalElements,newColElements,indexBase,comm));

    // Create the importer
    Importer_ = rcp (new import_type (oldDomainMap, colMap));
  }

  ++NumCompute_;
  Time_->stop ();
  ComputeTime_ += Time_->totalElapsedTime();
  IsComputed_ = true;
}


template<class MatrixType,class ContainerType>
void
BlockRelaxation<MatrixType,ContainerType>::
compute ()
{
  using Teuchos::rcp;
  typedef Tpetra::Vector<scalar_type,
    local_ordinal_type, global_ordinal_type, node_type> vector_type;
  typedef Tpetra::Import<local_ordinal_type,
    global_ordinal_type, node_type> import_type;

  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error, "Ifpack2::BlockRelaxation::compute: "
     "The matrix is null.  You must call setMatrix() with a nonnull matrix, "
     "then call initialize(), before you may call this method.");

  // We should have checked for this in setParameters(), so it's a
  // logic_error, not an invalid_argument or runtime_error.
  TEUCHOS_TEST_FOR_EXCEPTION
    (NumSweeps_ < 0, std::logic_error, "Ifpack2::BlockRelaxation::compute: "
     "NumSweeps_ = " << NumSweeps_ << " < 0.");

  if (! isInitialized ()) {
    initialize ();
  }

  if (hasBlockCrsMatrix_) {
    computeBlockCrs ();
    return;
  }

  Time_->start (true);

  // reset values
  IsComputed_ = false;

  // Extract the submatrices
  ExtractSubmatrices ();

  // Compute the weight vector if we're doing overlapped Jacobi (and
  // only if we're doing overlapped Jacobi).
  if (PrecType_ == Ifpack2::Details::JACOBI && OverlapLevel_ > 0) {
    // weight of each vertex
    W_ = rcp (new vector_type (A_->getRowMap ()));
    W_->putScalar (STS::zero ());
    Teuchos::ArrayRCP<scalar_type > w_ptr = W_->getDataNonConst(0);

    for (local_ordinal_type i = 0 ; i < NumLocalBlocks_ ; ++i) {
      for (size_t j = 0 ; j < Partitioner_->numRowsInPart(i) ; ++j) {
        // FIXME (mfh 12 Sep 2014) Should this really be int?
        // Perhaps it should be local_ordinal_type instead.
        int LID = (*Partitioner_)(i,j);
        w_ptr[LID]+= STS::one();
      }
    }
    W_->reciprocal (*W_);
  }

  // We need to import data from external processors. Here I create a
  // Tpetra::Import object if needed (stealing from A_ if possible)
  // Marzio's comment:
  // Note that I am doing some strange stuff to set the components of Y
  // from Y2 (to save some time).
  //
  if (IsParallel_ && (PrecType_ == Ifpack2::Details::GS ||
                      PrecType_ == Ifpack2::Details::SGS)) {
    Importer_ = A_->getGraph ()->getImporter ();
    if (Importer_.is_null ()) {
      Importer_ = rcp (new import_type (A_->getDomainMap (), A_->getColMap ()));
    }
  }

  ++NumCompute_;
  Time_->stop ();
  ComputeTime_ += Time_->totalElapsedTime();
  IsComputed_ = true;
}

template<class MatrixType,class ContainerType>
void
BlockRelaxation<MatrixType,ContainerType>::
ExtractSubmatrices ()
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (Partitioner_.is_null (), std::runtime_error, "Ifpack2::BlockRelaxation::"
     "ExtractSubmatrices: Partitioner object is null.");

  NumLocalBlocks_ = Partitioner_->numLocalParts ();
  Containers_.resize (NumLocalBlocks_);

  for (local_ordinal_type i = 0; i < NumLocalBlocks_; ++i) {
    const size_t numRows = Partitioner_->numRowsInPart (i);

    // Extract a list of the indices of each partitioner row.
    Teuchos::Array<local_ordinal_type> localRows (numRows);
    for (size_t j = 0; j < numRows; ++j) {
      localRows[j] = (*Partitioner_) (i,j);
    }
    if(numRows>1 || hasBlockCrsMatrix_) { // only do for non-singletons
      Containers_[i] = Teuchos::rcp (new ContainerType (A_, localRows ()));
      Containers_[i]->setParameters (List_);
      Containers_[i]->initialize ();
      Containers_[i]->compute ();
    }
  }
}

template<class MatrixType,class ContainerType>
void
BlockRelaxation<MatrixType,ContainerType>::
ApplyInverseJacobi (const MV& X, MV& Y) const
{
  const size_t NumVectors = X.getNumVectors ();
  MV AY (Y.getMap (), NumVectors);

  // Initial matvec not needed
  int starting_iteration = 0;
  if (ZeroStartingSolution_) {
    DoJacobi (X, Y);
    starting_iteration = 1;
  }

  const scalar_type ONE = STS::one ();
  for (int j = starting_iteration; j < NumSweeps_; ++j) {
    applyMat (Y, AY);
    AY.update (ONE, X, -ONE);
    DoJacobi (AY, Y);

    // Flops for matrix apply & update
    ApplyFlops_ += NumVectors * (2 * NumGlobalNonzeros_ + 2 * NumGlobalRows_);
  }

}

template<class MatrixType,class ContainerType>
void
BlockRelaxation<MatrixType,ContainerType>::
DoJacobi (const MV& X, MV& Y) const
{
  const size_t NumVectors = X.getNumVectors ();
  const scalar_type one = STS::one ();
  // Note: Flop counts copied naively from Ifpack.

  if (OverlapLevel_ == 0) {
    // Non-overlapping Jacobi
    for (local_ordinal_type i = 0; i < NumLocalBlocks_; ++i) {
      // may happen that a partition is empty
      if( Partitioner_->numRowsInPart (i) > 1 || hasBlockCrsMatrix_) {
        if(Containers_[i]->getNumRows () == 0 ) continue;
        Containers_[i]->apply (X, Y, Teuchos::NO_TRANS, DampingFactor_, one);
        ApplyFlops_ += NumVectors * 2 * NumGlobalRows_;
      }
      else if( Partitioner_->numRowsInPart (i) == 1) {
        // singleton, can't access Containers_[i] as it was never filled and may be null.
        local_ordinal_type LRID  = (*Partitioner_)(i,0);  // by definition, a singleton 1 row in block.
        getMatDiag();
        Teuchos::ArrayRCP< const scalar_type > Diag   = DiagRCP_->getData();
        scalar_type d = Diag[LRID];
        for(unsigned int nv = 0;nv < NumVectors ; ++nv ) {
          Teuchos::ArrayRCP< const scalar_type > xRCP = X.getData(nv);
          scalar_type x = xRCP[LRID];
          Teuchos::ArrayRCP<  scalar_type > yRCP = Y.getDataNonConst(nv);

          scalar_type newy= x/d;
          yRCP[LRID]= newy;
        }
      }
    }
  }
  else {
    // Overlapping Jacobi
    for (local_ordinal_type i = 0 ; i < NumLocalBlocks_ ; i++) {
      // may happen that a partition is empty
      if (Containers_[i]->getNumRows() == 0) continue;
      if ( Partitioner_->numRowsInPart (i)  > 1 ) {
        try {
          Containers_[i]->weightedApply(X,Y,*W_,Teuchos::NO_TRANS,DampingFactor_,one);
        } catch (std::exception& e) {
          std::cerr << "BlockRelaxation::DoJacobi: Containers_[" << i
                    << "]->weightedApply() threw an exception: " << e.what ()
                    << std::endl;
          throw;
        }
      } // end  Partitioner_->numRowsInPart (i)  != 1

      // NOTE: do not count (for simplicity) the flops due to overlapping rows
      ApplyFlops_ += NumVectors * 4 * NumGlobalRows_;
    }
  }
}

template<class MatrixType,class ContainerType>
void
BlockRelaxation<MatrixType,ContainerType>::
ApplyInverseGS (const MV& X, MV& Y) const
{
  MV Xcopy (X, Teuchos::Copy);
  for (int j = 0; j < NumSweeps_; ++j) {
    DoGaussSeidel (Xcopy, Y);
    if (j != NumSweeps_ - 1) {
      Tpetra::deep_copy (Xcopy, X);
    }
  }
}

template<class MatrixType,class ContainerType>
void
BlockRelaxation<MatrixType,ContainerType>::
DoGaussSeidel (MV& X, MV& Y) const
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;

  // Note: Flop counts copied naively from Ifpack.

  const scalar_type one = STS::one ();
  const size_t Length = A_->getNodeMaxNumRowEntries();
  const size_t NumVectors = X.getNumVectors();
  Array<scalar_type> Values;
  Array<local_ordinal_type> Indices;
  Indices.resize (Length);

  if(hasBlockCrsMatrix_)
  {
    Teuchos::RCP<const block_crs_matrix_type> A_bcrs =
      Teuchos::rcp_dynamic_cast<const block_crs_matrix_type> (A_);
    int bs = A_bcrs->getBlockSize();
    Values.resize (bs*bs*Length);
  }
  else
    Values.resize (Length);

  // an additonal vector is needed by parallel computations
  // (note that applications through Ifpack2_AdditiveSchwarz
  // are always seen are serial)
  RCP<MV> Y2;
  if (IsParallel_) {
    Y2 = rcp (new MV (Importer_->getTargetMap (), NumVectors));
  } else {
    Y2 = rcpFromRef (Y);
  }

  // I think I decided I need two extra vectors:
  // One to store the sum of the corrections (initialized to zero)
  // One to store the temporary residual (doesn't matter if it is zeroed or not)
  // My apologies for making the names clear and meaningful. (X=RHS, Y=guess?! Nice.)
  MV Residual (X.getMap (), NumVectors, false);

  ArrayRCP<ArrayRCP<scalar_type> >           x_ptr = X.get2dViewNonConst();
  ArrayRCP<ArrayRCP<scalar_type> >           y_ptr = Y.get2dViewNonConst();
  ArrayRCP<ArrayRCP<scalar_type> >          y2_ptr = Y2->get2dViewNonConst();
  ArrayRCP<ArrayRCP<scalar_type> >    residual_ptr = Residual.get2dViewNonConst();

  // data exchange is here, once per sweep
  if (IsParallel_)  Y2->doImport(Y,*Importer_,Tpetra::INSERT);

  for (local_ordinal_type i = 0; i < NumLocalBlocks_; ++i) {
    if( Partitioner_->numRowsInPart (i) > 1 || hasBlockCrsMatrix_) {
      if (Containers_[i]->getNumRows () == 0) continue;
      // update from previous block
      ArrayView<const local_ordinal_type> localRows =
        Containers_[i]->getLocalRows ();
      const size_t localNumRows = Containers_[i]->getNumRows ();
      for (size_t j = 0; j < localNumRows; ++j) {
        const local_ordinal_type LID = localRows[j]; // Containers_[i]->ID (j);
        size_t NumEntries;
        A_->getLocalRowCopy (LID, Indices (), Values (), NumEntries);

        for (size_t m = 0; m < NumVectors; ++m) {
          ArrayView<const scalar_type> x_local = (x_ptr())[m]();
          ArrayView<scalar_type>      y2_local = (y2_ptr())[m]();
          ArrayView<scalar_type>       r_local = (residual_ptr())[m]();

          if(hasBlockCrsMatrix_) {
            Teuchos::RCP<const block_crs_matrix_type> A_bcrs =
                  Teuchos::rcp_dynamic_cast<const block_crs_matrix_type> (A_);
            int bs = A_bcrs->getBlockSize();
            for (int localR = 0; localR < bs; localR++)
              r_local[LID*bs+localR] = x_local[LID*bs+localR];
            for (size_t k = 0; k < NumEntries; ++k) {
              const local_ordinal_type col = Indices[k];
              for (int localR = 0; localR < bs; localR++) {
                for(int localC = 0; localC < bs; localC++) {
//                  std::cerr << "r[" << LID*bs+localR << "] -= Values["
//                            << k*bs*bs+localR+localC*bs << "] * y2["
//                            << col*bs+localC << "]\n";
                  r_local[LID*bs+localR] -= Values[k*bs*bs+localR+localC*bs] * y2_local[col*bs+localC];
                }
              }
            }
          }
          else {
            r_local[LID] = x_local[LID];
            for (size_t k = 0; k < NumEntries; ++k) {
              const local_ordinal_type col = Indices[k];
              r_local[LID] -= Values[k] * y2_local[col];
            }
          }
        }
      }

//      Teuchos::RCP<Teuchos::FancyOStream> wrappedStream = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cout));
//      Residual.describe (*wrappedStream, Teuchos::VERB_EXTREME);

      // solve with this block
      //
      // Note: I'm abusing the ordering information, knowing that X/Y
      // and Y2 have the same ordering for on-proc unknowns.
      //
      // Note: Add flop counts for inverse apply
      Containers_[i]->apply (Residual, *Y2, Teuchos::NO_TRANS,
                             DampingFactor_,one);

//      Y2->describe (*wrappedStream, Teuchos::VERB_EXTREME);

      // operations for all getrow's
      ApplyFlops_ += NumVectors * (2 * NumGlobalNonzeros_ + 2 * NumGlobalRows_);
    }
    else if ( Partitioner_->numRowsInPart (i) == 1) {
      // singleton, can't access Containers_[i] as it was never filled and may be null.
      // a singleton calculation is exact, all residuals should be zero.
      local_ordinal_type LRID  = (*Partitioner_)(i,0);  // by definition, a singleton 1 row in block.
      getMatDiag();
      Teuchos::ArrayRCP< const scalar_type > Diag   = DiagRCP_->getData();
      scalar_type d = Diag[LRID];
      ArrayRCP<ArrayRCP<scalar_type> >          y2_ptr2 = Y2->get2dViewNonConst();
      for(unsigned int nv = 0;nv < NumVectors ; ++nv ) {
        Teuchos::ArrayRCP< const scalar_type > xRCP = X.getData(nv);
        scalar_type x = xRCP[LRID];
        ArrayView<scalar_type>      y2_local = (y2_ptr2())[nv]();
        scalar_type newy= x/d;
        y2_local[LRID]= newy;
      }
    } // end else
  } // end for NumLocalBlocks_
  // Attention: this is delicate... Not all combinations
  // of Y2 and Y will always work (tough for ML it should be ok)
  if (IsParallel_) {
    if(hasBlockCrsMatrix_) {
      Teuchos::RCP<const block_crs_matrix_type> A_bcrs =
            Teuchos::rcp_dynamic_cast<const block_crs_matrix_type> (A_);
      int bs = A_bcrs->getBlockSize();
      for (size_t m = 0; m < NumVectors; ++m) {
        ArrayView<scalar_type> y2_local = (y2_ptr())[m]();
        ArrayView<scalar_type> y_local = (y_ptr())[m]();
        for (size_t i = 0; i < NumMyRows_*bs; ++i) {
          y_local[i] = y2_local[i];
        }
      }
    }
    else {
      for (size_t m = 0; m < NumVectors; ++m) {
        ArrayView<scalar_type> y2_local = (y2_ptr())[m]();
        ArrayView<scalar_type> y_local = (y_ptr())[m]();
        for (size_t i = 0; i < NumMyRows_; ++i) {
          y_local[i] = y2_local[i];
        }
      }
    }
  }
}

template<class MatrixType,class ContainerType>
void
BlockRelaxation<MatrixType,ContainerType>::
ApplyInverseSGS (const MV& X, MV& Y) const
{
  MV Xcopy (X, Teuchos::Copy);
  for (int j = 0; j < NumSweeps_; ++j) {
    DoSGS (Xcopy, Y);
    if (j != NumSweeps_ - 1) {
      Tpetra::deep_copy (Xcopy, X);
    }
  }
}

template<class MatrixType,class ContainerType>
void
BlockRelaxation<MatrixType,ContainerType>::
DoSGS (MV& X, MV& Y) const
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;

  const scalar_type one = STS::one ();
  const size_t Length = A_->getNodeMaxNumRowEntries();
  const size_t NumVectors = X.getNumVectors();
  Array<scalar_type> Values;
  Array<local_ordinal_type> Indices;
  Indices.resize(Length);

  if(hasBlockCrsMatrix_)
  {
    Teuchos::RCP<const block_crs_matrix_type> A_bcrs =
      Teuchos::rcp_dynamic_cast<const block_crs_matrix_type> (A_);
    int bs = A_bcrs->getBlockSize();
    Values.resize (bs*bs*Length);
  }
  else
    Values.resize (Length);

  // an additonal vector is needed by parallel computations
  // (note that applications through Ifpack2_AdditiveSchwarz
  // are always seen are serial)
  RCP<MV> Y2;
  if (IsParallel_) {
    Y2 = rcp (new MV (Importer_->getTargetMap (), NumVectors));
  } else {
    Y2 = rcpFromRef (Y);
  }

  // I think I decided I need two extra vectors:
  // One to store the sum of the corrections (initialized to zero)
  // One to store the temporary residual (doesn't matter if it is zeroed or not)
  // My apologies for making the names clear and meaningful. (X=RHS, Y=guess?! Nice.)
  MV Residual (X.getMap (), NumVectors, false);

  ArrayRCP<ArrayRCP<scalar_type> >     x_ptr       = X.get2dViewNonConst();
  ArrayRCP<ArrayRCP<scalar_type> >     y_ptr       = Y.get2dViewNonConst();
  ArrayRCP<ArrayRCP<scalar_type> >     y2_ptr      = Y2->get2dViewNonConst();
  ArrayRCP<ArrayRCP<scalar_type> >    residual_ptr = Residual.get2dViewNonConst();

  // data exchange is here, once per sweep
  if (IsParallel_) {
    Y2->doImport (Y, *Importer_, Tpetra::INSERT);
  }

  // Forward Sweep
  for (local_ordinal_type i = 0; i < NumLocalBlocks_; ++i) {
    if( Partitioner_->numRowsInPart (i) > 1 || hasBlockCrsMatrix_) {
      if (Containers_[i]->getNumRows () == 0) {
        continue; // Skip empty partitions
      }
      // update from previous block
      ArrayView<const local_ordinal_type> localRows =
        Containers_[i]->getLocalRows ();
      for (size_t j = 0; j < Containers_[i]->getNumRows (); ++j) {
        const local_ordinal_type LID = localRows[j]; // Containers_[i]->ID (j);
        size_t NumEntries;
        A_->getLocalRowCopy (LID, Indices (), Values (), NumEntries);

        //set tmpresid = initresid - A*correction
        for (size_t m = 0; m < NumVectors; ++m) {
          ArrayView<const scalar_type> x_local = (x_ptr())[m]();
          ArrayView<scalar_type>      y2_local = (y2_ptr())[m]();
          ArrayView<scalar_type>       r_local = (residual_ptr())[m]();

          if(hasBlockCrsMatrix_) {
            Teuchos::RCP<const block_crs_matrix_type> A_bcrs =
                  Teuchos::rcp_dynamic_cast<const block_crs_matrix_type> (A_);
            int bs = A_bcrs->getBlockSize();
            for (int localR = 0; localR < bs; localR++)
              r_local[LID*bs+localR] = x_local[LID*bs+localR];
            for (size_t k = 0; k < NumEntries; ++k) {
              const local_ordinal_type col = Indices[k];
              for (int localR = 0; localR < bs; localR++) {
                for(int localC = 0; localC < bs; localC++) {
                  r_local[LID*bs+localR] -= Values[k*bs*bs+localR+localC*bs] * y2_local[col*bs+localC];
                }
              }
            }
          }
          else {
            r_local[LID] = x_local[LID];
            for (size_t k = 0 ; k < NumEntries ; k++) {
              local_ordinal_type col = Indices[k];
              r_local[LID] -= Values[k] * y2_local[col];
            }
          }
        }
      }
      // solve with this block
      //
      // Note: I'm abusing the ordering information, knowing that X/Y
      // and Y2 have the same ordering for on-proc unknowns.
      //
      // Note: Add flop counts for inverse apply
      Containers_[i]->apply (Residual, *Y2, Teuchos::NO_TRANS,
                           DampingFactor_, one);

      // operations for all getrow's
      ApplyFlops_ += NumVectors * (2 * NumGlobalNonzeros_ + 2 * NumGlobalRows_);

    }
    else if( Partitioner_->numRowsInPart (i) == 1) {
      // singleton, can't access Containers_[i] as it was never filled and may be null.
      local_ordinal_type LRID  = (*Partitioner_)(i,0);  // by definition, a singleton 1 row in block.
      getMatDiag();
      Teuchos::ArrayRCP< const scalar_type > Diag   = DiagRCP_->getData();
      scalar_type d = Diag[LRID];
      for(unsigned int nv = 0;nv < NumVectors ; ++nv ) {
        Teuchos::ArrayRCP< const scalar_type > xRCP = X.getData(nv);
        scalar_type x = xRCP[LRID];
        Teuchos::ArrayRCP<  scalar_type > yRCP = Y.getDataNonConst(nv);

        scalar_type newy= x/d;
        yRCP[LRID]= newy;
      }
    } // end else
  } // end forward sweep over NumLocalBlocks

  // Reverse Sweep
  //
  // mfh 12 July 2013: The unusual iteration bounds, and the use of
  // i-1 rather than i in the loop body, ensure correctness even if
  // local_ordinal_type is unsigned.  "i = NumLocalBlocks_-1; i >= 0;
  // i--" will loop forever if local_ordinal_type is unsigned, because
  // unsigned integers are (trivially) always nonnegative.
  for (local_ordinal_type i = NumLocalBlocks_; i > 0; --i) {
    //FIXME 7-June-2016 JJH should this have an "else" optimization, similar to the forward sweeps?
    if( hasBlockCrsMatrix_ || Partitioner_->numRowsInPart (i-1) > 1 ) {
      if (Containers_[i-1]->getNumRows () == 0) continue;

      // update from previous block
      ArrayView<const local_ordinal_type> localRows =
        Containers_[i-1]->getLocalRows ();
      for (size_t j = 0; j < Containers_[i-1]->getNumRows (); ++j) {
        const local_ordinal_type LID = localRows[j]; // Containers_[i-1]->ID (j);
        size_t NumEntries;
        A_->getLocalRowCopy (LID, Indices (), Values (), NumEntries);

        //set tmpresid = initresid - A*correction
        for (size_t m = 0; m < NumVectors; ++m) {
          ArrayView<const scalar_type> x_local = (x_ptr())[m]();
          ArrayView<scalar_type>      y2_local = (y2_ptr())[m]();
          ArrayView<scalar_type>       r_local = (residual_ptr())[m]();

          if(hasBlockCrsMatrix_) {
            Teuchos::RCP<const block_crs_matrix_type> A_bcrs =
                  Teuchos::rcp_dynamic_cast<const block_crs_matrix_type> (A_);
            int bs = A_bcrs->getBlockSize();
            for (int localR = 0; localR < bs; localR++)
              r_local[LID*bs+localR] = x_local[LID*bs+localR];
            for (size_t k = 0; k < NumEntries; ++k) {
              const local_ordinal_type col = Indices[k];
              for (int localR = 0; localR < bs; localR++) {
                for(int localC = 0; localC < bs; localC++) {
                  r_local[LID*bs+localR] -= Values[k*bs*bs+localR+localC*bs] * y2_local[col*bs+localC];
                }
              }
            }
          }
          else {
            r_local [LID] = x_local[LID];
            for (size_t k = 0; k < NumEntries; ++k)  {
              local_ordinal_type col = Indices[k];
              r_local[LID] -= Values[k] * y2_local[col];
            }
          }
        }
      }

      // solve with this block
      //
      // Note: I'm abusing the ordering information, knowing that X/Y
      // and Y2 have the same ordering for on-proc unknowns.
      //
      // Note: Add flop counts for inverse apply
      Containers_[i-1]->apply (Residual, *Y2, Teuchos::NO_TRANS,
                               DampingFactor_, one);

      // operations for all getrow's
      ApplyFlops_ += NumVectors * (2 * NumGlobalNonzeros_ + 2 * NumGlobalRows_);
    } // end  Partitioner_->numRowsInPart (i) != 1 ) {
    // else do nothing, as by definition with a singleton, the residuals are zero.
  } //end reverse sweep

  // Attention: this is delicate... Not all combinations
  // of Y2 and Y will always work (though for ML it should be ok)
  if (IsParallel_) {
    if(hasBlockCrsMatrix_) {
      Teuchos::RCP<const block_crs_matrix_type> A_bcrs =
            Teuchos::rcp_dynamic_cast<const block_crs_matrix_type> (A_);
      int bs = A_bcrs->getBlockSize();
      for (size_t m = 0; m < NumVectors; ++m) {
        ArrayView<scalar_type> y2_local = (y2_ptr())[m]();
        ArrayView<scalar_type> y_local = (y_ptr())[m]();
        for (size_t i = 0; i < NumMyRows_*bs; ++i) {
          y_local[i] = y2_local[i];
        }
      }
    }
    else {
      for (size_t m = 0; m < NumVectors; ++m) {
        ArrayView<scalar_type> y2_local = (y2_ptr())[m]();
        ArrayView<scalar_type> y_local = (y_ptr())[m]();
        for (size_t i = 0; i < NumMyRows_; ++i) {
          y_local[i] = y2_local[i];
        }
      }
    }
  }
}

template<class MatrixType, class ContainerType>
std::string
BlockRelaxation<MatrixType,ContainerType>::
description () const
{
  std::ostringstream out;

  // Output is a valid YAML dictionary in flow style.  If you don't
  // like everything on a single line, you should call describe()
  // instead.
  out << "\"Ifpack2::BlockRelaxation\": {";
  if (this->getObjectLabel () != "") {
    out << "Label: \"" << this->getObjectLabel () << "\", ";
  }
  out << "Initialized: " << (isInitialized () ? "true" : "false") << ", ";
  out << "Computed: " << (isComputed () ? "true" : "false") << ", ";

  if (A_.is_null ()) {
    out << "Matrix: null, ";
  }
  else {
    out << "Matrix: not null"
        << ", Global matrix dimensions: ["
        << A_->getGlobalNumRows () << ", " << A_->getGlobalNumCols () << "], ";
  }

  // It's useful to print this instance's relaxation method.  If you
  // want more info than that, call describe() instead.
  out << "\"relaxation: type\": ";
  if (PrecType_ == Ifpack2::Details::JACOBI) {
    out << "Block Jacobi";
  } else if (PrecType_ == Ifpack2::Details::GS) {
    out << "Block Gauss-Seidel";
  } else if (PrecType_ == Ifpack2::Details::SGS) {
    out << "Block Symmetric Gauss-Seidel";
  } else {
    out << "INVALID";
  }

  out  << ", " << "sweeps: " << NumSweeps_ << ", "
      << "damping factor: " << DampingFactor_ << ", ";

  out << "}";
  return out.str();
}

template<class MatrixType,class ContainerType>
void
BlockRelaxation<MatrixType,ContainerType>::
describe (Teuchos::FancyOStream& out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using std::endl;
  using std::setw;
  using Teuchos::TypeNameTraits;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;

  Teuchos::EVerbosityLevel vl = verbLevel;
  if (vl == VERB_DEFAULT) vl = VERB_LOW;
  const int myImageID = A_->getComm()->getRank();

  // Convention requires that describe() begin with a tab.
  Teuchos::OSTab tab (out);

  //    none: print nothing
  //     low: print O(1) info from node 0
  //  medium:
  //    high:
  // extreme:
  if (vl != VERB_NONE && myImageID == 0) {
    out << "Ifpack2::BlockRelaxation<"
        << TypeNameTraits<MatrixType>::name () << ", "
        << TypeNameTraits<ContainerType>::name () << " >:";
    Teuchos::OSTab tab1 (out);

    if (this->getObjectLabel () != "") {
      out << "label: \"" << this->getObjectLabel () << "\"" << endl;
    }
    out << "initialized: " << (isInitialized () ? "true" : "false") << endl
        << "computed: " << (isComputed () ? "true" : "false") << endl;

    std::string precType;
    if (PrecType_ == Ifpack2::Details::JACOBI) {
      precType = "Block Jacobi";
    } else if (PrecType_ == Ifpack2::Details::GS) {
      precType = "Block Gauss-Seidel";
    } else if (PrecType_ == Ifpack2::Details::SGS) {
      precType = "Block Symmetric Gauss-Seidel";
    }
    out << "type: " << precType << endl
        << "global number of rows: " << A_->getGlobalNumRows () << endl
        << "global number of columns: " << A_->getGlobalNumCols () << endl
        << "number of sweeps: " << NumSweeps_ << endl
        << "damping factor: " << DampingFactor_ << endl
        << "backwards mode: "
        << ((PrecType_ == Ifpack2::Details::GS && DoBackwardGS_) ? "true" : "false")
        << endl
        << "zero starting solution: "
        << (ZeroStartingSolution_ ? "true" : "false") << endl;

    out << "===============================================================================" << endl;
    out << "Phase           # calls    Total Time (s)     Total MFlops      MFlops/s       " << endl;
    out << "------------    -------    ---------------    ---------------   ---------------" << endl;
    out << setw(12) << "initialize()" << setw(5) << getNumInitialize() << "    " << setw(15) << getInitializeTime() << endl;
    out << setw(12) << "compute()" << setw(5) << getNumCompute()    << "    " << setw(15) << getComputeTime() << "    "
        << setw(15) << getComputeFlops() << "    "
        << setw(15) << (getComputeTime() != 0.0 ? getComputeFlops() / getComputeTime() * 1.0e-6 : 0.0) << endl;
    out << setw(12) << "apply()" << setw(5) << getNumApply()    << "    " << setw(15) << getApplyTime() << "    "
        << setw(15) << getApplyFlops() << "    "
        << setw(15) << (getApplyTime() != 0.0 ? getApplyFlops() / getApplyTime() * 1.0e-6 : 0.0) << endl;
    out << "===============================================================================" << endl;
    out << endl;
  }
}


template<class MatrixType,class ContainerType>
void
BlockRelaxation<MatrixType,ContainerType>::getMatDiag () const
{
  // TODO amk: Is this map correct for BlockCRSMatrix?
  if(DiagRCP_ == Teuchos::null) {
    DiagRCP_ = Teuchos::rcp(new vector_type(A_->getDomainMap ()));
    A_->getLocalDiagCopy (*DiagRCP_);
  }
}

}//namespace Ifpack2


#ifdef HAVE_IFPACK2_EXPLICIT_INSTANTIATION

// For ETI
#include "Ifpack2_DenseContainer_decl.hpp"
#include "Ifpack2_SparseContainer_decl.hpp"
#include "Ifpack2_TriDiContainer_decl.hpp"
#include "Ifpack2_BandedContainer_decl.hpp"
#include "Ifpack2_ILUT_decl.hpp"
#ifdef HAVE_IFPACK2_AMESOS2
#include "Ifpack2_Details_Amesos2Wrapper.hpp"
#endif

// There's no need to instantiate for CrsMatrix too.  All Ifpack2
// preconditioners can and should do dynamic casts if they need a type
// more specific than RowMatrix.

#ifdef HAVE_IFPACK2_AMESOS2
#define IFPACK2_BLOCKRELAXATION_AMESOS2_INSTANT(S,LO,GO,N) \
  template \
  class Ifpack2::BlockRelaxation<      \
    Tpetra::RowMatrix<S, LO, GO, N>, \
    Ifpack2::SparseContainer<       \
      Tpetra::RowMatrix<S, LO, GO, N>, \
      Ifpack2::Details::Amesos2Wrapper<Tpetra::RowMatrix<S,LO,GO,N> > > >;
#else
#define IFPACK2_BLOCKRELAXATION_AMESOS2_INSTANT(S,LO,GO,N) /* */
#endif

#define IFPACK2_BLOCKRELAXATION_INSTANT(S,LO,GO,N) \
  template \
  class Ifpack2::BlockRelaxation<      \
    Tpetra::RowMatrix<S, LO, GO, N>, \
    Ifpack2::SparseContainer<       \
      Tpetra::RowMatrix<S, LO, GO, N>, \
      Ifpack2::ILUT< ::Tpetra::RowMatrix<S,LO,GO,N> > > >; \
  template \
  class Ifpack2::BlockRelaxation<      \
    Tpetra::RowMatrix<S, LO, GO, N>, \
    Ifpack2::DenseContainer<        \
      Tpetra::RowMatrix<S, LO, GO, N>, \
      S > >; \
  template \
  class Ifpack2::BlockRelaxation<      \
    Tpetra::RowMatrix<S, LO, GO, N>, \
    Ifpack2::TriDiContainer<        \
      Tpetra::RowMatrix<S, LO, GO, N>, \
      S > >; \
  template \
  class Ifpack2::BlockRelaxation<      \
    Tpetra::RowMatrix<S, LO, GO, N>, \
    Ifpack2::BandedContainer<        \
      Tpetra::RowMatrix<S, LO, GO, N>, \
      S > >; \
  IFPACK2_BLOCKRELAXATION_AMESOS2_INSTANT(S,LO,GO,N)

#endif // HAVE_IFPACK2_EXPLICIT_INSTANTIATION

#endif // IFPACK2_BLOCKRELAXATION_DEF_HPP
