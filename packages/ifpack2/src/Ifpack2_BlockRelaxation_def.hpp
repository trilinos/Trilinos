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
    W_ = Teuchos::null;

    if (! A.is_null ()) {
      IsParallel_ = (A->getRowMap ()->getComm ()->getSize () != 1);
      Teuchos::RCP<const block_crs_matrix_type> A_bcrs =
        Teuchos::rcp_dynamic_cast<const block_crs_matrix_type> (A);
      hasBlockCrsMatrix_ = !A_bcrs.is_null();
    }
    if (! Container_.is_null ()) {
      Container_->clearBlocks();
    }
    NumLocalBlocks_ = 0;

    A_ = A;
    computeImporter();
  }
}

template<class MatrixType,class ContainerType>
BlockRelaxation<MatrixType,ContainerType>::
BlockRelaxation (const Teuchos::RCP<const row_matrix_type>& A)
:
  Time_ (Teuchos::rcp (new Teuchos::Time ("Ifpack2::BlockRelaxation"))),
  Container_ (Teuchos::null),
  Partitioner_ (Teuchos::null),
  PartitionerType_ ("linear"),
  NumSweeps_ (1),
  NumLocalBlocks_(0),
  containerType_ ("Dense"),
  PrecType_ (Ifpack2::Details::JACOBI),
  ZeroStartingSolution_ (true),
  hasBlockCrsMatrix_ (false),
  DoBackwardGS_ (false),
  OverlapLevel_ (0),
  DampingFactor_ (STS::one ()),
  IsInitialized_ (false),
  IsComputed_ (false),
  NumInitialize_ (0),
  NumCompute_ (0),
  NumApply_ (0),
  InitializeTime_ (0.0),
  ComputeTime_ (0.0),
  NumLocalRows_ (0),
  NumGlobalRows_ (0),
  NumGlobalNonzeros_ (0),
  W_ (Teuchos::null),
  Importer_ (Teuchos::null)
{
  setMatrix(A);
}

template<class MatrixType,class ContainerType>
BlockRelaxation<MatrixType,ContainerType>::
~BlockRelaxation ()
{}

template<class MatrixType,class ContainerType>
Teuchos::RCP<const Teuchos::ParameterList>
BlockRelaxation<MatrixType,ContainerType>::
getValidParameters () const
{
  Teuchos::RCP<Teuchos::ParameterList> validParams = Teuchos::parameterList ("Ifpack2::BlockRelaxation");

  validParams->set("relaxation: container", "TriDi");
  validParams->set("relaxation: backward mode",false);
  validParams->set("relaxation: type", "Jacobi");
  validParams->set("relaxation: sweeps", (int)1);
  validParams->set("relaxation: damping factor", STS::one());
  validParams->set("relaxation: zero starting solution", true);
  validParams->set("schwarz: compute condest", false); // mfh 24 Mar 2015: for backwards compatibility ONLY
  validParams->set("schwarz: combine mode", "ZERO"); // use string mode for this
  validParams->set("schwarz: use reordering", true);
  validParams->set("schwarz: filter singletons", false);
  validParams->set("schwarz: overlap level", (int)0);
  validParams->set("partitioner: type", "greedy");
  validParams->set("partitioner: local parts", (int)1);
  validParams->set("partitioner: overlap", (int)0);
  Teuchos::Array<Teuchos::ArrayRCP<int>> tmp0;
  validParams->set("partitioner: parts", tmp0);
  validParams->set("partitioner: maintain sparsity", false);

  Teuchos::ParameterList dummyList;
  validParams->set("Amesos2",dummyList);
  validParams->sublist("Amesos2").disableRecursiveValidation();
  validParams->set("Amesos2 solver name", "KLU2");

  Teuchos::ArrayRCP<int> tmp;
  validParams->set("partitioner: map", tmp);

  validParams->set("partitioner: line detection threshold",(double)0.0);
  validParams->set("partitioner: PDE equations",(int)1);
  Teuchos::RCP<Tpetra::MultiVector<double,
                                   typename MatrixType::local_ordinal_type,
                                   typename MatrixType::global_ordinal_type,
                                   typename MatrixType::node_type> > dummy;
  validParams->set("partitioner: coordinates",dummy);

  return validParams;
}

template<class MatrixType,class ContainerType>
void
BlockRelaxation<MatrixType,ContainerType>::
setParameters (const Teuchos::ParameterList& List)
{
  // Note that the validation process does not change List.
  Teuchos::RCP<const Teuchos::ParameterList> validparams;
  validparams = this->getValidParameters();
  List.validateParameters (*validparams);

  if (List.isParameter ("relaxation: container")) {
    // If the container type isn't a string, this will throw, but it
    // rightfully should.

    // If its value does not match the currently registered Container types,
    // the ContainerFactory will throw with an informative message.
    containerType_ = List.get<std::string> ("relaxation: container");
  }

  if (List.isParameter ("relaxation: type")) {
    const std::string relaxType = List.get<std::string> ("relaxation: type");

    if (relaxType == "Jacobi") {
      PrecType_ = Ifpack2::Details::JACOBI;
    }
    else if (relaxType == "MT Split Jacobi") {
      PrecType_ = Ifpack2::Details::MTSPLITJACOBI;
    }
    else if (relaxType == "Gauss-Seidel") {
      PrecType_ = Ifpack2::Details::GS;
    }
    else if (relaxType == "Symmetric Gauss-Seidel") {
      PrecType_ = Ifpack2::Details::SGS;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::invalid_argument, "Ifpack2::BlockRelaxation::"
         "setParameters: Invalid parameter value \"" << relaxType
         << "\" for parameter \"relaxation: type\".");
    }
  }

  if (List.isParameter ("relaxation: sweeps")) {
    NumSweeps_ = List.get<int> ("relaxation: sweeps");
  }

  // Users may specify this as a floating-point literal.  In that
  // case, it may have type double, even if scalar_type is something
  // else.  We take care to try the various types that make sense.
  if (List.isParameter ("relaxation: damping factor")) {
    if (List.isType<double> ("relaxation: damping factor")) {
      const double dampingFactor =
        List.get<double> ("relaxation: damping factor");
      DampingFactor_ = static_cast<scalar_type> (dampingFactor);
    }
    else if (List.isType<scalar_type> ("relaxation: damping factor")) {
      DampingFactor_ = List.get<scalar_type> ("relaxation: damping factor");
    }
    else if (List.isType<magnitude_type> ("relaxation: damping factor")) {
      const magnitude_type dampingFactor =
        List.get<magnitude_type> ("relaxation: damping factor");
      DampingFactor_ = static_cast<scalar_type> (dampingFactor);
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::invalid_argument, "Ifpack2::BlockRelaxation::"
         "setParameters: Parameter \"relaxation: damping factor\" "
         "has the wrong type.");
    }
  }

  if (List.isParameter ("relaxation: zero starting solution")) {
    ZeroStartingSolution_ = List.get<bool> ("relaxation: zero starting solution");
  }

  if (List.isParameter ("relaxation: backward mode")) {
    DoBackwardGS_ = List.get<bool> ("relaxation: backward mode");
  }

  if (List.isParameter ("partitioner: type")) {
    PartitionerType_ = List.get<std::string> ("partitioner: type");
  }

  // Users may specify this as an int literal, so we need to allow
  // both int and local_ordinal_type (not necessarily same as int).
  if (List.isParameter ("partitioner: local parts")) {
    if (List.isType<local_ordinal_type> ("partitioner: local parts")) {
      NumLocalBlocks_ = List.get<local_ordinal_type> ("partitioner: local parts");
    }
    else if (! std::is_same<int, local_ordinal_type>::value &&
             List.isType<int> ("partitioner: local parts")) {
      NumLocalBlocks_ = List.get<int> ("partitioner: local parts");
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::invalid_argument, "Ifpack2::BlockRelaxation::"
         "setParameters: Parameter \"partitioner: local parts\" "
         "has the wrong type.");
    }
  }

  if (List.isParameter ("partitioner: overlap level")) {
    if (List.isType<int> ("partitioner: overlap level")) {
      OverlapLevel_ = List.get<int> ("partitioner: overlap level");
    }
    else if (! std::is_same<int, local_ordinal_type>::value &&
             List.isType<local_ordinal_type> ("partitioner: overlap level")) {
      OverlapLevel_ = List.get<local_ordinal_type> ("partitioner: overlap level");
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::invalid_argument, "Ifpack2::BlockRelaxation::"
         "setParameters: Parameter \"partitioner: overlap level\" "
         "has the wrong type.");
    }
  }

  std::string defaultContainer = "Dense";
  if(std::is_same<ContainerType, Container<MatrixType> >::value)
  {
    //Generic container template parameter, container type specified in List
    Ifpack2::getParameter(List, "relaxation: container", defaultContainer);
  }
  // check parameters
  if (PrecType_ != Ifpack2::Details::JACOBI) {
    OverlapLevel_ = 0;
  }
  if (NumLocalBlocks_ < static_cast<local_ordinal_type> (0)) {
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
size_t BlockRelaxation<MatrixType,ContainerType>::getNodeSmootherComplexity() const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::BlockRelaxation::getNodeSmootherComplexity: "
    "The input matrix A is null.  Please call setMatrix() with a nonnull "
    "input matrix, then call compute(), before calling this method.");
  // Relaxation methods cost roughly one apply + one block-diagonal inverse per iteration
  // NOTE: This approximates all blocks as dense, which may overstate the cost if you have a sparse (or banded) container.
  size_t block_nnz = 0;
  for (local_ordinal_type i = 0; i < NumLocalBlocks_; ++i)
    block_nnz += Partitioner_->numRowsInPart(i) *Partitioner_->numRowsInPart(i);

  return block_nnz + A_->getNodeNumEntries();
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
    if (X_lcl_host.data () == Y_lcl_host.data ()) {
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
  case Ifpack2::Details::MTSPLITJACOBI:
    Container_->applyInverseJacobi(*X_copy, Y, ZeroStartingSolution_, NumSweeps_);
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
  hasBlockCrsMatrix_ = !A_bcrs.is_null();
  if (A_bcrs.is_null ()) {
    hasBlockCrsMatrix_ = false;
  }
  else {
    hasBlockCrsMatrix_ = true;
  }

  IsInitialized_ = false;
  Time_->start (true);

  NumLocalRows_      = A_->getNodeNumRows ();
  NumGlobalRows_     = A_->getGlobalNumRows ();
  NumGlobalNonzeros_ = A_->getGlobalNumEntries ();

  // NTS: Will need to add support for Zoltan2 partitions later Also,
  // will need a better way of handling the Graph typing issue.
  // Especially with ordinal types w.r.t the container.
  Partitioner_ = Teuchos::null;

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

  // We should have checked for this in setParameters(), so it's a
  // logic_error, not an invalid_argument or runtime_error.
  TEUCHOS_TEST_FOR_EXCEPTION
    (NumSweeps_ < 0, std::logic_error, "Ifpack2::BlockRelaxation::initialize: "
     "NumSweeps_ = " << NumSweeps_ << " < 0.");

  // Extract the submatrices
  ExtractSubmatricesStructure ();

  // Compute the weight vector if we're doing overlapped Jacobi (and
  // only if we're doing overlapped Jacobi).
  if (PrecType_ == Ifpack2::Details::JACOBI && OverlapLevel_ > 0) {
    TEUCHOS_TEST_FOR_EXCEPTION
      (hasBlockCrsMatrix_, std::runtime_error,
       "Ifpack2::BlockRelaxation::initialize: "
       "We do not support overlapped Jacobi yet for Tpetra::BlockCrsMatrix.  Sorry!");

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

  ++NumInitialize_;
  Time_->stop ();
  InitializeTime_ += Time_->totalElapsedTime ();
  IsInitialized_ = true;
}


template<class MatrixType,class ContainerType>
void
BlockRelaxation<MatrixType,ContainerType>::
compute ()
{
  using Teuchos::rcp;

  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error, "Ifpack2::BlockRelaxation::compute: "
     "The matrix is null.  You must call setMatrix() with a nonnull matrix, "
     "then call initialize(), before you may call this method.");

  if (! isInitialized ()) {
    initialize ();
  }

  Time_->start (true);

  // reset values
  IsComputed_ = false;

  Container_->compute();   // compute each block matrix

  ++NumCompute_;
  Time_->stop ();
  ComputeTime_ += Time_->totalElapsedTime();
  IsComputed_ = true;
}

template<class MatrixType, class ContainerType>
void
BlockRelaxation<MatrixType, ContainerType>::
ExtractSubmatricesStructure ()
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (Partitioner_.is_null (), std::runtime_error, "Ifpack2::BlockRelaxation::"
     "ExtractSubmatricesStructure: Partitioner object is null.");

  NumLocalBlocks_ = Partitioner_->numLocalParts ();
  std::string containerType = ContainerType::getName ();
  if (containerType == "Generic") {
    // ContainerType is Container<row_matrix_type>.  Thus, we need to
    // get the container name from the parameter list.  We use "Dense"
    // as the default container type.
    containerType = containerType_;
  }

  Teuchos::Array<Teuchos::Array<local_ordinal_type> > localRows(NumLocalBlocks_);
  for (local_ordinal_type i = 0; i < NumLocalBlocks_; ++i) {
    const size_t numRows = Partitioner_->numRowsInPart (i);

    localRows[i].resize(numRows);
    // Extract a list of the indices of each partitioner row.
    for (size_t j = 0; j < numRows; ++j) {
      localRows[i][j] = (*Partitioner_) (i,j);
    }
  }
  Container_ = ContainerFactory<MatrixType>::build(containerType, A_, localRows, Importer_, OverlapLevel_, DampingFactor_);
  Container_->setParameters(List_);
  Container_->initialize();
}

template<class MatrixType,class ContainerType>
void
BlockRelaxation<MatrixType,ContainerType>::
ApplyInverseJacobi (const MV& X, MV& Y) const
{
  const size_t NumVectors = X.getNumVectors ();
  typename ContainerType::HostView XView = X.template getLocalView<Kokkos::HostSpace>();
  typename ContainerType::HostView YView = Y.template getLocalView<Kokkos::HostSpace>();
  MV AY (Y.getMap (), NumVectors);

  typename ContainerType::HostView AYView = AY.template getLocalView<Kokkos::HostSpace>();
  // Initial matvec not needed
  int starting_iteration = 0;
  if (OverlapLevel_ > 0)
  {
    //Overlapping jacobi, with view of W_
    typename ContainerType::HostView WView = W_->template getLocalView<Kokkos::HostSpace>();
    if(ZeroStartingSolution_) {
      Container_->DoOverlappingJacobi(XView, YView, WView, X.getStride());
      starting_iteration = 1;
    }
    const scalar_type ONE = STS::one();
    for(int j = starting_iteration; j < NumSweeps_; j++)
    {
      applyMat (Y, AY);
      AY.update (ONE, X, -ONE);
      Container_->DoOverlappingJacobi (AYView, YView, WView, X.getStride());
    }
  }
  else
  {
    //Non-overlapping
    if(ZeroStartingSolution_)
    {
      Container_->DoJacobi (XView, YView, X.getStride());
      starting_iteration = 1;
    }
    const scalar_type ONE = STS::one();
    for(int j = starting_iteration; j < NumSweeps_; j++)
    {
      applyMat (Y, AY);
      AY.update (ONE, X, -ONE);
      Container_->DoJacobi (AYView, YView, X.getStride());
    }
  }
}

template<class MatrixType,class ContainerType>
void
BlockRelaxation<MatrixType,ContainerType>::
ApplyInverseGS (const MV& X, MV& Y) const
{
  using Teuchos::Ptr;
  using Teuchos::ptr;
  size_t numVecs = X.getNumVectors();
  //Get view of X (is never modified in this function)
  typename ContainerType::HostView XView = X.template getLocalView<Kokkos::HostSpace>();
  typename ContainerType::HostView YView = Y.template getLocalView<Kokkos::HostSpace>();
  //Pre-import Y, if parallel
  Ptr<MV> Y2;
  bool deleteY2 = false;
  if(IsParallel_)
  {
    Y2 = ptr(new MV(Importer_->getTargetMap(), numVecs));
    deleteY2 = true;
  }
  else
    Y2 = ptr(&Y);
  if(IsParallel_)
  {
    for(int j = 0; j < NumSweeps_; ++j)
    {
      //do import once per sweep
      Y2->doImport(Y, *Importer_, Tpetra::INSERT);
      typename ContainerType::HostView Y2View = Y2->template getLocalView<Kokkos::HostSpace>();
      Container_->DoGaussSeidel(XView, YView, Y2View, X.getStride());
    }
  }
  else
  {
    typename ContainerType::HostView Y2View = Y2->template getLocalView<Kokkos::HostSpace>();
    for(int j = 0; j < NumSweeps_; ++j)
    {
      Container_->DoGaussSeidel(XView, YView, Y2View, X.getStride());
    }
  }
  if(deleteY2)
    delete Y2.get();
}

template<class MatrixType,class ContainerType>
void
BlockRelaxation<MatrixType,ContainerType>::
ApplyInverseSGS (const MV& X, MV& Y) const
{
  using Teuchos::Ptr;
  using Teuchos::ptr;
  //Get view of X (is never modified in this function)
  typename ContainerType::HostView XView = X.template getLocalView<Kokkos::HostSpace>();
  typename ContainerType::HostView YView = Y.template getLocalView<Kokkos::HostSpace>();
  //Pre-import Y, if parallel
  Ptr<MV> Y2;
  bool deleteY2 = false;
  if(IsParallel_)
  {
    Y2 = ptr(new MV(Importer_->getTargetMap(), X.getNumVectors()));
    deleteY2 = true;
  }
  else
    Y2 = ptr(&Y);
  if(IsParallel_)
  {
    for(int j = 0; j < NumSweeps_; ++j)
    {
      //do import once per sweep
      Y2->doImport(Y, *Importer_, Tpetra::INSERT);
      typename ContainerType::HostView Y2View = Y2->template getLocalView<Kokkos::HostSpace>();
      Container_->DoSGS(XView, YView, Y2View, X.getStride());
    }
  }
  else
  {
    typename ContainerType::HostView Y2View = Y2->template getLocalView<Kokkos::HostSpace>();
    for(int j = 0; j < NumSweeps_; ++j)
    {
      Container_->DoSGS(XView, YView, Y2View, X.getStride());
    }
  }
  if(deleteY2)
    delete Y2.get();
}

template<class MatrixType,class ContainerType>
void BlockRelaxation<MatrixType,ContainerType>::computeImporter () const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Ptr;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::Comm;
  using Teuchos::rcp_dynamic_cast;
  if(IsParallel_)
  {
    if(hasBlockCrsMatrix_)
    {
      const RCP<const block_crs_matrix_type> bcrs =
        rcp_dynamic_cast<const block_crs_matrix_type>(A_);
      int bs_ = bcrs->getBlockSize();
      RCP<const map_type> oldDomainMap = A_->getDomainMap();
      RCP<const map_type> oldColMap = A_->getColMap();
      // Because A is a block CRS matrix, import will not do what you think it does
      // We have to construct the correct maps for it
      global_size_t numGlobalElements = oldColMap->getGlobalNumElements() * bs_;
      global_ordinal_type indexBase = oldColMap->getIndexBase();
      RCP<const Comm<int>> comm = oldColMap->getComm();
      ArrayView<const global_ordinal_type> oldColElements = oldColMap->getNodeElementList();
      Array<global_ordinal_type> newColElements(bs_ * oldColElements.size());
      for(int i = 0; i < oldColElements.size(); i++)
      {
        for(int j = 0; j < bs_; j++)
          newColElements[i * bs_ + j] = oldColElements[i] * bs_ + j;
      }
      RCP<map_type> colMap(new map_type(numGlobalElements, newColElements, indexBase, comm));
      // Create the importer
      Importer_ = rcp(new import_type(oldDomainMap, colMap));
    }
    else if(!A_.is_null())
    {
      Importer_ = A_->getGraph()->getImporter();
      if(Importer_.is_null())
        Importer_ = rcp(new import_type(A_->getDomainMap(), A_->getColMap()));
    }
  }
  //otherwise, Importer_ is not needed and is left NULL
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

  out << ", overlap: " << OverlapLevel_;

  out  << ", " << "sweeps: " << NumSweeps_ << ", "
      << "damping factor: " << DampingFactor_ << ", ";

  std::string containerType = ContainerType::getName();
  out << "block container: " << ((containerType == "Generic") ? containerType_ : containerType);

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
  }
}

}//namespace Ifpack2


// Macro that does explicit template instantiation (ETI) for
// Ifpack2::BlockRelaxation.  S, LO, GO, N correspond to the four
// template parameters of Ifpack2::Preconditioner and
// Tpetra::RowMatrix.
//
// We only instantiate for MatrixType = Tpetra::RowMatrix.  There's no
// need to instantiate for Tpetra::CrsMatrix too.  All Ifpack2
// preconditioners can and should do dynamic casts if they need a type
// more specific than RowMatrix.  This keeps build time short and
// library and executable sizes small.

#ifdef HAVE_IFPACK2_EXPLICIT_INSTANTIATION

#define IFPACK2_BLOCKRELAXATION_INSTANT(S,LO,GO,N) \
  template \
  class Ifpack2::BlockRelaxation<      \
    Tpetra::RowMatrix<S, LO, GO, N>,   \
    Ifpack2::Container<Tpetra::RowMatrix<S, LO, GO, N> > >;

#endif // HAVE_IFPACK2_EXPLICIT_INSTANTIATION

#endif // IFPACK2_BLOCKRELAXATION_DEF_HPP
