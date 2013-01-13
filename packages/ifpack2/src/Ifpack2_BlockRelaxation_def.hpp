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

#ifndef IFPACK2_BLOCKRELAXATION_DEF_HPP
#define IFPACK2_BLOCKRELAXATION_DEF_HPP

#include "Ifpack2_BlockRelaxation_decl.hpp"
#include "Ifpack2_LinearPartitioner_decl.hpp"

namespace Ifpack2 {

//==========================================================================
template<class MatrixType,class ContainerType>
BlockRelaxation<MatrixType,ContainerType>::BlockRelaxation(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A)
: A_(A),
  Time_( Teuchos::rcp( new Teuchos::Time("Ifpack2::BlockRelaxation") ) ),
  OverlapLevel_(0),
  PartitionerType_("linear"),
  NumSweeps_(1),
  PrecType_(Ifpack2::JACOBI),
  MinDiagonalValue_(0.0),
  DampingFactor_(1.0),
  IsParallel_(false),
  ZeroStartingSolution_(true),
  DoBackwardGS_(false),
  Condest_(-1.0),
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
  NumMyRows_(0),
  NumGlobalRows_(0),
  NumGlobalNonzeros_(0)
{
  TEUCHOS_TEST_FOR_EXCEPTION(A_ == Teuchos::null, std::runtime_error, 
      Teuchos::typeName(*this) << "::BlockRelaxation(): input matrix reference was null.");
}

//==========================================================================
template<class MatrixType,class ContainerType>
BlockRelaxation<MatrixType,ContainerType>::~BlockRelaxation() {
}

//==========================================================================
template<class MatrixType,class ContainerType>
void BlockRelaxation<MatrixType,ContainerType>::setParameters(const Teuchos::ParameterList& List)
{
  Teuchos::ParameterList validparams;
  Ifpack2::getValidParameters(validparams);
  List.validateParameters(validparams);

  std::string PT;
  if (PrecType_ == Ifpack2::JACOBI)
    PT = "Jacobi";
  else if (PrecType_ == Ifpack2::GS)
    PT = "Gauss-Seidel";
  else if (PrecType_ == Ifpack2::SGS)
    PT = "Symmetric Gauss-Seidel";

  Ifpack2::getParameter(List, "relaxation: type", PT);

  if (PT == "Jacobi")
    PrecType_ = Ifpack2::JACOBI;
  else if (PT == "Gauss-Seidel")
    PrecType_ = Ifpack2::GS;
  else if (PT == "Symmetric Gauss-Seidel")
    PrecType_ = Ifpack2::SGS;
  else {
    std::ostringstream osstr;
    osstr << "Ifpack2::BlockRelaxation::setParameters: unsupported parameter-value for 'relaxation: type' (" << PT << ")";
    std::string str = osstr.str();
    throw std::runtime_error(str);
  }

  Ifpack2::getParameter(List, "relaxation: sweeps",NumSweeps_);
  Ifpack2::getParameter(List, "relaxation: damping factor", DampingFactor_);
  Ifpack2::getParameter(List, "relaxation: min diagonal value", MinDiagonalValue_);
  Ifpack2::getParameter(List, "relaxation: zero starting solution", ZeroStartingSolution_);
  Ifpack2::getParameter(List, "relaxation: backward mode",DoBackwardGS_);
  Ifpack2::getParameter(List, "partitioner: type",PartitionerType_);
  Ifpack2::getParameter(List, "partitioner: local parts",NumLocalBlocks_);
  Ifpack2::getParameter(List, "partitioner: overlap",OverlapLevel_);

  // check parameters
  if (PrecType_ != Ifpack2::JACOBI)
    OverlapLevel_ = 0;
  if (NumLocalBlocks_ < 0)
    NumLocalBlocks_ = A_->getNodeNumRows() / (-NumLocalBlocks_);
  // other checks are performed in Partitioner_
  

  // NTS: Sanity check to be removed at a later date when Backward mode is enabled
  TEUCHOS_TEST_FOR_EXCEPTION(DoBackwardGS_ == true, std::runtime_error,
			     "Ifpack2::BlockRelaxation:setParameters ERROR: 'relaxation: backward mode' == true is not supported yet.");


  // copy the list as each subblock's constructor will
  // require it later
  List_ = List;
}

//==========================================================================
template<class MatrixType,class ContainerType>
const Teuchos::RCP<const Teuchos::Comm<int> > & 
BlockRelaxation<MatrixType,ContainerType>::getComm() const{
  return A_->getComm();
}

//==========================================================================
template<class MatrixType,class ContainerType>
Teuchos::RCP<const Tpetra::RowMatrix<typename MatrixType::scalar_type,
                                     typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >
BlockRelaxation<MatrixType,ContainerType>::getMatrix() const {
  return(A_);
}

//==========================================================================
template<class MatrixType,class ContainerType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
BlockRelaxation<MatrixType,ContainerType>::getDomainMap() const {
  return A_->getDomainMap();
}

//==========================================================================
template<class MatrixType,class ContainerType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >&
BlockRelaxation<MatrixType,ContainerType>::getRangeMap() const {
  return A_->getRangeMap();
}

//==========================================================================
template<class MatrixType,class ContainerType>
bool BlockRelaxation<MatrixType,ContainerType>::hasTransposeApply() const {
  return true;
}

//==========================================================================
template<class MatrixType,class ContainerType>
int BlockRelaxation<MatrixType,ContainerType>::getNumInitialize() const {
  return(NumInitialize_);
}

//==========================================================================
template<class MatrixType,class ContainerType>
int BlockRelaxation<MatrixType,ContainerType>::getNumCompute() const {
  return(NumCompute_);
}

//==========================================================================
template<class MatrixType,class ContainerType>
int BlockRelaxation<MatrixType,ContainerType>::getNumApply() const {
  return(NumApply_);
}

//==========================================================================
template<class MatrixType,class ContainerType>
double BlockRelaxation<MatrixType,ContainerType>::getInitializeTime() const {
  return(InitializeTime_);
}

//==========================================================================
template<class MatrixType,class ContainerType>
double BlockRelaxation<MatrixType,ContainerType>::getComputeTime() const {
  return(ComputeTime_);
}

//==========================================================================
template<class MatrixType,class ContainerType>
double BlockRelaxation<MatrixType,ContainerType>::getApplyTime() const {
  return(ApplyTime_);
}

//==========================================================================
template<class MatrixType,class ContainerType>
double BlockRelaxation<MatrixType,ContainerType>::getComputeFlops() const {
  return(ComputeFlops_);
}

//==========================================================================
template<class MatrixType,class ContainerType>
double BlockRelaxation<MatrixType,ContainerType>::getApplyFlops() const {
  return(ApplyFlops_);
}

//==========================================================================
template<class MatrixType,class ContainerType>
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType
BlockRelaxation<MatrixType,ContainerType>::getCondEst() const
{
  return(Condest_);
}

//==========================================================================
template<class MatrixType,class ContainerType>
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType
BlockRelaxation<MatrixType,ContainerType>::computeCondEst(
                     CondestType CT,
                     typename MatrixType::local_ordinal_type MaxIters, 
                     magnitudeType Tol,
     const Teuchos::Ptr<const Tpetra::RowMatrix<typename MatrixType::scalar_type,
                                                typename MatrixType::local_ordinal_type,
                                                typename MatrixType::global_ordinal_type,
                                                typename MatrixType::node_type> > &matrix)
{
  if (!isComputed()) // cannot compute right now
    return(-1.0);

  // always compute it. Call Condest() with no parameters to get
  // the previous estimate.
  Condest_ = Ifpack2::Condest(*this, CT, MaxIters, Tol, matrix);

  return(Condest_);
}

//==========================================================================
template<class MatrixType,class ContainerType>
void BlockRelaxation<MatrixType,ContainerType>::apply(
          const Tpetra::MultiVector<typename MatrixType::scalar_type,
                                    typename MatrixType::local_ordinal_type,
                                    typename MatrixType::global_ordinal_type,
                                    typename MatrixType::node_type>& X,
                Tpetra::MultiVector<typename MatrixType::scalar_type,
                                    typename MatrixType::local_ordinal_type,
                                    typename MatrixType::global_ordinal_type,
                                    typename MatrixType::node_type>& Y,
                Teuchos::ETransp mode,
                 Scalar alpha,
                 Scalar beta) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(isComputed() == false, std::runtime_error,
     "Ifpack2::BlockRelaxation::apply ERROR: isComputed() must be true prior to calling apply.");

  TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
     "Ifpack2::BlockRelaxation::apply ERROR: X.getNumVectors() != Y.getNumVectors().");

  TEUCHOS_TEST_FOR_EXCEPTION(mode != Teuchos::NO_TRANS, std::runtime_error,
			     "Ifpack2::BlockRelaxation::apply ERORR: transpose modes not supported.");

  Time_->start(true);

  // If X and Y are pointing to the same memory location,
  // we need to create an auxiliary vector, Xcopy
  Teuchos::RCP< const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Xcopy;
  if (X.getLocalMV().getValues() == Y.getLocalMV().getValues())
    Xcopy = Teuchos::rcp( new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(X) );
  else
    Xcopy = Teuchos::rcp( &X, false );

  if (ZeroStartingSolution_)
    Y.putScalar(0.0);

  // Flops are updated in each of the following.
  switch (PrecType_) {
  case Ifpack2::JACOBI:
    ApplyInverseJacobi(*Xcopy,Y);
    break;
  case Ifpack2::GS:
    ApplyInverseGS(*Xcopy,Y);
    break;
  case Ifpack2::SGS:
    ApplyInverseSGS(*Xcopy,Y);
    break;
  default:
    throw std::runtime_error("Ifpack2::BlockRelaxation::apply internal logic error.");
  }

  ++NumApply_;
  Time_->stop();
  ApplyTime_ += Time_->totalElapsedTime();
}

//==========================================================================
template<class MatrixType,class ContainerType>
void BlockRelaxation<MatrixType,ContainerType>::applyMat(
          const Tpetra::MultiVector<typename MatrixType::scalar_type,
                                    typename MatrixType::local_ordinal_type,
                                    typename MatrixType::global_ordinal_type,
                                    typename MatrixType::node_type>& X,
                Tpetra::MultiVector<typename MatrixType::scalar_type,
                                    typename MatrixType::local_ordinal_type,
                                    typename MatrixType::global_ordinal_type,
                                    typename MatrixType::node_type>& Y,
             Teuchos::ETransp mode) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(isComputed() == false, std::runtime_error,
     "Ifpack2::BlockRelaxation::applyMat() ERROR: isComputed() must be true prior to calling applyMat().");
  TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
     "Ifpack2::BlockRelaxation::applyMat() ERROR: X.getNumVectors() != Y.getNumVectors().");
  A_->apply(X, Y, mode);
}

//==========================================================================
template<class MatrixType,class ContainerType>
void BlockRelaxation<MatrixType,ContainerType>::initialize() {
  IsInitialized_ = false;

  TEUCHOS_TEST_FOR_EXCEPTION(A_ == Teuchos::null, std::runtime_error,
    "Ifpack2::BlockRelaxation::Initialize ERROR, Matrix is NULL");

  Time_->start(true);

  NumMyRows_         = A_->getNodeNumRows();
  NumGlobalRows_     = A_->getGlobalNumRows();
  NumGlobalNonzeros_ = A_->getGlobalNumEntries();

  // NTS: Will need to add support for Zoltan2 partitions later
  // Also, will need a better way of handling the Graph typing issue.  Especially with ordinal types w.r.t the container

  if (PartitionerType_ == "linear")
    Partitioner_ = Teuchos::rcp( new Ifpack2::LinearPartitioner<Tpetra::RowGraph<LocalOrdinal,GlobalOrdinal,Node> >(A_->getGraph()) );
  else
    TEUCHOS_TEST_FOR_EXCEPTION(0, std::runtime_error,"Ifpack2::BlockRelaxation::initialize, invalid partitioner type.");

  // need to partition the graph of A
  Partitioner_->setParameters(List_);
  Partitioner_->compute();

  // get actual number of partitions
  NumLocalBlocks_ = Partitioner_->numLocalParts();
  
  // Note: Unlike Ifpack, we'll punt on computing W_ until compute(), which is where
  // we assume that the type of relaxation has been chosen.
  
  if (A_->getComm()->getSize() != 1)
    IsParallel_ = true;
  else
    IsParallel_ = false;

  ++NumInitialize_;
  Time_->stop();
  InitializeTime_ += Time_->totalElapsedTime();
  IsInitialized_ = true;
}

//==========================================================================
template<class MatrixType,class ContainerType>
void BlockRelaxation<MatrixType,ContainerType>::compute()
{
  if (!isInitialized()) {
    initialize();
  }

  Time_->start(true);

  // reset values
  IsComputed_ = false;
  Condest_ = -1.0;

  TEUCHOS_TEST_FOR_EXCEPTION(NumSweeps_ < 0, std::runtime_error,
    "Ifpack2::BlockRelaxation::compute, NumSweeps_ must be >= 0");

  // Extract the submatrices
  ExtractSubmatrices();

  // Compute the weight vector if we're doing overlapped Jacobi (and only if we're doing overlapped Jacobi).
  if (PrecType_ == Ifpack2::JACOBI && OverlapLevel_ > 0) {
    // weight of each vertex
    W_ = Teuchos::rcp( new Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(A_->getRowMap()) );
    W_->putScalar(Teuchos::ScalarTraits<Scalar>::zero());
    Teuchos::ArrayRCP<Scalar > w_ptr = W_->getDataNonConst(0);

    for (LocalOrdinal i = 0 ; i < NumLocalBlocks_ ; ++i) {    
      for (size_t j = 0 ; j < Partitioner_->numRowsInPart(i) ; ++j) {
	int LID = (*Partitioner_)(i,j);
	w_ptr[LID]+= Teuchos::ScalarTraits<Scalar>::one();
      }
    }
    W_->reciprocal(*W_);
  }

  // We need to import data from external processors. Here I create a
  // Tpetra::Import object if needed (stealing from A_ if possible) 
  // Marzio's comment:
  // Note that I am doing some strange stuff to set the components of Y
  // from Y2 (to save some time).
  //
  if (IsParallel_ && ((PrecType_ == Ifpack2::GS) || (PrecType_ == Ifpack2::SGS))) {
    Importer_=A_->getGraph()->getImporter();
    if(Importer_==Teuchos::null)
      Importer_ = Teuchos::rcp( new Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node>(A_->getDomainMap(),
										    A_->getColMap()) );

    TEUCHOS_TEST_FOR_EXCEPTION(Importer_ == Teuchos::null, std::runtime_error,
      "Ifpack2::BlockRelaxation::compute ERROR failed to create Importer_");
  }

  ++NumCompute_;
  Time_->stop();
  ComputeTime_ += Time_->totalElapsedTime();
  IsComputed_ = true;
}

//==============================================================================
template<class MatrixType,class ContainerType>
void BlockRelaxation<MatrixType,ContainerType>::ExtractSubmatrices()
{
  TEUCHOS_TEST_FOR_EXCEPTION(Partitioner_==Teuchos::null, std::runtime_error,
			     "Ifpack2::BlockRelaxation::ExtractSubmatrices, partitioner is null.");

  NumLocalBlocks_ = Partitioner_->numLocalParts();

  Containers_.resize(NumLocalBlocks_);

  for (LocalOrdinal i = 0 ; i < NumLocalBlocks_ ; ++i) {
    size_t rows = Partitioner_->numRowsInPart(i);
    Containers_[i] = Teuchos::rcp( new ContainerType(rows) );
    TEUCHOS_TEST_FOR_EXCEPTION(Containers_[i]==Teuchos::null, std::runtime_error,
			     "Ifpack2::BlockRelaxation::ExtractSubmatrices, container consturctor failed.");
    
    Containers_[i]->setParameters(List_);
    Containers_[i]->initialize();
    // flops in initialize() will be computed on-the-fly in method initializeFlops().

    // set "global" ID of each partitioner row
    for (size_t j = 0 ; j < rows ; ++j) {
      Containers_[i]->ID(j) = (*Partitioner_)(i,j);
    }

    Containers_[i]->compute(A_);
    // flops in compute() will be computed on-the-fly in method computeFlops().
  }
}



//==========================================================================
template<class MatrixType,class ContainerType>
void BlockRelaxation<MatrixType,ContainerType>::ApplyInverseJacobi(
        const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
              Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  size_t NumVectors = X.getNumVectors();
  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> AY( Y.getMap(),NumVectors );
  
  // Initial matvec not needed
  int starting_iteration=0;
  if(ZeroStartingSolution_) {
    DoJacobi(X,Y);
    starting_iteration=1;
  }

  for (int j = starting_iteration; j < NumSweeps_ ; j++) {       
    applyMat(Y,AY);
    AY.update(1.0,X,-1.0);
    DoJacobi(AY,Y);

    // Flops for matrix apply & update
    ApplyFlops_ += NumVectors * (2 * NumGlobalNonzeros_ + 2 * NumGlobalRows_);
  }

}


//==========================================================================
template<class MatrixType,class ContainerType>
void BlockRelaxation<MatrixType,ContainerType>::DoJacobi(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  size_t NumVectors = X.getNumVectors();
  Scalar one=Teuchos::ScalarTraits<Scalar>::one();
  // Note: Flop counts copied naively from Ifpack.

  if (OverlapLevel_ == 0) {
    // Non-overlapping Jacobi
    for (LocalOrdinal i = 0 ; i < NumLocalBlocks_ ; i++) {     
      // may happen that a partition is empty
      if (Containers_[i]->getNumRows() == 0) continue;
      Containers_[i]->apply(X,Y,Teuchos::NO_TRANS,DampingFactor_,one);     
      ApplyFlops_ += NumVectors * 2 * NumGlobalRows_;
    }
  }
  else {
    // Overlapping Jacobi
    for (LocalOrdinal i = 0 ; i < NumLocalBlocks_ ; i++) {
      // may happen that a partition is empty
      if (Containers_[i]->getNumRows() == 0) continue;
      Containers_[i]->weightedApply(X,Y,*W_,Teuchos::NO_TRANS,DampingFactor_,one);
      // NOTE: do not count (for simplicity) the flops due to overlapping rows
      ApplyFlops_ += NumVectors * 4 * NumGlobalRows_;
    }
  }
}

//==========================================================================
template<class MatrixType,class ContainerType>
void BlockRelaxation<MatrixType,ContainerType>::ApplyInverseGS(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
							       Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{

  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>  Xcopy(X);
  for (int j = 0; j < NumSweeps_ ; j++) {
    DoGaussSeidel(Xcopy,Y);
    if(j != NumSweeps_-1)
      Xcopy=X;
  }
}


//==============================================================================
template<class MatrixType,class ContainerType>
void BlockRelaxation<MatrixType,ContainerType>::DoGaussSeidel(Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
							      Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  // Note: Flop counts copied naively from Ifpack.
  
  // allocations
  Scalar one=Teuchos::ScalarTraits<Scalar>::one();
  int Length = A_->getNodeMaxNumRowEntries(); 
  int NumVectors = X.getNumVectors();
  Teuchos::Array<Scalar>         Values;
  Teuchos::Array<LocalOrdinal>   Indices;
  Values.resize(Length);
  Indices.resize(Length);

  // an additonal vector is needed by parallel computations
  // (note that applications through Ifpack2_AdditiveSchwarz
  // are always seen are serial)
  Teuchos::RCP< Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node > > Y2;
  if (IsParallel_)
    Y2 = Teuchos::rcp( new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(Importer_->getTargetMap(), NumVectors) );
  else
    Y2 = Teuchos::rcp( &Y, false );


  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >       x_ptr   = X.get2dViewNonConst();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >       y_ptr   = Y.get2dViewNonConst();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >       y2_ptr  = Y2->get2dViewNonConst();

  // data exchange is here, once per sweep
  if (IsParallel_)  Y2->doImport(Y,*Importer_,Tpetra::INSERT);

  for (LocalOrdinal i = 0 ; i < NumLocalBlocks_ ; i++) {
    // may happen that a partition is empty
    if (Containers_[i]->getNumRows() == 0) continue;
    LocalOrdinal LID;

    // update from previous block
    for (size_t j = 0 ; j < Containers_[i]->getNumRows(); j++) {
      LID = Containers_[i]->ID(j);
      size_t NumEntries;
      A_->getLocalRowCopy(LID,Indices(),Values(),NumEntries);

      for (size_t k = 0 ; k < NumEntries ; k++) {
        LocalOrdinal col = Indices[k];
	for (int kk = 0 ; kk < NumVectors ; kk++)
	  x_ptr[kk][LID] -= Values[k] * y2_ptr[kk][col];	
      }
    }
    // solve with this block
    // Note: I'm abusing the ordering information, knowing that X/Y and Y2 have the same ordering for on-proc unknowns.
    // Note: Add flop counts for inverse apply
    Containers_[i]->apply(X,*Y2,Teuchos::NO_TRANS,DampingFactor_,one);
    
    // operations for all getrow's
    ApplyFlops_ += NumVectors * (2 * NumGlobalNonzeros_ + 2 * NumGlobalRows_);    
  } // end for NumLocalBlocks_

  // Attention: this is delicate... Not all combinations
  // of Y2 and Y will always work (tough for ML it should be ok)
  if (IsParallel_)
    for (int m = 0 ; m < NumVectors ; ++m) 
      for (size_t i = 0 ; i < NumMyRows_ ; ++i)
        y_ptr[m][i] = y2_ptr[m][i];

}

//==========================================================================
template<class MatrixType,class ContainerType>
void BlockRelaxation<MatrixType,ContainerType>::ApplyInverseSGS(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
								Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  
  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>  Xcopy(X);
  for (int j = 0; j < NumSweeps_ ; j++) {
    DoSGS(X,Xcopy,Y);
    if(j != NumSweeps_-1)
      Xcopy=X;
  }
}

//==========================================================================
template<class MatrixType,class ContainerType>
void BlockRelaxation<MatrixType,ContainerType>::DoSGS(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
						      Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Xcopy, 
						      Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  Scalar one=Teuchos::ScalarTraits<Scalar>::one();
  int Length = A_->getNodeMaxNumRowEntries(); 
  int NumVectors = X.getNumVectors();
  Teuchos::Array<Scalar>         Values;
  Teuchos::Array<LocalOrdinal>   Indices;
  Values.resize(Length);
  Indices.resize(Length);

  // an additonal vector is needed by parallel computations
  // (note that applications through Ifpack2_AdditiveSchwarz
  // are always seen are serial)
  Teuchos::RCP< Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node > > Y2;
  if (IsParallel_)
    Y2 = Teuchos::rcp( new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(Importer_->getTargetMap(), NumVectors) );
  else
    Y2 = Teuchos::rcp( &Y, false );

  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > x_ptr       = X.get2dView();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >       xcopy_ptr   = Xcopy.get2dViewNonConst();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >       y_ptr       = Y.get2dViewNonConst();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >       y2_ptr      = Y2->get2dViewNonConst();

  // data exchange is here, once per sweep
  if (IsParallel_)  Y2->doImport(Y,*Importer_,Tpetra::INSERT);

  // Forward Sweep
  for(LocalOrdinal i = 0 ; i < NumLocalBlocks_ ; i++) {
    // may happen that a partition is empty
    if (Containers_[i]->getNumRows() == 0) continue;
    LocalOrdinal LID;

    // update from previous block
    for(size_t j = 0 ; j < Containers_[i]->getNumRows(); j++) {
      LID = Containers_[i]->ID(j);
      size_t NumEntries;
      A_->getLocalRowCopy(LID,Indices(),Values(),NumEntries);

      for (size_t k = 0 ; k < NumEntries ; k++) {
        LocalOrdinal col = Indices[k];
        for (int kk = 0 ; kk < NumVectors ; kk++) 
          xcopy_ptr[kk][LID] -= Values[k] * y2_ptr[kk][col];
      }
    }
    // solve with this block
    // Note: I'm abusing the ordering information, knowing that X/Y and Y2 have the same ordering for on-proc unknowns.
    // Note: Add flop counts for inverse apply
    Containers_[i]->apply(Xcopy,*Y2,Teuchos::NO_TRANS,DampingFactor_,one);

    // operations for all getrow's
    ApplyFlops_ += NumVectors * (2 * NumGlobalNonzeros_ + 2 * NumGlobalRows_);
  }// end forward sweep

  // Reverse Sweep
  Xcopy = X;
  for(LocalOrdinal i = NumLocalBlocks_-1; i >=0 ; i--) {
    // may happen that a partition is empty
    if (Containers_[i]->getNumRows() == 0) continue;
    LocalOrdinal LID;

    // update from previous block
    for(size_t j = 0 ; j < Containers_[i]->getNumRows(); j++) {
      LID = Containers_[i]->ID(j);
      size_t NumEntries;
      A_->getLocalRowCopy(LID,Indices(),Values(),NumEntries);

      for (size_t k = 0 ; k < NumEntries ; k++) {
        LocalOrdinal col = Indices[k];
        for (int kk = 0 ; kk < NumVectors ; kk++) 
          xcopy_ptr[kk][LID] -= Values[k] * y2_ptr[kk][col];
      }
    }
    
    // solve with this block
    // Note: I'm abusing the ordering information, knowing that X/Y and Y2 have the same ordering for on-proc unknowns.
    // Note: Add flop counts for inverse apply
    Containers_[i]->apply(Xcopy,*Y2,Teuchos::NO_TRANS,DampingFactor_,one);

    // operations for all getrow's
    ApplyFlops_ += NumVectors * (2 * NumGlobalNonzeros_ + 2 * NumGlobalRows_);
  } //end reverse sweep


  // Attention: this is delicate... Not all combinations
  // of Y2 and Y will always work (tough for ML it should be ok)
  if (IsParallel_)
    for (int m = 0 ; m < NumVectors ; ++m) 
      for (size_t i = 0 ; i < NumMyRows_ ; ++i)
        y_ptr[m][i] = y2_ptr[m][i];
}

//==========================================================================
template<class MatrixType,class ContainerType>
std::string BlockRelaxation<MatrixType,ContainerType>::description() const {
  std::ostringstream oss;
  oss << Teuchos::Describable::description();
  if (isInitialized()) {
    if (isComputed()) {
      oss << "{status = initialized, computed";
    }
    else {
      oss << "{status = initialized, not computed";
    }
  }
  else {
    oss << "{status = not initialized, not computed";
  }
  //
  if (PrecType_ == Ifpack2::JACOBI)   oss << "Type = Block Jacobi, " << std::endl;
  else if (PrecType_ == Ifpack2::GS)  oss << "Type = Block Gauss-Seidel, " << std::endl;
  else if (PrecType_ == Ifpack2::SGS) oss << "Type = Block Sym. Gauss-Seidel, " << std::endl;
  //
  oss << ", global rows = " << A_->getGlobalNumRows()
      << ", global cols = " << A_->getGlobalNumCols();

  oss << "}";
  return oss.str();
}

//==========================================================================
template<class MatrixType,class ContainerType>
void BlockRelaxation<MatrixType,ContainerType>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
  using std::endl;
  using std::setw;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;
  Teuchos::EVerbosityLevel vl = verbLevel;
  if (vl == VERB_DEFAULT) vl = VERB_LOW;
  const int myImageID = A_->getComm()->getRank();
  Teuchos::OSTab tab(out);

  //    none: print nothing
  //     low: print O(1) info from node 0
  //  medium: 
  //    high: 
  // extreme: 
  if (vl != VERB_NONE && myImageID == 0) {
    out << this->description() << endl;
    out << endl;
    out << "===============================================================================" << endl;
    out << "Sweeps         = " << NumSweeps_ << endl;
    out << "damping factor = " << DampingFactor_ << endl;
    if (PrecType_ == Ifpack2::GS && DoBackwardGS_) {
      out << "Using backward mode (BGS only)" << endl;
    }
    if   (ZeroStartingSolution_) { out << "Using zero starting solution" << endl; }
    else                         { out << "Using input starting solution" << endl; }
    if   (Condest_ == -1.0) { out << "Condition number estimate       = N/A" << endl; }
    else                    { out << "Condition number estimate       = " << Condest_ << endl; }
    out << endl;
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

}//namespace Ifpack2

#endif // IFPACK2_BLOCKRELAXATION_DEF_HPP

