//@HEADER
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
// 1
// ***********************************************************************
//@HEADER

#ifndef IFPACK2_ADDITIVESCHWARZ_DEF_HPP
#define IFPACK2_ADDITIVESCHWARZ_DEF_HPP

#include "Ifpack2_AdditiveSchwarz_decl.hpp"

#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
#include "Xpetra_RowMatrix.hpp"
#include "Xpetra_TpetraRowMatrix.hpp"
#include "Zoltan2_XpetraRowMatrixInput.hpp"
#include "Zoltan2_OrderingProblem.hpp"
#endif

#include "Ifpack2_Condest.hpp"
#include "Ifpack2_OverlappingRowMatrix_def.hpp"
#include "Ifpack2_LocalFilter_def.hpp"
#include "Ifpack2_ReorderFilter_def.hpp"
#include "Ifpack2_SingletonFilter_def.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif

namespace Ifpack2 {

//==============================================================================
template<class MatrixType,class LocalInverseType>
AdditiveSchwarz<MatrixType,LocalInverseType>::AdditiveSchwarz(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > & Matrix_in, int OverlapLevel_in):
  Matrix_(Matrix_in),
  IsInitialized_(false),
  IsComputed_(false),
  IsOverlapping_(false),
  OverlapLevel_(OverlapLevel_in),
  CombineMode_(Tpetra::ADD),
  Condest_(-1.0),
  ComputeCondest_(true),
  UseReordering_(false),
  UseSubdomain_(false),
  FilterSingletons_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApply_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyTime_(0.0)
{
  
  // Reset Overlap parameters
  if (Matrix_->getComm()->getSize() == 1)
    OverlapLevel_ = 0;
  
  if ((OverlapLevel_ != 0) && (Matrix_->getComm()->getSize() > 1))
    IsOverlapping_ = true;

  // Sets parameters to default values
  Teuchos::ParameterList List_in;
  setParameters(List_in);
}

//==============================================================================
// Destructor
template<class MatrixType,class LocalInverseType>
AdditiveSchwarz<MatrixType,LocalInverseType>::~AdditiveSchwarz()
{
}

//==============================================================================
// Returns the Map associated with the domain of this operator, which must be compatible with X.getMap().
template<class MatrixType,class LocalInverseType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type > > & AdditiveSchwarz<MatrixType,LocalInverseType>::getDomainMap() const 
{ 
  return Matrix_->getDomainMap();
}
  
//==============================================================================
// Returns the Map associated with the range of this operator, which must be compatible with Y.getMap().
template<class MatrixType,class LocalInverseType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type> > & AdditiveSchwarz<MatrixType,LocalInverseType>::getRangeMap() const 
{
  return Matrix_->getRangeMap();
}

//==============================================================================
template<class MatrixType,class LocalInverseType>  
Teuchos::RCP<const Tpetra::RowMatrix<typename MatrixType::scalar_type, typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type> > AdditiveSchwarz<MatrixType,LocalInverseType>::getMatrix() const
{
  return Matrix_;
}
  
//==============================================================================
// Applies the effect of the preconditione.
template<class MatrixType,class LocalInverseType>
void AdditiveSchwarz<MatrixType,LocalInverseType>::apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
			    Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
			    Teuchos::ETransp mode,
			    Scalar alpha,
			    Scalar beta) const
{
  typedef typename Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MultiVectorType;

  TEUCHOS_TEST_FOR_EXCEPTION(IsComputed_ == false, std::runtime_error,
     "Ifpack2::AdditiveSchwarz::apply ERROR: isComputed() must be true prior to calling apply.");

  TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
     "Ifpack2::AdditiveSchwarz::apply ERROR: X.getNumVectors() != Y.getNumVectors().");

  size_t NumVectors = X.getNumVectors();

  Time_->start();

  Teuchos::RCP<MultiVectorType> OverlappingX,OverlappingY,Xtmp;
  
  if(IsOverlapping_){
    // Setup if we're overlapping
    OverlappingX = Teuchos::rcp( new MultiVectorType(OverlappingMatrix_->getRowMap(), X.getNumVectors()) );
    OverlappingY = Teuchos::rcp( new MultiVectorType(OverlappingMatrix_->getRowMap(), X.getNumVectors()) );
    OverlappingY->putScalar(0.0);
    OverlappingX->putScalar(0.0);    
    OverlappingMatrix_->importMultiVector(X,*OverlappingX,Tpetra::INSERT);
    // FIXME from Ifpack1: Will not work with non-zero starting solutions.
  }
  else{
    Xtmp=Teuchos::rcp(new MultiVectorType(X));
    OverlappingX=Xtmp;
    OverlappingY=Teuchos::rcp(&Y,false);		      
  }

  if (FilterSingletons_) {
    // process singleton filter
    MultiVectorType ReducedX(SingletonMatrix_->getRowMap(),NumVectors);
    MultiVectorType ReducedY(SingletonMatrix_->getRowMap(),NumVectors);
    SingletonMatrix_->SolveSingletons(*OverlappingX,*OverlappingY);
    SingletonMatrix_->CreateReducedRHS(*OverlappingY,*OverlappingX,ReducedX);

    // process reordering
    if (!UseReordering_) {
      Inverse_->apply(ReducedX,ReducedY);
    }
    else {
      MultiVectorType ReorderedX(ReducedX);
      MultiVectorType ReorderedY(ReducedY);
      ReorderedLocalizedMatrix_->permuteOriginalToReordered(ReducedX,ReorderedX);
      Inverse_->apply(ReorderedX,ReorderedY);
      ReorderedLocalizedMatrix_->permuteReorderedToOriginal(ReorderedY,ReducedY);
    }

    // finish up with singletons
    SingletonMatrix_->UpdateLHS(ReducedY,*OverlappingY);
  }
  else {

    // process reordering
    if (!UseReordering_) {
      Inverse_->apply(*OverlappingX,*OverlappingY);
    }
    else {
      MultiVectorType ReorderedX(*OverlappingX);
      MultiVectorType ReorderedY(*OverlappingY);
      ReorderedLocalizedMatrix_->permuteOriginalToReordered(*OverlappingX,ReorderedX);
      Inverse_->apply(ReorderedX,ReorderedY);
      ReorderedLocalizedMatrix_->permuteReorderedToOriginal(ReorderedY,*OverlappingY);
    }
  }

  if(IsOverlapping_)
    OverlappingMatrix_->exportMultiVector(*OverlappingY,Y,CombineMode_);

  NumApply_++;
  ApplyTime_ += Time_->stop();
}

//==============================================================================
// Sets all parameters for the preconditioner.
template<class MatrixType,class LocalInverseType>
void AdditiveSchwarz<MatrixType,LocalInverseType>::setParameters(const Teuchos::ParameterList& List_in)
{
  List_=List_in;

  // compute the condition number each time compute() is invoked.
  ComputeCondest_ = List_.get("schwarz: compute condest", ComputeCondest_);
  // combine mode
  if( Teuchos::ParameterEntry *combineModeEntry = List_.getEntryPtr("schwarz: combine mode") ) {
    if( typeid(std::string) == combineModeEntry->getAny().type() ) {
      std::string mode = List_.get("schwarz: combine mode", "Add");
      if (mode == "Add")
        CombineMode_ = Tpetra::ADD;
      else if (mode == "Insert")
        CombineMode_ = Tpetra::INSERT;
      else if (mode == "Replace")
        CombineMode_ = Tpetra::REPLACE;
      else if (mode == "AbsMax")
        CombineMode_ = Tpetra::ABSMAX;
      else {	
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error
				   ,"Error, The Tpetra combine mode of \""<<mode<<"\" is not valid!  Only the values"
				   " \"Add\", \"Insert\", \"Replace\", and \"AbsMax\" are accepted!"
				   );
      }
    }
    else if ( typeid(Tpetra::CombineMode) == combineModeEntry->getAny().type() ) {
      CombineMode_ = Teuchos::any_cast<Tpetra::CombineMode>(combineModeEntry->getAny());
    }
    else {
      // Throw exception with good error message!
      Teuchos::getParameter<std::string>(List_,"schwarz: combine mode");
    }
  }
  else {
    // Make the default be a string to be consistent with the valid parameters!
    List_.get("schwarz: combine mode","Add");
  }

  Ifpack2::getParameter(List_, "schwarz: overlap level",OverlapLevel_);  
  if(OverlapLevel_>0) {
    IsOverlapping_=true;
  }

  // Will we be doing reordering?
  // Note: Unlike Ifpack we'll use a "schwarz: reordering list" to give to Zoltan2...
  UseReordering_ = List_.get("schwarz: use reordering",false);

#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
  // If we don't have Zoltan2, we just turn the reordering off completely...
  UseReordering_=false;
#endif

  // Subdomain check
  if (List_.isParameter("schwarz: subdomain id") && List_.get("schwarz: subdomain id",-1) >  0)
    UseSubdomain_=true;

  // if true, filter singletons. NOTE: the filtered matrix can still have
  // singletons! A simple example: upper triangular matrix, if I remove
  // the lower node, I still get a matrix with a singleton! However, filter
  // singletons should help for PDE problems with Dirichlet BCs.
  FilterSingletons_ = List_.get("schwarz: filter singletons", FilterSingletons_);
}

  
//==============================================================================
// Computes all (graph-related) data necessary to initialize the preconditioner.
template<class MatrixType,class LocalInverseType>
void AdditiveSchwarz<MatrixType,LocalInverseType>::initialize()
{
  IsInitialized_ = false;
  IsComputed_ = false; // values required
  Condest_ = -1.0; // zero-out condest

  if (Time_ == Teuchos::null)
    Time_ = Teuchos::rcp( new Teuchos::Time("Ifpack2::AdditiveSchwarz"));

  Time_->start();
  
  // compute the overlapping matrix if necessary
  if (IsOverlapping_) {
    if(UseSubdomain_) {
      int sid = List_.get("subdomain id",-1);
      OverlappingMatrix_ = Teuchos::rcp( new Ifpack2::OverlappingRowMatrix<LocalMatrixType>(Matrix_, OverlapLevel_, sid) );
    }
    else
      OverlappingMatrix_ = Teuchos::rcp( new Ifpack2::OverlappingRowMatrix<LocalMatrixType>(Matrix_, OverlapLevel_) );

    TEUCHOS_TEST_FOR_EXCEPTION(OverlappingMatrix_ == Teuchos::null, std::runtime_error,
     "Ifpack2::AdditiveSchwarz::initialize ERROR: OverlappingRowMatrix constructor failed.");
  }

  // Setup
  setup();

  // Sanity Checks
  TEUCHOS_TEST_FOR_EXCEPTION(Inverse_ == Teuchos::null, std::runtime_error,
     "Ifpack2::AdditiveSchwarz::initialize ERROR: Setup() failed.");
 
  // Initialize inverse
  Inverse_->setParameters(List_);
  Inverse_->initialize();

  IsInitialized_ = true;
  NumInitialize_++;
  InitializeTime_ += Time_->stop();
}
  
//==============================================================================
// Returns true if the  preconditioner has been successfully initialized, false otherwise.
template<class MatrixType,class LocalInverseType>
bool AdditiveSchwarz<MatrixType,LocalInverseType>::isInitialized() const
{
  return IsInitialized_;
}
  
//==============================================================================
// Computes all (coefficient) data necessary to apply the preconditioner.
template<class MatrixType,class LocalInverseType>
void AdditiveSchwarz<MatrixType,LocalInverseType>::compute()
{
  if (!IsInitialized_) initialize();

  Time_->start();
  IsComputed_ = false;
  Condest_ = -1.0;

  Inverse_->compute();
  
  IsComputed_ = true; 
  ++NumCompute_;
  ComputeTime_ += Time_->stop();

}
  
//==============================================================================
// Returns true if the  preconditioner has been successfully computed, false otherwise.
template<class MatrixType,class LocalInverseType>
bool AdditiveSchwarz<MatrixType,LocalInverseType>::isComputed() const
{
  return IsComputed_;
}
  
//==============================================================================
// Computes the condition number estimate and returns its value.
template<class MatrixType,class LocalInverseType>
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType 
AdditiveSchwarz<MatrixType,LocalInverseType>::computeCondEst(CondestType CT,
							     LocalOrdinal MaxIters,
							     magnitudeType Tol,
							     const Teuchos::Ptr<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &Matrix_in)
{
  
  // If we haven't computed, we can't do a condest
  if (!isComputed()) return (-Teuchos::ScalarTraits<magnitudeType>::one());

  Condest_ = Ifpack2::Condest(*this, CT, MaxIters, Tol, Matrix_in);
  return(Condest_);
}


//==============================================================================
// Returns the computed condition number estimate, or -1.0 if not computed.
template<class MatrixType,class LocalInverseType>
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType 
AdditiveSchwarz<MatrixType,LocalInverseType>::getCondEst() const
{
  return Condest_;
}
  
//==============================================================================
// Returns the number of calls to initialize().
template<class MatrixType,class LocalInverseType>  
int AdditiveSchwarz<MatrixType,LocalInverseType>::getNumInitialize() const 
{
  return NumInitialize_;
}
  
//==============================================================================
// Returns the number of calls to compute().
template<class MatrixType,class LocalInverseType>
int AdditiveSchwarz<MatrixType,LocalInverseType>::getNumCompute() const
{
  return NumCompute_;
}
  
//==============================================================================
// Returns the number of calls to apply().
template<class MatrixType,class LocalInverseType>
int AdditiveSchwarz<MatrixType,LocalInverseType>::getNumApply() const 
{
  return NumApply_;
}
  
//==============================================================================
// Returns the time spent in initialize().
template<class MatrixType,class LocalInverseType>
double AdditiveSchwarz<MatrixType,LocalInverseType>::getInitializeTime() const
{
  return InitializeTime_;
}
  
//==============================================================================
// Returns the time spent in compute().
template<class MatrixType,class LocalInverseType>
double AdditiveSchwarz<MatrixType,LocalInverseType>::getComputeTime() const
{
  return ComputeTime_;
}

//==============================================================================
// Returns the time spent in Apply().
template<class MatrixType,class LocalInverseType>
double AdditiveSchwarz<MatrixType,LocalInverseType>::getApplyTime() const 
{
  return ApplyTime_;
}
//==============================================================================
/* Return a simple one-line description of this object. */
template<class MatrixType,class LocalInverseType>
std::string AdditiveSchwarz<MatrixType,LocalInverseType>::description() const
{
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
  oss<<"overlap level ="<<OverlapLevel_;
  oss << "}";
  return oss.str();
}

//==============================================================================
/* Print the object with some verbosity level to an FancyOStream object. */
template<class MatrixType,class LocalInverseType>
void AdditiveSchwarz<MatrixType,LocalInverseType>::describe(Teuchos::FancyOStream &os, const Teuchos::EVerbosityLevel verbLevel) const
{
  using std::endl;

  os << endl;
  os << "================================================================================" << endl;
  os << "Ifpack2::AdditiveSchwarz, overlap level = " << OverlapLevel_ << endl;
  if (CombineMode_ == Tpetra::INSERT)
    os << "Combine mode                          = Insert" << endl;
  else if (CombineMode_ == Tpetra::ADD)
    os << "Combine mode                          = Add" << endl;
  else if (CombineMode_ == Tpetra::REPLACE)
    os << "Combine mode                          = Replace" << endl;
  else if (CombineMode_ == Tpetra::ABSMAX)
    os << "Combine mode                          = AbsMax" << endl;

  os << "Condition number estimate             = " << Condest_ << endl;
  os << "Global number of rows                 = " << Matrix_->getGlobalNumRows() << endl;

  os << endl;
  os << "Phase           # calls   Total Time (s)       Total MFlops     MFlops/s" << endl;
  os << "-----           -------   --------------       ------------     --------" << endl;
  os << "Initialize()    "   << std::setw(5) << getNumInitialize()
     << "  " << std::setw(15) << getInitializeTime() << endl;
  os << "Compute()       "   << std::setw(5) << getNumCompute() 
     << "  " << std::setw(15) << getComputeTime() << endl;
  os << "Apply()         "   << std::setw(5) << getNumApply() 
     << "  " << std::setw(15) << getApplyTime() <<endl;
  os << "================================================================================" << endl;
  os << endl;
}

//==============================================================================
// Prints basic information on iostream. This function is used by operator<<.
template<class MatrixType,class LocalInverseType>
std::ostream& AdditiveSchwarz<MatrixType,LocalInverseType>::print(std::ostream& os) const
{
  Teuchos::FancyOStream fos(Teuchos::rcp(&os,false));
  fos.setOutputToRootOnly(0);
  describe(fos);
  return(os);
}

//==============================================================================
// Returns the level of overlap.
template<class MatrixType,class LocalInverseType>
int AdditiveSchwarz<MatrixType,LocalInverseType>::getOverlapLevel() const
{
  return OverlapLevel_;
}

//==============================================================================
// Sets up the localized matrix and the singleton filter.
template<class MatrixType,class LocalInverseType>
void AdditiveSchwarz<MatrixType,LocalInverseType>::setup() 
{  
#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
  typedef typename Xpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> XpetraMatrixType;
  typedef typename Xpetra::TpetraRowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> XpetraTpetraMatrixType;
#endif

  Teuchos::RCP<LocalMatrixType> ActiveMatrix;

  // Create Localized Matrix
  if (OverlappingMatrix_ != Teuchos::null) {
    if (UseSubdomain_) {
      //      int sid = List_.get("subdomain id",-1);
      throw std::runtime_error("Ifpack2::AdditiveSchwarz subdomain code not yet supported.");
      //Ifpack2_NodeFilter *tt = new Ifpack2_NodeFilter(OverlappingMatrix_,nodeID); //FIXME
      //LocalizedMatrix_ = Teuchos::rcp(tt);
    }
    else
      LocalizedMatrix_ = Teuchos::rcp( new Ifpack2::LocalFilter<LocalMatrixType>(OverlappingMatrix_) );
  }
  else {
    if (UseSubdomain_) {
      //      int sid = List_.get("subdomain id",-1);
      throw std::runtime_error("Ifpack2::AdditiveSchwarz subdomain code not yet supported.");
    }
    else
      LocalizedMatrix_ = Teuchos::rcp( new Ifpack2::LocalFilter<LocalMatrixType>(Matrix_) );
  }

  // Mark localized matrix as active
  ActiveMatrix = LocalizedMatrix_;
  TEUCHOS_TEST_FOR_EXCEPTION(LocalizedMatrix_ == Teuchos::null, std::runtime_error,
			     "Ifpack2::AdditiveSchwarz::Setup() ERROR: LocalFilter failed.");

  // Singleton Filtering
  if (FilterSingletons_) {
    SingletonMatrix_ = Teuchos::rcp( new Ifpack2::SingletonFilter<LocalMatrixType>(LocalizedMatrix_) );
    ActiveMatrix = SingletonMatrix_;
  }


  // Do reordering
  if (UseReordering_) {
#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
    // Unlike Ifpack, Zoltan2 does all the dirty work here.
    Teuchos::ParameterList zlist = List_.sublist("schwarz: reordering list");
    XpetraTpetraMatrixType XpetraMatrix(ActiveMatrix);
    Zoltan2::XpetraRowMatrixInput<XpetraMatrixType> Zoltan2Matrix(Teuchos::rcp<XpetraMatrixType>(&XpetraMatrix,false));

#ifdef HAVE_MPI
    // Grab the MPI Communicator and build the ordering problem with that
    MPI_Comm mycomm;

    Teuchos::RCP<const Teuchos::MpiComm<int> > mpicomm = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(ActiveMatrix->getComm());
    if(mpicomm == Teuchos::null) mycomm = MPI_COMM_SELF;
    else mycomm = *(mpicomm->getRawMpiComm());
    Zoltan2::OrderingProblem<Zoltan2::XpetraRowMatrixInput<XpetraMatrixType> > MyOrderingProblem(&Zoltan2Matrix, &zlist,mycomm);    
#else
    Zoltan2::OrderingProblem<Zoltan2::XpetraRowMatrixInput<XpetraMatrixType> > MyOrderingProblem(&Zoltan2Matrix, &zlist);    
#endif
    MyOrderingProblem.solve();

    // Now create the reordered matrix & mark it as active
    ReorderedLocalizedMatrix_ =  Teuchos::rcp(new Ifpack2::ReorderFilter<LocalMatrixType>(ActiveMatrix,Teuchos::rcp( new Zoltan2::OrderingSolution<GlobalOrdinal,LocalOrdinal>(*MyOrderingProblem.getSolution()))));

    ActiveMatrix = ReorderedLocalizedMatrix_;
#else
    throw std::runtime_error("Ifpack2::AdditiveSchwarz::Setup() You need to compile with Zoltan2 to support reordering.");
#endif
  }

  // Build the inverse
  Inverse_ = Teuchos::rcp(new LocalInverseType(ActiveMatrix));
  TEUCHOS_TEST_FOR_EXCEPTION(Inverse_ == Teuchos::null, std::runtime_error,
			     "Ifpack2::AdditiveSchwarz::Setup() ERROR: Inverse constructor failed.");
}


}// namespace Ifpack2

#endif // IFPACK2_ADDITIVESCHWARZ_DECL_HPP
