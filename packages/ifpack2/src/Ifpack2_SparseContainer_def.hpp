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

#ifndef IFPACK2_SPARSECONTAINER_DEF_HPP
#define IFPACK2_SPARSECONTAINER_DEF_HPP

#include "Ifpack2_SparseContainer_decl.hpp"
#ifdef HAVE_MPI
#include <mpi.h>
#include "Teuchos_DefaultMpiComm.hpp"
#else
#include "Teuchos_DefaultSerialComm.hpp"
#endif

namespace Ifpack2 {

//==============================================================================
template<class MatrixType, class InverseType>
SparseContainer<MatrixType,InverseType>::SparseContainer(const size_t NumRows, const size_t NumVectors) :
  NumRows_(NumRows),
  NumVectors_(NumVectors),
  IsInitialized_(false),
  IsComputed_(false)
{
#ifdef HAVE_MPI
  LocalComm_ = Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper((MPI_Comm)MPI_COMM_SELF)));
#else
  LocalComm_ = Teuchos::rcp( new Teuchos::SerialComm<int>() );
#endif
}

//==============================================================================
template<class MatrixType, class InverseType>
SparseContainer<MatrixType,InverseType>::SparseContainer(const SparseContainer<MatrixType,InverseType>& rhs)
{
  // nobody should ever call this  
}

//==============================================================================
template<class MatrixType, class InverseType>
SparseContainer<MatrixType,InverseType>::~SparseContainer()
{
  destroy();
}

//==============================================================================
template<class MatrixType, class InverseType>
size_t SparseContainer<MatrixType,InverseType>::getNumRows() const
{
  if (isInitialized()) return NumRows_;
  else return 0;
}

//==============================================================================
// Sets the number of vectors for X/Y
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::setNumVectors(const size_t NumVectors_in)
{
  TEUCHOS_TEST_FOR_EXCEPTION(NumVectors_in<=0, std::runtime_error, "Ifpack2::SparseContainer::setNumVectors NumVectors cannot be negative.");  

  // Do nothing if the vectors are the right size.  Always do something if we're not initialized.
  if(!IsInitialized_  || NumVectors_!=NumVectors_in) {
    NumVectors_=NumVectors_in;
    X_ = Teuchos::rcp( new Tpetra::MultiVector<InverseScalar,InverseLocalOrdinal,InverseGlobalOrdinal,InverseNode>(Map_,NumVectors_) );
    Y_ = Teuchos::rcp( new Tpetra::MultiVector<InverseScalar,InverseLocalOrdinal,InverseGlobalOrdinal,InverseNode>(Map_,NumVectors_) );
  }

}

//==============================================================================
// Returns the ID associated to local row i. 
template<class MatrixType, class InverseType>
typename MatrixType::local_ordinal_type & SparseContainer<MatrixType,InverseType>::ID(const size_t i)
{
  return GID_[i];
}

//==============================================================================
// Returns  true is the container has been successfully initialized.
template<class MatrixType, class InverseType>
bool SparseContainer<MatrixType,InverseType>::isInitialized() const
{
  return IsInitialized_;
}

//==============================================================================
// Returns  true is the container has been successfully computed.
template<class MatrixType, class InverseType>
bool SparseContainer<MatrixType,InverseType>::isComputed() const
{
  return IsComputed_;
}

//==============================================================================
// Sets all necessary parameters.
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::setParameters(const Teuchos::ParameterList& List)
{
  List_ = List;
}
 
//==============================================================================
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::initialize()
{
  if(IsInitialized_) destroy();
  IsInitialized_=false;

  // Allocations
  Map_ = Teuchos::rcp( new Tpetra::Map<InverseLocalOrdinal,InverseGlobalOrdinal,InverseNode>(NumRows_,0,LocalComm_) );
  Matrix_ = Teuchos::rcp( new Tpetra::CrsMatrix<InverseScalar,InverseLocalOrdinal,InverseGlobalOrdinal,InverseNode>(Map_,0) );
  GID_.resize(NumRows_);  
  setNumVectors(NumVectors_);

  // create the inverse
  Inverse_ = Teuchos::rcp( new InverseType(Matrix_) );
  TEUCHOS_TEST_FOR_EXCEPTION( Inverse_ == Teuchos::null, std::runtime_error, "Ifpack2::SparseContainer::initialize error in inverse constructor.");  
  Inverse_->setParameters(List_);

  // Note from IFPACK: Call Inverse_->Initialize() in Compute(). This saves
  // some time, because the diagonal blocks can be extracted faster,
  // and only once.
  IsInitialized_ = true;
}

//==============================================================================
// Finalizes the linear system matrix and prepares for the application of the inverse.
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::compute(const Teuchos::RCP<const Tpetra::RowMatrix<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode> >& Matrix)
{
  IsComputed_=false;
  TEUCHOS_TEST_FOR_EXCEPTION( !IsInitialized_, std::runtime_error, "Ifpack2::SparseContainer::compute please call initialize first.");  

  // extract the submatrices
  extract(Matrix);

  // initialize & compute the inverse operator
  Inverse_->initialize();
  Inverse_->compute();
  IsComputed_=true;
}

//==============================================================================
// Apply the inverse of the matrix to RHS, result is stored in LHS.
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::apply(const Tpetra::MultiVector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode>& X,
						    Tpetra::MultiVector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode>& Y,
						    Tpetra::CombineMode mode,
						    MatrixScalar alpha,
						    MatrixScalar beta)
{
 TEUCHOS_TEST_FOR_EXCEPTION( !IsComputed_, std::runtime_error, "Ifpack2::SparseContainer::apply compute has not been called.");  
 TEUCHOS_TEST_FOR_EXCEPTION( X.getNumVectors() != Y.getNumVectors(), std::runtime_error, "Ifpack2::SparseContainer::apply X/Y have different numbers of vectors.");   

 // NTS: Remove this test later
 TEUCHOS_TEST_FOR_EXCEPTION( mode!=Tpetra::REPLACE, std::runtime_error, "Ifpack2::SparseContainer::apply combine modes besides Tpetra:REPLACE not yet supported.");   

 // This is a noop if our internal vectors are already the right size.
 setNumVectors(X.getNumVectors());

 // Pointers, pointers, pointers!
 Teuchos::ArrayRCP<Teuchos::ArrayRCP<const MatrixScalar> >  global_x_ptr = X.get2dView();
 Teuchos::ArrayRCP<Teuchos::ArrayRCP<InverseScalar> >       local_x_ptr  = X_->get2dViewNonConst();
 Teuchos::ArrayRCP<Teuchos::ArrayRCP<const InverseScalar> > local_y_ptr  = Y_->get2dView();
 Teuchos::ArrayRCP<Teuchos::ArrayRCP<MatrixScalar> >        global_y_ptr = Y.get2dViewNonConst();

 // Copy in
 for(size_t j=0; j < NumRows_; j++){
   MatrixLocalOrdinal LID = ID(j);
   for (size_t k = 0 ; k < NumVectors_ ; k++) 
     local_x_ptr[k][j] = global_x_ptr[k][LID];   
 }

 // Apply inverse
 Inverse_->apply(*X_, *Y_);
 
 // copy out into solution vector Y
 for (size_t j = 0 ; j < NumRows_ ; j++) {
   MatrixLocalOrdinal LID = ID(j);
   for (size_t k = 0 ; k < NumVectors_ ; k++) 
     global_y_ptr[k][LID] = alpha * local_y_ptr[k][j] + beta * global_y_ptr[k][LID];
 }

}

//==============================================================================
// Destroys all data.
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::destroy()
{
  IsInitialized_ = false;
  IsComputed_ = false;
}

//============================================================================== 
// Prints basic information on iostream. This function is used by operator<<
template<class MatrixType, class InverseType>
std::ostream& SparseContainer<MatrixType,InverseType>::print(std::ostream& os) const
{
  Teuchos::FancyOStream fos(Teuchos::rcp(&os,false));
  fos.setOutputToRootOnly(0);
  describe(fos);
  return(os);
}


//==============================================================================
template<class MatrixType, class InverseType>
std::string SparseContainer<MatrixType,InverseType>::description() const
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

  oss << "}";
  return oss.str();
}

//==============================================================================
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::describe(Teuchos::FancyOStream &os, const Teuchos::EVerbosityLevel verbLevel) const
{
  using std::endl;
  if(verbLevel==Teuchos::VERB_NONE) return;
  os << "================================================================================" << endl;
  os << "Ifpack2_SparseContainer" << endl;
  os << "Number of rows          = " << NumRows_ << endl;
  os << "Number of vectors       = " << NumVectors_ << endl;
  os << "isInitialized()         = " << IsInitialized_ << endl;
  os << "isComputed()            = " << IsComputed_ << endl;
  os << "================================================================================" << endl;
  os << endl;
}

//============================================================================== 
// Extract the submatrices identified by the ID set int ID().
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::extract(const Teuchos::RCP<const Tpetra::RowMatrix<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode> >& Matrix_in) 
{
  size_t MatrixInNumRows= Matrix_in->getNodeNumRows();

  // Sanity checking
  for (size_t j = 0 ; j < NumRows_ ; ++j) {
    TEUCHOS_TEST_FOR_EXCEPTION( GID_[j] < 0 || (size_t) GID_[j] >= MatrixInNumRows, std::runtime_error, "Ifpack2::SparseContainer::applyInverse compute has not been called.");  
  }  

  int Length = Matrix_in->getNodeMaxNumRowEntries();
  Teuchos::Array<MatrixScalar>         Values;
  Teuchos::Array<MatrixLocalOrdinal>   Indices;
  Teuchos::Array<InverseScalar>        Values_insert;
  Teuchos::Array<InverseGlobalOrdinal> Indices_insert;

  Values.resize(Length);
  Indices.resize(Length);
  Values_insert.resize(Length);
  Indices_insert.resize(Length);

  for (size_t j = 0 ; j < NumRows_ ; ++j) {
    MatrixLocalOrdinal LRID = ID(j);
    size_t NumEntries;

    Matrix_in->getLocalRowCopy(LRID,Indices(),Values(),NumEntries);

    size_t num_entries_found=0;
    for (size_t k = 0 ; k < NumEntries ; ++k) {
      MatrixLocalOrdinal LCID = Indices[k];

      // skip off-processor elements
      if ((size_t)LCID >= MatrixInNumRows) 
	continue;

      // for local column IDs, look for each ID in the list
      // of columns hosted by this object
      // FIXME: use STL
      InverseLocalOrdinal jj = -1;
      for (size_t kk = 0 ; kk < NumRows_ ; ++kk)
	if (ID(kk) == LCID)
	  jj = kk;

      if (jj != -1) {
	Indices_insert[num_entries_found] = Map_->getGlobalElement(jj);
	Values_insert[num_entries_found]  = Values[k];
	num_entries_found++;
      }

    }
    Matrix_->insertGlobalValues(j,Indices_insert(0,num_entries_found),Values_insert(0,num_entries_found));
  }

  Matrix_->fillComplete();
}


} // namespace Ifpack2
#endif // IFPACK2_SPARSECONTAINER_HPP
