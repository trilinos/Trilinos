//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
//@HEADER


/*! \file BelosEpetraAdapter.cpp
    \brief Implementation of the interfaces between Belos virtual classes and Epetra concrete classes.
*/

#include "BelosEpetraAdapter.hpp"
#include "BelosInnerSolver.hpp"

using namespace Belos;

namespace Belos {

  bool 
  VectorSpaceTraits<Epetra_BlockMap>::
  compatible (const VectorSpaceTraits<Epetra_BlockMap>::vector_space_type& first, 
	      const VectorSpaceTraits<Epetra_BlockMap>::vector_space_type& second)
  {
    // FIXME (mfh 07 Mar 2011) Should I call PointSameAs() or
    // SameAs() here?
    //
    // FIXME (mfh 07 Mar 2011) Does it suffice to compare pointers?
    // If the two pointers point to the same object on one MPI
    // process, they should on all MPI processes, right?  Comparing
    // pointers avoids communication for the idiomatic Tpetra case of
    // multiple objects constructed from the same Map object, but
    // Epetra objects that inherit from Epetra_DistObject make a deep
    // copy of the Epetra_BlockMap, so comparing pointers in the
    // Epetra world may make less sense.
    return &first == &second || first.PointSameAs (second);
  }

  Teuchos::RCP<const VectorSpaceTraits<Epetra_BlockMap>::vector_space_type>
  VectorSpaceTraits<Epetra_BlockMap>::
  persistentView (const Teuchos::RCP<const VectorSpaceTraits<Epetra_BlockMap>::vector_space_type>& space)
  {
    using Teuchos::is_null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;
    using Teuchos::rcp_implicit_cast;

    // The argument might be a weak reference (Epetra objects'
    // methods typically return a "const Epetra_{Block}Map&").  If
    // it is, make a deep copy, using the constructor appropriate to
    // the (dynamic) type of space.  Otherwise, just return the RCP;
    // in that case, it's either a strong reference or null.
    if (space.strength() == Teuchos::RCP_WEAK)
      {
	// If space is an Epetra_Map, copy space using Epetra_Map's
	// constructor.  Otherwise, copy space using
	// Epetra_BlockMap's constructor.
	RCP<const Epetra_Map> mapMap = rcp_dynamic_cast<const Epetra_Map> (space);
	if (is_null (mapMap))
	  return rcp (new Epetra_BlockMap (*space));
	else
	  return rcp_implicit_cast<const Epetra_BlockMap> (rcp (new Epetra_Map (*mapMap)));
      }
    else // either a strong reference or null.
      return space;
  }

  Teuchos::RCP<const VectorSpaceTraits<Epetra_BlockMap>::vector_space_type>
  VectorSpaceTraits<Epetra_BlockMap>::
  persistentView (const VectorSpaceTraits<Epetra_BlockMap>::vector_space_type& space)
  {
    return persistentView (Teuchos::rcpFromRef (space));
  }

  Teuchos::RCP<const MultiVecTraits<double, Epetra_MultiVector>::vector_space_type> 
  MultiVecTraits<double, Epetra_MultiVector>::
  getRange (const Epetra_MultiVector& x) 
  {
    return VectorSpaceTraits<vector_space_type>::persistentView (x.Map());
  }

  Teuchos::RCP<const OperatorTraits<double, Epetra_MultiVector, Epetra_Operator>::vector_space_type>
  OperatorTraits<double, Epetra_MultiVector, Epetra_Operator>::
  getDomain (const Epetra_Operator& A) 
  {
    return VectorSpaceTraits<vector_space_type>::persistentView (A.OperatorDomainMap());
  }

  Teuchos::RCP<const OperatorTraits<double, Epetra_MultiVector, Epetra_Operator>::vector_space_type>
  OperatorTraits<double, Epetra_MultiVector, Epetra_Operator>::
  getRange (const Epetra_Operator& A) 
  {
    return VectorSpaceTraits<vector_space_type>::persistentView (A.OperatorDomainMap());
  }


  EpetraInnerSolver::
  EpetraInnerSolver (const Teuchos::RCP<inner_solver_type>& solver) :
    solver_ (solver), 
    domain_ (solver->getDomain()),
    range_ (solver->getRange())
  {}

  int
  EpetraInnerSolver::
  Apply (const Epetra_MultiVector& X,
	 Epetra_MultiVector& Y) const
  {
    using Teuchos::RCP;
    using Teuchos::rcpFromRef;
    typedef multivector_type MV;
    typedef Teuchos::ScalarTraits<scalar_type> STS;

    // FIXME (mfh 08 Mar 2011) There are many ways in which X and Y
    // could alias one another; for example, X could be a const column
    // view of Y.  We choose only to check for the simple case in
    // which the two pointers are identical.
    TEST_FOR_EXCEPTION(&X == &Y, std::invalid_argument,
		       "EpetraInnerSolver::Apply() does not allow aliasing of "
		       "X and Y.");
    // "X" is the right-hand side in this case, and "Y" is the
    // "left-hand side."
    RCP<const MV> X_ptr = rcpFromRef (X);
    RCP<MV> Y_ptr = rcpFromRef (Y);
    solver_->solve (Y_ptr, X_ptr);
    return 0;
  }

  int 
  EpetraInnerSolver::
  ApplyInverse (const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
  {
    return -1; // Never successful
  }

  int 
  EpetraInnerSolver::SetUseTranspose (bool UseTranspose) {
    (void) UseTranspose; // Silence compiler warning
    return -1; // Applying the transpose is not supported
  }

  double 
  EpetraInnerSolver::NormInf() const {
    throw std::logic_error("There is no reasonable way to compute the "
			   "infinity norm of an Epetra_Operator implemented "
			   "using InnerSolver.");
  }

  const char* 
  EpetraInnerSolver::Label() const {
    return "Epetra InnerSolver adaptor";
  }

  const Epetra_Map& 
  EpetraInnerSolver::OperatorDomainMap() const 
  {
    using Teuchos::rcp_dynamic_cast;
    // If it's not an Epetra_Map, we want the dynamic cast to fail
    // immediately with an informative error message.
    return *rcp_dynamic_cast<const Epetra_Map> (domain_, true);
  }

  const Epetra_Map& 
  EpetraInnerSolver::OperatorRangeMap() const 
  {
    using Teuchos::rcp_dynamic_cast;
    // If it's not an Epetra_Map, we want the dynamic cast to fail
    // immediately with an informative error message.
    return *rcp_dynamic_cast<const Epetra_Map> (range_, true);
  }

} // namespace Belos
  
  //-------------------------------------------------------------
  
  //////////////////////////////////////////////////////////////////////
  // Construction/Destruction
  //////////////////////////////////////////////////////////////////////
  
  
EpetraMultiVec::EpetraMultiVec(const Epetra_BlockMap& Map_in, double * array, 
			       const int numvecs, const int stride)
  : Epetra_MultiVector(Copy, Map_in, array, stride, numvecs) 
{
}


EpetraMultiVec::EpetraMultiVec(const Epetra_BlockMap& Map_in, const int numvecs, bool zeroOut)
  : Epetra_MultiVector(Map_in, numvecs, zeroOut) 
{
}


EpetraMultiVec::EpetraMultiVec(Epetra_DataAccess CV_in, const Epetra_MultiVector& P_vec, 				
			       const std::vector<int>& index )
  : Epetra_MultiVector(CV_in, P_vec, &(const_cast<std::vector<int> &>(index))[0], index.size())
{
}


EpetraMultiVec::EpetraMultiVec(const Epetra_MultiVector& P_vec)
  : Epetra_MultiVector(P_vec) 
{
}


EpetraMultiVec::~EpetraMultiVec() 
{
}
//
//  member functions inherited from Belos::MultiVec
//
//
//  Simulating a virtual copy constructor. If we could rely on the co-variance
//  of virtual functions, we could return a pointer to EpetraMultiVec
//  (the derived type) instead of a pointer to the pure virtual base class.
//

MultiVec<double>* EpetraMultiVec::Clone ( const int numvecs ) const
{
  EpetraMultiVec * ptr_apv = new EpetraMultiVec(Map(), numvecs, false);
  return ptr_apv; // safe upcast.
}
//
//  the following is a virtual copy constructor returning
//  a pointer to the pure virtual class. std::vector values are
//  copied.
//

MultiVec<double>* EpetraMultiVec::CloneCopy() const
{
  EpetraMultiVec *ptr_apv = new EpetraMultiVec(*this);
  return ptr_apv; // safe upcast
}


MultiVec<double>* EpetraMultiVec::CloneCopy ( const std::vector<int>& index ) const
{
  EpetraMultiVec * ptr_apv = new EpetraMultiVec(Copy, *this, index);
  return ptr_apv; // safe upcast.
}


MultiVec<double>* EpetraMultiVec::CloneViewNonConst ( const std::vector<int>& index ) 
{
  EpetraMultiVec * ptr_apv = new EpetraMultiVec(View, *this, index);
  return ptr_apv; // safe upcast.
}
  

const MultiVec<double>* EpetraMultiVec::CloneView ( const std::vector<int>& index ) const
{
  EpetraMultiVec * ptr_apv = new EpetraMultiVec(View, *this, index);
  return ptr_apv; // safe upcast.
}
  

void EpetraMultiVec::SetBlock( const MultiVec<double>& A, const std::vector<int>& index ) 
{	
  EpetraMultiVec temp_vec(View, *this, index);
  
  int numvecs = index.size();
  if ( A.GetNumberVecs() != numvecs ) {
    std::vector<int> index2( numvecs );
    for(int i=0; i<numvecs; i++)
      index2[i] = i;
    EpetraMultiVec *tmp_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); 
    TEST_FOR_EXCEPTION(tmp_vec==NULL, EpetraMultiVecFailure,
                       "Belos::EpetraMultiVec::SetBlock cast from Belos::MultiVec<> to Belos::EpetraMultiVec failed.");
    EpetraMultiVec A_vec(View, *tmp_vec, index2);
    temp_vec.MvAddMv( 1.0, A_vec, 0.0, A_vec );
  }
  else {
    temp_vec.MvAddMv( 1.0, A, 0.0, A );
  }
}								

//-------------------------------------------------------------
//
// *this <- alpha * A * B + beta * (*this)
//
//-------------------------------------------------------------

void EpetraMultiVec::MvTimesMatAddMv ( const double alpha, const MultiVec<double>& A, 
				       const Teuchos::SerialDenseMatrix<int,double>& B, const double beta ) 
{
  Epetra_LocalMap LocalMap(B.numRows(), 0, Map().Comm());
  Epetra_MultiVector B_Pvec(View, LocalMap, B.values(), B.stride(), B.numCols());
  
  EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); 
  TEST_FOR_EXCEPTION(A_vec==NULL, EpetraMultiVecFailure,
                     "Belos::EpetraMultiVec::MvTimesMatAddMv cast from Belos::MultiVec<> to Belos::EpetraMultiVec failed.");
  
  int info = Multiply( 'N', 'N', alpha, *A_vec, B_Pvec, beta );
  TEST_FOR_EXCEPTION(info!=0, EpetraMultiVecFailure, 
		     "Belos::EpetraMultiVec::MvTimesMatAddMv call to Multiply() returned a nonzero value.");

}

//-------------------------------------------------------------
//
// *this <- alpha * A + beta * B
//
//-------------------------------------------------------------
  
void EpetraMultiVec::MvAddMv ( const double alpha , const MultiVec<double>& A, 
			       const double beta, const MultiVec<double>& B) 
{
  EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A));
  TEST_FOR_EXCEPTION( A_vec==NULL, EpetraMultiVecFailure,
                     "Belos::EpetraMultiVec::MvAddMv cast from Belos::MultiVec<> to Belos::EpetraMultiVec failed.");
  EpetraMultiVec *B_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(B));
  TEST_FOR_EXCEPTION( B_vec==NULL, EpetraMultiVecFailure,
                     "Belos::EpetraMultiVec::MvAddMv cast from Belos::MultiVec<> to Belos::EpetraMultiVec failed.");
  
  int info = Update( alpha, *A_vec, beta, *B_vec, 0.0 );
  TEST_FOR_EXCEPTION(info!=0, EpetraMultiVecFailure, 
		     "Belos::EpetraMultiVec::MvAddMv call to Update() returned a nonzero value.");
}

//-------------------------------------------------------------
//
// this[i] = alpha[i] * this[i]
//
//-------------------------------------------------------------
void EpetraMultiVec::MvScale ( const std::vector<double>& alpha )
{
  // Check to make sure the vector is as long as the multivector has columns.
  int numvecs = this->NumVectors();
  TEST_FOR_EXCEPTION((int)alpha.size() != numvecs, EpetraMultiVecFailure, 
			 "Belos::MultiVecTraits<double,Epetra_MultiVec>::MvScale scaling vector (alpha) not same size as number of input vectors (mv).");
  int ret = 0;
  std::vector<int> tmp_index( 1, 0 );
  for (int i=0; i<numvecs; i++) {
    Epetra_MultiVector temp_vec(View, *this, &tmp_index[0], 1);
    ret = temp_vec.Scale( alpha[i] );
    TEST_FOR_EXCEPTION(ret!=0, EpetraMultiVecFailure, 
                      "Belos::MultiVecTraits<double,Epetra_MultiVec>::MvScale call to Scale() returned a nonzero value.");
    tmp_index[0]++;
  }
}

//-------------------------------------------------------------
//
// dense B <- alpha * A^T * (*this)
//
//-------------------------------------------------------------

void EpetraMultiVec::MvTransMv ( const double alpha, const MultiVec<double>& A,
				 Teuchos::SerialDenseMatrix<int,double>& B) const
{    
  EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A));
  
  if (A_vec) {
    Epetra_LocalMap LocalMap(B.numRows(), 0, Map().Comm());
    Epetra_MultiVector B_Pvec(View, LocalMap, B.values(), B.stride(), B.numCols());
    
    int info = B_Pvec.Multiply( 'T', 'N', alpha, *A_vec, *this, 0.0 );
    TEST_FOR_EXCEPTION(info!=0, EpetraMultiVecFailure, 
		       "Belos::EpetraMultiVec::MvTransMv call to Multiply() returned a nonzero value.");
  }
}

//-------------------------------------------------------------
//
// b[i] = A[i]^T * this[i]
// 
//-------------------------------------------------------------

void EpetraMultiVec::MvDot ( const MultiVec<double>& A, std::vector<double>& b ) const
{
  EpetraMultiVec *A_vec = dynamic_cast<EpetraMultiVec *>(&const_cast<MultiVec<double> &>(A)); 
  TEST_FOR_EXCEPTION(A_vec==NULL, EpetraMultiVecFailure,
                     "Belos::EpetraMultiVec::MvDot cast from Belos::MultiVec<> to Belos::EpetraMultiVec failed.");
  if (A_vec && ( (int)b.size() >= A_vec->NumVectors() ) ) {
     int info = this->Dot( *A_vec, &b[0] );
     TEST_FOR_EXCEPTION(info!=0, EpetraMultiVecFailure, 
			"Belos::EpetraMultiVec::MvDot call to Dot() returned a nonzero value.");   
  }
}

//-------------------------------------------------------------
//
// alpha[i] = norm of i-th column of (*this)
//
//-------------------------------------------------------------

void EpetraMultiVec::MvNorm ( std::vector<double>& normvec, NormType norm_type ) const {
  if ((int)normvec.size() >= GetNumberVecs()) {
    int info = 0;
    switch( norm_type ) {
    case ( OneNorm ) :
      info = Norm1(&normvec[0]);
      break;
    case ( TwoNorm ) :
      info = Norm2(&normvec[0]);
      break;
    case ( InfNorm ) :	
      info = NormInf(&normvec[0]);
      break;
    default:
      break;
    }
    TEST_FOR_EXCEPTION(info!=0, EpetraMultiVecFailure, 
		       "Belos::EpetraMultiVec::MvNorm call to Norm() returned a nonzero value.");
  }
}

///////////////////////////////////////////////////////////////
//
// implementation of the BelosEpetraOp class.
//
////////////////////////////////////////////////////////////////////
//
// BelosOperator constructors
//

EpetraOp::EpetraOp( const Teuchos::RCP<Epetra_Operator> &Op ) 
  : Epetra_Op(Op)
{
}

//
// BelosOperator applications
//
void EpetraOp::Apply ( const MultiVec<double>& x, 
                       MultiVec<double>& y, ETrans trans ) const {
  int info=0;
  MultiVec<double> & temp_x = const_cast<MultiVec<double> &>(x);
  Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
  Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);
  //
  TEST_FOR_EXCEPTION( vec_x==NULL || vec_y==NULL, EpetraOpFailure, 
		      "Belos::EpetraOp::Apply, x and/or y cannot be dynamic cast to an Epetra_MultiVector.");
  //
  // Set the operator to apply the transpose if that is requested.
  //
  if ( trans ) { 
    info = Epetra_Op->SetUseTranspose( true );
  }
  info = Epetra_Op->Apply( *vec_x, *vec_y );
  
  if ( trans ) { 
    info = Epetra_Op->SetUseTranspose( false );
  }
  
  TEST_FOR_EXCEPTION(info!=0, EpetraOpFailure, 
		     "Belos::EpetraOp::Apply call to Apply() returned a nonzero value."); 

}

///////////////////////////////////////////////////////////////
//
// implementation of the BelosEpetraPrecOp class.
//
////////////////////////////////////////////////////////////////////
//
// BelosOperator constructors
//

EpetraPrecOp::EpetraPrecOp( const Teuchos::RCP<Epetra_Operator> &Op ) 
  : Epetra_Op(Op)
{
}

//
// BelosOperator applications
//
void EpetraPrecOp::Apply ( const MultiVec<double>& x, 
			   MultiVec<double>& y, ETrans trans ) const {
  int info=0;
  MultiVec<double> & temp_x = const_cast<MultiVec<double> &>(x);
  Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
  Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);
  //
  TEST_FOR_EXCEPTION( vec_x==NULL || vec_y==NULL, EpetraOpFailure, 
		      "Belos::EpetraOp::Apply, x and/or y cannot be dynamic cast to an Epetra_MultiVector.");
  //
  // Set the operator to apply the transpose if that is requested.
  //
  if ( trans ) { 
    info = Epetra_Op->SetUseTranspose( true );
  }
  info = Epetra_Op->ApplyInverse( *vec_x, *vec_y );
  
  if ( trans ) { 
    info = Epetra_Op->SetUseTranspose( false );
  }
  
  TEST_FOR_EXCEPTION(info!=0, EpetraOpFailure, 
		     "Belos::EpetraOp::Apply call to Apply() returned a nonzero value."); 
}

int EpetraPrecOp::Apply ( const Epetra_MultiVector &X, Epetra_MultiVector &Y ) const 
{
  //
  // This operator applies Y = A^{-1}*X
  //
  int info = Epetra_Op->ApplyInverse( X, Y );
  return info;
}

int EpetraPrecOp::ApplyInverse( const Epetra_MultiVector &X, Epetra_MultiVector &Y ) const
{
  //
  // This operator applies Y = A*X
  //
  int info = Epetra_Op->Apply( X, Y );
  return info;
}

